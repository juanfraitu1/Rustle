#!/usr/bin/env python3
"""An annotation that carries the relatedness of its genes (human CHM13).

    related_gtf.py ANNOT.slim.gff3.gz NCBI_gene_group.gz OUT_PREFIX     (environment as famtruth)

No public GTF holds paralogy or pseudogene parents; this joins the tables onto one. Writes
OUT_PREFIX.relations.tsv  (gene_a, name_a, gene_b, name_b, relation, source), one row per related pair:
    paralog_recent   Ensembl Compara, last common ancestor Homo sapiens .. Catarrhini
    paralog_older    Ensembl Compara, any deeper level
    pseudogene_of    gene_a is a pseudogene of gene_b: NCBI Gene `Related functional gene` (GeneIDs through HGNC to
                     symbols to CHM13 genes), else the naming rule SYMBOL + P<number>
    same_source      two CHM13 genes lifted from the same GRCh38 gene (CAT/Liftoff `source_gene`): extra copies
    readthrough      NCBI Gene readthrough parent / child / sibling
OUT_PREFIX.gtf  gene, transcript and exon lines of the annotation with gene attributes
    related_group   connected component of paralog_recent + pseudogene_of + same_source (`.` when alone)
    paralogs_recent, pseudogene_parent, gene_group (HGNC), copy_of (source gene when it has several copies)
Old paralogues are in the table and NOT in related_group: their reads place uniquely, they are separate loci."""
import sys, os, re, gzip, csv, collections
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import famtruth as F
gff, ncbi, out = sys.argv[1:]
G = F.load_cat()
byname = collections.defaultdict(set)
for g, d in G.items(): byname[d['name']].add(g)
lv = {x.rstrip('\n').split('\t')[3] for x in F._open(F.PARALOGS) if x.count('\t') >= 3 and x.rstrip('\n').split('\t')[3] and not x.startswith('Gene')}
anyp, _ = F.truth_pairs(G, F.PARALOGS, lv); young, _ = F.truth_pairs(G, F.PARALOGS, F.LEVELS['recent'])
rel = {}                                              # (a, b, relation) -> source
for p in anyp: rel[(p[0], p[1], 'paralog_recent' if p in young else 'paralog_older')] = 'Ensembl Compara'
e2s, group = {}, {}
for r in csv.DictReader(open(os.environ['DENOVO_HGNC']), delimiter='\t'):
    if r.get('entrez_id'): e2s[r['entrez_id']] = r['symbol']
    if r.get('gene_group'): group[r['symbol']] = r['gene_group'].replace('"', '')
npg = 0
for l in gzip.open(ncbi, 'rt'):
    f = l.rstrip('\n').split('\t')
    if f[0] != '9606' or f[3] != '9606': continue
    a, b = e2s.get(f[1]), e2s.get(f[4])
    if not a or not b: continue
    if f[2] == 'Related functional gene':
        for x in byname.get(a, ()):
            for y in byname.get(b, ()):
                if x != y: rel[(x, y, 'pseudogene_of')] = 'NCBI Gene'; npg += 1
    elif f[2].startswith('Readthrough') or f[2] == 'Potential readthrough sibling':
        for x in byname.get(a, ()):
            for y in byname.get(b, ()):
                if x != y: rel[(min(x, y), max(x, y), 'readthrough')] = 'NCBI Gene'
have = {a for (a, b, r) in rel if r == 'pseudogene_of'}
for g, d in G.items():
    if 'pseudogene' in d['biotype'] and g not in have:
        m = re.match(r'^(.+?)P\d+$', d['name'])
        if m:
            for y in byname.get(m.group(1), ()):
                if G[y]['biotype'] == 'protein_coding': rel[(g, y, 'pseudogene_of')] = 'name'
src = collections.defaultdict(set)
for l in open(F.CAT_GFF):
    f = l.rstrip('\n').split('\t')
    if len(f) > 8 and f[2] == 'gene':
        a = dict(x.split('=', 1) for x in f[8].split(';') if '=' in x)
        s = a.get('source_gene', 'N/A')
        if s != 'N/A' and a.get('gene_id') in G: src[s].add(a['gene_id'])
copy_of = {}
for s, gs in src.items():
    if len(gs) > 1:
        gs = sorted(gs)
        for g in gs: copy_of[g] = s
        for i, x in enumerate(gs):
            for y in gs[i + 1:]: rel[(x, y, 'same_source')] = 'CAT/Liftoff'
parent = {}
def find(x):
    while parent.setdefault(x, x) != x: parent[x] = parent[parent[x]]; x = parent[x]
    return x
for (a, b, r) in rel:
    if r in ('paralog_recent', 'pseudogene_of', 'same_source'): parent[find(a)] = find(b)
comp = collections.defaultdict(list)
for g in list(parent): comp[find(g)].append(g)
gid = {}
for i, (_, gs) in enumerate(sorted(comp.items(), key=lambda kv: (-len(kv[1]), min(kv[1])))):
    for g in gs: gid[g] = f'RG{i + 1:05d}'
with open(out + '.relations.tsv', 'w') as w:
    w.write('gene_a\tname_a\tgene_b\tname_b\trelation\tsource\n')
    for (a, b, r), s in sorted(rel.items(), key=lambda kv: (kv[0][2], kv[0][0], kv[0][1])): w.write(f"{a}\t{G[a]['name']}\t{b}\t{G[b]['name']}\t{r}\t{s}\n")
rec = collections.defaultdict(set); par = collections.defaultdict(set)
for (a, b, r) in rel:
    if r == 'paralog_recent': rec[a].add(G[b]['name']); rec[b].add(G[a]['name'])
    elif r == 'pseudogene_of': par[a].add(G[b]['name'])
def attrs(g):
    d = G[g]; a = f'gene_id "{g}"; gene_name "{d["name"]}"; gene_biotype "{d["biotype"]}"; related_group "{gid.get(g, ".")}";'
    if rec[g]: a += f' paralogs_recent "{",".join(sorted(rec[g])[:50])}";'
    if par[g]: a += f' pseudogene_parent "{",".join(sorted(par[g]))}";'
    if d['name'] in group: a += f' gene_group "{group[d["name"]]}";'
    if g in copy_of: a += f' copy_of "{copy_of[g]}";'
    return a
t2g = {}
with open(out + '.gtf', 'w') as w:
    for l in gzip.open(gff, 'rt'):
        f = l.rstrip('\n').split('\t')
        if len(f) < 9: continue
        a = dict(x.split('=', 1) for x in f[8].split(';') if '=' in x)
        if f[2] == 'gene' and a.get('gene_id') in G: w.write('\t'.join(f[:8]) + '\t' + attrs(a['gene_id']) + '\n')
        elif f[2] == 'transcript' and a.get('Parent') in G:
            t2g[a['ID']] = a['Parent']; w.write('\t'.join(f[:8]) + f'\tgene_id "{a["Parent"]}"; transcript_id "{a["ID"]}"; related_group "{gid.get(a["Parent"], ".")}";\n')
        elif f[2] == 'exon' and a.get('Parent') in t2g: w.write('\t'.join(f[:8]) + f'\tgene_id "{t2g[a["Parent"]]}"; transcript_id "{a["Parent"]}";\n')
c = collections.Counter((r, s) for (a, b, r), s in rel.items())
print('relations:', ', '.join(f'{r} ({s}) {n}' for (r, s), n in sorted(c.items())))
sz = collections.Counter(len(v) for v in comp.values())
print(f'genes {len(G)}, in a related_group {len(gid)}, groups {len(comp)} (largest {max(sz)}; of 2 genes {sz[2]}, of >= 10 genes {sum(n for k, n in sz.items() if k >= 10)})')
