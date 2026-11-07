#!/usr/bin/env python3
"""Why O1 misses or invents Compara pairs (diagnostic, no bar; same truth and mapping as pairscore.py).

usage: pairdiag.py [--species human|hg002|...] RUN_PREFIX OUT_TSV   (env as pairscore: DENOVO_CAT_GFF/TSV, DENOVO_PARALOGS, DENOVO_HGNC)
Missed truth pairs (x, y):
  absorbed   a gene is expressed but no locus representative maps to it (its reads sit in a neighbour's / fused locus)
  unlinked   a gene has its own locus but that locus is in no family
  split      both genes are in families, different ones
  each with: any gated-PAF record between their loci (aligned: edge failed the gate or MCL cut it) or none (never aligned)
False pairs: older_paralog (Compara any level), hgnc_group, unrelated; unrelated ones also tagged fused when a copy of
either gene overlaps exons of >= 2 protein-coding genes.
"""
import collections, csv, itertools, os, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import famtruth as F

a = sys.argv[1:]
species = a[a.index('--species') + 1] if '--species' in a else 'human'
prefix, out = [x for i, x in enumerate(a) if not x.startswith('--') and (i == 0 or a[i - 1] not in F.SPECIES_FLAGS)][:2]
G = F.genes_from_argv(a, species)
idx = F.exon_index(G)
expressed = F.expressed_from_gtf(prefix + '.families.gtf', idx)
truth, universe = F.truth_pairs(G, F.PARALOGS, F.levels('simiiformes'))
truth = {tuple(sorted(p)) for p in truth if p[0] in expressed and p[1] in expressed}
fams, _, _ = F.load_catalog(prefix + '.fam.copies.tsv', G, idx)
fam_of = collections.defaultdict(set)
for k, v in fams.items():
    for g in v:
        fam_of[g].add(k)
pc = lambda g: G[g]['biotype'] == 'protein_coding' and g in universe
pred = {(x, y) for v in fams.values() for x, y in itertools.combinations(sorted(v), 2)
        if pc(x) and pc(y) and G[x]['name'] != G[y]['name']}

# loci: representative exons from the loci GFF -> gene; PAF key = chrom:start-end (1-based, as the GFF)
loci_of, genes_hit_by_locus = collections.defaultdict(set), {}
cur, ex = None, collections.defaultdict(list)
for line in open(prefix + '.fam.loci.gff3'):
    f = line.split('\t')
    if len(f) < 9:
        continue
    if f[2] == 'gene':
        cur = (f[0], f[6], f"{f[0]}:{f[3]}-{f[4]}")
    elif f[2] == 'exon':
        ex[cur].append((int(f[3]) - 1, int(f[4])))
for (chrom, strand, key), e in ex.items():
    e.sort()
    g = F.map_copy(idx, chrom, strand if strand in '+-' else '+', e)
    if g:
        loci_of[g].add(key)
    # every protein-coding gene whose exons this locus touches (fusion tag)
    hit = set()
    v = idx.get((chrom, strand if strand in '+-' else '+'))
    if v:
        import bisect
        lst, st, m = v
        for s, t in e:
            i = bisect.bisect_left(st, s - m)
            while i < len(lst) and lst[i][0] < t:
                if lst[i][1] > s and G[lst[i][2]]['biotype'] == 'protein_coding':
                    hit.add(lst[i][2])
                i += 1
    genes_hit_by_locus[key] = hit

missed = truth - pred
need = {k for p in missed for g in p for k in loci_of.get(g, ())}
ident = {}
for line in open(prefix + '.fam.loci.gated.paf'):
    f = line.split('\t', 12)
    if f[0] in need and f[5] in need and f[0] != f[5]:
        k = frozenset((f[0], f[5]))
        ident[k] = max(ident.get(k, 0.0), int(f[9]) / int(f[10]))

rows, cnt = [], collections.Counter()
for x, y in sorted(missed):
    if not loci_of.get(x) or not loci_of.get(y):
        cls = 'absorbed'
    elif not fam_of.get(x) or not fam_of.get(y):
        cls = 'unlinked'
    else:
        cls = 'split'
    ids = [ident[frozenset((a, b))] for a in loci_of.get(x, ()) for b in loci_of.get(y, ()) if frozenset((a, b)) in ident]
    aligned = f"aligned_{'hi' if max(ids) >= 0.90 else 'lo'}" if ids else 'never_aligned'
    cnt[(cls, aligned if cls != 'absorbed' else '-')] += 1
    rows.append(('missed', cls, aligned, x, y, G[x]['name'], G[y]['name'], f"{max(ids):.3f}" if ids else 'NA'))

anyp = set()
for r in csv.reader(F._open(F.PARALOGS), delimiter='\t'):
    if len(r) >= 2 and r[1] and not r[0].startswith('Gene'):
        a, b = r[0].split('.')[0], r[1].split('.')[0]
        anyp.add((min(a, b), max(a, b)))
groups = collections.defaultdict(set)
for r in csv.DictReader(open(os.environ['DENOVO_HGNC']), delimiter='\t'):
    for gg in r['gene_group_id'].split('|'):
        if gg and r['ensembl_gene_id']:
            groups[r['ensembl_gene_id']].add(gg)
ens = lambda g: {e.split('.')[0] for e in G[g]['ens']}
for x, y in sorted(pred - truth):
    if any((min(p, q), max(p, q)) in anyp for p in ens(x) for q in ens(y)):
        cls = 'older_paralog'
    elif any(groups[p] & groups[q] for p in ens(x) for q in ens(y)):
        cls = 'hgnc_group'
    else:
        fused = any(len(genes_hit_by_locus.get(k, ())) >= 2 for g in (x, y) for k in loci_of.get(g, ()))
        cls = 'unrelated_fused' if fused else 'unrelated'
    cnt[('false', cls)] += 1
    rows.append(('false', cls, '', x, y, G[x]['name'], G[y]['name'], ''))

with open(out, 'w') as fh:
    fh.write('kind\tclass\taligned\tgene_x\tgene_y\tname_x\tname_y\tmax_paf_identity\n')
    for r in rows:
        fh.write('\t'.join(r) + '\n')
print(f"truth {len(truth)} pred {len(pred)} hit {len(truth & pred)} missed {len(missed)} false {len(pred - truth)}")
for k in sorted(cnt):
    print(*k, cnt[k], sep='\t')
