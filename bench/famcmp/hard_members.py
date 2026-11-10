#!/usr/bin/env python3
"""PREREG r1601 (docs/PREREG_hard_members_2026-11-10.md): levels L1-L4 per hard member (missed by StringTie and FLAIR).

  hard_members.py GENES_TSV CAT_TSV CAT_SLIM_GFF3.gz OUT_TSV NAME=GTF ...
"""
import sys, re, gzip, collections

genes_tsv, cat_tsv, cat_gff, out = sys.argv[1:5]
arms = [a.split('=', 1) for a in sys.argv[5:]]
blocks = lambda s: [(int(a), int(b)) for a, b in (x.split('-') for x in s.split(','))]
ov = lambda a, b, x, y: a < y and x < b

with open(genes_tsv) as fh:
    hdr = fh.readline().rstrip('\n').split('\t')
    rows = [dict(zip(hdr, l.rstrip('\n').split('\t'))) for l in fh]
num = lambda r, k: float(r[k] or 0)
fam_of = {r['gene']: r['family'] for r in rows}
members = collections.defaultdict(list)
for r in rows: members[r['family']].append(r['gene'])
H = [r for r in rows if num(r, 'stringtie') == 0 and num(r, 'flair') == 0]
cls = lambda r: 'T-only' if num(r, 'primary_molecules') == 0 else ('mapq0' if num(r, 'mapq0_share') >= 0.5 else 'ordinary')

G = {}
for l in open(cat_tsv):
    f = l.rstrip('\n').split('\t')
    if f[0] in fam_of: G[f[0]] = dict(chrom=f[4], strand=f[7], exons=blocks(f[9]))
tx = collections.defaultdict(list); tx_gene = {}
for l in gzip.open(cat_gff, 'rt'):
    f = l.rstrip('\n').split('\t')
    if len(f) < 9: continue
    if f[2] == 'transcript':
        a = dict(x.split('=', 1) for x in f[8].split(';') if '=' in x)
        if a.get('gene_id') in fam_of: tx_gene[a['ID']] = a['gene_id']
    elif f[2] == 'exon':
        p = re.search(r'Parent=([^;]+)', f[8]).group(1)
        if p in tx_gene: tx[p].append((int(f[3]) - 1, int(f[4])))
introns = lambda ex: tuple((ex[i][1], ex[i + 1][0]) for i in range(len(ex) - 1))
chains = collections.defaultdict(set)
for t, ex in tx.items(): ex.sort(); chains[tx_gene[t]].add(introns(ex))
cat_introns = {g: set(i for c in cs for i in c) for g, cs in chains.items()}

def read_models(path, keys):
    ex = collections.defaultdict(list); st = {}
    op = gzip.open if path.endswith('.gz') else open
    for l in op(path, 'rt'):
        f = l.split('\t', 9)
        if len(f) < 9 or f[2] != 'exon' or (f[0], f[6]) not in keys: continue
        t = re.search(r'transcript_id "([^"]+)"', f[8]).group(1)
        ex[(f[0], t)].append((int(f[3]) - 1, int(f[4]))); st[(f[0], t)] = f[6]
    by = collections.defaultdict(list)
    for (c, t), e in ex.items(): e.sort(); by[(c, st[(c, t)])].append(e)
    return by

keys = {(v['chrom'], v['strand']) for v in G.values()}
res = {}
for name, path in arms:
    M = read_models(path, keys)
    for r in H:
        g = r['gene']; v = G[g]; gx = v['exons']; lo, hi = gx[0][0], gx[-1][1]
        others = [G[o]['exons'] for o in members[r['family']] if o != g and o in G and G[o]['chrom'] == v['chrom'] and G[o]['strand'] == v['strand']]
        lvl = 0
        for m in M.get((v['chrom'], v['strand']), []):
            if m[-1][1] <= lo or m[0][0] >= hi or not any(ov(a, b, x, y) for a, b in m for x, y in gx): continue
            lvl = max(lvl, 1)
            if any(ov(a, b, x, y) for oe in others for a, b in m for x, y in oe): continue
            lvl = max(lvl, 2)
            mi = introns(m)
            if any(chains.get(g, ())):  # spliced member: some CAT transcript has an intron
                if set(mi) & cat_introns.get(g, set()): lvl = max(lvl, 3)
                if mi and mi in chains.get(g, set()): lvl = 4
            else:
                inside = sum(max(0, min(b, y) - max(a, x)) for a, b in m for x, y in gx)
                if inside >= 0.5 * sum(b - a for a, b in m): lvl = max(lvl, 3)
        res[(name, g)] = lvl

with open(out, 'w') as w:
    w.write('family\tgene\tname\tclass\t' + '\t'.join(n for n, _ in arms) + '\n')
    for r in H: w.write(f"{r['family']}\t{r['gene']}\t{r['name']}\t{cls(r)}\t" + '\t'.join(str(res[(n, r['gene'])]) for n, _ in arms) + '\n')
print('arm\tn\tL1+\tL2+\tL3+\tL4')
for n, _ in arms:
    L = [res[(n, r['gene'])] for r in H]
    print(n, len(L), *(sum(x >= k for x in L) for k in (1, 2, 3)), sum(x == 4 for x in L), sep='\t')
for n, _ in arms:
    print(n, 'L3+ by class', dict(collections.Counter(cls(r) for r in H if res[(n, r['gene'])] >= 3)))
