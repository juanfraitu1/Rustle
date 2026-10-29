#!/usr/bin/env python3
"""HYBRID families (2026-10-30, mainstay diagnosis: the minimap2 base links repeat-rich pieces and fragments that exact k-mers cannot,
the sketch cut is more specific): graph = sketch edges (OUT.edges.tsv of sketchfam) UNION base alignment edges from the families
stage PAF (loci named chrom:start-end), weight = max(sketch containment, alignment matches / shorter locus length); families = the
same per-component UPGMA tree cut at 1 - root_cont as sketchfam --cluster tree.
  hybrid.py SKETCH_EDGES_TSV LOCI_PAF LOCI_GFF3 OUT [--root-cont 0.15] [--min-identity 0.8] [--min-cov 0.5] [--same-chrom]
                [--no-sketch]   (--no-sketch = alignment edges only, same cut)"""
import sys, re, collections, numpy as np
from scipy.cluster.hierarchy import linkage, fcluster
from scipy.spatial.distance import squareform
a = sys.argv[1:]; edges_f, paf, gff, out = a[:4]
opt = lambda k, d, f=float: f(a[a.index(k) + 1]) if k in a else d
ROOTC = opt('--root-cont', 0.15); MINID = opt('--min-identity', 0.8); MINCOV = opt('--min-cov', 0.5); SAME = '--same-chrom' in a; NOSK = '--no-sketch' in a
loci = {}; exons = collections.defaultdict(list); span2id = {}
for l in open(gff):
    if l[0] == '#': continue
    f = l.rstrip('\n').split('\t')
    if f[2] == 'gene':
        t = re.search(r'ID=gene-([^;]+)', f[8]).group(1); loci[t] = (f[0], f[6] if f[6] in '+-' else '+'); span2id[f"{f[0]}:{f[3]}-{f[4]}"] = t
    elif f[2] == 'exon':
        exons[re.search(r'gene=([^;]+)', f[8]).group(1)].append((int(f[3]), int(f[4])))
exlen = {t: sum(b - s + 1 for s, b in v) for t, v in exons.items()}
w = collections.defaultdict(float); n_sk = n_al = 0
if not NOSK:
    for l in open(edges_f):
        x, y, n, c, m = l.split('\t')[:5]; k = (min(x, y), max(x, y)); w[k] = max(w[k], float(c)); n_sk += 1
# alignment edges: per pair, best single alignment by matches / min exonic length, identity gate
best = {}
for l in open(paf):
    f = l.split('\t', 12)
    if len(f) < 12: continue
    q, t = span2id.get(f[0]), span2id.get(f[5])
    if not q or not t or q == t: continue
    if SAME and loci[q][0] != loci[t][0]: continue
    m, bl = int(f[9]), int(f[10])
    if bl == 0 or m / bl < MINID: continue
    k = (min(q, t), max(q, t)); c = m / max(1, min(exlen.get(q, 1), exlen.get(t, 1)))
    if c > best.get(k, 0): best[k] = c
for k, c in best.items():
    if c >= MINCOV:
        if c > w[k]: w[k] = min(1.0, c)
        n_al += 1
parent = {}
def find(x):
    parent.setdefault(x, x)
    while parent[x] != x: parent[x] = parent[parent[x]]; x = parent[x]
    return x
for x, y in w: parent[find(x)] = find(y)
comp = collections.defaultdict(list); cedges = collections.defaultdict(list)
for x, y in w: cedges[find(x)].append((x, y, w[(x, y)]))
for x in parent: comp[find(x)].append(x)
fams = []
for root, m in comp.items():
    if len(m) < 2: continue
    m = sorted(m); n = len(m); D = np.ones((n, n)); np.fill_diagonal(D, 0.0); ix = {t: i for i, t in enumerate(m)}
    for x, y, c in cedges[root]: D[ix[x], ix[y]] = D[ix[y], ix[x]] = 1.0 - c
    if n == 2: lab = [1, 1] if D[0, 1] <= 1.0 - ROOTC else [1, 2]
    else: lab = fcluster(linkage(squareform(D, checks=False), method='average'), t=1.0 - ROOTC, criterion='distance')
    g = collections.defaultdict(list)
    for t, l in zip(m, lab): g[l].append(t)
    fams += [v for v in g.values() if len(v) >= 2]
fams.sort(key=lambda v: -len(v))
with open(out + '.fam.copies.tsv', 'w') as o:
    o.write("family_id\tcopy_idx\ttid\tchrom\tstart\tend\tn_exon\tstrand\tn_reads\texons\n")
    for i, v in enumerate(fams):
        for j, t in enumerate(sorted(v)):
            c, s = loci[t]; ex = sorted(exons[t])
            o.write(f"HY{i}\t{j}\t{t}\t{c}\t{ex[0][0]}\t{ex[-1][1]}\t{len(ex)}\t{s}\t0\t{','.join(f'{p}-{q}' for p, q in ex)}\n")
print(f"sketch edges {n_sk}, alignment edges passing {n_al}, union {len(w)}, families {len(fams)}, members {sum(len(v) for v in fams)}, largest {[len(v) for v in fams[:5]]}", file=sys.stderr)
