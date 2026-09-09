#!/usr/bin/env python3
"""Examine families with a certified two-subfamily partition: components, chromosomes, identities, names.
usage: examine_partitions.py <clusters.tsv> <pairs.tsv> <partition.tsv> <gff> [--p 0.05]"""
import sys, csv, collections, re, bisect, statistics as st
import numpy as np
clusters_p, pairs_p, part_p, gff_p = sys.argv[1:5]
PT = float(sys.argv[sys.argv.index('--p') + 1]) if '--p' in sys.argv else 0.05
MIN_COMP = 3

# gene names: per-contig sorted intervals
genes = collections.defaultdict(list)
for line in open(gff_p):
    if line[0] == '#': continue
    f = line.rstrip('\n').split('\t')
    if len(f) < 9 or f[2] not in ('gene', 'pseudogene'): continue
    m = re.search(r'Name=([^;]+)', f[8])
    genes[f[0]].append((int(f[3]) - 1, int(f[4]), m.group(1) if m else '?'))
for c in genes: genes[c].sort()
gstart = {c: [g[0] for g in v] for c, v in genes.items()}
def name(c, s, e):
    v = genes.get(c)
    if not v: return '?'
    i = bisect.bisect_right(gstart[c], e); best, bo = '?', 0
    for gs, ge, n in v[max(0, i - 30):i]:
        o = min(e, ge) - max(s, gs)
        if o > bo: bo, best = o, n
    return best

# members
mem = collections.defaultdict(dict)       # fam -> "start-end" -> (chrom,start,end)
for r in csv.DictReader(open(clusters_p), delimiter='\t'):
    mem[r['cluster_id']][f"{r['start']}-{r['end']}"] = (r['chrom'], int(r['start']), int(r['end']))
pairs = collections.defaultdict(list)
for r in csv.DictReader(open(pairs_p), delimiter='\t'):
    pairs[r['family']].append((r['a'], r['b'], float(r['identity'])))
cert = [r for r in csv.DictReader(open(part_p), delimiter='\t') if r['partition_p'] != 'nan' and float(r['partition_p']) < PT]
cert.sort(key=lambda r: (float(r['partition_p']), -float(r['contrast'])))

def best_partition(pr):
    names = {}; E = []; I = []
    for a, b, i in pr:
        E.append((names.setdefault(a, len(names)), names.setdefault(b, len(names)))); I.append(i)
    E = np.array(E); I = np.array(I); n = len(names); inv_names = {v: k for k, v in names.items()}
    thr = np.unique(np.quantile(I, np.linspace(0.02, 0.98, 40)))[::-1]
    order = np.argsort(-I); parent = np.arange(n)
    def find(x):
        while parent[x] != x: parent[x] = parent[parent[x]]; x = parent[x]
        return x
    best, best_t, best_roots = -1.0, None, None; k = 0
    for t in thr:
        while k < I.size and I[order[k]] >= t:
            a, b = E[order[k]]; ra, rb = find(a), find(b)
            if ra != rb: parent[ra] = rb
            k += 1
        roots = np.array([find(i) for i in range(n)])
        _, inv, cnt = np.unique(roots, return_inverse=True, return_counts=True)
        big = cnt >= MIN_COMP
        if big.sum() < 2: continue
        c1, c2 = inv[E[:, 0]], inv[E[:, 1]]; keep = big[c1] & big[c2]
        if not keep.any(): continue
        same = c1[keep] == c2[keep]
        if same.sum() == 0 or (~same).sum() == 0: continue
        c = I[keep][same].mean() - I[keep][~same].mean()
        if c > best: best, best_t, best_roots = c, t, roots.copy()
    comps = collections.defaultdict(list)
    for i, r in enumerate(best_roots): comps[r].append(inv_names[i])
    comps = sorted([v for v in comps.values() if len(v) >= MIN_COMP], key=len, reverse=True)
    idx = {m: ci for ci, v in enumerate(comps) for m in v}
    within = [i for a, b, i in pr if a in idx and b in idx and idx[a] == idx[b]]
    between = [i for a, b, i in pr if a in idx and b in idx and idx[a] != idx[b]]
    return best_t, comps, st.median(within), st.median(between)

print(f"{len(cert)} families certified at p < {PT}\n")
for r in cert:
    fam = r['family']; t, comps, w, b = best_partition(pairs[fam])
    chroms = collections.Counter(mem[fam][m][0] for v in comps for m in v)
    flag = '⚠ OVER-MERGE?' if b < 0.85 else ('subfamilies' if b < 0.97 else 'shallow')
    print(f"== {fam}  members {len(mem[fam])}  p={r['partition_p']}  contrast {float(r['contrast']):.3f}  "
          f"cut {t:.4f}  within {w:.4f} / between {b:.4f}  [{flag}]  chroms {dict(chroms)}")
    for ci, v in enumerate(comps[:4]):
        nm = collections.Counter(name(*mem[fam][m]) for m in v)
        top = ', '.join(f"{k}×{n}" if n > 1 else k for k, n in nm.most_common(6))
        cc = collections.Counter(mem[fam][m][0] for m in v)
        print(f"     C{ci+1} n={len(v):3d}  {dict(cc)}  {top}")
