#!/usr/bin/env python3
"""Report a gene family as a HIERARCHY rather than as one thresholded cut.

Usage: family_hierarchy.py <all-vs-all.paf> [label_prefix_len]

WHY. A family and a subfamily are not different objects, they are different levels of one structure, and
choosing a single identity or coverage threshold asserts a level rather than measuring it. NPIP is the case
that forces the issue: at DNA level it is a single clique -- 168 of 171 annotated gene pairs clear the
shipped edge rule at median identity 0.978 -- yet the expected NPIPA/NPIPB division is invisible to every
signal tested (identity uniform at ~0.98, genomic-span coverage 0.86-1.00, annotated splice structure 0.71
between groups against 0.80-0.83 within). Producing four groups by tuning a threshold would be fitting to a
number, not finding structure.

So build the hierarchy once and let the levels speak. Two views:

  MINIMUM SPANNING TREE   the backbone: which copy is each copy's nearest relative, and how far. Cutting the
                          k-1 longest MST edges gives exactly k groups, so the MST edge lengths ARE the
                          thresholds -- a large gap between consecutive edges is a level with support, a
                          smooth run of edges means no natural cut exists at all.
  AVERAGE-LINKAGE TREE    the nesting, printed as an indented dendrogram with the merge distance shown, so
                          a reader can see which level a proposed subfamily corresponds to.

Distance is 1 - identity by default. Coverage is reported alongside but not mixed in: on genomic spans it
saturates (0.86-1.00) and would only add noise, while on exon-sums it is dominated by isoform choice rather
than by divergence.
"""
import sys
from collections import defaultdict

import numpy as np
from scipy.cluster.hierarchy import linkage, dendrogram
from scipy.sparse.csgraph import minimum_spanning_tree
from scipy.spatial.distance import squareform

PAF = sys.argv[1]
PFX = int(sys.argv[2]) if len(sys.argv) > 2 else 4

best = {}
names = set()
for ln in open(PAF):
    f = ln.rstrip("\n").split("\t")
    q, t = f[0], f[5]
    names.add(q)
    names.add(t)
    if q == t:
        continue
    ident = 1.0 - float(next((x[5:] for x in f[12:] if x.startswith("de:f:")), "1"))
    cov = (int(f[3]) - int(f[2])) / min(int(f[1]), int(f[6]))
    k = tuple(sorted((q, t)))
    if k not in best or ident > best[k][0]:
        best[k] = (ident, cov)
g = sorted(names)
n = len(g)
idx = {x: i for i, x in enumerate(g)}
# missing pairs are set to the worst observed distance rather than infinity, so the tree stays finite and
# an absent alignment is treated as "no closer than anything we did see" instead of silently dropping out
D = np.zeros((n, n))
seen = [1.0 - v[0] for v in best.values()]
worst = max(seen) * 1.5 if seen else 1.0
for i in range(n):
    for j in range(n):
        if i != j:
            v = best.get(tuple(sorted((g[i], g[j]))))
            D[i, j] = (1.0 - v[0]) if v else worst
print(f"{n} sequences, {len(best)} aligned pairs of {n*(n-1)//2}")
print(f"identity: min {1-max(seen):.4f}  median {1-sorted(seen)[len(seen)//2]:.4f}  max {1-min(seen):.4f}\n")

mst = minimum_spanning_tree(D).tocoo()
edges = sorted(zip(mst.data, mst.row, mst.col))
print("MINIMUM SPANNING TREE, edges shortest first (cutting the k-1 longest gives k groups)")
print(f"{'rank':>5}{'distance':>11}{'identity':>10}  pair")
prev = None
for r, (d, i, j) in enumerate(edges, 1):
    gap = "" if prev is None else f"   <-- gap x{d/prev:.1f}" if prev > 0 and d / prev >= 2.0 else ""
    print(f"{r:>5}{d:>11.5f}{1-d:>10.4f}  {g[i]} -- {g[j]}{gap}")
    prev = d

Z = linkage(squareform(D, checks=False), method="average")
print("\nAVERAGE-LINKAGE HIERARCHY (indent = nesting; distance = merge height)")


def show(node, depth):
    n_leaves = len(g)
    if node < n_leaves:
        print(f"{'  '*depth}{g[node]}")
        return
    row = Z[int(node) - n_leaves]
    print(f"{'  '*depth}+-- merge at d={row[2]:.5f} (identity {1-row[2]:.4f}), {int(row[3])} members")
    show(int(row[0]), depth + 1)
    show(int(row[1]), depth + 1)


show(len(Z) + len(g) - 1, 0)

# does any level correspond to the naming convention?
print("\ndoes any cut reproduce the naming groups?")
lab = {x: x[PFX] if len(x) > PFX else "?" for x in g}
groups = sorted(set(lab.values()))
from scipy.cluster.hierarchy import fcluster
for k in range(2, min(9, n)):
    cl = fcluster(Z, k, criterion="maxclust")
    part = defaultdict(set)
    for x, c in zip(g, cl):
        part[c].add(lab[x])
    pure = sum(1 for v in part.values() if len(v) == 1)
    print(f"  k={k}: {pure}/{k} clusters are label-pure  "
          f"{[''.join(sorted(v)) for v in part.values()]}")
