#!/usr/bin/env python3
"""Define gene families from genome self-alignment alone: NO annotation, NO reads, NO gene models.

Usage: dna_only_families.py <self_alignment.paf> <region_offset> [min_bp]

The earlier attempt failed because it MERGED overlapping duplication intervals into blocks. Any-overlap
merging chains transitively, so a 4 kb unit and a 50 kb unit sharing one base become one node, and the
result was 86 mega-blocks at 6.3x gene size whose best NPIP component was 12% NPIP and 88% other genes.
That was an artifact of the aggregation, not a limit of self-alignment: individual alignment records in the
same region are gene-scale (median 4,271 bp where they overlap an NPIP gene).

So this keeps the records and fixes the node definition instead. Two intervals become the same node only if
they RECIPROCALLY overlap -- each covering at least `RECIP` of the other -- which prevents a short unit from
being swallowed by a long one and stops the transitive chaining that produced the blocks.

Nodes are then clustered by the identity of the alignments between them, using the same average-linkage
hierarchy as `family_hierarchy.py`, so "family" and "subfamily" stay levels of one tree rather than
separate thresholded objects.

Genes are read ONLY at the end, to score what the blind procedure produced.
"""
import gzip
import re
import sys
from collections import defaultdict

import numpy as np
from scipy.cluster.hierarchy import linkage, fcluster
from scipy.spatial.distance import squareform

PAF = sys.argv[1]
OFF = int(sys.argv[2])
MIN_BP = int(sys.argv[3]) if len(sys.argv) > 3 else 1000
RECIP = 0.50

recs = []
for ln in open(PAF):
    f = ln.split("\t")
    qs, qe, ts, te = int(f[2]), int(f[3]), int(f[7]), int(f[8])
    if qe - qs < MIN_BP or te - ts < MIN_BP:
        continue
    de = next((float(x[5:]) for x in f[12:] if x.startswith("de:f:")), None)
    if de is None:
        continue
    recs.append((qs, qe, ts, te, 1.0 - de))
print(f"{len(recs)} self-alignment records >= {MIN_BP} bp")

# ---- nodes: reciprocal-overlap clusters of intervals (NOT transitive any-overlap merging) --------
iv = sorted({(a, b) for (qs, qe, ts, te, _) in recs for (a, b) in ((qs, qe), (ts, te))})
parent = list(range(len(iv)))


def find(x):
    while parent[x] != x:
        parent[x] = parent[parent[x]]
        x = parent[x]
    return x


for i in range(len(iv)):
    a, b = iv[i]
    for j in range(i + 1, len(iv)):
        c, d = iv[j]
        if c >= b:
            break
        o = min(b, d) - max(a, c)
        if o >= RECIP * (b - a) and o >= RECIP * (d - c):
            ri, rj = find(i), find(j)
            if ri != rj:
                parent[ri] = rj
grp = defaultdict(list)
for i in range(len(iv)):
    grp[find(i)].append(i)
nodes = []
for members in grp.values():
    lo = min(iv[i][0] for i in members)
    hi = max(iv[i][1] for i in members)
    nodes.append((lo, hi))
nodes.sort()
print(f"{len(iv)} distinct intervals -> {len(nodes)} nodes "
      f"(median {int(np.median([b-a for a,b in nodes])):,} bp, max {max(b-a for a,b in nodes):,} bp)")


def node_of(a, b):
    best, bo = None, 0
    for k, (x, y) in enumerate(nodes):
        o = min(b, y) - max(a, x)
        if o > bo:
            best, bo = k, o
    return best


w = defaultdict(list)
for (qs, qe, ts, te, ident) in recs:
    u, v = node_of(qs, qe), node_of(ts, te)
    if u is not None and v is not None and u != v:
        w[(min(u, v), max(u, v))].append(ident)
print(f"{len(w)} node pairs carry an alignment\n")

n = len(nodes)
seen = [max(v) for v in w.values()]
worst = (1.0 - min(seen)) * 1.5
D = np.full((n, n), worst)
np.fill_diagonal(D, 0.0)
for (u, v), vals in w.items():
    D[u, v] = D[v, u] = 1.0 - max(vals)
Z = linkage(squareform(D, checks=False), method="average")

# ---- only now, genes ------------------------------------------------------------------------------
genes = []
with gzip.open("/mnt/c/Users/jfris/Desktop/Reference/chm13v2.0_RefSeq_full.gff.gz", "rt") as fh:
    for ln in fh:
        if ln[0] == "#":
            continue
        f = ln.rstrip("\n").split("\t")
        if len(f) < 9 or f[2] != "gene" or f[0] != "chr16":
            continue
        m = re.search(r"Name=([^;]+)", f[8])
        s, e = int(f[3]) - 1 - OFF, int(f[4]) - OFF
        if m and 0 <= s:
            genes.append((s, e, m.group(1)))


def genes_at(a, b):
    return {g for (x, y, g) in genes if a < y and b > x}


npip = {g for _, _, g in genes if g.upper().startswith("NPIP")}
print(f"{'k':>4}{'best NPIP cluster':>20}{'NPIP in it':>12}{'other genes':>13}{'purity':>8}")
for k in (2, 4, 6, 8, 12, 16, 24, 32):
    if k >= n:
        break
    cl = fcluster(Z, k, criterion="maxclust")
    best = None
    for c in set(cl):
        gs = set()
        for i in range(n):
            if cl[i] == c:
                gs |= genes_at(*nodes[i])
        hit = len(gs & npip)
        if best is None or hit > best[0]:
            best = (hit, c, gs)
    hit, c, gs = best
    print(f"{k:>4}{('cluster '+str(c)):>20}{f'{hit}/{len(npip)}':>12}"
          f"{len(gs)-hit:>13}{hit/max(len(gs),1):>8.2f}")
