#!/usr/bin/env python3
"""Find subfamily boundaries as GAPS in a family's within-family identity distribution.

Unbiased in the sense that matters here: no threshold is chosen by the operator and no external
truth is consulted. The largest gap in the sorted pairwise identities is located, then tested
against a null of "no gap" — identities resampled from a smooth unimodal fit with the same mean,
variance and n. A family with a real subfamily boundary has a gap the null does not produce; a
family homogenised by gene conversion does not.

usage: identity_gap.py <pairs.tsv>   with columns  a  b  identity
"""
import sys, math, random, statistics as st

rows = [l.split('\t') for l in open(sys.argv[1]) if l.strip() and not l.startswith('#')]
pairs = [(r[0], r[1], float(r[2])) for r in rows]
v = sorted(p[2] for p in pairs)
n = len(v)
if n < 6:
    sys.exit(f'need >= 6 pairs, got {n}')

def largest_gap(xs):
    """(gap, position) of the widest jump between consecutive sorted values, ignoring the
    outermost 10% on each side so a single outlier cannot masquerade as a boundary."""
    lo, hi = max(1, int(0.1 * len(xs))), min(len(xs) - 1, int(0.9 * len(xs)))
    best, at = 0.0, None
    for i in range(lo, hi):
        d = xs[i] - xs[i - 1]
        if d > best:
            best, at = d, (xs[i - 1] + xs[i]) / 2
    return best, at

obs, cut = largest_gap(v)
mu, sd = st.mean(v), st.pstdev(v)

# Null: same n, same mean/sd, drawn from a single smooth mode (no boundary).
random.seed(0)
null = []
for _ in range(10000):
    s = sorted(min(1.0, max(0.0, random.gauss(mu, sd))) for _ in range(n))
    null.append(largest_gap(s)[0])
null.sort()
p = sum(1 for x in null if x >= obs) / len(null)

print(f'pairs {n}   identity  min {v[0]:.4f}  median {st.median(v):.4f}  max {v[-1]:.4f}')
print(f'largest interior gap  {obs:.4f}  at identity {cut:.4f}')
print(f'null (unimodal, same mean/sd)  median gap {st.median(null):.4f}  95th {null[int(.95*len(null))]:.4f}')
print(f'p = {p:.4f}   -> {"SPLIT: a boundary the null does not produce" if p < 0.05 else "NO SPLIT: gap is what a single mode gives"}')

if p < 0.05:
    # The partition is the CONNECTED COMPONENTS of the subgraph above the gap — not a per-member
    # vote (a first attempt used majority-of-own-edges and put every member on one side, because
    # most pairs sit below the cut whatever the structure).
    members = sorted({x for a, b, _ in pairs for x in (a, b)})
    idx = {m: i for i, m in enumerate(members)}
    parent = list(range(len(members)))
    def find(x):
        while parent[x] != x:
            parent[x] = parent[parent[x]]; x = parent[x]
        return x
    def union(x, y):
        rx, ry = find(x), find(y)
        if rx != ry: parent[rx] = ry
    for a, b, i in pairs:
        if i >= cut: union(idx[a], idx[b])
    comp = {}
    for m in members: comp.setdefault(find(idx[m]), []).append(m)
    groups = sorted(comp.values(), key=len, reverse=True)
    print(f'  components of the subgraph above {cut:.4f}: {len(groups)}')
    for k, mem in enumerate(groups, 1):
        if len(mem) == 1: continue
        print(f'  G{k} ({len(mem)}): {sorted(mem)}')
    single = [m for g in groups if len(g) == 1 for m in g]
    if single: print(f'  singletons ({len(single)}): {sorted(single)}')
