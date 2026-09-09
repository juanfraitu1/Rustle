#!/usr/bin/env python3
"""Find subfamily boundaries as GAPS in a family's within-family identity distribution.

Unbiased in the sense that matters here: no threshold is chosen by the operator and no external
truth is consulted. The largest gap in the sorted pairwise identities is located, then tested
against a null of "no gap" — identities resampled from a smooth unimodal fit with the same mean,
variance and n. A family with a real subfamily boundary has a gap the null does not produce; a
family homogenised by gene conversion does not.

usage: identity_gap.py <pairs.tsv>   with columns  a  b  identity
"""
import sys, math, random, statistics as st, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import modality as M

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

# ---------------------------------------------------------------- nulls
# Three nulls, because the verdict must not rest on a distribution choice.
#  gauss  — the first form. ⚠ POOR for identities pressed against 1.0 (mass beyond the boundary).
#  beta   — moment-matched to the data; the natural family on [0,1].
#  smooth — smoothed bootstrap: resample the observed values and add Gaussian noise at Silverman's
#           bandwidth. Non-parametric, keeps the observed shape, and is the standard way to ask
#           "is this more clustered than one smooth mode?" without naming a family.
def nulls(n, v, mu, sd, reps=10000, seed=0):
    rnd = random.Random(seed)
    out = {}
    out['gauss'] = [largest_gap(sorted(min(1.0, max(0.0, rnd.gauss(mu, sd))) for _ in range(n)))[0]
                    for _ in range(reps)]
    res = {}
    if 0 < sd and 0 < mu < 1:
        t = mu * (1 - mu) / (sd * sd) - 1
        if t > 0:
            al, be = mu * t, (1 - mu) * t
            res['beta'] = [largest_gap(sorted(rnd.betavariate(al, be) for _ in range(n)))[0]
                           for _ in range(reps)]
    out.update(res)
    h = 0.9 * min(sd, (st.quantiles(v, n=4)[2] - st.quantiles(v, n=4)[0]) / 1.34) * n ** -0.2 if sd > 0 else 0.0
    if h > 0:
        out['smooth'] = [largest_gap(sorted(min(1.0, max(0.0, rnd.choice(v) + rnd.gauss(0, h)))
                                            for _ in range(n)))[0] for _ in range(reps)]
    return out

NUL = nulls(n, v, mu, sd)
pv = {k: sum(1 for x in d if x >= obs) / len(d) for k, d in NUL.items()}
p = max(pv.values())          # the CONSERVATIVE verdict: the least favourable null governs

print(f'pairs {n}   identity  min {v[0]:.4f}  median {st.median(v):.4f}  max {v[-1]:.4f}')
print(f'largest interior gap  {obs:.4f}  at identity {cut:.4f}')
for k in sorted(NUL):
    d = sorted(NUL[k])
    print(f'  null {k:7s} median {st.median(d):.4f}  95th {d[int(.95*len(d))]:.4f}  p = {pv[k]:.4f}')
print(f'p = {p:.4f} (worst null governs)   -> '
      f'{"SPLIT: a boundary no smooth mode produces" if p < 0.05 else "NO SPLIT: gap is what a single mode gives"}')

# ---------------------------------------------------------------- whole-distribution tests (PREREG_modality)
sp, hc = M.silverman_p(v, B=1000)
dbic = M.mixture_dbic(v)
print(f"silverman  p = {sp:.4f}  (h_crit {hc:.4f})   -> {'MULTIMODAL' if sp < 0.05 else 'unimodal'}")
print(f"mixture    dBIC = {dbic:+.1f}                 -> {'TWO components (very strong)' if dbic >= 10 else 'one component'}")

pp, contrast = M.partition_perm_p(pairs, B=1000)
if contrast != contrast or contrast < 0:
    print("partition  no threshold yields two components of >= 3 members -> NO SUBFAMILY PARTITION EXISTS")
else:
    print(f"partition  contrast {contrast:.4f}  permutation p = {pp:.4f}   -> "
          f"{'PARTITION CERTIFIED (identities are not exchangeable across edges)' if pp < 0.05 else 'not certified'}")

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
