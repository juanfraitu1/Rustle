#!/usr/bin/env python3
"""Whole-distribution tests of subfamily structure in a family's pairwise identities.

Two statistics, both using every pair rather than the single largest gap:

  silverman_p  — Silverman's (1981) critical-bandwidth test of unimodality. Find the smallest
                 Gaussian-KDE bandwidth h_crit at which the sample has ONE mode (mode count is
                 monotone in h for a Gaussian kernel, so bisection is exact); bootstrap from that
                 KDE with the variance correction and report how often a unimodal sample needs a
                 bandwidth >= h_crit. Small p = the data need more smoothing than one mode explains.
  mixture_dbic — 1- vs 2-component Gaussian mixture by EM; dBIC = BIC1 - BIC2. Kass & Raftery:
                 >= 10 is "very strong" evidence for two components.

No free parameters are tuned on any family; the only knobs are grid size and bootstrap count.
"""
import math
import numpy as np

GRID = 512

def _kde_modes(x, h, grid):
    """Number of local maxima of the Gaussian KDE of x (bandwidth h) evaluated on grid."""
    z = (grid[:, None] - x[None, :]) / h
    f = np.exp(-0.5 * z * z).sum(axis=1)
    d = np.diff(f)
    # a maximum is a + to - sign change of the derivative
    return int(np.sum((d[:-1] > 0) & (d[1:] <= 0)))

def h_crit(x, k=1, tol=1e-3):
    """Smallest bandwidth at which the KDE has <= k modes (bisection; monotone in h)."""
    x = np.asarray(x, dtype=float)
    lo_g, hi_g = x.min(), x.max()
    rng = hi_g - lo_g
    if rng <= 0:
        return 0.0
    grid = np.linspace(lo_g, hi_g, GRID)
    hi = rng                      # one mode guaranteed here
    lo = rng * 1e-4
    if _kde_modes(x, lo, grid) <= k:
        return lo
    while (hi - lo) / hi > tol:
        mid = math.sqrt(lo * hi)  # geometric bisection: h spans orders of magnitude
        if _kde_modes(x, mid, grid) <= k:
            hi = mid
        else:
            lo = mid
    return hi

def silverman_p(x, B=200, seed=0):
    """P(a unimodal population needs bandwidth >= h_crit(x)). Returns (p, h_crit)."""
    x = np.asarray(x, dtype=float)
    n = x.size
    if n < 8 or x.std() == 0:
        return float('nan'), float('nan')
    hc = h_crit(x)
    if hc <= 0:
        return float('nan'), hc
    rnd = np.random.default_rng(seed)
    mu, var = x.mean(), x.var()
    shrink = 1.0 / math.sqrt(1.0 + hc * hc / var)   # Silverman's variance correction
    ge = 0
    for _ in range(B):
        ys = x[rnd.integers(0, n, n)] + hc * rnd.standard_normal(n)
        ys = mu + (ys - mu) * shrink
        if h_crit(ys) >= hc:
            ge += 1
    return ge / B, hc

def mixture_dbic(x, iters=300):
    """dBIC = BIC(1 component) - BIC(2 components) for a 1-D Gaussian mixture fitted by EM."""
    x = np.asarray(x, dtype=float)
    n = x.size
    if n < 8 or x.std() == 0:
        return float('nan')
    # one component
    mu1, v1 = x.mean(), x.var() + 1e-12
    ll1 = -0.5 * n * (math.log(2 * math.pi * v1) + 1.0)
    bic1 = -2 * ll1 + 2 * math.log(n)
    # two components, initialised on the two halves of the sorted sample
    xs = np.sort(x)
    m = n // 2
    mu = np.array([xs[:m].mean(), xs[m:].mean()])
    var = np.array([xs[:m].var(), xs[m:].var()]) + 1e-10
    w = np.array([0.5, 0.5])
    for _ in range(iters):
        # E
        lp = -0.5 * ((x[:, None] - mu[None, :]) ** 2 / var[None, :] + np.log(2 * math.pi * var[None, :])) + np.log(w[None, :])
        mx = lp.max(axis=1, keepdims=True)
        r = np.exp(lp - mx)
        r /= r.sum(axis=1, keepdims=True)
        # M
        nk = r.sum(axis=0) + 1e-12
        w = nk / n
        mu = (r * x[:, None]).sum(axis=0) / nk
        var = (r * (x[:, None] - mu[None, :]) ** 2).sum(axis=0) / nk + 1e-10
    lp = -0.5 * ((x[:, None] - mu[None, :]) ** 2 / var[None, :] + np.log(2 * math.pi * var[None, :])) + np.log(w[None, :])
    mx = lp.max(axis=1, keepdims=True)
    ll2 = float((mx[:, 0] + np.log(np.exp(lp - mx).sum(axis=1))).sum())
    bic2 = -2 * ll2 + 5 * math.log(n)
    return bic1 - bic2


# ---------------------------------------------------------------------------------------------
# A PARTITION test, because a clean subfamily split need not make the 1-D marginal multimodal:
# on human NPIP every A–A pair beats every A–B pair, yet within-B pairs alone span the whole
# range, so the pooled histogram is one smeared hump and Silverman says "unimodal" (p = 0.125).
# The structure is in WHICH pairs are high — the graph — and a 1-D test discards that.
#
# Statistic: over threshold-induced partitions (components of the subgraph with identity >= t),
# the largest contrast  mean(within-component identity) - mean(between-component identity),
# counting only components with >= MIN_COMP members so sister pairs cannot masquerade as a
# subfamily. Null: permute identities across the edges. The null has EXACTLY the observed
# marginal, so only the assignment of identities to pairs is tested — and the same maximisation
# over t runs inside every permutation, so choosing the best t costs nothing unearned.
# ---------------------------------------------------------------------------------------------
MIN_COMP = 3

def _partition_contrast(edges, idn, n_nodes, thresholds):
    """Best within-minus-between contrast over threshold-induced component partitions."""
    order = np.argsort(-idn)
    best = -1.0
    parent = np.arange(n_nodes)
    def find(x):
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x
    k = 0
    ei = edges[order]; ii = idn[order]
    for t in thresholds:                      # descending
        while k < ii.size and ii[k] >= t:
            a, b = ei[k]
            ra, rb = find(a), find(b)
            if ra != rb:
                parent[ra] = rb
            k += 1
        roots = np.array([find(i) for i in range(n_nodes)])
        _, inv, cnt = np.unique(roots, return_inverse=True, return_counts=True)
        big = cnt >= MIN_COMP
        if big.sum() < 2:
            continue
        comp = inv[edges[:, 0]]; comp2 = inv[edges[:, 1]]
        keep = big[comp] & big[comp2]
        if keep.sum() == 0:
            continue
        same = comp[keep] == comp2[keep]
        if same.sum() == 0 or (~same).sum() == 0:
            continue
        c = idn[keep][same].mean() - idn[keep][~same].mean()
        if c > best:
            best = c
    return best

def partition_perm_p(pairs, B=200, seed=0, n_thr=40):
    """pairs: iterable of (a, b, identity). Returns (p, observed contrast)."""
    names = {}
    E = []; I = []
    for a, b, i in pairs:
        E.append((names.setdefault(a, len(names)), names.setdefault(b, len(names)))); I.append(i)
    E = np.array(E, dtype=int); I = np.array(I, dtype=float)
    n = len(names)
    if E.shape[0] < 6 or n < 2 * MIN_COMP:
        return float('nan'), float('nan')
    thr = np.unique(np.quantile(I, np.linspace(0.02, 0.98, n_thr)))[::-1]
    obs = _partition_contrast(E, I, n, thr)
    if obs < 0:
        return float('nan'), obs
    rnd = np.random.default_rng(seed)
    ge = 0
    for _ in range(B):
        if _partition_contrast(E, rnd.permutation(I), n, thr) >= obs:
            ge += 1
    return ge / B, obs
