#!/usr/bin/env python3
"""Port of `annotation_families::mcl` (Rust) for graphs built outside `mcl_families`.

Self-loops of weight 1, column normalisation, expansion M*M, inflation, absolute prune, convergence at 1e-7, at most 100
iterations; clusters = columns joined to their heaviest row (union-find). Verified against the Rust clusters on the
development RefSeq E1 graph: pairwise 0.9999 / 1.000 on shared members (§6km).
"""
import collections

import numpy as np
import scipy.sparse as sp


def mcl(edges, inflation=2.8, prune=1e-9, max_iter=100):
    """edges: {(a, b): weight} over hashable node ids (undirected). Returns clusters (lists of node ids), size >= 1."""
    nodes = sorted({x for e in edges for x in e})
    ix = {n: i for i, n in enumerate(nodes)}
    n = len(nodes)
    if n == 0:
        return []
    rows, cols, vals = list(range(n)), list(range(n)), [1.0] * n
    for (a, b), w in edges.items():
        if a == b:
            continue
        rows += [ix[a], ix[b]]
        cols += [ix[b], ix[a]]
        vals += [w, w]
    M = sp.csc_matrix((vals, (rows, cols)), shape=(n, n))

    def norm(X):
        s = np.asarray(X.sum(axis=0)).ravel()
        s[s == 0] = 1
        return (X @ sp.diags(1 / s)).tocsc()

    M = norm(M)
    for _ in range(max_iter):
        N = (M @ M).tocsc()
        N.data = N.data ** inflation
        N.data[N.data < prune] = 0
        N.eliminate_zeros()
        N = norm(N)
        d = abs(N - M)
        done = d.nnz == 0 or d.max() < 1e-7
        M = N
        if done:
            break
    parent = list(range(n))

    def find(x):
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x
    for j in range(n):
        s, e = M.indptr[j], M.indptr[j + 1]
        if s == e:
            continue
        r = M.indices[s:e][np.argmax(M.data[s:e])]
        ra, rb = find(j), find(int(r))
        if ra != rb:
            parent[max(ra, rb)] = min(ra, rb)
    groups = collections.defaultdict(list)
    for i in range(n):
        groups[find(i)].append(nodes[i])
    return list(groups.values())
