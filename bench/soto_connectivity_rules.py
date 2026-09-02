#!/usr/bin/env python3
"""Combinatorial co-membership criteria, scored against Soto 2025's families (ledger §6bz).

§6by shipped `distance`, which is a *reachability* statement. This asks whether a criterion with a
cleaner combinatorial character does better — the shape the advisor asks for: a named quantity on a
graph with a theorem behind it, not a tuned threshold.

The candidate:

    LOCAL EDGE CONNECTIVITY  lambda(u,v) = the minimum number of E_r edges whose removal separates
    u from v = (Menger) the maximum number of EDGE-DISJOINT paths joining them.

    lambda >= 1 : connected — what the catalog asserts today (transitive closure)
    lambda >= 2 : NO SINGLE ALIGNMENT RECORD'S LOSS SEPARATES THEM

`lambda(u,v) >= 2` is the per-pair form of the family-level `cut_certified` the catalog already
reports, and it says something an identity threshold cannot: the assertion does not rest on one
alignment. Computed by unit-capacity max-flow (BFS augmenting paths); families are <= 54 copies so
this is trivial.

⛔ **Identity-keyed rules are UNADJUDICABLE on this truth set** — Soto is SD98 by construction, so
"raise the floor" wins for the same reason the truth set exists (§6by). Nothing below uses identity.
⚠ Precision is UNDERSTATED throughout (Soto is CAT-bounded); the RATIO between strata is the robust
quantity, never the level.

Usage: soto_connectivity_rules.py <bed> <arm_dir>
"""
import csv
import sys
import os
import math
import collections
import itertools


def wilson(k, n, z=1.96):
    if n == 0:
        return (float("nan"), float("nan"))
    p = k / n
    d = 1 + z * z / n
    c = (p + z * z / (2 * n)) / d
    h = z * math.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / d
    return (max(0.0, c - h), min(1.0, c + h))


def max_flow_unit(adj, s, t, cap_limit=4):
    """Edge-disjoint paths between s and t, capped at `cap_limit` (we only need >=2 vs <2).

    Unit-capacity Edmonds-Karp on the undirected graph: each undirected edge becomes a residual
    pair, and each augmenting BFS path consumes one unit.
    """
    if s == t:
        return 0
    res = collections.defaultdict(int)
    for u in adj:
        for v in adj[u]:
            res[(u, v)] = 1
    flow = 0
    while flow < cap_limit:
        prev = {s: None}
        q = collections.deque([s])
        while q and t not in prev:
            u = q.popleft()
            for v in adj[u]:
                if v not in prev and res[(u, v)] > 0:
                    prev[v] = u
                    q.append(v)
        if t not in prev:
            break
        v = t
        while prev[v] is not None:
            u = prev[v]
            res[(u, v)] -= 1
            res[(v, u)] += 1
            v = u
        flow += 1
    return flow


def main():
    bed, arm = sys.argv[1], sys.argv[2]
    genes, fam = {}, {}
    for line in open(bed):
        f = line.rstrip("\n").split("\t")
        genes[f[3]] = (f[0], int(f[1]), int(f[2]))
        fam[f[3]] = f[3].split("|")[1]

    copies = list(csv.DictReader(open(os.path.join(arm, "cat.copies.tsv")), delimiter="\t"))
    a = {}
    for i, r in enumerate(copies):
        c, s, e = r["chrom"], int(r["start"]), int(r["end"])
        best = None
        for n, (gc, gs, ge) in genes.items():
            if gc != c:
                continue
            ov = min(ge, e) - max(gs, s)
            if ov <= 0:
                continue
            fr = ov / max(1, ge - gs)
            if best is None or fr > best[0]:
                best = (fr, n)
        if best:
            a[i] = best[1]

    byfam = collections.defaultdict(list)
    for g, f in fam.items():
        byfam[f].append(g)
    truth = {tuple(sorted(p)) for gs in byfam.values() for p in itertools.combinations(sorted(gs), 2)}

    key = lambda r: (r["chrom"], int(r["start"]), int(r["end"]))
    adj_all = collections.defaultdict(set)
    for r in csv.DictReader(open(os.path.join(arm, "dump", "e.edges.tsv")), delimiter="\t"):
        ki = (r["chrom_i"], int(r["start_i"]), int(r["end_i"]))
        kj = (r["chrom_j"], int(r["start_j"]), int(r["end_j"]))
        adj_all[ki].add(kj)
        adj_all[kj].add(ki)

    fams = collections.defaultdict(list)
    for i, r in enumerate(copies):
        fams[r["family_id"]].append((key(r), a.get(i)))

    rows = []   # (pair, dist, lam, common_nbrs, famsize)
    for fid, mem in fams.items():
        S = {k for k, _ in mem}
        adj = {k: (adj_all[k] & S) for k in S}
        for (ku, gu), (kv, gv) in itertools.combinations(mem, 2):
            if gu is None or gv is None or gu == gv:
                continue
            # BFS distance
            seen = {ku: 0}
            q = collections.deque([ku])
            while q:
                x = q.popleft()
                for y in adj[x]:
                    if y not in seen:
                        seen[y] = seen[x] + 1
                        q.append(y)
            d = seen.get(kv)
            if d is None:
                continue
            lam = max_flow_unit(adj, ku, kv)
            cn = len(adj[ku] & adj[kv])
            rows.append((tuple(sorted((gu, gv))), d, lam, cn, len(mem)))

    def rep(label, sub, base_tp=None):
        tp = sum(1 for r in sub if r[0] in truth)
        n = len(sub)
        lo, hi = wilson(tp, n)
        kept = f"  keeps {tp}/{base_tp} = {tp/base_tp:6.1%} of TP" if base_tp else ""
        print(f"  {label:46s} n={n:5d} TP={tp:5d} P={tp/n if n else float('nan'):.4f} "
              f"[{lo:.4f},{hi:.4f}]{kept}")
        return tp

    print(f"arm={os.path.basename(arm.rstrip('/'))}   judgeable pairs={len(rows)}")
    base = rep("BASELINE: every pair the catalog asserts", rows)
    print()
    rep("d = 1 (direct edge)  [§6by, shipped]", [r for r in rows if r[1] == 1], base)
    print()
    print("  --- LOCAL EDGE CONNECTIVITY (Menger: max edge-disjoint paths) ---")
    for k in (1, 2, 3):
        rep(f"lambda(u,v) >= {k}", [r for r in rows if r[2] >= k], base)
    print()
    print("  --- decomposition: is lambda>=2 doing more than d=1? ---")
    rep("d = 1 AND lambda >= 2", [r for r in rows if r[1] == 1 and r[2] >= 2], base)
    rep("d = 1 AND lambda == 1  (edge is a BRIDGE)", [r for r in rows if r[1] == 1 and r[2] == 1], base)
    rep("d >= 2 AND lambda >= 2", [r for r in rows if r[1] >= 2 and r[2] >= 2], base)
    print()
    print("  --- triangle / common neighbours ---")
    rep("common neighbours >= 1", [r for r in rows if r[3] >= 1], base)
    rep("d = 1 AND common neighbours >= 1", [r for r in rows if r[1] == 1 and r[3] >= 1], base)
    print()
    print("  --- CONTROL: does lambda>=2 survive family-size stratification? ---")
    for lo_, hi_, lab in ((2, 5, "2-5"), (6, 15, "6-15"), (16, 10 ** 9, "16+")):
        s2 = [r for r in rows if lo_ <= r[4] <= hi_ and r[2] >= 2]
        s1 = [r for r in rows if lo_ <= r[4] <= hi_ and r[2] < 2]
        def P(s):
            tp = sum(1 for r in s if r[0] in truth)
            return (tp / len(s) if s else float("nan")), len(s)
        p2, n2 = P(s2)
        p1, n1 = P(s1)
        print(f"  family {lab:>5s}:  lambda>=2 n={n2:5d} P={p2:.4f}   lambda<2 n={n1:5d} P={p1:.4f}")


if __name__ == "__main__":
    main()
