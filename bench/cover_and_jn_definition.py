#!/usr/bin/env python3
"""Two candidate improvements to the family definition, per
`docs/PREREG_cover_and_jn_weight_2026-09-21.md` (md5 `eed4b0fc`).

**TEST 1 — the prediction may emit a COVER.** O1 emits a strict partition (0 of 2,670 loci in >1
cluster), so a fusion gene is inexpressible. Register 845 refuted dual membership on the TRUTH side and
left an explicit re-open condition, now met (chr10 AGAP 40% chimeras). After MCL, node v joins cluster C
iff v has >= k edges into C -- purely combinatorial, k=2 being triangle support (§6kd).

**TEST 2 — neighbourhood Jaccard as an MCL edge WEIGHT, size-gated.** J_N is the project's strongest
separator (AUC 0.924 at component size >= 10) and the only one to survive size-residualisation, but was
only ever judged as a post-clustering merge (register 933, refuted). Register 916/917: a scalar cannot
be judged at the operator's weakest setting. w' = w * (1 + J_N) inside components >= gate, else w.
The gate is part of the rule: J_N is AUC 0.500 -- exact chance -- at component size 2.

⚠TRUTH IS THE COVER. Soto's S1C has an explicit `No. Assigned Families` column and 149/2,334 gene IDs
are multi-family; the derived set size matches that column for 2,334/2,334. Every prior scorer dropped
those genes, which also drops whole families below the >= 3 floor (chr16: 15 fams/70 members as a cover
vs 8/43 as a partition). Scoring a cover prediction against a partition-restricted truth is rigged.

⚠The comparator is mcl_port-MCL ON THE SAME GRAPH, never the shipped Rust F (register 917: mcl_port is
not bit-identical to the Rust MCL).
"""
import argparse
import collections
import csv
import os
import re
import sys

import numpy as np
from scipy.optimize import linear_sum_assignment

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import mcl_port  # noqa: E402


def read_graph(path):
    adj = collections.defaultdict(dict)
    nodes = set()
    for line in open(path):
        f = line.rstrip('\n').split('\t')
        if len(f) == 3 and f[0] != f[1]:
            a, b, w = f[0], f[1], float(f[2])
            adj[a][b] = max(adj[a].get(b, 0.0), w)
            adj[b][a] = max(adj[b].get(a, 0.0), w)
            nodes.add(a); nodes.add(b)
        elif f:
            nodes.add(f[0])
    return adj, nodes


def components(adj, nodes):
    seen, comp = set(), {}
    for n0 in nodes:
        if n0 in seen:
            continue
        stack, mem = [n0], []
        while stack:
            x = stack.pop()
            if x in seen:
                continue
            seen.add(x); mem.append(x)
            stack.extend(k for k in adj.get(x, ()) if k not in seen)
        for m in mem:
            comp[m] = len(mem)
    return comp


def gene_names(gff, chrom):
    out = {}
    for line in open(gff):
        if line.startswith('#'):
            continue
        f = line.rstrip('\n').split('\t')
        if len(f) < 9 or f[0] != chrom or f[2] not in ('gene', 'pseudogene'):
            continue
        m = re.search(r'Name=([^;]+)', f[8])
        if m:
            out[f'{f[0]}:{f[3]}-{f[4]}'] = m.group(1)
    return out


def soto_cover(s1c):
    """name -> set(Family ID)  (the COVER), plus the set of names that are ambiguous by gene ID."""
    byname = collections.defaultdict(set)
    name_ids = collections.defaultdict(set)
    byid = collections.defaultdict(set)
    for r in csv.DictReader(open(s1c), delimiter='\t'):
        fid = (r.get('Family ID') or '').strip()
        nm = (r.get('Gene Name') or '').strip()
        gid = (r.get('Gene ID') or '').strip()
        if not (fid and fid != 'N/A' and nm and gid):
            continue
        byname[nm].add(fid); name_ids[nm].add(gid); byid[gid].add(fid)
    ambiguous = {n for n, g in name_ids.items()
                 if len(g) > 1 and len({frozenset(byid[x]) for x in g}) > 1}
    return byname, ambiguous


def truth_families(cover, on_chrom, drop):
    fam = collections.defaultdict(set)
    for nm, fids in cover.items():
        if nm not in on_chrom or nm in drop:
            continue
        for f in fids:
            fam[f].add(nm)
    return {f: sorted(v) for f, v in fam.items() if len(v) >= 3}


def score(truth, clusters):
    """one-to-one bipartite on overlap; unmatched truth families score 0 and stay in the mean."""
    troots, cids = sorted(truth), sorted(clusters)
    if not troots:
        return []
    ov = np.zeros((len(troots), len(cids)), dtype=int)
    for i, r in enumerate(troots):
        t = set(truth[r])
        for j, c in enumerate(cids):
            ov[i, j] = len(t & set(clusters[c]))
    rows, cols = linear_sum_assignment(-ov) if cids else ([], [])
    matched = {troots[i]: cids[j] for i, j in zip(rows, cols) if ov[i, j] > 0}
    out = []
    for r in troots:
        t = set(truth[r]); c = matched.get(r)
        if c is None:
            out.append((r, len(t), 0.0, 0.0, 0.0)); continue
        p = set(clusters[c]); hit = len(t & p)
        sens = hit / len(t); prec = hit / len(p) if p else 0.0
        f = 0.0 if sens + prec == 0 else 2 * sens * prec / (sens + prec)
        out.append((r, len(t), sens, prec, f))
    return out


def mcl_clusters(adj, nodes):
    edges = {}
    for u in adj:
        for v, w in adj[u].items():
            if u < v:
                edges[(u, v)] = w
    groups = mcl_port.mcl(edges, inflation=2.8) if edges else []
    out = {}
    for i, g in enumerate(sorted((set(x) for x in groups), key=lambda s: -len(s))):
        if len(g) >= 2:
            out[f'C{i}'] = sorted(g)
    return out


def jn_weighted(adj, nodes, gate):
    comp = components(adj, nodes)
    new = collections.defaultdict(dict)
    for u in adj:
        for v, w in adj[u].items():
            if u >= v:
                continue
            if comp.get(u, 1) >= gate:
                nu = set(adj[u]) - {v}; nv = set(adj[v]) - {u}
                uni = nu | nv
                j = len(nu & nv) / len(uni) if uni else 0.0
                w = w * (1.0 + j)
            new[u][v] = w; new[v][u] = w
    return new


def add_cover(clusters, adj, k):
    """node v joins cluster C (v not in C) iff v has >= k edges into C."""
    member = collections.defaultdict(set)
    for cid, mem in clusters.items():
        for m in mem:
            member[m].add(cid)
    out = {cid: set(mem) for cid, mem in clusters.items()}
    added = 0
    for v in list(adj):
        cnt = collections.Counter()
        for u in adj[v]:
            for cid in member.get(u, ()):
                cnt[cid] += 1
        for cid, n in cnt.items():
            if cid not in member.get(v, ()) and n >= k:
                out[cid].add(v); added += 1
    return {c: sorted(m) for c, m in out.items()}, added


def run(a):
    cover, ambiguous = soto_cover(a.soto)
    drop = ambiguous if a.drop_ambiguous else set()
    arms = collections.defaultdict(list)
    meta = collections.defaultdict(lambda: [0, 0, 0])   # arm -> [nodes, covered, 2-member groups]
    for c in a.chroms.split(','):
        adj, nodes = read_graph(f'{a.graphs}/{c}.graph.tsv')
        names = gene_names(a.gff, c)
        on = {names[n] for n in nodes if n in names}
        truth = truth_families(cover, set(names.values()), drop)
        if not truth:
            print(f'  {c}: no truth family >= 3 members, skipped', file=sys.stderr); continue

        def emit(label, cl):
            named = {cid: [names.get(x, x) for x in mem] for cid, mem in cl.items()}
            arms[label] += score(truth, named)
            m = meta[label]
            m[0] += len(nodes)
            m[1] += len({x for mem in cl.values() for x in mem})
            m[2] += sum(1 for mem in cl.values() if len(mem) == 2)

        base = mcl_clusters(adj, nodes)
        emit('baseline MCL', base)
        if a.test in ('cover', 'both'):
            for k in (2, 3, 4):
                cl, _ = add_cover(base, adj, k)
                emit(f'cover k>={k}', cl)
        if a.test in ('jn', 'both'):
            for g in (3, 5, 8):
                emit(f'J_N gate>={g}', mcl_clusters(jn_weighted(adj, nodes, g), nodes))
        print(f'  {c}: {len(nodes)} nodes, {len(truth)} truth families '
              f'({sum(len(v) for v in truth.values())} members)', file=sys.stderr)

    print(f"\ncover truth | chroms {a.chroms} | ambiguous names "
          f"{'DROPPED' if a.drop_ambiguous else 'kept'} ({len(ambiguous)})")
    print(f"  {'arm':16s} {'fams':>5} {'sens':>7} {'prec':>7} {'F':>7} {'dF':>8} "
          f"{'matched':>8} {'3-mem hit':>10} {'cov%':>6} {'2-grp':>6}")
    bF = None
    for label in ['baseline MCL'] + [k for k in arms if k != 'baseline MCL']:
        rows = arms[label]
        if not rows:
            continue
        F = sum(r[4] for r in rows) / len(rows)
        S = sum(r[2] for r in rows) / len(rows)
        P = sum(r[3] for r in rows) / len(rows)
        if bF is None:
            bF = F
        mt = sum(1 for r in rows if r[4] > 0)
        m3 = sum(1 for r in rows if r[1] == 3 and r[4] > 0)
        n, cov, two = meta[label]
        print(f"  {label:16s} {len(rows):>5} {S:>7.4f} {P:>7.4f} {F:>7.4f} {F-bF:>+8.4f} "
              f"{mt:>8} {m3:>10} {100*cov/n if n else 0:>5.1f}% {two:>6}")
    if a.per_family:
        base = {r[0]: r[4] for r in arms['baseline MCL']}
        print("\n  per-family movement vs baseline (register 917: a pooled gain can be 2 families)")
        for label in [k for k in arms if k != 'baseline MCL']:
            up = dn = 0; tot = 0.0
            worst = None
            for r in arms[label]:
                d = r[4] - base.get(r[0], 0.0)
                tot += d
                if d > 1e-9: up += 1
                elif d < -1e-9:
                    dn += 1
                    if worst is None or d < worst[1]:
                        worst = (r[0], d)
            print(f"    {label:16s} better {up:>3} | worse {dn:>3} | unchanged "
                  f"{len(arms[label])-up-dn:>3} | sum dF {tot:>+7.4f}"
                  + (f" | worst {worst[0]} {worst[1]:+.3f}" if worst else ""))


def main():
    ap = argparse.ArgumentParser()
    for x in ('--graphs', '--gff', '--soto', '--chroms'):
        ap.add_argument(x, required=True)
    ap.add_argument('--test', default='both', choices=['cover', 'jn', 'both'])
    ap.add_argument('--per-family', action='store_true',
                    help='per-family movement vs baseline -- register 917: a pooled gain can be 2 families')
    ap.add_argument('--drop-ambiguous', action='store_true',
                    help='sensitivity arm: drop the names that are multi-family only by gene-ID collision')
    run(ap.parse_args())


if __name__ == '__main__':
    main()
