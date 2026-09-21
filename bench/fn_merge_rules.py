#!/usr/bin/env python3
"""Post-clustering merge rules against FALSE NEGATIVES, per
`docs/PREREG_false_negative_rules_2026-09-21.md`.

The FN decomposition on held-out chr2/chr8/chr10 (263 truth pairs, 88 false negatives) found the
largest bucket is a GROUPING loss, not an edge loss: for **40 of 88** the edge survives every filter
and MCL then separates the two members. Those are recoverable without touching edge construction.

Three one-line merges over pairs of MCL clusters in the same component:

    R1 edge-count     >= m shipped edges run between A and B
    R2 edge-fraction  edges(A,B) / min(|A|,|B|) >= f
    R3 nbr-jaccard    mean J_N over the A-B edges >= j      (§6u3's metric)

⚠ R3 abstains below component size 5: §6u3 measured J_N at AUC 0.500 — exact chance — at component
size 2, and 60% of truth families are pairs. The restriction is part of the rule, not a knob.

Usage: fn_merge_rules.py --graphs DIR --clusters DIR --chrom chrN --rule r1 --param 3 --out PREFIX
"""
import argparse
import collections


def load(graphs, clusters, chrom):
    adj = collections.defaultdict(set)
    for line in open(f'{graphs}/{chrom}.graph.tsv'):
        f = line.rstrip('\n').split('\t')
        if len(f) == 3 and f[0] != f[1]:
            adj[f[0]].add(f[1]); adj[f[1]].add(f[0])
    node_cl, members = {}, collections.defaultdict(list)
    for line in open(f'{clusters}/{chrom}_fam.clusters.tsv'):
        if line.startswith('cluster_id'):
            continue
        q = line.rstrip('\n').split('\t')
        key = f'{q[5]}:{q[6]}-{q[7]}'
        node_cl[key] = q[0]; members[q[0]].append(key)
    # connected components of the shipped graph, for the size restriction
    seen, comp = set(), {}
    for n0 in adj:
        if n0 in seen:
            continue
        st, mem = [n0], []
        while st:
            x = st.pop()
            if x in seen:
                continue
            seen.add(x); mem.append(x); st.extend(adj[x] - seen)
        for m in mem:
            comp[m] = len(mem)
    return adj, node_cl, members, comp


def main():
    ap = argparse.ArgumentParser()
    for x in ('--graphs', '--clusters', '--chrom', '--rule', '--out'):
        ap.add_argument(x, required=True)
    ap.add_argument('--param', type=float, required=True)
    a = ap.parse_args()

    adj, node_cl, members, comp = load(a.graphs, a.clusters, a.chrom)

    # cross-cluster edges
    cross = collections.defaultdict(list)
    for u in adj:
        cu = node_cl.get(u)
        if not cu:
            continue
        for v in adj[u]:
            cv = node_cl.get(v)
            if not cv or cv == cu:
                continue
            cross[tuple(sorted((cu, cv)))].append((u, v))

    parent = {}

    def find(x):
        parent.setdefault(x, x)
        while parent[x] != x:
            parent[x] = parent[parent[x]]; x = parent[x]
        return x

    merged = 0
    for (ca, cb), es in cross.items():
        n = len(es) // 2 or len(es)          # each edge appears from both endpoints
        ok = False
        if a.rule == 'r1':
            ok = n >= a.param
        elif a.rule == 'r2':
            ok = n / max(1, min(len(members[ca]), len(members[cb]))) >= a.param
        elif a.rule == 'r3':
            sz = max(comp.get(es[0][0], 1), comp.get(es[0][1], 1))
            if sz >= 5:
                js = []
                for u, v in es:
                    nu, nv = adj[u] - {v}, adj[v] - {u}
                    uni = nu | nv
                    js.append(len(nu & nv) / len(uni) if uni else 0.0)
                ok = (sum(js) / len(js)) >= a.param if js else False
        if ok:
            ra, rb = find(ca), find(cb)
            if ra != rb:
                parent[ra] = rb; merged += 1

    out = collections.defaultdict(list)
    for cl, mem in members.items():
        out[find(cl)].extend(mem)
    with open(a.out + '.clusters.tsv', 'w') as fh:
        fh.write('cluster_id\tsize\tdensity\tfrac_in\tcorroborated\tchrom\tstart\tend\n')
        for i, (_, mem) in enumerate(sorted(out.items(), key=lambda kv: -len(kv[1]))):
            if len(mem) < 2:
                continue
            for key in sorted(set(mem)):
                ch, rng = key.rsplit(':', 1); s, e = rng.split('-')
                fh.write(f'FM{i}\t{len(mem)}\tNA\tNA\tNA\t{ch}\t{s}\t{e}\n')
    print(f'{a.chrom} {a.rule}={a.param}: {len(members)} clusters -> '
          f'{sum(1 for v in out.values() if len(v) >= 2)} ({merged} merges)')


if __name__ == '__main__':
    main()
