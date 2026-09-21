#!/usr/bin/env python3
"""Clustering-operator bakeoff on HELD-OUT substrate, per
`docs/PREREG_operator_heldout_2026-09-20.md`.

`bench/CLUSTERING_OPERATOR_BAKEOFF.md` (§6n5) ran these operators on NPIP and TBC1D3 and warned that
picking a winner there is the dev-set selection trap. This runs the same question on the SHIPPED DNA
gene-body graph (`mcl_families --dump-graph`, so the edge set and exon conjunct are the shipped ones)
for chromosomes with zero ledger exposure, scored against Soto's published families.

§6n5's lesson was that the F column hides the cost — every triangle-based operator dissolved all 13
two-copy components, and 2 is the modal family size — so node coverage and two-member groups are
reported beside F, never instead of it.

Usage: operator_bakeoff_heldout.py --graphs DIR --chroms chr2,chr8,chr10 --out DIR
"""
import argparse
import collections
import os
import sys

import networkx as nx

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import mcl_port  # noqa: E402


def read_graph(path):
    g = nx.Graph()
    for line in open(path):
        f = line.rstrip('\n').split('\t')
        if len(f) == 3:
            a, b, w = f[0], f[1], float(f[2])
            if a == b:
                g.add_node(a)
            else:
                g.add_edge(a, b, weight=w)
        elif len(f) == 2:
            g.add_node(f[0])
    return g


def groups(op, g):
    """-> list of node sets, each of size >= 1."""
    if op == 'components':
        return [set(c) for c in nx.connected_components(g)]
    if op.endswith('-truss'):
        k = int(op[0])
        h = nx.k_truss(g, k)
        return [set(c) for c in nx.connected_components(h)]
    if op.endswith('-clique'):
        k = int(op[0])
        return [set(c) for c in nx.algorithms.community.k_clique_communities(g, k)]
    if op == 'louvain':
        return [set(c) for c in nx.algorithms.community.louvain_communities(g, seed=0)]
    if op == 'greedy-modularity':
        return [set(c) for c in nx.algorithms.community.greedy_modularity_communities(g)]
    if op == 'label-prop':
        return [set(c) for c in nx.algorithms.community.label_propagation_communities(g)]
    if op == 'mcl':
        edges = {(u, v): d.get('weight', 1.0) for u, v, d in g.edges(data=True)}
        return [set(c) for c in mcl_port.mcl(edges, inflation=2.8)]
    raise SystemExit(f'unknown operator {op}')


OPS = ['components', '3-truss', '4-truss', '3-clique', '4-clique',
       'louvain', 'greedy-modularity', 'label-prop', 'mcl']


def main():
    ap = argparse.ArgumentParser()
    for x in ('--graphs', '--chroms', '--out'):
        ap.add_argument(x, required=True)
    a = ap.parse_args()
    os.makedirs(a.out, exist_ok=True)
    chroms = a.chroms.split(',')

    stats = collections.defaultdict(lambda: [0, 0, 0, 0])   # op -> [nodes, covered, groups>=2, pairs]
    for c in chroms:
        g = read_graph(f'{a.graphs}/{c}.graph.tsv')
        for op in OPS:
            try:
                gs = [s for s in groups(op, g)]
            except Exception as e:                           # noqa: BLE001
                print(f'  {c} {op}: FAILED ({e})'); continue
            keep = [s for s in gs if len(s) >= 2]
            st = stats[op]
            st[0] += g.number_of_nodes()
            st[1] += sum(len(s) for s in keep)
            st[2] += len(keep)
            st[3] += sum(1 for s in keep if len(s) == 2)
            with open(f'{a.out}/{c}_{op}.clusters.tsv', 'w') as fh:
                fh.write('cluster_id\tsize\tdensity\tfrac_in\tcorroborated\tchrom\tstart\tend\n')
                for i, s in enumerate(sorted(keep, key=lambda x: -len(x))):
                    for name in sorted(s):
                        try:
                            ch, rng = name.rsplit(':', 1); st_, e_ = rng.split('-')
                        except ValueError:
                            continue
                        fh.write(f'OP{i}\t{len(s)}\tNA\tNA\tNA\t{ch}\t{int(st_)-1}\t{int(e_)}\n')
    print(f"{'operator':20s} {'groups>=2':>10} {'nodes covered':>14} {'coverage':>9} {'2-member groups':>16}")
    for op in OPS:
        n, cov, grp, two = stats[op]
        print(f'{op:20s} {grp:>10} {cov:>14} {100*cov/n if n else 0:>8.1f}% {two:>16}')


if __name__ == '__main__':
    main()
