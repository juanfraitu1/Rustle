#!/usr/bin/env python3
"""Row 1662 (docs/PREREG_family_join_2026-11-13.md): join MCL families through strong exon-sharing edges.

  famjoin.py BASE.fam.copies.tsv LOCI.tsv STRONG_GRAPH.tsv OUT.fam.copies.tsv
  famjoin.py --selftest
STRONG_GRAPH = `mcl_families --min-exon-share 0.75 --dump-graph` (u, v, weight; nodes = locus representatives
`chrom:start1-end`). Two families with an edge between their loci become one (transitive); other rows are unchanged.
"""
import sys


def join(fam_of_node, edges):
    """{node: family}, [(u, v)] -> {family: joined family} (union-find; the smallest id names the join)."""
    par = {f: f for f in set(fam_of_node.values())}
    def find(x):
        while par[x] != x:
            par[x] = par[par[x]]; x = par[x]
        return x
    for u, v in edges:
        a, b = fam_of_node.get(u), fam_of_node.get(v)
        if a is None or b is None:
            continue
        a, b = find(a), find(b)
        if a != b:
            par[max(a, b)] = min(a, b)
    return {f: find(f) for f in par}


def main(base, loci, graph, out):
    rep = dict(l.rstrip('\n').split('\t')[:2] for l in open(loci))
    rows = [l.rstrip('\n').split('\t') for l in open(base)]
    h = rows[0]; c, s, e = h.index('chrom'), h.index('locus_start'), h.index('locus_end')
    node = lambda r: rep.get(f'{r[c]}:{int(r[s]) + 1}-{r[e]}', f'{r[c]}:{int(r[s]) + 1}-{r[e]}')
    fam_of_node = {node(r): r[0] for r in rows[1:]}
    edges = [l.split('\t')[:2] for l in open(graph)]
    m = join(fam_of_node, edges)
    with open(out, 'w') as w:
        w.write('\t'.join(h) + '\n')
        for r in rows[1:]:
            w.write('\t'.join([m[r[0]]] + r[1:]) + '\n')
    size = lambda fams: max(fams.count(f) for f in set(fams))
    b, a = [r[0] for r in rows[1:]], [m[r[0]] for r in rows[1:]]
    print(f'{base}: families {len(set(b))} -> {len(set(a))}, largest {size(b)} -> {size(a)}, '
          f'strong edges between two families {sum(1 for u, v in edges if u in fam_of_node and v in fam_of_node and fam_of_node[u] != fam_of_node[v])}')


def selftest():
    m = join({'a': 'F1', 'b': 'F2', 'c': 'F3', 'd': 'F4'}, [('a', 'b'), ('b', 'c'), ('d', 'x'), ('y', 'z')])
    assert m == {'F1': 'F1', 'F2': 'F1', 'F3': 'F1', 'F4': 'F4'}, m   # transitive; a locus in no family joins nothing
    print('selftest ok')


if __name__ == '__main__':
    selftest() if '--selftest' in sys.argv else main(*sys.argv[1:5])
