#!/usr/bin/env python3
"""Pool graph and clusters (docs/PREREG_unmapped_rescue_2026-10-08.md section 1, steps 1-2)."""
import collections


def edges(paf_lines, delta, min_frac):
    """read-read edges of an all-vs-all PAF (minimap2 -x ava-pb -c): the overlap's gap-compressed divergence de <= delta and the alignment block
    covers >= min_frac of the SHORTER read. Self hits and records without a de tag are dropped. Yields (a, b) with a < b (duplicates possible)."""
    for ln in paf_lines:
        f = ln.rstrip("\n").split("\t")
        if len(f) < 12 or f[0] == f[5]:
            continue
        de = next((float(x[5:]) for x in f[12:] if x.startswith("de:f:")), None)
        if de is None or de > delta:
            continue
        if int(f[10]) >= min_frac * min(int(f[1]), int(f[6])):
            yield (f[0], f[5]) if f[0] < f[5] else (f[5], f[0])


def components(edge_list, nodes):
    """{node: component id}: union-find over the edges; nodes without an edge are singletons"""
    parent = {n: n for n in nodes}

    def find(x):
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x
    for a, b in edge_list:
        if a in parent and b in parent:
            ra, rb = find(a), find(b)
            if ra != rb:
                parent[ra] = rb
    return {n: find(n) for n in nodes}


def clusters(comp, min_size=3):
    """{component id: [reads]} for components with >= min_size reads, reads sorted"""
    by = collections.defaultdict(list)
    for n, c in comp.items():
        by[c].append(n)
    return {c: sorted(v) for c, v in by.items() if len(v) >= min_size}
