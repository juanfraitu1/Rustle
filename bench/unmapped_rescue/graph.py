#!/usr/bin/env python3
"""Pool graph and clusters (docs/PREREG_unmapped_rescue_2026-10-08.md section 1, steps 1-2)."""
import collections


def proper_overlap(qlen, qs, qe, strand, tlen, ts, te, blk, delta):
    """Amendment 22: at each end of the alignment at least one of the two reads is reached (unaligned remainder <= delta x block): a containment or a
    suffix-prefix overlap, not a shared middle. Reverse strand: the query start pairs with the target end."""
    tol = delta * blk
    if strand == "+":
        e1, e2 = min(qs, ts), min(qlen - qe, tlen - te)
    else:
        e1, e2 = min(qs, tlen - te), min(qlen - qe, ts)
    return e1 <= tol and e2 <= tol


def edges(paf_lines, delta, min_frac, proper=False, edit=False):
    """read-read edges of an all-vs-all PAF (minimap2 -x ava-pb -c): the overlap's gap-compressed divergence de <= delta and the alignment block
    covers >= min_frac of the SHORTER read (with proper=True also a proper overlap, Amendment 22; with edit=True the divergence is NM / block, Amendment 24). Self hits and records without a de tag are dropped. Yields (a, b) with a < b (duplicates possible)."""
    for ln in paf_lines:
        f = ln.rstrip("\n").split("\t")
        if len(f) < 12 or f[0] == f[5]:
            continue
        de = next((float(x[5:]) for x in f[12:] if x.startswith("de:f:")), None)
        if de is None or de > delta:
            continue
        if edit:        # Amendment 24: every gap counts, NM / block (the gap-compressed de hides an exon-scale deletion)
            nm = next((int(x[5:]) for x in f[12:] if x.startswith("NM:i:")), None)
            if nm is None or nm > delta * int(f[10]):
                continue
        if int(f[10]) >= min_frac * min(int(f[1]), int(f[6])):
            if proper and not proper_overlap(int(f[1]), int(f[2]), int(f[3]), f[4], int(f[6]), int(f[7]), int(f[8]), int(f[10]), delta):
                continue
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
