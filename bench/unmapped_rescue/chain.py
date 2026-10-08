#!/usr/bin/env python3
"""Chain-aware refinement of the components of the pool graph (docs/PREREG_unmapped_rescue_2026-10-08.md Amendment 23).

A component is a chain when a short read is CONTAINED in two long reads that disagree with each other: both containments are proper overlaps, so single-linkage
joins the long reads. Star clustering inside the component fixes it: longest-first, a read joins only a cluster whose representative it is compatible with."""
import subprocess

import graph as G

MAX_COMPONENT = 60


def star_clusters(lens, compat, min_size=3):
    """lens {read: length}; compat = set of sorted (a, b) pairs. Longest-first: the longest unassigned read is the representative, its cluster = itself + the
    unassigned reads compatible with it. Clusters below min_size are dropped. -> [[representative, members...]]"""
    order = sorted(lens, key=lambda r: (-lens[r], r))
    left, out = set(order), []
    for rep in order:
        if rep not in left:
            continue
        members = [r for r in order if r in left and r != rep and ((rep, r) if rep < r else (r, rep)) in compat]
        cluster = [rep] + members
        left -= set(cluster)
        if len(cluster) >= min_size:
            out.append(cluster)
    return out


def refine(clusters, seqs, allvsall_fn, delta, min_size=3, max_component=MAX_COMPONENT):
    """clusters {id: [reads]} -> {id: [reads]}: components of at most max_component reads are re-clustered by star clustering over the pairs that pass the
    Amendment 22 edge (allvsall_fn({read: seq}) -> PAF lines); larger components are kept. A refined cluster is named after its representative."""
    out = {}
    for cid, reads in clusters.items():
        if len(reads) > max_component:
            out[cid] = list(reads)
            continue
        sub = {r: seqs[r] for r in reads}
        compat = set(G.edges(allvsall_fn(sub), delta, 0.5, proper=True))
        for cl in star_clusters({r: len(s) for r, s in sub.items()}, compat, min_size):
            out[cl[0]] = cl
    return out


def minimap_allvsall(work, threads=2):
    import os
    os.makedirs(work, exist_ok=True)

    def f(reads):
        p = f"{work}/set.fa"
        with open(p, "w") as o:
            for n, s in reads.items():
                o.write(f">{n}\n{s}\n")
        return subprocess.run(f"minimap2 -x map-hifi -c -N 20 -t {threads} {p} {p}", shell=True, stdout=subprocess.PIPE, stderr=subprocess.DEVNULL, text=True).stdout.splitlines()
    return f
