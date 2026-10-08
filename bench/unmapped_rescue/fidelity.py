#!/usr/bin/env python3
"""Does a cluster's consensus reproduce the erased copy? (prereg section 4, candidate fidelity). The consensus is aligned to the UNMASKED genome (splice:hq);
truth = the best hit overlaps the erased interval of the cluster's majority family. Reported: identity, query coverage and the share at identity >= 0.999."""
import csv
import json
import sys

import score as S


def best_hits(paf_lines):
    """{query: (identity, coverage, ref, start, end)} the best record per query by matches"""
    b = {}
    for ln in paf_lines:
        f = ln.rstrip("\n").split("\t")
        q = f[0].split("|")[0]
        m = int(f[9])
        if q not in b or m > b[q][5]:
            b[q] = (m / max(1, int(f[10])), (int(f[3]) - int(f[2])) / int(f[1]), f[5], int(f[7]), int(f[8]), m)
    return {q: v[:5] for q, v in b.items()}


def fidelity(clusters, truth, hits, erased, only=None):
    """erased: {family: (chrom, start, end)}; only: restrict to these cluster ids (e.g. the attributed ones).
    -> clusters_with_family (truth majority), on_erased_copy (best hit overlaps the erased interval; a property of the truth labels, not a result by itself),
    identity_ge_0_999, identity_x_coverage_ge_0_999 (the registered metric, among those on the erased copy)"""
    n = on = hi = reg = 0
    for c, rs in clusters.items():
        if only is not None and c not in only:
            continue
        fam = S.majority(rs, truth)
        if fam is None or fam not in erased or c not in hits:
            continue
        n += 1
        ident, cov, ref, a, b = hits[c]
        e = erased[fam]
        if ref == e[0] and a < e[2] and e[1] < b:
            on += 1
            hi += ident >= 0.999
            reg += ident * cov >= 0.999
    return dict(clusters_with_family=n, on_erased_copy=on, identity_ge_0_999=hi, identity_x_coverage_ge_0_999=reg)
