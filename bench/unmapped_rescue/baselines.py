#!/usr/bin/env python3
"""Single-read comparators of docs/PREREG_unmapped_rescue_2026-10-08.md section 1: B0 (A13's rule) and B1 (single-read dc-megablast cover score)."""
import collections
import random

try:
    from attribute import family_of
except ImportError:                      # pragma: no cover
    from .attribute import family_of

COV_MIN = 0.5      # A13 read coverage floor
DE_MAX = 0.20      # A13 divergence ceiling


def sample(reads, k, seed=1):
    """seeded sample of k reads (sorted first so the draw does not depend on input order)"""
    rs = sorted(reads)
    return sorted(random.Random(seed).sample(rs, k)) if len(rs) > k else rs


def single_read_nucleotide(paf_lines):
    """B0: minimap2 map-ont of each read to the family copies, read coverage >= 0.5 and de <= 0.20 (A13's floors); the read goes to the family of its
    hit with the most matching bases, None when two families tie. -> {read: family | None}; a read without a usable hit is absent"""
    best = collections.defaultdict(dict)
    for ln in paf_lines:
        f = ln.rstrip("\n").split("\t")
        de = next((float(x[5:]) for x in f[12:] if x.startswith("de:f:")), None)
        if de is None or de > DE_MAX or (int(f[3]) - int(f[2])) / int(f[1]) < COV_MIN:
            continue
        fam = family_of(f[5])
        best[f[0]][fam] = max(best[f[0]].get(fam, 0), int(f[9]))
    out = {}
    for q, fams in best.items():
        ranked = sorted(fams.items(), key=lambda kv: -kv[1])
        out[q] = ranked[0][0] if len(ranked) == 1 or ranked[0][1] > ranked[1][1] else None
    return out
