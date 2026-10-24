#!/usr/bin/env python3
"""Liftoff copy pairs against matched decoys: the statistics, as pure functions.

Prereg docs/archive/2026-10/merged/PREREG_liftoff_pairs_null_2026-10-14.md. Loci, pairs and the coverage rule are figure 8's (`figures/_liftoff.py`;
prereg docs/archive/2026-09/PREREG_liftoff_loci_2026-09-25.md): a Liftoff locus is covered by a family member when the
member covers >= 50% of its exon bases.

A pair (a, c) = a source gene's annotated placement a and one of its extra copies c (`_liftoff.copy_pairs`, read-supported).
It is SHARED when one of our families covers both. The null replaces c by a decoy x: a read-supported Liftoff locus of another
source gene, not overlapping a, matched to c on (1) whether a family member covers it and (2) its distance bin from a (another
contig; or the same contig at < 100 kb, 100 kb-1 Mb, 1-10 Mb, >= 10 Mb). When no decoy matches, the match is widened: level 1 =
the same contig relation at any distance, level 2 = any contig. A pair without a decoy at level 2 leaves the test, on both sides.
"""
import bisect
import collections
import math
import os
import random
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "..", "figures"))
import _liftoff as LF  # noqa: E402

BINS = (100_000, 1_000_000, 10_000_000)
BIN_NAMES = ("lt100k", "100k_1M", "1M_10M", "ge10M")
OTHER = "other_contig"
ELIGIBLE = ("in_place", "moved", "extra_copy")


def gap(a, b):
    """Bases between two loci's spans on one contig; 0 when the spans overlap or touch."""
    return max(0, max(a["iv"][0][0], b["iv"][0][0]) - min(a["iv"][-1][1], b["iv"][-1][1]))


def dist_bin(a, b):
    if a["contig"] != b["contig"]:
        return OTHER
    return BIN_NAMES[bisect.bisect_right(BINS, gap(a, b))]


def cover_families(rows, fam_loci):
    """{locus_key: frozenset(family)}: the families with a member covering >= MATCH_COV of the locus's exon bases (the rule
    of `_liftoff.pair_families`)."""
    idx = LF.Index(fam_loci)
    out = {}
    for g in rows:
        n = LF.iv_len(g["iv"])
        out[LF.locus_key(g)] = frozenset(fam_loci[i]["family"] for i in idx.near(g["contig"], g["iv"])
                                         if n and LF.iv_inter(g["iv"], fam_loci[i]["iv"]) / n >= LF.MATCH_COV)
    return out


def decoy_pool(loci, sup, keep=None):
    """The decoy candidates: read-supported Liftoff loci (in place, moved or extra copy; exon union >= SHORT_BP;
    >= SUPPORT_READS reads), on kept contigs."""
    return [r for r in loci if r["cls"] in ELIGIBLE and not r["short"] and r["iv"]
            and sup.get(LF.locus_key(r), 0) >= LF.SUPPORT_READS and (keep is None or keep(r["contig"]))]


class Pool:
    """The decoy pool as arrays, for vectorised candidate search."""

    def __init__(self, rows, fams):
        self.rows = rows
        self.code = {c: i for i, c in enumerate(sorted({r["contig"] for r in rows}))}
        self.scode = {s: i for i, s in enumerate(sorted({r["source_id"] for r in rows}))}
        self.contig = np.array([self.code[r["contig"]] for r in rows], dtype=np.int64)
        self.source = np.array([self.scode[r["source_id"]] for r in rows], dtype=np.int64)
        self.start = np.array([r["iv"][0][0] for r in rows], dtype=np.int64)
        self.end = np.array([r["iv"][-1][1] for r in rows], dtype=np.int64)
        self.covered = np.array([bool(fams[LF.locus_key(r)]) for r in rows], dtype=bool)
        self.fams = [fams[LF.locus_key(r)] for r in rows]

    def candidates(self, a, c, c_covered):
        """-> (pool indices of c's decoys, widening level 0 / 1 / 2), or (empty, 3) when no decoy exists."""
        a0, a1 = a["iv"][0][0], a["iv"][-1][1]
        same = self.contig == self.code.get(a["contig"], -1)
        base = (self.covered == c_covered) & (self.source != self.scode.get(a["source_id"], -1))
        base &= ~(same & (self.start < a1) & (self.end > a0))
        b = dist_bin(a, c)
        if b == OTHER:
            levels = (base & ~same, base & ~same, base)
        else:
            g = np.maximum(0, np.maximum(self.start, a0) - np.minimum(self.end, a1))
            in_bin = same & (np.searchsorted(BINS, g, side="right") == BIN_NAMES.index(b))
            levels = (base & in_bin, base & same, base)
        for level, m in enumerate(levels):
            idx = np.flatnonzero(m)
            if idx.size:
                return idx, level
        return np.empty(0, dtype=np.int64), 3


def p_upper(null, observed):
    """(1 + draws >= observed) / (1 + draws); NaN without an observation."""
    if observed is None or (isinstance(observed, float) and math.isnan(observed)):
        return float("nan")
    return (1 + sum(1 for v in null if v >= observed)) / (1 + len(null))


def run_test(pairs, pool, fams, n_null, seed):
    """The decoy test on one list of pairs (a, c). -> dict with n, k, rate, ci95 (2,000 bootstrap resamples of source genes,
    figure 8's seed), the null's mean and max, p, the widening levels and the coverage of the pairs."""
    levels = collections.Counter()
    kept = []
    for a, c in pairs:
        idx, level = pool.candidates(a, c, bool(fams[LF.locus_key(c)]))
        levels[level] += 1
        if level < 3:
            kept.append((a, c, idx))
    hits = [int(bool(fams[LF.locus_key(a)] & fams[LF.locus_key(c)])) for a, c, _ in kept]
    k, n = sum(hits), len(kept)
    afams = [fams[LF.locus_key(a)] for a, _, _ in kept]
    rng = random.Random(seed)
    null = [sum(1 for fa, (_, _, idx) in zip(afams, kept) if fa & pool.fams[idx[rng.randrange(idx.size)]])
            for _ in range(n_null)]
    lo, hi = LF.boot_ci([a["source_id"] for a, _, _ in kept], hits)
    covered = collections.Counter(int(bool(fams[LF.locus_key(a)])) + int(bool(fams[LF.locus_key(c)]))
                                  for a, c, _ in kept)
    null_mean = sum(null) / n_null if null else float("nan")
    return dict(n=n, k=k, rate=k / n if n else float("nan"), ci95=[lo, hi],
                null_n=n_null, null_seed=seed, null_mean=null_mean, null_rate=null_mean / n if n else float("nan"),
                null_max=max(null) if null else None, p=p_upper(null, k if n else None),
                levels={"0": levels[0], "1": levels[1], "2": levels[2], "none": levels[3]},
                covered={"both": covered[2], "one": covered[1], "none": covered[0]})
