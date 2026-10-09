#!/usr/bin/env python3
"""Amendment 43b: variant k-mers of a consensus (absent from the reference segments it aligns to) and their DNA label from the individual's whole-genome reads.

Pure (tested in test_varkmers.py): kmer_codes (canonical 21-mers, the kmerhit encoding A0 C1 G2 T3), query_junctions, not_across, label."""
import re

import numpy as np

K = 21
LUT = np.full(256, 255, dtype=np.uint8)
for i, c in enumerate("ACGT"):
    LUT[ord(c)] = i
    LUT[ord(c.lower())] = i
CIG = re.compile(r"(\d+)([MIDNSHP=X])")


def kmer_codes(seq, k=K):
    """-> (start positions, canonical codes) of the k-mers without N"""
    a = LUT[np.frombuffer(seq.encode(), dtype=np.uint8)]
    if len(a) < k:
        return np.zeros(0, dtype=np.int64), np.zeros(0, dtype=np.uint64)
    win = np.lib.stride_tricks.sliding_window_view(a, k)
    ok = (win != 255).all(axis=1)
    w = win.astype(np.uint64)
    sh = (2 * np.arange(k - 1, -1, -1)).astype(np.uint64)
    fw = np.bitwise_or.reduce(w << sh, axis=1)
    rv = np.bitwise_or.reduce((np.uint64(3) - w) << (2 * np.arange(k, dtype=np.uint64)), axis=1)
    can = np.minimum(fw, rv)
    pos = np.nonzero(ok)[0]
    return pos.astype(np.int64), can[ok]


def query_junctions(cigar):
    """query positions of the introns (N) of an alignment"""
    q, out = 0, []
    for n, op in CIG.findall(cigar):
        n = int(n)
        if op in "MI=XS":
            q += n
        elif op == "N":
            out.append(q)
    return out


def not_across(pos, junctions, k=K):
    """boolean mask of the k-mers (start positions) that do not contain a junction strictly inside them"""
    keep = np.ones(len(pos), dtype=bool)
    for j in junctions:
        keep &= ~((pos < j) & (pos + k > j))
    return keep


def label(counts, min_n=10, present=2, hi=0.5, lo=0.1):
    """DNA-SUPPORTED: >= min_n variant k-mers and >= hi of them seen >= present times; RNA-ONLY: <= lo of them; else UNDECIDED"""
    n = len(counts)
    if n < min_n:
        return "UNDECIDED"
    f = sum(c >= present for c in counts) / n
    return "DNA-SUPPORTED" if f >= hi else "RNA-ONLY" if f <= lo else "UNDECIDED"
