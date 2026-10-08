#!/usr/bin/env python3
"""The 5' artifact read from the library (docs/PREREG_unmapped_rescue_2026-10-08.md Amendment 15): the 5' soft clips of cleanly aligned reads of surviving copies.

five_prime_clip: the 5' end of the ORIGINAL read (left soft clip of a forward read; right soft clip of a reverse read, reverse-complemented).
signature: counts of pure single-base clips of 1 to 3 bases among clean primary reads; gate: is G over-represented (binomial, one-sided)."""
import re

import polish as P

CIG = re.compile(r"(\d+)([MIDNSHP=X])")
RC = str.maketrans("ACGTacgt", "TGCAtgca")
MAX_CLIP = 3
ALPHA = 1e-6
MIN_CLIPS = 20
MAX_ALL = 8
FLOOR = 1e-4


def five_prime_clip(flag, cigar, seq):
    ops = CIG.findall(cigar)
    ops = [(int(n), o) for n, o in ops if o != "H"]
    if not ops:
        return 0, ""
    if not flag & 16:
        n, o = ops[0]
        return (n, seq[:n]) if o == "S" else (0, "")
    n, o = ops[-1]
    return (n, seq[len(seq) - n:].translate(RC)[::-1]) if o == "S" else (0, "")


def signature(sam_lines, keep=lambda name: True, max_de=0.0096, min_mapq=10, max_len=MAX_CLIP):
    pure = {b: 0 for b in "ACGT"}
    g_len, g_all, reads, no_clip, other = {}, {}, 0, 0, 0
    for ln in sam_lines:
        if ln[0] == "@":
            continue
        f = ln.rstrip("\n").split("\t")
        if int(f[1]) & 2308 or f[2] == "*" or int(f[4]) < min_mapq or not keep(f[0]):
            continue
        de = next((float(t[5:]) for t in f[11:] if t.startswith("de:f:")), 1.0)
        if de > max_de:
            continue
        reads += 1
        n, clip = five_prime_clip(int(f[1]), f[5], f[9])
        if 0 < n <= MAX_ALL and len(set(clip)) == 1 and clip[0] == "G":
            g_all[n] = g_all.get(n, 0) + 1
        if n == 0:
            no_clip += 1
        elif n <= max_len and len(set(clip)) == 1 and clip[0] in pure:
            pure[clip[0]] += 1
            if clip[0] == "G":
                g_len[n] = g_len.get(n, 0) + 1
        else:
            other += 1
    return dict(reads=reads, no_clip=no_clip, pure=pure, g_lengths=g_len, g_lengths_all=g_all, other_clip=other)


def gate(sig, alpha=ALPHA, min_clips=MIN_CLIPS):
    """-> (artifact present, p). p = P(Binomial(n, 0.5) >= G count), n = pure-clip count; false if fewer than min_clips pure clips"""
    n = sum(sig["pure"].values())
    if n < min_clips:
        return False, 1.0
    p = P.binom_tail(sig["pure"]["G"], n, 0.5)
    return p < alpha, p


def artifact_distribution(sig, max_len=MAX_ALL, floor=FLOOR):
    """a(l), l = 0..max_len: the share of clean survivor reads whose 5' clip is a pure G run of l bases (l = 0: no clip); a floor for lengths never seen (Amendment 16)"""
    n = max(1, sig["reads"])
    a = {0: max(floor, sig["no_clip"] / n)}
    for l in range(1, max_len + 1):
        a[l] = max(floor, sig["g_lengths_all"].get(l, sig["g_lengths_all"].get(str(l), 0)) / n)
    return a
