#!/usr/bin/env python3
"""Shared full-length HiFi transcript-read simulator (deterministic). IsoSeq reads are full-length,
high-accuracy transcripts, so each read = one transcript sequence + HiFi errors (substitutions +
rare short indels), with optional 5'/3' degradation truncation. A KNOWN error rate is what the
identifiability theorem needs (PSVs must clear the error floor). No wall-clock RNG.
"""
import random

BASES = "ACGT"


def simulate_reads(seq, n, err=0.003, indel=0.001, seed=0, trunc_frac=0.0):
    """n reads from `seq` with per-base substitution rate `err` and indel rate `indel` each.
    trunc_frac: max fraction trimmed from a random end (models IsoSeq 5'/3' degradation).
    Returns list of (read_seq, qual_str). Deterministic given (seq, n, seed)."""
    out = []
    for i in range(n):
        rng = random.Random((seed * 1_000_003) ^ (i * 2_654_435_761) ^ len(seq))
        s = seq
        if trunc_frac > 0:
            t = int(rng.random() * trunc_frac * len(s))
            if t:
                if rng.random() < 0.5:
                    s = s[t:]
                else:
                    s = s[:len(s) - t]
        buf = []
        for ch in s:
            r = rng.random()
            if r < err:
                buf.append(rng.choice([b for b in BASES if b != ch]))   # substitution
            elif r < err + indel:
                continue                                                # deletion
            elif r < err + 2 * indel:
                buf.append(ch); buf.append(rng.choice(BASES))           # insertion
            else:
                buf.append(ch)
        rd = "".join(buf)
        out.append((rd, "~" * len(rd)))   # '~' = Q93 placeholder (HiFi-grade)
    return out


def write_fastq(fh, name, read_qual):
    rd, q = read_qual
    fh.write(f"@{name}\n{rd}\n+\n{q}\n")
