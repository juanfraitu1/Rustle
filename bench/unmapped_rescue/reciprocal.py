#!/usr/bin/env python3
"""Reciprocal best match for an O3 COPY call (docs/PREREG_unmapped_rescue_2026-10-08.md Amendment 27).

A consensus that is >= 0.999 on the truth haplotype and more than the allele cutoff from every reference locus may still be a DIVERGED ORTHOLOG (a locus both haplotypes
carry, 1 to 7% apart). It is an absent copy only if the reference locus it matches belongs to another truth locus: map the reference locus' transcript back to the truth
haplotype; if it lands on the truth hit of the consensus, the two are one-to-one orthologs."""
import re

CIG = re.compile(r"(\d+)([MIDNSHP=X])")
_RC = str.maketrans("ACGTacgt", "TGCAtgca")


def target_segments(cigar, tstart):
    """reference intervals (0-based, half-open) covered by an alignment: matches, mismatches and bases deleted from the query are kept; N (an intron) and I are not;
    adjacent intervals are merged"""
    segs, pos = [], tstart
    for n, op in CIG.findall(cigar):
        n = int(n)
        if op in "M=XD":
            if segs and segs[-1][1] == pos:
                segs[-1] = (segs[-1][0], pos + n)
            else:
                segs.append((pos, pos + n))
            pos += n
        elif op == "N":
            pos += n
    return segs


def same_locus(a, b, min_frac=0.5):
    """a, b = (chrom, start, end) or None; same chromosome and an overlap of at least min_frac of the shorter span"""
    if a is None or b is None or a[0] != b[0]:
        return False
    ov = min(a[2], b[2]) - max(a[1], b[1])
    return ov >= min_frac * min(a[2] - a[1], b[2] - b[1]) and ov > 0


def orient(seq, strand):
    return seq if strand == "+" else seq.translate(_RC)[::-1]
