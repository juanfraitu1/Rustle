#!/usr/bin/env python3
"""DNA verification of consensus sequences by k-mer lookup in the individual's whole-genome reads (docs/PREREG_unmapped_rescue_2026-10-08.md Amendment 32).

covered_bases: a base is covered if at least one of the k-mers that contain it has count >= min_count (robust to the k-mers that span exon junctions of a cDNA)."""
import statistics


def covered_bases(counts, k=21, min_count=2):
    n = len(counts)
    if n == 0:
        return 0
    L = n + k - 1
    diff = [0] * (L + 1)
    for i, c in enumerate(counts):
        if c >= min_count:
            diff[i] += 1
            diff[i + k] -= 1
    cov, run = 0, 0
    for j in range(L):
        run += diff[j]
        cov += run > 0
    return cov


def summarize(counts, length, k=21, min_count=2):
    present = [c for c in counts if c >= min_count]
    return dict(kmers=len(counts), frac_ge1=sum(c >= 1 for c in counts) / max(1, len(counts)), frac_ge2=len(present) / max(1, len(counts)),
                median_ge2=statistics.median(present) if present else None, covered=covered_bases(counts, k, min_count) / max(1, length))
