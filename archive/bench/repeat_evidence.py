#!/usr/bin/env python3
"""Evidence-agnostic families (spec 2026-09-14 §1a): repeat evidence R1-R4 as one BED format
(chrom, start, end, class, source); 0-based half-open. R1 RepeatMasker .out (class), R2 lowercase runs of a soft-masked
assembly, R3 WindowMasker + DustMasker intervals, R4 meryl high-copy runs (dup_evidence.py)."""
import bisect
import collections
import gzip
import re

import pysam


def _open(path):
    path = str(path)
    return gzip.open(path, "rt") if path.endswith(".gz") else open(path)


def parse_rmsk(path, contigs=None):
    out = []
    for line in _open(path):
        f = line.split()
        if len(f) < 11 or not f[0].isdigit():
            continue
        if contigs is not None and f[4] not in contigs:
            continue
        out.append((f[4], int(f[5]) - 1, int(f[6]), f[10]))
    return out


def lowercase_runs(fasta, contigs=None):
    g = pysam.FastaFile(str(fasta))
    out = []
    for c in g.references:
        if contigs is not None and c not in contigs:
            continue
        for m in re.finditer(r"[a-z]+", g.fetch(c)):
            out.append((c, m.start(), m.end(), "."))
    return out


def parse_masker_intervals(path):
    out, chrom = [], None
    for line in open(path):
        line = line.strip()
        if not line:
            continue
        if line.startswith(">"):
            chrom = line[1:].split()[0]
            continue
        a, b = line.split(" - ")
        out.append((chrom, int(a), int(b) + 1, "."))
    return out


def merge(ivs):
    by = collections.defaultdict(list)
    for c, s, e, _ in ivs:
        by[c].append((s, e))
    out = {}
    for c, v in by.items():
        m = []
        for s, e in sorted(v):
            if m and s <= m[-1][1]:
                m[-1] = (m[-1][0], max(m[-1][1], e))
            else:
                m.append((s, e))
        out[c] = m
    return out


def masked_bases(merged, chrom, s, e):
    v = merged.get(chrom)
    if not v or e <= s:
        return 0
    lo, hi = 0, len(v)
    while lo < hi:  # first interval ending after s
        mid = (lo + hi) // 2
        if v[mid][1] <= s:
            lo = mid + 1
        else:
            hi = mid
    tot = 0
    for a, b in v[lo:]:
        if a >= e:
            break
        tot += max(0, min(b, e) - max(a, s))
    return tot


def write_bed(ivs, path, source):
    with open(path, "w") as fh:
        for c, s, e, cls in ivs:
            fh.write(f"{c}\t{s}\t{e}\t{cls}\t{source}\n")


def read_bed(path):
    return [(f[0], int(f[1]), int(f[2]), f[3]) for f in (l.rstrip("\n").split("\t") for l in open(path)) if len(f) >= 4]
