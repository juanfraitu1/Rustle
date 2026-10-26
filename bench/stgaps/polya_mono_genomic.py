#!/usr/bin/env python3
"""Variant of polya_mono.py: a read is primed when the GENOME 20 bp past its aligned end is >= 60% A or holds A6
(r1064/r1135 instrument, copy_assign TES_PRIME_*), on + (after reference_end) or, mirrored, on - (before
reference_start, T-rich). Strand of an unspliced read is unknown, so either side counts.
usage: polya_mono_genomic.py polya_mono.tsv FASTA BAM"""
import sys
from collections import Counter
import pysam

tsv, fa_path, bam_path = sys.argv[1:4]
fa, bam = pysam.FastaFile(fa_path), pysam.AlignmentFile(bam_path, "rb")


def rich(s, b):
    return bool(s) and (s.count(b) >= 0.6 * len(s) or b * 6 in s)


rows = [l.rstrip("\n").split("\t") for l in open(tsv)][1:]
out = []
for t, ctg, s, e, nr, n, p, c in rows:
    s, e = int(s), int(e)
    tot = pr = 0
    for r in bam.fetch(ctg, s, e):
        if r.flag & 2820 or r.cigartuples is None or any(op == 3 for op, _ in r.cigartuples):
            continue
        if min(e, r.reference_end) - max(s, r.reference_start) < 0.5 * (e - s):
            continue
        tot += 1
        dn = fa.fetch(ctg, r.reference_end, r.reference_end + 20).upper()
        up = fa.fetch(ctg, max(0, r.reference_start - 20), r.reference_start).upper()
        pr += rich(dn, "A") or rich(up, "T")
    out.append((t, tot, pr, c))
print(f"single-exon {len(out)}")
for thr in (0.5, 0.8, 1.0):
    hit = [r for r in out if r[1] and r[2] / r[1] >= thr]
    print(f"genomic-primed share >= {thr}: {len(hit)}; classes {dict(Counter(r[3] for r in hit))}")
