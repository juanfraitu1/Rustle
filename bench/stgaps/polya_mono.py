#!/usr/bin/env python3
"""StringTie's unspliced-read poly-A filter (rlink.cpp:486 genome_polyA_terminated, applied at :605 to unspliced
long reads only), measured on the single-exon transcripts of a polished GTF.

usage: polya_mono.py POLISHED.gtf TMAP BAM OUT.tsv
Per single-exon transcript: unspliced primary reads (flag & 2820 == 0, no N) overlapping >= 50% of the transcript,
how many are genome-polyA-terminated, and the gffcompare class. Prints the class mix of transcripts whose
non-primed reads fall below the transcript's own `reads` * (1 - primed share) floor proxy (see summary lines).
"""
import re
import sys
from collections import Counter

import pysam

gtf, tmap, bam_path, out = sys.argv[1:5]
RE_TID = re.compile(r'transcript_id "([^"]+)"')
RE_READS = re.compile(r'reads "(\d+)"')


def polya_terminated(seq, clip_l, clip_r):
    """Exact port of genome_polyA_terminated: A-run at the aligned 3' end or T-run at the aligned 5' start."""
    n = len(seq) - clip_r
    a = e = 0
    p = n - 1
    while p >= 0 and e < 25:
        e += 1
        a += seq[p] == "A"
        if e - a > 5:
            break
        if a >= 10 and a / e >= 0.8:
            return True
        p -= 1
    a = e = 0
    p = clip_l
    while p < n and e < 25:
        e += 1
        a += seq[p] == "T"
        if e - a > 5:
            break
        if a >= 10 and a / e >= 0.8:
            return True
        p += 1
    return False


cls = {}
for line in open(tmap):
    f = line.rstrip("\n").split("\t")
    if f[0] != "ref_gene_id" and f[5] == "1":
        cls[f[4]] = f[2]
mono = {}
for line in open(gtf):
    f = line.split("\t")
    if len(f) >= 9 and f[2] == "transcript":
        t = RE_TID.search(f[8]).group(1)
        if t in cls:
            r = RE_READS.search(f[8])
            mono[t] = (f[0], int(f[3]) - 1, int(f[4]), int(r.group(1)) if r else 0)

bam = pysam.AlignmentFile(bam_path, "rb")
rows = []
for t, (ctg, s, e, nr) in sorted(mono.items(), key=lambda x: x[1]):
    n = p = 0
    for r in bam.fetch(ctg, s, e):
        if r.flag & 2820 or r.cigartuples is None or any(op == 3 for op, _ in r.cigartuples):
            continue
        ov = min(e, r.reference_end) - max(s, r.reference_start)
        if ov < 0.5 * (e - s):
            continue
        ct = r.cigartuples
        cl = ct[0][1] if ct[0][0] == 4 else 0
        cr = ct[-1][1] if ct[-1][0] == 4 else 0
        n += 1
        p += polya_terminated(r.query_sequence, cl, cr)
    rows.append((t, ctg, s, e, nr, n, p, cls[t]))

with open(out, "w") as fo:
    fo.write("tid\tcontig\tstart\tend\treads_attr\tunspliced\tprimed\tclass\n")
    for r in rows:
        fo.write("\t".join(map(str, r)) + "\n")
print(f"single-exon transcripts {len(rows)}; classes {dict(Counter(r[7] for r in rows))}")
for thr in (0.5, 0.8, 1.0):
    hit = [r for r in rows if r[5] and r[6] / r[5] >= thr]
    print(f"primed share >= {thr}: {len(hit)} transcripts; classes {dict(Counter(r[7] for r in hit))}")
