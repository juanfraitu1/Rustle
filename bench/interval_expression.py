#!/usr/bin/env python3
"""Prereg Addendum K1: per-interval read support from aligned blocks.

usage: interval_expression.py <bam> <intervals.tsv> <out.tsv> [contig,contig,...]
<intervals.tsv> has a header with at least `key chrom start end` (0-based half-open). For each interval, counts primary
reads (not unmapped / secondary / supplementary) with >= 1 aligned base (M/=/X block) inside it: `u` = MAPQ >= 1,
`m0` = MAPQ 0. A read counts once per interval however many of its blocks fall inside.
"""
import collections
import csv
import sys

import pysam

BIN = 100_000
bam_path, iv_path, out_path = sys.argv[1], sys.argv[2], sys.argv[3]
only = set(sys.argv[4].split(",")) if len(sys.argv) > 4 else None

ivs = list(csv.DictReader(open(iv_path), delimiter="\t"))
bins = collections.defaultdict(list)  # (chrom, bin) -> interval indices
for k, r in enumerate(ivs):
    s, e = int(r["start"]), int(r["end"])
    for b in range(s // BIN, (e - 1) // BIN + 1):
        bins[(r["chrom"], b)].append(k)
u = [0] * len(ivs)
m0 = [0] * len(ivs)
chroms = sorted({r["chrom"] for r in ivs if only is None or r["chrom"] in only})
bam = pysam.AlignmentFile(bam_path)
n_reads = 0
for chrom in chroms:
    if chrom not in bam.references:
        continue
    for read in bam.fetch(chrom):
        if read.is_unmapped or read.is_secondary or read.is_supplementary:
            continue
        n_reads += 1
        hit = set()
        for bs, be in read.get_blocks():
            for b in range(bs // BIN, (be - 1) // BIN + 1):
                for k in bins.get((chrom, b), ()):
                    if k in hit:
                        continue
                    r = ivs[k]
                    if int(r["start"]) < be and bs < int(r["end"]):
                        hit.add(k)
        for k in hit:
            if read.mapping_quality >= 1:
                u[k] += 1
            else:
                m0[k] += 1
with open(out_path, "w") as fh:
    fh.write("key\tchrom\tstart\tend\tu\tm0\n")
    for k, r in enumerate(ivs):
        fh.write(f"{r['key']}\t{r['chrom']}\t{r['start']}\t{r['end']}\t{u[k]}\t{m0[k]}\n")
print(f"[expr] {n_reads} primary reads over {len(chroms)} contigs -> {len(ivs)} intervals -> {out_path}", file=sys.stderr)
