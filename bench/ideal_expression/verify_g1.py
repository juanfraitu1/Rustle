#!/usr/bin/env python3
"""G1 of docs/PREREG_ideal_expression_2026-10-06.md, written separately from the generator: every kept intron of every simulated transcript must be >= 50 bp and canonical
(GT-AG, GC-AG, AT-AC on the transcript's strand), re-read from the genome with this script's own motif table; no simulated molecule may exceed 30 kb or keep an intron > 200 kb.

    verify_g1.py TRUTHPREFIX
"""
import csv
import sys

import pysam

GENOME = "/mnt/linuxdisk/home/juanfraitu/winloci_data/chm13v2.0.fa"
COMP = str.maketrans("ACGT", "TGCA")
OK = {"GTAG", "GCAG", "ATAC"}

fa = pysam.FastaFile(GENOME)
n_tx = n_int = bad = 0
for r in csv.DictReader(open(sys.argv[1] + ".transcripts.tsv"), delimiter="\t"):
    if not int(r["simulated"]):
        continue
    n_tx += 1
    if int(r["spliced_len"]) > 30_000 or int(r["spliced_len"]) < 120:
        bad += 1; print("length", r["transcript"], r["spliced_len"])
    for x in (r["chain"].split(",") if r["chain"] else []):
        d, a = map(int, x.split("-"))
        n_int += 1
        donor, acc = fa.fetch(r["chrom"], d, d + 2).upper(), fa.fetch(r["chrom"], a - 2, a).upper()
        # `donor` = the 2 bp at the genomic start of the intron, `acc` = the 2 bp at its genomic end. Plus strand: donor+acc. Minus strand: the transcript reads the intron
        # backwards, so its donor is revcomp(acc) and its acceptor is revcomp(donor).
        motif = donor + acc if r["strand"] == "+" else acc[::-1].translate(COMP) + donor[::-1].translate(COMP)
        if motif not in OK or a - d < 50 or a - d > 200_000:
            bad += 1; print("intron", r["transcript"], d, a, r["strand"], motif)
print(f"G1: {n_tx} simulated transcripts, {n_int} kept introns, {bad} violations")
sys.exit(1 if bad else 0)
