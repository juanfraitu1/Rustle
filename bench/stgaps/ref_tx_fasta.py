#!/usr/bin/env python3
"""Spliced sequences of every multi-exon reference transcript (FAMILY_LAYER_SPEC Part A).
usage: ref_tx_fasta.py REF.gtf[.gz] GENOME.fa OUT.fa OUT.tsv   (tsv: tid gene chrom strand start end)"""
import gzip, re, sys
from collections import defaultdict
import pysam
gtf, fa_path, out_fa, out_tsv = sys.argv[1:5]
RT, RG = re.compile(r'transcript_id "([^"]+)"'), re.compile(r'gene_id "([^"]*)"')
ex, meta = defaultdict(list), {}
for line in (gzip.open(gtf, "rt") if gtf.endswith(".gz") else open(gtf)):
    f = line.split("\t")
    if len(f) < 9 or f[2] != "exon":
        continue
    t = RT.search(f[8]).group(1)
    ex[t].append((int(f[3]) - 1, int(f[4])))
    meta.setdefault(t, (RG.search(f[8]).group(1) or t, f[0], f[6]))  # empty gene_id: the transcript is its own gene
fa = pysam.FastaFile(fa_path)
comp = str.maketrans("ACGTNacgtn", "TGCANtgcan")
with open(out_fa, "w") as fo, open(out_tsv, "w") as ft:
    for t, e in ex.items():
        if len(e) < 2:
            continue
        g, c, s = meta[t]
        if c not in fa.references:
            continue
        e.sort()
        seq = "".join(fa.fetch(c, a, b) for a, b in e).upper()
        if s == "-":
            seq = seq[::-1].translate(comp)
        fo.write(f">{t}\n{seq}\n")
        ft.write(f"{t}\t{g}\t{c}\t{s}\t{e[0][0]}\t{e[-1][1]}\n")
