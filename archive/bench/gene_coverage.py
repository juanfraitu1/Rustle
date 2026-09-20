#!/usr/bin/env python3
"""Count primary reads per gene in the real GGO IsoSeq BAM -> defines the under-covered ('lacking')
set for the ideal-coverage top-up dataset, and sizes it."""
import pysam

META = "/tmp/gene_reps_gw.meta.tsv"
BAM = "/home/juanfra/winloci_scratch/GGO.bam"
OUT = "/home/juanfra/winloci_scratch/gene_cov.tsv"
TARGETS = (20, 30, 40, 50)

bam = pysam.AlignmentFile(BAM, "rb")
rows = []
with open(META) as fh:
    fh.readline()
    for line in fh:
        g, c, s, e, strand, ln = line.rstrip("\n").split("\t")
        s, e, ln = int(s), int(e), int(ln)
        n = 0
        for a in bam.fetch(c, s, e):
            if not (a.is_unmapped or a.is_secondary or a.is_supplementary):
                n += 1
        rows.append((g, c, s, e, ln, n))

with open(OUT, "w") as fh:
    fh.write("gene\tchrom\tstart\tend\ttx_len\tn_reads\n")
    for r in rows:
        fh.write("\t".join(map(str, r)) + "\n")

import numpy as np
nr = np.array([r[5] for r in rows])
lens = np.array([r[4] for r in rows])
print(f"genes={len(rows)}  reads: median={int(np.median(nr))} mean={nr.mean():.1f} "
      f"zero={int((nr==0).sum())}")
for t in TARGETS:
    under = nr < t
    topup = np.clip(t - nr, 0, None)[under].sum()
    # rough FASTQ size: reads * avg tx_len * ~2.1 bytes/base (seq+qual)
    avglen = lens[under].mean() if under.any() else 0
    gb = topup * avglen * 2.1 / 1e9
    print(f"target {t}x: genes under={int(under.sum())}  reads to simulate={int(topup)}  "
          f"~FASTQ {gb:.1f} GB (avg tx_len {avglen:.0f})")
print(f"[wrote {OUT}]")
