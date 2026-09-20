#!/usr/bin/env python3
"""Pass 1 the READ-COHERENCE way, GENOME-WIDE, annotation-free: one pass over the BAM, group PRIMARY
reads by their exact intron chain -> de-novo transcript SKELETONS (direct transfrags, no annotation, no
flow). This is the annotation-free Pass-1 input that the family + copy-assignment machinery (Pass 2)
would consume instead of RefSeq gene reps. Reports the de-novo skeleton set genome-wide.
Run with /home/juanfra/miniforge3/bin/python.
"""
import pysam
from collections import defaultdict

BAM = "/home/juanfra/winloci_scratch/GGO_mm.bam"  # the COMPLETE -N50 -p0.1 multimapping BAM (the old
                                                   # GGO.bam was the undercount; recompute against this)
OUT = "/home/juanfra/winloci_scratch/denovo_skeletons.tsv"
MIN_READS = 2


def main():
    bam = pysam.AlignmentFile(BAM, "rb")
    groups = defaultdict(lambda: [0, 10**12, 0])   # key -> [nreads, min_start, max_end]
    n_primary = 0
    for aln in bam.fetch():
        if aln.is_unmapped or aln.is_secondary or aln.is_supplementary:
            continue
        n_primary += 1
        rpos = aln.reference_start
        introns = []
        for op, length in aln.cigartuples:
            if op == 3:
                introns.append((rpos, rpos + length)); rpos += length
            elif op in (0, 2, 7, 8):
                rpos += length
        key = (aln.reference_name, tuple(introns))   # intron chain = the read-coherence skeleton key
        g = groups[key]
        g[0] += 1
        g[1] = min(g[1], aln.reference_start)
        g[2] = max(g[2], aln.reference_end)

    multi = {k: v for k, v in groups.items() if v[0] >= MIN_READS}
    spliced = {k: v for k, v in multi.items() if len(k[1]) >= 1}   # multi-exon skeletons
    import numpy as np
    sizes = np.array([v[0] for v in multi.values()])
    nex = np.array([len(k[1]) + 1 for k in multi])  # exon count = introns + 1
    with open(OUT, "w") as fh:
        fh.write("chrom\tstart\tend\tn_exon\tn_reads\tintrons\n")
        for (c, intr), (nr, s, e) in sorted(multi.items(), key=lambda kv: -kv[1][0]):
            fh.write(f"{c}\t{s}\t{e}\t{len(intr)+1}\t{nr}\t{';'.join(f'{a}-{b}' for a,b in intr)}\n")

    print(f"primary reads scanned: {n_primary:,}")
    print(f"de-novo read-coherence skeletons (>= {MIN_READS} reads): {len(multi):,} "
          f"(multi-exon: {len(spliced):,}; single-exon: {len(multi)-len(spliced):,})")
    print(f"  reads/skeleton: median={int(np.median(sizes))} mean={sizes.mean():.1f} max={sizes.max()}")
    print(f"  exon count: median={int(np.median(nex))} max={nex.max()}")
    print(f"  (compare: RefSeq annotated transcripts ~84,709; annotation gene reps 22,983)")
    print(f"[wrote {OUT}]  -> this is the annotation-free Pass-1 skeleton set; Pass 2 (minimizer-LSH +")
    print("   POA family graph + de-novo PSV split) runs the SAME machinery on these instead of gene reps.")


if __name__ == "__main__":
    main()
