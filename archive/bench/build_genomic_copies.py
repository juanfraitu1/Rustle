#!/usr/bin/env python3
"""Genomic self-alignment pass: find copies the annotation AND protein-clustering miss — pseudogene
copies (no CDS -> no protein) and UNANNOTATED copies. Map all gene-rep transcripts back to the genome
(minimap2 -cx splice, secondary alns), then classify each strong homology hit by what annotated feature
(if any) it lands on. New copies = hits in annotated PSEUDOGENES + hits OUTSIDE all annotation.
RNA overlay: are the new copies transcribed? Run with /home/juanfra/miniforge3/bin/python.
"""
import bisect
import os
import subprocess
from collections import defaultdict

import pysam

OUT = "/home/juanfra/winloci_scratch/gcopies"
FASTA = "/home/juanfra/winloci_scratch/GGO.fasta"
QUERY = "/tmp/gene_reps_gw.fa"
META = "/tmp/gene_reps_gw.meta.tsv"
ANNOT = "/home/juanfra/winloci_scratch/annot_intervals.tsv"
BAM = "/home/juanfra/winloci_scratch/GGO.bam"
MM2 = "/home/juanfra/miniforge3/bin/minimap2"
PAF = f"{OUT}/gene_reps_vs_genome.paf"
QCOV_MIN, ID_MIN = 0.50, 0.80      # a real copy: >=50% of the transcript aligns at >=80% identity
T_EXPR = 5


def main():
    os.makedirs(OUT, exist_ok=True)
    # gene own-locus (to call 'self')
    own = {}
    with open(META) as fh:
        fh.readline()
        for line in fh:
            g, c, s, e, st, ln = line.rstrip("\n").split("\t")
            own[g] = (c, int(s), int(e))
    # annotated intervals per chrom (sorted starts) + biotype/gene
    iv = defaultdict(list)
    with open(ANNOT) as fh:
        fh.readline()
        for line in fh:
            c, s, e, bt, gn = line.rstrip("\n").split("\t")
            iv[c].append((int(s), int(e), bt, gn))
    for c in iv:
        iv[c].sort()
    starts = {c: [x[0] for x in lst] for c, lst in iv.items()}

    def classify(qgene, chrom, ts, te):
        lst = iv.get(chrom)
        if not lst:
            return "unannotated"
        i = bisect.bisect_right(starts[chrom], te)
        hits = []
        for j in range(i - 1, -1, -1):
            s, e, bt, gn = lst[j]
            if e <= ts:
                if ts - s > 3_000_000:
                    break
                continue
            if s < te and e > ts:
                hits.append((bt, gn))
        if any(gn == qgene for _, gn in hits):
            return "self"
        nonpseudo = [bt for bt, gn in hits if "pseudogene" not in bt]
        if nonpseudo:
            return "annotated_gene"
        if hits:
            return "pseudogene_copy"
        return "unannotated"

    if not os.path.exists(PAF):
        print("minimap2 -cx splice (gene reps vs genome, secondary)...")
        with open(PAF, "w") as fh:
            subprocess.run([MM2, "-cx", "splice", "-N", "20", "-p", "0.4", "--secondary=yes",
                            "-t", "5", FASTA, QUERY], stdout=fh, stderr=subprocess.DEVNULL, check=True)
    # parse PAF
    new_copies = []     # (chrom, ts, te, qgene, ident, kind)
    counts = defaultdict(int)
    with open(PAF) as fh:
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) < 12:
                continue
            q, qlen, qs, qe, strand, t, tlen, ts, te, nm, aln, mapq = f[:12]
            qlen, qs, qe, ts, te, nm, aln = int(qlen), int(qs), int(qe), int(ts), int(te), int(nm), int(aln)
            qcov = (qe - qs) / qlen if qlen else 0
            ident = nm / aln if aln else 0
            if qcov < QCOV_MIN or ident < ID_MIN:
                continue
            kind = classify(q, t, ts, te)
            counts[kind] += 1
            if kind in ("pseudogene_copy", "unannotated"):
                new_copies.append((t, ts, te, q, round(ident, 3), kind))

    # dedup new-copy loci (merge overlapping across queries); keep best (highest ident) gene label
    by_chrom = defaultdict(list)
    for c, ts, te, q, idt, kind in new_copies:
        by_chrom[c].append((ts, te, q, idt, kind))
    merged = []
    for c, lst in by_chrom.items():
        lst.sort()
        cs, ce, cq, ci, ck = lst[0]
        for ts, te, q, idt, kind in lst[1:]:
            if ts <= ce:
                ce = max(ce, te)
                if idt > ci:
                    cq, ci, ck = q, idt, kind
            else:
                merged.append((c, cs, ce, cq, ci, ck)); cs, ce, cq, ci, ck = ts, te, q, idt, kind
        merged.append((c, cs, ce, cq, ci, ck))

    # RNA overlay on the new-copy loci
    bam = pysam.AlignmentFile(BAM, "rb")
    n_tx = 0
    rows = []
    for c, s, e, q, idt, kind in merged:
        n = sum(1 for a in bam.fetch(c, s, e)
                if not (a.is_unmapped or a.is_secondary or a.is_supplementary))
        if n >= T_EXPR:
            n_tx += 1
        rows.append((c, s, e, q, idt, kind, n))
    n_pseudo = sum(1 for r in merged if r[5] == "pseudogene_copy")
    n_unann = sum(1 for r in merged if r[5] == "unannotated")

    with open(f"{OUT}/new_copies.tsv", "w") as fh:
        fh.write("chrom\tstart\tend\tbest_gene\tidentity\tkind\tn_reads\ttranscribed\n")
        for c, s, e, q, idt, kind, n in sorted(rows, key=lambda r: -r[6]):
            fh.write(f"{c}\t{s}\t{e}\t{q}\t{idt}\t{kind}\t{n}\t{int(n >= T_EXPR)}\n")

    print(f"PAF hit classes (>= {QCOV_MIN} qcov, >= {ID_MIN} id): " +
          " ".join(f"{k}={v}" for k, v in sorted(counts.items())))
    print(f"distinct NEW-copy loci (merged): {len(merged)}  "
          f"(pseudogene_copy={n_pseudo}, unannotated={n_unann})")
    print(f"  of new-copy loci, transcribed (>= {T_EXPR} reads): {n_tx} "
          f"({100*n_tx/max(len(merged),1):.0f}%)")
    print(f"[wrote {OUT}/new_copies.tsv]")


if __name__ == "__main__":
    main()
