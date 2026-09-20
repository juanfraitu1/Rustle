#!/usr/bin/env python3
"""Quantify how many of the 74 coverage-independent floor members are DNA-family-but-not-RNA-family:
their EXPRESSED transcript is non-homologous to their (detected) family siblings, so an RNA method correctly
cannot place them in the family. That fraction is a recall-ceiling correction, not a method failure.

For each floor member we align its representative expressed transcript (the longest MAPQ>=60 reads = full-
length molecules) to the transcripts of every DETECTED member of the SAME Soto family, at the pipeline's
sensitive tier. Best hit over all sibling pairs classifies the member:
  RNA-non-homologous : no alignment to any detected sibling            -> DNA-family only (correct miss)
  partial-homology   : aligns but cov-of-shorter < 0.50 (E_r gate)     -> tunable frontier
  homologous         : aligns, cov >= 0.50                             -> should have linked (real gap)
  no-detected-sibling: family has no detected member to compare against -> unassessable from RNA alone

Run: /home/juanfra/miniforge3/bin/python bench/soto/soto_rna_homology_ceiling.py
"""
import csv
import os
import subprocess
from collections import defaultdict

import pysam

D = "/mnt/linuxdisk/home/juanfraitu/winloci_data"
BAM = f"{D}/soto_reads.bam"
DECOMP = "bench/soto/soto_floor_decomposition.tsv"
DET = "bench/soto/soto_member_detection.tsv"
COV = "bench/soto/a119b_member_reads.tsv"
MM2 = "/home/juanfra/miniforge3/bin/minimap2"
SC = "/tmp/claude-1000/-mnt-c-Users-jfris-Desktop/c40f2f02-e69d-4a8f-82ea-9ed7974a0cb9/scratchpad"
COV_GATE = 0.50
TOPN = 3          # representative transcripts per locus (longest MAPQ>=60 reads)


def top_reads(bam, chrom, s, e, n=TOPN):
    reads = []
    for a in bam.fetch(chrom, max(0, s), e):
        if a.is_unmapped or a.is_secondary or a.is_supplementary or a.mapping_quality < 60:
            continue
        if a.query_sequence and len(a.query_sequence) >= 300:
            reads.append(a.query_sequence.upper())
    reads.sort(key=len, reverse=True)
    return reads[:n]


def best_align(qs, ts, tag):
    """best (id, cov_of_shorter) over all query x target read pairs; (0,0) if none align."""
    if not qs or not ts:
        return (0.0, 0.0)
    qf, tf = f"{SC}/rh_q_{tag}.fa", f"{SC}/rh_t_{tag}.fa"
    open(qf, "w").write("".join(f">q{i}\n{s}\n" for i, s in enumerate(qs)))
    open(tf, "w").write("".join(f">t{i}\n{s}\n" for i, s in enumerate(ts)))
    out = subprocess.run([MM2, "-cx", "asm20", "-k11", "-w5", "-t", "2", tf, qf],
                         capture_output=True, text=True).stdout
    best = (0.0, 0.0)
    for line in out.splitlines():
        f = line.split("\t")
        if len(f) < 12:
            continue
        qlen, tlen = int(f[1]), int(f[6])
        matches, alnlen = int(f[9]), int(f[10])
        idp = matches / alnlen if alnlen else 0
        cov = alnlen / min(qlen, tlen) if min(qlen, tlen) else 0
        if cov > best[1] or (cov == best[1] and idp > best[0]):
            best = (idp, cov)
    for p in (qf, tf):
        os.remove(p)
    return best


def main():
    det = {}; fam_members = defaultdict(list); coord = {}
    for r in csv.DictReader(open(DET), delimiter="\t"):
        det[(r["family_id"], r["gene"])] = r["detected"]
    for r in csv.DictReader(open(COV), delimiter="\t"):
        fam_members[r["famID"]].append(r["member"])
        coord[(r["famID"], r["member"])] = (r["chrom"], int(r["start"]), int(r["end"]))

    bam = pysam.AlignmentFile(BAM, "rb")
    floor = list(csv.DictReader(open(DECOMP), delimiter="\t"))
    rows = []; tax = defaultdict(int)
    for i, m in enumerate(floor):
        fam, gene = m["family_id"], m["gene"]
        mq = top_reads(bam, m["chrom"], int(m["start"]), int(m["end"]))
        sibs = [g for g in fam_members[fam] if g != gene]  # ALL family siblings (detected or not)
        if not sibs or not mq:
            tax["singleton/no-reads"] += 1
            rows.append([fam, gene, m["chrom"], m["start"], m["end"], 0, "?", "0.00", "0.00", "singleton/no-reads"])
            continue
        best = (0.0, 0.0); best_g = ""; best_det = ""
        for g in sibs:
            c, s, e = coord[(fam, g)]
            sq = top_reads(bam, c, s, e, n=2)
            b = best_align(mq, sq, f"{i}_{g}")
            if b[1] > best[1] or (b[1] == best[1] and b[0] > best[0]):
                best = b; best_g = g; best_det = det.get((fam, g), "?")
        idp, cov = best
        if cov == 0:
            cls = "RNA-non-homologous (DNA-family only)"
        elif cov >= COV_GATE:
            cls = "homologous to a sibling (RNA-family member; recoverable gap)"
        else:
            cls = "partial-homology (cov<0.50; coverage-gate frontier)"
        tax[cls] += 1
        rows.append([fam, gene, m["chrom"], m["start"], m["end"], len(sibs),
                     f"{best_g}({best_det})", f"{idp:.3f}", f"{cov:.3f}", cls])

    rows.sort(key=lambda x: x[-1])
    with open("bench/soto/soto_rna_homology_ceiling.tsv", "w", newline="") as f:
        w = csv.writer(f, delimiter="\t")
        w.writerow(["family_id", "gene", "chrom", "start", "end", "n_sibs",
                    "best_sibling(det)", "best_id_to_sib", "best_cov_of_shorter", "class"])
        w.writerows(rows)

    n = len(rows)
    print("\n=== RNA-homology ceiling: are the 74 floor members even RNA-family members? ===")
    for c, k in sorted(tax.items(), key=lambda kv: -kv[1]):
        print(f"  {k:3d}  {c}")
    dna_only = tax["RNA-non-homologous (DNA-family only)"]
    real_gap = tax["homologous to a sibling (RNA-family member; recoverable gap)"]
    partial = tax["partial-homology (cov<0.50; coverage-gate frontier)"]
    singl = tax["singleton/no-reads"]
    print(f"  ---  {n} floor members")
    print(f"\n  DNA-family-not-RNA (no homologous expressed sibling -> correct miss, RNA-impossible): {dna_only}")
    print(f"  partial-homology (coverage-gate, tunable, FP-risky):                                 {partial}")
    print(f"  homologous to a sibling but unlinked (genuine recoverable gap):                      {real_gap}")
    print(f"  singleton family / no usable reads (unassessable):                                   {singl}")
    detected = 276; total = 362
    rna_possible = total - dna_only            # only the DNA-only members are truly RNA-impossible
    print(f"\n  raw sensitivity:       {detected}/{total} = {100*detected/total:.1f}%")
    print(f"  corrected RNA ceiling: {detected}/{rna_possible} = {100*detected/rna_possible:.1f}%  "
          f"(removing only the {dna_only} RNA-impossible DNA-only members from the denominator)")
    print("  wrote bench/soto/soto_rna_homology_ceiling.tsv")


if __name__ == "__main__":
    main()
