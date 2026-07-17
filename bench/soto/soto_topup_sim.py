#!/usr/bin/env python3
"""Coverage top-up for the Soto benchmark: give every UNDETECTED member IDEAL coverage of its OWN observed
transcript, then re-run detection — so families lost to coverage recover while the K=0 floor stays missing.

Method (real-read replication, not annotation-based): for each undetected member we fetch its own PRIMARY
reads from the scoped Soto BAM and resample them with HiFi sequencing noise up to TARGET coverage. This is
faithful (uses the member's actual spliced transcript, not a guessed annotation — 62/86 undetected members
are AC*/AL* clones with no gene-level exon model) AND makes K=0 self-enforcing: a member whose expressed
sequence is identical to a sibling gets identical top-up reads that still map MAPQ-0 / primary-at-sibling, so
it will still not seed its own locus. Only members with real coverage give a template — a member with 0
primary reads (silent / secondary-sink) has nothing to resample and is reported, not fabricated.

Disk-light: simulate -> map -> delete FASTQ/SAM; keep BAMs. Big outputs to /mnt/linuxdisk (not the VHDX).
Run: /home/juanfra/miniforge3/bin/python bench/soto/soto_topup_sim.py
"""
import os
import subprocess
import sys

import pysam

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
from sim_reads import simulate_reads, write_fastq

D = "/mnt/linuxdisk/home/juanfraitu/winloci_data"
REAL_BAM = f"{D}/soto_reads.bam"
FASTA = f"{D}/chm13v2.0.fa"
OUT = f"{D}/topup"
DET = "bench/soto/soto_member_detection.tsv"
COV = "bench/soto/a119b_member_reads.tsv"
MM2 = "/home/juanfra/miniforge3/bin/minimap2"
SAM = "/home/juanfra/miniforge3/bin/samtools"
TARGET = 40          # top every undetected member up to this many reads (build_topup precedent = "ideal coverage")
MAX_TEMPLATES = 60   # cap templates fetched per member
ERR, INDEL, TRUNC = 0.003, 0.001, 0.10


def load_undetected():
    """[(fam, gene, chrom, start, end, primary_reads), ...] for members detected == N."""
    det = {}
    for line in open(DET):
        p = line.rstrip("\n").split("\t")
        if p[0] == "family_id" or len(p) < 6:
            continue
        det[(p[0], p[1])] = p[5]
    out = []
    for line in open(COV):
        p = line.rstrip("\n").split("\t")
        if p[0] == "chrom" or len(p) < 6:
            continue
        chrom, start, end, gene, fam, nr = p[0], int(p[1]), int(p[2]), p[3], p[4], int(p[5])
        if det.get((fam, gene)) == "N":
            out.append((fam, gene, chrom, start, end, nr))
    return out


def fetch_primary_templates(bam, chrom, start, end):
    """Full query sequences of PRIMARY reads overlapping the member span (the member's own transcripts)."""
    seqs = []
    for aln in bam.fetch(chrom, max(0, start), end):
        if aln.is_unmapped or aln.is_secondary or aln.is_supplementary:
            continue
        if aln.query_sequence and len(aln.query_sequence) >= 100:
            seqs.append(aln.query_sequence.upper())
            if len(seqs) >= MAX_TEMPLATES:
                break
    return seqs


def main():
    os.makedirs(OUT, exist_ok=True)
    undet = load_undetected()
    bam = pysam.AlignmentFile(REAL_BAM, "rb")
    fq = f"{OUT}/topup.fq"
    truth = open(f"{OUT}/topup_truth.tsv", "w")
    truth.write("family_id\tgene\tchrom\tstart\tend\tprimary_reads\tn_templates\tn_added\tstatus\n")
    n_reads = n_topped = n_notemplate = 0
    with open(fq, "w") as fh:
        for (fam, gene, chrom, start, end, nr) in undet:
            add = max(0, TARGET - nr)
            if add == 0:
                truth.write(f"{fam}\t{gene}\t{chrom}\t{start}\t{end}\t{nr}\t0\t0\talready>=TARGET\n")
                continue
            templates = fetch_primary_templates(bam, chrom, start, end)
            if not templates:
                truth.write(f"{fam}\t{gene}\t{chrom}\t{start}\t{end}\t{nr}\t0\t0\tno-primary-template\n")
                n_notemplate += 1
                continue
            seed0 = abs(hash((fam, gene))) % 2_000_000
            for i in range(add):
                tmpl = templates[i % len(templates)]
                rd = simulate_reads(tmpl, 1, err=ERR, indel=INDEL, seed=seed0 + i, trunc_frac=TRUNC)[0]
                write_fastq(fh, f"SIMTOPUP|{fam}|{gene}|{i}", rd)
                n_reads += 1
            truth.write(f"{fam}\t{gene}\t{chrom}\t{start}\t{end}\t{nr}\t{len(templates)}\t{add}\ttopped-up\n")
            n_topped += 1
    truth.close()
    print(f"topped up {n_topped} members (+{n_reads} reads); {n_notemplate} had no primary template (silent); "
          f"{len(undet)-n_topped-n_notemplate} already >= {TARGET} reads")

    print("mapping top-up reads (minimap2 -ax splice:hq -N 50 -p 0.1 -> CHM13, same flags as real Soto)...")
    sam = f"{OUT}/topup.sam"
    with open(sam, "w") as sf:
        subprocess.run([MM2, "-ax", "splice:hq", "-N", "50", "-p", "0.1", "-t", "6", FASTA, fq],
                       stdout=sf, stderr=subprocess.DEVNULL, check=True)
    synth = f"{OUT}/topup.bam"
    subprocess.run([SAM, "sort", "-@", "4", "-o", synth, sam], check=True, stderr=subprocess.DEVNULL)
    subprocess.run([SAM, "index", synth], check=True)
    os.remove(sam); os.remove(fq)

    merged = f"{OUT}/soto_reads_topup.bam"
    print("merging real scoped BAM + top-up -> soto_reads_topup.bam ...")
    subprocess.run([SAM, "merge", "-f", "-@", "4", merged, REAL_BAM, synth], check=True, stderr=subprocess.DEVNULL)
    subprocess.run([SAM, "index", merged], check=True)

    mapped = sum(1 for a in pysam.AlignmentFile(synth, "rb")
                 if not (a.is_unmapped or a.is_secondary or a.is_supplementary))
    print(f"[done] top-up primary mapped={mapped}/{n_reads}; merged={os.path.getsize(merged)/1e9:.1f}GB")
    print(f"  outputs: {synth} (topup only), {merged} (real+topup), {OUT}/topup_truth.tsv")


if __name__ == "__main__":
    main()
