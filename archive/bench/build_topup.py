#!/usr/bin/env python3
"""Synthetic benchmark #1: 'ideal coverage' GGO. Simulate ONLY what the real IsoSeq lacks — top up
every under-covered gene to TARGET reads from its representative transcript — then map to the genome
and merge with the real BAM, so every gene is well-covered (coverage never the limiting factor for
recall / copy / allele analyses). Disk-light: simulate -> map -> delete FASTQ; keep only BAMs.
Run with /home/juanfra/miniforge3/bin/python.
"""
import os
import subprocess
import sys

import pysam

sys.path.insert(0, os.path.dirname(__file__))
from sim_reads import simulate_reads, write_fastq
from copy_specific_junctions import load_gene_exons

OUT = "/home/juanfra/winloci_scratch/topup"
FASTA = "/home/juanfra/winloci_scratch/GGO.fasta"
REAL_BAM = "/home/juanfra/winloci_scratch/GGO.bam"
COV = "/home/juanfra/winloci_scratch/gene_cov.tsv"
MM2 = "/home/juanfra/miniforge3/bin/minimap2"
SAM = "/home/juanfra/miniforge3/bin/samtools"
TARGET = 40
ERR, INDEL = 0.003, 0.001
_RC = bytes.maketrans(b"ACGTacgtNn", b"TGCAtgcaNn")


def main():
    os.makedirs(OUT, exist_ok=True)
    genes = load_gene_exons()
    fa = pysam.FastaFile(FASTA)
    need = {}
    for line in open(COV):
        if line.startswith("gene\t"):
            continue
        f = line.rstrip("\n").split("\t"); g = f[0]; nr = int(f[5])
        if nr < TARGET:
            need[g] = TARGET - nr
    fq = f"{OUT}/topup.fq"
    n_reads = n_genes = 0
    truth = open(f"{OUT}/topup_truth.tsv", "w"); truth.write("read\tgene\tchrom\n")
    with open(fq, "w") as fh:
        for g in sorted(need):
            if g not in genes:
                continue
            chrom, strand, exons = genes[g]
            seq = "".join(fa.fetch(chrom, s, e).upper() for s, e in sorted(exons))
            if strand == "-":
                seq = seq.encode().translate(_RC)[::-1].decode()
            if len(seq) < 100:
                continue
            for i, rq in enumerate(simulate_reads(seq, need[g], err=ERR, indel=INDEL,
                                                  seed=abs(hash(g)) % 2_000_000, trunc_frac=0.10)):
                name = f"SIMTOPUP|{g}|{i}"
                write_fastq(fh, name, rq); truth.write(f"{name}\t{g}\t{chrom}\n")
                n_reads += 1
            n_genes += 1
    truth.close()
    print(f"simulated {n_reads} top-up reads for {n_genes} under-covered genes -> {fq}")

    synth_bam = f"{OUT}/topup.bam"
    print("mapping to genome (minimap2 -ax splice:hq)...")
    with open(f"{OUT}/topup.sam", "w") as sf:
        subprocess.run([MM2, "-ax", "splice:hq", "-t", "5", FASTA, fq], stdout=sf,
                       stderr=subprocess.DEVNULL, check=True)
    subprocess.run([SAM, "sort", "-@", "4", "-o", synth_bam, f"{OUT}/topup.sam"], check=True,
                   stderr=subprocess.DEVNULL)
    subprocess.run([SAM, "index", synth_bam], check=True)
    os.remove(f"{OUT}/topup.sam"); os.remove(fq)   # disk-light: drop SAM + FASTQ

    ideal = f"{OUT}/GGO_ideal.bam"
    print("merging real + synthetic -> GGO_ideal.bam ...")
    subprocess.run([SAM, "merge", "-f", "-@", "4", ideal, REAL_BAM, synth_bam], check=True,
                   stderr=subprocess.DEVNULL)
    subprocess.run([SAM, "index", ideal], check=True)

    mapped = sum(1 for a in pysam.AlignmentFile(synth_bam, "rb")
                 if not (a.is_unmapped or a.is_secondary or a.is_supplementary))
    sz = os.path.getsize(ideal) / 1e9
    print(f"[done] synthetic primary mapped={mapped}/{n_reads} "
          f"({100*mapped/max(n_reads,1):.1f}%); GGO_ideal.bam={sz:.1f}GB")
    print(f"  outputs: {synth_bam} (synthetic only), {ideal} (real+synthetic), {OUT}/topup_truth.tsv")


if __name__ == "__main__":
    main()
