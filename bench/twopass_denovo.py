#!/usr/bin/env python3
"""ANNOTATION-FREE two-pass prototype (read-coherence -> family + copy assignment), on the collapsed
5-copy regime where read-coherence's value over flow is starkest:

  PASS 1 (read-coherence): map all copies' reads to a SINGLE-copy reference (collapsed) and group reads
    by exact intron chain -> transcript SKELETONS, keeping each read's sequence (so PSVs survive). Flow
    would emit ~1 transcript here and LOSE the copies; read-coherence keeps one skeleton + its reads.
  PASS 2: call PSVs DE NOVO from the skeleton's read pileup (no annotation), split reads by PSV
    allele-vector into copies (copy_split), and declare the multi-copy family. Assign each read (incl.
    hard ones) to a copy; accuracy vs ground truth, across the PSV ladder K (the identifiability boundary).

No annotation, no StringTie. Run with /home/juanfra/miniforge3/bin/python. Deterministic.
"""
import math
import os
import sys
import subprocess
from collections import defaultdict, Counter

import pysam

sys.path.insert(0, os.path.dirname(__file__))
from sim_reads import simulate_reads
from copy_specific_junctions import load_gene_exons

OUT = "/home/juanfra/winloci_scratch/twopass"
FASTA = "/home/juanfra/winloci_scratch/GGO.fasta"
MM2 = "/home/juanfra/miniforge3/bin/minimap2"
SAM = "/home/juanfra/miniforge3/bin/samtools"
BASES = "ACGT"
N_COPIES, READS_PER_COPY, ERR = 5, 40, 0.003
K_LADDER = [0, 1, 2, 4, 8]
D_MIN, C_MIN, F_MIN = 8, 3, 0.15
MIN_READS = 3


def introns_of(aln):
    out = []; rpos = aln.reference_start
    for op, length in aln.cigartuples:
        if op == 3:
            out.append((rpos, rpos + length)); rpos += length
        elif op in (0, 2, 7, 8):
            rpos += length
    return tuple(out)


def call_psvs(bam, chrom, lo, hi):
    """De-novo PSV columns from the collapsed pileup + per-read alleles."""
    col_reads = {}
    for pc in bam.pileup(chrom, lo, hi, truncate=True, stepper="nofilter", min_base_quality=0):
        tally = Counter(); per = {}
        for pr in pc.pileups:
            if pr.is_del or pr.is_refskip or pr.query_position is None:
                continue
            seq = pr.alignment.query_sequence
            if seq is None:
                continue
            b = seq[pr.query_position]
            if b in BASES:
                tally[b] += 1; per[pr.alignment.query_name] = b
        d = sum(tally.values())
        if d < D_MIN or len(tally) < 2:
            continue
        minor = tally.most_common(2)[1][1]
        if minor >= C_MIN and minor / d >= F_MIN:
            col_reads[pc.reference_pos] = per
    cols = sorted(col_reads)
    read_vec = defaultdict(dict)
    for ci, c in enumerate(cols):
        for rn, b in col_reads[c].items():
            read_vec[rn][ci] = b
    return cols, read_vec


def main():
    os.makedirs(OUT, exist_ok=True)
    genes = load_gene_exons()
    chrom, strand, exons = genes["AASDHPPT"]
    gstart = min(s for s, _ in exons); gend = max(e for _, e in exons)
    fa = pysam.FastaFile(FASTA)
    gseq = fa.fetch(chrom, gstart, gend).upper()                 # SINGLE-copy genomic reference
    mrna = "".join(gseq[s - gstart:e - gstart] for s, e in exons)
    ref_path = f"{OUT}/ref_single.fa"
    with open(ref_path, "w") as fh:
        fh.write(">COLLAPSED\n")
        for i in range(0, len(gseq), 80):
            fh.write(gseq[i:i+80] + "\n")

    def allele(c, j):
        return BASES[(c // (4 ** j)) % 4]

    print(f"{'K':>3} {'skeletons':>10} {'denovo_PSV':>11} {'copies_rec':>11} {'assign_acc':>11}  note")
    rows = []
    for K in K_LADDER:
        psv_idx = [(j + 1) * len(mrna) // (K + 1) for j in range(K)]
        # 5 copy mRNAs; simulate reads; write FASTQ (read name carries TRUE copy)
        fq = f"{OUT}/reads_K{K}.fq"
        truth = {}
        with open(fq, "w") as fh:
            for c in range(N_COPIES):
                ms = list(mrna)
                for j in range(K):
                    ms[psv_idx[j]] = allele(c, j)
                cmrna = "".join(ms)
                for i, (rd, q) in enumerate(simulate_reads(cmrna, READS_PER_COPY, err=ERR, indel=0.0,
                                                           seed=4000 * K + c)):
                    nm = f"K{K}_c{c}_r{i}"; truth[nm] = c
                    fh.write(f"@{nm}\n{rd}\n+\n{q}\n")
        bam = f"{OUT}/K{K}.bam"
        with open(f"{OUT}/K{K}.sam", "w") as sf:
            subprocess.run([MM2, "-ax", "splice:hq", "-t", "4", ref_path, fq],
                           stdout=sf, stderr=subprocess.DEVNULL, check=True)
        subprocess.run([SAM, "sort", "-o", bam, f"{OUT}/K{K}.sam"], check=True, stderr=subprocess.DEVNULL)
        subprocess.run([SAM, "index", bam], check=True)
        os.remove(f"{OUT}/K{K}.sam"); os.remove(fq)

        b = pysam.AlignmentFile(bam, "rb")
        # PASS 1: read-coherence — group reads by intron chain (skeletons)
        chains = defaultdict(list)
        for aln in b.fetch():
            if aln.is_secondary or aln.is_supplementary or aln.is_unmapped:
                continue
            chains[introns_of(aln)].append(aln.query_name)
        skeletons = [c for c in chains.values() if len(c) >= MIN_READS]
        dom = max(skeletons, key=len)               # the collapsed transcript skeleton
        dom_set = set(dom)

        # PASS 2: de-novo PSVs from the skeleton's collapsed pileup + copy split
        # (reads are mapped to the single-copy ref contig "COLLAPSED", coords 0..len(gseq))
        cols, read_vec = call_psvs(b, "COLLAPSED", 0, len(gseq))
        # recovered copy haplotypes = distinct allele vectors (over covered cols) with >= MIN_READS
        vec_counts = Counter()
        read_full = {}
        for rn in dom_set:
            v = tuple(read_vec[rn].get(i) for i in range(len(cols)))
            read_full[rn] = v
            if any(x is not None for x in v):
                vec_counts[v] += 1
        haplos = [v for v, n in vec_counts.items() if n >= MIN_READS]
        # map each haplotype to the true copy its reads mostly come from; assign reads to nearest haplo
        def dist(a, b):
            return sum(1 for x, y in zip(a, b) if x is not None and y is not None and x != y)
        # consensus per haplotype already = the vector; assign read to nearest haplo
        correct = total = 0
        haplo_truth = {}
        # majority truth per haplotype
        for h in haplos:
            members = [rn for rn in dom_set if read_full[rn] == h]
            tc = Counter(truth[rn] for rn in members).most_common(1)
            haplo_truth[h] = tc[0][0] if tc else -1
        n_rec = len(set(haplo_truth.values())) if haplos else 0
        for rn in dom_set:
            v = read_full[rn]
            if not haplos or all(x is None for x in v):
                continue
            nearest = min(haplos, key=lambda h: dist(v, h))
            total += 1
            correct += (haplo_truth[nearest] == truth[rn])
        acc = (correct / total) if total else 0.0
        note = "identical->unassignable" if K == 0 else ("flow would give 1 transcript; we recover copies")
        rows.append((K, len(skeletons), len(cols), n_rec, round(acc, 3)))
        print(f"{K:>3} {len(skeletons):>10} {len(cols):>11} {n_rec:>11} {acc:>11.3f}  {note}")

    import json
    json.dump({"rows": rows}, open(f"{OUT}/twopass_summary.json", "w"), indent=2)
    print("\nPASS 1 collapses all copies' reads into ONE skeleton (flow would stop here = copies lost).")
    print("PASS 2 calls PSVs de novo and splits the skeleton into copies; recovery + assignment follow")
    print("the identifiability boundary (K=0 identical = unassignable; resolves as PSVs accrue). NO annotation.")


if __name__ == "__main__":
    main()
