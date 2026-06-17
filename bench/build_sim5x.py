#!/usr/bin/env python3
"""Synthetic benchmark #2: a gene with 5 near-identical copies in TANDEM in the reference — so a read
has '5 equally good places to put it' (the regime the advisor asks about). A PSV divergence ladder
(K = 0 identical .. 8 private exonic PSVs per copy) maps the identifiability boundary:

  COORDINATES can't assign reads (5 equal places => minimap2 MAPQ 0) almost regardless of K, but
  the PSV allele each read carries CAN assign it to its copy once K PSVs clear the HiFi error floor.

Builds, per K: a genomic contig (5 copies + spacers), spliced HiFi reads per copy (ground-truth copy
in the name), minimap2 -ax splice:hq mapping, and the MAPQ-vs-PSV-assignment summary. Deterministic.
Run with /home/juanfra/miniforge3/bin/python.
"""
import os
import subprocess
import sys
from collections import defaultdict

import pysam

sys.path.insert(0, os.path.dirname(__file__))
from sim_reads import simulate_reads, write_fastq
from copy_specific_junctions import load_gene_exons

OUT = "/home/juanfra/winloci_scratch/sim5x"
FASTA = "/home/juanfra/winloci_scratch/GGO.fasta"
MM2 = "/home/juanfra/miniforge3/bin/minimap2"
SAM = "/home/juanfra/miniforge3/bin/samtools"
BASES = "ACGT"
N_COPIES = 5
READS_PER_COPY = 40
ERR = 0.003          # HiFi substitution error (indel=0 so PSV index is stable = theorem-faithful)
K_LADDER = [0, 1, 2, 4, 8]
SPACER = 2000


def pick_base_gene():
    """Deterministic: a + strand gene, 6-9 exons, mRNA 2000-3500bp, well-expressed (realistic)."""
    nexon = {}
    for line in open("/tmp/gene_chains.tsv"):
        if line.startswith("gene\t"):
            continue
        g, c, st, ne, exl, inl = line.rstrip("\n").split("\t")
        nexon[g] = (int(ne), st, sum(int(x) for x in exl.split(",")) if exl else 0)
    cov = {}
    for line in open("/home/juanfra/winloci_scratch/gene_cov.tsv"):
        if line.startswith("gene\t"):
            continue
        f = line.rstrip("\n").split("\t"); cov[f[0]] = int(f[5])
    cands = [g for g, (ne, st, ml) in nexon.items()
             if st == "+" and 6 <= ne <= 9 and 2000 <= ml <= 3500 and cov.get(g, 0) >= 100]
    return sorted(cands)[0]


def main():
    os.makedirs(OUT, exist_ok=True)
    genes = load_gene_exons()
    base = pick_base_gene()
    chrom, strand, exons = genes[base]           # exons transcription order (+ strand => ascending)
    gstart = min(s for s, _ in exons); gend = max(e for _, e in exons)
    fa = pysam.FastaFile(FASTA)
    gseq = fa.fetch(chrom, gstart, gend).upper()
    ex_rel = [(s - gstart, e - gstart) for s, e in exons]     # exon intervals within gseq
    # mRNA + map mrna_index -> gene offset
    mrna = []; idx2off = []
    for s, e in ex_rel:
        for off in range(s, e):
            mrna.append(gseq[off]); idx2off.append(off)
    mrna = "".join(mrna)
    print(f"base gene {base} {chrom}:{gstart}-{gend} +  exons={len(exons)} mRNA={len(mrna)}bp")

    # K PSV mRNA indices, spread across the transcript
    def psv_indices(K):
        return [ (j + 1) * len(mrna) // (K + 1) for j in range(K) ]

    def allele(ref, c, j):
        # base-4 digits of the copy index across columns: each copy gets a DISTINCT allele vector
        # once K >= ceil(log4(5)) = 2. (K=1 has only 4 bases for 5 copies -> copies 0 & 4 collide.)
        return BASES[(c // (4 ** j)) % 4]

    summary = []
    for K in K_LADDER:
        psv_mrna = psv_indices(K)
        psv_off = [idx2off[i] for i in psv_mrna]      # gene offsets of the PSVs
        # build per-copy genomic + mRNA
        copy_gseq = []; copy_mrna = []
        for c in range(N_COPIES):
            gs = list(gseq); ms = list(mrna)
            for j in range(K):
                gs[psv_off[j]] = allele(gseq[psv_off[j]], c, j)
                ms[psv_mrna[j]] = allele(mrna[psv_mrna[j]], c, j)
            copy_gseq.append("".join(gs)); copy_mrna.append("".join(ms))
        if K == 4:
            saved = dict(copy_mrna=list(copy_mrna), psv_mrna=list(psv_mrna))
        # contig = 5 copies + spacers
        spacer = "A" * SPACER
        contig = spacer.join(copy_gseq)
        unit = len(gseq) + SPACER
        ref_path = f"{OUT}/sim5x_K{K}.ref.fa"
        with open(ref_path, "w") as fh:
            fh.write(f">SIM5X_K{K}\n")
            for i in range(0, len(contig), 80):
                fh.write(contig[i:i+80] + "\n")
        # reads (ground-truth copy in name) + record each read's PSV alleles from the read seq
        fq_path = f"{OUT}/sim5x_K{K}.reads.fq"
        read_psv = {}   # name -> tuple alleles at the K psv mrna indices (with error)
        with open(fq_path, "w") as fh:
            for c in range(N_COPIES):
                reads = simulate_reads(copy_mrna[c], READS_PER_COPY, err=ERR, indel=0.0, seed=1000*K + c)
                for i, rq in enumerate(reads):
                    name = f"K{K}_c{c}_r{i}"
                    write_fastq(fh, name, rq)
                    rd = rq[0]
                    read_psv[name] = tuple(rd[m] if m < len(rd) else "N" for m in psv_mrna)
        # map (spliced HiFi)
        bam = f"{OUT}/sim5x_K{K}.bam"
        p1 = subprocess.run([MM2, "-ax", "splice:hq", "-t", "4", ref_path, fq_path],
                            capture_output=True, text=True)
        with open(f"{OUT}/sim5x_K{K}.sam", "w") as fh:
            fh.write(p1.stdout)
        subprocess.run([SAM, "sort", "-o", bam, f"{OUT}/sim5x_K{K}.sam"], check=True,
                       stderr=subprocess.DEVNULL)
        subprocess.run([SAM, "index", bam], check=True)
        os.remove(f"{OUT}/sim5x_K{K}.sam"); os.remove(fq_path)

        # --- analysis: mapping (coordinates) vs PSV ---
        b = pysam.AlignmentFile(bam, "rb")
        n = mq0 = mismap = 0
        for aln in b.fetch():
            if aln.is_secondary or aln.is_supplementary or aln.is_unmapped:
                continue
            n += 1
            if aln.mapping_quality == 0:
                mq0 += 1
            true_c = int(aln.query_name.split("_c")[1].split("_")[0])
            mapped_c = aln.reference_start // unit
            if aln.mapping_quality > 0 and mapped_c != true_c:
                mismap += 1
        # PSV assignment: assign each read to the copy whose designed allele-vector it best matches
        true_vec = {c: tuple(allele(mrna[psv_mrna[j]], c, j) for j in range(K)) for c in range(N_COPIES)}
        # identifiable copies = # copies pairwise-distinct in designed vectors
        identifiable = sum(1 for c in range(N_COPIES)
                           if all(true_vec[c] != true_vec[d] for d in range(N_COPIES) if d != c))
        correct = total = 0
        for name, obs in read_psv.items():
            tc = int(name.split("_c")[1].split("_")[0]); total += 1
            if K == 0:
                continue  # no PSVs -> unassignable (baseline 1/N)
            # nearest designed vector by Hamming
            best = min(range(N_COPIES), key=lambda c: sum(a != b for a, b in zip(obs, true_vec[c])))
            if best == tc:
                correct += 1
        psv_acc = (correct / total) if (K > 0 and total) else (1.0 / N_COPIES)
        summary.append(dict(K=K, n=n, pct_mq0=round(100*mq0/n, 1) if n else 0,
                            mismap=mismap, identifiable=identifiable, psv_acc=round(psv_acc, 3)))
        print(f"K={K}: reads={n} %MAPQ0={summary[-1]['pct_mq0']} mismapped={mismap} "
              f"identifiable_copies={identifiable}/5 PSV_assign_acc={summary[-1]['psv_acc']}")

    # ---- error-floor sweep at fixed K=4 (PSV-only, no remap): the theorem's 2nd axis ----
    # copies stay identifiable, but PSV assignment degrades once per-base error overwhelms the K columns.
    err_curve = []
    Kf = 4
    pm = saved["psv_mrna"]
    tvec = {c: tuple(allele(mrna[pm[j]], c, j) for j in range(Kf)) for c in range(N_COPIES)}
    for e in (0.003, 0.01, 0.03, 0.06, 0.10, 0.15):
        correct = total = 0
        for c in range(N_COPIES):
            for i, rq in enumerate(simulate_reads(saved["copy_mrna"][c], READS_PER_COPY,
                                                  err=e, indel=0.0, seed=77000 + c)):
                rd = rq[0]; obs = tuple(rd[m] if m < len(rd) else "N" for m in pm)
                best = min(range(N_COPIES), key=lambda cc: sum(a != b for a, b in zip(obs, tvec[cc])))
                total += 1; correct += (best == c)
        err_curve.append((e, round(correct / total, 3)))
        print(f"  err-sweep K=4 e={e}: PSV_assign_acc={err_curve[-1][1]}")

    # write summary
    import json
    with open(f"{OUT}/sim5x_summary.tsv", "w") as fh:
        fh.write("K\treads\tpct_MAPQ0\tmismapped\tidentifiable_copies\tPSV_assign_acc\n")
        for s in summary:
            fh.write(f"{s['K']}\t{s['n']}\t{s['pct_mq0']}\t{s['mismap']}\t{s['identifiable']}\t{s['psv_acc']}\n")
    with open(f"{OUT}/sim5x_errsweep.tsv", "w") as fh:
        fh.write("error\tPSV_assign_acc_K4\n")
        for e, a in err_curve:
            fh.write(f"{e}\t{a}\n")
    print(f"\n[wrote {OUT}/ : ref.fa + bam per K + sim5x_summary.tsv + sim5x_errsweep.tsv]  base_gene={base}")
    json.dump({"base_gene": base, "summary": summary, "err_curve": err_curve},
              open(f"{OUT}/sim5x_summary.json", "w"), indent=2)


if __name__ == "__main__":
    main()
