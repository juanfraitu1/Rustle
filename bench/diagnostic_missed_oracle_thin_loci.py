#!/usr/bin/env python3
"""Diagnostic: which diploid-CN oracle multi-copy genes are missed by gw_xcbase,
and do their loci have thin-locus (support 1-2) read evidence?

Outputs a TSV of missed genes with their annotation coordinates and the count of
primary-read intron-chain groups at support 1, 2, and >=3 in a +/-1 Mb window.

Run: PYTHONHASHSEED=0 /home/juanfra/miniforge3/bin/python bench/diagnostic_missed_oracle_thin_loci.py
"""
import csv
import os
import sys
from collections import defaultdict

if os.environ.get("PYTHONHASHSEED") != "0":
    os.environ["PYTHONHASHSEED"] = "0"
    os.execv(sys.executable, [sys.executable] + sys.argv)

import pysam

SCRATCH = "/home/juanfra/winloci_scratch"
CATALOG_COPY = os.path.join(SCRATCH, "gw_xcbase.copies.tsv")
ANNOT = os.path.join(SCRATCH, "annot_intervals.tsv")
ORACLE = os.path.join("bench", "diploid_cn_oracle.tsv")
BAM = os.path.join(SCRATCH, "GGO.bam")
OUT = os.path.join("bench", "diagnostic_missed_oracle_thin_loci.tsv")


def load_annot():
    gene_to_coords = {}
    with open(ANNOT) as fh:
        next(fh)
        for ln in fh:
            f = ln.rstrip("\n").split("\t")
            c, s, e, g = f[0], int(f[1]), int(f[2]), f[4]
            # keep the longest span for a gene symbol
            if g not in gene_to_coords or (e - s) > (gene_to_coords[g][2] - gene_to_coords[g][1]):
                gene_to_coords[g] = (c, s, e)
    return gene_to_coords


def load_oracle_multicopy_genes():
    genes = []
    with open(ORACLE) as fh:
        next(fh)
        for ln in fh:
            f = ln.rstrip("\n").split("\t")
            if len(f) < 13:
                continue
            gene, cls = f[1], f[12]
            if gene and gene != "NA" and cls == "multi_copy":
                genes.append(gene)
    return sorted(set(genes))


def load_covered_genes(copies_tsv, oracle_set, gene_coords):
    covered = set()
    with open(copies_tsv) as fh:
        r = csv.DictReader(fh, delimiter="\t")
        for row in r:
            chrom, s, e = row["chrom"], int(row["start"]), int(row["end"])
            for g in oracle_set:
                c, gs, ge = gene_coords.get(g, (None, None, None))
                if c is None:
                    continue
                if chrom == c and e >= gs and s <= ge:
                    covered.add(g)
    return covered


def count_thin_loci(bam, chrom, start, end, win=1_000_000):
    lo, hi = max(0, start - win), end + win
    chains = defaultdict(int)
    for aln in bam.fetch(chrom, lo, hi):
        if aln.is_unmapped or aln.is_secondary or aln.is_supplementary:
            continue
        rp = aln.reference_start
        intr = []
        for op, length in aln.cigartuples:
            if op == 3:
                intr.append((rp, rp + length))
                rp += length
            elif op in (0, 2, 7, 8):
                rp += length
        if not intr:
            continue
        chains[tuple(intr)] += 1
    sup1 = sum(1 for v in chains.values() if v == 1)
    sup2 = sum(1 for v in chains.values() if v == 2)
    sup3plus = sum(1 for v in chains.values() if v >= 3)
    return sup1, sup2, sup3plus, len(chains)


def main():
    gene_coords = load_annot()
    oracle_genes = load_oracle_multicopy_genes()
    oracle_set = set(oracle_genes)
    covered = load_covered_genes(CATALOG_COPY, oracle_set, gene_coords)
    missed = sorted(oracle_set - covered)
    print(f"Oracle multi-copy genes: {len(oracle_genes)}")
    print(f"Covered by gw_xcbase: {len(covered)}")
    print(f"Missed: {len(missed)}")
    print("Missed genes:", missed)

    bam = pysam.AlignmentFile(BAM, "rb")
    with open(OUT, "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(["gene", "chrom", "start", "end", "sup1", "sup2", "sup3plus", "n_chains"])
        for g in missed:
            c, s, e = gene_coords.get(g, (None, None, None))
            if c is None:
                w.writerow([g, "NA", "NA", "NA", "NA", "NA", "NA", "NA"])
                continue
            sup1, sup2, sup3, n = count_thin_loci(bam, c, s, e)
            w.writerow([g, c, s, e, sup1, sup2, sup3, n])
            print(f"{g}: {c}:{s}-{e}  thin(1/2/3+)=({sup1}/{sup2}/{sup3})  chains={n}")
    bam.close()
    print(f"\nWrote {OUT}")


if __name__ == "__main__":
    main()
