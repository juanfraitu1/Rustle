#!/usr/bin/env python3
"""eval_psv_catalog.py — evaluate PSV-filtered catalog against known multi-copy gene families.

Ground truth is a curated set of classic multi-copy gene families (GSTM, RABL2, IFITM, etc.)
plus any gene symbol that appears at >=2 distinct RefSeq gene loci.
"""
import csv
import os
import re
import sys
from collections import defaultdict

HERE = os.path.dirname(os.path.abspath(__file__))
REFINE_TSV = os.path.join(HERE, "family_rna_refine.tsv")
FILTERED_TSV = os.path.join(HERE, "family_rna_refine_psv_filtered.tsv")
REJECTED_TSV = os.path.join(HERE, "family_rna_refine_psv_rejected.tsv")
GFF = "/home/juanfra/winloci_scratch/GGO_tx.gff"

# Curated classic multi-copy / paralog families in primates
CURATED_FAMILIES = {
    "GSTM": ["GSTM1", "GSTM2", "GSTM3", "GSTM4", "GSTM5"],
    "GSTT": ["GSTT1", "GSTT2", "GSTT2B", "GSTT3", "GSTT4"],
    "GSTP": ["GSTP1"],
    "RABL2": ["RABL2A", "RABL2B"],
    "IFITM": ["IFITM1", "IFITM2", "IFITM3", "IFITM5", "IFITM6", "IFITM7", "IFITM10"],
    "HLA": ["HLA-A", "HLA-B", "HLA-C", "HLA-DMA", "HLA-DMB", "HLA-DOA", "HLA-DOB",
            "HLA-DPA1", "HLA-DPB1", "HLA-DQA1", "HLA-DQB1", "HLA-DRA", "HLA-DRB1",
            "HLA-E", "HLA-F", "HLA-G"],
    "PSG": ["PSG1", "PSG2", "PSG3", "PSG4", "PSG5", "PSG6", "PSG7", "PSG8", "PSG9", "PSG10", "PSG11"],
    "UGT": ["UGT1A1", "UGT1A3", "UGT1A4", "UGT1A5", "UGT1A6", "UGT1A7", "UGT1A8",
            "UGT1A9", "UGT1A10", "UGT2A1", "UGT2A2", "UGT2A3", "UGT2B4", "UGT2B7",
            "UGT2B10", "UGT2B11", "UGT2B15", "UGT2B17", "UGT2B28", "UGT3A1", "UGT3A2", "UGT8"],
    "KRT": ["KRT" + x for x in ["1", "2", "3", "4", "5", "6A", "6B", "6C", "7", "8", "9",
            "10", "12", "13", "14", "15", "16", "17", "18", "19", "20", "23", "24",
            "25", "26", "27", "28", "31", "32", "33A", "33B", "34", "35", "36", "37",
            "38", "71", "72", "73", "74", "75", "76", "77", "78", "79", "80", "81",
            "82", "83", "84", "85", "86", "87", "88", "89", "90", "222"]],
    "HSP": ["HSPA1A", "HSPA1B", "HSPA1L", "HSPA2", "HSPA4", "HSPA4L", "HSPA5",
            "HSPA6", "HSPA8", "HSPA9", "HSPA12A", "HSPA12B", "HSPA13", "HSPA14"],
    "TRIM": ["TRIM" + str(i) for i in range(1, 80)],
    "TBC1D3": ["TBC1D3" + s for s in ["", "B", "C", "D", "E", "F", "G", "H", "I", "K", "L"]],
    "DEF": ["DEFA1", "DEFA3", "DEFA4", "DEFA5", "DEFA6", "DEFB1", "DEFB4A", "DEFB4B"],
    "PRAME": ["PRAME"],
    "CCL": ["CCL" + str(i) for i in range(1, 29)],
    "CXCL": ["CXCL" + str(i) for i in range(1, 18)],
}


def load_ground_truth_symbols():
    """Return set of known multi-copy gene symbols: curated + >=2 RefSeq loci."""
    symbols = set()
    for fam, genes in CURATED_FAMILIES.items():
        for g in genes:
            symbols.add(g)

    # also add any named symbol with >=2 gene loci in RefSeq
    counts = defaultdict(int)
    with open(GFF) as f:
        for line in f:
            if line.startswith("#"):
                continue
            p = line.rstrip().split("\t")
            if len(p) < 9 or p[2] != "gene":
                continue
            m = re.search(r'gene=([^;]+)', p[8])
            if not m:
                continue
            g = m.group(1)
            if g.startswith("LOC") or g.startswith("TRNA") or g.startswith("RNA"):
                continue
            counts[g] += 1
    for g, c in counts.items():
        if c >= 2:
            symbols.add(g)
    return symbols


def load_catalog_families(path):
    """Return {family_id: set(gene_symbols)}."""
    fam_genes = defaultdict(set)
    with open(path) as f:
        for r in csv.DictReader(f, delimiter="\t"):
            fam_genes[int(r["family_id"])].add(r["member_gene"])
            fam_genes[int(r["family_id"])].add(r["dominant_gene"])
    return fam_genes


def evaluate(path, gt_symbols, name):
    fam_genes = load_catalog_families(path)
    n_fams = len(fam_genes)
    tp = 0
    fp = 0
    found_symbols = set()
    for fid, genes in fam_genes.items():
        if genes & gt_symbols:
            tp += 1
            found_symbols.update(genes & gt_symbols)
        else:
            fp += 1
    precision = tp / n_fams if n_fams else 0.0
    recall = len(found_symbols) / len(gt_symbols) if gt_symbols else 0.0
    f1 = 2 * precision * recall / (precision + recall) if (precision + recall) else 0.0
    print(f"\n{name}")
    print(f"  families: {n_fams}")
    print(f"  TP (contains known symbol): {tp}")
    print(f"  FP (no known symbol):       {fp}")
    print(f"  Precision: {precision:.3f}")
    print(f"  Known symbols found: {len(found_symbols)} / {len(gt_symbols)}  Recall: {recall:.3f}")
    print(f"  F1: {f1:.3f}")
    return tp, fp, found_symbols


def main():
    gt_symbols = load_ground_truth_symbols()
    print(f"Ground-truth multi-copy gene symbols: {len(gt_symbols)}")

    evaluate(FILTERED_TSV, gt_symbols, "PSV-filtered catalog (402 families)")
    all_fam_genes = load_catalog_families(REFINE_TSV)
    # for original catalog, use same gt_symbols
    evaluate(REFINE_TSV, gt_symbols, "Original catalog (551 families)")

    # show some FP examples from filtered catalog
    print("\nFiltered-catalog families with NO known multi-copy symbol (sample FP):")
    filtered = load_catalog_families(FILTERED_TSV)
    fp_examples = []
    for fid, genes in filtered.items():
        if not (genes & gt_symbols):
            fp_examples.append((fid, sorted(genes - {"", "NA"})))
    for fid, genes in fp_examples[:15]:
        print(f"  fam{fid}: {genes}")


if __name__ == "__main__":
    main()
