#!/usr/bin/env python3
"""Evaluate a gw_* catalog vs baseline on oracle genes and known-family windows."""
import csv
import os
import sys

if os.environ.get("PYTHONHASHSEED") != "0":
    os.environ["PYTHONHASHSEED"] = "0"
    os.execv(sys.executable, [sys.executable] + sys.argv)

SCRATCH = "/home/juanfra/winloci_scratch"
ANNOT = os.path.join(SCRATCH, "annot_intervals.tsv")
ORACLE = os.path.join("bench", "diploid_cn_oracle.tsv")
SHOWCASE = os.path.join("bench", "KNOWN_FAMILY_SHOWCASE.md")


def load_annot():
    gene_to_coords = {}
    with open(ANNOT) as fh:
        next(fh)
        for ln in fh:
            f = ln.rstrip("\n").split("\t")
            c, s, e, g = f[0], int(f[1]), int(f[2]), f[4]
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


def covered_oracle_genes(copies_tsv, oracle_set, gene_coords):
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


def catalog_stats(copies_tsv, families_tsv):
    with open(families_tsv) as fh:
        n_fam = sum(1 for _ in fh) - 1
    with open(copies_tsv) as fh:
        r = csv.DictReader(fh, delimiter="\t")
        copies = list(r)
    n_copies = len(copies)
    n_xchrom = len(set(c["family_id"] for c in copies))  # rough
    return n_fam, n_copies


def main():
    import argparse
    ap = argparse.ArgumentParser()
    ap.add_argument("--baseline-copies", default=os.path.join(SCRATCH, "gw_xcbase.copies.tsv"))
    ap.add_argument("--baseline-families", default=os.path.join(SCRATCH, "gw_xcbase.families.tsv"))
    ap.add_argument("--test-copies", default=os.path.join(SCRATCH, "gw_rescue.copies.tsv"))
    ap.add_argument("--test-families", default=os.path.join(SCRATCH, "gw_rescue.families.tsv"))
    args = ap.parse_args()

    gene_coords = load_annot()
    oracle_genes = load_oracle_multicopy_genes()
    oracle_set = set(oracle_genes)

    base_cov = covered_oracle_genes(args.baseline_copies, oracle_set, gene_coords)
    test_cov = covered_oracle_genes(args.test_copies, oracle_set, gene_coords)

    base_fam, base_copy = catalog_stats(args.baseline_copies, args.baseline_families)
    test_fam, test_copy = catalog_stats(args.test_copies, args.test_families)

    print(f"Oracle multi-copy genes: {len(oracle_set)}")
    print(f"Baseline covered: {len(base_cov)} ({len(base_cov)/len(oracle_set):.3f})")
    print(f"Test covered:     {len(test_cov)} ({len(test_cov)/len(oracle_set):.3f})")
    print(f"Newly covered:    {sorted(test_cov - base_cov)}")
    print(f"Lost:             {sorted(base_cov - test_cov)}")
    print(f"Still missed:     {sorted(oracle_set - test_cov)}")
    print(f"Families:         {base_fam} -> {test_fam} ({test_fam - base_fam:+d})")
    print(f"Copies:           {base_copy} -> {test_copy} ({test_copy - base_copy:+d})")


if __name__ == "__main__":
    main()
