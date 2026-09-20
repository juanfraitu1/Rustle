#!/usr/bin/env python3
"""Diagnostic: compute exon-colinearity scores between known-family siblings that are
split in the current catalog (MAGEA, HERC2, ZNF92).  Uses the same machinery as the
proposed family_merge_colinear step."""
import os
import sys

if os.environ.get("PYTHONHASHSEED") != "0":
    os.environ["PYTHONHASHSEED"] = "0"
    os.execv(sys.executable, [sys.executable] + sys.argv)

import csv
from collections import defaultdict

import pysam

BENCH = os.path.dirname(os.path.abspath(__file__))
SCRATCH = "/home/juanfra/winloci_scratch"
sys.path.insert(0, BENCH)
import family_er_pr as FP
import recombination_bridge_detector as R
import colinear_multiexon_gate as CM

CATALOG = os.path.join(BENCH, "family_rna_refine.tsv")
GENOME = os.path.join(SCRATCH, "GGO.fasta")
META = os.path.join(SCRATCH, "denovo_transcripts.meta.tsv")


def load_catalog():
    fam = defaultdict(list)
    with open(CATALOG) as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            fam[row["family_id"]].append(dict(
                dn=row["member_dn"], gene=row["member_gene"],
                chrom=row["chrom"], start=int(row["start"]), end=int(row["end"])))
    return fam


def best_colinear_between_blocks(block_a, block_b, skel, strand, fa, id_thresh=R.ID_THRESH):
    """Best strict-LIS colinear shared-exon count and adjacent-junction count between any locus in A and any in B."""
    recs_a = R.family_exons(block_a, skel, strand, fa)
    recs_b = R.family_exons(block_b, skel, strand, fa)
    if not recs_a or not recs_b:
        return 0, 0, len(recs_a), len(recs_b)
    recs = recs_a + recs_b
    best = R.exon_match_tensor(recs)
    best_count = 0
    best_junc = 0
    for ia, ra in enumerate(recs_a):
        for ib_off, rb in enumerate(recs_b):
            ib = len(recs_a) + ib_off
            cnt = CM.colinear_count(best, ia, ib, len(ra["exons"]), strict=True, thresh=id_thresh)
            best_count = max(best_count, cnt)
            # adjacent junctions preserved
            pairs = []
            for i in range(len(ra["exons"])):
                idv, j = best.get((ia, i, ib), (0.0, -1))
                if idv >= id_thresh and j >= 0:
                    pairs.append((i, j))
            pairs.sort()
            junc = 0
            for k in range(1, len(pairs)):
                if pairs[k][0] == pairs[k - 1][0] + 1 and pairs[k][1] == pairs[k - 1][1] + 1:
                    junc += 1
            best_junc = max(best_junc, junc)
    return best_count, best_junc, len(recs_a), len(recs_b)


def main():
    fa = pysam.FastaFile(GENOME)
    skel = R.load_skeletons()
    strand = R.load_strand()
    fam = load_catalog()

    pairs = [
        ("MAGEA", "510", "508", "MAGEA9 sibling"),
        ("HERC2", "385", "384", "HERC2 duplicon"),
        ("ZNF92", "42", "191", "LOC101133668"),
        ("ZNF92", "42", "83", "ZNF626"),
        ("ZNF92", "42", "22", "ZNF208"),
    ]

    thresholds = [0.70, 0.60, 0.55, 0.50]
    print("family\tfid_A\tfid_B\tlabel\t" + "\t".join(f"col_{t:.2f}" for t in thresholds) + "\tjunc_0.70\tn_loci_A\tn_loci_B")
    for name_a, fid_a, fid_b, label in pairs:
        block_a = fam.get(fid_a, [])
        block_b = fam.get(fid_b, [])
        cols = []
        for t in thresholds:
            c, _, _, _ = best_colinear_between_blocks(block_a, block_b, skel, strand, fa, id_thresh=t)
            cols.append(c)
        _, junc, na, nb = best_colinear_between_blocks(block_a, block_b, skel, strand, fa, id_thresh=0.70)
        print(f"{name_a}\t{fid_a}\t{fid_b}\t{label}\t" + "\t".join(map(str, cols)) + f"\t{junc}\t{na}\t{nb}")


if __name__ == "__main__":
    main()
