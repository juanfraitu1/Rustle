#!/usr/bin/env python3
"""Post hoc (docs/SEED_POOL_REAL_READS_2026-10-07.md section 5): at each copy that is expressed under Amendment A, how many PRIMARY and how many
SECONDARY alignments carry >= k of the copy's annotated introns? The answer says which pool can see the copy: a copy whose evidence is a handful of
primaries and dozens of secondaries is only buildable from the secondary pool.

    secondary_support.py --copies copies.tsv --truth truth.gtf --bam reads.bam --support support.copies.tsv --family NPIP --out secondary_support.tsv

For every copy with `ann_expressed` == 1 in --support: k = min(2, the copy's longest annotated intron chain); an alignment counts when it is on the
copy's strand, overlaps the copy's exon union, is not supplementary and carries >= k introns of the copy's annotated set (exact donor / acceptor,
introns >= 50 bp). Prints one row per copy and the medians; writes the TSV (columns: name, k, primary, secondary, molecules).
"""
import argparse
import csv
import os
import statistics
import sys

import pysam

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))
import copy_support as cs  # noqa: E402


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--copies", required=True)
    ap.add_argument("--truth", required=True)
    ap.add_argument("--bam", required=True)
    ap.add_argument("--support", required=True, help="copy_support.py's PREFIX.copies.tsv (its ann_expressed column picks the copies)")
    ap.add_argument("--family", required=True)
    ap.add_argument("--out", required=True)
    a = ap.parse_args()
    copies = [r for r in csv.DictReader(open(a.copies), delimiter="\t") if r["family"] == a.family]
    tx = cs.gtf_transcripts(a.truth, {r["isoform_gene"] for r in copies} | {r["cid"] for r in copies})
    expressed = {r["name"] for r in csv.DictReader(open(a.support), delimiter="\t") if r["ann_expressed"] == "1"}
    bam = pysam.AlignmentFile(a.bam)
    rows = []
    for c in copies:
        name = c.get("refseq_name") or c.get("cat_name") or c["name"]
        if name not in expressed:
            continue
        texons = tx.get(c["cid"]) or tx.get(c["isoform_gene"], {})
        ann_all = set().union(*[set(cs.introns_of(ex)) for ex in texons.values()])
        k = min(2, max(len(cs.introns_of(ex)) for ex in texons.values()))
        union = cs.merge([e for ex in texons.values() for e in ex])
        lo, hi = union[0][0], union[-1][1]
        prim = sec = 0
        mols = set()
        for rd in bam.fetch(c["chrom"], lo, hi):
            if rd.is_unmapped or rd.is_supplementary or ("-" if rd.is_reverse else "+") != c["strand"]:
                continue
            blocks, juncs = cs.read_blocks_junctions(rd)
            if cs.inter(blocks, union) == 0:
                continue
            if sum(1 for j in juncs if j[0] >= lo and j[1] <= hi and j in ann_all) >= k:
                if rd.is_secondary:
                    sec += 1
                else:
                    prim += 1
                mols.add(rd.query_name)
        rows.append((name, k, prim, sec, len(mols)))
    with open(a.out, "w") as fh:
        fh.write("name\tk\tprimary\tsecondary\tmolecules\n")
        for r in rows:
            fh.write("\t".join(str(x) for x in r) + "\n")
    print(f"{'copy':14}{'k':>2}{'primary':>9}{'secondary':>11}{'molecules':>11}")
    for r in rows:
        print(f"{r[0]:14}{r[1]:>2}{r[2]:>9}{r[3]:>11}{r[4]:>11}")
    if rows:
        print(f"{a.family}: {len(rows)} copies; median primary {statistics.median(r[2] for r in rows):g}, median secondary {statistics.median(r[3] for r in rows):g}; "
              f"copies with more secondary than primary alignments: {sum(1 for r in rows if r[3] > r[2])}")


if __name__ == "__main__":
    main()
