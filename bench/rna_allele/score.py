#!/usr/bin/env python3
"""Score the RNA-only allele calls against the frozen truth (docs/archive/2026-10/PREREG_rna_allele_haplotype_count_2026-10-01.md, H1, H1b, H2 and the
blind spot; H3/H4 need every chromosome holding a family and are reported when those are present).

    python3 score.py --truth truth.tsv --calls-dir calls [--chroms NC_073244.2,NC_073247.2]
"""
import argparse
import collections
import csv
import glob
import os


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--truth", required=True)
    ap.add_argument("--calls-dir", required=True)
    ap.add_argument("--chroms", default="")
    a = ap.parse_args(argv)
    truth = {r["gene_id"]: r for r in csv.DictReader(open(a.truth), delimiter="\t")}
    calls = {}
    for path in sorted(glob.glob(f"{a.calls_dir}/*.calls.tsv")):
        c = os.path.basename(path).split(".calls")[0]
        if a.chroms and c not in a.chroms.split(","):
            continue
        for r in csv.DictReader(open(path), delimiter="\t"):
            calls[r["gene_id"]] = r
    print(f"genes with calls: {len(calls):,}")
    by = collections.defaultdict(collections.Counter)
    for g, c in calls.items():
        t = truth.get(g)
        if not t:
            continue
        by[c["set"]][(t["class"], c["call"])] += 1
    for st in ("S_X", "S_single", "S_multi", "S_fam"):
        k = by.get(st)
        if not k:
            continue
        expressed = {cl: sum(v for (tc, call), v in k.items() if tc == cl and call != "NA") for cl in ("T2d", "T2i", "T1", "T?")}
        two = {cl: k[(cl, "2")] for cl in ("T2d", "T2i", "T1", "T?")}
        inc = sum(v for (tc, call), v in k.items() if call == "inconsistent")
        n2 = sum(two[cl] for cl in ("T2d", "T2i", "T1"))
        print(f"\n== {st}: expressed {sum(expressed.values()):,} (T2d {expressed['T2d']}, T2i {expressed['T2i']}, T1 {expressed['T1']}, "
              f"T? {expressed['T?']}); called 2: {sum(two.values())} (T2d {two['T2d']}, T2i {two['T2i']}, T1 {two['T1']}, T? {two['T?']}); "
              f"inconsistent {inc}")
        if st == "S_X":
            fx = two["T1"] / expressed["T1"] if expressed["T1"] else 0
            print(f"   H1 f_X = {two['T1']}/{expressed['T1']} = {fx:.4f} -> {'PASS' if fx <= 0.02 else 'FAIL' if fx > 0.05 else 'marginal'}")
            continue
        prec = two["T2d"] / n2 if n2 else float("nan")
        f1 = two["T1"] / expressed["T1"] if expressed["T1"] else float("nan")
        fi = two["T2i"] / expressed["T2i"] if expressed["T2i"] else float("nan")
        rec = two["T2d"] / expressed["T2d"] if expressed["T2d"] else float("nan")
        print(f"   precision of '2' = {two['T2d']}/{n2} = {prec:.3f}; false-2 on T1 = {two['T1']}/{expressed['T1']} = {f1:.3f}; "
              f"false-2 on T2i = {two['T2i']}/{expressed['T2i']} = {fi:.3f}; recall on T2d = {rec:.3f}")
        if st in ("S_multi", "S_fam"):
            v = ("HOLDS" if prec >= 0.90 and (f1 != f1 or f1 <= 0.10) else
                 "PARTIAL" if prec >= 0.75 and (f1 != f1 or f1 <= 0.25) else "FAILS")
            print(f"   H2 ({st}) -> {v}")
        blind_t1 = sum(v for (tc, call), v in k.items() if tc == "T1" and call == "1+")
        blind_t2i = sum(v for (tc, call), v in k.items() if tc == "T2i" and call == "1+")
        print(f"   blind spot: expressed T1 called 1+ = {blind_t1}; expressed T2i called 1+ = {blind_t2i}")


if __name__ == "__main__":
    main()
