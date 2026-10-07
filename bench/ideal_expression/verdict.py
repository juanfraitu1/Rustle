#!/usr/bin/env python3
"""Combine the two replicates of an arm (docs/PREREG_ideal_expression_2026-10-06.md): per family the verdict is the LOWER of the two (NO < PARTLY < YES), UNSTABLE when they differ,
CEILING-LIMITED when either replicate is; lists the copies whose IDEAL-FOUND flag differs between the replicates.

    verdict.py --rep1 SCOREPREFIX1 --rep2 SCOREPREFIX2 --out OUT.json
"""
import argparse
import csv
import json

ORDER = {"NO": 0, "PARTLY": 1, "YES": 2}


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--rep1", required=True)
    ap.add_argument("--rep2", required=True)
    ap.add_argument("--out", required=True)
    a = ap.parse_args()
    s1, s2 = json.load(open(a.rep1 + ".summary.json")), json.load(open(a.rep2 + ".summary.json"))
    out = {}
    for fam in s1:
        if fam == "G7":
            out[fam] = dict(rep1=s1[fam], rep2=s2.get(fam))
            continue
        v1, v2 = s1[fam]["rule"], s2[fam]["rule"]
        if "CEILING-LIMITED" in (v1, v2):
            v = "CEILING-LIMITED"
        elif v1 == v2:
            v = v1
        else:
            v = "UNSTABLE (" + v1 + " / " + v2 + "; lower = " + (v1 if ORDER[v1] <= ORDER[v2] else v2) + ")"
        f1 = {r["cid"]: r for r in csv.DictReader(open(f"{a.rep1}.{fam}.copies.tsv"), delimiter="\t")}
        f2 = {r["cid"]: r for r in csv.DictReader(open(f"{a.rep2}.{fam}.copies.tsv"), delimiter="\t")}
        flips = sorted(c for c in f1 if c in f2 and f1[c]["IDEAL_FOUND"] != f2[c]["IDEAL_FOUND"])
        out[fam] = dict(verdict=v, rep1=dict(rule=v1, R=s1[fam]["R"], found_R=s1[fam]["R_counts"]["IDEAL_FOUND"], found_ALL=s1[fam]["ALL_counts"]["IDEAL_FOUND"], N=s1[fam]["N"], CP=s1[fam]["CP"]),
                        rep2=dict(rule=v2, R=s2[fam]["R"], found_R=s2[fam]["R_counts"]["IDEAL_FOUND"], found_ALL=s2[fam]["ALL_counts"]["IDEAL_FOUND"], N=s2[fam]["N"], CP=s2[fam]["CP"]),
                        copies_flipping_between_replicates=flips)
        print(f"{fam}: {v} | rep1 {v1} (R {s1[fam]['R']}, found {s1[fam]['R_counts']['IDEAL_FOUND']}, all {s1[fam]['ALL_counts']['IDEAL_FOUND']}/{s1[fam]['N']}) | "
              f"rep2 {v2} (R {s2[fam]['R']}, found {s2[fam]['R_counts']['IDEAL_FOUND']}, all {s2[fam]['ALL_counts']['IDEAL_FOUND']}/{s2[fam]['N']}) | flips {flips}")
    json.dump(out, open(a.out, "w"), indent=1)


if __name__ == "__main__":
    main()
