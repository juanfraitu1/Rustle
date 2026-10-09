#!/usr/bin/env python3
"""Amendment 39 report: end-support flags of every cluster of the 16 control runs and the 4 real runs (local-mode consensus). Miniforge python."""
import collections
import csv
import json
import os
import statistics as st
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
sys.path.insert(0, os.path.dirname(HERE))
import discover as D  # noqa: E402
import endtrim as E  # noqa: E402
import run_polish as RP  # noqa: E402
import seeds as SD  # noqa: E402

GOOD = ("PRESENT", "ALLELE-LIKE")


def main():
    animals = ("ggo_testis", "a119b", "ptr", "ppy")
    runs = [a + s for a in animals for s in ("", "_control", "_control_s6", "_control_s7")]
    rows = []
    for run in runs:
        d = f"{D.W}/discover_{run}"
        cl = collections.defaultdict(list)
        for r in csv.DictReader(open(f"{d}/clusters.tsv"), delimiter="\t"):
            cl["cl" + r["cluster"]].append(r["read"])
        reads = SD.read_fa(open(f"{d}/reads_path.txt").read().strip())
        cons = {n.split("|")[0]: s for n, s in RP.read_cons(D.variant_paths(d, "l")["cons"]).items()}
        cg = {r["k"]: r for r in json.load(open(D.variant_paths(d, "")["rescored"]))}
        cls = {r["k"]: r for r in json.load(open(D.variant_paths(d, "l")["rescored"]))}
        for k in cls:
            names = sorted(cl[k], key=lambda r: -len(reads[r]))[:100]
            f = E.end_support(cons[k], [reads[r] for r in names], f"{d}/tmp_ends")
            rows.append(dict(run=run, k=k, real="control" not in run, cls=cls[k]["cls"], cls_global=cg[k]["cls"], length=cls[k]["length"], reads=cls[k]["reads"], **f))
    json.dump(rows, open(f"{D.W}/end_flags.json", "w"))
    for label, sel in (("control runs", [r for r in rows if not r["real"]]), ("real runs", [r for r in rows if r["real"]])):
        print(f"== {label}: {len(sel)} clusters, flag undefined (< 3 aligned reads) in {sum(r['flag5'] is None for r in sel)}")
        for c in ("PRESENT", "ALLELE-LIKE", "DIVERGED", "NOVEL", "ELSEWHERE", "UNSUPPORTED"):
            x = [r for r in sel if r["cls"] == c and r["flag5"] is not None]
            if x:
                print(f"  {c:12s} n {len(x):4d}  5' flagged {sum(r['flag5'] for r in x) / len(x):.3f}  3' flagged {sum(r['flag3'] for r in x) / len(x):.3f}  either {sum(r['flag5'] or r['flag3'] for r in x) / len(x):.3f}  median flank 5'/3' {st.median(r['gap5'] for r in x)}/{st.median(r['gap3'] for r in x)}")
    ctl = [r for r in rows if not r["real"] and r["flag5"] is not None]
    changed = [r for r in ctl if (r["cls_global"] in GOOD) != (r["cls"] in GOOD)]
    fl = lambda r: r["flag5"] or r["flag3"]
    print(f"control clusters whose good/not-good status differs between the global and the local consensus: {len(changed)}; flagged {sum(map(fl, changed))}; "
          f"base rate of flagged among all control clusters {sum(map(fl, ctl)) / len(ctl):.3f}")
    reg = [r for r in ctl if r["cls_global"] in GOOD and r["cls"] not in GOOD]
    print("local-mode regressions (good -> not good):", [(r["run"], r["length"], r["reads"], r["cls"], dict(f5=r["flag5"], f3=r["flag3"], gap5=r["gap5"], gap3=r["gap3"])) for r in reg])
    print("the two Amendment-37 tail clusters:", [(r["run"], r["length"], r["reads"], r["cls"], dict(f5=r["flag5"], f3=r["flag3"], gap5=r["gap5"], gap3=r["gap3"])) for r in ctl if (r["run"], r["length"], r["reads"]) in (("ptr_control_s7", 2650, 5), ("ppy_control", 2170, 5))])
    big = [r for r in rows if r["real"] and r["cls"] in ("ELSEWHERE", "NOVEL") and r["reads"] >= 20]
    print("real-run candidates with >= 20 reads:", [(r["run"], r["reads"], r["length"], r["cls"][:4], r["flag5"], r["flag3"], r["gap5"], r["gap3"]) for r in sorted(big, key=lambda r: -r["reads"])])


if __name__ == "__main__":
    main()
