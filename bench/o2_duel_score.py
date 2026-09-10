#!/usr/bin/env python3
"""Score PREREG_best_by_duel (819c1615): --best-by-duel run vs the flag-off base, over the contested set.
  python3 bench/o2_duel_score.py --base ours_final.assignments.tsv --new ours_duel.assignments.tsv
"""
import argparse, csv, statistics
from collections import Counter


def load(p):
    with open(p) as fh:
        return {r["read_name"]: r for r in csv.DictReader(fh, delimiter="\t")}


def med(xs):
    return statistics.median(xs) if xs else float("nan")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--base", required=True)
    ap.add_argument("--new", required=True)
    a = ap.parse_args()
    b, n = load(a.base), load(a.new)
    cb = {k for k, r in b.items() if r["origin_rejected"] == "0" and int(r["n_candidates"]) >= 2}
    cn = {k for k, r in n.items() if r["origin_rejected"] == "0" and int(r["n_candidates"]) >= 2}
    print(f"contested base {len(cb)} / new {len(cn)} / symmetric diff {len(cb ^ cn)}; status base {dict(Counter(b[k]['status'] for k in cb))} new {dict(Counter(n[k]['status'] for k in cn))}")
    asg = [k for k in cb if b[k]["status"] == "assigned"]
    ch = [k for k in asg if n[k]["status"] != "assigned" or n[k]["catalog_copy_idx"] != b[k]["catalog_copy_idx"]]
    print(f"P1 base assigned {len(asg)}: changed status/copy = {len(ch)}  (pass iff 0)")
    for k in ch[:5]:
        print(f"   {k}: {b[k]['status']} {b[k]['catalog_copy_idx']} -> {n[k]['status']} {n[k]['catalog_copy_idx']}")
    neg = [k for k in cb if float(b[k]["margin"]) < 0]
    newa = [k for k in neg if n[k]["status"] == "assigned"]
    pos_newa = [k for k in cb if float(b[k]["margin"]) >= 0 and b[k]["status"] != "assigned" and n[k]["status"] == "assigned"]
    print(f"P2 rows with margin<0: {len(neg)}; now assigned: {len(newa)} (pass 10..65); new assignments from margin>=0 rows: {len(pos_newa)} (pass iff 0)")
    print(f"   margin<0 rows now: {dict(Counter(n[k]['status'] for k in neg))}; bk changed in {sum(1 for k in neg if n[k]['catalog_copy_idx'] != b[k]['catalog_copy_idx'])}")
    pos = [k for k in cb if float(b[k]["margin"]) >= 0]
    kept = sum(1 for k in pos if n[k]["status"] == b[k]["status"] and n[k]["catalog_copy_idx"] == b[k]["catalog_copy_idx"])
    print(f"P3 rows with margin>=0 keeping bk and status: {kept}/{len(pos)} = {100*kept/max(1,len(pos)):.2f}%  (pass >=99.5)")
    if newa:
        print(f"P4 new assignments: n_decisive med {med(int(n[k]['n_decisive']) for k in newa)}, margin med {med(float(n[k]['margin']) for k in newa):.1f}; copies {Counter(n[k]['catalog_copy_idx'] for k in newa).most_common()}")
        print(f"   base bk -> new bk pairs: {Counter((b[k]['catalog_copy_idx'], n[k]['catalog_copy_idx']) for k in newa).most_common(8)}")
    other = [k for k in b if k not in cb and (k not in n or n[k]["status"] != b[k]["status"])]
    print(f"sanity: non-contested rows with a changed status: {len(other)}")


if __name__ == "__main__":
    main()
