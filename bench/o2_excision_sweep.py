#!/usr/bin/env python3
"""PREREG_excision_sweep (4c832450): for every copy k, follow the molecules `--base` assigned to k into the
run without k (`<dir>/no{k}.assignments.tsv`, remap `<dir>/remap_no{k}.tsv` new->old).

  python3 bench/o2_excision_sweep.py --base ours_final.assignments.tsv --dir excise_all --copies copies16.tsv \
      --paf human_gspans.paf --n 26
"""
import argparse, csv, re, statistics
from collections import Counter, defaultdict


def contested(rows):
    return {n: r for n, r in rows.items() if r["origin_rejected"] == "0" and int(r["n_candidates"]) >= 2}


def load(path):
    with open(path) as fh:
        return {r["read_name"]: r for r in csv.DictReader(fh, delimiter="\t")}


def identity_matrix(paf, copies):
    """1 - X/aligned over the genomic unit spans, all fragments summed, keyed by (copy_idx, copy_idx)."""
    span2idx = {}
    with open(copies) as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            span2idx[f"{r['chrom']}:{int(r['start'])+1}-{r['end']}"] = r["copy_idx"]
    acc = defaultdict(lambda: [0, 0])
    for l in open(paf):
        f = l.rstrip().split("\t")
        if f[0] == f[5] or f[0] not in span2idx or f[5] not in span2idx:
            continue
        k = tuple(sorted((span2idx[f[0]], span2idx[f[5]])))
        cg = next(x for x in f[12:] if x.startswith("cg:Z:"))[5:]
        acc[k][0] += int(f[10])
        acc[k][1] += sum(int(n) for n, op in re.findall(r"(\d+)([=XID])", cg) if op == "X")
    return {k: (1 - x / a if a else 0.0) for k, (a, x) in acc.items()}


def med(xs):
    return statistics.median(xs) if xs else float("nan")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--base", required=True)
    ap.add_argument("--dir", required=True)
    ap.add_argument("--copies", required=True)
    ap.add_argument("--paf", required=True)
    ap.add_argument("--n", type=int, default=26)
    a = ap.parse_args()
    base = contested(load(a.base))
    ident = identity_matrix(a.paf, a.copies)
    by_copy = defaultdict(list)
    for n, r in base.items():
        if r["status"] == "assigned":
            by_copy[r["catalog_copy_idx"]].append(n)
    tot_abst = tot = 0
    moved_all = []
    p2_fail, p4_fail, p5_fail = [], [], []
    print(f"{'k':>2} {'assigned':>8} {'abstain':>8} {'moved':>5}  {'moved->(copy,identity,margin)':<40} {'others kept':>12} {'contested status kept':>22}")
    for k in range(a.n):
        kk = str(k)
        remap = {}
        with open(f"{a.dir}/remap_no{k}.tsv") as fh:
            for l in fh:
                new, old = l.split()
                remap[new] = old
        ex = load(f"{a.dir}/no{k}.assignments.tsv")
        for r in ex.values():
            r["old_copy"] = remap.get(r["catalog_copy_idx"], "?")
        mine = by_copy.get(kk, [])
        abst = [n for n in mine if ex.get(n, {}).get("status") != "assigned"]
        moved = [(n, ex[n]["old_copy"], float(ex[n]["margin"])) for n in mine if n in ex and ex[n]["status"] == "assigned"]
        tot += len(mine); tot_abst += len(abst)
        moved_all += [(kk, c, m, ident.get(tuple(sorted((kk, c))), float("nan"))) for _, c, m in moved]
        # P4: molecules assigned to other copies keep status + copy
        others = [n for c, ns in by_copy.items() if c != kk for n in ns]
        kept = sum(1 for n in others if n in ex and ex[n]["status"] == "assigned" and ex[n]["old_copy"] == base[n]["catalog_copy_idx"])
        # P5-style: contested statuses of molecules not assigned to k
        rest = [n for n in base if base[n].get("catalog_copy_idx") != kk or base[n]["status"] != "assigned"]
        kept_st = sum(1 for n in rest if n in ex and ex[n]["status"] == base[n]["status"])
        frac_o = kept / len(others) if others else 1.0
        frac_s = kept_st / len(rest) if rest else 1.0
        if len(mine) >= 5 and len(abst) / len(mine) < 0.8: p2_fail.append((kk, len(abst), len(mine)))
        if frac_o < 0.99: p4_fail.append((kk, kept, len(others)))
        if not mine and 1 - frac_s > 0.01: p5_fail.append((kk, len(rest) - kept_st, len(rest)))
        mv = ", ".join(f"({c},{ident.get(tuple(sorted((kk, c))), float('nan')):.3f},{m:.0f})" for _, c, m in moved[:4])
        print(f"{kk:>2} {len(mine):>8} {len(abst):>8} {len(moved):>5}  {mv:<40} {kept:>5}/{len(others):<6} {kept_st:>9}/{len(rest):<6} {100*frac_s:5.1f}%")
    print(f"\nP1 pooled abstention: {tot_abst}/{tot} = {100*tot_abst/max(1,tot):.1f}%  (pass >=95, refuted <90)")
    print(f"P2 copies with >=5 assigned below 80% abstention: {p2_fail or 'none'}  (refuted if any <70%)")
    if moved_all:
        low = [m for m in moved_all if not (m[3] >= 0.97)]
        conf = [m for m in moved_all if m[2] >= 40]
        print(f"P3 moved {len(moved_all)}: identity med {med([m[3] for m in moved_all]):.4f}, margin med {med([m[2] for m in moved_all]):.1f}; to a copy <0.97: {len(low)} ({100*len(low)/len(moved_all):.0f}%); silent-confident (margin>=40): {len(conf)}  (refuted: >30% low-identity or >=3 silent-confident)")
        print("   moved detail (from, to, margin, identity):", sorted(Counter((m[0], m[1]) for m in moved_all).items()))
    else:
        print("P3 moved 0")
    print(f"P4 runs where others' assignments kept <99%: {p4_fail or 'none'}  (refuted if any <97%)")
    print(f"P5 zero-assignment copies changing >1% of contested statuses: {p5_fail or 'none'}  (refuted if any >3%)")


if __name__ == "__main__":
    main()
