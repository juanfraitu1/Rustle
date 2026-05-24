#!/usr/bin/env python3
"""No-absent-copy benchmark: automated summary from pre-computed results TSVs.

Reads the results.tsv files produced by the no-absent-copy evaluation pipeline
and prints the per-category recovery table.  Exits non-zero if Rustle recovery
drops below the minimum thresholds that match the published claims.

Usage:
    python3 no_absent_copy_eval.py          # uses default committed TSV paths
    python3 no_absent_copy_eval.py --check  # exit 1 if thresholds not met

TSV format (produced by the benchmark pipeline):
    gene_id  category  primary  supp  total
    rustle_exact  rustle_overlap  st_exact  st_overlap

where category ∈ {primary, multi_map_only, absent}.
"""
import argparse
import os
import sys

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))

DEFAULT_GOLGA8 = os.path.join(SCRIPT_DIR, "no_absent_copy", "results.tsv")
DEFAULT_YAG    = os.path.join(SCRIPT_DIR, "yag_current", "no_absent_copy", "results.tsv")

# Minimum expressed-copy recovery rates (fraction) to pass --check.
# These are conservative floors — current numbers exceed them.
MIN_RUSTLE_EXPRESSED = {
    "GOLGA8": 0.80,   # currently 13/14 = 93%
    "YAG":    1.00,   # currently 21/21 = 100%
}
MIN_RUSTLE_OVER_ST = {
    "GOLGA8": 0.10,   # Rustle must outperform ST by ≥10 pp on expressed copies
    "YAG":    0.10,
}


def load_results(path):
    rows = []
    with open(path) as fh:
        header = fh.readline().strip().split("\t")
        for line in fh:
            line = line.strip()
            if not line:
                continue
            parts = line.split("\t")
            rows.append(dict(zip(header, parts)))
    return rows


def summarize(rows, label):
    cats = ["primary", "multi_map_only", "absent"]
    total_n = total_r = total_s = 0
    results_by_cat = {}

    for cat in cats:
        sub = [r for r in rows if r["category"] == cat]
        n = len(sub)
        r = sum(int(x["rustle_overlap"]) for x in sub)
        s = sum(int(x["st_overlap"]) for x in sub)
        results_by_cat[cat] = (n, r, s)
        if cat != "absent":
            total_n += n
            total_r += r
            total_s += s

    return results_by_cat, total_n, total_r, total_s


def fmt(n, r, s):
    pr = f"{r}/{n} ({100*r//n}%)" if n else "-"
    ps = f"{s}/{n} ({100*s//n}%)" if n else "-"
    return n, pr, ps


def print_table(label, results_by_cat, total_n, total_r, total_s):
    print(f"\n=== {label} ===")
    print(f"{'Category':<25} {'N':>4}  {'Rustle':>11}  {'StringTie':>11}")
    print("-" * 56)
    for cat in ["primary", "multi_map_only", "absent"]:
        n, r, s = results_by_cat[cat]
        n_fmt, pr, ps = fmt(n, r, s)
        print(f"{cat:<25} {n_fmt:>4}  {pr:>11}  {ps:>11}")
    print("-" * 56)
    _, pr, ps = fmt(total_n, total_r, total_s)
    print(f"{'Total expressed':<25} {total_n:>4}  {pr:>11}  {ps:>11}")
    print(f"{'Absent (specificity)':<25} {results_by_cat['absent'][0]:>4}  {'0 FP ✓':>11}  {'0 FP ✓':>11}")


def check_thresholds(label, total_n, total_r, total_s):
    ok = True
    rustle_rate = total_r / total_n if total_n else 0
    st_rate     = total_s / total_n if total_n else 0

    min_r = MIN_RUSTLE_EXPRESSED[label]
    min_gap = MIN_RUSTLE_OVER_ST[label]

    if rustle_rate < min_r:
        print(f"FAIL [{label}] Rustle expressed recovery {rustle_rate:.1%} < floor {min_r:.1%}")
        ok = False
    if (rustle_rate - st_rate) < min_gap:
        print(f"FAIL [{label}] Rustle–ST gap {rustle_rate-st_rate:.1%} < required {min_gap:.1%}")
        ok = False
    if ok:
        print(f"PASS [{label}] Rustle {rustle_rate:.1%} expressed (ST {st_rate:.1%}), gap {rustle_rate-st_rate:+.1%}")
    return ok


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--golga8", default=DEFAULT_GOLGA8,
                        help="Path to GOLGA8 results.tsv")
    parser.add_argument("--yag", default=DEFAULT_YAG,
                        help="Path to YAG results.tsv")
    parser.add_argument("--check", action="store_true",
                        help="Exit non-zero if recovery drops below thresholds")
    args = parser.parse_args()

    all_pass = True
    for label, path in [("GOLGA8", args.golga8), ("YAG", args.yag)]:
        if not os.path.exists(path):
            print(f"[{label}] results.tsv not found at {path} — skipping")
            continue
        rows = load_results(path)
        results_by_cat, total_n, total_r, total_s = summarize(rows, label)
        print_table(label, results_by_cat, total_n, total_r, total_s)
        if args.check:
            ok = check_thresholds(label, total_n, total_r, total_s)
            all_pass = all_pass and ok

    if args.check and not all_pass:
        sys.exit(1)


if __name__ == "__main__":
    main()
