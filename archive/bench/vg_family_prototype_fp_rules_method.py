#!/usr/bin/env python3
"""vg_family_prototype_fp_rules_method.py — grid-search empirical FP gates using
ONLY features computable from the method itself (no annotation-derived labels).

Method-computable features used here:
  n_members, n_chrom, same_chrom_pairs, mean_recip_overlap, max_recip_overlap,
  strand_majority, mean_n_exons, mean_pair_mult, median_pair_mult, max_pair_mult,
  pair_hub_frac, mean_repeat_frac
"""
import os
import sys

if os.environ.get("PYTHONHASHSEED") != "0":
    os.environ["PYTHONHASHSEED"] = "0"
    os.execv(sys.executable, [sys.executable] + sys.argv)

import csv
import json
from itertools import combinations

BENCH = os.path.dirname(os.path.abspath(__file__))
FEATURES_TSV = os.path.join(BENCH, "vg_family_prototype_fp_features.tsv")
OUT_JSON = os.path.join(BENCH, "vg_family_prototype_fp_rules_method.json")
OUT_TXT = os.path.join(BENCH, "vg_family_prototype_fp_rules_method.txt")

METHOD_FEATURES = [
    "n_members", "n_chrom", "same_chrom_pairs",
    "mean_recip_overlap", "max_recip_overlap", "strand_majority",
    "mean_n_exons", "mean_pair_mult", "median_pair_mult", "max_pair_mult",
    "pair_hub_frac", "mean_repeat_frac",
]


def load_features(path):
    rows = []
    with open(path) as fh:
        r = csv.DictReader(fh, delimiter="\t")
        for row in r:
            rows.append({k: (float(v) if "." in v else int(v)) for k, v in row.items()})
    return rows


def apply_rule(rows, pred):
    flagged = [r for r in rows if pred(r)]
    fp = sum(1 for r in flagged if r["any_dna_fp"])
    ep = sum(1 for r in flagged if r["ep_impure"])
    real = len(flagged) - fp - ep
    return {"flagged": len(flagged), "removed_dna_fp": fp,
            "removed_ep_impure": ep, "collateral_real": real}


def score(res, real_penalty=1.0):
    bad = res["removed_dna_fp"] + res["removed_ep_impure"]
    return bad - real_penalty * res["collateral_real"]


def make_pred(name, thresh):
    if name == "strand_majority<1.0_AND_n_chrom==1":
        return lambda r: r["strand_majority"] < 1.0 and r["n_chrom"] == 1
    return lambda r: r[name] >= thresh


def main():
    rows = load_features(FEATURES_TSV)
    bad_total = sum(1 for r in rows if r["any_dna_fp"] or r["ep_impure"])
    real_total = len(rows) - bad_total

    # threshold grids per feature
    grids = {
        "n_members": [30, 40, 50, 60, 70, 80, 100],
        "n_chrom": [2, 3, 4, 5],
        "same_chrom_pairs": [200, 400, 600, 800, 1000, 1200, 1400],
        "mean_recip_overlap": [0.3, 0.4, 0.5, 0.6, 0.7],
        "max_recip_overlap": [0.8, 0.9, 0.95, 0.99, 1.0],
        "strand_majority": [0.5, 0.6, 0.7, 0.8, 0.9],
        "mean_n_exons": [8, 10, 12, 14, 16, 18, 20],
        "mean_pair_mult": [8, 10, 12, 15, 18, 20, 25, 30],
        "median_pair_mult": [10, 15, 20, 25, 30, 40],
        "max_pair_mult": [20, 25, 30, 40, 50, 60, 80],
        "pair_hub_frac": [0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.40, 0.50],
        "mean_repeat_frac": [0.3, 0.35, 0.4, 0.45, 0.5, 0.55],
    }

    all_results = []

    # single rules
    for feat in METHOD_FEATURES:
        for t in grids.get(feat, []):
            res = apply_rule(rows, make_pred(feat, t))
            all_results.append({"rule": f"{feat}>={t}", **res,
                                "score": score(res)})

    # architecture rule
    res = apply_rule(rows, make_pred("strand_majority<1.0_AND_n_chrom==1", None))
    all_results.append({"rule": "strand_majority<1.0 AND n_chrom==1", **res, "score": score(res)})

    # pick top singles, then build 2- and 3-feature ANDs with distinct atoms
    top_single = sorted(all_results, key=lambda x: x["score"], reverse=True)[:25]
    atoms = []
    seen_atoms = set()
    for r in top_single:
        if "AND" in r["rule"]:
            continue
        name, thresh = r["rule"].rsplit(">=", 1)
        if name not in METHOD_FEATURES:
            continue
        thresh = float(thresh) if "." in thresh else int(thresh)
        key = (name, thresh)
        if key not in seen_atoms:
            seen_atoms.add(key)
            atoms.append(key)

    two_results = []
    seen2 = set()
    for (n1, t1), (n2, t2) in combinations(atoms, 2):
        if n1 == n2:
            continue
        key = tuple(sorted([(n1, t1), (n2, t2)]))
        if key in seen2:
            continue
        seen2.add(key)
        pred = lambda r, n1=n1, t1=t1, n2=n2, t2=t2: r[n1] >= t1 and r[n2] >= t2
        res = apply_rule(rows, pred)
        two_results.append({"rule": f"{n1}>={t1} AND {n2}>={t2}", **res, "score": score(res)})

    top_two = sorted(two_results, key=lambda x: x["score"], reverse=True)[:20]
    three_results = []
    seen3 = set()
    for r2a, r2b in combinations(top_two, 2):
        atoms_a = [(a.rsplit(">=", 1)[0], a.rsplit(">=", 1)[1]) for a in r2a["rule"].split(" AND ")]
        atoms_b = [(a.rsplit(">=", 1)[0], a.rsplit(">=", 1)[1]) for a in r2b["rule"].split(" AND ")]
        merged = set(atoms_a) | set(atoms_b)
        if len(merged) != 3:
            continue
        key = tuple(sorted(merged))
        if key in seen3:
            continue
        seen3.add(key)
        def _pred(r, merged=merged):
            for name, thresh in merged:
                thresh = float(thresh) if "." in thresh else int(thresh)
                if r[name] < thresh:
                    return False
            return True
        res = apply_rule(rows, _pred)
        rule_str = " AND ".join(f"{n}>={t}" for n, t in sorted(merged))
        three_results.append({"rule": rule_str, **res, "score": score(res)})

    all_results += two_results + three_results

    def precision(r):
        flagged = r["flagged"]
        return (r["removed_dna_fp"] + r["removed_ep_impure"]) / flagged if flagged else 0.0

    def recall_bad(r):
        return (r["removed_dna_fp"] + r["removed_ep_impure"]) / bad_total if bad_total else 0.0

    def real_recall(r):
        return 1 - r["collateral_real"] / real_total if real_total else 1.0

    for r in all_results:
        r["precision"] = round(precision(r), 4)
        r["recall_bad"] = round(recall_bad(r), 4)
        r["real_recall"] = round(real_recall(r), 4)
        r["f1_like"] = round(2 * precision(r) * recall_bad(r) / (precision(r) + recall_bad(r)), 4) if (precision(r) + recall_bad(r)) > 0 else 0.0

    all_results.sort(key=lambda x: x["score"], reverse=True)

    with open(OUT_JSON, "w") as fh:
        json.dump({
            "total": len(rows),
            "bad_total": bad_total,
            "real_total": real_total,
            "method_features": METHOD_FEATURES,
            "top_rules": all_results[:100],
        }, fh, indent=2)
    print(f"[write] {OUT_JSON}", flush=True)

    lines = []
    lines.append(f"Total families: {len(rows)}  bad={bad_total}  real={real_total}")
    lines.append("")
    lines.append("=== Top 30 method-computable rules by score ===")
    lines.append(f"{'rank':>4} {'rule':<55} {'flag':>4} {'bad':>4} {'FP':>3} {'EP':>3} {'real$':>4} {'prec':>6} {'rec':>6} {'score':>6}")
    for i, r in enumerate(all_results[:30], 1):
        bad = r["removed_dna_fp"] + r["removed_ep_impure"]
        lines.append(f"{i:>4} {r['rule']:<55} {r['flagged']:>4} {bad:>4} {r['removed_dna_fp']:>3} {r['removed_ep_impure']:>3} {r['collateral_real']:>4} {r['precision']:>6.3f} {r['recall_bad']:>6.3f} {r['score']:>6.1f}")

    lines.append("")
    lines.append("=== Best precision at fixed bad-recall thresholds ===")
    for target in [0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.40, 0.50]:
        elig = [r for r in all_results if r["recall_bad"] >= target]
        if elig:
            best = max(elig, key=lambda x: x["precision"])
            lines.append(f"recall>={target:.2f}: {best['rule']}  prec={best['precision']:.3f}  flagged={best['flagged']}  real_lost={best['collateral_real']}")
        else:
            lines.append(f"recall>={target:.2f}: no rule reaches this recall")

    lines.append("")
    lines.append("=== Zero-collateral rules that catch the most bad families ===")
    zero = sorted([r for r in all_results if r["collateral_real"] == 0], key=lambda x: x["score"], reverse=True)[:20]
    for r in zero:
        bad = r["removed_dna_fp"] + r["removed_ep_impure"]
        lines.append(f"{r['rule']:<55} bad={bad:>3} FP={r['removed_dna_fp']} EP={r['removed_ep_impure']} flagged={r['flagged']:>3}")

    txt = "\n".join(lines)
    print(txt)
    with open(OUT_TXT, "w") as fh:
        fh.write(txt + "\n")
    print(f"[write] {OUT_TXT}", flush=True)


if __name__ == "__main__":
    main()
