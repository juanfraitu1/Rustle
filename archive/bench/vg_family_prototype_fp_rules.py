#!/usr/bin/env python3
"""vg_family_prototype_fp_rules.py — grid-search empirical FP gates.

Loads the per-family feature TSV produced by vg_family_prototype_fp_characterize.py
and tests combinations of structural/graph features for separating false-positive
families (DNA-confirmed FP and/or E_p-impure) from real families with minimal
collateral damage.

Run: PYTHONHASHSEED=0 /home/juanfra/miniforge3/bin/python bench/vg_family_prototype_fp_rules.py
"""
import os
import sys

if os.environ.get("PYTHONHASHSEED") != "0":
    os.environ["PYTHONHASHSEED"] = "0"
    os.execv(sys.executable, [sys.executable] + sys.argv)

import csv
import json
from itertools import combinations, product
from collections import Counter

BENCH = os.path.dirname(os.path.abspath(__file__))
FEATURES_TSV = os.path.join(BENCH, "vg_family_prototype_fp_features.tsv")
OUT_JSON = os.path.join(BENCH, "vg_family_prototype_fp_rules.json")
OUT_TXT = os.path.join(BENCH, "vg_family_prototype_fp_rules.txt")


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
    return {
        "flagged": len(flagged),
        "removed_dna_fp": fp,
        "removed_ep_impure": ep,
        "collateral_real": real,
    }


def score(res, bad_total, real_total, real_penalty=1.0):
    """Higher is better. Penalize collateral real families."""
    bad_removed = res["removed_dna_fp"] + res["removed_ep_impure"]
    return bad_removed - real_penalty * res["collateral_real"]


def main():
    rows = load_features(FEATURES_TSV)
    total = len(rows)
    bad_rows = [r for r in rows if r["any_dna_fp"] or r["ep_impure"]]
    real_rows = [r for r in rows if not (r["any_dna_fp"] or r["ep_impure"])]
    bad_total = len(bad_rows)
    real_total = len(real_rows)

    print(f"[*] loaded {total} families: bad={bad_total} real={real_total}", flush=True)

    atomic = {
        "mean_pair_mult>=X":     [("mean_pair_mult", t) for t in [8, 10, 12, 15, 18, 20, 25, 30]],
        "max_pair_mult>=X":      [("max_pair_mult", t) for t in [15, 20, 25, 30, 35, 40, 50, 60, 80]],
        "pair_hub_frac>=X":      [("pair_hub_frac", t) for t in [0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.40, 0.50]],
        "n_distinct_loci>=X":    [("n_distinct_loci", t) for t in [1, 2, 3, 4, 5]],
        "n_protein_families>=X": [("n_protein_families", t) for t in [1, 2, 3]],
        "mean_n_exons>=X":       [("mean_n_exons", t) for t in [8, 10, 12, 14, 16, 18, 20]],
        "max_recip_overlap>=X":  [("max_recip_overlap", t) for t in [0.5, 0.7, 0.8, 0.9, 0.95, 0.99]],
        "mean_recip_overlap>=X": [("mean_recip_overlap", t) for t in [0.2, 0.3, 0.4, 0.5]],
        "n_members>=X":          [("n_members", t) for t in [30, 40, 50, 60, 80, 100]],
        "n_chrom>=X":            [("n_chrom", t) for t in [2, 3, 4, 5]],
    }

    def make_pred(name, threshold):
        if name == "strand_majority<1.0 AND n_chrom==1":
            return lambda r: r["strand_majority"] < 1.0 and r["n_chrom"] == 1
        return lambda r: r[name] >= threshold

    # single-feature rules
    single_results = []
    for label, pairs in atomic.items():
        for name, thresh in pairs:
            pred = make_pred(name, thresh)
            res = apply_rule(rows, pred)
            single_results.append({
                "rule": f"{name}>={thresh}",
                **res,
                "score": score(res, bad_total, real_total),
            })

    # add the strand architecture rule
    res = apply_rule(rows, make_pred("strand_majority<1.0 AND n_chrom==1", None))
    single_results.append({"rule": "strand_majority<1.0 AND n_chrom==1", **res, "score": score(res, bad_total, real_total)})

    # two-feature AND rules (conjunctive) — keep search tractable by using top single features
    top_single = sorted(single_results, key=lambda x: x["score"], reverse=True)[:20]
    top_atoms = []
    for r in top_single:
        rule = r["rule"]
        if "AND" in rule:
            continue
        name, thresh = rule.rsplit(">=", 1)
        top_atoms.append((name, float(thresh) if "." in thresh else int(thresh)))

    two_results = []
    seen = set()
    for (n1, t1), (n2, t2) in combinations(top_atoms, 2):
        key = tuple(sorted([(n1, t1), (n2, t2)]))
        if key in seen:
            continue
        seen.add(key)
        pred = lambda r, n1=n1, t1=t1, n2=n2, t2=t2: r[n1] >= t1 and r[n2] >= t2
        res = apply_rule(rows, pred)
        two_results.append({"rule": f"{n1}>={t1} AND {n2}>={t2}", **res, "score": score(res, bad_total, real_total)})

    # three-feature AND rules from top two-feature rules
    top_two = sorted(two_results, key=lambda x: x["score"], reverse=True)[:15]
    three_results = []
    seen3 = set()
    for r2a, r2b in combinations(top_two, 2):
        atoms_a = set(r2a["rule"].split(" AND "))
        atoms_b = set(r2b["rule"].split(" AND "))
        merged = sorted(atoms_a | atoms_b)
        if len(merged) != 3:
            continue
        key = tuple(merged)
        if key in seen3:
            continue
        seen3.add(key)
        # build predicate
        def _pred(r, merged=merged):
            for atom in merged:
                name, thresh = atom.rsplit(">=", 1)
                thresh = float(thresh) if "." in thresh else int(thresh)
                if r[name] < thresh:
                    return False
            return True
        res = apply_rule(rows, _pred)
        three_results.append({"rule": " AND ".join(merged), **res, "score": score(res, bad_total, real_total)})

    all_results = single_results + two_results + three_results

    # rank by score, also by precision-at-recall bins
    def precision(res):
        flagged = res["flagged"]
        return (res["removed_dna_fp"] + res["removed_ep_impure"]) / flagged if flagged else 0.0

    def recall_bad(res):
        return (res["removed_dna_fp"] + res["removed_ep_impure"]) / bad_total if bad_total else 0.0

    def real_recall(res):
        return 1 - res["collateral_real"] / real_total if real_total else 1.0

    for r in all_results:
        r["precision"] = round(precision(r), 4)
        r["recall_bad"] = round(recall_bad(r), 4)
        r["real_recall"] = round(real_recall(r), 4)
        r["f1_like"] = round(2 * precision(r) * recall_bad(r) / (precision(r) + recall_bad(r)), 4) if (precision(r) + recall_bad(r)) > 0 else 0.0

    all_results.sort(key=lambda x: x["score"], reverse=True)

    with open(OUT_JSON, "w") as fh:
        json.dump({
            "total": total,
            "bad_total": bad_total,
            "real_total": real_total,
            "top_rules": all_results[:100],
        }, fh, indent=2, sort_keys=False)
    print(f"[write] {OUT_JSON}", flush=True)

    # human-readable summary
    lines = []
    lines.append(f"Total families: {total}")
    lines.append(f"Bad families (DNA-FP OR E_p-impure): {bad_total}")
    lines.append(f"Real families: {real_total}")
    lines.append("")
    lines.append("=== Top 30 rules by score (bad_removed - collateral_real) ===")
    lines.append(f"{'rank':>4} {'rule':<55} {'flag':>4} {'bad':>4} {'FP':>3} {'EP':>3} {'real$':>4} {'prec':>6} {'rec':>6} {'score':>6}")
    for i, r in enumerate(all_results[:30], 1):
        bad = r["removed_dna_fp"] + r["removed_ep_impure"]
        lines.append(f"{i:>4} {r['rule']:<55} {r['flagged']:>4} {bad:>4} {r['removed_dna_fp']:>3} {r['removed_ep_impure']:>3} {r['collateral_real']:>4} {r['precision']:>6.3f} {r['recall_bad']:>6.3f} {r['score']:>6.1f}")

    # best precision at given recall levels
    lines.append("")
    lines.append("=== Best precision at fixed bad-recall thresholds ===")
    for target_recall in [0.10, 0.20, 0.30, 0.40, 0.50]:
        eligible = [r for r in all_results if r["recall_bad"] >= target_recall]
        if eligible:
            best = max(eligible, key=lambda x: x["precision"])
            lines.append(f"recall>={target_recall:.2f}: {best['rule']}  prec={best['precision']:.3f}  flagged={best['flagged']}  real_lost={best['collateral_real']}")
        else:
            lines.append(f"recall>={target_recall:.2f}: no rule reaches this recall")

    # best rules with zero collateral
    lines.append("")
    lines.append("=== Zero-collateral rules that catch the most bad families ===")
    zero_coll = sorted([r for r in all_results if r["collateral_real"] == 0], key=lambda x: x["score"], reverse=True)[:15]
    for r in zero_coll:
        bad = r["removed_dna_fp"] + r["removed_ep_impure"]
        lines.append(f"{r['rule']:<55} bad={bad:>3} FP={r['removed_dna_fp']} EP={r['removed_ep_impure']} flagged={r['flagged']:>3}")

    txt = "\n".join(lines)
    print(txt)
    with open(OUT_TXT, "w") as fh:
        fh.write(txt + "\n")
    print(f"[write] {OUT_TXT}", flush=True)


if __name__ == "__main__":
    main()
