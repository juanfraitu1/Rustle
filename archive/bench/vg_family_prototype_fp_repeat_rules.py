#!/usr/bin/env python3
"""vg_family_prototype_fp_repeat_rules.py — grid-search repeat-aware FP gates.

Loads the enriched repeat-feature table and tests many combinations of
repeat/TE signals with VG topology.
"""
import os
import sys

if os.environ.get("PYTHONHASHSEED") != "0":
    os.environ["PYTHONHASHSEED"] = "0"
    os.execv(sys.executable, [sys.executable] + sys.argv)

import csv
import json
from itertools import product, combinations

BENCH = os.path.dirname(os.path.abspath(__file__))
FEATURES_TSV = os.path.join(BENCH, "vg_family_prototype_fp_repeat_features.tsv")
OUT_JSON = os.path.join(BENCH, "vg_family_prototype_fp_repeat_rules.json")
OUT_TXT = os.path.join(BENCH, "vg_family_prototype_fp_repeat_rules.txt")


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


def score(res):
    return (res["removed_dna_fp"] + res["removed_ep_impure"]) - res["collateral_real"]


def main():
    rows = load_features(FEATURES_TSV)
    bad_total = sum(1 for r in rows if r["any_dna_fp"] or r["ep_impure"])

    # features to combine
    repeat_feats = ["mean_repeat_frac", "mean_exon_repeat_frac", "max_exon_repeat_frac",
                    "mean_node_repeat_frac", "max_node_repeat_frac", "repeat_hub_frac",
                    "repeat_rich_frac_0.3", "repeat_rich_frac_0.5",
                    "repeat_weighted_mult", "high_repeat_hub_proxy"]
    topo_feats = ["n_members", "mean_pair_mult", "pair_hub_frac", "max_pair_mult"]

    thresholds = {
        "mean_repeat_frac": [0.20, 0.25, 0.30, 0.35, 0.40, 0.45, 0.50],
        "mean_exon_repeat_frac": [0.05, 0.08, 0.10, 0.12, 0.15, 0.20, 0.25, 0.30],
        "max_exon_repeat_frac": [0.10, 0.15, 0.20, 0.25, 0.30, 0.40, 0.50, 0.60],
        "mean_node_repeat_frac": [0.05, 0.08, 0.10, 0.12, 0.15, 0.20, 0.25, 0.30],
        "max_node_repeat_frac": [0.10, 0.15, 0.20, 0.25, 0.30, 0.40, 0.50, 0.60],
        "repeat_hub_frac": [0.02, 0.05, 0.10, 0.15, 0.20, 0.25, 0.30],
        "repeat_rich_frac_0.3": [0.05, 0.10, 0.20, 0.30, 0.40, 0.50],
        "repeat_rich_frac_0.5": [0.05, 0.10, 0.20, 0.30, 0.40, 0.50],
        "repeat_weighted_mult": [1, 2, 3, 5, 8, 10, 15, 20],
        "high_repeat_hub_proxy": [0.02, 0.05, 0.10, 0.15, 0.20],
        "n_members": [30, 40, 50, 60, 70, 80],
        "mean_pair_mult": [10, 15, 20, 25, 30],
        "pair_hub_frac": [0.10, 0.20, 0.30, 0.40, 0.50],
        "max_pair_mult": [20, 30, 40, 50, 60],
    }

    all_results = []

    # single-feature repeat rules
    for feat in repeat_feats + topo_feats:
        for t in thresholds.get(feat, []):
            res = apply_rule(rows, lambda r, f=feat, t=t: r[f] >= t)
            all_results.append({"rule": f"{feat}>={t}", **res, "score": score(res)})

    # two-feature AND: one repeat + one topo
    for rfeat in repeat_feats:
        for t1 in thresholds.get(rfeat, []):
            for tfeat in topo_feats:
                for t2 in thresholds.get(tfeat, []):
                    rule = f"{rfeat}>={t1} AND {tfeat}>={t2}"
                    pred = lambda r, rf=rfeat, t1=t1, tf=tfeat, t2=t2: r[rf] >= t1 and r[tf] >= t2
                    res = apply_rule(rows, pred)
                    all_results.append({"rule": rule, **res, "score": score(res)})

    # a few three-feature combos from top two-feature
    top = sorted(all_results, key=lambda x: x["score"], reverse=True)[:50]
    atoms = set()
    for r in top:
        for atom in r["rule"].split(" AND "):
            atoms.add(atom)
    atoms = sorted(atoms)
    seen3 = set()
    for a, b, c in combinations(atoms, 3):
        n1, t1 = a.rsplit(">=", 1); n2, t2 = b.rsplit(">=", 1); n3, t3 = c.rsplit(">=", 1)
        if len({n1, n2, n3}) != 3:
            continue
        key = tuple(sorted([a, b, c]))
        if key in seen3:
            continue
        seen3.add(key)
        t1 = float(t1) if "." in t1 else int(t1)
        t2 = float(t2) if "." in t2 else int(t2)
        t3 = float(t3) if "." in t3 else int(t3)
        pred = lambda r, n1=n1, t1=t1, n2=n2, t2=t2, n3=n3, t3=t3: r[n1] >= t1 and r[n2] >= t2 and r[n3] >= t3
        res = apply_rule(rows, pred)
        all_results.append({"rule": " AND ".join(key), **res, "score": score(res)})

    def precision(r):
        flagged = r["flagged"]
        return (r["removed_dna_fp"] + r["removed_ep_impure"]) / flagged if flagged else 0.0

    def recall_bad(r):
        return (r["removed_dna_fp"] + r["removed_ep_impure"]) / bad_total if bad_total else 0.0

    for r in all_results:
        r["precision"] = round(precision(r), 4)
        r["recall_bad"] = round(recall_bad(r), 4)

    all_results.sort(key=lambda x: x["score"], reverse=True)

    with open(OUT_JSON, "w") as fh:
        json.dump({"total": len(rows), "bad_total": bad_total,
                   "top_rules": all_results[:200]}, fh, indent=2)
    print(f"[write] {OUT_JSON}", flush=True)

    lines = []
    lines.append(f"Total families: {len(rows)}  bad: {bad_total}")
    lines.append("")
    lines.append("=== Top 40 repeat-aware rules by score ===")
    lines.append(f"{'rank':>4} {'rule':<60} {'flag':>5} {'bad':>5} {'FP':>4} {'EP':>4} {'real$':>6} {'prec':>6} {'rec':>6} {'score':>6}")
    for i, r in enumerate(all_results[:40], 1):
        bad = r["removed_dna_fp"] + r["removed_ep_impure"]
        lines.append(f"{i:>4} {r['rule']:<60} {r['flagged']:>5} {bad:>5} {r['removed_dna_fp']:>4} {r['removed_ep_impure']:>4} {r['collateral_real']:>6} {r['precision']:>6.3f} {r['recall_bad']:>6.3f} {r['score']:>6.1f}")

    lines.append("")
    lines.append("=== Zero-collateral rules catching the most bad families ===")
    zero = sorted([r for r in all_results if r["collateral_real"] == 0], key=lambda x: x["score"], reverse=True)[:20]
    for r in zero:
        bad = r["removed_dna_fp"] + r["removed_ep_impure"]
        lines.append(f"{r['rule']:<60} bad={bad:>3} FP={r['removed_dna_fp']} EP={r['removed_ep_impure']}")

    txt = "\n".join(lines)
    print(txt)
    with open(OUT_TXT, "w") as fh:
        fh.write(txt + "\n")
    print(f"[write] {OUT_TXT}", flush=True)


if __name__ == "__main__":
    main()
