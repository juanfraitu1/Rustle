#!/usr/bin/env python3
"""vg_family_prototype_fp_diagnose.py — diagnose which method-computable features
separate FP subtypes from real families.
"""
import os
import sys

if os.environ.get("PYTHONHASHSEED") != "0":
    os.environ["PYTHONHASHSEED"] = "0"
    os.execv(sys.executable, [sys.executable] + sys.argv)

import csv
import statistics
from collections import defaultdict

BENCH = os.path.dirname(os.path.abspath(__file__))
FEATURES_TSV = os.path.join(BENCH, "vg_family_prototype_fp_features.tsv")


def load_features(path):
    rows = []
    with open(path) as fh:
        r = csv.DictReader(fh, delimiter="\t")
        for row in r:
            rows.append({k: (float(v) if "." in v else int(v)) for k, v in row.items()})
    return rows


def main():
    rows = load_features(FEATURES_TSV)

    groups = {
        "real": [r for r in rows if not (r["any_dna_fp"] or r["ep_impure"])],
        "ep_impure": [r for r in rows if r["ep_impure"]],
        "dna_fp": [r for r in rows if r["any_dna_fp"]],
        "fp_multifam": [r for r in rows if r["fp_multifam"]],
        "fp_oversize": [r for r in rows if r["fp_oversize"]],
        "fp_allele": [r for r in rows if r["fp_allele"]],
    }

    # method-computable features
    features = [
        "n_members", "n_chrom", "same_chrom_pairs",
        "mean_recip_overlap", "max_recip_overlap", "strand_majority",
        "mean_n_exons", "mean_pair_mult", "median_pair_mult", "max_pair_mult",
        "pair_hub_frac", "mean_repeat_frac",
    ]

    print("Feature distributions by group (mean ± std):")
    print(f"{'feature':<22}", end="")
    for g in groups:
        print(f"{g:>12}", end="")
    print()
    print("-" * (22 + 12 * len(groups)))
    for feat in features:
        print(f"{feat:<22}", end="")
        for g, rs in groups.items():
            vals = [r[feat] for r in rs]
            if vals:
                mu = statistics.mean(vals)
                sd = statistics.stdev(vals) if len(vals) > 1 else 0.0
                print(f"{mu:>6.2f}±{sd:<4.2f}", end="")
            else:
                print(f"{'—':>12}", end="")
        print()

    print("\n\nThreshold-separation tables:")
    for feat in features:
        print(f"\n{feat}:")
        print(f"{'threshold':>12} {'flagged':>8} {'real_lost':>10} {'ep_lost':>8} {'fp_lost':>8} {'precision':>10}")
        vals = sorted(set(r[feat] for r in rows))
        # sample percentiles
        for thresh in [statistics.quantiles(vals, n=20)[i] for i in [3, 5, 7, 10, 13, 15, 17]] + [max(vals)]:
            flagged = [r for r in rows if r[feat] >= thresh]
            if not flagged:
                continue
            real_lost = sum(1 for r in flagged if r in groups["real"])
            ep_lost = sum(1 for r in flagged if r in groups["ep_impure"])
            fp_lost = sum(1 for r in flagged if r in groups["dna_fp"])
            prec = (ep_lost + fp_lost) / len(flagged)
            print(f"{thresh:>12.4f} {len(flagged):>8} {real_lost:>10} {ep_lost:>8} {fp_lost:>8} {prec:>10.3f}")


if __name__ == "__main__":
    main()
