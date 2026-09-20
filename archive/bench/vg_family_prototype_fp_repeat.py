#!/usr/bin/env python3
"""vg_family_prototype_fp_repeat.py — add repeat/TE-aware features and gates.

Loads the existing feature table, computes exon-only soft-mask fractions from
GGO.fasta ( RepeatMasker soft-masked ), and tests empirical gates that combine
repeat content with VG topology.

Run: PYTHONHASHSEED=0 /home/juanfra/miniforge3/bin/python bench/vg_family_prototype_fp_repeat.py
"""
import os
import sys

if os.environ.get("PYTHONHASHSEED") != "0":
    os.environ["PYTHONHASHSEED"] = "0"
    os.execv(sys.executable, [sys.executable] + sys.argv)

import csv
import json
import statistics
from collections import defaultdict

import pysam

BENCH = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, BENCH)

from vg_family_prototype import load_meta, load_skeletons, FASTA

FEATURES_TSV = os.path.join(BENCH, "vg_family_prototype_fp_features.tsv")
OUT_TSV = os.path.join(BENCH, "vg_family_prototype_fp_repeat_features.tsv")
OUT_JSON = os.path.join(BENCH, "vg_family_prototype_fp_repeat.json")
OUT_TXT = os.path.join(BENCH, "vg_family_prototype_fp_repeat.txt")


def exon_intervals(start, end, introns):
    """Return list of (exon_start, exon_end) within [start,end] given introns."""
    pts = [start]
    for a, b in introns:
        if start < a and b < end:
            pts.append(a)
            pts.append(b)
    pts.append(end)
    ivs = []
    for i in range(0, len(pts) - 1, 2):
        s, e = pts[i], pts[i + 1]
        if e > s:
            ivs.append((s, e))
    return ivs


def repeat_frac_seq(seq):
    if not seq:
        return 0.0
    return sum(1 for c in seq if c.islower()) / len(seq)


def compute_exon_repeat(fa, loc):
    """Mean soft-mask fraction over spliced exons of a locus."""
    ivs = exon_intervals(loc["start"], loc["end"], loc["introns"])
    if not ivs:
        return 0.0
    fracs = []
    for s, e in ivs:
        seq = fa.fetch(loc["chrom"], s, e)
        if seq:
            fracs.append(repeat_frac_seq(seq))
    return statistics.mean(fracs) if fracs else 0.0


def main():
    print("[*] loading existing features ...", flush=True)
    rows = {}
    with open(FEATURES_TSV) as fh:
        r = csv.DictReader(fh, delimiter="\t")
        for row in r:
            rows[int(row["family_id"])] = {k: (float(v) if "." in v else int(v)) for k, v in row.items()}

    print("[*] loading meta + skeletons + genome ...", flush=True)
    meta = load_meta()
    loci = load_skeletons(meta)
    fa = pysam.FastaFile(FASTA)

    # map lid -> loc
    loc_by_id = {loc["lid"]: loc for loc in loci}

    # load catalog to know members per family
    catalog = defaultdict(list)
    with open(os.path.join(BENCH, "vg_family_prototype.tsv")) as fh:
        r = csv.DictReader(fh, delimiter="\t")
        for row in r:
            catalog[int(row["fam_id"])].append(row["member"])

    print("[*] computing exon-only repeat fractions ...", flush=True)
    exon_repeat = {}
    for loc in loci:
        exon_repeat[loc["lid"]] = compute_exon_repeat(fa, loc)

    # enrich features
    enriched = []
    for fid, feat in rows.items():
        members = catalog.get(fid, [])
        span_fracs = [feat["mean_repeat_frac"]]  # existing whole-span mean
        exon_fracs = [exon_repeat.get(m, 0.0) for m in members]
        if not exon_fracs:
            exon_fracs = [0.0]

        new = dict(feat)
        new["mean_exon_repeat_frac"] = round(statistics.mean(exon_fracs), 4)
        new["max_exon_repeat_frac"] = round(max(exon_fracs), 4)
        new["median_exon_repeat_frac"] = round(statistics.median(exon_fracs), 4)
        new["repeat_rich_frac_0.5"] = round(sum(1 for x in exon_fracs if x >= 0.5) / len(exon_fracs), 4)
        new["repeat_rich_frac_0.3"] = round(sum(1 for x in exon_fracs if x >= 0.3) / len(exon_fracs), 4)

        # repeat-weighted pair multiplicity: average member exon repeat * mean_pair_mult
        new["repeat_weighted_mult"] = round(new["mean_exon_repeat_frac"] * feat["mean_pair_mult"], 4)

        # high-repeat hub: pair nodes whose multiplicity >= 30 and member exon repeat >= 0.5
        # (approximate: we don't recompute node-level repeat here)
        new["high_repeat_hub_proxy"] = round(
            feat["pair_hub_frac"] * new["repeat_rich_frac_0.5"], 4)

        enriched.append(new)

    # write enriched TSV
    cols = list(enriched[0].keys())
    with open(OUT_TSV, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=cols, delimiter="\t")
        w.writeheader()
        for row in enriched:
            w.writerow(row)
    print(f"[write] {OUT_TSV}", flush=True)

    # summarize distributions
    groups = {
        "real": [r for r in enriched if not (r["any_dna_fp"] or r["ep_impure"])],
        "ep_impure": [r for r in enriched if r["ep_impure"]],
        "dna_fp": [r for r in enriched if r["any_dna_fp"]],
    }
    repeat_features = [
        "mean_repeat_frac", "mean_exon_repeat_frac", "max_exon_repeat_frac",
        "median_exon_repeat_frac", "repeat_rich_frac_0.5", "repeat_rich_frac_0.3",
        "repeat_weighted_mult", "high_repeat_hub_proxy",
    ]

    summary = {}
    for gname, grp in groups.items():
        summary[gname] = {}
        for f in repeat_features:
            vals = [r[f] for r in grp]
            summary[gname][f] = dict(
                n=len(vals),
                mean=round(statistics.mean(vals), 4) if vals else 0.0,
                median=round(statistics.median(vals), 4) if vals else 0.0,
                stdev=round(statistics.stdev(vals), 4) if len(vals) > 1 else 0.0,
            )

    # rule search on repeat-aware combinations
    def apply_rule(pred):
        flagged = [r for r in enriched if pred(r)]
        fp = sum(1 for r in flagged if r["any_dna_fp"])
        ep = sum(1 for r in flagged if r["ep_impure"])
        real = len(flagged) - fp - ep
        return dict(flagged=len(flagged), removed_dna_fp=fp,
                    removed_ep_impure=ep, collateral_real=real)

    rules = {
        "mean_exon_repeat_frac>=0.5": apply_rule(lambda r: r["mean_exon_repeat_frac"] >= 0.5),
        "max_exon_repeat_frac>=0.8": apply_rule(lambda r: r["max_exon_repeat_frac"] >= 0.8),
        "repeat_rich_frac_0.5>=0.5": apply_rule(lambda r: r["repeat_rich_frac_0.5"] >= 0.5),
        "repeat_weighted_mult>=10": apply_rule(lambda r: r["repeat_weighted_mult"] >= 10),
        "high_repeat_hub_proxy>=0.1": apply_rule(lambda r: r["high_repeat_hub_proxy"] >= 0.1),
        "mean_pair_mult>=20 AND mean_exon_repeat_frac>=0.3": apply_rule(
            lambda r: r["mean_pair_mult"] >= 20 and r["mean_exon_repeat_frac"] >= 0.3),
        "pair_hub_frac>=0.3 AND repeat_rich_frac_0.5>=0.3": apply_rule(
            lambda r: r["pair_hub_frac"] >= 0.3 and r["repeat_rich_frac_0.5"] >= 0.3),
        "n_members>=60 AND mean_exon_repeat_frac>=0.3": apply_rule(
            lambda r: r["n_members"] >= 60 and r["mean_exon_repeat_frac"] >= 0.3),
        "n_members>=50 AND repeat_weighted_mult>=5": apply_rule(
            lambda r: r["n_members"] >= 50 and r["repeat_weighted_mult"] >= 5),
    }

    summary["rules"] = rules

    with open(OUT_JSON, "w") as fh:
        json.dump(summary, fh, indent=2, sort_keys=True)
    print(f"[write] {OUT_JSON}", flush=True)

    # human-readable summary
    lines = []
    lines.append("Repeat/TE-aware feature distributions")
    lines.append(f"{'feature':<28} {'real_mean':>10} {'real_med':>10} {'ep_mean':>10} {'ep_med':>10} {'fp_mean':>10} {'fp_med':>10}")
    for f in repeat_features:
        vals = []
        for g in ["real", "ep_impure", "dna_fp"]:
            s = summary.get(g, {}).get(f, {})
            vals.extend([s.get("mean", 0), s.get("median", 0)])
        lines.append(f"{f:<28} {vals[0]:>10.4f} {vals[1]:>10.4f} {vals[2]:>10.4f} {vals[3]:>10.4f} {vals[4]:>10.4f} {vals[5]:>10.4f}")

    lines.append("")
    lines.append("Repeat-aware rule sweeps")
    lines.append(f"{'rule':<55} {'flag':>5} {'bad':>5} {'FP':>4} {'EP':>4} {'real$':>6}")
    for name, res in rules.items():
        bad = res["removed_dna_fp"] + res["removed_ep_impure"]
        lines.append(f"{name:<55} {res['flagged']:>5} {bad:>5} {res['removed_dna_fp']:>4} {res['removed_ep_impure']:>4} {res['collateral_real']:>6}")

    txt = "\n".join(lines)
    print(txt)
    with open(OUT_TXT, "w") as fh:
        fh.write(txt + "\n")
    print(f"[write] {OUT_TXT}", flush=True)


if __name__ == "__main__":
    main()
