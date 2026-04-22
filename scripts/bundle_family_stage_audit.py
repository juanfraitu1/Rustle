#!/usr/bin/env python3
"""Summarize current bundle trace stages by family key.

Input: normalized current trace TSV from trace_decision_extract.py
Output:
- <prefix>.summary.tsv
- <prefix>.families.tsv
- <prefix>.corrected_groups.tsv
"""

from __future__ import annotations

import argparse
import csv
from collections import Counter, defaultdict
from pathlib import Path


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("input", type=Path)
    parser.add_argument("--output-prefix", required=True)
    parser.add_argument("--top-n", type=int, default=40)
    return parser.parse_args()


def parse_reason(reason: str) -> dict[str, str]:
    out: dict[str, str] = {}
    for part in reason.split(";"):
        if "=" not in part:
            continue
        key, value = part.split("=", 1)
        out[key] = value
    return out


def as_int(value: str) -> int:
    try:
        return int(value)
    except (TypeError, ValueError):
        return 0


def update_span(stats: dict, start: int, end: int) -> None:
    if start <= 0 or end <= 0:
        return
    if stats["span_start"] == 0 or start < stats["span_start"]:
        stats["span_start"] = start
    if end > stats["span_end"]:
        stats["span_end"] = end


def main() -> int:
    args = parse_args()
    in_path = args.input
    out_prefix = Path(args.output_prefix)
    families: dict[str, dict] = defaultdict(
        lambda: {
            "bundle_chrom": "",
            "bundle_strand": "",
            "span_start": 0,
            "span_end": 0,
            "corrected_key": "",
            "family_input_rows": 0,
            "family_input_support_sum": 0,
            "mode_selected_rows": 0,
            "mode_selected_support_sum": 0,
            "mode_member_rows": 0,
            "protected_raw_split_rows": 0,
            "protected_raw_split_support_sum": 0,
            "corrected_key_rows": 0,
            "corrected_key_support_sum": 0,
            "corrected_same_as_raw_rows": 0,
            "corrected_same_as_raw_support_sum": 0,
            "other_reason_rows": 0,
        }
    )
    totals = Counter()

    with in_path.open() as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            family_key = row.get("family_key", "")
            if not family_key:
                continue
            stats = families[family_key]
            stats["bundle_chrom"] = row.get("bundle_chrom", "") or stats["bundle_chrom"]
            stats["bundle_strand"] = row.get("bundle_strand", "") or stats["bundle_strand"]
            update_span(stats, as_int(row.get("span_start", "")), as_int(row.get("span_end", "")))
            support = as_int(row.get("support", ""))
            stage = row.get("event", "")
            decision = row.get("decision", "")
            reason_bits = parse_reason(row.get("reason", ""))
            family_select = reason_bits.get("family_select", "")

            totals[f"{stage}:{decision or '-'}"] += 1

            if stage == "family_input":
                stats["family_input_rows"] += 1
                stats["family_input_support_sum"] += support
                if family_select == "protected_raw_split":
                    stats["corrected_key"] = reason_bits.get("corrected_key", "") or stats["corrected_key"]
                    stats["protected_raw_split_rows"] += 1
                    stats["protected_raw_split_support_sum"] += support
                elif family_select == "corrected_key":
                    stats["corrected_key_rows"] += 1
                    stats["corrected_key_support_sum"] += support
                elif family_select == "corrected_same_as_raw":
                    stats["corrected_same_as_raw_rows"] += 1
                    stats["corrected_same_as_raw_support_sum"] += support
                else:
                    stats["other_reason_rows"] += 1
            elif stage == "mode_eval":
                if decision == "selected":
                    stats["mode_selected_rows"] += 1
                    stats["mode_selected_support_sum"] += support
                elif decision == "mode_member":
                    stats["mode_member_rows"] += 1

    family_rows = []
    for family_key, stats in families.items():
        row = dict(stats)
        row["family_key"] = family_key
        row["protected_split_fraction"] = (
            f"{row['protected_raw_split_rows'] / row['family_input_rows']:.4f}"
            if row["family_input_rows"]
            else "0.0000"
        )
        row["selected_per_input"] = (
            f"{row['mode_selected_rows'] / row['family_input_rows']:.4f}"
            if row["family_input_rows"]
            else "0.0000"
        )
        family_rows.append(row)

    family_rows.sort(
        key=lambda row: (
            row["protected_raw_split_rows"],
            row["mode_selected_rows"],
            row["family_input_support_sum"],
            row["span_end"] - row["span_start"],
        ),
        reverse=True,
    )

    summary_path = out_prefix.with_suffix(".summary.tsv")
    families_path = out_prefix.with_suffix(".families.tsv")
    corrected_groups_path = out_prefix.with_suffix(".corrected_groups.tsv")

    with summary_path.open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(["metric", "value"])
        writer.writerow(["families_total", len(family_rows)])
        writer.writerow(["family_input_rows_total", sum(row["family_input_rows"] for row in family_rows)])
        writer.writerow(["mode_selected_rows_total", sum(row["mode_selected_rows"] for row in family_rows)])
        writer.writerow(
            ["protected_raw_split_rows_total", sum(row["protected_raw_split_rows"] for row in family_rows)]
        )
        writer.writerow(
            [
                "families_with_any_protected_raw_split",
                sum(1 for row in family_rows if row["protected_raw_split_rows"] > 0),
            ]
        )
        writer.writerow(
            [
                "families_with_selected_and_protected_raw_split",
                sum(
                    1
                    for row in family_rows
                    if row["protected_raw_split_rows"] > 0 and row["mode_selected_rows"] > 0
                ),
            ]
        )
        for key in sorted(totals):
            writer.writerow([f"stage_count:{key}", totals[key]])

    fields = [
        "bundle_chrom",
        "bundle_strand",
        "span_start",
        "span_end",
        "family_key",
        "family_input_rows",
        "family_input_support_sum",
        "mode_selected_rows",
        "mode_selected_support_sum",
        "mode_member_rows",
        "protected_raw_split_rows",
        "protected_raw_split_support_sum",
        "corrected_key_rows",
        "corrected_key_support_sum",
        "corrected_same_as_raw_rows",
        "corrected_same_as_raw_support_sum",
        "other_reason_rows",
        "protected_split_fraction",
        "selected_per_input",
    ]
    with families_path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        for row in family_rows[: args.top_n]:
            writer.writerow({field: row[field] for field in fields})

    corrected_groups: dict[str, dict] = defaultdict(
        lambda: {
            "bundle_chrom": "",
            "bundle_strand": "",
            "span_start": 0,
            "span_end": 0,
            "families": 0,
            "family_input_rows": 0,
            "mode_selected_rows": 0,
            "protected_raw_split_rows": 0,
            "protected_raw_split_support_sum": 0,
        }
    )
    for row in family_rows:
        corrected_key = row["corrected_key"]
        if not corrected_key:
            continue
        group = corrected_groups[corrected_key]
        group["bundle_chrom"] = row["bundle_chrom"] or group["bundle_chrom"]
        group["bundle_strand"] = row["bundle_strand"] or group["bundle_strand"]
        if group["span_start"] == 0 or row["span_start"] < group["span_start"]:
            group["span_start"] = row["span_start"]
        if row["span_end"] > group["span_end"]:
            group["span_end"] = row["span_end"]
        group["families"] += 1
        group["family_input_rows"] += row["family_input_rows"]
        group["mode_selected_rows"] += row["mode_selected_rows"]
        group["protected_raw_split_rows"] += row["protected_raw_split_rows"]
        group["protected_raw_split_support_sum"] += row["protected_raw_split_support_sum"]

    corrected_rows = []
    for corrected_key, stats in corrected_groups.items():
        row = dict(stats)
        row["corrected_key"] = corrected_key
        row["selected_per_family"] = (
            f"{row['mode_selected_rows'] / row['families']:.4f}" if row["families"] else "0.0000"
        )
        corrected_rows.append(row)
    corrected_rows.sort(
        key=lambda row: (
            row["protected_raw_split_rows"],
            row["families"],
            row["mode_selected_rows"],
            row["span_end"] - row["span_start"],
        ),
        reverse=True,
    )

    corrected_fields = [
        "bundle_chrom",
        "bundle_strand",
        "span_start",
        "span_end",
        "corrected_key",
        "families",
        "family_input_rows",
        "mode_selected_rows",
        "protected_raw_split_rows",
        "protected_raw_split_support_sum",
        "selected_per_family",
    ]
    with corrected_groups_path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=corrected_fields, delimiter="\t")
        writer.writeheader()
        for row in corrected_rows[: args.top_n]:
            writer.writerow({field: row[field] for field in corrected_fields})

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
