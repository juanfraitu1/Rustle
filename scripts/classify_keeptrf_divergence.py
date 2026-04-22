#!/usr/bin/env python3
"""Classify normalized keeptrf vs native divergence relations."""

from __future__ import annotations

import argparse
import csv
from collections import Counter
from pathlib import Path


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("details_tsv", type=Path)
    parser.add_argument("--output-prefix", required=True)
    return parser.parse_args()


def parse_chain(chain: str) -> list[tuple[int, int]]:
    out = []
    for part in chain.split(","):
        part = part.strip()
        if not part or "-" not in part:
            continue
        a, b = part.split("-", 1)
        try:
            out.append((int(a), int(b)))
        except ValueError:
            continue
    return out


def classify_relation(lhs: str, rhs: str) -> str:
    a = parse_chain(lhs)
    b = parse_chain(rhs)
    if not a or not b:
        return "empty_or_unparsed"
    if a == b:
        return "exact"
    if len(a) == len(b):
        diff = [i for i, (x, y) in enumerate(zip(a, b)) if x != y]
        if len(diff) == 1:
            idx = diff[0]
            if idx == len(a) - 1:
                return "terminal_single_junction_substitution"
            return "internal_single_junction_substitution"
        if len(diff) == 2 and diff == [len(a) - 2, len(a) - 1]:
            return "terminal_two_junction_substitution"
        return "same_length_other"
    shorter, longer = (a, b) if len(a) < len(b) else (b, a)
    if longer[: len(shorter)] == shorter:
        return "terminal_extension_or_truncation"
    if longer[-len(shorter) :] == shorter:
        return "leading_extension_or_truncation"
    return "different_length_other"


def main() -> int:
    args = parse_args()
    rows = list(csv.DictReader(args.details_tsv.open(), delimiter="\t"))

    summary_path = Path(f"{args.output_prefix}.summary.tsv")
    details_path = Path(f"{args.output_prefix}.details.tsv")

    counts = Counter()
    with details_path.open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(
            [
                "rep_idx",
                "first_divergence",
                "normalized_chain",
                "comparison_chain",
                "relation_class",
            ]
        )
        for row in rows:
            divergence = row["first_divergence"]
            comparison = ""
            if divergence == "mode_member_only":
                comparison = row.get("selected_competitor_chains", "").split(" || ")[0]
            elif divergence == "boundary_only_pre_family":
                comparison = row.get("boundary_emit_replacements", "").split(" || ")[0]
            if not comparison:
                continue
            relation = classify_relation(row["normalized_chain"], comparison)
            counts[(divergence, relation)] += 1
            writer.writerow(
                [
                    row["rep_idx"],
                    divergence,
                    row["normalized_chain"],
                    comparison,
                    relation,
                ]
            )

    with summary_path.open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(["first_divergence", "relation_class", "count"])
        for (divergence, relation), count in sorted(counts.items()):
            writer.writerow([divergence, relation, count])
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
