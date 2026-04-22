#!/usr/bin/env python3
"""Classify exact-chain boundary blockers into coarse failure families.

Inputs:
- exact-chain summary TSV from `trace_exact_chain_boundary_candidates.py`
- backward-oracle refs TSV

This is intentionally simple. It helps separate:
- path pair exists but emission shifts it later
- raw annotated pair exists but path never keeps it
- no annotated pair exists even at raw exact-chain level
- non-boundary exact-chain cases that should not drive broad rules
"""

from __future__ import annotations

import argparse
import csv
from collections import Counter
from pathlib import Path


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--exact-summary", required=True, type=Path)
    parser.add_argument("--oracle-refs", required=True, type=Path)
    parser.add_argument("--output-stem", required=True)
    return parser.parse_args()


def as_int(text: str) -> int:
    try:
        return int(text)
    except (TypeError, ValueError):
        return 0


def classify_row(summary: dict[str, str], oracle: dict[str, str] | None) -> tuple[str, str]:
    oracle_verdict = (oracle or {}).get("oracle_verdict", "")
    current_classification = (oracle or {}).get("current_classification", "")
    raw_has_ref = as_int(summary.get("raw_has_ref_boundary", "0"))
    path_has_ref = as_int(summary.get("path_has_ref_boundary", "0"))
    emitted_has_ref = as_int(summary.get("emitted_has_ref_boundary", "0"))
    path_count = as_int(summary.get("path_exact_candidate_count", "0"))
    emitted_count = as_int(summary.get("emitted_same_chain_count", "0"))
    best_path = summary.get("best_path_boundary", "")
    best_emitted = summary.get("best_emitted_boundary", "")

    if current_classification not in {"exact_chain", ""}:
        return "non_exact_chain_driver", "current best query is not a pure exact-chain case"

    if oracle_verdict == "boundary_only_current_miss":
        if path_has_ref and not emitted_has_ref:
            return "path_has_ref_pair_but_emission_loses_it", "same-chain ref boundary exists as path candidate but not in emitted output"
        if raw_has_ref and not path_has_ref:
            return "raw_has_ref_pair_but_path_never_keeps_it", "exact-chain raw reads contain ref boundary but the path layer does not"
        if path_count > 0 and emitted_count > 0 and best_path and best_emitted and best_path != best_emitted:
            return "path_pair_present_but_emission_shifted", "path family exists but final emission lands on a different same-chain boundary pair"
        if path_count > 0 and emitted_count > 0 and best_path == best_emitted:
            return "stable_nonref_same_chain_pair", "path and emission agree, but on a non-reference same-chain boundary pair"
        if path_count > 0 and emitted_count == 0:
            return "path_pair_present_but_not_emitted", "exact-chain path family exists but does not survive emission"
        if raw_has_ref or as_int(summary.get("raw_exact_boundary_count", "0")) > 0:
            return "raw_exact_chain_only", "boundary evidence exists in exact-chain reads but not in path/emission"
        return "no_ref_pair_even_at_raw_exact_chain", "no exact-chain raw read pair lands on the annotated boundary"

    if oracle_verdict in {"backward_keeps_exact_span", "canonical_ref_family_preferred_backward"}:
        return "backward_prefers_other_family_or_tail", oracle_verdict

    return "other", oracle_verdict or "unclassified"


def main() -> int:
    args = parse_args()
    oracle_rows: dict[str, dict[str, str]] = {}
    with args.oracle_refs.open() as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            oracle_rows[row["ref_id"]] = row

    out_rows: list[dict[str, str]] = []
    counts = Counter()
    with args.exact_summary.open() as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            ref_id = row["ref_id"]
            if not ref_id.startswith("STRG."):
                continue
            oracle = oracle_rows.get(ref_id)
            category, note = classify_row(row, oracle)
            out = {
                "ref_id": ref_id,
                "category": category,
                "note": note,
                "oracle_verdict": (oracle or {}).get("oracle_verdict", ""),
                "current_classification": (oracle or {}).get("current_classification", ""),
                "raw_has_ref_boundary": row.get("raw_has_ref_boundary", ""),
                "path_has_ref_boundary": row.get("path_has_ref_boundary", ""),
                "emitted_has_ref_boundary": row.get("emitted_has_ref_boundary", ""),
                "raw_exact_boundary_count": row.get("raw_exact_boundary_count", ""),
                "path_exact_candidate_count": row.get("path_exact_candidate_count", ""),
                "emitted_same_chain_count": row.get("emitted_same_chain_count", ""),
                "best_path_boundary": row.get("best_path_boundary", ""),
                "best_emitted_boundary": row.get("best_emitted_boundary", ""),
            }
            out_rows.append(out)
            counts[category] += 1

    out_rows.sort(key=lambda row: row["ref_id"])
    out_path = Path(f"{args.output_stem}.refs.tsv")
    with out_path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, delimiter="\t", fieldnames=list(out_rows[0].keys()))
        writer.writeheader()
        writer.writerows(out_rows)

    summary_path = Path(f"{args.output_stem}.summary.tsv")
    with summary_path.open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(["category", "count"])
        for category, count in sorted(counts.items()):
            writer.writerow([category, count])

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
