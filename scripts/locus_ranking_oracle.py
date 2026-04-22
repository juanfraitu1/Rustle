#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
from collections import Counter, defaultdict
from pathlib import Path


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Join backward verdicts, exact-chain candidates, and emitted denovo paths into a ranking oracle table."
    )
    parser.add_argument("--backward-oracle-refs", required=True, type=Path)
    parser.add_argument("--exact-candidates", required=True, type=Path)
    parser.add_argument("--path-emission", required=True, type=Path)
    parser.add_argument("--output-stem", required=True, type=Path)
    return parser.parse_args()


def load_tsv(path: Path) -> list[dict[str, str]]:
    with path.open() as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def span_key(start: str, end: str) -> str:
    return f"{start}-{end}"


def as_int(value: str | None) -> int:
    if not value:
        return 0
    return int(float(value))


def dedupe_candidates(rows: list[dict[str, str]]) -> list[dict[str, str]]:
    best: dict[tuple[str, str, str, str, str], dict[str, str]] = {}
    for row in rows:
        key = (
            row.get("ref_id", ""),
            row.get("source", ""),
            span_key(row.get("start", ""), row.get("end", "")),
            row.get("kind", ""),
            row.get("tx_id", ""),
        )
        prev = best.get(key)
        if prev is None or as_int(row.get("support")) > as_int(prev.get("support")):
            best[key] = row
    return list(best.values())


def main() -> int:
    args = parse_args()

    backward_rows = load_tsv(args.backward_oracle_refs)
    candidate_rows = dedupe_candidates(load_tsv(args.exact_candidates))
    path_rows = load_tsv(args.path_emission)

    emitted_by_tx: dict[str, dict[str, str]] = {}
    emitted_by_span: dict[str, list[dict[str, str]]] = defaultdict(list)
    path_by_span: dict[str, list[dict[str, str]]] = defaultdict(list)
    for row in path_rows:
        derived_span = span_key(row.get("derived_start", ""), row.get("derived_end", ""))
        path_by_span[derived_span].append(row)
        if row.get("emitted") == "1":
            emitted_by_span[derived_span].append(row)
            tx_id = row.get("emitted_tx_id", "")
            if tx_id:
                emitted_by_tx[tx_id] = row

    rows_out: list[dict[str, str]] = []
    summary = Counter()
    heuristic_prefer_refs: set[str] = set()
    ref_boundary_exact_full_refs: set[str] = set()
    ref_boundary_raw_refs: set[str] = set()

    candidates_by_ref: dict[str, list[dict[str, str]]] = defaultdict(list)
    for row in candidate_rows:
        candidates_by_ref[row["ref_id"]].append(row)

    for brow in backward_rows:
        ref_id = brow["ref_id"]
        ref_span = brow["ref_span"]
        ref_class = brow["current_classification"]
        oracle_verdict = brow["oracle_verdict"]
        current_best_tx = brow["current_best_query_id"]
        current_best_span = brow["current_best_query_span"]

        current_path = emitted_by_tx.get(current_best_tx, {})
        current_path_id = current_path.get("path_id", "")
        current_exact_count = current_path.get("exact_count", "")
        current_total_count = current_path.get("total_count", "")
        current_derived_chain = current_path.get("derived_junction_chain", "")

        ref_candidates = candidates_by_ref.get(ref_id, [])
        if not ref_candidates:
            rows_out.append(
                {
                    "ref_id": ref_id,
                    "ref_span": ref_span,
                    "oracle_verdict": oracle_verdict,
                    "current_classification": ref_class,
                    "current_best_query_id": current_best_tx,
                    "current_best_query_span": current_best_span,
                    "current_best_path_id": current_path_id,
                    "candidate_rank": "0",
                    "candidate_source": "",
                    "candidate_support": "0",
                    "candidate_span": "",
                    "candidate_is_ref_boundary": "0",
                    "candidate_tx_id": "",
                    "candidate_path_ids": "",
                    "candidate_exact_counts": "",
                    "candidate_total_counts": "",
                    "candidate_emitted_tx_ids": "",
                    "candidate_same_as_current_best": "0",
                    "candidate_has_backward_exact_span": brow["backward_exact_span_kept"],
                    "candidate_should_prefer": "0",
                    "candidate_note": "no_exact_chain_candidates",
                }
            )
            summary["refs_no_exact_chain_candidates"] += 1
            continue

        # Rank per-ref candidates by support desc, then source preference, then span.
        source_rank = {"exact_full": 0, "DENOVO:VG_CONSTRAINED_FLOW": 1, "total_support": 2, "bam_exact_chain": 3}
        ref_candidates = sorted(
            ref_candidates,
            key=lambda row: (
                -as_int(row.get("support")),
                source_rank.get(row.get("source", ""), 9),
                row.get("start", ""),
                row.get("end", ""),
            ),
        )

        best_ref_boundary_exact_full_support = 0
        for cand in ref_candidates:
            if cand.get("source") == "exact_full" and cand.get("is_ref_boundary") == "1":
                best_ref_boundary_exact_full_support = max(best_ref_boundary_exact_full_support, as_int(cand.get("support")))

        for idx, cand in enumerate(ref_candidates, start=1):
            candidate_span = span_key(cand.get("start", ""), cand.get("end", ""))
            linked_paths = path_by_span.get(candidate_span, [])
            linked_path_ids = ",".join(row.get("path_id", "") for row in linked_paths)
            linked_exact_counts = ",".join(row.get("exact_count", "") for row in linked_paths)
            linked_total_counts = ",".join(row.get("total_count", "") for row in linked_paths)
            linked_emitted_tx_ids = ",".join(
                row.get("emitted_tx_id", "") for row in linked_paths if row.get("emitted_tx_id", "")
            )

            same_as_current = "1" if (
                candidate_span == current_best_span
                or (cand.get("tx_id") and cand.get("tx_id") == current_best_tx)
            ) else "0"

            should_prefer = "0"
            note_parts: list[str] = []

            if cand.get("is_ref_boundary") == "1":
                note_parts.append("ref_boundary")
            if candidate_span == ref_span:
                note_parts.append("ref_span")
            if linked_paths:
                note_parts.append("linked_path")
            if cand.get("source") == "exact_full":
                note_parts.append("exact_full")
            if same_as_current == "1":
                note_parts.append("current_best")

            # Conservative research heuristic:
            # prefer only when the candidate is an exact_full ref-boundary family and the
            # current best is a non-ref boundary same-chain miss.
            if (
                oracle_verdict == "boundary_only_current_miss"
                and cand.get("source") == "exact_full"
                and cand.get("is_ref_boundary") == "1"
                and current_best_span != ref_span
                and candidate_span == ref_span
            ):
                should_prefer = "1"
                note_parts.append("heuristic_prefer")
            elif (
                oracle_verdict == "boundary_only_current_miss"
                and cand.get("source") == "bam_exact_chain"
                and cand.get("is_ref_boundary") == "1"
                and best_ref_boundary_exact_full_support == 0
            ):
                note_parts.append("raw_only_ref_boundary")

            rows_out.append(
                {
                    "ref_id": ref_id,
                    "ref_span": ref_span,
                    "oracle_verdict": oracle_verdict,
                    "current_classification": ref_class,
                    "current_best_query_id": current_best_tx,
                    "current_best_query_span": current_best_span,
                    "current_best_path_id": current_path_id,
                    "current_best_exact_count": current_exact_count,
                    "current_best_total_count": current_total_count,
                    "current_best_chain": current_derived_chain,
                    "candidate_rank": str(idx),
                    "candidate_source": cand.get("source", ""),
                    "candidate_support": cand.get("support", ""),
                    "candidate_span": candidate_span,
                    "candidate_is_ref_boundary": cand.get("is_ref_boundary", "0"),
                    "candidate_kind": cand.get("kind", ""),
                    "candidate_tx_id": cand.get("tx_id", ""),
                    "candidate_path_ids": linked_path_ids,
                    "candidate_exact_counts": linked_exact_counts,
                    "candidate_total_counts": linked_total_counts,
                    "candidate_emitted_tx_ids": linked_emitted_tx_ids,
                    "candidate_same_as_current_best": same_as_current,
                    "candidate_has_backward_exact_span": brow["backward_exact_span_kept"],
                    "candidate_should_prefer": should_prefer,
                    "candidate_note": ",".join(note_parts),
                }
            )

            if idx == 1:
                summary["refs_with_ranked_candidates"] += 1
            if cand.get("source") == "exact_full" and cand.get("is_ref_boundary") == "1":
                summary["ref_boundary_exact_full_rows"] += 1
                ref_boundary_exact_full_refs.add(ref_id)
            if cand.get("source") == "bam_exact_chain" and cand.get("is_ref_boundary") == "1":
                summary["ref_boundary_raw_rows"] += 1
                ref_boundary_raw_refs.add(ref_id)
            if should_prefer == "1":
                summary["heuristic_prefer_rows"] += 1
                heuristic_prefer_refs.add(ref_id)

    summary["heuristic_prefer_refs"] = len(heuristic_prefer_refs)
    summary["ref_boundary_exact_full_refs"] = len(ref_boundary_exact_full_refs)
    summary["ref_boundary_raw_refs"] = len(ref_boundary_raw_refs)

    out_refs = args.output_stem.with_suffix(".refs.tsv")
    with out_refs.open("w", newline="") as handle:
        fieldnames = [
            "ref_id",
            "ref_span",
            "oracle_verdict",
            "current_classification",
            "current_best_query_id",
            "current_best_query_span",
            "current_best_path_id",
            "current_best_exact_count",
            "current_best_total_count",
            "current_best_chain",
            "candidate_rank",
            "candidate_source",
            "candidate_support",
            "candidate_span",
            "candidate_is_ref_boundary",
            "candidate_kind",
            "candidate_tx_id",
            "candidate_path_ids",
            "candidate_exact_counts",
            "candidate_total_counts",
            "candidate_emitted_tx_ids",
            "candidate_same_as_current_best",
            "candidate_has_backward_exact_span",
            "candidate_should_prefer",
            "candidate_note",
        ]
        writer = csv.DictWriter(handle, delimiter="\t", fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows_out)

    out_summary = args.output_stem.with_suffix(".summary.tsv")
    with out_summary.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, delimiter="\t", fieldnames=["metric", "count"])
        writer.writeheader()
        for metric in sorted(summary):
            writer.writerow({"metric": metric, "count": summary[metric]})

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
