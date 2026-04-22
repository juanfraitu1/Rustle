#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import re
from collections import defaultdict
from pathlib import Path


EXON_PAIR_RE = re.compile(r"(\d+)-(\d+)")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Compare current VG endpoint-family partitioning against backward prediction families."
    )
    parser.add_argument("--current-path-emission", required=True, type=Path)
    parser.add_argument("--current-read-path-map", required=True, type=Path)
    parser.add_argument("--backward-predictions", required=True, type=Path)
    parser.add_argument("--target-start", required=True, type=int)
    parser.add_argument("--target-end", required=True, type=int)
    parser.add_argument("--output-stem", required=True, type=Path)
    return parser.parse_args()


def load_tsv(path: Path) -> list[dict[str, str]]:
    with path.open() as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def first_donor_from_chain(chain: str) -> int | None:
    if not chain:
        return None
    token = chain.split(",")[0]
    donor, _acceptor = token.split("-")
    return int(donor)


def parse_pred_exons(raw: str, span_start: int, span_end: int) -> list[tuple[int, int]]:
    pairs = [(int(a), int(b)) for a, b in EXON_PAIR_RE.findall(raw)]
    if len(pairs) >= 2 and pairs[0] == (span_start, span_end):
        pairs = pairs[1:]
    return pairs


def main() -> int:
    args = parse_args()
    current_paths = load_tsv(args.current_path_emission)
    current_reads = load_tsv(args.current_read_path_map)
    backward_preds = load_tsv(args.backward_predictions)

    selected_paths: dict[str, dict[str, str]] = {}
    current_summary_rows: list[dict[str, str]] = []
    current_group_summary: dict[tuple[int, int], dict[str, object]] = {}

    for row in current_paths:
        if int(row["derived_start"]) != args.target_start:
            continue
        first_donor = first_donor_from_chain(row.get("derived_junction_chain", ""))
        if first_donor is None:
            continue
        derived_end = int(row["derived_end"])
        key = (first_donor, derived_end)
        group = current_group_summary.setdefault(
            key,
            {
                "first_donor": first_donor,
                "derived_end": derived_end,
                "n_paths": 0,
                "sum_exact_count": 0,
                "sum_total_count": 0,
                "path_ids": [],
                "emitted_tx_ids": [],
            },
        )
        group["n_paths"] += 1
        group["sum_exact_count"] += int(row["exact_count"])
        group["sum_total_count"] += int(row["total_count"])
        group["path_ids"].append(row["path_id"])
        if row.get("emitted_tx_id"):
            group["emitted_tx_ids"].append(row["emitted_tx_id"])
        selected_paths[row["path_id"]] = row
        current_summary_rows.append(
            {
                "source": "current_path",
                "path_id": row["path_id"],
                "first_donor": str(first_donor),
                "derived_end": row["derived_end"],
                "exact_count": row["exact_count"],
                "total_count": row["total_count"],
                "emitted_tx_id": row.get("emitted_tx_id", ""),
                "emitted_end": row.get("emitted_end", ""),
            }
        )

    current_group_rows: list[dict[str, str]] = []
    for (_first_donor, _end), row in sorted(
        current_group_summary.items(),
        key=lambda item: (-int(item[1]["sum_exact_count"]), item[0][0], item[0][1]),
    ):
        current_group_rows.append(
            {
                "first_donor": str(row["first_donor"]),
                "derived_end": str(row["derived_end"]),
                "n_paths": str(row["n_paths"]),
                "sum_exact_count": str(row["sum_exact_count"]),
                "sum_total_count": str(row["sum_total_count"]),
                "path_ids": ",".join(row["path_ids"]),
                "emitted_tx_ids": ",".join(row["emitted_tx_ids"]),
            }
        )

    current_read_rows: list[dict[str, str]] = []
    current_read_group: dict[tuple[int, int, int, int, int], int] = defaultdict(int)
    for row in current_reads:
        path_id = row["path_id"]
        path = selected_paths.get(path_id)
        if not path:
            continue
        first_donor = first_donor_from_chain(path.get("derived_junction_chain", ""))
        if first_donor is None:
            continue
        key = (
            first_donor,
            int(path["derived_end"]),
            int(row["observed_end"]),
            int(row["normalized_end"]),
            int(row["tx_end"]),
        )
        current_read_group[key] += int(row["read_support"])
        current_read_rows.append(
            {
                "path_id": path_id,
                "first_donor": str(first_donor),
                "derived_end": path["derived_end"],
                "read_support": row["read_support"],
                "observed_end": row["observed_end"],
                "normalized_end": row["normalized_end"],
                "tx_end": row["tx_end"],
            }
        )

    current_read_group_rows: list[dict[str, str]] = []
    for key, support in sorted(current_read_group.items(), key=lambda item: (-item[1], *item[0])):
        first_donor, derived_end, observed_end, normalized_end, tx_end = key
        current_read_group_rows.append(
            {
                "first_donor": str(first_donor),
                "derived_end": str(derived_end),
                "observed_end": str(observed_end),
                "normalized_end": str(normalized_end),
                "tx_end": str(tx_end),
                "read_support_sum": str(support),
            }
        )

    backward_pred_rows: list[dict[str, str]] = []
    backward_group: dict[tuple[int, int], dict[str, object]] = {}
    final_fates: dict[str, str] = {}
    for row in backward_preds:
        if row.get("event") == "FINAL_KEEP":
            final_fates[row.get("pred_idx", "")] = "FINAL_KEEP"
        elif row.get("event") == "FILTER":
            final_fates[row.get("pred_idx", "")] = row.get("reason", "FILTER")

    for row in backward_preds:
        if row.get("event") != "PRED":
            continue
        if int(row["span_start"]) != args.target_start or int(row["span_end"]) != args.target_end:
            continue
        exons = parse_pred_exons(row.get("raw", ""), int(row["span_start"]), int(row["span_end"]))
        if len(exons) < 2:
            continue
        first_donor = exons[0][1]
        exons_count = len(exons)
        pred_idx = row["pred_idx"]
        fate = final_fates.get(pred_idx, "")
        key = (first_donor, int(row["span_end"]))
        group = backward_group.setdefault(
            key,
            {
                "first_donor": first_donor,
                "span_end": int(row["span_end"]),
                "n_preds": 0,
                "cov_sum": 0.0,
                "pred_indices": [],
                "fates": [],
            },
        )
        group["n_preds"] += 1
        group["cov_sum"] += float(row["cov"] or 0.0)
        group["pred_indices"].append(pred_idx)
        if fate:
            group["fates"].append(f"{pred_idx}:{fate}")
        backward_pred_rows.append(
            {
                "pred_idx": pred_idx,
                "first_donor": str(first_donor),
                "span_end": row["span_end"],
                "cov": row["cov"],
                "exons": str(exons_count),
                "final_fate": fate,
            }
        )

    backward_group_rows: list[dict[str, str]] = []
    for (_first_donor, _end), row in sorted(
        backward_group.items(),
        key=lambda item: (-float(item[1]["cov_sum"]), item[0][0], item[0][1]),
    ):
        backward_group_rows.append(
            {
                "first_donor": str(row["first_donor"]),
                "span_end": str(row["span_end"]),
                "n_preds": str(row["n_preds"]),
                "cov_sum": f"{row['cov_sum']:.4f}",
                "pred_indices": ",".join(row["pred_indices"]),
                "fates": ",".join(row["fates"]),
            }
        )

    def write_rows(path: Path, fieldnames: list[str], rows: list[dict[str, str]]) -> None:
        with path.open("w", newline="") as handle:
            writer = csv.DictWriter(handle, delimiter="\t", fieldnames=fieldnames)
            writer.writeheader()
            writer.writerows(rows)

    write_rows(
        args.output_stem.with_name(args.output_stem.name + ".current_paths.tsv"),
        ["source", "path_id", "first_donor", "derived_end", "exact_count", "total_count", "emitted_tx_id", "emitted_end"],
        current_summary_rows,
    )
    write_rows(
        args.output_stem.with_name(args.output_stem.name + ".current_groups.tsv"),
        ["first_donor", "derived_end", "n_paths", "sum_exact_count", "sum_total_count", "path_ids", "emitted_tx_ids"],
        current_group_rows,
    )
    write_rows(
        args.output_stem.with_name(args.output_stem.name + ".current_read_groups.tsv"),
        ["first_donor", "derived_end", "observed_end", "normalized_end", "tx_end", "read_support_sum"],
        current_read_group_rows,
    )
    write_rows(
        args.output_stem.with_name(args.output_stem.name + ".backward_preds.tsv"),
        ["pred_idx", "first_donor", "span_end", "cov", "exons", "final_fate"],
        backward_pred_rows,
    )
    write_rows(
        args.output_stem.with_name(args.output_stem.name + ".backward_groups.tsv"),
        ["first_donor", "span_end", "n_preds", "cov_sum", "pred_indices", "fates"],
        backward_group_rows,
    )

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
