#!/usr/bin/env python3

from __future__ import annotations

import argparse
import csv
from pathlib import Path


def parse_int(row: dict[str, str], key: str) -> int:
    value = row.get(key, "")
    return int(value) if value else 0


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--parity-tsv", type=Path, required=True)
    ap.add_argument("--out-tsv", type=Path, required=True)
    ap.add_argument("--out-md", type=Path, required=True)
    ap.add_argument("--top-n", type=int, default=25)
    args = ap.parse_args()

    rows: list[dict[str, object]] = []
    with args.parity_tsv.open() as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            stage_id = row["stage_id"]
            if stage_id not in {
                "seed_flow",
                "checktrf",
                "pairwise_overlap_filter",
                "isofrac",
                "collapse_single_exon_runoff",
                "polymerase_runoff_filter",
                "readthr_gate",
                "junction_support_filter",
                "final_family_selection",
            }:
                continue
            seed_fail = (
                parse_int(row, "seed_unwitnessed")
                + parse_int(row, "seed_zero_flux")
                + parse_int(row, "seed_low_cov")
                + parse_int(row, "seed_back_fail")
                + parse_int(row, "seed_fwd_fail")
            )
            unresolved_checktrf = max(
                0,
                parse_int(row, "checktrf_total") - parse_int(row, "seed_rescued"),
            )
            predcluster_removed = parse_int(row, "predcluster_removed")
            junction_support_removed = parse_int(row, "junction_support_removed")
            rows.append(
                {
                    "stage_id": stage_id,
                    "chrom": row["chrom"],
                    "start": parse_int(row, "start"),
                    "end": parse_int(row, "end"),
                    "strand": row["strand"],
                    "n_transcripts": parse_int(row, "n_transcripts"),
                    "seed_total": parse_int(row, "seed_total"),
                    "seed_stored": parse_int(row, "seed_stored"),
                    "seed_rescued": parse_int(row, "seed_rescued"),
                    "seed_unwitnessed": parse_int(row, "seed_unwitnessed"),
                    "seed_zero_flux": parse_int(row, "seed_zero_flux"),
                    "seed_low_cov": parse_int(row, "seed_low_cov"),
                    "seed_back_fail": parse_int(row, "seed_back_fail"),
                    "seed_fwd_fail": parse_int(row, "seed_fwd_fail"),
                    "checktrf_total": parse_int(row, "checktrf_total"),
                    "checktrf_pre_filter": parse_int(row, "checktrf_pre_filter"),
                    "predcluster_before": parse_int(row, "predcluster_before"),
                    "predcluster_removed": predcluster_removed,
                    "junction_support_removed": junction_support_removed,
                    "seed_fail_sum": seed_fail,
                    "unresolved_checktrf": unresolved_checktrf,
                    "priority_score": (
                        seed_fail * 10
                        + unresolved_checktrf * 5
                        + predcluster_removed * 3
                        + junction_support_removed * 4
                    ),
                    "note": row["note"],
                }
            )

    rows.sort(
        key=lambda row: (
            -int(row["priority_score"]),
            -int(row["seed_fail_sum"]),
            -int(row["unresolved_checktrf"]),
            str(row["chrom"]),
            int(row["start"]),
            str(row["stage_id"]),
        )
    )

    with args.out_tsv.open("w", newline="") as handle:
        fieldnames = [
            "stage_id",
            "chrom",
            "start",
            "end",
            "strand",
            "n_transcripts",
            "seed_total",
            "seed_stored",
            "seed_rescued",
            "seed_unwitnessed",
            "seed_zero_flux",
            "seed_low_cov",
            "seed_back_fail",
            "seed_fwd_fail",
            "checktrf_total",
            "checktrf_pre_filter",
            "predcluster_before",
            "predcluster_removed",
            "junction_support_removed",
            "seed_fail_sum",
            "unresolved_checktrf",
            "priority_score",
            "note",
        ]
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        for row in rows[: args.top_n]:
            writer.writerow(row)

    with args.out_md.open("w") as handle:
        handle.write("# Live Parity Failure Ranking\n\n")
        handle.write(f"Source parity TSV: `{args.parity_tsv}`\n\n")
        handle.write(f"Top rows: `{min(args.top_n, len(rows))}`\n\n")
        for idx, row in enumerate(rows[: args.top_n], start=1):
            handle.write(
                f"{idx}. `{row['stage_id']}` `{row['chrom']}:{row['start']}-{row['end']}({row['strand']})` "
                f"priority `{row['priority_score']}` seed_fail `{row['seed_fail_sum']}` "
                f"unresolved_checktrf `{row['unresolved_checktrf']}` "
                f"pred_removed `{row['predcluster_removed']}` junction_removed `{row['junction_support_removed']}` "
                f"stored `{row['seed_stored']}/{row['seed_total']}` rescued `{row['seed_rescued']}` "
                f"detail `{row['note']}`\n"
            )


if __name__ == "__main__":
    main()
