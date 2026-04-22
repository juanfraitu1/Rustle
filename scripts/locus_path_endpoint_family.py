#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
from collections import defaultdict
from pathlib import Path


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description="Summarize endpoint-family support from denovo graph path tables.")
    p.add_argument("--paths", required=True, type=Path)
    p.add_argument("--output", required=True, type=Path)
    return p.parse_args()


def main() -> int:
    args = parse_args()

    pair_rows: dict[tuple[int, int], dict[str, object]] = {}
    start_rows: dict[int, dict[str, object]] = {}
    end_rows: dict[int, dict[str, object]] = {}
    total_exact = 0

    with args.paths.open() as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            if row.get("source") != "exact_full":
                continue
            nodes = [int(x) for x in row.get("nodes", "").split(",") if x]
            if len(nodes) < 2:
                continue
            start = nodes[0]
            end = nodes[-1]
            count = int(float(row.get("count", "0") or 0))
            weight = float(row.get("weight", "0") or 0.0)
            total_exact += count

            pair_key = (start, end)
            pair = pair_rows.setdefault(
                pair_key,
                {
                    "start_node": start,
                    "end_node": end,
                    "count": 0,
                    "weight": 0.0,
                    "n_paths": 0,
                    "junction_chains": set(),
                },
            )
            pair["count"] += count
            pair["weight"] += weight
            pair["n_paths"] += 1
            pair["junction_chains"].add(row.get("junction_chain", ""))

            srow = start_rows.setdefault(
                start,
                {"start_node": start, "count": 0, "weight": 0.0, "distinct_end_nodes": set()},
            )
            srow["count"] += count
            srow["weight"] += weight
            srow["distinct_end_nodes"].add(end)

            erow = end_rows.setdefault(
                end,
                {"end_node": end, "count": 0, "weight": 0.0, "distinct_start_nodes": set()},
            )
            erow["count"] += count
            erow["weight"] += weight
            erow["distinct_start_nodes"].add(start)

    out_rows: list[dict[str, str]] = []
    for (start, end), row in sorted(
        pair_rows.items(),
        key=lambda item: (-int(item[1]["count"]), -float(item[1]["weight"]), item[0][0], item[0][1]),
    ):
        out_rows.append(
            {
                "start_node": str(start),
                "end_node": str(end),
                "count": str(row["count"]),
                "weight": f"{row['weight']:.3f}",
                "count_frac": f"{(row['count'] / total_exact):.4f}" if total_exact else "0.0000",
                "n_paths": str(row["n_paths"]),
                "junction_chain_count": str(len(row["junction_chains"])),
            }
        )

    with args.output.open("w", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            delimiter="\t",
            fieldnames=[
                "start_node",
                "end_node",
                "count",
                "weight",
                "count_frac",
                "n_paths",
                "junction_chain_count",
            ],
        )
        writer.writeheader()
        writer.writerows(out_rows)

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
