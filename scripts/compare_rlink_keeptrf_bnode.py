#!/usr/bin/env python3
"""Compare rlink keeptrf/usepath representatives against native BNODE trace stages."""

from __future__ import annotations

import argparse
import csv
from collections import defaultdict
from dataclasses import dataclass
from pathlib import Path
import re


@dataclass
class CurrentRow:
    event: str
    family_key: str
    mode_idx: str
    path_idx: str
    support: int
    mode_total_reads: int
    gate: int
    decision: str
    reason: str
    emit_chain: str
    boundary_chain: str
    selected_path_idx: str


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--rlink-keeptrf", required=True, type=Path)
    parser.add_argument("--current-trace", required=True, type=Path)
    parser.add_argument("--output-prefix", required=True)
    return parser.parse_args()


def read_keeptrf(path: Path) -> list[dict[str, str]]:
    with path.open() as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def read_current_trace(path: Path) -> list[CurrentRow]:
    out: list[CurrentRow] = []
    with path.open() as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            out.append(
                CurrentRow(
                    event=row.get("event", ""),
                    family_key=row.get("family_key", ""),
                    mode_idx=row.get("mode_idx", ""),
                    path_idx=row.get("path_idx", ""),
                    support=int(row.get("support") or 0),
                    mode_total_reads=int(row.get("mode_total_reads") or 0),
                    gate=int(row.get("gate") or 0),
                    decision=row.get("decision", ""),
                    reason=row.get("reason", ""),
                    emit_chain=row.get("emit_chain", ""),
                    boundary_chain=row.get("boundary_chain", ""),
                    selected_path_idx=row.get("selected_path_idx", ""),
                )
            )
    return out


CHAIN_PART_RE = re.compile(r"^(\d+)-(\d+)$")


def normalize_rlink_chain(chain: str) -> str:
    parts = []
    for part in chain.split(","):
        part = part.strip()
        if not part:
            continue
        match = CHAIN_PART_RE.match(part)
        if not match:
            parts.append(part)
            continue
        donor, acceptor = int(match.group(1)), int(match.group(2))
        parts.append(f"{donor}-{acceptor + 1}")
    return ",".join(parts)


def main() -> int:
    args = parse_args()
    keeptrf = read_keeptrf(args.rlink_keeptrf)
    current = read_current_trace(args.current_trace)

    family_input_by_chain: dict[str, list[CurrentRow]] = defaultdict(list)
    family_input_by_boundary_chain: dict[str, list[CurrentRow]] = defaultdict(list)
    mode_rows_by_chain: dict[str, list[CurrentRow]] = defaultdict(list)
    selected_by_family_mode: dict[tuple[str, str], CurrentRow] = {}

    for row in current:
        if row.event == "family_input":
            family_input_by_chain[row.emit_chain].append(row)
            family_input_by_boundary_chain[row.boundary_chain].append(row)
        elif row.event == "mode_eval":
            mode_rows_by_chain[row.emit_chain].append(row)
            if row.decision == "selected":
                selected_by_family_mode[(row.family_key, row.mode_idx)] = row

    summary_path = Path(f"{args.output_prefix}.summary.tsv")
    details_path = Path(f"{args.output_prefix}.details.tsv")

    counts = defaultdict(int)
    with details_path.open("w", newline="") as detail_handle:
        writer = csv.writer(detail_handle, delimiter="\t")
        writer.writerow(
            [
                "rep_idx",
                "usepath",
                "group_cov",
                "seed_complete",
                "junction_chain",
                "normalized_chain",
                "family_input_count",
                "family_input_boundary_count",
                "mode_eval_count",
                "selected_count",
                "first_divergence",
                "current_families",
                "candidate_mode_keys",
                "boundary_emit_replacements",
                "selected_competitor_chains",
                "selected_competitor_supports",
            ]
        )

        for row in keeptrf:
            chain = row.get("junction_chain", "")
            normalized_chain = normalize_rlink_chain(chain)
            family_rows = family_input_by_chain.get(normalized_chain, [])
            boundary_rows = family_input_by_boundary_chain.get(normalized_chain, [])
            mode_rows = mode_rows_by_chain.get(normalized_chain, [])
            selected_rows = [r for r in mode_rows if r.decision == "selected"]

            if not family_rows:
                if boundary_rows:
                    divergence = "boundary_only_pre_family"
                else:
                    divergence = "missing_pre_family"
            elif not mode_rows:
                divergence = "missing_mode_eval"
            elif selected_rows:
                divergence = "selected_exact"
            else:
                divergence = "mode_member_only"

            competitor_chains: list[str] = []
            competitor_supports: list[str] = []
            mode_keys: list[str] = []
            boundary_emit_replacements = sorted({r.emit_chain for r in boundary_rows if r.emit_chain})
            for cand in mode_rows:
                key = (cand.family_key, cand.mode_idx)
                mode_keys.append(f"{cand.family_key}:{cand.mode_idx}")
                selected = selected_by_family_mode.get(key)
                if selected and selected.emit_chain != chain:
                    competitor_chains.append(selected.emit_chain)
                    competitor_supports.append(str(selected.support))

            counts[divergence] += 1
            writer.writerow(
                [
                    row.get("rep_idx", ""),
                    row.get("usepath", ""),
                    row.get("group_cov", ""),
                    row.get("seed_complete", ""),
                    chain,
                    normalized_chain,
                    len(family_rows),
                    len(boundary_rows),
                    len(mode_rows),
                    len(selected_rows),
                    divergence,
                    ",".join(sorted({r.family_key for r in family_rows})),
                    " || ".join(sorted(set(mode_keys))),
                    " || ".join(boundary_emit_replacements),
                    " || ".join(sorted(set(competitor_chains))),
                    " || ".join(sorted(set(competitor_supports))),
                ]
            )

    with summary_path.open("w", newline="") as summary_handle:
        writer = csv.writer(summary_handle, delimiter="\t")
        writer.writerow(["metric", "value"])
        writer.writerow(["rlink_keeptrf_rows", len(keeptrf)])
        writer.writerow(["current_trace_rows", len(current)])
        for key in sorted(counts):
            writer.writerow([key, counts[key]])

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
