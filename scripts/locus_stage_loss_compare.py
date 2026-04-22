#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import re
from collections import Counter, defaultdict
from pathlib import Path


PAIR_RE = re.compile(r"(\d+)-(\d+)")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Build a stage-by-stage loss/replacement table for one locus family."
    )
    parser.add_argument("--backward-read-pattern-paths", type=Path)
    parser.add_argument("--backward-read-graph-paths", required=True, type=Path)
    parser.add_argument("--backward-predictions", required=True, type=Path)
    parser.add_argument("--current-read-path-map", required=True, type=Path)
    parser.add_argument("--current-path-emission", required=True, type=Path)
    parser.add_argument("--current-derived-start", required=True, type=int)
    parser.add_argument("--backward-pred-span-start", type=int)
    parser.add_argument("--backward-node0-start", required=True, type=int)
    parser.add_argument("--target-end", required=True, type=int)
    parser.add_argument("--output-stem", required=True, type=Path)
    return parser.parse_args()


def load_tsv(path: Path) -> list[dict[str, str]]:
    with path.open() as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def parse_pred_exons(raw: str, span_start: int, span_end: int) -> list[tuple[int, int]]:
    pairs = [(int(a), int(b)) for a, b in PAIR_RE.findall(raw)]
    if len(pairs) >= 2 and pairs[0] == (span_start, span_end):
        pairs = pairs[1:]
    return pairs


def first_donor_from_chain(chain: str) -> int | None:
    if not chain:
        return None
    token = chain.split(",", 1)[0]
    donor, _acceptor = token.split("-", 1)
    return int(donor)


def parse_node_spans(text: str) -> list[tuple[int, int]]:
    spans: list[tuple[int, int]] = []
    for part in filter(None, text.split(";")):
        if "-" not in part:
            continue
        start, end = part.split("-", 1)
        spans.append((int(start), int(end)))
    return spans


def pred_label(row: dict[str, str]) -> str:
    pred_idx = row.get("pred_idx", "")
    if pred_idx:
        return f"pred[{pred_idx}]"
    return row.get("pred_id", "") or row.get("raw", "").split(":", 1)[0]


def main() -> int:
    args = parse_args()
    backward_read_paths = load_tsv(args.backward_read_graph_paths)
    backward_read_pattern_paths = load_tsv(args.backward_read_pattern_paths) if args.backward_read_pattern_paths else []
    backward_predictions = load_tsv(args.backward_predictions)
    current_read_path_map = load_tsv(args.current_read_path_map)
    current_path_emission = load_tsv(args.current_path_emission)
    backward_pred_span_start = args.backward_pred_span_start or args.current_derived_start

    # Stage 0: backward PATH_read_pattern terminal states
    backward_pattern_terminal_counts = Counter()
    for row in backward_read_pattern_paths:
        spans = parse_node_spans(row.get("node_spans", ""))
        if len(spans) < 2:
            continue
        first_donor = spans[0][1]
        last_end = spans[-1][1]
        backward_pattern_terminal_counts[(first_donor, last_end)] += 1

    # Stage 1: backward read-path terminal states
    backward_terminal_counts = Counter()
    for row in backward_read_paths:
        if row.get("node0_start") != str(args.backward_node0_start):
            continue
        end = row.get("nodeLast_end", "")
        if end:
            backward_terminal_counts[int(end)] += 1

    # Stage 2: backward exact-end prediction families
    backward_pred_families: dict[int, dict[str, object]] = {}
    for row in backward_predictions:
        if row.get("event") != "PRED":
            continue
        if row.get("span_start") != str(backward_pred_span_start):
            continue
        if row.get("span_end") != str(args.target_end):
            continue
        exons = parse_pred_exons(row.get("raw", ""), backward_pred_span_start, args.target_end)
        if len(exons) < 2:
            continue
        first_donor = exons[0][1]
        entry = backward_pred_families.setdefault(
            first_donor,
            {"n_preds": 0, "cov_sum": 0.0, "pred_idxs": []},
        )
        entry["n_preds"] += 1
        entry["cov_sum"] += float(row.get("cov") or 0.0)
        entry["pred_idxs"].append(pred_label(row))

    # Stage 3: current raw-read family states before emission
    selected_paths: dict[str, dict[str, str]] = {}
    current_emission_families: dict[tuple[int, int], dict[str, object]] = {}
    for row in current_path_emission:
        if row.get("derived_start") != str(args.current_derived_start):
            continue
        donor = first_donor_from_chain(row.get("derived_junction_chain", ""))
        if donor is None:
            continue
        key = (donor, int(row["derived_end"]))
        entry = current_emission_families.setdefault(
            key,
            {
                "exact_count": 0,
                "total_count": 0,
                "path_ids": [],
                "emitted_tx_ids": [],
            },
        )
        entry["exact_count"] += int(row.get("exact_count") or 0)
        entry["total_count"] += int(row.get("total_count") or 0)
        entry["path_ids"].append(row["path_id"])
        if row.get("emitted_tx_id"):
            entry["emitted_tx_ids"].append(row["emitted_tx_id"])
        selected_paths[row["path_id"]] = row

    current_read_families: dict[tuple[int, int, int, int], float] = defaultdict(float)
    for row in current_read_path_map:
        path = selected_paths.get(row.get("path_id", ""))
        if not path:
            continue
        donor = first_donor_from_chain(path.get("derived_junction_chain", ""))
        if donor is None:
            continue
        key = (
            donor,
            int(row["observed_end"]),
            int(row["normalized_end"]),
            int(row["tx_end"]),
        )
        current_read_families[key] += float(row.get("read_support") or 0.0)

    donors = sorted(
        set(backward_pred_families)
        | {donor for donor, _ in current_emission_families}
        | {donor for donor, *_ in current_read_families},
    )

    rows: list[dict[str, str]] = []
    for donor in donors:
        backward_pred = backward_pred_families.get(donor)
        read_variants = sorted(
            [
                (obs_end, norm_end, tx_end, support)
                for (d, obs_end, norm_end, tx_end), support in current_read_families.items()
                if d == donor
            ],
            key=lambda item: (-item[3], item[0], item[1], item[2]),
        )
        emit_variants = sorted(
            [
                (derived_end, info["exact_count"], info["total_count"], ",".join(info["path_ids"]))
                for (d, derived_end), info in current_emission_families.items()
                if d == donor
            ],
            key=lambda item: (-item[1], -item[2], item[0]),
        )
        rows.append(
            {
                "first_donor": str(donor),
                "backward_pattern_exact_end_rows": str(backward_pattern_terminal_counts.get((donor, args.target_end), 0)),
                "backward_locus_exact_end_readpath_rows": str(backward_terminal_counts.get(args.target_end, 0)),
                "backward_exact_pred_cov_sum": f"{(backward_pred or {}).get('cov_sum', 0.0):.4f}",
                "backward_exact_pred_n": str((backward_pred or {}).get("n_preds", 0)),
                "backward_exact_pred_ids": ",".join((backward_pred or {}).get("pred_idxs", [])),
                "current_read_variants": "; ".join(
                    f"obs={obs} norm={norm} tx={tx} support={support:.0f}"
                    for obs, norm, tx, support in read_variants
                ),
                "current_emitted_variants": "; ".join(
                    f"derived={derived} exact={exact} total={total} paths={paths}"
                    for derived, exact, total, paths in emit_variants
                ),
            }
        )

    tsv_path = args.output_stem.with_suffix(".tsv")
    with tsv_path.open("w", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=[
                "first_donor",
                "backward_pattern_exact_end_rows",
                "backward_locus_exact_end_readpath_rows",
                "backward_exact_pred_cov_sum",
                "backward_exact_pred_n",
                "backward_exact_pred_ids",
                "current_read_variants",
                "current_emitted_variants",
            ],
            delimiter="\t",
        )
        writer.writeheader()
        writer.writerows(rows)

    md_path = args.output_stem.with_suffix(".md")
    with md_path.open("w") as handle:
        handle.write("## Locus stage-loss compare\n\n")
        handle.write(
            f"Target family: backward read-path `node0_start={args.backward_node0_start}`, backward pred `span_start={backward_pred_span_start}` vs current emitted `derived_start={args.current_derived_start}` with exact end `{args.target_end}`.\n\n"
        )
        handle.write("Main pattern:\n")
        handle.write("- old backward keeps strong exact-end support at the read-pattern, read-path, and prediction stages when the family is truly preserved early\n")
        handle.write("- current VG still carries matching donor families, but support drifts first inside read-family `observed_end -> normalized_end/tx_end`, then again at `tx_end -> derived_end`\n\n")
        for row in rows:
            handle.write(f"- donor `{row['first_donor']}`:\n")
            handle.write(
                f"  old exact-end read-pattern rows: `{row['backward_pattern_exact_end_rows']}`\n"
            )
            handle.write(
                f"  old exact-end preds: cov `{row['backward_exact_pred_cov_sum']}` across `{row['backward_exact_pred_n']}` preds ({row['backward_exact_pred_ids']})\n"
            )
            handle.write(f"  current read variants: {row['current_read_variants'] or '(none)'}\n")
            handle.write(f"  current emitted variants: {row['current_emitted_variants'] or '(none)'}\n")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
