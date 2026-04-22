#!/usr/bin/env python3
"""Build a bundle-scoped backward-vs-current terminal-variant oracle table.

This script is intentionally focused on the current parity workflow:

- backward oracle evidence from the deep StringTie-style trace log
- backward final outputs from the old-path harness
- current family/mode/selected evidence from the native audit trace

Rows are keyed by normalized `last_junction` for one terminal donor.
"""

from __future__ import annotations

import argparse
import csv
from collections import defaultdict
from dataclasses import dataclass
from pathlib import Path
from typing import DefaultDict, Dict, Iterable, List, Set


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--trace-log", required=True, help="Backward deep trace log.")
    parser.add_argument(
        "--current-trace",
        required=True,
        help="Current native terminal audit trace TSV.",
    )
    parser.add_argument(
        "--current-summary",
        required=True,
        help="Current terminal audit summary TSV.",
    )
    parser.add_argument(
        "--old-transcripts",
        required=True,
        help="Old-path harness transcripts TSV for the focused bundle.",
    )
    parser.add_argument(
        "--terminal-donor",
        required=True,
        help="Terminal donor coordinate to audit, e.g. 16621027.",
    )
    parser.add_argument(
        "--output-stem",
        required=True,
        help="Prefix for <stem>.summary.tsv and <stem>.variants.tsv",
    )
    return parser.parse_args()


def parse_int(text: str) -> int:
    try:
        return int(text)
    except (TypeError, ValueError):
        return 0


def parse_float(text: str) -> float:
    try:
        return float(text)
    except (TypeError, ValueError):
        return 0.0


def split_csv_field(text: str) -> List[str]:
    if not text:
        return []
    return [part.strip() for part in text.split(",") if part.strip()]


def normalize_old_junction(junc: str) -> str:
    donor, acceptor = junc.split("-")
    return f"{donor}-{int(acceptor) + 1}"


def last_junction_from_chain(chain: str, old_half_open: bool = False) -> str:
    parts = split_csv_field(chain)
    if not parts:
        return ""
    last = parts[-1]
    return normalize_old_junction(last) if old_half_open else last


def span_label(start: str, end: str) -> str:
    if not start or not end:
        return ""
    return f"{start}-{end}"


@dataclass
class CurrentAgg:
    family_input_rows: int = 0
    family_input_support_sum: int = 0
    mode_member_rows: int = 0
    mode_member_support_sum: int = 0
    selected_rows: int = 0
    selected_support_sum: int = 0
    selected_mode_total_reads_max: int = 0
    selected_reasons: Set[str] | None = None
    selected_boundary_chains: Set[str] | None = None

    def __post_init__(self) -> None:
        if self.selected_reasons is None:
            self.selected_reasons = set()
        if self.selected_boundary_chains is None:
            self.selected_boundary_chains = set()


def load_current_trace(path: Path) -> Dict[str, CurrentAgg]:
    out: Dict[str, CurrentAgg] = {}
    with path.open() as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            key = row["last_junction"]
            agg = out.setdefault(key, CurrentAgg())
            agg.family_input_rows += parse_int(row.get("family_input_rows", "0"))
            agg.family_input_support_sum += parse_int(row.get("family_input_support_sum", "0"))
            agg.mode_member_rows += parse_int(row.get("mode_member_rows", "0"))
            agg.mode_member_support_sum += parse_int(row.get("mode_member_support_sum", "0"))
            agg.selected_rows += parse_int(row.get("selected_rows", "0"))
            agg.selected_support_sum += parse_int(row.get("selected_support_sum", "0"))
            agg.selected_mode_total_reads_max = max(
                agg.selected_mode_total_reads_max,
                parse_int(row.get("selected_mode_total_reads", "0")),
            )
            reason = (row.get("selected_reason", "") or "").strip()
            if reason:
                agg.selected_reasons.add(reason)
            chain = (row.get("selected_boundary_chain", "") or "").strip()
            if chain:
                agg.selected_boundary_chains.add(chain)
    return out


def load_current_summary(path: Path) -> Dict[str, dict]:
    out: Dict[str, dict] = {}
    with path.open() as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            out[row["last_junction"]] = {
                "native_selected_ids": split_csv_field(row.get("native_selected_ids", "")),
                "native_selected_spans": split_csv_field(row.get("native_selected_spans", "")),
                "rlink_selected_ids": split_csv_field(row.get("rlink_selected_ids", "")),
                "rlink_selected_spans": split_csv_field(row.get("rlink_selected_spans", "")),
                "rlink_trace_final_fates": split_csv_field(row.get("rlink_trace_final_fates", "")),
            }
    return out


def load_old_transcripts(path: Path, donor: str) -> Dict[str, dict]:
    tx_ids: DefaultDict[str, List[str]] = defaultdict(list)
    spans: DefaultDict[str, List[str]] = defaultdict(list)
    coverages: DefaultDict[str, List[str]] = defaultdict(list)
    with path.open() as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            last = last_junction_from_chain(row.get("junction_chain", ""), old_half_open=True)
            if not last or not last.startswith(f"{donor}-"):
                continue
            tx_ids[last].append(row["tx_idx"])
            spans[last].append(span_label(row.get("start", ""), row.get("end", "")))
            cov = row.get("coverage", "")
            if cov:
                coverages[last].append(cov)
    out: Dict[str, dict] = {}
    for key in set(tx_ids) | set(spans) | set(coverages):
        out[key] = {
            "tx_ids": tx_ids[key],
            "spans": spans[key],
            "coverages": coverages[key],
        }
    return out


def scan_trace_log(path: Path, last_junctions: Iterable[str]) -> Dict[str, dict]:
    metrics: Dict[str, dict] = {
        key: {
            "build_bad_junc": 0,
            "build_ejunc_delete": 0,
            "good_junc_accept": 0,
            "jdec_reason_good": 0,
            "example_bad_junc_line": "",
            "example_accept_line": "",
        }
        for key in last_junctions
    }
    needles = list(metrics.keys())
    with path.open(errors="replace") as handle:
        for line_no, line in enumerate(handle, 1):
            text = line.rstrip("\n")
            for needle in needles:
                if needle not in text:
                    continue
                entry = metrics[needle]
                if "BAD_JUNC" in text:
                    entry["build_bad_junc"] += 1
                    if not entry["example_bad_junc_line"]:
                        entry["example_bad_junc_line"] = f"{line_no}:{text}"
                if "EJUNC_DELETE" in text:
                    entry["build_ejunc_delete"] += 1
                if "good_junc: ACCEPTED" in text:
                    entry["good_junc_accept"] += 1
                    if not entry["example_accept_line"]:
                        entry["example_accept_line"] = f"{line_no}:{text}"
                if "reason=GOOD" in text and "jdec[" in text:
                    entry["jdec_reason_good"] += 1
    return metrics


def determine_status(
    backward_supported: bool,
    current_selected_spans: Set[str],
    backward_selected_spans: Set[str],
    current: CurrentAgg,
) -> str:
    exact_overlap = bool(current_selected_spans & backward_selected_spans)
    if exact_overlap:
        return "selected_exact_overlap"
    if backward_supported and current.selected_rows == 0 and current.mode_member_rows > 0:
        return "missing_after_mode_selection"
    if backward_supported and current.selected_rows == 0 and current.family_input_rows > 0:
        return "missing_before_mode_selection"
    if backward_supported and current.family_input_rows == 0:
        return "absent_from_current_family_input"
    if current.selected_rows > 0 and not backward_selected_spans:
        return "current_only_selected"
    if current.selected_rows > 0 and backward_selected_spans and not exact_overlap:
        return "selected_but_not_exact_backward"
    if backward_supported:
        return "backward_supported_other"
    return "no_backward_oracle_signal"


def main() -> int:
    args = parse_args()
    donor = args.terminal_donor

    current_trace = load_current_trace(Path(args.current_trace))
    current_summary = load_current_summary(Path(args.current_summary))
    old_transcripts = load_old_transcripts(Path(args.old_transcripts), donor)

    last_junctions = sorted(
        {
            key
            for key in set(current_summary) | set(current_trace) | set(old_transcripts)
            if key.startswith(f"{donor}-")
        },
        key=lambda item: int(item.split("-")[1]),
    )

    trace_metrics = scan_trace_log(Path(args.trace_log), last_junctions)

    variant_rows: List[dict] = []
    status_counts: DefaultDict[str, int] = defaultdict(int)

    for key in last_junctions:
        current = current_trace.get(key, CurrentAgg())
        summary = current_summary.get(
            key,
            {
                "native_selected_ids": [],
                "native_selected_spans": [],
                "rlink_selected_ids": [],
                "rlink_selected_spans": [],
                "rlink_trace_final_fates": [],
            },
        )
        old = old_transcripts.get(key, {"tx_ids": [], "spans": [], "coverages": []})
        trace = trace_metrics.get(
            key,
            {
                "build_bad_junc": 0,
                "build_ejunc_delete": 0,
                "good_junc_accept": 0,
                "jdec_reason_good": 0,
                "example_bad_junc_line": "",
                "example_accept_line": "",
            },
        )

        current_selected_spans = set(summary["native_selected_spans"])
        backward_selected_spans = set(summary["rlink_selected_spans"])
        backward_supported = bool(
            backward_selected_spans
            or old["tx_ids"]
            or trace["good_junc_accept"]
            or trace["jdec_reason_good"]
        )
        status = determine_status(
            backward_supported,
            current_selected_spans,
            backward_selected_spans,
            current,
        )
        status_counts[status] += 1

        variant_rows.append(
            {
                "last_junction": key,
                "legacy_bad_junc_lines": trace["build_bad_junc"],
                "legacy_ejunc_delete_lines": trace["build_ejunc_delete"],
                "legacy_good_junc_accept_lines": trace["good_junc_accept"],
                "legacy_jdec_reason_good_lines": trace["jdec_reason_good"],
                "legacy_bad_junc_example": trace["example_bad_junc_line"],
                "legacy_accept_example": trace["example_accept_line"],
                "oldrust_tx_count": len(old["tx_ids"]),
                "oldrust_tx_ids": ",".join(old["tx_ids"]),
                "oldrust_tx_spans": ",".join(old["spans"]),
                "oldrust_coverages": ",".join(old["coverages"]),
                "rlink_selected_count": len(summary["rlink_selected_spans"]),
                "rlink_selected_ids": ",".join(summary["rlink_selected_ids"]),
                "rlink_selected_spans": ",".join(summary["rlink_selected_spans"]),
                "current_family_input_rows": current.family_input_rows,
                "current_family_input_support_sum": current.family_input_support_sum,
                "current_mode_member_rows": current.mode_member_rows,
                "current_mode_member_support_sum": current.mode_member_support_sum,
                "current_selected_rows": current.selected_rows,
                "current_selected_support_sum": current.selected_support_sum,
                "current_selected_mode_total_reads_max": current.selected_mode_total_reads_max,
                "current_selected_reasons": ",".join(sorted(current.selected_reasons)),
                "current_selected_boundary_chains": "|".join(
                    sorted(current.selected_boundary_chains)
                ),
                "current_selected_spans": ",".join(summary["native_selected_spans"]),
                "exact_span_overlap_count": len(current_selected_spans & backward_selected_spans),
                "backward_supported": int(backward_supported),
                "status": status,
            }
        )

    summary_rows = [
        {"metric": "terminal_donor", "value": donor},
        {"metric": "variants_total", "value": len(variant_rows)},
        {
            "metric": "backward_supported_variants",
            "value": sum(row["backward_supported"] for row in variant_rows),
        },
        {
            "metric": "exact_current_backward_overlap_variants",
            "value": sum(1 for row in variant_rows if row["exact_span_overlap_count"] > 0),
        },
        {
            "metric": "missing_after_mode_selection",
            "value": status_counts["missing_after_mode_selection"],
        },
        {
            "metric": "missing_before_mode_selection",
            "value": status_counts["missing_before_mode_selection"],
        },
        {
            "metric": "selected_but_not_exact_backward",
            "value": status_counts["selected_but_not_exact_backward"],
        },
        {
            "metric": "selected_exact_overlap",
            "value": status_counts["selected_exact_overlap"],
        },
    ]

    out_variants = Path(f"{args.output_stem}.variants.tsv")
    out_summary = Path(f"{args.output_stem}.summary.tsv")

    with out_variants.open("w", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=list(variant_rows[0].keys()) if variant_rows else ["last_junction"],
            delimiter="\t",
        )
        writer.writeheader()
        for row in variant_rows:
            writer.writerow(row)

    with out_summary.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=["metric", "value"], delimiter="\t")
        writer.writeheader()
        for row in summary_rows:
            writer.writerow(row)

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
