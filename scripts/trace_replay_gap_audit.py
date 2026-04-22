#!/usr/bin/env python3
"""Audit how close the deep trace is to a full backward-process replay.

This is intentionally conservative. It inventories major event families in the
deep trace, records whether the current staged explainer already reconstructs
them, and highlights which remaining families still block a "replay from the
trace alone" workflow.
"""

from __future__ import annotations

import argparse
import csv
import re
from collections import Counter
from pathlib import Path


EVENT_SPECS = [
    {
        "name": "merge_read_to_group",
        "pattern": re.compile(r"--- merge_read_to_group:"),
        "coverage": "partial",
        "criticality": "high",
        "note": "Deep trace has read->group lineage, but current explainer only summarizes graph-build rows rather than replaying group merges exactly.",
    },
    {
        "name": "create_bundle",
        "pattern": re.compile(r"--- create_bundle:"),
        "coverage": "covered",
        "criticality": "high",
        "note": "Bundle/bnode construction is represented in graph-build stage tables.",
    },
    {
        "name": "build_graphs_group_to_bnode",
        "pattern": re.compile(r"GROUP_TO_BNODE|ALLOC_BUNDLE2GRAPH|BUNDLE2GRAPH"),
        "coverage": "covered",
        "criticality": "high",
        "note": "Bundle->graph projection is represented in graph-build stage tables.",
    },
    {
        "name": "build_graphs_junction_decision",
        "pattern": re.compile(r"\bjdec\[|good_junc: ACCEPTED|JUNC_COLOR_BREAK"),
        "coverage": "covered",
        "criticality": "high",
        "note": "Current explainer writes structured junction-trace tables.",
    },
    {
        "name": "build_graphs_junction_mutation",
        "pattern": re.compile(r"JUNC_DELETE|JUNC_DEMOTE"),
        "coverage": "covered",
        "criticality": "high",
        "note": "Junction deletion/demotion is staged explicitly into a bundle-scoped junction-mutation table.",
    },
    {
        "name": "build_graphs_routing",
        "pattern": re.compile(r"ROUTE_STRANDED|ROUTE_UNSTRANDED"),
        "coverage": "covered",
        "criticality": "medium",
        "note": "Read routing by strand/group is staged explicitly into bundle-scoped route tables.",
    },
    {
        "name": "path_update_abund_core",
        "pattern": re.compile(r"PATH_update_abund: (LR_ENTRY|LR_FLAGS|FINAL_NODES|NEW_TRF|ADD_TO_TRF|SET_LONGSTART|SET_LONGEND|SKIP_SINGLE_NODE)"),
        "coverage": "covered",
        "criticality": "high",
        "note": "Read->path->raw-transfrag lineage is already reconstructed into staged lineage tables.",
    },
    {
        "name": "path_update_abund_boundary_logic",
        "pattern": re.compile(r"PATH_update_abund: (SINK_CHECK_ENTER|SOURCE_CHECK_ENTER|SINK_PROXIMITY|SOURCE_PROXIMITY|SINK_TRIM|SOURCE_TRIM)"),
        "coverage": "partial",
        "criticality": "high",
        "note": "Boundary-adjustment events exist in the deep trace, but current lineage only preserves a reduced summary of terminal evidence.",
    },
    {
        "name": "path_read_pattern",
        "pattern": re.compile(r"PATH_read_pattern:"),
        "coverage": "partial",
        "criticality": "high",
        "note": "Current explainer stages raw PATH_read_pattern rows and aggregated read-pattern path summaries, but it does not yet replay them into mutable graph state.",
    },
    {
        "name": "longtrim",
        "pattern": re.compile(r"LONGTRIM_(BOUND|BNODE|SPLIT)"),
        "coverage": "partial",
        "criticality": "medium",
        "note": "Longtrim rows are counted and partly staged, but not replayed as a mutable graph operation sequence.",
    },
    {
        "name": "create_graph_nodes",
        "pattern": re.compile(r"--- create_graph: (NEW_NODE|NODE_SHRINK_JSTART|NODE_SHRINK_JEND|NODE_FINAL_END|SRCSINK_COVBUILD_SOURCE)"),
        "coverage": "partial",
        "criticality": "high",
        "note": "Graph mutation rows are staged, but not yet replayed into a faithful per-graph object history.",
    },
    {
        "name": "process_transfrags",
        "pattern": re.compile(r"--- process_transfrags:"),
        "coverage": "partial",
        "criticality": "high",
        "note": "SRCSINK decisions are extracted, but eliminate/transform causality is not yet fully replayable.",
    },
    {
        "name": "parse_trflong_seed_core",
        "pattern": re.compile(r"--- parse_trflong: (SEED |SEED_PATHPAT|back_to_source|fwd_to_sink|FLOW_PATHPAT|long_max_flow|DEPLETION_BEFORE|DEPLETION_AFTER|SEED_DECISION)"),
        "coverage": "covered",
        "criticality": "high",
        "note": "Seed lineage and decision tables already reconstruct most of this stage.",
    },
    {
        "name": "parse_trflong_boundary_tiebreak",
        "pattern": re.compile(
            r"--- parse_trflong: HARD_BOUNDARY_CHECK|--- (fwd_to_sink|back_to_source): step|--- flow: (PATHPAT_OR|TIEBREAK)|TIEBREAK_(BACK|FWD)"
        ),
        "coverage": "partial",
        "criticality": "high",
        "note": "Low-level step traces are now staged explicitly, and newer parity traces can also supply `TIEBREAK_BACK` / `TIEBREAK_FWD`, but exact HARD_BOUNDARY_CHECK / PATHPAT_OR / TIEBREAK replay still depends on what is present in the trace set.",
    },
    {
        "name": "checktrf",
        "pattern": re.compile(r"--- checktrf:"),
        "coverage": "missing",
        "criticality": "high",
        "note": "Deep trace API has explicit checktrf gate/rescue/complete/reject events, but current replay still relies on side traces rather than deep-log reconstruction.",
    },
    {
        "name": "predcluster_pred_dump",
        "pattern": re.compile(r"--- print_predcluster: pred\["),
        "coverage": "covered",
        "criticality": "high",
        "note": "Prediction tables are reconstructed and linked into final fate edges.",
    },
    {
        "name": "predcluster_filter_fate",
        "pattern": re.compile(r"--- print_predcluster: (ENTER|BEFORE_PAIRWISE|FILTER|PRED_FATE|FINAL_KEEP|FINAL_DROP|FINAL_FATE|FINAL_SURVIVORS|FINAL_LEDGER|PRED_FATE_SUMMARY)"),
        "coverage": "covered",
        "criticality": "high",
        "note": "Prediction edge/fate layers are already reconstructed well.",
    },
    {
        "name": "guide_explicit",
        "pattern": re.compile(r"\bguide="),
        "coverage": "partial",
        "criticality": "medium",
        "note": "Guide metadata is carried through many rows, but guide-driven causal decisions are not yet isolated as their own replay layer.",
    },
]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--trace-log", required=True, type=Path)
    parser.add_argument("--parity-trace-log", type=Path)
    parser.add_argument("--output-stem", required=True, type=Path)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    text = args.trace_log.read_text(errors="replace")
    lines = text.splitlines()
    parity_text = args.parity_trace_log.read_text(errors="replace") if args.parity_trace_log else ""
    parity_lines = parity_text.splitlines() if parity_text else []

    rows: list[dict[str, str]] = []
    coverage_counts = Counter()
    critical_missing: list[str] = []
    deep_trace_only_blockers: list[str] = []

    for spec in EVENT_SPECS:
        match = spec["pattern"].search(text)
        count = sum(1 for line in lines if spec["pattern"].search(line))
        coverage = spec["coverage"]
        note = spec["note"]
        if spec["name"] == "checktrf" and args.parity_trace_log:
            parity_count = sum(
                1
                for line in parity_lines
                if "PARITY_SEED_PROC " in line
                or "PARITY_CHECKTRF_START" in line
                or "PARITY_CHK " in line
                or "PARITY_CHECKTRF_OUTCOME" in line
            )
            if parity_count > 0:
                count += parity_count
                match = match or True
                coverage = "partial"
                note = (
                    "Deep trace still lacks explicit `--- checktrf:` rows, but `trace_full.log` provides "
                    "`PARITY_SEED_PROC` plus `PARITY_CHK`"
                    + (" and legacy batch/outcome markers" if any("PARITY_CHECKTRF_START" in line or "PARITY_CHECKTRF_OUTCOME" in line for line in parity_lines) else "")
                    + ", which is enough to reconstruct a first-class parity checktrf layer."
                )
        example = ""
        if match:
            for idx, line in enumerate(lines, 1):
                if spec["pattern"].search(line):
                    example = f"{idx}:{line}"
                    break
            if not example and spec["name"] == "checktrf" and parity_lines:
                for idx, line in enumerate(parity_lines, 1):
                    if (
                        "PARITY_SEED_PROC " in line
                        or "PARITY_CHECKTRF_START" in line
                        or "PARITY_CHK " in line
                        or "PARITY_CHECKTRF_OUTCOME" in line
                    ):
                        example = f"{idx}:{line}"
                        break
        coverage_counts[coverage] += 1
        if coverage in {"missing", "partial"} and spec["criticality"] == "high":
            critical_missing.append(spec["name"])
        if coverage == "missing" and spec["criticality"] == "high":
            deep_trace_only_blockers.append(spec["name"])
        rows.append(
            {
                "event_family": spec["name"],
                "coverage": coverage,
                "criticality": spec["criticality"],
                "count": str(count),
                "present_in_trace": "1" if count > 0 else "0",
                "example": example,
                "note": note,
            }
        )

    tsv_path = args.output_stem.with_suffix(".tsv")
    with tsv_path.open("w", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=["event_family", "coverage", "criticality", "count", "present_in_trace", "example", "note"],
            delimiter="\t",
        )
        writer.writeheader()
        writer.writerows(rows)

    md_path = args.output_stem.with_suffix(".md")
    with md_path.open("w") as handle:
        handle.write("# Deep Trace Replay Gap Audit\n\n")
        handle.write(f"- Trace: `{args.trace_log.name}`\n")
        if args.parity_trace_log:
            handle.write(f"- Parity trace: `{args.parity_trace_log.name}`\n")
        handle.write(f"- Event families audited: `{len(EVENT_SPECS)}`\n")
        handle.write(f"- Covered: `{coverage_counts['covered']}`\n")
        handle.write(f"- Partial: `{coverage_counts['partial']}`\n")
        handle.write(f"- Missing: `{coverage_counts['missing']}`\n\n")
        handle.write("## Main Result\n\n")
        if deep_trace_only_blockers:
            handle.write(
                "- A replay from the deep trace alone is **not yet supported** by the current machinery.\n"
            )
            handle.write(
                "- The highest-signal missing high-criticality families are: "
                + ", ".join(f"`{name}`" for name in deep_trace_only_blockers)
                + ".\n"
            )
        else:
            handle.write("- No high-criticality missing families were detected in this audit.\n")
        if critical_missing:
            handle.write(
                "- High-criticality families that still need full replay coverage or stronger normalization: "
                + ", ".join(f"`{name}`" for name in critical_missing)
                + ".\n"
            )
        handle.write("\n## Recommendations\n\n")
        handle.write("- Add explicit reconstruction for `PATH_read_pattern`.\n")
        handle.write("- Add low-level flow replay for `HARD_BOUNDARY_CHECK`, step traces, `PATHPAT_OR`, and both legacy `TIEBREAK` plus newer `TIEBREAK_BACK` / `TIEBREAK_FWD` rows.\n")
        handle.write("- Reconstruct `checktrf` directly from the deep trace rather than side traces.\n")
        handle.write("- Promote `JUNC_DELETE` / `JUNC_DEMOTE` / route events into structured stage tables.\n")
        handle.write("- Convert graph mutation rows into a real mutable replay object if exact graph recreation is the goal.\n")
        handle.write("- Keep intentional Iso-Seq-style representative collapse separate from trace-missing blockers when interpreting late-emission losses.\n")


if __name__ == "__main__":
    main()
