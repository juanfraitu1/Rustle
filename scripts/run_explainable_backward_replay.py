#!/usr/bin/env python3
"""Run backward replay only on bundles with currently explainable lineage.

This intentionally excludes known guard bundles whose current trace slices do
not recover enough read->path lineage for a useful backward replay package.
At the moment that means:

- include: bundle1014 / STRG.171
- include: bundle1074 / STRG.212
- exclude for this replay pass: bundle70 / STRG.88, bundle904 / STRG.15
"""

from __future__ import annotations

import argparse
import csv
import subprocess
import sys
from pathlib import Path


EXPLAINABLE_BUNDLES = [
    {
        "bundle_id": "1014",
        "gene_id": "STRG.171",
        "bundle_summary": "GGO_19_bundle1014_familytrace.bundles.tsv",
        "legacy_trace_norm": "GGO_19_bundle1014_legacy_trace_norm.tsv",
        "rust_trace_norm": "GGO_19_bundle1014_oldrust_trace_strg171_norm.tsv",
        "old_transcripts": "GGO_19_bundle1014_oldpath_trace_strg171.transcripts.tsv",
    },
    {
        "bundle_id": "1074",
        "gene_id": "STRG.212",
        "bundle_summary": "GGO_19_bundle1074_familytrace.bundles.tsv",
        "legacy_trace_norm": "GGO_19_bundle1074_legacy_trace_norm.tsv",
        "rust_trace_norm": "GGO_19_bundle1074_oldrust_trace_strg212_norm.tsv",
        "old_transcripts": "GGO_19_bundle1074_oldpath_trace_strg212.transcripts.tsv",
    },
]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=Path("/tmp/explainable_backward_replay"),
        help="Directory for regenerated explainable-only backward replay outputs.",
    )
    return parser.parse_args()


def read_stage_counts(path: Path) -> dict[str, int]:
    out: dict[str, int] = {}
    with path.open() as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            try:
                out[row["stage"]] = int(row["count"])
            except (TypeError, ValueError):
                out[row["stage"]] = 0
    return out


def read_layer_counts(path: Path) -> dict[str, int]:
    out: dict[str, int] = {}
    with path.open() as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            try:
                out[row["layer_id"]] = int(row["count"])
            except (TypeError, ValueError):
                out[row["layer_id"]] = 0
    return out


def main() -> int:
    args = parse_args()
    repo_root = Path(__file__).resolve().parent.parent
    output_dir = args.output_dir.resolve()
    output_dir.mkdir(parents=True, exist_ok=True)

    explain_script = repo_root / "scripts" / "explain_bundle_stage_graphs.py"
    bundle_log = repo_root / "GGO_19.log"
    reference_gtf = repo_root / "GGO_19.gtf"
    legacy_trace_log = repo_root / "trace_GGO_19_deep.log"
    parity_trace_log = repo_root / "trace_full.log"

    summary_rows: list[dict[str, str]] = []

    for spec in EXPLAINABLE_BUNDLES:
        stem = output_dir / f"bundle{spec['bundle_id']}_stage_explain"
        cmd = [
            sys.executable,
            str(explain_script),
            "--bundle-summary",
            str(repo_root / spec["bundle_summary"]),
            "--bundle-id",
            spec["bundle_id"],
            "--bundle-log",
            str(bundle_log),
            "--reference-gtf",
            str(reference_gtf),
            "--legacy-trace-log",
            str(legacy_trace_log),
            "--parity-trace-log",
            str(parity_trace_log),
            "--legacy-trace-norm",
            str(repo_root / spec["legacy_trace_norm"]),
            "--rust-trace-norm",
            str(repo_root / spec["rust_trace_norm"]),
            "--old-transcripts",
            str(repo_root / spec["old_transcripts"]),
            "--output-stem",
            str(stem),
        ]
        subprocess.run(cmd, cwd=repo_root, check=True)

        stages = read_stage_counts(stem.with_suffix(".stages.tsv"))
        layers = read_layer_counts(stem.with_suffix(".layers.tsv"))
        summary_rows.append(
            {
                "bundle_id": spec["bundle_id"],
                "gene_id": spec["gene_id"],
                "status": "included_explainable",
                "focused_output_stem": str(stem),
                "reference_overlap": str(stages.get("reference_overlap", 0)),
                "rawgraph_transfrags": str(stages.get("rawgraph_transfrags", 0)),
                "checktrf_seeds": str(stages.get("checktrf_seeds", 0)),
                "checktrf_outputs": str(stages.get("checktrf_outputs", 0)),
                "lineage_read_events": str(stages.get("lineage_read_events", 0)),
                "lineage_read_graph_paths": str(stages.get("lineage_read_graph_paths", 0)),
                "lineage_read_to_rawtrf": str(stages.get("lineage_read_to_rawtrf", 0)),
                "lineage_rawtrf_provenance": str(stages.get("lineage_rawtrf_provenance", 0)),
                "lineage_seed_decisions": str(stages.get("lineage_seed_decisions", 0)),
                "lineage_gate_rows": str(stages.get("lineage_gate_rows", 0)),
                "lineage_junction_rows": str(stages.get("lineage_junction_rows", 0)),
                "lineage_parity_graph_rows": str(stages.get("lineage_parity_graph_rows", 0)),
                "lineage_parity_boundary_check_rows": str(stages.get("lineage_parity_boundary_check_rows", 0)),
                "lineage_parity_checktrf_event_rows": str(stages.get("lineage_parity_checktrf_event_rows", 0)),
                "L0_bundle_graph_context": str(layers.get("L0_bundle_graph_context", 0)),
                "L1_junction_policy": str(layers.get("L1_junction_policy", 0)),
                "L2_path_lineage": str(layers.get("L2_path_lineage", 0)),
                "L3_raw_transfrag_state": str(layers.get("L3_raw_transfrag_state", 0)),
                "L4_seed_flow": str(layers.get("L4_seed_flow", 0)),
                "L5_checktrf_rescue": str(layers.get("L5_checktrf_rescue", 0)),
                "L6_prediction_fate": str(layers.get("L6_prediction_fate", 0)),
                "L7_transcript_emission": str(layers.get("L7_transcript_emission", 0)),
            }
        )

    summary_path = output_dir / "summary.tsv"
    with summary_path.open("w", newline="") as handle:
        fieldnames = list(summary_rows[0].keys()) if summary_rows else []
        writer = csv.DictWriter(handle, delimiter="\t", fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(summary_rows)

    md_path = output_dir / "summary.md"
    with md_path.open("w") as handle:
        handle.write("# Explainable Backward Replay\n\n")
        handle.write("- Included bundles: `1014` / `STRG.171`, `1074` / `STRG.212`\n")
        handle.write("- Guard bundles intentionally excluded for now: `70` / `STRG.88`, `904` / `STRG.15`\n")
        handle.write("- Reason: current data already expose usable read->graph->rawtrf->seed provenance on `1014` and `1074`, while `70` and `904` remain weaker drivers for a clean backward replay pass.\n\n")
        handle.write("## Bundle Summary\n\n")
        for row in summary_rows:
            handle.write(
                f"- bundle `{row['bundle_id']}` / `{row['gene_id']}`:"
                f" path lineage `{row['L2_path_lineage']}`, raw-transfrag state `{row['L3_raw_transfrag_state']}`,"
                f" seed flow `{row['L4_seed_flow']}`, checktrf rescue `{row['L5_checktrf_rescue']}`,"
                f" transcript emission `{row['L7_transcript_emission']}`,"
                f" parity graph `{row['lineage_parity_graph_rows']}`,"
                f" parity boundary checks `{row['lineage_parity_boundary_check_rows']}`.\n"
            )

    print(summary_path)
    print(md_path)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
