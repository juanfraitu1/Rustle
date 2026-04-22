#!/usr/bin/env python3
"""Aggregate blocker analysis for the current explainable backward-replay subset."""

from __future__ import annotations

import argparse
import csv
from collections import Counter, defaultdict
from pathlib import Path


FOCUSED = [
    {
        "gene_id": "STRG.171",
        "oracle_refs": "GGO_19_STRG171_backward_oracle_live.refs.tsv",
        "exact_classes": "GGO_19_STRG171_exact_boundary_classes.refs.tsv",
    },
    {
        "gene_id": "STRG.212",
        "oracle_refs": "GGO_19_STRG212_backward_oracle_live.refs.tsv",
        "exact_classes": "GGO_19_STRG212_exact_boundary_classes.refs.tsv",
        "ranking_summary": "GGO_19_STRG212_ranking_oracle.summary.tsv",
    },
]

LATE_EMISSION_CATEGORIES = {
    "path_pair_present_but_emission_shifted",
    "path_has_ref_pair_but_emission_loses_it",
}

EARLY_PATH_CATEGORIES = {
    "raw_has_ref_pair_but_path_never_keeps_it",
    "raw_exact_chain_only",
}

REPRESENTATIVE_COLLAPSE_CATEGORIES = {
    "stable_nonref_same_chain_pair",
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=Path("/tmp/explainable_backward_blockers"),
        help="Directory for explainable-only blocker summary outputs.",
    )
    parser.add_argument("--strg171-exact-classes", type=Path, help="Override exact-boundary classes TSV for STRG.171.")
    parser.add_argument("--strg212-exact-classes", type=Path, help="Override exact-boundary classes TSV for STRG.212.")
    parser.add_argument("--strg212-ranking-summary", type=Path, help="Override ranking-summary TSV for STRG.212.")
    parser.add_argument(
        "--replay-dir",
        type=Path,
        help="Optional explainable backward replay output dir with per-bundle stage-explain tables.",
    )
    return parser.parse_args()


def load_tsv(path: Path) -> list[dict[str, str]]:
    with path.open() as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def load_metric_map(path: Path) -> dict[str, str]:
    out: dict[str, str] = {}
    with path.open() as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            out[row["metric"]] = row["count"]
    return out


def parse_boundary_pair(text: str) -> tuple[int, int]:
    try:
        start, end = text.split("-", 1)
        return int(start), int(end)
    except (AttributeError, TypeError, ValueError):
        return 0, 0


def boundary_pairs_match(left: str, right: str, tolerance: int = 10) -> bool:
    left_start, left_end = parse_boundary_pair(left)
    right_start, right_end = parse_boundary_pair(right)
    if not left_start or not left_end or not right_start or not right_end:
        return False
    return abs(left_start - right_start) <= tolerance and abs(left_end - right_end) <= tolerance


def load_replay_bundle_tables(replay_dir: Path | None) -> dict[str, dict[str, list[dict[str, str]]]]:
    if replay_dir is None:
        return {}

    summary_path = replay_dir / "summary.tsv"
    if not summary_path.exists():
        return {}

    out: dict[str, dict[str, list[dict[str, str]]]] = {}
    with summary_path.open() as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            gene_id = row["gene_id"]
            stem = Path(row["focused_output_stem"])
            parity_boundary_path = stem.with_suffix(".lineage.parity_boundary_checks.tsv")
            flow_trace_path = stem.with_suffix(".lineage.flow_trace.tsv")
            out[gene_id] = {
                "parity_boundary_rows": load_tsv(parity_boundary_path) if parity_boundary_path.exists() else [],
                "flow_trace_rows": load_tsv(flow_trace_path) if flow_trace_path.exists() else [],
            }
    return out


def classify_late_emission_subtype(
    category: str,
    best_path_boundary: str,
    best_emitted_boundary: str,
    parity_boundary_rows: list[dict[str, str]],
    flow_trace_rows: list[dict[str, str]],
) -> tuple[str, dict[str, str]]:
    path_boundary_matches = [
        row for row in parity_boundary_rows if boundary_pairs_match(best_path_boundary, f"{row.get('span_start', '')}-{row.get('span_end', '')}")
    ]
    emitted_boundary_matches = [
        row for row in parity_boundary_rows if boundary_pairs_match(best_emitted_boundary, f"{row.get('span_start', '')}-{row.get('span_end', '')}")
    ]
    tiebreak_rows = [row for row in flow_trace_rows if row.get("event") in {"TIEBREAK_BACK", "TIEBREAK_FWD"}]
    nodecov_tiebreak_rows = [row for row in tiebreak_rows if row.get("reason") == "nodecov_tiebreak"]

    if not path_boundary_matches:
        subtype = "boundary_gate_or_seed_flow"
    elif path_boundary_matches and emitted_boundary_matches and best_path_boundary != best_emitted_boundary:
        subtype = "final_emission_choice_after_boundary_gate"
    elif category == "path_has_ref_pair_but_emission_loses_it":
        subtype = "final_emission_loss_after_boundary_gate"
    elif tiebreak_rows:
        subtype = "final_tiebreak_or_ranking"
    else:
        subtype = "final_emission_after_boundary_gate"

    context = {
        "late_subtype": subtype,
        "path_parity_boundary_match_count": str(len(path_boundary_matches)),
        "emitted_parity_boundary_match_count": str(len(emitted_boundary_matches)),
        "path_parity_hardstart_rows": str(sum(1 for row in path_boundary_matches if row.get("hardstart") == "1")),
        "path_parity_hardend_rows": str(sum(1 for row in path_boundary_matches if row.get("hardend") == "1")),
        "emitted_parity_hardstart_rows": str(sum(1 for row in emitted_boundary_matches if row.get("hardstart") == "1")),
        "emitted_parity_hardend_rows": str(sum(1 for row in emitted_boundary_matches if row.get("hardend") == "1")),
        "bundle_tiebreak_rows": str(len(tiebreak_rows)),
        "bundle_nodecov_tiebreak_rows": str(len(nodecov_tiebreak_rows)),
    }
    return subtype, context


def main() -> int:
    args = parse_args()
    repo_root = Path(__file__).resolve().parent.parent
    output_dir = args.output_dir.resolve()
    output_dir.mkdir(parents=True, exist_ok=True)

    exact_class_overrides = {
        "STRG.171": args.strg171_exact_classes,
        "STRG.212": args.strg212_exact_classes,
    }
    ranking_overrides = {
        "STRG.212": args.strg212_ranking_summary,
    }
    replay_tables = load_replay_bundle_tables(args.replay_dir)

    per_ref_rows: list[dict[str, str]] = []
    per_gene_rows: list[dict[str, str]] = []

    verdict_counts = Counter()
    category_counts = Counter()
    stage_bucket_counts = Counter()
    late_subtype_counts = Counter()

    for spec in FOCUSED:
        oracle_rows = load_tsv(repo_root / spec["oracle_refs"])
        exact_classes_path = exact_class_overrides.get(spec["gene_id"]) or (repo_root / spec["exact_classes"])
        class_rows = {
            row["ref_id"]: row
            for row in load_tsv(exact_classes_path)
            if row["ref_id"].startswith(f"{spec['gene_id']}.")
        }
        ranking_path = ranking_overrides.get(spec["gene_id"]) or (repo_root / spec["ranking_summary"]) if spec.get("ranking_summary") else None
        ranking_metrics = load_metric_map(ranking_path) if ranking_path else {}

        gene_verdicts = Counter()
        gene_categories = Counter()
        gene_stage_buckets = Counter()
        gene_late_subtypes = Counter()
        gene_replay = replay_tables.get(spec["gene_id"], {})
        parity_boundary_rows = gene_replay.get("parity_boundary_rows", [])
        flow_trace_rows = gene_replay.get("flow_trace_rows", [])

        for row in oracle_rows:
            ref_id = row["ref_id"]
            cls = class_rows.get(ref_id, {})
            verdict = row.get("oracle_verdict", "")
            category = cls.get("category", "")
            note = cls.get("note", "")

            if category in LATE_EMISSION_CATEGORIES:
                stage_bucket = "late_emission"
            elif category in EARLY_PATH_CATEGORIES:
                stage_bucket = "early_path_or_partition"
            elif category in REPRESENTATIVE_COLLAPSE_CATEGORIES:
                stage_bucket = "representative_collapse"
            elif category == "non_exact_chain_driver":
                stage_bucket = "non_exact_chain_driver"
            else:
                stage_bucket = "unclassified"

            late_subtype = ""
            late_context = {
                "path_parity_boundary_match_count": "",
                "emitted_parity_boundary_match_count": "",
                "path_parity_hardstart_rows": "",
                "path_parity_hardend_rows": "",
                "emitted_parity_hardstart_rows": "",
                "emitted_parity_hardend_rows": "",
                "bundle_tiebreak_rows": "",
                "bundle_nodecov_tiebreak_rows": "",
            }
            if stage_bucket == "late_emission":
                late_subtype, late_context = classify_late_emission_subtype(
                    category,
                    cls.get("best_path_boundary", ""),
                    cls.get("best_emitted_boundary", ""),
                    parity_boundary_rows,
                    flow_trace_rows,
                )

            per_ref_rows.append(
                {
                    "gene_id": spec["gene_id"],
                    "ref_id": ref_id,
                    "oracle_verdict": verdict,
                    "current_classification": row.get("current_classification", ""),
                    "exact_boundary_category": category,
                    "exact_boundary_note": note,
                    "stage_bucket": stage_bucket,
                    "late_subtype": late_subtype,
                    "raw_has_ref_boundary": cls.get("raw_has_ref_boundary", ""),
                    "path_has_ref_boundary": cls.get("path_has_ref_boundary", ""),
                    "emitted_has_ref_boundary": cls.get("emitted_has_ref_boundary", ""),
                    "raw_exact_boundary_count": cls.get("raw_exact_boundary_count", ""),
                    "path_exact_candidate_count": cls.get("path_exact_candidate_count", ""),
                    "emitted_same_chain_count": cls.get("emitted_same_chain_count", ""),
                    "best_path_boundary": cls.get("best_path_boundary", ""),
                    "best_emitted_boundary": cls.get("best_emitted_boundary", ""),
                    **late_context,
                }
            )

            verdict_counts[verdict] += 1
            gene_verdicts[verdict] += 1
            category_counts[category] += 1
            gene_categories[category] += 1
            stage_bucket_counts[stage_bucket] += 1
            gene_stage_buckets[stage_bucket] += 1
            if late_subtype:
                late_subtype_counts[late_subtype] += 1
                gene_late_subtypes[late_subtype] += 1

        per_gene_rows.append(
            {
                "gene_id": spec["gene_id"],
                "focused_refs": str(len(oracle_rows)),
                "boundary_only_current_miss": str(gene_verdicts.get("boundary_only_current_miss", 0)),
                "canonical_ref_family_preferred_backward": str(gene_verdicts.get("canonical_ref_family_preferred_backward", 0)),
                "backward_keeps_exact_span": str(gene_verdicts.get("backward_keeps_exact_span", 0)),
                "late_emission": str(gene_stage_buckets.get("late_emission", 0)),
                "early_path_or_partition": str(gene_stage_buckets.get("early_path_or_partition", 0)),
                "representative_collapse": str(gene_stage_buckets.get("representative_collapse", 0)),
                "non_exact_chain_driver": str(gene_stage_buckets.get("non_exact_chain_driver", 0)),
                "late_boundary_gate_or_seed_flow": str(gene_late_subtypes.get("boundary_gate_or_seed_flow", 0)),
                "late_final_emission_choice_after_boundary_gate": str(gene_late_subtypes.get("final_emission_choice_after_boundary_gate", 0)),
                "late_final_emission_loss_after_boundary_gate": str(gene_late_subtypes.get("final_emission_loss_after_boundary_gate", 0)),
                "late_final_tiebreak_or_ranking": str(gene_late_subtypes.get("final_tiebreak_or_ranking", 0)),
                "late_final_emission_after_boundary_gate": str(gene_late_subtypes.get("final_emission_after_boundary_gate", 0)),
                "heuristic_prefer_rows": ranking_metrics.get("heuristic_prefer_rows", ""),
                "heuristic_prefer_refs": ranking_metrics.get("heuristic_prefer_refs", ""),
                "ref_boundary_exact_full_rows": ranking_metrics.get("ref_boundary_exact_full_rows", ""),
                "ref_boundary_exact_full_refs": ranking_metrics.get("ref_boundary_exact_full_refs", ""),
                "ref_boundary_raw_rows": ranking_metrics.get("ref_boundary_raw_rows", ""),
                "ref_boundary_raw_refs": ranking_metrics.get("ref_boundary_raw_refs", ""),
            }
        )

    summary_rows = [
        {"metric": "focused_genes", "value": str(len(FOCUSED))},
        {"metric": "focused_refs_total", "value": str(sum(int(row["focused_refs"]) for row in per_gene_rows))},
    ]
    for key, value in sorted(verdict_counts.items()):
        summary_rows.append({"metric": f"oracle_verdict:{key}", "value": str(value)})
    for key, value in sorted(category_counts.items()):
        summary_rows.append({"metric": f"exact_boundary_category:{key}", "value": str(value)})
    for key, value in sorted(stage_bucket_counts.items()):
        summary_rows.append({"metric": f"stage_bucket:{key}", "value": str(value)})
    for key, value in sorted(late_subtype_counts.items()):
        summary_rows.append({"metric": f"late_subtype:{key}", "value": str(value)})

    per_ref_path = output_dir / "per_ref.tsv"
    with per_ref_path.open("w", newline="") as handle:
        fieldnames = list(per_ref_rows[0].keys()) if per_ref_rows else []
        writer = csv.DictWriter(handle, delimiter="\t", fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(per_ref_rows)

    per_gene_path = output_dir / "per_gene.tsv"
    with per_gene_path.open("w", newline="") as handle:
        fieldnames = list(per_gene_rows[0].keys()) if per_gene_rows else []
        writer = csv.DictWriter(handle, delimiter="\t", fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(per_gene_rows)

    summary_path = output_dir / "summary.tsv"
    with summary_path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, delimiter="\t", fieldnames=["metric", "value"])
        writer.writeheader()
        writer.writerows(summary_rows)

    md_path = output_dir / "summary.md"
    with md_path.open("w") as handle:
        handle.write("# Explainable-Only Backward Blockers\n\n")
        handle.write("- Focused loci: `STRG.171`, `STRG.212`\n")
        handle.write("- These are the current clean backward-replay calibration loci.\n")
        handle.write("- Guard loci such as `STRG.88` / bundle `70` and `STRG.15` / bundle `904` are intentionally excluded.\n\n")
        handle.write("## Aggregate\n\n")
        handle.write(f"- Focused refs: `{sum(int(row['focused_refs']) for row in per_gene_rows)}`\n")
        handle.write(f"- Late-emission blockers: `{stage_bucket_counts.get('late_emission', 0)}`\n")
        handle.write(f"- Early path / partition blockers: `{stage_bucket_counts.get('early_path_or_partition', 0)}`\n")
        handle.write(f"- Representative-collapse blockers: `{stage_bucket_counts.get('representative_collapse', 0)}`\n")
        handle.write(f"- Non-exact-chain drivers: `{stage_bucket_counts.get('non_exact_chain_driver', 0)}`\n")
        if late_subtype_counts:
            handle.write("\n## Late-Emission Split\n\n")
            for key, value in sorted(late_subtype_counts.items()):
                handle.write(f"- `{key}`: `{value}`\n")
        handle.write("\n## Oracle Verdicts\n\n")
        for key, value in sorted(verdict_counts.items()):
            handle.write(f"- `{key}`: `{value}`\n")
        handle.write("\n## Per Gene\n\n")
        for row in per_gene_rows:
            handle.write(
                f"- `{row['gene_id']}`: refs `{row['focused_refs']}`,"
                f" late-emission `{row['late_emission']}`,"
                f" early path/partition `{row['early_path_or_partition']}`,"
                f" representative-collapse `{row['representative_collapse']}`,"
                f" non-exact-chain `{row['non_exact_chain_driver']}`"
            )
            if row["late_boundary_gate_or_seed_flow"] or row["late_final_emission_choice_after_boundary_gate"] or row["late_final_emission_loss_after_boundary_gate"] or row["late_final_tiebreak_or_ranking"] or row["late_final_emission_after_boundary_gate"]:
                handle.write(
                    f"; late split gate/seed `{row['late_boundary_gate_or_seed_flow'] or '0'}`,"
                    f" final-choice `{row['late_final_emission_choice_after_boundary_gate'] or '0'}`,"
                    f" final-loss `{row['late_final_emission_loss_after_boundary_gate'] or '0'}`,"
                    f" tiebreak/ranking `{row['late_final_tiebreak_or_ranking'] or '0'}`,"
                    f" final-after-gate `{row['late_final_emission_after_boundary_gate'] or '0'}`"
                )
            if row["heuristic_prefer_rows"]:
                handle.write(
                    f", heuristic-prefer exact-boundary rows `{row['heuristic_prefer_rows']}`"
                )
            if row["heuristic_prefer_refs"]:
                handle.write(
                    f", heuristic-prefer refs `{row['heuristic_prefer_refs']}`"
                )
            handle.write(".\n")
        handle.write("\n## Interpretation\n\n")
        handle.write(
            "- The explainable subset is dominated by stable same-chain representative selection, not by generic late-emission loss.\n"
        )
        handle.write(
            "- `STRG.171` is now a clean representative-collapse locus: the path family and emitted family usually agree on the same non-reference boundary pair.\n"
        )
        handle.write(
            "- `STRG.212` still mixes representative collapse with a smaller earlier path/partition failure set where exact-boundary evidence exists in raw support but does not survive as a kept path.\n"
        )
        if late_subtype_counts:
            handle.write(
                "- With parity boundary-check rows joined in, late-emission blockers can now be split into cases that already survive the boundary-gate stage versus cases that still look boundary-gate or seed-flow limited.\n"
            )
            handle.write(
                "- In practice this means the remaining late cases are no longer all equivalent: some now look like final selection/tiebreak problems rather than missing boundary-gate evidence.\n"
            )
        handle.write(
            "- The current ranking oracle still shows only a narrow subset of `STRG.212` is fixable by a pure ranking heuristic, so the next patch should stay narrow and target the parity-backed late-emission cases directly.\n"
        )

    print(summary_path)
    print(per_gene_path)
    print(per_ref_path)
    print(md_path)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
