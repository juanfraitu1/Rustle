#!/usr/bin/env python3

from __future__ import annotations

import argparse
import csv
import re
from collections import defaultdict
from dataclasses import dataclass
from pathlib import Path


NOTE_INT_RE = re.compile(r"([A-Za-z_]+)=([0-9]+)")


@dataclass
class GeneLocus:
    gene_id: str
    chrom: str
    start: int
    end: int
    strand: str
    transcript_count: int = 0


def parse_attrs(attr_field: str) -> dict[str, str]:
    out: dict[str, str] = {}
    for part in attr_field.split(";"):
        part = part.strip()
        if not part or " " not in part:
            continue
        key, value = part.split(" ", 1)
        out[key] = value.strip().strip('"')
    return out


def load_gene_loci(gtf_path: Path, gene_ids: list[str]) -> dict[str, GeneLocus]:
    wanted = set(gene_ids)
    loci: dict[str, GeneLocus] = {}
    with gtf_path.open() as handle:
        for line in handle:
            if not line or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 9 or fields[2] != "transcript":
                continue
            attrs = parse_attrs(fields[8])
            gene_id = attrs.get("gene_id")
            if gene_id not in wanted:
                continue
            chrom = fields[0]
            start = int(fields[3]) - 1
            end = int(fields[4])
            strand = fields[6]
            if gene_id not in loci:
                loci[gene_id] = GeneLocus(
                    gene_id=gene_id,
                    chrom=chrom,
                    start=start,
                    end=end,
                    strand=strand,
                    transcript_count=1,
                )
            else:
                locus = loci[gene_id]
                locus.start = min(locus.start, start)
                locus.end = max(locus.end, end)
                locus.transcript_count += 1
    return loci


def overlaps(chrom: str, start: int, end: int, locus: GeneLocus) -> bool:
    return chrom == locus.chrom and not (end <= locus.start or start >= locus.end)


def parse_note(note: str) -> dict[str, int]:
    return {k: int(v) for k, v in NOTE_INT_RE.findall(note)}


def aggregate_live_stage_rows(parity_tsv: Path, loci: dict[str, GeneLocus]) -> dict[str, dict[str, dict[str, int]]]:
    out: dict[str, dict[str, dict[str, int]]] = {
        gene_id: defaultdict(lambda: defaultdict(int)) for gene_id in loci
    }
    with parity_tsv.open() as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            chrom = row["chrom"]
            start = int(row["start"])
            end = int(row["end"])
            for gene_id, locus in loci.items():
                if not overlaps(chrom, start, end, locus):
                    continue
                stage = row["stage_id"]
                stage_out = out[gene_id][stage]
                stage_out["rows"] += 1
                for key in (
                    "n_reads",
                    "n_junctions",
                    "n_nodes",
                    "n_edges",
                    "n_transfrags",
                    "n_transcripts",
                    "seed_total",
                    "seed_stored",
                    "seed_rescued",
                    "bundle_bnodes",
                    "bundle_read_bnode_links",
                    "bundle_color_count",
                    "struct_hash",
                    "graph_active_transfrags",
                    "graph_long_transfrags",
                    "graph_pattern_size",
                    "graph_hard_boundary_nodes",
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
                ):
                    value = row.get(key)
                    if value:
                        if key == "struct_hash":
                            stage_out[key] ^= int(value)
                        else:
                            stage_out[key] += int(value)
                for key, value in parse_note(row["note"]).items():
                    stage_out.setdefault(key, 0)
                    stage_out[key] += value
    return out


def parse_stage_int(row: dict[str, str], key: str) -> int:
    value = row.get(key, "")
    return int(value) if value else 0


def load_backward_stage_counts(backward_dir: Path, gene_ids: list[str]) -> dict[str, dict[str, dict[str, int]]]:
    out: dict[str, dict[str, dict[str, int]]] = {}
    for gene_id in gene_ids:
        candidates = [
            backward_dir / f"isoseq_backward_trace_{gene_id}_cachecold.backward_trace_stages.tsv",
            backward_dir / f"isoseq_backward_trace_{gene_id.replace('.', '')}_cachecold.backward_trace_stages.tsv",
            backward_dir / f"backward_{gene_id}.backward_trace_stages.tsv",
            backward_dir / f"backward_{gene_id.replace('.', '')}.backward_trace_stages.tsv",
        ]
        stages: dict[str, dict[str, int]] = {}
        stage_path = next((path for path in candidates if path.exists()), None)
        if stage_path is not None:
            with stage_path.open() as handle:
                reader = csv.DictReader(handle, delimiter="\t")
                for row in reader:
                    stages[row["stage_id"]] = {
                        "locus_trace_rows": parse_stage_int(row, "locus_trace_rows"),
                        "bundle_create_rows": parse_stage_int(row, "bundle_create_rows"),
                        "bundle_process_rows": parse_stage_int(row, "bundle_process_rows"),
                        "bundle_struct_hash": parse_stage_int(row, "bundle_struct_hash"),
                        "graph_distinct_nodes": parse_stage_int(row, "graph_distinct_nodes"),
                        "graph_distinct_bids": parse_stage_int(row, "graph_distinct_bids"),
                        "graph_distinct_ops": parse_stage_int(row, "graph_distinct_ops"),
                        "graph_struct_hash": parse_stage_int(row, "graph_struct_hash"),
                    }
        out[gene_id] = stages
    return out


def earliest_divergence(live: dict[str, dict[str, int]]) -> tuple[str, str]:
    seed = live.get("seed_flow", {})
    seed_fail = (
        seed.get("seed_unwitnessed", seed.get("unwitnessed", 0))
        + seed.get("seed_zero_flux", seed.get("zero_flux", 0))
        + seed.get("seed_low_cov", seed.get("low_cov", 0))
        + seed.get("seed_back_fail", seed.get("back_fail", 0))
        + seed.get("seed_fwd_fail", seed.get("fwd_fail", 0))
    )
    if seed_fail > 0:
        return "seed_flow", (
            f"unwitnessed={seed.get('seed_unwitnessed', seed.get('unwitnessed', 0))} "
            f"zero_flux={seed.get('seed_zero_flux', seed.get('zero_flux', 0))} "
            f"low_cov={seed.get('seed_low_cov', seed.get('low_cov', 0))} "
            f"back_fail={seed.get('seed_back_fail', seed.get('back_fail', 0))} "
            f"fwd_fail={seed.get('seed_fwd_fail', seed.get('fwd_fail', 0))}"
        )
    checktrf = live.get("checktrf", {})
    unresolved = checktrf.get("checktrf_total", 0) - checktrf.get("seed_rescued", checktrf.get("rescued", 0))
    if unresolved > 0:
        return "checktrf", (
            f"checktrf_total={checktrf.get('checktrf_total', 0)} rescued={checktrf.get('seed_rescued', checktrf.get('rescued', 0))}"
        )
    for stage_id in (
        "pairwise_overlap_filter",
        "isofrac",
        "collapse_single_exon_runoff",
        "polymerase_runoff_filter",
        "readthr_gate",
        "junction_support_filter",
    ):
        stage = live.get(stage_id, {})
        removed = stage.get("predcluster_removed", stage.get("removed", 0))
        if removed > 0:
            return stage_id, (
                f"before={stage.get('predcluster_before', stage.get('before', 0))} "
                f"removed={removed} survivors={stage.get('n_transcripts', 0)}"
            )
    final_stage = live.get("final_family_selection", {})
    final_nonstored = max(0, final_stage.get("seed_total", 0) - final_stage.get("seed_stored", 0))
    if final_nonstored > 0:
        return "final_family_selection", (
            f"seed_total={final_stage.get('seed_total', 0)} seed_stored={final_stage.get('seed_stored', 0)} "
            f"seed_rescued={final_stage.get('seed_rescued', 0)}"
        )
    return "none", "no live divergence signal in compact stages"


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--reference-gtf", type=Path, default=Path("GGO_19.gtf"))
    ap.add_argument("--live-parity-tsv", type=Path, default=Path("/tmp/rustle_compat_fullrun2.parity.tsv"))
    ap.add_argument("--backward-dir", type=Path, default=Path("/tmp"))
    ap.add_argument(
        "--genes",
        nargs="+",
        default=["STRG.15", "STRG.88", "STRG.171", "STRG.212"],
    )
    ap.add_argument("--out-tsv", type=Path, default=Path("/tmp/live_vs_backward_calibration.tsv"))
    ap.add_argument("--out-md", type=Path, default=Path("/tmp/live_vs_backward_calibration.md"))
    args = ap.parse_args()

    loci = load_gene_loci(args.reference_gtf, args.genes)
    live = aggregate_live_stage_rows(args.live_parity_tsv, loci)
    backward = load_backward_stage_counts(args.backward_dir, args.genes)

    with args.out_tsv.open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(
            [
                "gene_id",
                "chrom",
                "start",
                "end",
                "strand",
                "ref_transcripts",
                "live_bundle_rows",
                "live_bundle_bnodes",
                "live_bundle_read_bnode_links",
                "live_bundle_color_count",
                "live_bundle_hash",
                "live_graph_rows",
                "live_graph_active_transfrags",
                "live_graph_long_transfrags",
                "live_graph_pattern_size",
                "live_graph_hard_boundary_nodes",
                "live_graph_hash",
                "backward_bundle_rows",
                "backward_bundle_create_rows",
                "backward_bundle_process_rows",
                "backward_bundle_hash",
                "backward_graph_rows",
                "backward_graph_distinct_nodes",
                "backward_graph_distinct_bids",
                "backward_graph_distinct_ops",
                "backward_graph_hash",
                "live_seed_rows",
                "live_seed_fail_sum",
                "live_checktrf_rows",
                "live_checktrf_total",
                "live_checktrf_rescued",
                "live_final_rows",
                "live_final_seed_total",
                "live_final_seed_stored",
                "live_final_seed_rescued",
                "backward_seed_rows",
                "backward_checktrf_rows",
                "backward_final_rows",
                "earliest_live_divergence",
                "divergence_detail",
            ]
        )
        for gene_id in args.genes:
            locus = loci[gene_id]
            live_gene = live[gene_id]
            bundle_stage = live_gene.get("bundle_context", {})
            graph_stage = live_gene.get("graph_evolution", {})
            backward_bundle = backward[gene_id].get("bundle_context", {})
            backward_graph = backward[gene_id].get("graph_evolution", {})
            seed = live_gene.get("seed_flow", {})
            checktrf = live_gene.get("checktrf", {})
            final = live_gene.get("final_family_selection", {})
            seed_fail_sum = (
                seed.get("seed_unwitnessed", seed.get("unwitnessed", 0))
                + seed.get("seed_zero_flux", seed.get("zero_flux", 0))
                + seed.get("seed_low_cov", seed.get("low_cov", 0))
                + seed.get("seed_back_fail", seed.get("back_fail", 0))
                + seed.get("seed_fwd_fail", seed.get("fwd_fail", 0))
            )
            earliest, detail = earliest_divergence(live_gene)
            writer.writerow(
                [
                    gene_id,
                    locus.chrom,
                    locus.start,
                    locus.end,
                    locus.strand,
                    locus.transcript_count,
                    bundle_stage.get("rows", 0),
                    bundle_stage.get("bundle_bnodes", 0),
                    bundle_stage.get("bundle_read_bnode_links", 0),
                    bundle_stage.get("bundle_color_count", 0),
                    bundle_stage.get("struct_hash", 0),
                    graph_stage.get("rows", 0),
                    graph_stage.get("graph_active_transfrags", 0),
                    graph_stage.get("graph_long_transfrags", 0),
                    graph_stage.get("graph_pattern_size", 0),
                    graph_stage.get("graph_hard_boundary_nodes", 0),
                    graph_stage.get("struct_hash", 0),
                    backward_bundle.get("locus_trace_rows", 0),
                    backward_bundle.get("bundle_create_rows", 0),
                    backward_bundle.get("bundle_process_rows", 0),
                    backward_bundle.get("bundle_struct_hash", 0),
                    backward_graph.get("locus_trace_rows", 0),
                    backward_graph.get("graph_distinct_nodes", 0),
                    backward_graph.get("graph_distinct_bids", 0),
                    backward_graph.get("graph_distinct_ops", 0),
                    backward_graph.get("graph_struct_hash", 0),
                    seed.get("rows", 0),
                    seed_fail_sum,
                    checktrf.get("rows", 0),
                    checktrf.get("checktrf_total", 0),
                    checktrf.get("seed_rescued", checktrf.get("rescued", 0)),
                    final.get("rows", 0),
                    final.get("seed_total", 0),
                    final.get("seed_stored", 0),
                    final.get("seed_rescued", 0),
                    backward[gene_id].get("seed_flow", {}).get("locus_trace_rows", 0),
                    backward[gene_id].get("checktrf", {}).get("locus_trace_rows", 0),
                    backward[gene_id].get("final_family_selection", {}).get("locus_trace_rows", 0),
                    earliest,
                    detail,
                ]
            )

    with args.out_md.open("w") as handle:
        handle.write("# Live vs Backward Calibration\n\n")
        handle.write(
            f"Live parity source: `{args.live_parity_tsv}`\n\n"
            f"Backward stage source dir: `{args.backward_dir}`\n\n"
        )
        for gene_id in args.genes:
            locus = loci[gene_id]
            live_gene = live[gene_id]
            bundle_stage = live_gene.get("bundle_context", {})
            graph_stage = live_gene.get("graph_evolution", {})
            backward_bundle = backward[gene_id].get("bundle_context", {})
            backward_graph = backward[gene_id].get("graph_evolution", {})
            seed = live_gene.get("seed_flow", {})
            checktrf = live_gene.get("checktrf", {})
            final = live_gene.get("final_family_selection", {})
            earliest, detail = earliest_divergence(live_gene)
            handle.write(f"## {gene_id}\n\n")
            handle.write(
                f"- locus: `{locus.chrom}:{locus.start}-{locus.end}({locus.strand})`, ref transcripts `{locus.transcript_count}`\n"
            )
            handle.write(
                f"- live rows: bundle `{bundle_stage.get('rows', 0)}`, "
                f"graph `{graph_stage.get('rows', 0)}`, "
                f"seed `{seed.get('rows', 0)}`, checktrf `{checktrf.get('rows', 0)}`, "
                f"final `{final.get('rows', 0)}`\n"
            )
            handle.write(
                f"- live bundle detail: bnodes `{bundle_stage.get('bundle_bnodes', 0)}`, "
                f"read_bnode_links `{bundle_stage.get('bundle_read_bnode_links', 0)}`, "
                f"colors `{bundle_stage.get('bundle_color_count', 0)}`, "
                f"hash `{bundle_stage.get('struct_hash', 0)}`\n"
            )
            handle.write(
                f"- live graph detail: active_trf `{graph_stage.get('graph_active_transfrags', 0)}`, "
                f"long_trf `{graph_stage.get('graph_long_transfrags', 0)}`, "
                f"pattern_size `{graph_stage.get('graph_pattern_size', 0)}`, "
                f"hard_boundary_nodes `{graph_stage.get('graph_hard_boundary_nodes', 0)}`, "
                f"hash `{graph_stage.get('struct_hash', 0)}`\n"
            )
            handle.write(
                f"- backward bundle detail: rows `{backward_bundle.get('locus_trace_rows', 0)}`, "
                f"create_rows `{backward_bundle.get('bundle_create_rows', 0)}`, "
                f"process_rows `{backward_bundle.get('bundle_process_rows', 0)}`, "
                f"hash `{backward_bundle.get('bundle_struct_hash', 0)}`\n"
            )
            handle.write(
                f"- backward graph detail: rows `{backward_graph.get('locus_trace_rows', 0)}`, "
                f"distinct_nodes `{backward_graph.get('graph_distinct_nodes', 0)}`, "
                f"distinct_bids `{backward_graph.get('graph_distinct_bids', 0)}`, "
                f"distinct_ops `{backward_graph.get('graph_distinct_ops', 0)}`, "
                f"hash `{backward_graph.get('graph_struct_hash', 0)}`\n"
            )
            handle.write(
                f"- live seed failures: unwitnessed `{seed.get('seed_unwitnessed', seed.get('unwitnessed', 0))}`, "
                f"zero_flux `{seed.get('seed_zero_flux', seed.get('zero_flux', 0))}`, "
                f"low_cov `{seed.get('seed_low_cov', seed.get('low_cov', 0))}`, "
                f"back_fail `{seed.get('seed_back_fail', seed.get('back_fail', 0))}`, "
                f"fwd_fail `{seed.get('seed_fwd_fail', seed.get('fwd_fail', 0))}`\n"
            )
            handle.write(
                f"- live checktrf: total `{checktrf.get('checktrf_total', 0)}`, rescued `{checktrf.get('seed_rescued', checktrf.get('rescued', 0))}`\n"
            )
            handle.write(
                f"- live final: seed_total `{final.get('seed_total', 0)}`, seed_stored `{final.get('seed_stored', 0)}`, seed_rescued `{final.get('seed_rescued', 0)}`\n"
            )
            handle.write(
                f"- backward rows: seed `{backward[gene_id].get('seed_flow', {}).get('locus_trace_rows', 0)}`, "
                f"checktrf `{backward[gene_id].get('checktrf', {}).get('locus_trace_rows', 0)}`, "
                f"final `{backward[gene_id].get('final_family_selection', {}).get('locus_trace_rows', 0)}`\n"
            )
            handle.write(f"- earliest live divergence: `{earliest}`\n")
            handle.write(f"- detail: `{detail}`\n\n")


if __name__ == "__main__":
    main()
