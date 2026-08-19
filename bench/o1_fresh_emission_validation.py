#!/usr/bin/env python3
"""Freshly emit the expanded O1 panel and compare loci only after prediction.

The expanded purity audit reads frozen catalog nodes.  This script takes one step
upstream: it builds a predeclared BED from those loci, extracts records from the
original full-genome-aligned BAMs, runs the current ``gw_family_catalog`` binary,
and only then compares fresh copy coordinates with the frozen audit nodes.

The experiment proves reproducible node construction on the regional substrate.  It
does not prove whole-genome partition identity because the panel BAM contains only
records overlapping the predeclared loci.
"""

from __future__ import annotations

import argparse
import csv
import os
import shutil
import subprocess
from collections import Counter, defaultdict
from pathlib import Path

from o1_expanded_family_audit import CASES
from o1_gene_family_audit import write_tsv


INPUTS = {
    "GGO": {
        "bam": Path("/home/juanfra/winloci_scratch/GGO_mm.bam"),
        "fasta": Path("/home/juanfra/winloci_scratch/GGO.fasta"),
    },
    "HSA": {
        "bam": Path("/mnt/linuxdisk/home/juanfraitu/winloci_data/A119b.t2t.bam"),
        "fasta": Path("/mnt/linuxdisk/home/juanfraitu/winloci_data/chm13v2.0.fa"),
    },
}

FRESH_SNAPSHOTS = (
    ("fresh.copies.tsv", "copies.tsv"),
    ("fresh.families.tsv", "families.tsv"),
    ("fresh.log", "run.log"),
    ("fresh_er.nodes.tsv", "er.nodes.tsv"),
    ("fresh_er.edges.tsv", "er.edges.tsv"),
    ("fresh_er.rule.tsv", "er.rule.tsv"),
    ("fresh_er.params.tsv", "er.params.tsv"),
)

FRESH_COLOURS = {
    "KNOWN_TARGET": "#34a853",
    "RNA_CANDIDATE": "#4285f4",
    "RELATED_OUTGROUP": "#a142f4",
    "CONFLICT_IN_TARGET": "#d93025",
    "CONFLICT_SEPARATED": "#ff6d01",
    "REVIEW_UNTYPED": "#9aa0a6",
    "UNSCORED": "#dadce0",
}


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open() as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def parse_locus(value: str) -> tuple[str, int, int]:
    chrom, coords = value.rsplit(":", 1)
    start, end = map(int, coords.split("-"))
    return chrom, start, end


def merge_intervals(intervals: list[tuple[str, int, int]]) -> list[tuple[str, int, int]]:
    answer: list[tuple[str, int, int]] = []
    for chrom, start, end in sorted(intervals):
        if answer and answer[-1][0] == chrom and start <= answer[-1][2]:
            old_chrom, old_start, old_end = answer[-1]
            answer[-1] = (old_chrom, old_start, max(old_end, end))
        else:
            answer.append((chrom, start, end))
    return answer


def write_panel_beds(nodes: list[dict[str, str]], work: Path, padding: int) -> dict[str, Path]:
    beds = {}
    for species in INPUTS:
        intervals = []
        for row in nodes:
            if row["species"] != species:
                continue
            chrom, start, end = parse_locus(row["locus"])
            intervals.append((chrom, max(0, start - padding), end + padding))
        path = work / f"{species}.panel.bed"
        with path.open("w") as fh:
            for chrom, start, end in merge_intervals(intervals):
                fh.write(f"{chrom}\t{start}\t{end}\n")
        beds[species] = path
    return beds


def run_fresh(binary: Path, nodes: list[dict[str, str]], work: Path, padding: int, threads: int) -> None:
    work.mkdir(parents=True, exist_ok=True)
    beds = write_panel_beds(nodes, work, padding)
    for species, paths in INPUTS.items():
        bam = work / f"{species}.panel.bam"
        prefix = work / f"{species}.fresh"
        log = work / f"{species}.fresh.log"
        subprocess.run(
            ["samtools", "view", "-b", "-M", "-L", str(beds[species]),
             str(paths["bam"]), "-o", str(bam)],
            check=True,
        )
        subprocess.run(["samtools", "index", "-@", str(threads), str(bam)], check=True)
        env = os.environ.copy()
        env["RUSTLE_ER_EDGE_DUMP"] = str(work / f"{species}.fresh_er")
        with log.open("w") as stderr:
            subprocess.run(
                [str(binary), "--bam", str(bam), "--fasta", str(paths["fasta"]),
                 "--threads", str(threads), "--rna-forward-only", "--out", str(prefix)],
                check=True,
                stdout=stderr,
                stderr=subprocess.STDOUT,
                env=env,
            )


def best_match(
    expected: tuple[str, int, int], fresh: list[dict[str, str]]
) -> tuple[dict[str, str] | None, float, float]:
    chrom, start, end = expected
    best: tuple[float, float, float, dict[str, str]] | None = None
    for row in fresh:
        if row["chrom"] != chrom:
            continue
        other_start, other_end = int(row["start"]), int(row["end"])
        overlap = max(0, min(end, other_end) - max(start, other_start))
        if overlap == 0:
            continue
        containment = overlap / max(min(end - start, other_end - other_start), 1)
        jaccard = overlap / max(max(end, other_end) - min(start, other_start), 1)
        centre_delta = abs((start + end) - (other_start + other_end))
        candidate = (containment, jaccard, -centre_delta, row)
        if best is None or candidate[:3] > best[:3]:
            best = candidate
    if best is None:
        return None, 0.0, 0.0
    return best[3], best[0], best[1]


def compare(nodes: list[dict[str, str]], work: Path, out_dir: Path) -> tuple[list[dict[str, object]], list[dict[str, object]]]:
    fresh_by_species = {
        species: read_tsv(work / f"{species}.fresh.copies.tsv")
        for species in INPUTS
    }
    matches: list[dict[str, object]] = []
    matched_by_case: dict[tuple[str, str], list[dict[str, object]]] = defaultdict(list)
    for row in nodes:
        expected = parse_locus(row["locus"])
        fresh, containment, jaccard = best_match(expected, fresh_by_species[row["species"]])
        emitted = fresh is not None and containment >= 0.50
        result: dict[str, object] = {
            **row,
            "fresh_emitted": int(emitted),
            "overlap_of_shorter": f"{containment:.4f}",
            "jaccard": f"{jaccard:.4f}",
            "fresh_family": fresh["family_id"] if emitted and fresh else ".",
            "fresh_copy": fresh["copy_idx"] if emitted and fresh else ".",
            "fresh_locus": (
                f'{fresh["chrom"]}:{fresh["start"]}-{fresh["end"]}'
                if emitted and fresh else "."
            ),
        }
        matches.append(result)
        matched_by_case[(row["species"], row["case"])].append(result)

    summaries: list[dict[str, object]] = []
    for species, label, old_family, _ in CASES:
        subset = matched_by_case[(species, label)]
        named_matches = [
            row for row in subset
            if row["disposition"].startswith("KEEP_KNOWN") and row["fresh_emitted"] == 1
        ]
        votes = Counter(str(row["fresh_family"]) for row in named_matches)
        target_family = sorted(votes, key=lambda value: (-votes[value], value))[0] if votes else "."
        for row in subset:
            row["fresh_target_family"] = target_family
            row["in_fresh_target_family"] = int(
                row["fresh_emitted"] == 1 and row["fresh_family"] == target_family
            )
        conflicts = [row for row in subset if row["disposition"] == "REVIEW_CONFLICTING_GENE"]
        related = [row for row in subset if row["disposition"] == "REVIEW_RELATED_OUTGROUP"]
        primary = [row for row in subset if row["disposition"].startswith("KEEP_")]
        named = [row for row in subset if row["disposition"].startswith("KEEP_KNOWN")]
        summaries.append(
            {
                "species": species,
                "case": label,
                "frozen_family": old_family,
                "fresh_target_family": target_family,
                "frozen_nodes": len(subset),
                "fresh_matched_nodes": sum(int(row["fresh_emitted"]) for row in subset),
                "fresh_unique_nodes": len({
                    (row["fresh_family"], row["fresh_copy"])
                    for row in subset if row["fresh_emitted"] == 1
                }),
                "named_expected": len(named),
                "named_emitted": sum(int(row["fresh_emitted"]) for row in named),
                "named_in_target_family": sum(int(row["in_fresh_target_family"]) for row in named),
                "primary_expected": len(primary),
                "primary_in_target_family": sum(int(row["in_fresh_target_family"]) for row in primary),
                "conflicts_expected": len(conflicts),
                "conflicts_reemitted": sum(int(row["fresh_emitted"]) for row in conflicts),
                "conflicts_in_target_family": sum(int(row["in_fresh_target_family"]) for row in conflicts),
                "related_outgroups_expected": len(related),
                "related_outgroups_reemitted": sum(int(row["fresh_emitted"]) for row in related),
                "related_outgroups_in_target_family": sum(int(row["in_fresh_target_family"]) for row in related),
                "fresh_family_ids_for_named": ",".join(sorted({
                    str(row["fresh_family"]) for row in named if row["fresh_emitted"] == 1
                })) or ".",
            }
        )

    out_dir.mkdir(parents=True, exist_ok=True)
    write_tsv(out_dir / "fresh_node_matches.tsv", matches, list(matches[0]))
    write_tsv(out_dir / "fresh_family_summary.tsv", summaries, list(summaries[0]))
    return summaries, matches


def snapshot_fresh_outputs(work: Path, out_dir: Path) -> None:
    """Keep the small, auditable Rustle outputs; omit BAM, FASTA, and minimap scratch."""
    for species in INPUTS:
        for source_suffix, destination_suffix in FRESH_SNAPSHOTS:
            source = work / f"{species}.{source_suffix}"
            if not source.exists():
                raise RuntimeError(f"missing fresh output {source}")
            shutil.copyfile(source, out_dir / f"{species}.fresh.{destination_suffix}")


def fresh_node_status(rows: list[dict[str, object]]) -> tuple[str, str, str]:
    """Collapse post-hoc expected-node matches to one visual label per fresh node."""
    genes = sorted({str(row["genes"]) for row in rows if row["genes"] != "."})
    cases = sorted({str(row["case"]) for row in rows})
    if any(
        row["disposition"] == "REVIEW_CONFLICTING_GENE"
        and row["in_fresh_target_family"] == 1
        for row in rows
    ):
        status = "CONFLICT_IN_TARGET"
    elif any(row["disposition"] == "REVIEW_CONFLICTING_GENE" for row in rows):
        status = "CONFLICT_SEPARATED"
    elif any(row["disposition"] == "REVIEW_RELATED_OUTGROUP" for row in rows):
        status = "RELATED_OUTGROUP"
    elif any(str(row["disposition"]).startswith("KEEP_KNOWN") for row in rows):
        status = "KNOWN_TARGET"
    elif any(row["disposition"] == "KEEP_RNA_CANDIDATE" for row in rows):
        status = "RNA_CANDIDATE"
    else:
        status = "REVIEW_UNTYPED"
    return status, ",".join(genes) or ".", ",".join(cases) or "."


def write_fresh_graphs(
    work: Path, out_dir: Path, matches: list[dict[str, object]]
) -> list[dict[str, object]]:
    """Emit GFA views of actual fresh family representatives and actual E_r edges."""
    graph_dir = out_dir / "graphs"
    graph_dir.mkdir(exist_ok=True)
    stats: list[dict[str, object]] = []
    for species in INPUTS:
        copies = read_tsv(work / f"{species}.fresh.copies.tsv")
        er_nodes = read_tsv(work / f"{species}.fresh_er.nodes.tsv")
        er_edges = read_tsv(work / f"{species}.fresh_er.edges.tsv")

        copy_by_key: dict[str, dict[str, str]] = {}
        copy_key_by_id: dict[tuple[str, str], str] = {}
        for copy in copies:
            er_node, containment, _ = best_match(
                (copy["chrom"], int(copy["start"]), int(copy["end"])), er_nodes
            )
            if er_node is None or containment < 0.50:
                continue
            key = er_node["node_key"]
            copy_by_key[key] = copy
            copy_key_by_id[(copy["family_id"], copy["copy_idx"])] = key

        evidence: dict[str, list[dict[str, object]]] = defaultdict(list)
        for row in matches:
            if row["species"] != species or row["fresh_emitted"] != 1:
                continue
            key = copy_key_by_id.get((str(row["fresh_family"]), str(row["fresh_copy"])))
            if key is not None:
                evidence[key].append(row)

        graph_path = graph_dir / f"{species}.fresh_emitted.gfa"
        colour_path = graph_dir / f"{species}.fresh_emitted.colours.csv"
        with graph_path.open("w") as graph, colour_path.open("w") as colours:
            graph.write("H\tVN:Z:1.0\tTS:Z:fresh_emitted_Er_family_graph\n")
            colours.write("Node,Colour\n")
            for key, copy in sorted(copy_by_key.items()):
                if key in evidence:
                    status, genes, cases = fresh_node_status(evidence[key])
                else:
                    status, genes, cases = "UNSCORED", ".", "."
                span = int(copy["end"]) - int(copy["start"])
                graph.write(
                    f'S\t{key}\t*\tLN:i:{span}\tFF:Z:{copy["family_id"]}\t'
                    f'CI:i:{copy["copy_idx"]}\tDP:Z:{status}\tGN:Z:{genes}\tCA:Z:{cases}\t'
                    f'LO:Z:{copy["chrom"]}:{copy["start"]}-{copy["end"]}\n'
                )
                colours.write(f"{key},{FRESH_COLOURS[status]}\n")
            graph_edges = 0
            for edge in er_edges:
                a, b = edge["node_key_i"], edge["node_key_j"]
                if a not in copy_by_key or b not in copy_by_key:
                    continue
                graph.write(
                    f'L\t{a}\t+\t{b}\t+\t0M\tID:f:{edge["identity"]}\t'
                    f'CV:f:{edge["coverage"]}\tMT:Z:{edge["metric_tier"]}\n'
                )
                graph_edges += 1
        stats.append(
            {
                "species": species,
                "fresh_family_copies": len(copies),
                "copies_mapped_to_er_nodes": len(copy_by_key),
                "fresh_er_edges_in_graph": graph_edges,
            }
        )
    write_tsv(out_dir / "fresh_graph_summary.tsv", stats, list(stats[0]))
    return stats


def o2_dependency_warnings(work: Path) -> list[str]:
    warnings = set()
    for species in INPUTS:
        log = work / f"{species}.fresh.log"
        for line in log.read_text().splitlines():
            if "[o1-perp-o2] WARNING:" in line:
                warnings.add(f"{species}: {line.split('WARNING:', 1)[1].strip()}")
    return sorted(warnings)


def write_report(
    out_dir: Path,
    summaries: list[dict[str, object]],
    matches: list[dict[str, object]],
    graph_stats: list[dict[str, object]],
    o2_warnings: list[str],
    padding: int,
) -> None:
    total = sum(int(row["frozen_nodes"]) for row in summaries)
    matched = sum(int(row["fresh_matched_nodes"]) for row in summaries)
    named = sum(int(row["named_expected"]) for row in summaries)
    named_emitted = sum(int(row["named_emitted"]) for row in summaries)
    named_in_target = sum(int(row["named_in_target_family"]) for row in summaries)
    conflicts = sum(int(row["conflicts_expected"]) for row in summaries)
    conflicts_emitted = sum(int(row["conflicts_reemitted"]) for row in summaries)
    conflicts_same = sum(int(row["conflicts_in_target_family"]) for row in summaries)
    related = sum(int(row["related_outgroups_expected"]) for row in summaries)
    related_emitted = sum(int(row["related_outgroups_reemitted"]) for row in summaries)
    related_same = sum(int(row["related_outgroups_in_target_family"]) for row in summaries)
    lines = [
        "# O1 fresh-emission validation",
        "",
        "The current Rustle binary rebuilt nodes from regional BAMs extracted from the original",
        "full-genome alignments. Frozen node ids, family ids, and dispositions were not inputs to Rustle;",
        "they were used only after emission for coordinate matching.",
        "",
        f"Intervals were fixed from the 19-family audit with ±{padding:,} bp padding. Fresh loci match",
        f"**{matched}/{total}** frozen nodes at ≥0.50 overlap of the shorter interval, including",
        f"**{named_emitted}/{named}** independently named target nodes; **{named_in_target}/{named}** land in",
        f"the modal fresh family for their test case. Of **{conflicts}** previously",
        f"flagged conflicting-gene nodes, **{conflicts_emitted}** are independently re-emitted and",
        f"**{conflicts_same}** re-enter the same fresh component as the named target family.",
        f"Separately, **{related_emitted}/{related}** documented broad-family outgroups are re-emitted,",
        f"and **{related_same}/{related}** remain in the broad RNA family.",
        "",
        "| species | family | frozen→fresh matched | named emitted/in target | related emitted/in target | conflicts emitted/in target | fresh named-family ids |",
        "|---|---|---:|---:|---:|---:|---|",
    ]
    for row in summaries:
        lines.append(
            f'| {row["species"]} | {row["case"]} | {row["frozen_nodes"]}→{row["fresh_matched_nodes"]} '
            f'({row["fresh_unique_nodes"]} unique) | {row["named_emitted"]}/{row["named_in_target_family"]} | '
            f'{row["related_outgroups_reemitted"]}/{row["related_outgroups_in_target_family"]} | '
            f'{row["conflicts_reemitted"]}/{row["conflicts_in_target_family"]} | '
            f'{row["fresh_family_ids_for_named"]} |'
        )
    lines += [
        "",
        "## Emitted evidence",
        "",
        "The `*.fresh.copies.tsv`, `*.fresh.families.tsv`, and `*.fresh.er.*.tsv` files are direct",
        "snapshots from the new Rustle runs. The two `graphs/*.fresh_emitted.gfa` files contain only",
        "freshly emitted family copies and the actual fresh `E_r` edges between representative nodes;",
        "no audit edge is synthesized. Colours are post-hoc labels: green named target, blue RNA",
        "candidate, purple documented broad-family outgroup, red conflicting gene still in the target",
        "family, orange conflicting gene emitted in another family, grey review/unscored.",
        "",
        "| species | fresh family copies | mapped to E_r nodes | fresh E_r edges in GFA |",
        "|---|---:|---:|---:|",
    ]
    for row in graph_stats:
        lines.append(
            f'| {row["species"]} | {row["fresh_family_copies"]} | '
            f'{row["copies_mapped_to_er_nodes"]} | {row["fresh_er_edges_in_graph"]} |'
        )
    if o2_warnings:
        lines += [
            "",
            "Run-certificate disclosure: " + " ".join(o2_warnings),
            "Consequently, this fresh HSA node set is not a function of sequence evidence alone.",
        ]
    lines += [
        "",
        "## Interpretation boundary",
        "",
        "A re-emitted conflicting node is not an audit artifact: current upstream locus construction",
        "produced it again. If it also enters the target's fresh component, the false merge is an O1",
        "family-edge/typing problem. Failure to re-emit is evidence of representative-selection or",
        "regional-substrate instability, not proof that the old node was biologically false.",
        "A related outgroup is different: its broad-family edge is biologically expected, while its",
        "exclusion applies only to the recent-copy subfamily view.",
        "",
        "This remains a regional reconstruction. `samtools -M -L` preserves records from the original",
        "whole-genome alignment but removes records outside the panel intervals. Therefore the experiment",
        "validates node reproducibility and local regrouping, not byte-identical whole-genome partitioning.",
        "",
    ]
    (out_dir / "README.md").write_text("\n".join(lines))


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--audit-nodes", type=Path,
        default=Path("bench/o1_expanded_family_audit/expanded_family_nodes.tsv"),
    )
    parser.add_argument(
        "--binary", type=Path,
        default=Path("/tmp/rustle_o1_target/release/gw_family_catalog"),
    )
    parser.add_argument("--work", type=Path, default=Path("/tmp/rustle_o1_fresh_emission"))
    parser.add_argument(
        "--out-dir", type=Path,
        default=Path("bench/o1_fresh_emission_validation"),
    )
    parser.add_argument("--padding", type=int, default=10_000)
    parser.add_argument("--threads", type=int, default=4)
    parser.add_argument(
        "--reuse", action="store_true",
        help="reuse existing fresh Rustle outputs and only rerun comparison/reporting",
    )
    args = parser.parse_args()
    nodes = read_tsv(args.audit_nodes)
    if not args.reuse:
        run_fresh(args.binary, nodes, args.work, args.padding, args.threads)
    summaries, matches = compare(nodes, args.work, args.out_dir)
    snapshot_fresh_outputs(args.work, args.out_dir)
    graph_stats = write_fresh_graphs(args.work, args.out_dir, matches)
    warnings = o2_dependency_warnings(args.work)
    write_report(args.out_dir, summaries, matches, graph_stats, warnings, args.padding)
    print(f"wrote fresh-emission validation to {args.out_dir}")


if __name__ == "__main__":
    main()
