#!/usr/bin/env python3
"""Expanded post-hoc O1 purity audit over frozen whole-catalog RNA graphs.

This complements ``o1_gene_family_audit.py``.  The original seven-family panel has
purpose-built DNA/RNA evidence and can test member recovery.  This larger panel uses
already-emitted GGO/HSA catalogs and therefore tests only a narrower question:
given a recognizable family-bearing graph, which nodes can be typed as members,
which remain plausible unnamed candidates, and which are independently named as a
different gene?  It is not a discovery-recall denominator.
"""

from __future__ import annotations

import argparse
import csv
from collections import defaultdict, deque
from pathlib import Path

from o1_gene_family_audit import components, paf_edges, pair, write_tsv


# Deliberately includes clean, difficult, and known repeat-bridge cases.  Catalog
# ids merely locate the frozen graph; gene symbols provide the independent typing.
CASES = (
    ("GGO", "C4", "GWFAM179", ("C4A", "C4B")),
    ("GGO", "TSPYL", "GWFAM196", ("TSPYL",)),
    ("GGO", "GSTM", "GWFAM34", ("GSTM",)),
    ("GGO", "RGPD", "GWFAM332", ("RGPD",)),
    ("GGO", "ANKRD18", "GWFAM338", ("ANKRD18",)),
    ("GGO", "HERC2", "GWFAM387", ("HERC2",)),
    ("GGO", "MAGEA", "GWFAM455", ("MAGEA",)),
    ("GGO", "DAZ", "GWFAM479", ("DAZ",)),
    ("GGO", "APOBEC3", "GWFAM491", ("APOBEC3",)),
    ("HSA", "GOLGA6_8", "GWFAM141", ("GOLGA6", "GOLGA8")),
    ("HSA", "TBC1D3", "GWFAM235", ("TBC1D3",)),
    ("HSA", "RANBP2_RGPD", "GWFAM280", ("RANBP2", "RGPD")),
    ("HSA", "RABL2", "GWFAM288", ("RABL2",)),
    ("HSA", "AMY1", "GWFAM32", ("AMY1",)),
    ("HSA", "PCDHB", "GWFAM341", ("PCDHB",)),
    ("HSA", "MAGEA", "GWFAM375", ("MAGEA",)),
    ("HSA", "RBMY1", "GWFAM393", ("RBMY1",)),
    ("HSA", "NBPF_core", "GWFAM5", ("NBPF",)),
    ("HSA", "NBPF_repeat_bridge", "GWFAM33", ("NBPF",)),
)

# A target subfamily can have an independently documented ancestral homolog that is
# a valid member of the broad homology family but not of the recent-copy subfamily.
# Keep this narrow and evidence-backed: it is an evaluation truth-layer, never an
# input to Rustle's discovery graph.
RELATED_OUTGROUPS = {
    ("HSA", "GOLGA6_8"): {"GOLGA2"},
}

COLOURS = {
    "KEEP_KNOWN": "#34a853",
    "KEEP_KNOWN_OUTLIER": "#f9ab00",
    "KEEP_RNA_CANDIDATE": "#4285f4",
    "REVIEW_RELATED_OUTGROUP": "#a142f4",
    "REVIEW_CONFLICTING_GENE": "#ff6d01",
    "REVIEW_UNTYPED": "#9aa0a6",
}


def genes(row: dict[str, str]) -> set[str]:
    return {g for g in row["all_pc_genes"].split(",") if g and g != "."}


def best_gene(row: dict[str, str]) -> str:
    value = row.get("best_pc_gene", ".")
    return value if value and value != "." else "."


def matches_family(gene: str, stems: tuple[str, ...]) -> bool:
    return any(gene.upper().startswith(stem) for stem in stems)


def informative_other_genes(row: dict[str, str], stems: tuple[str, ...]) -> set[str]:
    return {
        gene
        for gene in genes(row)
        if not gene.upper().startswith("LOC") and not matches_family(gene, stems)
    }


def reachable(seeds: set[str], edges: set[tuple[str, str]]) -> set[str]:
    adjacent: dict[str, set[str]] = defaultdict(set)
    for a, b in edges:
        adjacent[a].add(b)
        adjacent[b].add(a)
    seen = set(seeds)
    todo = deque(seeds)
    while todo:
        node = todo.popleft()
        for other in adjacent[node] - seen:
            seen.add(other)
            todo.append(other)
    return seen


def component_sizes(nodes: set[str], edges: set[tuple[str, str]]) -> str:
    return "+".join(str(len(c)) for c in components(nodes, edges))


def load_catalog(path: Path) -> dict[str, list[dict[str, str]]]:
    answer: dict[str, list[dict[str, str]]] = defaultdict(list)
    with path.open() as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            answer[row["family_id"]].append(row)
    return answer


def write_graph(
    path: Path,
    rows: dict[str, dict[str, str]],
    nodes: set[str],
    edges: set[tuple[str, str]],
    dispositions: dict[str, tuple[str, str]],
    all_forward: set[tuple[str, str]],
    audit: bool,
) -> None:
    with path.open("w") as fh:
        graph_type = "expanded_gene_family_audit" if audit else "expanded_gene_family_primary"
        fh.write(f"H\tVN:Z:1.0\tTS:Z:{graph_type}\n")
        for node in sorted(nodes):
            row = rows[node]
            status, _ = dispositions[node]
            symbol = ",".join(sorted(genes(row))) or "."
            fh.write(
                f'S\t{node}\t*\tLN:i:{row["span"]}\tDP:Z:{status}\tGN:Z:{symbol}\t'
                f'LO:Z:{row["chrom"]}:{row["start"]}-{row["end"]}\n'
            )
        for a, b in sorted(edges):
            orientation = "FORWARD" if (a, b) in all_forward else "REVERSE_ONLY"
            member = int(dispositions[a][0].startswith("KEEP_") and dispositions[b][0].startswith("KEEP_"))
            fh.write(f"L\t{a}\t+\t{b}\t+\t0M\tOR:Z:{orientation}\tMB:i:{member}\n")


def run(catalog_root: Path, out_dir: Path) -> tuple[list[dict[str, object]], list[dict[str, object]]]:
    out_dir.mkdir(parents=True, exist_ok=True)
    graph_dir = out_dir / "graphs"
    graph_dir.mkdir(exist_ok=True)
    catalogs = {
        species: load_catalog(catalog_root / f"{species}_copy_annot.tsv")
        for species in {case[0] for case in CASES}
    }
    edge_sets = {
        species: paf_edges(catalog_root / f"{species}_copies.er.paf")
        for species in catalogs
    }
    family_rows: list[dict[str, object]] = []
    node_rows: list[dict[str, object]] = []

    for species, label, family_id, stems in CASES:
        source_rows = catalogs[species].get(family_id)
        if not source_rows:
            raise RuntimeError(f"missing frozen family {species}/{family_id}")
        by_node = {f'{family_id}~{row["copy_idx"]}': row for row in source_rows}
        nodes = set(by_node)
        all_edges, forward_edges = edge_sets[species]
        all_local = {edge for edge in all_edges if set(edge) <= nodes}
        forward_local = {edge for edge in forward_edges if set(edge) <= nodes}
        known = {
            node for node, row in by_node.items()
            if any(matches_family(gene, stems) for gene in genes(row))
        }
        related_names = RELATED_OUTGROUPS.get((species, label), set())
        related = {
            node for node, row in by_node.items()
            if node not in known and best_gene(row) in related_names
        }
        conflicting = {
            node for node, row in by_node.items()
            if node not in known and node not in related and informative_other_genes(row, stems)
        }
        non_subfamily = conflicting | related
        clean_forward = {
            edge for edge in forward_local if not (set(edge) & non_subfamily)
        }
        reached = reachable(known, clean_forward)
        dispositions: dict[str, tuple[str, str]] = {}
        for node, row in by_node.items():
            other = informative_other_genes(row, stems)
            if node in known:
                supported = bool(reachable(known - {node}, clean_forward) & {node})
                status = "KEEP_KNOWN" if supported else "KEEP_KNOWN_OUTLIER"
                reason = "independent same-family protein-coding annotation"
            elif node in related:
                status = "REVIEW_RELATED_OUTGROUP"
                reason = (
                    f"best annotation {best_gene(row)} is a documented ancestral/broad-family "
                    "homolog, but not a recent target-subfamily copy"
                )
            elif other:
                status = "REVIEW_CONFLICTING_GENE"
                reason = "independent annotation names another gene: " + ",".join(sorted(other))
            elif node in reached:
                status = "KEEP_RNA_CANDIDATE"
                reason = "forward RNA homology reaches a named family member; no conflicting named gene"
            else:
                status = "REVIEW_UNTYPED"
                reason = "no same-family annotation or forward path to a named member"
            dispositions[node] = (status, reason)

        primary = {node for node, (status, _) in dispositions.items() if status.startswith("KEEP_")}
        primary_edges = {edge for edge in clean_forward if set(edge) <= primary}
        known_components = components(known, {edge for edge in clean_forward if set(edge) <= known})
        if len(known) < 2:
            status = "ANNOTATION_LIMITED"
        elif len(known_components) == 1:
            status = "PASS_NAMED_CONNECTED"
        else:
            status = "PARTIAL_NAMED_SPLIT"
        possible = len(primary) * (len(primary) - 1) // 2
        family_rows.append(
            {
                "species": species,
                "case": label,
                "catalog_family": family_id,
                "audit_nodes": len(nodes),
                "named_members": len(known),
                "primary_nodes": len(primary),
                "rna_candidates": sum(s == "KEEP_RNA_CANDIDATE" for s, _ in dispositions.values()),
                "related_outgroups": sum(s == "REVIEW_RELATED_OUTGROUP" for s, _ in dispositions.values()),
                "conflicting_nodes": sum(s == "REVIEW_CONFLICTING_GENE" for s, _ in dispositions.values()),
                "untyped_review": sum(s == "REVIEW_UNTYPED" for s, _ in dispositions.values()),
                "edges_r0": len(all_local),
                "edges_forward": len(forward_local),
                "reverse_only_removed": len(all_local - forward_local),
                "components_r0": component_sizes(nodes, all_local),
                "components_forward": component_sizes(nodes, forward_local),
                "primary_components": component_sizes(primary, primary_edges),
                "primary_density": f"{len(primary_edges) / max(possible, 1):.4f}",
                "status": status,
            }
        )
        for node, row in sorted(by_node.items()):
            disposition, reason = dispositions[node]
            node_rows.append(
                {
                    "species": species,
                    "case": label,
                    "catalog_family": family_id,
                    "node": node,
                    "locus": f'{row["chrom"]}:{row["start"]}-{row["end"]}',
                    "genes": ",".join(sorted(genes(row))) or ".",
                    "disposition": disposition,
                    "reason": reason,
                }
            )

        prefix = f"{species}.{label}"
        write_graph(
            graph_dir / f"{prefix}.gene_family.gfa", by_node, primary, primary_edges,
            dispositions, forward_local, audit=False,
        )
        write_graph(
            graph_dir / f"{prefix}.audit.gfa", by_node, nodes, all_local,
            dispositions, forward_local, audit=True,
        )
        for suffix, graph_nodes in (("", primary), (".audit", nodes)):
            with (graph_dir / f"{prefix}{suffix}.colours.csv").open("w") as fh:
                fh.write("Node,Colour\n")
                for node in sorted(graph_nodes):
                    fh.write(f"{node},{COLOURS[dispositions[node][0]]}\n")

    write_tsv(out_dir / "expanded_family_certificates.tsv", family_rows, list(family_rows[0]))
    write_tsv(out_dir / "expanded_family_nodes.tsv", node_rows, list(node_rows[0]))
    return family_rows, node_rows


def write_report(out_dir: Path, families: list[dict[str, object]], nodes: list[dict[str, object]]) -> None:
    total_audit = sum(int(row["audit_nodes"]) for row in families)
    total_primary = sum(int(row["primary_nodes"]) for row in families)
    total_named = sum(int(row["named_members"]) for row in families)
    total_conflict = sum(int(row["conflicting_nodes"]) for row in families)
    total_related = sum(int(row["related_outgroups"]) for row in families)
    total_reverse = sum(int(row["reverse_only_removed"]) for row in families)
    passed = sum(row["status"] == "PASS_NAMED_CONNECTED" for row in families)
    testable = sum(int(row["named_members"]) >= 2 for row in families)
    annotation_limited = sum(row["status"] == "ANNOTATION_LIMITED" for row in families)
    lines = [
        "# Expanded O1 known-family purity audit",
        "",
        "**Scope:** 19 post-hoc family-bearing graphs from frozen gorilla and human catalogs. This",
        "measures within-emitted-graph typing purity and orientation sensitivity. It is not discovery",
        "recall: choosing a catalog family after emission conditions on Rustle having found it.",
        "",
        f"Primary graphs retain **{total_primary}/{total_audit}** audit nodes and all **{total_named}**",
        f"independently named target-family members. They withhold **{total_conflict}** nodes named as",
        f"unrelated genes and report **{total_related}** broad-family/recent-subfamily outgroup separately.",
        f"The orientation guard removes **{total_reverse}** within-family edges. Named",
        f"members remain connected in **{passed}/{testable} testable** cases; the remaining",
        f"**{annotation_limited}** cases have only one independently named member and are reported as",
        "annotation-limited rather than forced passes.",
        "",
        "| species | family | audit→primary nodes | named | related outgroups | conflicts | R0→forward edges | primary components | status |",
        "|---|---|---:|---:|---:|---:|---:|---:|---|",
    ]
    for row in families:
        lines.append(
            f'| {row["species"]} | {row["case"]} | {row["audit_nodes"]}→{row["primary_nodes"]} | '
            f'{row["named_members"]} | {row["related_outgroups"]} | {row["conflicting_nodes"]} | '
            f'{row["edges_r0"]}→{row["edges_forward"]} | {row["primary_components"]} | {row["status"]} |'
        )
    lines += [
        "",
        "`graphs/<species>.<family>.gene_family.gfa` is the purified view. Its sibling `.audit.gfa`",
        "retains every emitted node and reverse-only edge. Purple nodes are documented broad-family",
        "homologs outside the recent-copy subfamily; orange nodes are independently named unrelated",
        "genes; grey nodes are untyped; blue nodes are unnamed candidates connected to a named member.",
        "",
        "The panel deliberately contains difficult mixtures (GOLGA6/8, RANBP2/RGPD, MAGEA, NBPF) and",
        "the known NBPF/TTC6-DNAH14 repeat-bridge component. These are more informative for purity than",
        "adding only clean two-copy pairs. Family stems are an evaluation typing device, not a production",
        "coordinate blacklist or a substitute for O1 discovery.",
        "",
    ]
    (out_dir / "README.md").write_text("\n".join(lines))


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--catalog-root", type=Path,
        default=Path("/home/juanfra/winloci_scratch/o1_errorcensus"),
    )
    parser.add_argument(
        "--out-dir", type=Path,
        default=Path("bench/o1_expanded_family_audit"),
    )
    args = parser.parse_args()
    families, nodes = run(args.catalog_root, args.out_dir)
    write_report(args.out_dir, families, nodes)
    print(f"wrote expanded O1 audit to {args.out_dir}")


if __name__ == "__main__":
    main()
