#!/usr/bin/env python3
"""Export a typed, block-witness provenance prototype for selected fresh O1 cases.

This is deliberately an evidence graph, not yet the final block-class algorithm.
Every RNA/DNA block node is one real passing minimap2 record. Overlapping pairwise
witnesses are retained rather than transitively collapsed into a fabricated block
class. All within-human duplication relations are UNROOTED.
"""

from __future__ import annotations

import argparse
from collections import defaultdict
from dataclasses import dataclass
from pathlib import Path

from o1_fresh_emission_validation import best_match, read_tsv
from o1_gene_family_audit import write_tsv


DEFAULT_CASES = ("GOLGA6_8", "MAGEA", "RABL2", "NBPF_core", "NBPF_repeat_bridge")

COLOURS = {
    "KEEP_KNOWN": "#34a853",
    "KEEP_KNOWN_OUTLIER": "#f9ab00",
    "KEEP_RNA_CANDIDATE": "#4285f4",
    "REVIEW_RELATED_OUTGROUP": "#a142f4",
    "REVIEW_CONFLICTING_GENE": "#ff6d01",
    "REVIEW_UNTYPED": "#9aa0a6",
    "RNA_BLOCK": "#f06292",
    "DNA_BLOCK": "#00acc1",
}


@dataclass(frozen=True)
class PafRecord:
    query: int
    query_len: int
    query_start: int
    query_end: int
    strand: str
    target: int
    target_len: int
    target_start: int
    target_end: int
    identity: float
    coverage: float


def parse_paf_line(line: str) -> PafRecord:
    fields = line.rstrip("\n").split("\t")
    if len(fields) < 12:
        raise ValueError("PAF record has fewer than 12 fields")
    tags = {}
    for field in fields[12:]:
        parts = field.split(":", 2)
        if len(parts) == 3:
            tags[parts[0]] = parts[2]
    query_len, target_len = int(fields[1]), int(fields[6])
    query_start, query_end = int(fields[2]), int(fields[3])
    target_start, target_end = int(fields[7]), int(fields[8])
    if "de" in tags:
        identity = 1.0 - float(tags["de"])
    else:
        identity = int(fields[9]) / max(int(fields[10]), 1)
    if query_len <= target_len:
        coverage = (query_end - query_start) / max(query_len, 1)
    else:
        coverage = (target_end - target_start) / max(target_len, 1)
    return PafRecord(
        int(fields[0]), query_len, query_start, query_end, fields[4],
        int(fields[5]), target_len, target_start, target_end, identity, coverage,
    )


def load_paf(path: Path) -> dict[tuple[int, int], list[PafRecord]]:
    records: dict[tuple[int, int], list[PafRecord]] = defaultdict(list)
    with path.open() as handle:
        for line in handle:
            if not line.strip():
                continue
            record = parse_paf_line(line)
            records[tuple(sorted((record.query, record.target)))].append(record)
    return records


def passing_witness(
    records: list[PafRecord],
    min_identity: float,
    min_coverage: float,
    require_forward: bool,
    expected_identity: float | None = None,
    expected_coverage: float | None = None,
) -> PafRecord | None:
    passing = [
        record for record in records
        if record.identity >= min_identity
        and record.coverage >= min_coverage
        and (not require_forward or record.strand == "+")
    ]
    if not passing:
        return None
    if expected_identity is not None and expected_coverage is not None:
        return min(
            passing,
            key=lambda record: (
                abs(record.identity - expected_identity)
                + abs(record.coverage - expected_coverage),
                -record.coverage,
                -record.identity,
            ),
        )
    return max(passing, key=lambda record: (record.coverage, record.identity))


def pair_id(a: tuple[str, str], b: tuple[str, str]) -> tuple[tuple[str, str], tuple[str, str]]:
    return (a, b) if a <= b else (b, a)


def edge_pair(a: str, b: str) -> tuple[str, str]:
    return (a, b) if a <= b else (b, a)


def component_sizes(nodes: set[str], edges: set[tuple[str, str]]) -> list[int]:
    """Return descending component sizes, including isolated nodes."""
    neighbours: dict[str, set[str]] = {node: set() for node in nodes}
    for left, right in edges:
        neighbours[left].add(right)
        neighbours[right].add(left)
    remaining = set(nodes)
    sizes: list[int] = []
    while remaining:
        seed = remaining.pop()
        stack = [seed]
        size = 0
        while stack:
            node = stack.pop()
            size += 1
            unseen = neighbours[node] & remaining
            remaining.difference_update(unseen)
            stack.extend(unseen)
        sizes.append(size)
    return sorted(sizes, reverse=True)


def membership_role(status: str) -> str:
    return {
        "KEEP_KNOWN": "FAMILY_MEMBER",
        "KEEP_KNOWN_OUTLIER": "FAMILY_MEMBER_OUTLIER",
        "KEEP_RNA_CANDIDATE": "CANDIDATE_MEMBER",
        "REVIEW_RELATED_OUTGROUP": "CONTEXT_OUTGROUP",
        "REVIEW_CONFLICTING_GENE": "EXCLUDE_CONFOUND",
        "REVIEW_UNTYPED": "REVIEW",
    }[status]


def disposition(rows: list[dict[str, str]]) -> tuple[str, str]:
    priority = (
        "REVIEW_CONFLICTING_GENE",
        "REVIEW_RELATED_OUTGROUP",
        "KEEP_KNOWN",
        "KEEP_KNOWN_OUTLIER",
        "KEEP_RNA_CANDIDATE",
        "REVIEW_UNTYPED",
    )
    values = {row["disposition"] for row in rows}
    status = next((value for value in priority if value in values), "REVIEW_UNTYPED")
    genes = sorted({gene for row in rows for gene in row["genes"].split(",") if gene != "."})
    return status, ",".join(genes) or "."


def interval_on_node(record: PafRecord, node_idx: int) -> tuple[int, int]:
    if record.query == node_idx:
        return record.query_start, record.query_end
    if record.target == node_idx:
        return record.target_start, record.target_end
    raise ValueError(f"node {node_idx} is not in PAF record")


def relation_row(
    case: str,
    source: str,
    target: str,
    relation: str,
    layer: str,
    witness: PafRecord | None,
    source_idx: int | None,
    target_idx: int | None,
    evidence: str,
) -> dict[str, object]:
    if witness is None:
        source_interval = target_interval = "."
        identity = coverage = "."
        orientation = "."
    else:
        assert source_idx is not None and target_idx is not None
        ss, se = interval_on_node(witness, source_idx)
        ts, te = interval_on_node(witness, target_idx)
        source_interval = f"{ss}-{se}"
        target_interval = f"{ts}-{te}"
        identity = f"{witness.identity:.4f}"
        coverage = f"{witness.coverage:.4f}"
        orientation = witness.strand
    return {
        "case": case,
        "source_id": source,
        "target_id": target,
        "relation": relation,
        "layer": layer,
        "identity": identity,
        "coverage": coverage,
        "orientation": orientation,
        "source_interval": source_interval,
        "target_interval": target_interval,
        "evidence": evidence,
        "direction_status": "UNROOTED",
    }


def write_case_graph(
    graph_dir: Path,
    case: str,
    entities: list[dict[str, object]],
    relations: list[dict[str, object]],
) -> None:
    graph = graph_dir / f"HSA.{case}.provenance_witness.gfa"
    colours = graph_dir / f"HSA.{case}.provenance_witness.colours.csv"
    with graph.open("w") as out, colours.open("w") as colour_out:
        out.write("H\tVN:Z:1.0\tTS:Z:typed_provenance_witness_prototype\n")
        colour_out.write("Node,Colour\n")
        for entity in entities:
            node = str(entity["entity_id"])
            kind = str(entity["entity_type"])
            status = str(entity["disposition"])
            length = max(int(entity["end"]) - int(entity["start"]), 1)
            out.write(
                f'S\t{node}\t*\tLN:i:{length}\tTY:Z:{kind}\tDP:Z:{status}\t'
                f'GN:Z:{entity["annotation"]}\tFF:Z:{entity["fresh_family"]}\n'
            )
            colour_out.write(f"{node},{COLOURS[status]}\n")
        for relation in relations:
            out.write(
                f'L\t{relation["source_id"]}\t+\t{relation["target_id"]}\t+\t0M\t'
                f'RT:Z:{relation["relation"]}\tLY:Z:{relation["layer"]}\t'
                f'DS:Z:{relation["direction_status"]}\n'
            )


def run(
    evidence_dir: Path,
    matches_path: Path,
    rna_paf_path: Path,
    joint_edges_path: Path,
    joint_dna_paf_path: Path,
    out_dir: Path,
    cases: tuple[str, ...],
) -> tuple[list[dict[str, object]], list[dict[str, object]], list[dict[str, object]]]:
    copies = read_tsv(evidence_dir / "HSA.fresh.copies.tsv")
    er_nodes = read_tsv(evidence_dir / "HSA.fresh.er.nodes.tsv")
    er_edges = read_tsv(evidence_dir / "HSA.fresh.er.edges.tsv")
    matches = [
        row for row in read_tsv(matches_path)
        if row["species"] == "HSA" and row["case"] in cases and row["fresh_emitted"] == "1"
    ]
    joint_edges = read_tsv(joint_edges_path)
    rna_paf = load_paf(rna_paf_path)
    dna_paf = load_paf(joint_dna_paf_path)

    node_key_by_copy: dict[tuple[str, str], str] = {}
    copy_by_key: dict[str, dict[str, str]] = {}
    for copy in copies:
        node, containment, _ = best_match(
            (copy["chrom"], int(copy["start"]), int(copy["end"])), er_nodes
        )
        if node is None or containment < 0.50:
            continue
        key = node["node_key"]
        node_key_by_copy[(copy["family_id"], copy["copy_idx"])] = key
        copy_by_key[key] = copy

    expected_by_key: dict[str, list[dict[str, str]]] = defaultdict(list)
    keys_by_case: dict[str, set[str]] = defaultdict(set)
    for row in matches:
        key = node_key_by_copy.get((row["fresh_family"], row["fresh_copy"]))
        if key is not None:
            expected_by_key[key].append(row)
            keys_by_case[row["case"]].add(key)

    idx_by_key = {row["node_key"]: int(row["idx"]) for row in er_nodes}
    # The joint PAF indexes the coordinate-sorted emitted copies in file order.
    global_joint_idx = {
        (copy["family_id"], copy["copy_idx"]): index for index, copy in enumerate(copies)
    }
    out_dir.mkdir(parents=True, exist_ok=True)
    graph_dir = out_dir / "graphs"
    graph_dir.mkdir(exist_ok=True)
    all_entities: list[dict[str, object]] = []
    all_relations: list[dict[str, object]] = []
    summaries: list[dict[str, object]] = []

    for case in cases:
        selected = keys_by_case[case]
        entities: list[dict[str, object]] = []
        relations: list[dict[str, object]] = []
        contained: dict[str, list[tuple[int, str, str]]] = defaultdict(list)
        for key in sorted(selected):
            copy = copy_by_key[key]
            status, annotation = disposition(expected_by_key[key])
            entities.append(
                {
                    "case": case,
                    "entity_id": key,
                    "entity_type": "LOCUS",
                    "chrom": copy["chrom"],
                    "start": copy["start"],
                    "end": copy["end"],
                    "strand": copy["strand"],
                    "fresh_family": copy["family_id"],
                    "fresh_copy": copy["copy_idx"],
                    "annotation": annotation,
                    "disposition": status,
                    "membership_role": membership_role(status),
                    "coordinate_space": "genome",
                }
            )

        rna_count = 0
        rna_block_count = 0
        rna_pairs: set[tuple[str, str]] = set()
        for edge in er_edges:
            a, b = edge["node_key_i"], edge["node_key_j"]
            if a not in selected or b not in selected:
                continue
            ai, bi = int(edge["rep_i"]), int(edge["rep_j"])
            witness = passing_witness(
                rna_paf.get(tuple(sorted((ai, bi))), []),
                0.60,
                0.50,
                True,
                float(edge["identity"]),
                float(edge["coverage"]),
            )
            if witness is None:
                raise RuntimeError(f"no RNA witness for {case}: {a} -- {b}")
            block_id = f"RB~{case}~{rna_block_count}"
            qspan = witness.query_end - witness.query_start
            tspan = witness.target_end - witness.target_start
            entities.append(
                {
                    "case": case,
                    "entity_id": block_id,
                    "entity_type": "BLOCK",
                    "chrom": ".",
                    "start": 0,
                    "end": min(qspan, tspan),
                    "strand": witness.strand,
                    "fresh_family": ".",
                    "fresh_copy": ".",
                    "annotation": ".",
                    "disposition": "RNA_BLOCK",
                    "membership_role": "EVIDENCE_WITNESS",
                    "coordinate_space": "exon_sum_pairwise_witness",
                }
            )
            relations.append(relation_row(case, a, b, "RNA_HOMOLOGY", "RNA", witness, ai, bi, "fresh_Er"))
            for node, idx in ((a, ai), (b, bi)):
                start, _ = interval_on_node(witness, idx)
                contained[node].append((start, block_id, "RNA"))
                relations.append(
                    relation_row(case, node, block_id, "CONTAINS", "RNA", witness, idx, bi if idx == ai else ai, "minimap2_PAF")
                )
            rna_count += 1
            rna_block_count += 1
            rna_pairs.add(edge_pair(a, b))

        selected_copies = {
            (copy_by_key[key]["family_id"], copy_by_key[key]["copy_idx"]): key
            for key in selected
        }
        dna_count = 0
        dna_block_count = 0
        dna_pairs: set[tuple[str, str]] = set()
        for edge in joint_edges:
            left = (edge["family_i"], edge["copy_i"])
            right = (edge["family_j"], edge["copy_j"])
            if left not in selected_copies or right not in selected_copies or edge["dna_edge"] != "true":
                continue
            a, b = selected_copies[left], selected_copies[right]
            ai, bi = global_joint_idx[left], global_joint_idx[right]
            witness = passing_witness(
                dna_paf.get(tuple(sorted((ai, bi))), []), 0.60, 0.50, False
            )
            if witness is None:
                raise RuntimeError(f"no DNA witness for {case}: {a} -- {b}")
            block_id = f"DB~{case}~{dna_block_count}"
            qspan = witness.query_end - witness.query_start
            tspan = witness.target_end - witness.target_start
            entities.append(
                {
                    "case": case,
                    "entity_id": block_id,
                    "entity_type": "BLOCK",
                    "chrom": ".",
                    "start": 0,
                    "end": min(qspan, tspan),
                    "strand": witness.strand,
                    "fresh_family": ".",
                    "fresh_copy": ".",
                    "annotation": ".",
                    "disposition": "DNA_BLOCK",
                    "membership_role": "EVIDENCE_WITNESS",
                    "coordinate_space": "transcript_normalised_genomic_span_pairwise_witness",
                }
            )
            relations.append(relation_row(case, a, b, "DNA_DUPLICATION", "DNA", witness, ai, bi, edge["evidence_class"]))
            for node, idx in ((a, ai), (b, bi)):
                start, _ = interval_on_node(witness, idx)
                contained[node].append((start, block_id, "DNA"))
                relations.append(
                    relation_row(case, node, block_id, "CONTAINS", "DNA", witness, idx, bi if idx == ai else ai, "minimap2_PAF")
                )
            dna_count += 1
            dna_block_count += 1
            dna_pairs.add(edge_pair(a, b))

        paths: list[dict[str, object]] = []
        for key in sorted(selected):
            ordered = sorted(contained[key])
            paths.append(
                {
                    "case": case,
                    "locus_id": key,
                    "path": ",".join(f"{block}:{layer}@{start}" for start, block, layer in ordered) or ".",
                    "path_status": (
                        "PAIRWISE_WITNESSES_NOT_YET_CONSOLIDATED"
                        if ordered else "NO_SELECTED_PAIRWISE_WITNESS"
                    ),
                }
            )
        write_tsv(out_dir / f"HSA.{case}.paths.tsv", paths, list(paths[0]))

        both_pairs = rna_pairs & dna_pairs
        entity_by_id = {str(entity["entity_id"]): entity for entity in entities}
        locus_rows: list[dict[str, object]] = []
        for key in sorted(selected):
            rna_neighbours = {
                other for pair in rna_pairs if key in pair for other in pair if other != key
            }
            dna_neighbours = {
                other for pair in dna_pairs if key in pair for other in pair if other != key
            }
            both_neighbours = {
                other for pair in both_pairs if key in pair for other in pair if other != key
            }
            entity = entity_by_id[key]
            locus_rows.append(
                {
                    "case": case,
                    "locus_id": key,
                    "annotation": entity["annotation"],
                    "disposition": entity["disposition"],
                    "membership_role": entity["membership_role"],
                    "fresh_family": entity["fresh_family"],
                    "rna_degree": len(rna_neighbours),
                    "dna_degree": len(dna_neighbours),
                    "both_layer_degree": len(both_neighbours),
                    "rna_only_degree": len(rna_neighbours - dna_neighbours),
                    "dna_only_degree": len(dna_neighbours - rna_neighbours),
                }
            )
        write_tsv(out_dir / f"HSA.{case}.loci.tsv", locus_rows, list(locus_rows[0]))

        write_case_graph(graph_dir, case, entities, relations)
        all_entities.extend(entities)
        all_relations.extend(relations)
        possible_pairs = len(selected) * (len(selected) - 1) // 2
        union_pairs = rna_pairs | dna_pairs
        rna_components = component_sizes(selected, rna_pairs)
        dna_components = component_sizes(selected, dna_pairs)
        family_by_key = {key: copy_by_key[key]["family_id"] for key in selected}
        summaries.append(
            {
                "case": case,
                "emitted_loci": len(selected),
                "fresh_family_ids": ",".join(sorted({copy_by_key[key]["family_id"] for key in selected})),
                "rna_homology_edges": rna_count,
                "rna_block_witnesses": rna_block_count,
                "dna_duplication_edges": dna_count,
                "dna_block_witnesses": dna_block_count,
                "both_layer_edges": len(both_pairs),
                "rna_only_edges": len(rna_pairs - dna_pairs),
                "dna_only_edges": len(dna_pairs - rna_pairs),
                "layer_edge_jaccard": f"{len(both_pairs) / len(union_pairs):.4f}" if union_pairs else ".",
                "rna_edge_density": f"{len(rna_pairs) / possible_pairs:.4f}" if possible_pairs else ".",
                "dna_edge_density": f"{len(dna_pairs) / possible_pairs:.4f}" if possible_pairs else ".",
                "rna_components": ",".join(map(str, rna_components)),
                "dna_components": ",".join(map(str, dna_components)),
                "rna_cross_fresh_family_edges": sum(
                    family_by_key[left] != family_by_key[right] for left, right in rna_pairs
                ),
                "dna_cross_fresh_family_edges": sum(
                    family_by_key[left] != family_by_key[right] for left, right in dna_pairs
                ),
                "rooted_provenance_edges": 0,
                "direction_status": "UNROOTED",
            }
        )

    write_tsv(out_dir / "entities.tsv", all_entities, list(all_entities[0]))
    write_tsv(out_dir / "relations.tsv", all_relations, list(all_relations[0]))
    write_tsv(out_dir / "case_summary.tsv", summaries, list(summaries[0]))
    return summaries, all_entities, all_relations


def write_report(out_dir: Path, summaries: list[dict[str, object]]) -> None:
    lines = [
        "# O1 duplication-provenance witness prototype",
        "",
        "This experiment instantiates the typed model on fresh HSA families already on disk.",
        "Every block entity is one real passing minimap2 PAF witness. RNA and DNA relations remain",
        "separate; all directions are `UNROOTED`. Pairwise witnesses deliberately remain overlapping",
        "and unmerged, so this is not yet the final non-transitive block-class decomposition.",
        "",
        "| case | loci | fresh families | RNA/DNA/both | layer Jaccard | RNA comps | DNA comps | direction |",
        "|---|---:|---|---:|---:|---|---|---|",
    ]
    for row in summaries:
        lines.append(
            f'| {row["case"]} | {row["emitted_loci"]} | {row["fresh_family_ids"]} | '
            f'{row["rna_homology_edges"]}/{row["dna_duplication_edges"]}/{row["both_layer_edges"]} | '
            f'{row["layer_edge_jaccard"]} | {row["rna_components"]} | {row["dna_components"]} | '
            f'{row["direction_status"]} |'
        )
    lines += [
        "",
        "## Interpretation",
        "",
        "A connected component in one layer is not promoted to a family in another layer. In",
        "particular, a DNA block witness records shared duplicated sequence but does not create an",
        "RNA family edge. The GFA files are visual projections; `entities.tsv`, `relations.tsv`, and",
        "the per-case locus/path tables are the normative evidence. Component-size lists include",
        "isolated loci, so `15,1,1` means one 15-locus component and two singletons.",
        "",
        "The next implementation step is to consolidate overlapping pairwise witnesses into stable",
        "multi-locus block classes using reciprocal overlap plus quasi-clique cohesion. Until then,",
        "the path tables explicitly say `NOT_YET_CONSOLIDATED`, and no ancestry is inferred.",
        "",
    ]
    (out_dir / "README.md").write_text("\n".join(lines))


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--evidence-dir", type=Path,
        default=Path("bench/o1_fresh_emission_validation"),
    )
    parser.add_argument(
        "--matches", type=Path,
        default=Path("bench/o1_fresh_emission_validation/fresh_node_matches.tsv"),
    )
    parser.add_argument(
        "--rna-paf", type=Path,
        default=Path("/tmp/rustle_o1_fresh_emission/HSA.fresh_er.er._k11_w5.0.paf"),
    )
    parser.add_argument(
        "--joint-edges", type=Path,
        default=Path("/tmp/rustle_o1_fresh_emission/HSA.joint.joint_edges.tsv"),
    )
    parser.add_argument(
        "--joint-dna-paf", type=Path,
        default=Path("/tmp/rustle_o1_fresh_emission/HSA.joint_er.joint_dna._k11_w5.2.paf"),
    )
    parser.add_argument(
        "--out-dir", type=Path,
        default=Path("bench/o1_provenance_witness_prototype"),
    )
    parser.add_argument("--cases", nargs="+", default=list(DEFAULT_CASES))
    args = parser.parse_args()
    summaries, _, _ = run(
        args.evidence_dir,
        args.matches,
        args.rna_paf,
        args.joint_edges,
        args.joint_dna_paf,
        args.out_dir,
        tuple(args.cases),
    )
    write_report(args.out_dir, summaries)
    print(f"wrote provenance-witness prototype to {args.out_dir}")


if __name__ == "__main__":
    main()
