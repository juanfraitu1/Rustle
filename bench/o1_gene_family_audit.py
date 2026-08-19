#!/usr/bin/env python3
"""Audit O1 as a *gene-family* graph, not merely an SD-similarity graph.

The seven-family panel was originally assembled from segmental-duplication (SD)
intervals.  An SD edge is useful discovery evidence, but it is not sufficient to
call two loci members of the same gene family.  This audit therefore keeps three
objects separate:

* transcript-oriented RNA E_r edges (the O1 membership evidence);
* genomic DNA E_r edges (structural corroboration / RNA-null rescue evidence);
* independent gene annotations (the type check from SD family to gene family).

It writes a compact GFA for each known family, a node-level disposition table,
and a separately typed audit of the 18 previously reviewed Soto-missing loci.
No result from this script mutates a Rustle catalog or blacklists a coordinate.
"""

from __future__ import annotations

import argparse
import csv
import re
import urllib.parse
from collections import defaultdict, deque
from pathlib import Path


FAMILIES = ("NPIP", "TBC1D3", "RABL2", "APOBEC3", "MAGEA", "GSTM", "HERC2")
GAMMA = 0.20

# These are family-symbol tests, not claims that every overlapping LOC record is
# a member.  Unnamed nodes must still connect to a named anchor through forward
# transcript homology before they enter the primary RNA graph.
PANEL_STEMS = {
    "NPIP": ("NPIP",),
    "TBC1D3": ("TBC1D3",),
    "RABL2": ("RABL2",),
    "APOBEC3": ("APOBEC3",),
    "MAGEA": ("MAGEA",),
    "GSTM": ("GSTM",),
    "HERC2": ("HERC2",),
}

# Adjudicated panel contaminants.  The rule is deliberately conjunctive:
# independent other-gene annotation + no family seed + no qualifying forward
# RNA path to a named family anchor.  Coordinates are panel certificates, not a
# production blacklist.
OBVIOUS_PANEL_FALSE_POSITIVES = {
    "L~chr16_2117540_2119871": "PKD1 interval; RNA-null and no NPIP seed",
    "L~chr16_15105563_15115512": (
        "PDXDC1/LOC100505915 interval; NPIP connection is reverse-complement-only"
    ),
    "L~chr17_48453760_48455808": "NPEPPS interval; RNA-null and no TBC1D3 seed",
}

DISPOSITION_COLOURS = {
    "KEEP_KNOWN": "#34a853",
    "KEEP_KNOWN_OUTLIER": "#f9ab00",
    "KEEP_RNA_CANDIDATE": "#4285f4",
    "REVIEW_CONFLICTING_GENE": "#ff6d01",
    "REVIEW_SD_ONLY": "#a142f4",
    "REVIEW_DNA_ONLY": "#9aa0a6",
    "REMOVE_OBVIOUS_FP": "#ea4335",
}

# Soto IDs can contain one coherent gene family, several overlapping gene
# families, or anonymous SD sequence.  Empty tuples intentionally mean that the
# SD cluster cannot type a gene family without additional adjudication.
CANDIDATE_FAMILY_STEMS = {
    "ID_393": ("MST1",),
    "ID_400": ("NBPF", "NOTCH2NL"),
    "ID_212": (),
    "ID_116": ("GOLGA6",),
    "ID_71": ("HERC2",),
    "ID_179": ("ULK4",),
    "ID_163": ("GUSB",),
    "ID_35": (),
    "ID_280": ("ANKRD20A", "ANKYRIN REPEAT DOMAIN-CONTAINING PROTEIN 20A"),
}


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open() as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def write_tsv(path: Path, rows: list[dict[str, object]], fields: list[str]) -> None:
    with path.open("w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def pair(a: str, b: str) -> tuple[str, str]:
    return tuple(sorted((a, b)))


def components(nodes: set[str], edges: set[tuple[str, str]]) -> list[set[str]]:
    adjacent = {node: set() for node in nodes}
    for a, b in edges:
        if a in adjacent and b in adjacent:
            adjacent[a].add(b)
            adjacent[b].add(a)
    unseen = set(nodes)
    answer = []
    while unseen:
        seed = min(unseen)
        todo = [seed]
        component = set()
        while todo:
            node = todo.pop()
            if node not in unseen:
                continue
            unseen.remove(node)
            component.add(node)
            todo.extend(adjacent[node] & unseen)
        answer.append(component)
    return sorted(answer, key=lambda c: (-len(c), min(c)))


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


def tsv_edges(path: Path) -> set[tuple[str, str]]:
    return {pair(row["node_a"], row["node_b"]) for row in read_tsv(path)}


def paf_edges(path: Path) -> tuple[set[tuple[str, str]], set[tuple[str, str]]]:
    """Return (all R0 edges, forward R0 edges) using Rustle's shipped floor."""
    all_edges: set[tuple[str, str]] = set()
    forward_edges: set[tuple[str, str]] = set()
    with path.open() as fh:
        for line in fh:
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 12 or fields[0] == fields[5]:
                continue
            qlen, qstart, qend = map(int, (fields[1], fields[2], fields[3]))
            tlen, tstart, tend = map(int, (fields[6], fields[7], fields[8]))
            de = next((float(tag[5:]) for tag in fields[12:] if tag.startswith("de:f:")), None)
            identity = 1.0 - de if de is not None else int(fields[9]) / max(int(fields[10]), 1)
            # This is the current production single-record predicate.  The span
            # is measured on whichever PAF axis is the shorter full sequence.
            span = qend - qstart if qlen <= tlen else tend - tstart
            coverage = span / max(min(qlen, tlen), 1)
            if identity < 0.60 or coverage < 0.50:
                continue
            edge = pair(fields[0], fields[5])
            all_edges.add(edge)
            if fields[4] == "+":
                forward_edges.add(edge)
    return all_edges, forward_edges


def symbols(row: dict[str, str]) -> set[str]:
    out = set()
    for column in ("seed_genes", "overlapping_refseq_genes"):
        out.update(x for x in row[column].split(";") if x and x != ".")
    return out


def has_family_symbol(row: dict[str, str], family: str) -> bool:
    return any(any(symbol.upper().startswith(stem) for stem in PANEL_STEMS[family]) for symbol in symbols(row))


def conflicting_gene_symbols(row: dict[str, str], family: str) -> set[str]:
    """Return independently named, non-family genes overlapping an untyped candidate.

    LOC-style placeholders do not establish a conflicting gene identity.  A named
    gene such as PDXDC1 does: forward homology can keep it reviewable, but cannot by
    itself relabel that annotated gene as an NPIP copy.  This is a typing guard, not
    a sequence threshold and not a coordinate blacklist.
    """
    family_stems = PANEL_STEMS[family]
    return {
        symbol
        for symbol in symbols(row)
        if not symbol.upper().startswith("LOC")
        and not any(symbol.upper().startswith(stem) for stem in family_stems)
    }


def format_components(nodes: set[str], edges: set[tuple[str, str]]) -> str:
    return "+".join(str(len(component)) for component in components(nodes, edges))


def panel_audit(panel_root: Path, orientation_root: Path, out_dir: Path) -> list[dict[str, object]]:
    node_rows_out: list[dict[str, object]] = []
    family_rows: list[dict[str, object]] = []
    graph_dir = out_dir / "known_family_graphs"
    graph_dir.mkdir(parents=True, exist_ok=True)

    for family in FAMILIES:
        rows = read_tsv(panel_root / "nodes" / f"{family}.tsv")
        by_node = {row["node"]: row for row in rows}
        nodes = set(by_node)
        named = {node for node, row in by_node.items() if has_family_symbol(row, family)}
        rna_eligible = {node for node, row in by_node.items() if int(row["n_exon_blocks"]) > 0}
        r0_edges, rna_edges = paf_edges(orientation_root / f"{family}.transcript_oriented.paf")
        r0_edges = {edge for edge in r0_edges if set(edge) <= nodes}
        rna_edges = {edge for edge in rna_edges if set(edge) <= nodes}
        dna_edges = tsv_edges(panel_root / "edges" / f"{family}.dna.tsv")
        dna_edges = {edge for edge in dna_edges if set(edge) <= nodes}

        # RNA may recruit an unnamed copy, while DNA may rescue an independently
        # named RNA-null member.  A node annotated as another named gene is removed
        # before reachability so it cannot act as a hidden bridge that recruits more
        # anonymous nodes after that conflicting node itself is withheld.
        conflicting_nodes = {
            node for node, row in by_node.items()
            if node not in named and conflicting_gene_symbols(row, family)
        }
        raw_rna_reached = reachable(named, rna_edges)
        clean_rna_edges = {
            edge for edge in rna_edges if not (set(edge) & conflicting_nodes)
        }
        rna_reached = reachable(named, clean_rna_edges)
        typed_dna_edges = {
            edge for edge in dna_edges if edge[0] in named and edge[1] in named
        }
        membership_edges = clean_rna_edges | typed_dna_edges
        retained = reachable(named, membership_edges)

        dispositions: dict[str, tuple[str, str]] = {}
        for node in sorted(nodes):
            if node in OBVIOUS_PANEL_FALSE_POSITIVES and node not in rna_reached:
                dispositions[node] = ("REMOVE_OBVIOUS_FP", OBVIOUS_PANEL_FALSE_POSITIVES[node])
            elif node in named:
                if node in reachable(named - {node}, membership_edges):
                    dispositions[node] = ("KEEP_KNOWN", "independent family symbol and graph support")
                else:
                    dispositions[node] = (
                        "KEEP_KNOWN_OUTLIER",
                        "independent family symbol but no passing family edge; report as a false negative",
                    )
            elif node in conflicting_nodes and node in raw_rna_reached:
                conflicts = conflicting_gene_symbols(by_node[node], family)
                dispositions[node] = (
                    "REVIEW_CONFLICTING_GENE",
                    "forward transcript homology reaches the family, but independent annotation "
                    f'names a different gene: {",".join(sorted(conflicts))}',
                )
            elif node in rna_reached:
                dispositions[node] = (
                    "KEEP_RNA_CANDIDATE",
                    "forward transcript homology reaches an independently named family member",
                )
            elif node in rna_eligible:
                dispositions[node] = (
                    "REVIEW_SD_ONLY",
                    "RNA evidence exists but no forward transcript path reaches a named member",
                )
            else:
                dispositions[node] = (
                    "REVIEW_DNA_ONLY",
                    "DNA SD homology without RNA or independent gene-family typing",
                )

        primary = {node for node, (status, _) in dispositions.items() if status.startswith("KEEP_")}
        primary_edges = {edge for edge in membership_edges if set(edge) <= primary}
        known_supported = sum(any(node in edge for edge in membership_edges) for node in named)
        known_outliers = len(named) - known_supported
        joint = dna_edges | rna_edges
        both = dna_edges & rna_edges
        primary_component_sizes = format_components(primary, primary_edges)
        primary_connected = len(components(primary, primary_edges)) == 1
        possible_primary_edges = len(primary) * (len(primary) - 1) // 2
        primary_density = len(primary_edges) / max(possible_primary_edges, 1)
        if len(named) < 2:
            family_status = "ANNOTATION_LIMITED"
        elif known_outliers == 0 and primary_connected and primary_density >= GAMMA:
            family_status = "PASS"
        else:
            family_status = "PARTIAL_KNOWN_OUTLIER"

        family_rows.append(
            {
                "family": family,
                "panel_nodes": len(nodes),
                "named_gene_members": len(named),
                "primary_nodes": len(primary),
                "primary_edges": len(primary_edges),
                "primary_components": primary_component_sizes,
                "primary_density": f"{primary_density:.4f}",
                "gamma": f"{GAMMA:.2f}",
                "rna_eligible": len(rna_eligible),
                "rna_forward_edges": len(rna_edges),
                "rna_components": format_components(nodes, rna_edges),
                "dna_edges": len(dna_edges),
                "dna_components": format_components(nodes, dna_edges),
                "joint_edges": len(joint),
                "joint_components": format_components(nodes, joint),
                "edges_both": len(both),
                "edges_rna_only": len(rna_edges - dna_edges),
                "edges_dna_only": len(dna_edges - rna_edges),
                "known_supported": known_supported,
                "known_outliers": known_outliers,
                "obvious_fp_removed": sum(s == "REMOVE_OBVIOUS_FP" for s, _ in dispositions.values()),
                "conflicting_gene_review": sum(
                    s == "REVIEW_CONFLICTING_GENE" for s, _ in dispositions.values()
                ),
                "status": family_status,
            }
        )

        for node in sorted(nodes):
            status, reason = dispositions[node]
            row = by_node[node]
            node_rows_out.append(
                {
                    "family": family,
                    "node": node,
                    "locus": f'{row["chrom"]}:{row["start0"]}-{row["end"]}',
                    "symbols": ";".join(sorted(symbols(row))) or ".",
                    "rna_eligible": int(node in rna_eligible),
                    "named_gene_member": int(node in named),
                    "rna_reaches_named": int(node in rna_reached),
                    "dna_degree": sum(node in edge for edge in dna_edges),
                    "rna_forward_degree": sum(node in edge for edge in rna_edges),
                    "disposition": status,
                    "reason": reason,
                }
            )

        # The default gene-family graph is deliberately pure: REVIEW/REMOVE nodes
        # are not rendered as family members.  The companion audit graph retains
        # the complete SD discovery universe and marks which links were admitted.
        with (graph_dir / f"{family}.gene_family.gfa").open("w") as fh:
            fh.write("H\tVN:Z:1.0\tTS:Z:typed_gene_family_primary\n")
            for node in sorted(primary):
                row = by_node[node]
                status, _ = dispositions[node]
                fh.write(
                    f'S\t{node}\t*\tLN:i:{row["span"]}\tGF:Z:{family}\t'
                    f'DP:Z:{status}\tLO:Z:{row["chrom"]}:{row["start0"]}-{row["end"]}\n'
                )
            for a, b in sorted(primary_edges):
                evidence = "RNA_DNA" if (a, b) in both else "RNA_ONLY" if (a, b) in rna_edges else "DNA_ONLY"
                fh.write(f"L\t{a}\t+\t{b}\t+\t0M\tEV:Z:{evidence}\tMB:i:1\n")
        with (graph_dir / f"{family}.audit.gfa").open("w") as fh:
            fh.write("H\tVN:Z:1.0\tTS:Z:typed_gene_family_audit\n")
            for node in sorted(nodes):
                row = by_node[node]
                status, _ = dispositions[node]
                fh.write(
                    f'S\t{node}\t*\tLN:i:{row["span"]}\tGF:Z:{family}\t'
                    f'DP:Z:{status}\tLO:Z:{row["chrom"]}:{row["start0"]}-{row["end"]}\n'
                )
            for a, b in sorted(joint):
                evidence = "RNA_DNA" if (a, b) in both else "RNA_ONLY" if (a, b) in rna_edges else "DNA_ONLY"
                membership = int((a, b) in primary_edges)
                fh.write(f"L\t{a}\t+\t{b}\t+\t0M\tEV:Z:{evidence}\tMB:i:{membership}\n")
        with (graph_dir / f"{family}.colours.csv").open("w") as fh:
            fh.write("Node,Colour\n")
            for node in sorted(primary):
                status, _ = dispositions[node]
                fh.write(f"{node},{DISPOSITION_COLOURS[status]}\n")
        with (graph_dir / f"{family}.audit.colours.csv").open("w") as fh:
            fh.write("Node,Colour\n")
            for node in sorted(nodes):
                status, _ = dispositions[node]
                fh.write(f"{node},{DISPOSITION_COLOURS[status]}\n")

    write_tsv(
        out_dir / "known_family_certificates.tsv",
        family_rows,
        list(family_rows[0]),
    )
    write_tsv(
        out_dir / "known_family_nodes.tsv",
        node_rows_out,
        list(node_rows_out[0]),
    )
    return family_rows


def parse_attributes(value: str) -> dict[str, str]:
    out = {}
    for item in value.split(";"):
        if "=" in item:
            key, val = item.split("=", 1)
            out[key] = urllib.parse.unquote(val)
    return out


def load_genes(gff: Path, wanted_chroms: set[str]) -> dict[str, list[dict[str, object]]]:
    genes: dict[str, list[dict[str, object]]] = defaultdict(list)
    with gff.open() as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 9 or fields[0] not in wanted_chroms or fields[2] != "gene":
                continue
            attrs = parse_attributes(fields[8])
            genes[fields[0]].append(
                {
                    "start": int(fields[3]) - 1,
                    "end": int(fields[4]),
                    "name": attrs.get("Name", attrs.get("gene", ".")),
                    "description": attrs.get("description", "."),
                    "biotype": attrs.get("gene_biotype", "."),
                }
            )
    return genes


def candidate_audit(candidate_tsv: Path, annotation: Path, out_dir: Path) -> list[dict[str, object]]:
    rows = read_tsv(candidate_tsv)
    parsed = []
    for row in rows:
        chrom, rest = row["locus"].split(":", 1)
        start, end = map(int, rest.split("-"))
        parsed.append((row, chrom, start, end))
    genes = load_genes(annotation, {chrom for _, chrom, _, _ in parsed})
    answer: list[dict[str, object]] = []

    for row, chrom, start, end in parsed:
        family = row["best_soto_family"]
        stems = CANDIDATE_FAMILY_STEMS.get(family, ())
        overlaps = []
        for gene in genes.get(chrom, []):
            overlap = max(0, min(end, int(gene["end"])) - max(start, int(gene["start"])))
            if overlap:
                overlaps.append((overlap / max(end - start, 1), gene))
        overlaps.sort(key=lambda value: -value[0])
        annotation_text = " ".join(
            f'{gene["name"]} {gene["description"]}'.upper() for _, gene in overlaps
        )
        matching = [
            (coverage, gene)
            for coverage, gene in overlaps
            if coverage >= 0.50 and any(stem in f'{gene["name"]} {gene["description"]}'.upper() for stem in stems)
        ]
        conflicting = [
            (coverage, gene)
            for coverage, gene in overlaps
            if coverage >= 0.50 and gene["biotype"] == "protein_coding" and not any(
                stem in f'{gene["name"]} {gene["description"]}'.upper() for stem in stems
            )
        ]
        identity = float(row["best_soto_id"])
        coverage = float(row["best_soto_cov"])
        paralogs = int(row["n_genomic_paralog_loci"])

        if paralogs < 2 or "single-copy" in row["class"]:
            disposition = "REJECT_NOT_MULTICOPY"
            reason = "whole-genome check found only one locus"
        elif coverage > 1.0:
            disposition = "REVIEW_INVALID_LEGACY_COVERAGE"
            reason = "legacy coverage exceeds 1.0 and must be recomputed before assignment"
        elif identity < 0.60 or coverage < 0.50:
            disposition = "REVIEW_PARTIAL_SD_HOMOLOGY"
            reason = "real duplicated locus, but it does not pass the O1 family edge floor"
        elif not stems:
            disposition = "SD_EXTENSION_ONLY"
            reason = "the Soto SD cluster is gene-family-heterogeneous or untyped"
        elif matching and any(
            gene["biotype"] in {"protein_coding", "pseudogene", "transcribed_pseudogene"}
            for _, gene in matching
        ):
            disposition = "SUPPORTED_SOTO_MISSING_GENE_COPY"
            reason = "SD homology, multi-copy genome evidence, and independent gene-family annotation agree"
        elif matching:
            disposition = "SUPPORTED_NONCODING_FAMILY_LOCUS"
            reason = "same-family annotation agrees, but the record is noncoding and is not counted as a gene copy"
        elif conflicting:
            disposition = "REJECT_GENE_FAMILY_ASSIGNMENT"
            reason = "the locus is independently annotated as a different protein-coding gene"
        elif not overlaps:
            disposition = "NOVEL_GENE_COPY_CANDIDATE"
            reason = "coherent-family homology and genomic paralogs, but no independent gene annotation"
        else:
            disposition = "SD_EXTENSION_ONLY"
            reason = "SD homology lacks an independent same-family gene annotation"

        best_annotation = overlaps[0][1] if overlaps else None
        matched_annotation = matching[0][1] if matching else None
        copy_key = (
            str(matched_annotation["name"])
            if matched_annotation is not None
            else row["locus"]
        )
        answer.append(
            {
                "locus": row["locus"],
                "soto_sd_family": family,
                "gene_family_stems": ",".join(stems) or ".",
                "genomic_paralog_loci": paralogs,
                "identity": f"{identity:.3f}",
                "coverage_of_shorter": f"{coverage:.3f}",
                "best_annotation": best_annotation["name"] if best_annotation else ".",
                "best_annotation_biotype": best_annotation["biotype"] if best_annotation else ".",
                "independent_family_annotation": matched_annotation["name"] if matched_annotation else ".",
                "copy_key": copy_key,
                "disposition": disposition,
                "reason": reason,
            }
        )

    write_tsv(out_dir / "soto_missing_candidate_audit.tsv", answer, list(answer[0]))
    return answer


def write_report(
    out_dir: Path,
    family_rows: list[dict[str, object]],
    candidate_rows: list[dict[str, object]],
) -> None:
    supported_copy_keys = {
        row["copy_key"]
        for row in candidate_rows
        if row["disposition"] == "SUPPORTED_SOTO_MISSING_GENE_COPY"
    }
    supported_rows = []
    seen_copy_keys = set()
    for row in candidate_rows:
        if row["disposition"] != "SUPPORTED_SOTO_MISSING_GENE_COPY" or row["copy_key"] in seen_copy_keys:
            continue
        seen_copy_keys.add(row["copy_key"])
        supported_rows.append(row)
    novel = [row for row in candidate_rows if row["disposition"] == "NOVEL_GENE_COPY_CANDIDATE"]
    removed = sum(int(row["obvious_fp_removed"]) for row in family_rows)
    conflicts = sum(int(row["conflicting_gene_review"]) for row in family_rows)
    lines = [
        "# O1 typed gene-family audit",
        "",
        "This audit treats Soto families as segmental-duplication discovery strata, not as gene-family truth.",
        "RNA forward homology supplies primary membership; DNA supplies corroboration and can rescue an",
        "independently named RNA-null member, but DNA-only SD homology cannot recruit an anonymous locus.",
        "",
        "## Known multi-copy families",
        "",
        "| family | named members | supported | primary graph (nodes/edges/components) | RNA edges | DNA edges | RNA-only / both / DNA-only | status |",
        "|---|---:|---:|---:|---:|---:|---:|---|",
    ]
    for row in family_rows:
        lines.append(
            f'| {row["family"]} | {row["named_gene_members"]} | {row["known_supported"]} | '
            f'{row["primary_nodes"]}/{row["primary_edges"]}/{row["primary_components"]} | '
            f'{row["rna_forward_edges"]} | {row["dna_edges"]} | '
            f'{row["edges_rna_only"]} / {row["edges_both"]} / {row["edges_dna_only"]} | {row["status"]} |'
        )
    lines += [
        "",
        f"The panel removes {removed} adjudicated contaminants and withholds {conflicts} independently",
        "conflicting gene annotations from primary membership. GSTM3 and APOBEC3H remain named",
        "members but are reported as graph false negatives rather than being relabeled or deleted.",
        "HERC2 has a connected ten-node candidate graph but only one independently named member in",
        "this annotation, so it is annotation-limited rather than counted as a validated pass.",
        "",
        "## Soto-missing loci",
        "",
        f"Independent annotation supports {len(supported_copy_keys)} unique Soto-missing gene copies.",
        f"A further {len(novel)} loci are novel gene-copy candidates, not confirmed copies, because they",
        "lack independent same-family annotation.",
        "",
        "| locus | SD family | independent annotation | disposition |",
        "|---|---|---|---|",
    ]
    highlighted = supported_rows + [
        row for row in candidate_rows if row["disposition"] in {
            "SUPPORTED_NONCODING_FAMILY_LOCUS",
            "NOVEL_GENE_COPY_CANDIDATE",
            "REJECT_GENE_FAMILY_ASSIGNMENT",
            "REJECT_NOT_MULTICOPY",
        }
    ]
    for row in highlighted:
            lines.append(
                f'| {row["locus"]} | {row["soto_sd_family"]} | '
                f'{row["independent_family_annotation"]} | {row["disposition"]} |'
            )
    lines += [
        "",
        "The full TSV retains SD-only and partial-homology rows so they remain reviewable without",
        "inflating the gene-family copy count. `<family>.gene_family.gfa` contains only primary members",
        "and membership edges; `<family>.audit.gfa` retains every SD-discovery node for inspection. GFA",
        "link tag `EV` records RNA/DNA evidence and audit tag `MB:i:1` marks admitted edges. The companion",
        "Bandage colour CSV uses green for known members, blue for RNA-recruited candidates, orange for",
        "known false negatives, grey/",
        "purple for review nodes, and red for adjudicated contaminants.",
        "",
    ]
    (out_dir / "README.md").write_text("\n".join(lines))


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--panel-root",
        type=Path,
        default=Path("/home/juanfra/winloci_scratch/o1_closure"),
    )
    parser.add_argument(
        "--orientation-root",
        type=Path,
        default=Path("/tmp/rustle_o1_orientation"),
        help="corrected panel PAFs produced by bench/o1_orientation_gate.py",
    )
    parser.add_argument(
        "--candidates",
        type=Path,
        default=Path("bench/soto/candidate18_classification.tsv"),
    )
    parser.add_argument(
        "--annotation",
        type=Path,
        default=Path("../Reference/HSA_genomic.gff"),
    )
    parser.add_argument(
        "--out-dir",
        type=Path,
        default=Path("bench/o1_gene_family_audit"),
    )
    args = parser.parse_args()
    args.out_dir.mkdir(parents=True, exist_ok=True)
    family_rows = panel_audit(args.panel_root, args.orientation_root, args.out_dir)
    candidate_rows = candidate_audit(args.candidates, args.annotation, args.out_dir)
    write_report(args.out_dir, family_rows, candidate_rows)
    print(f"wrote typed O1 audit to {args.out_dir}")


if __name__ == "__main__":
    main()
