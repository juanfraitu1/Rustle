"""Shared utilities for demo visualizers.

All renderers consume the same TSV schemas. This module centralizes the schemas,
palette, and graphviz wrapping so the renderers can stay focused.
"""
from __future__ import annotations

import csv
from dataclasses import dataclass
from pathlib import Path
from typing import Iterator, Sequence

# Color palette — keep accessible (deuteranopia-safe).
PALETTE = {
    "node_default":     "#E8F1F8",
    "node_source_sink": "#F8E5C8",
    "edge_collinear":   "#1F4E79",
    "edge_junction":    "#C0504D",
    "edge_residual":    "#9D9D9D",
    "edge_augmenting":  "#E37222",
    "iso_pair":         ["#1F77B4", "#2CA02C", "#D62728", "#9467BD", "#FF7F0E"],
}


@dataclass(frozen=True)
class GraphNode:
    bdstart: int
    bdend: int
    strand: str
    chrom: str
    node_idx: int
    start: int
    end: int
    cov: float
    hardstart: bool
    hardend: bool
    nchildren: int
    nparents: int


@dataclass(frozen=True)
class GraphEdge:
    bdstart: int
    bdend: int
    strand: str
    chrom: str
    from_idx: int
    to_idx: int
    from_start: int
    from_end: int
    to_start: int
    to_end: int
    edge_kind: str   # "collinear" or "junction"
    from_cov: float
    to_cov: float
    bottleneck_cov: float


def _bundle_match(row: dict, locus: tuple[str, int, int] | None) -> bool:
    if locus is None:
        return True
    chrom, start, end = locus
    if row["chrom"] != chrom:
        return False
    bs, be = int(row["bdstart"]), int(row["bdend"])
    return be >= start and bs <= end


def read_nodes_tsv(path: Path, locus: tuple[str, int, int] | None = None) -> list[GraphNode]:
    out: list[GraphNode] = []
    with Path(path).open() as f:
        for row in csv.DictReader(f, delimiter="\t"):
            if not _bundle_match(row, locus):
                continue
            out.append(GraphNode(
                bdstart=int(row["bdstart"]),
                bdend=int(row["bdend"]),
                strand=row["strand"],
                chrom=row["chrom"],
                node_idx=int(row["node_idx"]),
                start=int(row["start"]),
                end=int(row["end"]),
                cov=float(row["cov"]),
                hardstart=row["hardstart"] in ("1", "true", "True"),
                hardend=row["hardend"] in ("1", "true", "True"),
                nchildren=int(row["nchildren"]),
                nparents=int(row["nparents"]),
            ))
    return out


def read_edges_tsv(path: Path, locus: tuple[str, int, int] | None = None) -> list[GraphEdge]:
    out: list[GraphEdge] = []
    with Path(path).open() as f:
        for row in csv.DictReader(f, delimiter="\t"):
            if not _bundle_match(row, locus):
                continue
            out.append(GraphEdge(
                bdstart=int(row["bdstart"]),
                bdend=int(row["bdend"]),
                strand=row["strand"],
                chrom=row["chrom"],
                from_idx=int(row["from_idx"]),
                to_idx=int(row["to_idx"]),
                from_start=int(row["from_start"]),
                from_end=int(row["from_end"]),
                to_start=int(row["to_start"]),
                to_end=int(row["to_end"]),
                edge_kind=row["edge_kind"],
                from_cov=float(row["from_cov"]),
                to_cov=float(row["to_cov"]),
                bottleneck_cov=float(row["bottleneck_cov"]),
            ))
    return out


def group_by_bundle(nodes: Sequence[GraphNode]) -> Iterator[tuple[tuple[str, int, int, str], list[GraphNode]]]:
    """Yield ((chrom, bdstart, bdend, strand), nodes_in_bundle)."""
    keyed: dict[tuple, list[GraphNode]] = {}
    for n in nodes:
        keyed.setdefault((n.chrom, n.bdstart, n.bdend, n.strand), []).append(n)
    for k, v in keyed.items():
        yield k, sorted(v, key=lambda n: n.node_idx)
