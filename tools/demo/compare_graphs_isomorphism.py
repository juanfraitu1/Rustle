#!/usr/bin/env python3
"""Three-way splice-graph isomorphism comparator.

For each --locus given (>=2), build a NetworkX DiGraph from the matching nodes
and edges, then test pairwise isomorphism using VF2.

"Similar" definition: graph isomorphism preserving (a) node count, (b) edge
multiset by edge_kind, (c) source/sink topology (hardstart/hardend nodes).
Genomic coordinates normalized by offset within the locus before comparison.

Usage:
    compare_graphs_isomorphism.py \\
        --nodes-tsv ... --edges-tsv ... \\
        --locus chr19:104789647-104796276 \\
        --locus chr19:104830536-104837094 \\
        --locus chr19:104871356-104877901 \\
        --output-json out.json --output-png out.png
"""
from __future__ import annotations

import argparse
import json
import subprocess
import sys
from pathlib import Path
from tempfile import TemporaryDirectory

import networkx as nx
from networkx.algorithms.isomorphism import DiGraphMatcher

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from demo.common import (
    GraphEdge,
    GraphNode,
    group_by_bundle,
    read_edges_tsv,
    read_nodes_tsv,
)


def _parse_locus(s: str) -> tuple[str, int, int]:
    chrom, rng = s.split(":")
    start, end = rng.split("-")
    return (chrom, int(start), int(end))


def _build_nx(nodes: list[GraphNode], edges: list[GraphEdge]) -> nx.DiGraph:
    g = nx.DiGraph()
    if not nodes:
        return g
    bundles = list(group_by_bundle(nodes))
    (key, bnodes) = bundles[0]
    chrom, bs, be, _ = key
    bedges = [e for e in edges if e.chrom == chrom and e.bdstart == bs and e.bdend == be]
    # Normalize coordinates by subtracting bdstart so the same shape at different
    # genomic positions matches.
    for n in bnodes:
        g.add_node(
            n.node_idx,
            type="source" if n.hardstart else ("sink" if n.hardend else "internal"),
            length=n.end - n.start,
        )
    for e in bedges:
        g.add_edge(e.from_idx, e.to_idx, kind=e.edge_kind)
    return g


def _node_match(a: dict, b: dict) -> bool:
    return a.get("type") == b.get("type")


def _edge_match(a: dict, b: dict) -> bool:
    return a.get("kind") == b.get("kind")


def _kind_counts(g: nx.DiGraph) -> dict:
    c: dict[str, int] = {}
    for _, _, d in g.edges(data=True):
        c[d.get("kind", "?")] = c.get(d.get("kind", "?"), 0) + 1
    return c


def _compare_pair(g1: nx.DiGraph, g2: nx.DiGraph) -> dict:
    matcher = DiGraphMatcher(g1, g2, node_match=_node_match, edge_match=_edge_match)
    is_iso = matcher.is_isomorphic()
    out = {"isomorphic": is_iso}
    if is_iso:
        out["bijection"] = {str(k): v for k, v in matcher.mapping.items()}
    else:
        out["diff"] = {
            "node_counts": (g1.number_of_nodes(), g2.number_of_nodes()),
            "edge_counts": (g1.number_of_edges(), g2.number_of_edges()),
            "edge_kind_counts": (
                _kind_counts(g1), _kind_counts(g2),
            ),
        }
    return out


def _to_dot_for_graph(nodes: list[GraphNode], edges: list[GraphEdge], title: str) -> str:
    # Reuse render_splice_graph.to_dot for visual consistency.
    from demo.render_splice_graph import to_dot
    bundles = list(group_by_bundle(nodes))
    (key, bnodes) = bundles[0]
    chrom, bs, be, _ = key
    bedges = [e for e in edges if e.chrom == chrom and e.bdstart == bs and e.bdend == be]
    return to_dot(bnodes, bedges, title=title)


def _render_panel(loci, dots: list[str], output_png: Path, mappings: dict | None) -> None:
    """Compose N graphviz dot images into a horizontal panel via matplotlib."""
    import matplotlib.pyplot as plt
    from matplotlib.image import imread

    with TemporaryDirectory() as td:
        td_path = Path(td)
        pngs = []
        for i, dot in enumerate(dots):
            png = td_path / f"copy{i}.png"
            p = subprocess.run(
                ["dot", "-Tpng", "-o", str(png)],
                input=dot, text=True, capture_output=True,
            )
            if p.returncode != 0:
                raise SystemExit(f"graphviz failed for copy {i}: {p.stderr}")
            pngs.append(png)

        fig, axes = plt.subplots(1, len(pngs), figsize=(6 * len(pngs), 6))
        if len(pngs) == 1:
            axes = [axes]
        for ax, png, locus in zip(axes, pngs, loci):
            ax.imshow(imread(png))
            ax.set_title(f"{locus[0]}:{locus[1]}-{locus[2]}", fontsize=10)
            ax.axis("off")
        if mappings:
            note = "Bijection: " + "  /  ".join(
                f"{k}: " + ", ".join(f"{a}↔{b}" for a, b in sorted(v.items())[:5])
                + ("…" if len(v) > 5 else "")
                for k, v in mappings.items()
            )
            fig.suptitle(note, fontsize=10)
        fig.tight_layout()
        fig.savefig(output_png, dpi=150, bbox_inches="tight")
        plt.close(fig)


def main(argv: list[str]) -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--nodes-tsv", required=True, type=Path)
    ap.add_argument("--edges-tsv", required=True, type=Path)
    ap.add_argument("--locus", action="append", required=True,
                    help="chrom:start-end (specify >=2 times)")
    ap.add_argument("--output-json", required=True, type=Path)
    ap.add_argument("--output-png", required=True, type=Path)
    args = ap.parse_args(argv)

    if len(args.locus) < 2:
        raise SystemExit("need >=2 --locus arguments")

    loci = [_parse_locus(s) for s in args.locus]
    graphs = []
    dots = []
    for locus in loci:
        nodes = read_nodes_tsv(args.nodes_tsv, locus)
        edges = read_edges_tsv(args.edges_tsv, locus)
        if not nodes:
            raise SystemExit(f"no nodes for {locus}")
        g = _build_nx(nodes, edges)
        graphs.append(g)
        dots.append(_to_dot_for_graph(nodes, edges, title=f"{locus[0]}:{locus[1]}-{locus[2]}"))

    # Compare all pairs (1<->2, 1<->3, etc.)
    pairwise = {}
    bijections = {}
    for i in range(len(graphs)):
        for j in range(i + 1, len(graphs)):
            label = f"{i+1}↔{j+1}"
            res = _compare_pair(graphs[i], graphs[j])
            pairwise[label] = res
            if res.get("isomorphic"):
                bijections[label] = res["bijection"]

    all_iso = all(p["isomorphic"] for p in pairwise.values())
    out = {
        "isomorphic": all_iso,
        "pairwise": pairwise,
        "bijection": bijections if all_iso else {},
        "graph_summary": [
            {"locus": f"{c}:{s}-{e}", "n_nodes": g.number_of_nodes(),
             "n_edges": g.number_of_edges(),
             "edge_kinds": _kind_counts(g)}
            for (c, s, e), g in zip(loci, graphs)
        ],
    }
    if not all_iso:
        out["diff"] = {k: v.get("diff") for k, v in pairwise.items() if not v["isomorphic"]}

    args.output_json.write_text(json.dumps(out, indent=2))
    _render_panel(loci, dots, args.output_png, bijections if all_iso else None)
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
