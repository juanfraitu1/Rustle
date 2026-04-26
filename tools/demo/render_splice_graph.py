#!/usr/bin/env python3
"""Render a single splice graph from RUSTLE_PARITY_GRAPH_TSV + RUSTLE_PARITY_GRAPH_EDGES_TSV.

Mechanical layout via graphviz dot. No hand-coded positioning. Output is PNG
plus the underlying .dot file (for reproducibility — paste into a graphviz
viewer if PNG is contested).

Usage:
    render_splice_graph.py --nodes-tsv ... --edges-tsv ... \\
        --output-png out.png --output-dot out.dot \\
        [--locus chrom:start-end] [--title "Title"]
"""
from __future__ import annotations

import argparse
import subprocess
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from demo.common import (
    PALETTE,
    GraphEdge,
    GraphNode,
    group_by_bundle,
    read_edges_tsv,
    read_nodes_tsv,
)


def _parse_locus(s: str | None) -> tuple[str, int, int] | None:
    if not s:
        return None
    chrom, rng = s.split(":")
    start, end = rng.split("-")
    return (chrom, int(start), int(end))


def _node_label(n: GraphNode) -> str:
    badge = ""
    if n.hardstart:
        badge += "▶"
    if n.hardend:
        badge += "◀"
    return f"#{n.node_idx}{badge}\\n[{n.start}-{n.end}]\\ncov={n.cov:.1f}"


def _edge_color(e: GraphEdge) -> str:
    return PALETTE["edge_junction"] if e.edge_kind == "junction" else PALETTE["edge_collinear"]


def _edge_label(e: GraphEdge) -> str:
    return f"{e.edge_kind}\\nbn={e.bottleneck_cov:.1f}"


def to_dot(nodes: list[GraphNode], edges: list[GraphEdge], title: str = "") -> str:
    out = ['digraph G {']
    out.append('  rankdir=LR;')
    out.append('  node [fontname="Helvetica" fontsize=10];')
    out.append('  edge [fontname="Helvetica" fontsize=9];')
    if title:
        out.append(f'  labelloc="t"; label="{title}"; fontsize=14;')

    for n in nodes:
        fill = PALETTE["node_source_sink"] if (n.hardstart or n.hardend) else PALETTE["node_default"]
        label = _node_label(n).replace('"', '\\"')
        out.append(
            f'  n{n.node_idx} [shape=box style="rounded,filled" fillcolor="{fill}" label="{label}"];'
        )
    for e in edges:
        out.append(
            f'  n{e.from_idx} -> n{e.to_idx} '
            f'[color="{_edge_color(e)}" label="{_edge_label(e)}" '
            f'penwidth={1.0 + min(4.0, e.bottleneck_cov / 5.0):.1f}];'
        )
    out.append('}')
    return "\n".join(out)


def render(dot_text: str, output_png: Path, output_dot: Path | None) -> None:
    if output_dot:
        output_dot.write_text(dot_text)
    p = subprocess.run(
        ["dot", "-Tpng", "-o", str(output_png)],
        input=dot_text, text=True, capture_output=True,
    )
    if p.returncode != 0:
        raise SystemExit(f"graphviz failed: {p.stderr}")


def main(argv: list[str]) -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--nodes-tsv", required=True, type=Path)
    ap.add_argument("--edges-tsv", required=True, type=Path)
    ap.add_argument("--output-png", required=True, type=Path)
    ap.add_argument("--output-dot", type=Path, default=None)
    ap.add_argument("--locus", default=None)
    ap.add_argument("--title", default="")
    args = ap.parse_args(argv)

    locus = _parse_locus(args.locus)
    nodes = read_nodes_tsv(args.nodes_tsv, locus)
    edges = read_edges_tsv(args.edges_tsv, locus)

    if not nodes:
        raise SystemExit(f"No nodes match locus {args.locus} in {args.nodes_tsv}")

    bundles = list(group_by_bundle(nodes))
    if len(bundles) > 1:
        # If multiple bundles match, render only the first and warn.
        print(f"WARN: locus matched {len(bundles)} bundles; rendering first only", file=sys.stderr)
    (key, bnodes) = bundles[0]
    chrom, bs, be, _ = key
    bedges = [e for e in edges if e.chrom == chrom and e.bdstart == bs and e.bdend == be]

    title = args.title or f"{chrom}:{bs}-{be} ({len(bnodes)} nodes, {len(bedges)} edges)"
    dot = to_dot(bnodes, bedges, title=title)
    render(dot, args.output_png, args.output_dot)
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
