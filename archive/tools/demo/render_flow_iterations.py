#!/usr/bin/env python3
"""Per-iteration flow decomposition renderer.

For each iteration row in the flow TSV, render the splice graph with the
augmenting path highlighted in orange, capacity-after-iteration shown on
edges. Compose into a horizontal grid.

Usage:
    render_flow_iterations.py --flow-tsv ... --nodes-tsv ... --edges-tsv ...
        --locus chrom:start-end --output-png out.png
"""
from __future__ import annotations

import argparse
import csv
import subprocess
import sys
from dataclasses import dataclass
from pathlib import Path
from tempfile import TemporaryDirectory

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from demo.common import GraphEdge, GraphNode, group_by_bundle, read_edges_tsv, read_nodes_tsv, PALETTE
from demo.render_splice_graph import to_dot


@dataclass(frozen=True)
class FlowIter:
    chrom: str
    bdstart: int
    bdend: int
    bundle_run: int
    iter_idx: int
    path_nodes: list[str]
    bottleneck: float
    total_flow_after: float


def _parse_locus(s: str) -> tuple[str, int, int]:
    chrom, rng = s.split(":")
    start, end = rng.split("-")
    return (chrom, int(start), int(end))


def read_flow_tsv(path: Path, locus: tuple[str, int, int]) -> list[FlowIter]:
    out = []
    chrom_f, lstart, lend = locus
    with Path(path).open() as f:
        for row in csv.DictReader(f, delimiter="\t"):
            if row["chrom"] != chrom_f:
                continue
            bs, be = int(row["bdstart"]), int(row["bdend"])
            if be < lstart or bs > lend:
                continue
            out.append(FlowIter(
                chrom=row["chrom"], bdstart=bs, bdend=be,
                bundle_run=int(row["bundle_run"]),
                iter_idx=int(row["iter_idx"]),
                path_nodes=row["path_nodes"].split(","),
                bottleneck=float(row["bottleneck"]),
                total_flow_after=float(row["total_flow_after"]),
            ))
    return sorted(out, key=lambda r: (r.bundle_run, r.iter_idx))


def _highlighted_dot(nodes: list[GraphNode], edges: list[GraphEdge], path_node_strs: list[str], iter_idx: int, bottleneck: float, title_extra: str) -> str:
    """Same as render_splice_graph.to_dot but highlights edges along the augmenting path in orange."""
    # Identify path edge pairs: consecutive node strs (skipping S/T sentinels).
    path_idx_pairs = set()
    cleaned = [p for p in path_node_strs if p not in ("S", "T")]
    for a, b in zip(cleaned, cleaned[1:]):
        try:
            path_idx_pairs.add((int(a), int(b)))
        except ValueError:
            continue

    # Reuse to_dot but tag highlighted edges by intercepting:
    base = to_dot(nodes, edges, title=f"iter {iter_idx} — bottleneck={bottleneck:.2f}{title_extra}")
    out_lines = []
    for line in base.splitlines():
        emitted = False
        for (a, b) in path_idx_pairs:
            tag = f"n{a} -> n{b}"
            if line.lstrip().startswith(tag):
                # rewrite color + penwidth
                line = line.replace(
                    'color="' + PALETTE["edge_collinear"] + '"',
                    f'color="{PALETTE["edge_augmenting"]}"'
                ).replace(
                    'color="' + PALETTE["edge_junction"] + '"',
                    f'color="{PALETTE["edge_augmenting"]}"'
                )
                # Bump pen width
                if "penwidth=" in line:
                    pre, _, post = line.rpartition("penwidth=")
                    val_str = post.split("]")[0]
                    line = pre + "penwidth=4.0" + post[len(val_str):]
                emitted = True
                break
        out_lines.append(line)
    return "\n".join(out_lines)


def main(argv: list[str]) -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--flow-tsv", required=True, type=Path)
    ap.add_argument("--nodes-tsv", required=True, type=Path)
    ap.add_argument("--edges-tsv", required=True, type=Path)
    ap.add_argument("--locus", required=True)
    ap.add_argument("--output-png", required=True, type=Path)
    args = ap.parse_args(argv)

    locus = _parse_locus(args.locus)
    iters = read_flow_tsv(args.flow_tsv, locus)
    nodes = read_nodes_tsv(args.nodes_tsv, locus)
    edges = read_edges_tsv(args.edges_tsv, locus)
    if not iters:
        raise SystemExit(f"no flow iterations for locus {args.locus}")
    if not nodes:
        raise SystemExit(f"no nodes for locus {args.locus}")

    # Render one PNG per iteration, then stitch.
    import matplotlib.pyplot as plt
    from matplotlib.image import imread

    with TemporaryDirectory() as td:
        td_path = Path(td)
        pngs = []
        for it in iters:
            dot = _highlighted_dot(
                nodes, edges, it.path_nodes, it.iter_idx, it.bottleneck,
                title_extra=f"  total flow={it.total_flow_after:.2f}"
            )
            png = td_path / f"iter_{it.bundle_run:02d}_{it.iter_idx:02d}.png"
            p = subprocess.run(
                ["dot", "-Tpng", "-o", str(png)],
                input=dot, text=True, capture_output=True,
            )
            if p.returncode != 0:
                raise SystemExit(f"dot failed: {p.stderr}")
            pngs.append((it, png))

        cols = min(3, len(pngs))
        rows = (len(pngs) + cols - 1) // cols
        fig, axes = plt.subplots(rows, cols, figsize=(7 * cols, 5 * rows))
        if rows == 1 and cols == 1:
            axes = [[axes]]
        elif rows == 1:
            axes = [axes]
        elif cols == 1:
            axes = [[a] for a in axes]
        for k, (it, png) in enumerate(pngs):
            r, c = k // cols, k % cols
            ax = axes[r][c]
            ax.imshow(imread(png))
            ax.set_title(f"iter {it.iter_idx}: bn={it.bottleneck:.2f}, total={it.total_flow_after:.2f}", fontsize=10)
            ax.axis("off")
        # blank out unused subplots
        for k in range(len(pngs), rows * cols):
            r, c = k // cols, k % cols
            axes[r][c].axis("off")
        fig.suptitle(f"Edmonds-Karp iterations at {args.locus}", fontsize=12)
        fig.tight_layout()
        fig.savefig(args.output_png, dpi=150, bbox_inches="tight")
        plt.close(fig)
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
