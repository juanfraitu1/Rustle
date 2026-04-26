#!/usr/bin/env python3
"""Phase 3 conceptual sketch: variation graph for a 2-copy gene family.

Hand-coded illustration of where the project is going. Two copies share a
linear "spine" with bubbles at copy-discriminating SNP positions. A single
read carrying SNP variants is shown aligned to the graph, settling on
copy A's path.

Usage:
    render_phase3_vg_sketch.py --output-png out.png
"""
from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import matplotlib.patches as patches


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--output-png", required=True, type=Path)
    args = ap.parse_args()

    fig, ax = plt.subplots(figsize=(11, 4.5))
    ax.set_xlim(0, 12)
    ax.set_ylim(-2.5, 3)
    ax.axis("off")

    # Spine nodes (shared)
    spine_y = 0
    spine_xs = [1, 3, 5.5, 8, 10.5]
    spine_labels = ["E1", "E2-pre", "E2-post", "E3", "E4"]
    for x, lbl in zip(spine_xs, spine_labels):
        c = patches.FancyBboxPatch((x - 0.45, spine_y - 0.3), 0.9, 0.6,
                                   boxstyle="round,pad=0.02", linewidth=1.2,
                                   edgecolor="#1F4E79", facecolor="#E8F1F8")
        ax.add_patch(c)
        ax.text(x, spine_y, lbl, ha="center", va="center", fontsize=9)

    # Spine edges
    for x1, x2 in zip(spine_xs, spine_xs[1:]):
        ax.annotate("", xy=(x2 - 0.5, 0), xytext=(x1 + 0.5, 0),
                    arrowprops=dict(arrowstyle="->", color="#1F4E79", lw=1.2))

    # Bubble at midpoint between E2-pre and E2-post: SNP G/T
    bx_mid = (spine_xs[1] + spine_xs[2]) / 2
    for dy, base, color, lbl in [(1.1, "G", "#2CA02C", "copy A"), (-1.1, "T", "#D62728", "copy B")]:
        c = patches.Circle((bx_mid, dy), 0.32, edgecolor=color, facecolor="white", linewidth=2)
        ax.add_patch(c)
        ax.text(bx_mid, dy, base, ha="center", va="center", fontsize=10, color=color, fontweight="bold")
        ax.text(bx_mid + 0.6, dy, lbl, ha="left", va="center", fontsize=8, color=color)
    # Bubble edges
    for dy in (1.1, -1.1):
        ax.annotate("", xy=(bx_mid - 0.32, dy * 0.7), xytext=(spine_xs[1] + 0.5, 0),
                    arrowprops=dict(arrowstyle="->", color="#888", lw=0.9))
        ax.annotate("", xy=(spine_xs[2] - 0.5, 0), xytext=(bx_mid + 0.32, dy * 0.7),
                    arrowprops=dict(arrowstyle="->", color="#888", lw=0.9))

    # Bubble at midpoint between E3 and E4: SNP C/A
    bx2 = (spine_xs[3] + spine_xs[4]) / 2
    for dy, base, color in [(1.1, "C", "#2CA02C"), (-1.1, "A", "#D62728")]:
        c = patches.Circle((bx2, dy), 0.32, edgecolor=color, facecolor="white", linewidth=2)
        ax.add_patch(c)
        ax.text(bx2, dy, base, ha="center", va="center", fontsize=10, color=color, fontweight="bold")
    for dy in (1.1, -1.1):
        ax.annotate("", xy=(bx2 - 0.32, dy * 0.7), xytext=(spine_xs[3] + 0.5, 0),
                    arrowprops=dict(arrowstyle="->", color="#888", lw=0.9))
        ax.annotate("", xy=(spine_xs[4] - 0.5, 0), xytext=(bx2 + 0.32, dy * 0.7),
                    arrowprops=dict(arrowstyle="->", color="#888", lw=0.9))

    # Read aligned on copy A path (top)
    read_y = 2.3
    ax.plot([spine_xs[0], spine_xs[-1]], [read_y, read_y], color="#E37222", lw=4, solid_capstyle="round")
    ax.text(spine_xs[0] - 0.6, read_y, "read", ha="right", va="center", fontsize=9, color="#E37222")
    # Mark read SNP calls
    ax.text(bx_mid, read_y - 0.3, "G", ha="center", color="#2CA02C", fontsize=9, fontweight="bold")
    ax.text(bx2, read_y - 0.3, "C", ha="center", color="#2CA02C", fontsize=9, fontweight="bold")
    # Annotate
    ax.text(6, read_y + 0.45, "SNPs G, C → copy A path", ha="center", fontsize=10, color="#2CA02C")

    ax.text(6, -2.2, "Future: align reads to family variation graph; resolve to specific copy by SNP support",
            ha="center", fontsize=10, style="italic", color="#555")

    fig.tight_layout()
    fig.savefig(args.output_png, dpi=150, bbox_inches="tight")


if __name__ == "__main__":
    main()
