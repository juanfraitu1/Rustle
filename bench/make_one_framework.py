#!/usr/bin/env python
"""One-framework slide: the family graph, the variation graph, copies/isoforms, and the conflict graph are
ONE object seen at different scopes/views. Hero = the variation graph; satellites = its scope (family graph)
and its read-view (conflict graph); the two path-axes (bubbles=copies, junctions=isoforms) labelled on it.
Output: bench/slides/one_framework.png
Run: /home/juanfra/miniforge3/bin/python bench/make_one_framework.py
"""
import os
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Circle, FancyBboxPatch, FancyArrowPatch, Rectangle

OUT = os.path.join(os.path.dirname(__file__), "slides")
os.makedirs(OUT, exist_ok=True)
plt.rcParams.update({"font.family": "DejaVu Sans"})
NAVY = "#14285a"; BLUE = "#377eb8"; RED = "#e41a1c"; GREEN = "#4daf4a"
GREY = "#8a8f98"; ORANGE = "#e08a00"; DARK = "#222222"; SPINE = "#c9c9c9"


def nd(ax, x, y, label, fc, r=0.24, tc="white", fs=10):
    ax.add_patch(Circle((x, y), r, fc=fc, ec="black", lw=1.2, zorder=5))
    ax.text(x, y, label, ha="center", va="center", fontsize=fs, weight="bold", color=tc, zorder=6)


def build():
    fig, ax = plt.subplots(figsize=(14, 7.6))
    ax.axis("off"); ax.set_xlim(0, 14); ax.set_ylim(0, 7.6)

    # ============== HERO: the variation graph ==============
    ax.text(7.0, 7.15, "TWO graphs, one arrow  —  loci  →  sequence", ha="center", fontsize=18, weight="bold", color=NAVY)
    ax.text(7.0, 6.75, "GRAPH B  ·  the variation graph over SEQUENCE (one per family)", ha="center", fontsize=10.5, style="italic", color="#666")
    sx0, sx1, sy = 4.7, 9.4, 4.25
    ax.add_patch(Rectangle((sx0, sy - 0.09), sx1 - sx0, 0.18, fc=SPINE, ec="none", zorder=1))
    ax.text((sx0 + sx1) / 2, sy - 0.42, "shared backbone (spine)", ha="center", fontsize=8.5, color="#999")
    # PSV bubbles  (copy axis)
    bx = [5.6, 6.9]
    top = [("A", BLUE), ("C", BLUE)]; bot = [("G", RED), ("T", RED)]
    for i, x in enumerate(bx):
        ax.plot([x, x], [sy + 0.09, 5.2], color="#c4c4c4", lw=1, zorder=1)
        ax.plot([x, x], [sy - 0.09, 3.3], color="#c4c4c4", lw=1, zorder=1)
        nd(ax, x, 5.2, top[i][0], top[i][1]); nd(ax, x, 3.3, bot[i][0], bot[i][1])
    # two copy-paths threading the bubbles
    ax.plot([sx0 + 0.05, 5.6, 6.9, 8.0], [sy, 5.2, 5.2, sy], color=BLUE, lw=2.6, alpha=.85, zorder=3)
    ax.plot([sx0 + 0.05, 5.6, 6.9, 8.0], [sy, 3.3, 3.3, sy], color=RED, lw=2.6, alpha=.85, zorder=3)
    ax.text(6.3, 5.95, "PSV bubbles  →  COPY", ha="center", fontsize=12, weight="bold", color=NAVY)
    ax.text(6.3, 5.62, "(which paralog — a path through the bubbles)", ha="center", fontsize=9, style="italic", color="#666")
    # junction (isoform axis): an exon that a read either includes or skips
    ax.add_patch(Rectangle((8.05, sy - 0.05), 0.85, 0.10, fc="#8fbf8f", ec="black", lw=0.6, zorder=2))
    ax.text(8.47, sy + 0.22, "exon", ha="center", fontsize=8, color="#2e7d32")
    ax.add_patch(FancyArrowPatch((8.05, sy - 0.02), (8.9, sy - 0.02), connectionstyle="arc3,rad=0.9",
                                 arrowstyle="-", lw=2.0, ls=(0, (4, 2)), color=ORANGE, zorder=2))
    ax.text(8.47, 2.95, "junction  →  ISOFORM", ha="center", fontsize=11, weight="bold", color=ORANGE)
    ax.text(8.47, 2.63, "(which transcript —\na path through the junctions)", ha="center", fontsize=8.6, style="italic", color="#8a6d00")

    # ============== LEFT satellite: FAMILY GRAPH (scope) ==============
    loci = {"a": (1.0, 4.9), "b": (2.0, 5.2), "c": (2.15, 4.1), "d": (1.05, 3.75)}
    for p, q in [("a", "b"), ("b", "c"), ("c", "d"), ("a", "c")]:
        ax.plot([loci[p][0], loci[q][0]], [loci[p][1], loci[q][1]], color=GREEN, lw=1.8, zorder=2)
    for k, (x, y) in loci.items():
        ax.add_patch(Circle((x, y), 0.2, fc=GREEN, ec="black", lw=1, zorder=3))
    ax.text(1.55, 5.75, "GRAPH A · over LOCI", ha="center", fontsize=11, weight="bold", color=GREEN)
    ax.text(1.55, 3.25, "homology: which loci\nare one family\n(γ-quasi-clique, O1)", ha="center", va="top", fontsize=8.6, color="#555")
    ax.add_patch(FancyArrowPatch((2.6, 4.5), (4.5, 4.4), arrowstyle="-|>", mutation_scale=18, lw=2, color=GREY))
    ax.text(3.55, 4.9, "which loci share", ha="center", fontsize=9, weight="bold", color=DARK)
    ax.text(3.55, 4.66, "a variation graph", ha="center", fontsize=9, weight="bold", color=DARK)

    # ====== RIGHT satellite: COUNTING COPIES = colour the read-conflict graph (a VIEW of graph B) ======
    ax.add_patch(FancyArrowPatch((9.5, 4.4), (11.05, 4.4), arrowstyle="-|>", mutation_scale=18, lw=2, color=GREY))
    ax.text(10.28, 4.9, "a VIEW of", ha="center", fontsize=9, weight="bold", color=DARK)
    ax.text(10.28, 4.66, "graph B", ha="center", fontsize=9, weight="bold", color=DARK)
    # reads coloured by the copy (colour class) they belong to
    ax.plot([11.45, 13.65], [4.95, 4.95], color=BLUE, lw=1.1, ls=(0, (2, 2)), alpha=.5, zorder=2)
    ax.plot([11.45, 13.65], [3.95, 3.95], color=RED, lw=1.1, ls=(0, (2, 2)), alpha=.5, zorder=2)
    for x0, x1, y, c in [(11.6, 12.4, 4.95, BLUE), (12.75, 13.55, 4.95, BLUE), (11.6, 12.4, 3.95, RED), (12.75, 13.55, 3.95, RED)]:
        ax.plot([x0, x1], [y, y], color=c, lw=5, alpha=.7, solid_capstyle="round", zorder=3)
    ax.text(12.55, 5.75, "COUNTING COPIES (χ_H)", ha="center", fontsize=10.5, weight="bold", color=NAVY)
    ax.text(12.55, 3.35, "colour the read-conflict graph:\ncopy number = χ_H = fewest colours\n(MCC = χ(H); a lower bound)", ha="center", va="top", fontsize=8.4, color="#555")

    # ============== unifying caption ==============
    ax.add_patch(FancyBboxPatch((0.5, 0.55), 13.0, 1.35, boxstyle="round,pad=0.02,rounding_size=0.05",
                                fc="#f4f7fb", ec=NAVY, lw=1.4))
    ax.text(7.0, 1.5, "TWO graphs, one arrow.", ha="center", fontsize=13, weight="bold", color=NAVY)
    ax.text(7.0, 1.08, "GRAPH A over LOCI (homology → which loci are one family)     →     GRAPH B over SEQUENCE (the variation graph → how the copies differ)",
            ha="center", fontsize=9.4, color=DARK)
    ax.text(7.0, 0.72, "The PSV matrix, χ_H, and copies-as-paths are all VIEWS of graph B — not separate graphs.", ha="center", fontsize=10.5, weight="bold", style="italic", color=GREEN)

    p = os.path.join(OUT, "one_framework.png")
    fig.savefig(p, bbox_inches="tight", facecolor="white", dpi=150)
    plt.close()
    return p


if __name__ == "__main__":
    print("wrote", build())
