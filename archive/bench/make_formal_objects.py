#!/usr/bin/env python
"""Backup slide figure — the two formal objects: family = gamma-quasi-clique, copy number = chi(H).
Output: bench/slides/formal_objects.png
Run: /home/juanfra/miniforge3/bin/python bench/make_formal_objects.py
"""
import os
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Circle, FancyBboxPatch, FancyArrowPatch

OUT = os.path.join(os.path.dirname(__file__), "slides")
os.makedirs(OUT, exist_ok=True)
plt.rcParams.update({"font.family": "DejaVu Sans"})
GREEN, BLUE, RED, GREY = "#4daf4a", "#377eb8", "#e41a1c", "#999999"
NAVY = "#14285a"; DARK = "#222222"; YELLOW = "#ffe680"


def node(ax, x, y, label, fc, r=0.34, tc="white"):
    ax.add_patch(Circle((x, y), r, fc=fc, ec="black", lw=1.4, zorder=4))
    ax.text(x, y, label, ha="center", va="center", fontsize=11, weight="bold", color=tc, zorder=5)


def edge(ax, p, q, color="#555", lw=2.2, ls="-"):
    ax.plot([p[0], q[0]], [p[1], q[1]], color=color, lw=lw, ls=ls, zorder=2)


def build():
    fig = plt.figure(figsize=(13.2, 5.7))
    axL = fig.add_subplot(1, 2, 1); axR = fig.add_subplot(1, 2, 2)
    for ax in (axL, axR):
        ax.axis("off"); ax.set_xlim(0, 10); ax.set_ylim(0, 6)

    # ================= LEFT: FAMILY = gamma-quasi-clique =================
    axL.add_patch(FancyBboxPatch((0.6, 5.4), 8.8, 0.55, boxstyle="round,pad=0.02,rounding_size=0.05", fc=GREEN, ec=GREEN))
    axL.text(5, 5.67, "FAMILY  =  γ-quasi-clique", ha="center", va="center", fontsize=14, weight="bold", color="white")
    axL.text(5, 5.22, "a dense subgraph of the locus homology graph", ha="center", fontsize=9.7, style="italic", color="#666")
    # -- the family: a dense cluster (5 of 6 edges; A-D missing) --
    fam = {"A": (1.55, 4.3), "B": (3.15, 4.55), "C": (3.35, 3.15), "D": (1.4, 2.85)}
    for a, b in [("A", "B"), ("B", "C"), ("C", "D"), ("A", "C"), ("B", "D")]:
        edge(axL, fam[a], fam[b], color=BLUE, lw=2.4)
    axL.add_patch(FancyBboxPatch((0.85, 2.4), 3.2, 2.5, boxstyle="round,pad=0.02,rounding_size=0.04",
                                 fc="none", ec=GREEN, lw=2.2, ls=(0, (4, 2))))
    for k, (x, y) in fam.items():
        node(axL, x, y, k, BLUE)
    axL.text(2.45, 1.95, "one FAMILY: dense —\n≥ γ of ties present\n(A–D missing, still one)",
             ha="center", va="top", fontsize=9.6, color=DARK, weight="bold")
    # -- the over-merge case: a chain through a shared repeat --
    ch = {"X": (5.75, 4.35), "M": (7.15, 3.55), "Z": (8.55, 4.35)}
    edge(axL, ch["X"], ch["M"], color=GREY, lw=2.4); edge(axL, ch["M"], ch["Z"], color=GREY, lw=2.4)
    node(axL, *ch["X"], "X", GREY); node(axL, *ch["Z"], "Z", GREY)
    axL.add_patch(Circle(ch["M"], 0.34, fc=YELLOW, ec="black", lw=1.4, zorder=4))
    axL.text(*ch["M"], "M", ha="center", va="center", fontsize=11, weight="bold", color="black", zorder=5)
    axL.text(7.15, 3.0, "shared exon / repeat", ha="center", fontsize=8.4, color="#8a6d00", style="italic")
    axL.text(7.15, 1.95, "single-linkage chains X–Z\nvia M → over-merge\n(density ≥ γ blocks it)",
             ha="center", va="top", fontsize=9.6, color=RED, weight="bold")
    axL.text(5, 0.35, "γ = 1: strict clique   ·   γ → 0: component   ·   we use γ ≈ 0.4",
             ha="center", fontsize=9.0, style="italic", color="#777")

    # ============ RIGHT: COPY NUMBER = chromatic number of the read-conflict graph ============
    axR.add_patch(FancyBboxPatch((0.6, 5.4), 8.8, 0.55, boxstyle="round,pad=0.02,rounding_size=0.05", fc=NAVY, ec=NAVY))
    axR.text(5, 5.67, "COPY NUMBER  =  χ_H", ha="center", va="center", fontsize=14, weight="bold", color="white")
    axR.text(5, 5.22, "chromatic number of the read-conflict graph  (Minimum Copy Cover; Lemma 1)", ha="center", fontsize=8.6, style="italic", color="#666")
    # a copy = one COLOUR class of reads; reads are coloured by the copy they belong to
    for y, c, lab in [(4.35, BLUE, "copy 1 (colour 1)"), (3.2, RED, "copy 2 (colour 2)")]:
        axR.plot([1.3, 8.2], [y, y], color=c, lw=1.3, ls=(0, (2, 2)), alpha=.45, zorder=2)
        axR.text(8.5, y, lab, ha="left", va="center", fontsize=8.6, weight="bold", color=c)
    for x0, x1, y, c in [(1.5, 3.1, 4.35, BLUE), (3.6, 5.2, 4.35, BLUE), (5.9, 8.0, 4.35, BLUE),
                         (1.8, 3.4, 3.2, RED), (4.0, 5.6, 3.2, RED), (6.2, 8.0, 3.2, RED)]:
        axR.plot([x0, x1], [y, y], color=c, lw=6, alpha=.75, solid_capstyle="round", zorder=3)
    axR.text(5.0, 2.35, "a copy = one COLOUR class (reads that never disagree at a PSV)", ha="center", va="center",
             fontsize=9.0, weight="bold", color=DARK)
    axR.text(5.0, 1.5, "copy number = fewest colours = χ_H   (MCC = χ(H))", ha="center", fontsize=12, weight="bold", color=NAVY)
    axR.text(5, 0.55, "greedy colouring gives a LOWER bound  ·  K = 0 identical copies carry no PSV → collapse to one.",
             ha="center", va="center", fontsize=8.4, style="italic", color="#777")

    p = os.path.join(OUT, "formal_objects.png")
    fig.savefig(p, bbox_inches="tight", facecolor="white", dpi=150)
    plt.close()
    return p


if __name__ == "__main__":
    print("wrote", build())
