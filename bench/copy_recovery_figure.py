#!/usr/bin/env python3
"""Explanatory figure: how the paralog-family splice graph looks with PRIMARY reads only
vs. with SECONDARY (multimapper) alignments added — illustrating why primary-only assembly
(StringTie / long-read) starves a paralog copy, and how the secondary alignments + PSVs
resolve it. Schematic, annotated with the real GGO result (XM_055380753.2)."""
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, Rectangle
from matplotlib.lines import Line2D
import numpy as np

PRIM = "#222222"     # primary alignments (black)
SEC  = "#e8730c"     # secondary / multimapper alignments (orange)
EXON = "#cfe3f7"     # exon fill
EXON_B_off = "#eeeeee"
PSV_A = "#1a7f37"    # PSV allele for copy A (green)
PSV_B = "#e8730c"    # PSV allele for copy B (orange)

# exon x-spans (schematic genomic coords), shared structure between the two copies
EXONS = [(4, 14), (26, 36), (50, 60), (74, 90)]
PSV_X = [31, 55]     # two PSV positions (paralog-distinguishing bases)

def draw_copy(ax, y, label, exon_color, faded=False, psv_color=None, edge=PRIM):
    a = 0.35 if faded else 1.0
    # exon blocks
    for (x0, x1) in EXONS:
        ax.add_patch(Rectangle((x0, y-0.6), x1-x0, 1.2, facecolor=exon_color,
                               edgecolor=edge, lw=1.6, alpha=a, zorder=3))
    # introns (junction arcs) connecting consecutive exons
    for (x0, x1), (x2, x3) in zip(EXONS[:-1], EXONS[1:]):
        xm = (x1 + x2) / 2
        ax.plot([x1, xm, x2], [y+0.55, y+1.15, y+0.55], color=edge, lw=1.4, alpha=a, zorder=2)
    # PSV ticks
    if psv_color:
        for px in PSV_X:
            ax.plot([px, px], [y-0.6, y+0.6], color=psv_color, lw=2.4, zorder=5)
            ax.scatter([px], [y+0.6], marker="v", s=36, color=psv_color, zorder=6)
    ax.text(0.5, y, label, ha="right", va="center", fontsize=10, fontweight="bold")

def draw_reads(ax, y, spans, color, n_label=None, alpha=0.9):
    for i, (x0, x1, yo) in enumerate(spans):
        ax.plot([x0, x1], [y-yo, y-yo], color=color, lw=2.2, solid_capstyle="round",
                alpha=alpha, zorder=4)
    if n_label:
        ax.text(106, y-0.4, n_label, color=color, fontsize=8.0, va="center", ha="right")

fig, (axL, axR) = plt.subplots(1, 2, figsize=(13.5, 5.6))
fig.suptitle("Multimapper-resolved paralog-copy recovery: the splice graph with vs. without secondary alignments",
             fontsize=12.5, fontweight="bold", y=0.98)

# read spans (x0, x1, y-offset) — schematic reads spanning exon pairs
reads = [(5, 35, 1.6), (27, 59, 2.0), (51, 89, 2.4), (5, 59, 2.8),
         (27, 89, 3.2), (5, 35, 3.6), (51, 89, 4.0)]

for ax in (axL, axR):
    ax.set_xlim(-17, 108); ax.set_ylim(-5.6, 6.2); ax.axis("off")

# ---- LEFT: primary alignments only ----
axL.set_title("A.  Primary alignments only  (StringTie / long-read)", fontsize=11, loc="left")
draw_copy(axL, 4.4, "Copy A\n(expressed)", EXON, edge=PRIM, psv_color=PSV_A)
draw_copy(axL, 0.2, "Copy B\n(paralog)", EXON_B_off, faded=True, edge="#999999")
# all multimapper reads land their PRIMARY on Copy A (arbitrary MAPQ-0 tie-break)
draw_reads(axL, 4.0, reads, PRIM, n_label="primary reads → A")
axL.text(47, -2.4, "Copy B has ~no primary reads\n→ starved → NOT assembled\n(StringTie misses XM_055380753.2)",
         ha="center", va="center", fontsize=9.5, color="#888888",
         bbox=dict(boxstyle="round", fc="#f6f6f6", ec="#cccccc"))

# ---- RIGHT: + secondary alignments ----
axR.set_title("B.  + Secondary (multimapper) alignments", fontsize=11, loc="left")
draw_copy(axR, 4.4, "Copy A\n(expressed)", EXON, edge=PRIM, psv_color=PSV_A)
draw_copy(axR, 0.2, "Copy B\n(recovered)", "#fbe3cc", edge=SEC, psv_color=PSV_B)
draw_reads(axR, 4.0, reads, PRIM, n_label="primary (same reads)")
# the SAME molecules' secondary alignments reveal copy B; PSV allele assigns them there
secreads = [(5, 35, 1.6), (27, 59, 2.0), (51, 89, 2.4), (27, 89, 2.8), (5, 59, 3.2)]
draw_reads(axR, 0.0, [(x0, x1, -yo) for (x0, x1, yo) in secreads], SEC,
           n_label="secondary → Copy B")
# dotted "same molecule" links between a primary on A and its secondary on B
for (x0, x1, yo) in secreads[:3]:
    axR.plot([(x0+x1)/2, (x0+x1)/2], [4.0-yo, 0.0+yo], ls=":", color=SEC, lw=1.0, alpha=0.7, zorder=1)
axR.text(47, -2.7, "PSV alleles (▼) assign each multimapper to its copy of origin\n"
                   "→ Copy B's graph is built from its rightful reads\n"
                   "→ recovered: XM_055380753.2 (class c) + 21 more isoforms,\n"
                   "   0 chain-overlap with StringTie",
         ha="center", va="center", fontsize=9.5, color="#8a4500",
         bbox=dict(boxstyle="round", fc="#fff6ec", ec=SEC))

# legend
leg = [Line2D([0],[0], color=PRIM, lw=2.6, label="primary alignment"),
       Line2D([0],[0], color=SEC, lw=2.6, label="secondary (multimapper) alignment"),
       Line2D([0],[0], marker="v", color="w", markerfacecolor=PSV_A, markersize=9, label="PSV (paralog-sequence variant)")]
fig.legend(handles=leg, loc="lower center", ncol=3, frameon=False, fontsize=9.5, bbox_to_anchor=(0.5, 0.005))

fig.text(0.5, 0.075, "Same reads, same graph topology — the secondary alignments + PSVs let the family graph "
                     "partition the reads into the correct copy, recovering a copy primary-only assembly starves.",
         ha="center", fontsize=9, style="italic", color="#444444")
fig.subplots_adjust(top=0.90, bottom=0.16, left=0.04, right=0.985, wspace=0.08)
out = "/mnt/c/Users/jfris/Desktop/Rustle/bench/copy_recovery_figure.png"
fig.savefig(out, dpi=160)
print(f"wrote {out}")
