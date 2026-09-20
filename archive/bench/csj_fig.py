#!/usr/bin/env python3
"""Figure: copy-specific junctions = the unification (find copies -> assign reads to copies ->
compare junction usage between copies). Flagship: DAZ1 vs DAZL differential splicing."""
import csv
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch

BASE = os.path.dirname(__file__)
NAVY = "#22313f"; ORANGE = "#e8590c"; GREEN = "#1e8449"; SLATE = "#8898aa"; RED = "#c0392b"
plt.rcParams["font.family"] = "DejaVu Sans"

fig, (axL, axR) = plt.subplots(1, 2, figsize=(13.4, 5.6), gridspec_kw={"width_ratios": [1.05, 1.0]})

# LEFT: unification flow + counts
axL.set_xlim(0, 10); axL.set_ylim(0, 10); axL.axis("off")
axL.set_title("Copy-specific junctions = the unification", fontsize=12.5, color=NAVY, fontweight="bold")
steps = [
    ("1. find copies", "cross-chromosome recent-dup pairs\n(248 with BOTH copies expressed)", GREEN),
    ("2. assign reads to copies", "by genomic locus (cross-chrom);\nPSV copy_split in the collapsed regime", ORANGE),
    ("3. compare junctions", "align structures → per-copy PSI of\n696 homologous junctions (Fisher)", NAVY),
]
y = 8.4
for t, d, col in steps:
    axL.add_patch(FancyBboxPatch((0.6, y-0.72), 8.8, 1.44, boxstyle="round,pad=0.02,rounding_size=0.1",
                                 facecolor=col, alpha=0.13, edgecolor=col, lw=1.5))
    axL.text(0.95, y+0.28, t, fontsize=10.5, color=col, fontweight="bold", va="center")
    axL.text(0.95, y-0.32, d, fontsize=8.3, color=NAVY, va="center")
    if y > 3:
        axL.annotate("", xy=(5, y-0.74-0.42), xytext=(5, y-0.74-0.02),
                     arrowprops=dict(arrowstyle="-|>", color=SLATE, lw=1.6))
    y -= 2.15
axL.text(5, 1.55, "146 DIFFERENTIAL copy-specific junctions / 66 paralog pairs",
         fontsize=10, color=RED, ha="center", fontweight="bold")
axL.text(5, 0.95, "(both copies splice but differ; retrocopy/unspliced cases separated)",
         fontsize=8, color=SLATE, ha="center", style="italic")
axL.text(5, 0.35, "+ 1,408 copy-private exon junctions (a copy gained/lost an exon)",
         fontsize=8.4, color=NAVY, ha="center")

# RIGHT: flagship DAZ1 vs DAZL + two more
axR.set_xlim(0, 10); axR.set_ylim(0, 10); axR.axis("off")
axR.set_title("Paralog copies splice differently — e.g. DAZ1 vs DAZL", fontsize=12, color=NAVY, fontweight="bold")
examples = [
    ("DAZ1",  0.074, "DAZL",        0.992, "q=3e-151"),
    ("HERC2", 1.000, "LOC101149126", 0.000, "q=2e-3"),
    ("FRG1",  0.955, "LOC115933219", 0.000, "q=2e-32"),
    ("SORL1", 0.995, "LOC115932779", 0.000, "q=1e-29"),
]
y = 8.2
for a, pa, b, pb, q in examples:
    axR.text(0.3, y+0.32, f"{a}  vs  {b}", fontsize=9.5, color=NAVY, fontweight="bold")
    axR.text(9.7, y+0.32, q, fontsize=8, color=SLATE, ha="right")
    for (name, p, yy, col) in [(a, pa, y, GREEN), (b, pb, y-0.55, ORANGE)]:
        axR.add_patch(FancyBboxPatch((2.6, yy-0.2), 5.5, 0.4, boxstyle="square,pad=0",
                                     facecolor="#eceff1", edgecolor="none"))
        axR.add_patch(FancyBboxPatch((2.6, yy-0.2), 5.5*p, 0.4, boxstyle="square,pad=0",
                                     facecolor=col, edgecolor="none", alpha=0.9))
        axR.text(2.45, yy, name, fontsize=7.5, color=NAVY, ha="right", va="center")
        axR.text(8.25, yy, f"PSI={p:.2f}", fontsize=7.8, color=col, va="center", fontweight="bold")
    y -= 1.9
axR.text(5, 0.4, "homologous junction used by one copy, ~not by the other",
         fontsize=8.5, color=SLATE, ha="center", style="italic")

fig.suptitle("Copy-specific junctions: do paralog copies splice differently?",
             fontsize=13.5, fontweight="bold", color=NAVY, y=1.0)
fig.tight_layout(rect=[0, 0, 1, 0.93])
out = os.path.join(BASE, "copy_specific_junctions.png")
fig.savefig(out, dpi=160, facecolor="white")
print("saved", out)
