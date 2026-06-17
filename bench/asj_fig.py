#!/usr/bin/env python3
"""Figure: allele-specific junctions. Left = funnel (tested -> ASJ -> genetic -> splice-proximal) +
dPSI distribution; right = a mechanistic example (SNP at the splice site flips junction usage per
molecule). Numbers from asj_aggregate.py / asj_verify.py / asj_evidence.py.
"""
import csv
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.patches import FancyBboxPatch

BASE = os.path.dirname(__file__)
NAVY = "#22313f"; ORANGE = "#e8590c"; GREEN = "#1e8449"; SLATE = "#8898aa"; RED = "#c0392b"
plt.rcParams["font.family"] = "DejaVu Sans"

calls = list(csv.DictReader(open(os.path.join(BASE, "asj_calls_verified.tsv")), delimiter="\t"))
dpsi = np.array([float(r["dPSI"]) for r in calls])

fig = plt.figure(figsize=(13.4, 6.0))
axL = fig.add_axes([0.04, 0.10, 0.40, 0.78])
axR = fig.add_axes([0.54, 0.12, 0.43, 0.76])

# ---- LEFT funnel ----
axL.set_xlim(0, 10); axL.set_ylim(0, 10); axL.axis("off")
axL.text(5, 9.6, "Allele-specific junctions (genome-wide)", fontsize=12.5, fontweight="bold",
         color=NAVY, ha="center")
stages = [
    ("74,674", "alternatively-spliced junctions tested\n(7,898 phaseable het genes)", NAVY, 1.0),
    ("475", "allele-specific (BH-FDR q<0.05, |dPSI|≥0.3)\nacross 235 genes", ORANGE, 0.66),
    ("120", "genetic (transversion anchor, not RNA-edit)\n59 genes · 0 collapsed-paralog", GREEN, 0.42),
    ("20", "splice-proximal (SNP ≤100bp from junction)\n= candidate splice-disrupting variants", RED, 0.22),
]
y = 8.2; cx = 5.0; maxw = 8.6
for i, (val, lab, col, frac) in enumerate(stages):
    w = maxw * (0.18 + 0.82 * frac); h = 1.25
    axL.add_patch(FancyBboxPatch((cx - w/2, y - h/2), w, h, boxstyle="round,pad=0.02,rounding_size=0.07",
                                 facecolor=col, edgecolor="none", alpha=0.92))
    axL.text(cx, y + 0.22, val, fontsize=17, fontweight="bold", color="white", ha="center", va="center")
    axL.text(cx, y - 0.36, lab, fontsize=7.6, color="white", ha="center", va="center")
    if i < len(stages) - 1:
        axL.annotate("", xy=(cx, y - h/2 - 0.55), xytext=(cx, y - h/2 - 0.02),
                     arrowprops=dict(arrowstyle="-|>", color=SLATE, lw=1.8))
    y -= 2.0
axL.text(5, 0.2, f"strong effects: median |dPSI|=0.64 · {int((dpsi>=0.7).sum())} ≥0.7 · "
                 f"{int((dpsi>=0.999).sum())} full switches (1.0)",
         fontsize=8.4, color=NAVY, ha="center", style="italic")

# ---- RIGHT mechanistic example: PSMD2 (SNP at splice donor) ----
axR.set_xlim(0, 10); axR.set_ylim(0, 10); axR.axis("off")
axR.set_title("Mechanism: the molecule's allele decides if it splices", fontsize=12, color=NAVY,
              fontweight="bold")
axR.text(5, 8.9, "PSMD2 — SNP at the splice donor (canonical GT-AG)", fontsize=10, color=NAVY, ha="center")
# exon-intron cartoon
ex_y = 7.3
axR.add_patch(FancyBboxPatch((0.8, ex_y-0.3), 3.0, 0.6, boxstyle="square,pad=0", facecolor=SLATE, alpha=0.5, edgecolor="none"))
axR.add_patch(FancyBboxPatch((6.2, ex_y-0.3), 3.0, 0.6, boxstyle="square,pad=0", facecolor=SLATE, alpha=0.5, edgecolor="none"))
axR.plot([3.8, 6.2], [ex_y, ex_y], color=SLATE, lw=1, ls=(0,(2,2)))
axR.text(2.3, ex_y, "exon", fontsize=8, color="white", ha="center", va="center")
axR.text(7.7, ex_y, "exon", fontsize=8, color="white", ha="center", va="center")
axR.plot([3.95, 3.95], [ex_y-0.55, ex_y+0.55], color=NAVY, lw=1.5)
axR.text(3.95, ex_y+0.85, "SNP at donor\n(T vs G)", fontsize=8, color=NAVY, ha="center")
axR.text(4.4, ex_y-0.7, "GT...AG", fontsize=8, color=GREEN, ha="left")
# per-allele PSI bars
for (al, psi, used, span, col, yy) in [("G  (donor present)", 1.0, 14, 14, GREEN, 5.0),
                                       ("T  (donor disrupted)", 0.0, 0, 18, RED, 3.6)]:
    axR.text(0.3, yy, al, fontsize=9.5, color=NAVY, va="center")
    axR.add_patch(FancyBboxPatch((4.0, yy-0.32), 5.0, 0.64, boxstyle="square,pad=0",
                                 facecolor="#eceff1", edgecolor="none"))
    axR.add_patch(FancyBboxPatch((4.0, yy-0.32), 5.0*psi, 0.64, boxstyle="square,pad=0",
                                 facecolor=col, edgecolor="none", alpha=0.9))
    axR.text(9.15, yy, f"PSI={psi:.0f}  ({used}/{span})", fontsize=8.5, color=col, va="center", fontweight="bold")
axR.text(5, 2.4, "dPSI = 1.00,  q = 5.8e-7,  uniquely mapped,  transversion (genetic)",
         fontsize=8.6, color=NAVY, ha="center", fontweight="bold")
axR.text(5, 1.5, "Other textbook hits: DAXX (SNP at acceptor, 0/10 vs 9/9), KCNAB2, CSNK1D.\n"
                 "Long reads observe allele + junctions on the SAME molecule — no statistical phasing.",
         fontsize=8.0, color=SLATE, ha="center", va="center", style="italic")

fig.suptitle("Quantifying allele-specific junctions (het = phasing signal; the task-c confound, reused)",
             fontsize=13.5, fontweight="bold", color=NAVY, y=0.99)
out = os.path.join(BASE, "asj_findings.png")
fig.savefig(out, dpi=160, facecolor="white")
print("saved", out)
