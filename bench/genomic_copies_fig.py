#!/usr/bin/env python3
"""Figure: the COMPLETE copy roster across three tiers (annotated coding / pseudogene / unannotated)
with the RNA 'what transcribes' overlay on each."""
import os
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

BASE = os.path.dirname(__file__)
NAVY = "#22313f"; GREEN = "#1e8449"; SLATE = "#b0bec5"; ORANGE = "#e8590c"
plt.rcParams["font.family"] = "DejaVu Sans"

# from dna_rna_overlay (coding) + genomic self-alignment (new copies)
CODING_TOT, CODING_TX = 14545, 10490
PS_TOT, PS_TX = 1313, 239
UN_TOT, UN_TX = 665, 104

fig, (axA, axB) = plt.subplots(1, 2, figsize=(13.2, 5.4), gridspec_kw={"width_ratios": [1.15, 1.0]})

labels = ["annotated\ncoding family\n(protein cluster)", "pseudogene\ncopy\n(no CDS)", "unannotated\ncopy\n(no annotation)"]
tot = [CODING_TOT, PS_TOT, UN_TOT]; tx = [CODING_TX, PS_TX, UN_TX]
sil = [t - x for t, x in zip(tot, tx)]
x = np.arange(3)
axA.bar(x, tx, color=GREEN, label="transcribed (≥5 reads)")
axA.bar(x, sil, bottom=tx, color=SLATE, label="silent")
for i in range(3):
    axA.text(i, tot[i] + 200, f"{tot[i]:,}", ha="center", fontsize=10, fontweight="bold", color=NAVY)
    axA.text(i, tx[i] / 2, f"{100*tx[i]/tot[i]:.0f}%\ntx", ha="center", va="center", fontsize=8.5,
             color="white", fontweight="bold")
axA.set_xticks(x); axA.set_xticklabels(labels, fontsize=8.5)
axA.set_ylabel("gene copies"); axA.legend(fontsize=8.5, loc="upper right")
axA.set_title("Complete copy roster: 3 tiers, with 'what transcribes' overlay", fontsize=11, color=NAVY)
axA.text(1.0, CODING_TOT*0.62,
         "genomic self-alignment adds\n+1,978 copies the annotation &\nprotein clustering missed\n"
         "(mostly silent: ~17% transcribed)",
         fontsize=8, color=ORANGE, ha="center",
         bbox=dict(boxstyle="round,pad=0.4", facecolor="#fff4ec", edgecolor=ORANGE))
for s in axA.spines.values():
    s.set_color("#cccccc")

axB.axis("off"); axB.set_xlim(0, 10); axB.set_ylim(0, 10)
axB.set_title("The complete two-tier picture", fontsize=11.5, color=NAVY, fontweight="bold")
lines = [
    ("DNA roster", "annotated coding (protein families) + pseudogene + unannotated copies", NAVY),
    ("coding-family copies", f"{CODING_TOT:,} ({100*CODING_TX/CODING_TOT:.0f}% transcribed)", GREEN),
    ("pseudogene copies", f"{PS_TOT:,} (18% transcribed) — protein-clustering blind spot", ORANGE),
    ("unannotated copies", f"{UN_TOT:,} (16% transcribed) — annotation blind spot", ORANGE),
    ("RNA layer", "shows which copies are live (transcribed) per copy", GREEN),
    ("gradient", "coding 72% >> pseudogene/unannotated ~17% transcribed", NAVY),
    ("example", "SDHA retrocopy: unannotated/pseudo, 335 reads, 94.8% id", NAVY),
]
yy = 8.4
for k, v, c in lines:
    axB.text(0.3, yy, f"{k}:", fontsize=9.5, color=NAVY, fontweight="bold")
    axB.text(3.5, yy, v, fontsize=8.6, color=c)
    yy -= 1.05
axB.text(5, 0.4,
         "DNA enumerates every copy (incl. pseudogene + unannotated);\nRNA shows what is actually transcribed",
         fontsize=8.8, color=NAVY, ha="center", style="italic")

fig.suptitle("Genomic self-alignment completes the copy roster + RNA overlay",
             fontsize=13, fontweight="bold", color=NAVY, y=1.0)
fig.tight_layout(rect=[0, 0, 1, 0.94])
out = os.path.join(BASE, "genomic_copies.png")
fig.savefig(out, dpi=160, facecolor="white")
print("saved", out)
