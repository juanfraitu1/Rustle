#!/usr/bin/env python3
"""Figure: annotation-free two-pass (read-coherence -> family + copy assignment) vs flow."""
import json
import os
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch

BASE = os.path.dirname(__file__)
NAVY = "#22313f"; GREEN = "#1e8449"; ORANGE = "#e8590c"; RED = "#c0392b"; SLATE = "#8898aa"
plt.rcParams["font.family"] = "DejaVu Sans"
rows = json.load(open("/home/juanfra/winloci_scratch/twopass/twopass_summary.json"))["rows"]

fig, (axA, axB) = plt.subplots(1, 2, figsize=(13.4, 5.4), gridspec_kw={"width_ratios": [1.15, 1.0]})

# Panel A: architecture (read-coherence vs flow)
axA.set_xlim(0, 10); axA.set_ylim(0, 10); axA.axis("off")
axA.set_title("Annotation-free two-pass (collapsed copies)", fontsize=11.5, color=NAVY, fontweight="bold")
def box(x, y, w, h, txt, col, fs=8.5):
    axA.add_patch(FancyBboxPatch((x, y), w, h, boxstyle="round,pad=0.03,rounding_size=0.1",
                                 facecolor=col, alpha=0.16, edgecolor=col, lw=1.6))
    axA.text(x + w/2, y + h/2, txt, ha="center", va="center", fontsize=fs, color=NAVY)
box(3.2, 8.6, 3.6, 1.0, "reads from 5 copies\n(collapsed onto 1 locus)", SLATE)
axA.add_patch(FancyArrowPatch((5, 8.6), (5, 7.7), arrowstyle="-|>", color=SLATE, lw=1.6, mutation_scale=12))
box(0.6, 6.4, 4.0, 1.2, "PASS 1: read-coherence\ngroup by intron chain →\nONE skeleton (+ its reads)", GREEN)
box(5.4, 6.4, 4.0, 1.2, "FLOW (StringTie):\n→ ~1 transcript\nCOPIES LOST", RED, fs=8.2)
axA.add_patch(FancyArrowPatch((2.6, 6.4), (2.6, 5.5), arrowstyle="-|>", color=GREEN, lw=1.6, mutation_scale=12))
axA.text(7.4, 6.2, "✗ per-read PSVs\ncollapsed away", fontsize=7.5, color=RED, ha="center", va="top")
box(0.4, 3.9, 4.4, 1.4, "PASS 2 (no annotation):\nde-novo PSVs from the\ncollapsed pileup → split\nskeleton into N copies", ORANGE, fs=8.2)
axA.add_patch(FancyArrowPatch((2.6, 3.9), (2.6, 3.1), arrowstyle="-|>", color=ORANGE, lw=1.6, mutation_scale=12))
box(0.8, 1.7, 3.6, 1.2, "multi-copy FAMILY\n+ reads assigned\nto copies", NAVY)
axA.text(5.2, 2.3, "read-coherence PRESERVES per-read\nPSVs → Pass 2 recovers the copies\nflow would have lost",
         fontsize=8.3, color=GREEN, va="center", style="italic")

# Panel B: K-ladder results
K = [r[0] for r in rows]; rec = [r[3] for r in rows]; acc = [r[4]*100 for r in rows]; psv = [r[2] for r in rows]
x = range(len(K))
axB.bar([i-0.2 for i in x], rec, width=0.4, color=ORANGE, label="copies recovered (of 5)")
ax2 = axB.twinx()
ax2.plot(x, acc, "-o", color=GREEN, lw=2.4, label="read→copy assignment acc")
axB.set_xticks(list(x)); axB.set_xticklabels(K)
axB.set_xlabel("PSVs per copy (K)"); axB.set_ylabel("copies recovered (of 5)", color=ORANGE)
ax2.set_ylabel("assignment accuracy (%)", color=GREEN); ax2.set_ylim(-5, 108); axB.set_ylim(0, 5.5)
axB.set_title("Recovers copies + assigns reads, by identifiability\n(Pass 1 = 1 skeleton at every K)", fontsize=10.5, color=NAVY)
for i, p in enumerate(psv):
    axB.text(i-0.2, rec[i]+0.1, f"{p}psv", ha="center", fontsize=7, color=SLATE)
axB.annotate("K=0 identical:\nunassignable", (0, 0.1), fontsize=7.5, color=RED, xytext=(0.4, 2.5),
             arrowprops=dict(arrowstyle="->", color=RED, lw=1))
axB.legend(fontsize=7.5, loc="center left"); ax2.legend(fontsize=7.5, loc="lower right")

fig.suptitle("Read-coherence as Pass 1 preserves per-read PSVs → recovers copies flow loses (annotation-free)",
             fontsize=12.5, fontweight="bold", color=NAVY, y=1.0)
fig.tight_layout(rect=[0, 0, 1, 0.94])
out = os.path.join(BASE, "twopass_denovo.png")
fig.savefig(out, dpi=160, facecolor="white")
print("saved", out)
