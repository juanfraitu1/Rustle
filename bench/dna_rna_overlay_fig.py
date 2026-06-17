#!/usr/bin/env python3
"""Figure: two-tier overlay — DNA/protein families + RNA 'what transcribes' layer."""
import os
from collections import defaultdict
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

BASE = os.path.dirname(__file__)
NAVY = "#22313f"; GREEN = "#1e8449"; ORANGE = "#e8590c"; SLATE = "#8898aa"; RED = "#c0392b"
plt.rcParams["font.family"] = "DejaVu Sans"

# read overlay tsv
copies = tx = well = 0
for line in open(os.path.join(BASE, "dna_families_overlay.tsv")):
    if line.startswith("dna_family\t"):
        continue
    f = line.rstrip("\n").split("\t")
    copies += int(f[1]); tx += int(f[2]); well += int(f[3])
silent = copies - tx
mid = tx - well

fig, (axA, axB) = plt.subplots(1, 2, figsize=(13.0, 5.2), gridspec_kw={"width_ratios": [1.0, 1.1]})

# Panel A: DNA vs RNA family count + complementarity
axA.bar(["DNA / protein\n(mmseqs)", "RNA / POA\n(transcript)"], [3587, 1337],
        color=[NAVY, GREEN], alpha=0.9, width=0.6)
axA.set_ylabel("multi-copy families (size ≥ 2)")
axA.set_title("DNA tier finds 2.7× more families (catches ancient ones)", fontsize=11, color=NAVY)
for i, v in enumerate([3587, 1337]):
    axA.text(i, v + 50, f"{v:,}", ha="center", fontsize=11, fontweight="bold", color=NAVY)
axA.text(0.5, 2600,
         "complementary, both ways:\n• DNA recovers RNA-missed ANCIENT families\n   (DEFB, SIGLEC…)\n"
         "• RNA recovers DNA-split families\n   (DAZ1/DAZL — repeat structure)\n• both over-merge superfamilies",
         fontsize=8.2, color=SLATE, ha="center", va="top",
         bbox=dict(boxstyle="round,pad=0.4", facecolor="#f4f6f8", edgecolor=SLATE))
for s in axA.spines.values():
    s.set_color("#cccccc")

# Panel B: genome vs transcriptome (what actually transcribes)
axB.axis("equal")
wedges, _ = axB.pie([well, mid, silent], colors=[GREEN, ORANGE, RED], startangle=90,
                    wedgeprops=dict(width=0.42, edgecolor="white"))
axB.text(0, 0.12, f"{copies:,}", ha="center", fontsize=15, fontweight="bold", color=NAVY)
axB.text(0, -0.18, "DNA copies", ha="center", fontsize=9, color=SLATE)
axB.set_title("What actually transcribes\n(per DNA-defined copy, real GGO IsoSeq)", fontsize=11, color=NAVY)
axB.legend(wedges,
           [f"well-expressed ≥40   {well:,} ({100*well/copies:.0f}%)",
            f"transcribed 5–40     {mid:,} ({100*mid/copies:.0f}%)",
            f"SILENT (<5)          {silent:,} ({100*silent/copies:.0f}%)"],
           loc="lower center", bbox_to_anchor=(0.5, -0.28), fontsize=9, frameon=False)

fig.suptitle("Two-tier overlay: DNA defines the family roster, RNA shows what's live\n"
             "28% of paralog copies are present in the genome but transcriptionally silent",
             fontsize=12.5, fontweight="bold", color=NAVY, y=1.04)
fig.tight_layout(rect=[0, 0, 1, 0.9])
out = os.path.join(BASE, "dna_rna_overlay.png")
fig.savefig(out, dpi=160, facecolor="white", bbox_inches="tight")
print("saved", out)
