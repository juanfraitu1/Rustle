#!/usr/bin/env python3
"""Figure for the minimizer-free intron-chain (structural) copy discovery.
Left: structure vs sequence concordance (independent confirmation). Right: candidate-gen tradeoff.
"""
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Circle

NAVY = "#22313f"; ORANGE = "#e8590c"; GREEN = "#1e8449"; SLATE = "#8898aa"
plt.rcParams["font.family"] = "DejaVu Sans"

fig, (axA, axB) = plt.subplots(1, 2, figsize=(13.2, 5.2), gridspec_kw={"width_ratios": [1.0, 1.15]})

# Panel A: two overlapping circles — structure vs sequence (cross-chrom copies)
axA.set_xlim(0, 10); axA.set_ylim(0, 10); axA.axis("off")
axA.set_title("Independent confirmation: structure vs sequence\n(cross-chromosome copies)",
              fontsize=11.5, color=NAVY, fontweight="bold")
c_seq = Circle((6.4, 5.2), 3.2, facecolor=NAVY, alpha=0.14, edgecolor=NAVY, lw=2)
c_str = Circle((4.1, 5.2), 1.7, facecolor=GREEN, alpha=0.18, edgecolor=GREEN, lw=2)
axA.add_patch(c_seq); axA.add_patch(c_str)
axA.text(7.6, 6.7, "sequence (POA)", fontsize=10, color=NAVY, fontweight="bold", ha="center")
axA.text(7.7, 5.0, "7,964", fontsize=15, color=NAVY, fontweight="bold", ha="center")
axA.text(7.7, 4.4, "sequence-only\n(retrocopies, single-exon)", fontsize=7.6, color=SLATE, ha="center")
axA.text(3.0, 7.4, "intron-chain\n(structure)", fontsize=10, color=GREEN, fontweight="bold", ha="center")
axA.text(3.0, 5.2, "20", fontsize=13, color=GREEN, fontweight="bold", ha="center")
axA.text(3.0, 4.7, "structure-\nonly", fontsize=7.6, color=SLATE, ha="center")
axA.text(5.05, 5.2, "340", fontsize=15, color=ORANGE, fontweight="bold", ha="center")
axA.text(5.05, 4.6, "BOTH", fontsize=8.5, color=ORANGE, fontweight="bold", ha="center")
axA.text(5.0, 1.0, "94% of structural cross-chrom copies are also found by sequence\n"
                   "→ minimizer-free structure is a strong independent confirmation axis",
         fontsize=8.8, color=NAVY, ha="center", style="italic")

# Panel B: candidate-gen tradeoff (cross-chrom confirmed, structure-only, recall)
axB.axis("off")
axB.set_title("Candidate-generation tradeoff (no sequence minimizers)", fontsize=11.5,
              color=NAVY, fontweight="bold")
rows = [
    ("config", "cross-chrom", "struct-only", "univ recall", "runtime"),
    ("length-window, matched≥4  (default)", "360", "20", "1/5", "~5 min"),
    ("length-window, matched≥3", "993", "527", "1/5", "~5 min"),
    ("2-intron-shingle, matched≥3", "38,588", "37,138", "3/5", "~50 min"),
]
y = 8.6
for ri, r in enumerate(rows):
    fw = "bold" if ri == 0 else "normal"
    col = NAVY if ri == 0 else (ORANGE if ri == 1 else SLATE)
    axB.text(0.2, y, r[0], fontsize=8.6, color=col, fontweight=fw, family="DejaVu Sans")
    for ci, x in enumerate((5.9, 7.2, 8.5, 9.5)):
        axB.text(x, y, r[ci + 1], fontsize=8.6, color=col, fontweight=fw, ha="center")
    y -= 0.85
    if ri == 0:
        axB.plot([0.15, 10], [y + 0.35, y + 0.35], color=SLATE, lw=0.6)
axB.set_xlim(0, 10.5); axB.set_ylim(0, 9.5)
axB.text(0.2, y - 0.1,
         "• single intron length not discriminative (common sizes) → 2-intron shingles needed\n"
         "• exon length alone → 173k chance matches; coupling exon+intron fixes precision\n"
         "• length-window precision ↔ misses partial copies (inherent tradeoff)\n"
         "• RABL2A/B recovered (matched 6/8, indel-robust); domain-sharers rejected by construction\n"
         "• blind spot: retrocopies / single-exon (no intron chain) — sequence axis covers these",
         fontsize=8.0, color=NAVY, va="top")

fig.suptitle("Minimizer-free copy discovery: intron-chain alignment (structural axis)",
             fontsize=13.5, fontweight="bold", color=NAVY, y=1.0)
fig.tight_layout(rect=[0, 0, 1, 0.94])
out = "/mnt/c/Users/jfris/Desktop/Rustle/bench/intronchain_discovery.png"
fig.savefig(out, dpi=160)
print("saved", out)
