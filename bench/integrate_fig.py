#!/usr/bin/env python3
"""Figure: integrated family-detection -> read-to-copy assignment, hard-multimapper focus."""
import json
import os
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

BASE = os.path.dirname(__file__)
NAVY = "#22313f"; GREEN = "#1e8449"; ORANGE = "#e8590c"; RED = "#c0392b"; SLATE = "#8898aa"
plt.rcParams["font.family"] = "DejaVu Sans"
d = json.load(open("/home/juanfra/winloci_scratch/sim5x/integrate_summary.json"))
S = d["summary"]   # (K, fscore, detected, npsv, hard, acc_hard, acc_easy)

fig, (axA, axB) = plt.subplots(1, 2, figsize=(13.2, 5.3), gridspec_kw={"width_ratios": [1.1, 1.0]})

# Panel A: synthetic end-to-end (truth) — hard-MM assignment by K
K = [s[0] for s in S]
hard = [s[4] for s in S]
acc = [(s[5] if s[5] is not None else 0) for s in S]
x = range(len(K))
ax2 = axA.twinx()
b = axA.bar([i - 0.2 for i in x], hard, width=0.4, color=SLATE, label="hard multimappers (MAPQ-0)")
ax2.plot(x, [a * 100 for a in acc], "-o", color=GREEN, lw=2.4, label="PSV assignment acc (hard MM)")
axA.set_xticks(list(x)); axA.set_xticklabels(K)
axA.set_xlabel("PSVs per copy (K) — synthetic 5-copy, WITH ground truth")
axA.set_ylabel("hard multimappers (MAPQ-0)", color=SLATE)
ax2.set_ylabel("PSV assignment accuracy on hard MM (%)", color=GREEN)
ax2.set_ylim(-5, 108)
axA.set_title("End-to-end (truth): family detected at every K (✓);\n"
              "PSV resolves hard MM up to identifiability", fontsize=10.5, color=NAVY)
axA.annotate("K=0 identical:\n200 MAPQ-0, UNASSIGNABLE\n(no PSV info — impossible)", (0, 200),
             fontsize=7.5, color=RED, xytext=(0.5, 150),
             arrowprops=dict(arrowstyle="->", color=RED, lw=1))
axA.annotate("K≥2: minimap2 itself\nresolves (0 hard MM)\nfor full-length reads", (3, 5),
             fontsize=7.5, color=NAVY, xytext=(2.3, 120),
             arrowprops=dict(arrowstyle="->", color=SLATE, lw=1))
axA.legend(fontsize=7.5, loc="center left"); ax2.legend(fontsize=7.5, loc="center right")

# Panel B: real GGO census + bottom line
axB.axis("off"); axB.set_xlim(0, 10); axB.set_ylim(0, 10)
axB.set_title("Real GGO: hard multimappers are RARE", fontsize=11.5, color=NAVY, fontweight="bold")
axB.text(0.3, 8.7,
         "of 174,459 reads at multi-copy family loci:\n"
         "  hard multimappers (MAPQ-0) = 1,156  (0.7%)\n"
         "  families with ≥5 hard MM = 26 / 400  (6%)\n"
         "  (a few families carry most: FAM156 = 358)",
         fontsize=9.2, color=NAVY, va="top",
         bbox=dict(boxstyle="round,pad=0.5", facecolor="#f4f6f8", edgecolor=SLATE))
axB.text(0.3, 5.2, "Bottom line", fontsize=10, color=NAVY, fontweight="bold")
for i, (t, c) in enumerate([
    ("Family identification: YES — graph def, validated, fixes over-merge.", GREEN),
    ("Read→copy assignment: YES, validated WITH TRUTH (synthetic),", GREEN),
    ("   correct up to the identifiability boundary; declares the unassignable.", GREEN),
    ("Hard MM: PSV's edge over minimap2 is the MARGINAL (single-PSV) regime;", ORANGE),
    ("   full-length reads let the aligner resolve once ≥2 PSVs are spanned.", ORANGE),
    ("Real GGO under-exercises the hard regime (0.7%) — needs deep", RED),
    ("   co-located data (testis HiFi DAZ/RBMY) for a rich showcase.", RED),
]):
    axB.text(0.3, 4.5 - i * 0.62, t, fontsize=8.2, color=c)

fig.suptitle("Integrated: identify families → assign reads to copies (hard multimappers)",
             fontsize=13, fontweight="bold", color=NAVY, y=1.0)
fig.tight_layout(rect=[0, 0, 1, 0.94])
out = os.path.join(BASE, "integrate_end_to_end.png")
fig.savefig(out, dpi=160, facecolor="white")
print("saved", out)
