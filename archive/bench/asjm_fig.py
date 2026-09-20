#!/usr/bin/env python3
"""Figure: multi-SNP haplotype phasing gain over single-anchor for allele-specific junctions."""
import csv
import os
from collections import defaultdict

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Circle

BASE = os.path.dirname(__file__)
NAVY = "#22313f"; ORANGE = "#e8590c"; GREEN = "#1e8449"; SLATE = "#8898aa"
plt.rcParams["font.family"] = "DejaVu Sans"

calls = list(csv.DictReader(open(os.path.join(BASE, "asjm_calls.tsv")), delimiter="\t"))
nd = defaultdict(int)
for r in calls:
    nd[min(int(r["n_snp"]), 10)] += 1

# numbers from asjm_aggregate
SINGLE, MULTI, SHARED, MO, SO = 475, 630, 448, 182, 27

fig, (axA, axB) = plt.subplots(1, 2, figsize=(13.2, 5.2), gridspec_kw={"width_ratios": [1.0, 1.1]})

# Panel A: Venn single vs multi
axA.set_xlim(0, 10); axA.set_ylim(0, 10); axA.axis("off")
axA.set_title("Multi-SNP phasing vs single-anchor\n(allele-specific junctions)", fontsize=12,
              color=NAVY, fontweight="bold")
axA.add_patch(Circle((4.0, 5.0), 2.7, facecolor=SLATE, alpha=0.18, edgecolor=SLATE, lw=2))
axA.add_patch(Circle((6.2, 5.0), 3.0, facecolor=GREEN, alpha=0.16, edgecolor=GREEN, lw=2))
axA.text(2.6, 7.2, "single-anchor", fontsize=10, color=SLATE, fontweight="bold", ha="center")
axA.text(2.6, 6.6, f"{SINGLE}", fontsize=12, color=SLATE, fontweight="bold", ha="center")
axA.text(7.4, 7.7, "multi-SNP", fontsize=10.5, color=GREEN, fontweight="bold", ha="center")
axA.text(7.4, 7.1, f"{MULTI}", fontsize=13, color=GREEN, fontweight="bold", ha="center")
axA.text(2.7, 5.0, f"{SO}\nsingle-only", fontsize=8.5, color=SLATE, ha="center", va="center")
axA.text(4.9, 5.0, f"{SHARED}\nshared", fontsize=9.5, color=NAVY, ha="center", va="center", fontweight="bold")
axA.text(7.6, 5.0, f"{MO}\nmulti-only\n(reach gain)", fontsize=9, color=GREEN, ha="center", va="center", fontweight="bold")
axA.text(5, 1.0, f"+{MULTI-SINGLE} ASJ (+{round(100*(MULTI-SINGLE)/SINGLE)}%) by phasing reads the single anchor couldn't reach",
         fontsize=9, color=GREEN, ha="center", style="italic", fontweight="bold")

# Panel B: ASJ by n_snp
axB.set_title("ASJ by # het SNPs used to phase the gene", fontsize=12, color=NAVY, fontweight="bold")
ks = sorted(nd)
vals = [nd[k] for k in ks]
labels = [f"{k}" if k < 10 else "10+" for k in ks]
cols = [SLATE if k == 1 else GREEN for k in ks]
axB.bar(range(len(ks)), vals, color=cols, alpha=0.9)
axB.set_xticks(range(len(ks))); axB.set_xticklabels(labels)
axB.set_xlabel("# het SNPs (n_snp);  n_snp=1 ≈ single-anchor (grey)")
axB.set_ylabel("allele-specific junctions")
for i, v in enumerate(vals):
    axB.text(i, v + 3, str(v), ha="center", fontsize=8, color=NAVY)
axB.text(0.97, 0.92, f"{sum(v for k,v in nd.items() if k>=2)} ASJ from\nmulti-SNP genes (>=2)",
         transform=axB.transAxes, ha="right", va="top", fontsize=9, color=GREEN, fontweight="bold")
for s in axB.spines.values():
    s.set_color("#cccccc")

fig.suptitle("Multi-SNP haplotype phasing for allele-specific junctions",
             fontsize=13.5, fontweight="bold", color=NAVY, y=1.0)
fig.tight_layout(rect=[0, 0, 1, 0.93])
out = os.path.join(BASE, "asjm_findings.png")
fig.savefig(out, dpi=160, facecolor="white")
print("saved", out)
