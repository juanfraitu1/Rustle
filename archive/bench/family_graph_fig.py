#!/usr/bin/env python3
"""Figure: graph-based (POA variation graph) family definition — clean separation + over-merge fix."""
import json
import os
import random

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

BASE = os.path.dirname(__file__)
NAVY = "#22313f"; GREEN = "#1e8449"; ORANGE = "#e8590c"; RED = "#c0392b"; SLATE = "#8898aa"
plt.rcParams["font.family"] = "DejaVu Sans"
T = 0.045

rows = []
for line in open(os.path.join(BASE, "family_graph_scores.tsv")):
    if line.startswith("family_id"):
        continue
    f = line.rstrip("\n").split("\t")
    rows.append((int(f[1]), float(f[2])))
lab = json.load(open(os.path.join(BASE, "family_graph_labeled.json")))
cur = {k: v[0] for k, v in lab["curated"].items()}
dom = {k: v[0] for k, v in lab["domain_sharers"].items()}
big = [s for n, s in rows if n >= 25]
small = [s for n, s in rows if n <= 5]

fig, (axA, axB) = plt.subplots(1, 2, figsize=(13.2, 5.4), gridspec_kw={"width_ratios": [1.1, 1.0]})

rng = random.Random(1)
def strip(ax, yi, vals, color, label):
    ys = [yi + (rng.random() - 0.5) * 0.3 for _ in vals]
    ax.scatter(vals, ys, s=55, color=color, alpha=0.85, edgecolor="white", linewidth=0.5, zorder=3, label=label)

strip(axA, 3, list(cur.values()), GREEN, f"curated real families (n={len(cur)})")
strip(axA, 2, list(dom.values()), RED, f"domain-sharers (n={len(dom)})")
strip(axA, 1, big, ORANGE, f"big components N≥25 (n={len(big)})")
strip(axA, 0, [s for n, s in rows if 2 <= n <= 5][:300], SLATE, "small families N≤5 (sample)")
for nm, sc in cur.items():
    if nm in ("RABL2", "APOBEC3"):
        axA.annotate(nm, (sc, 3), fontsize=8, color=GREEN, xytext=(sc, 3.3))
axA.axvline(T, color="black", ls="--", lw=1.4)
axA.text(T + 0.005, 3.6, f"T = {T}", fontsize=10)
axA.set_yticks([3, 2, 1, 0]); axA.set_yticklabels(["real\nfamilies", "domain-\nsharers", "big\ncomponents", "small\nfamilies"])
axA.set_xlabel("graph core-conservation score  |R| / median|s|")
axA.set_xlim(-0.02, 1.02); axA.set_ylim(-0.6, 4.0)
axA.set_title("Clean separation on the labeled set\n(real ≥0.062 > T > domain-sharers ≤0.029)", fontsize=10.5, color=NAVY)
axA.legend(fontsize=7.5, loc="lower right"); axA.grid(axis="x", alpha=0.2)

axB.axis("off"); axB.set_xlim(0, 10); axB.set_ylim(0, 10)
axB.set_title("Formal graph-based family definition", fontsize=11.5, color=NAVY, fontweight="bold")
axB.text(0.2, 8.9,
         "Members S={s₁..s_N}; POA graph G, columns c₁..c_M.\n"
         "supp(c)=#{i: sᵢ non-gap at c}.\n"
         "CORE R = longest contiguous run with\n"
         "    supp(c) ≥ max(2, ⌈N/2⌉)   (majority spine).\n"
         "FAMILY at level T  ⟺  |R| ≥ T·median|sᵢ|.\n"
         "score(S) = |R| / median|sᵢ|  ∈ [0,1].",
         fontsize=9.2, color=NAVY, va="top", family="DejaVu Sans",
         bbox=dict(boxstyle="round,pad=0.5", facecolor="#f4f6f8", edgecolor=SLATE))
lines = [
    ("reduces to pairwise", "at N=2, ⌈N/2⌉→2 ⇒ both copies ⇒ the validated pairwise core (RABL2/DAZ)", GREEN),
    ("separates", "real ≥0.062 vs domain-sharers ≤0.029 (gap at T=0.045): perfect on labeled set", GREEN),
    ("fixes over-merge", "7/14 mega-components scored <T = NOT single families (pairwise CC couldn't)", ORANGE),
    ("one statistic", "from one graph, no O(N²) pairwise / transitive closure", NAVY),
    ("graded", "recent dups high (0.16), diverged families moderate (0.06–0.08): tightness measure", SLATE),
]
yy = 3.9
for k, v, c in lines:
    axB.text(0.2, yy, f"• {k}:", fontsize=9, color=c, fontweight="bold")
    axB.text(0.5, yy - 0.42, v, fontsize=8, color=NAVY)
    yy -= 0.95

fig.suptitle("Upgrade #1 (option): per-family POA variation graph → a clean formal family definition",
             fontsize=13, fontweight="bold", color=NAVY, y=1.0)
fig.tight_layout(rect=[0, 0, 1, 0.94])
out = os.path.join(BASE, "family_graph_definition.png")
fig.savefig(out, dpi=160, facecolor="white")
print("saved", out)
