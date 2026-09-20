#!/usr/bin/env python3
"""Two-panel figure for the DNA-truth precision/recall validation:
   A. recall decomposition funnel (why raw recall is low = by design)
   B. domain-sharer non-separability (no cDNA-coverage threshold separates real
      UTR-divergent paralogs from domain-sharers -> the DNA definition is ill-posed).
Reads bench/family_def_dna_pr_edges.tsv + family_def_dna_pr.json.
Run: /home/juanfra/miniforge3/bin/python bench/family_def_dna_pr_fig.py
"""
import json
import os
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
DUMP = os.path.join(HERE, "family_def_dna_pr_edges.tsv")
JS = json.load(open(os.path.join(HERE, "family_def_dna_pr.json")))
OUT = os.path.join(HERE, "family_def_dna_pr_fig.png")

NAVY = "#1F2D5A"; TEAL = "#127A6E"; ORANGE = "#D96B27"; GREY = "#9aa0a6"

# ---- load edges ----
rows = []
with open(DUMP) as f:
    h = f.readline().rstrip("\n").split("\t")
    idx = {k: i for i, k in enumerate(h)}
    for line in f:
        c = line.rstrip("\n").split("\t")
        rows.append(c)
def col(c, k): return c[idx[k]]

fig, ax = plt.subplots(1, 2, figsize=(14, 5.4))

# ================= Panel A: recall funnel =================
fn = JS["funnel"]
tp = JS["results"]["loose"]["tp"]
segs = [
    ("recovered by RNA (TP)", tp, TEAL),
    ("sub-quorum (1-2 tied reads)", fn["4_sub_quorum (tied 1-2)"], "#E0A040"),
    ("resolvable (reads decidable)", fn["3_resolvable (cross-map, de-decidable)"], "#C98A6A"),
    ("reads place uniquely (no cross-map)", fn["2_no_cross_map (reads place uniquely)"], "#B0B6BD"),
    ("transcriptionally SILENT (unexpressed)", fn["1_unexpressed (>=1 copy no RNA)"], "#6B7280"),
]
total = sum(s[1] for s in segs)
left = 0
for name, val, color in segs:
    ax[0].barh(0, val, left=left, color=color, edgecolor="white", height=0.55)
    if val / total > 0.04:
        ax[0].text(left + val / 2, 0, f"{val:,}\n{100*val/total:.0f}%", ha="center", va="center",
                   fontsize=8.5, color="white", weight="bold")
    left += val
ax[0].set_xlim(0, total); ax[0].set_ylim(-1.1, 1.1); ax[0].set_yticks([])
ax[0].set_xlabel(f"DNA paralog edges (LOOSE truth: {total:,} pairs, id≥.90, cov≥.30)")
ax[0].set_title("A. Recall is low because the copies are SILENT, not missed", fontsize=12, weight="bold", color=NAVY)
# recall brackets
def bracket(x0, x1, y, txt, color):
    ax[0].annotate("", xy=(x1, y), xytext=(x0, y), arrowprops=dict(arrowstyle="<->", color=color, lw=1.4))
    ax[0].text((x0 + x1) / 2, y + 0.10, txt, ha="center", fontsize=8.5, color=color, weight="bold")
bracket(0, tp, 0.62, f"raw recall = {tp:,}/{total:,} = {tp/total:.2f}", NAVY)
crossmap = tp + fn["3_resolvable (cross-map, de-decidable)"] + fn["4_sub_quorum (tied 1-2)"]
bracket(0, crossmap, -0.62, f"recall | cross-mapping universe = {JS['recall_crossmap']:.2f}  (honest single figure)", TEAL)
ax[0].text(total*0.50, -0.95,
           "423 sub-quorum (1–2 tied reads) = the only genuine near-misses; "
           "524 'resolvable' = reads decide → excluded by choice",
           ha="center", fontsize=8.0, color=ORANGE, style="italic")
ax[0].legend(handles=[plt.Rectangle((0,0),1,1,color=c) for _,_,c in segs],
             labels=[n for n,_,_ in segs], fontsize=7.4, loc="upper center",
             bbox_to_anchor=(0.5, -0.16), ncol=2, frameon=False)

# ================= Panel B: domain-sharer non-separability =================
# scatter LOOSE-DNA edges in (min_cov, max_cov); color by RNA-confirmation.
xs_y, ys_y, xs_n, ys_n = [], [], [], []
for c in rows:
    if col(c, "in_dna_loose") != "1":
        continue
    ca = float(col(c, "cov_a")); cb = float(col(c, "cov_b"))
    mn, mx = min(ca, cb), max(ca, cb)
    if col(c, "in_rna") == "1":
        xs_y.append(mn); ys_y.append(mx)
    else:
        xs_n.append(mn); ys_n.append(mx)
ax[1].scatter(xs_n, ys_n, s=5, c=GREY, alpha=0.30, label=f"DNA-only (RNA says no)  n={len(xs_n):,}", rasterized=True)
ax[1].scatter(xs_y, ys_y, s=10, c=TEAL, alpha=0.55, label=f"RNA-confirmed family  n={len(xs_y):,}", rasterized=True)
# reciprocal (whole-gene) threshold: keeps points with min_cov>=0.5
ax[1].axvline(0.50, color=ORANGE, ls="--", lw=1.6)
ax[1].text(0.515, 0.04, "reciprocal-cov ≥0.50\n(‘whole-gene’ bar)", color=ORANGE, fontsize=8)
# named panel pairs
PANEL = {
    "RABL2A~RABL2B": (0.34, 0.38, "real paralog", TEAL),
    "CREB1~METTL21A": (0.54, 0.62, "domain-sharer", ORANGE),
    "ASDURF~ASNSD1": (0.13, 0.39, "domain-sharer", ORANGE),
    "GCA~KCNH7": (0.19, 0.70, "domain-sharer", ORANGE),
    "RFPL2~RFPL3": (0.50, 0.70, "real paralog", TEAL),
}
for name, (mn, mx, kind, color) in PANEL.items():
    ax[1].scatter([mn], [mx], s=90, marker="*", c=color, edgecolor="black", zorder=5, linewidths=0.6)
    dy = 0.05 if "RABL2" not in name else -0.07
    ax[1].annotate(name, (mn, mx), textcoords="offset points", xytext=(6, 6 if dy>0 else -12),
                   fontsize=7.6, color="black")
ax[1].set_xlabel("min coverage (the more-diverged copy)")
ax[1].set_ylabel("max coverage")
ax[1].set_xlim(0, 1.02); ax[1].set_ylim(0.25, 1.02)
ax[1].set_title("B. No CLEAN coverage cut: the distributions overlap",
                fontsize=12.5, weight="bold", color=NAVY)
ax[1].legend(fontsize=8, loc="lower right")
ax[1].text(0.035, 0.30, "RNA-confirmed median recip-cov 0.75\nvs DNA-only 0.27 — shifted but overlapping",
           fontsize=7.8, color=TEAL, style="italic")
ax[1].annotate("‘whole-gene’ bar DROPS RABL2 (real)\nyet KEEPS CREB1~METTL21A (domain-sharer)",
               xy=(0.34, 0.38), xytext=(0.04, 0.90), fontsize=8.0, color=NAVY,
               arrowprops=dict(arrowstyle="->", color=NAVY, lw=1))

fig.suptitle("Precision/recall vs an independent DNA-sequence ground truth "
             "(all-vs-all cDNA homology, same 34k gene vertices)",
             fontsize=12.5, weight="bold", y=1.0)
fig.tight_layout(rect=[0, 0, 1, 0.97])
fig.savefig(OUT, dpi=150, bbox_inches="tight")
print(f"[+] wrote {OUT}")
