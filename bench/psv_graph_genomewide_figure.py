#!/usr/bin/env python3
"""Figure for the genome-wide PSV resolution scan (VERIFIED numbers).
Reads psv_graph_verify.json (cross-family-deduped, no_psv-triaged).
Panel A: refined family-class split. Panel B: identifiability frontier (#copies vs K).
Panel C: read assignment over unique families.
Run: /home/juanfra/miniforge3/bin/python bench/psv_graph_genomewide_figure.py
"""
import collections
import json
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
NAVY, TEAL, ORANGE, GREY, RED = "#1F2D5A", "#127A6E", "#D96B27", "#9aa0a6", "#B23B3B"
v = json.load(open(os.path.join(HERE, "psv_graph_verify.json")))
fams = v["families_refined"]
N = v["unique_families"]
nres = v["resolvable"]
ngk0 = v["no_psv_triage"].get("genuine_K0", 0)
naf = v["no_psv_triage"].get("align_fail", 0)
pct_res = round(100 * nres / N, 1)
pct_k0 = round(100 * ngk0 / N, 1)

CLS_COL = {"fully_resolvable": TEAL, "partial": ORANGE, "genuine_K0": GREY, "align_fail": "#d9c3b0"}

fig = plt.figure(figsize=(14.5, 4.9))
gs = fig.add_gridspec(1, 3, width_ratios=[1.05, 1.25, 0.95], wspace=0.34)

# ---- Panel A: refined classes ----
axA = fig.add_subplot(gs[0])
cnt = collections.Counter(f["cls_refined"] for f in fams)
order = ["fully_resolvable", "partial", "genuine_K0", "align_fail"]
labs = ["fully\nresolvable", "partial", "genuine\nK=0", "indet.\n(aln fail)"]
cols = [CLS_COL[k] for k in order]
vals = [cnt.get(k, 0) for k in order]
bars = axA.bar(range(len(order)), vals, color=cols, edgecolor="white")
for b, val in zip(bars, vals):
    axA.text(b.get_x() + b.get_width() / 2, val + 0.6, str(val),
             ha="center", va="bottom", fontsize=10.5, weight="bold", color="#333")
axA.set_xticks(range(len(order)))
axA.set_xticklabels(labs, fontsize=9)
axA.set_ylabel("families", fontsize=10)
axA.set_ylim(0, max(vals) * 1.16)
axA.set_title(f"A  {N} unique multi-copy families\n"
              f"{pct_res}% read-resolvable · {pct_k0}% genuine K-frontier",
              fontsize=11, weight="bold", color=NAVY, loc="left")
for sp in ("top", "right"):
    axA.spines[sp].set_visible(False)

# ---- Panel B: identifiability frontier (#copies vs K) ----
axB = fig.add_subplot(gs[1])
jit = collections.Counter()
for f in fams:
    x, y = f["n_copies"], f["K"]
    o = jit[(x, y)] * 0.07
    jit[(x, y)] += 1
    axB.scatter(x + o - 0.14, y + o - 0.14, s=20 + 0.02 * f["reads"],
                color=CLS_COL.get(f["cls_refined"], GREY), alpha=0.62,
                edgecolor="none", zorder=3)
mx = max(f["n_copies"] for f in fams)
axB.plot([1, mx], [1, mx], "--", color=NAVY, lw=1.2, zorder=1,
         label="K = #copies  (fully resolvable)")
axB.axhline(1, ls=":", color=RED, lw=1.3, zorder=1, label="K = 1  (collapsed / K-frontier)")
axB.set_xlabel("# copies in family", fontsize=10)
axB.set_ylabel("K = resolvable haplotypes", fontsize=10)
axB.set_title("B  the identifiability frontier\n(each point a family; size ∝ reads)",
              fontsize=11, weight="bold", color=NAVY, loc="left")
axB.legend(fontsize=8, loc="upper left", frameon=False)
axB.set_xlim(1.4, mx + 0.6)
axB.set_ylim(0.4, mx + 0.6)
for sp in ("top", "right"):
    axB.spines[sp].set_visible(False)

# ---- Panel C: read assignment ----
axC = fig.add_subplot(gs[2])
rlabs = ["→ ONE\ncopy", "→ collapsed\ngroup", "no PSV\ncovered"]
rvals = [v["reads_to_one_copy"], v["reads_to_collapsed_group"], v["reads_unassignable"]]
bars = axC.bar(range(3), rvals, color=[TEAL, GREY, "#e3e3e3"], edgecolor="white")
for b, val in zip(bars, rvals):
    axC.text(b.get_x() + b.get_width() / 2, val + max(rvals) * 0.015, f"{val:,}",
             ha="center", va="bottom", fontsize=9.5, weight="bold", color="#333")
axC.set_xticks(range(3))
axC.set_xticklabels(rlabs, fontsize=9)
axC.set_ylabel("reads", fontsize=10)
axC.set_ylim(0, max(rvals) * 1.16)
axC.set_title(f"C  {v['total_reads']:,} reads threaded\n"
              f"single-copy calls {v['single_copy_validation_pct']}% agree w/ best map",
              fontsize=11, weight="bold", color=NAVY, loc="left")
for sp in ("top", "right"):
    axC.spines[sp].set_visible(False)

fig.suptitle("PSV-aware variation graphs, genome-wide: where great-ape gene-family copies are "
             "read-resolvable — and where they hit the K-frontier",
             fontsize=12.5, weight="bold", color=NAVY, y=1.05)
fig.savefig(os.path.join(HERE, "psv_graph_genomewide.png"), dpi=150, bbox_inches="tight")
print("[+] wrote psv_graph_genomewide.png")
