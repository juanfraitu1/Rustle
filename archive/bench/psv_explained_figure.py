#!/usr/bin/env python3
"""Explanatory figure: what a PSV is, how read->copy assignment works, and how a PSV
is distinguished from a sequencing error. For the paralog-copy-recovery method."""
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle, FancyBboxPatch
import numpy as np

A_COL = "#1a7f37"   # matches Copy A allele (green)
B_COL = "#e8730c"   # matches Copy B allele (orange)
ERR   = "#d11"      # sequencing error (red)
SHARE = "#b9b9b9"   # shared (non-PSV) base

# 12 columns; PSVs at these columns (copies differ); rest shared
NCOL = 12
PSV = [2, 5, 8, 10]
shared = list("CTACGTACGTGA")          # shared backbone bases
copyA  = dict(zip(PSV, "GACT"))        # Copy A alleles at PSVs
copyB  = dict(zip(PSV, "AGTC"))        # Copy B alleles at PSVs

def base_at(copy, col):
    if col in PSV: return (copyA if copy=="A" else copyB)[col]
    return shared[col]

# reads: (label, origin, dict col->base override for errors, covered PSV cols)
reads = [
    ("read 1", "A", {}, PSV),
    ("read 2", "B", {}, PSV),
    ("read 3", "B", {5:"C"}, PSV),       # one sequencing error at PSV col 5 (C matches NEITHER copy: A=A, B=G)
    ("read 4", "A", {}, PSV),
    ("read 5", "B", {}, [2]),            # only spans 1 PSV -> below K -> unassigned
]
K = 3

fig, ax = plt.subplots(figsize=(12.5, 5.6))
ax.set_xlim(-3.5, NCOL+6.0); ax.set_ylim(-6.0, 2.9); ax.axis("off")
fig.suptitle("Paralog-sequence variants (PSVs): assigning reads to copies, and PSV vs. sequencing error",
             fontsize=12.5, fontweight="bold", y=0.97)

def cell(x, y, ch, color, bold=False, bg=None):
    if bg: ax.add_patch(Rectangle((x-0.45, y-0.42), 0.9, 0.84, facecolor=bg, edgecolor="none", zorder=1))
    ax.text(x, y, ch, ha="center", va="center", fontsize=12,
            color=color, fontweight=("bold" if bold else "normal"), family="monospace", zorder=3)

# PSV column highlight bands (centered on each PSV column; span Copy A row down to last read row)
band_top = 1.85
band_bot = -1.1 - (len(reads)-1)*0.95 - 0.5
for c in PSV:
    ax.add_patch(Rectangle((c-0.5, band_bot), 1.0, band_top-band_bot, facecolor="#fff3cf", edgecolor="none", zorder=0))
ax.text(NCOL/2-0.5, 2.4, "genomic position →", ha="center", fontsize=9, color="#666", style="italic")

# Copy A / Copy B reference rows
for c in range(NCOL):
    cell(c, 1.4, base_at("A", c), A_COL if c in PSV else SHARE, bold=(c in PSV))
    cell(c, 0.5, base_at("B", c), B_COL if c in PSV else SHARE, bold=(c in PSV))
ax.text(-0.9, 1.4, "Copy A", ha="right", va="center", fontweight="bold", fontsize=10)
ax.text(-0.9, 0.5, "Copy B", ha="right", va="center", fontweight="bold", fontsize=10)
for px in PSV: ax.text(px, 2.05, "PSV", ha="center", fontsize=7.5, color="#a9791b", fontweight="bold")

# divider
ax.plot([-3.3, NCOL-0.4], [-0.2, -0.2], color="#ccc", lw=1)

# reads
for i,(lab, origin, errs, covered) in enumerate(reads):
    y = -1.1 - i*0.95
    votesA=votesB=0; err_here=False
    for c in range(NCOL):
        if c not in covered:                       # read doesn't cover this column
            cell(c, y, "·", "#ddd"); continue
        b = errs.get(c, base_at(origin, c))
        if c in PSV:
            mA = (b==copyA[c]); mB = (b==copyB[c])
            if mA and not mB: col=A_COL; votesA+=1
            elif mB and not mA: col=B_COL; votesB+=1
            else: col=ERR; err_here=True
            cell(c, y, b, col, bold=True, bg=("#ffe2e2" if col==ERR else None))
        else:
            cell(c, y, b, SHARE)
    ax.text(-0.9, y, lab, ha="right", va="center", fontsize=9.5)
    # assignment verdict
    vmax, vmin = max(votesA,votesB), min(votesA,votesB)
    if vmax>=K and (vmax-vmin)>=1:
        who = "A" if votesA>votesB else "B"
        col = A_COL if who=="A" else B_COL
        txt = f"→ Copy {who}   ({votesA} A / {votesB} B votes)"
    else:
        col="#888"; txt=f"→ unassigned   ({votesA} A / {votesB} B; < K={K})"
    ax.text(NCOL+0.2, y, txt, ha="left", va="center", fontsize=9.5, color=col,
            fontweight=("bold" if col!="#888" else "normal"))

# one concise caption — figure-relative so it never clips
fig.text(0.5, 0.115,
         "PSV = fixed copy difference, consistent across a copy's reads.   Error = one random base "
         "(read 3, red) → outvoted by the other PSVs.   read 5 (< K=3 PSVs) → unassigned.",
         ha="center", fontsize=9.5, bbox=dict(boxstyle="round", fc="#f7f7f7", ec="#cccccc"))

# legend
from matplotlib.lines import Line2D
leg=[Line2D([0],[0],marker="s",color="w",markerfacecolor=A_COL,markersize=11,label="matches Copy A allele"),
     Line2D([0],[0],marker="s",color="w",markerfacecolor=B_COL,markersize=11,label="matches Copy B allele"),
     Line2D([0],[0],marker="s",color="w",markerfacecolor=ERR,markersize=11,label="sequencing error (matches neither)"),
     Line2D([0],[0],marker="s",color="w",markerfacecolor=SHARE,markersize=11,label="shared (non-PSV) base")]
fig.legend(handles=leg, loc="lower center", ncol=4, frameon=False, fontsize=9, bbox_to_anchor=(0.5,0.02))
fig.subplots_adjust(top=0.90, bottom=0.20, left=0.02, right=0.98)
out="/mnt/c/Users/jfris/Desktop/Rustle/bench/psv_explained_figure.png"
fig.savefig(out, dpi=160); print(f"wrote {out}")
