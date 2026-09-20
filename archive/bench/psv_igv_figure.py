#!/usr/bin/env python3
"""IGV-style mock pileup of BOTH paralog copy loci, so the read->copy assignment can be
corroborated by eye: a read is 'clean' (no PSV mismatches) at the copy it truly belongs to,
and shows a COLUMN of colored PSV mismatches at the sister copy. A lone scattered mismatch
is a sequencing error (random), not a PSV (a column)."""
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle, FancyBboxPatch
import numpy as np

# IGV base colors (read's base at a mismatch): A green, C blue, G brown, T red
IGV = {"A":"#4aa64a", "C":"#2b6cd4", "G":"#b8860b", "T":"#e23b3b"}
READ = "#c9d4e0"; READ_E = "#9fb0c2"; EXON = "#5b7aa8"

# gene model (3 exons) + 4 PSVs inside exons
EXONS = [(8,24),(40,60),(74,94)]
PSV   = [16, 46, 54, 84]
refA  = dict(zip(PSV, "GACT"))   # Copy A alleles
refB  = dict(zip(PSV, "AGTC"))   # Copy B alleles
ERR_READ, ERR_X, ERR_BASE = 4, 33, "T"   # one read carries a sequencing error at a non-PSV pos

# 12 reads: origin = the molecule's true copy (= where it should be assigned)
origins = ["A","B","A","B","B","A","B","A","A","B","A","B"]

def draw_locus(ax, title, ref, sister, gene_y):
    ax.set_xlim(-13, 100); ax.set_ylim(-0.4, gene_y+1.6); ax.axis("off")
    ax.text(-12, gene_y+1.0, title, ha="left", va="center", fontsize=10.5, fontweight="bold")
    # PSV guides
    for px in PSV:
        ax.plot([px,px],[ -0.2, gene_y+0.4], color="#e7c46b", lw=8, alpha=0.30, zorder=0)
        ax.text(px, gene_y+0.45, "PSV", ha="center", fontsize=6.5, color="#a9791b")
    # gene model: introns line + exon boxes
    ax.plot([EXONS[0][0], EXONS[-1][1]], [gene_y, gene_y], color=EXON, lw=1.2, zorder=1)
    for (x0,x1) in EXONS:
        ax.add_patch(Rectangle((x0, gene_y-0.22), x1-x0, 0.44, facecolor=EXON, edgecolor="none", zorder=2))
    # reads pileup
    for r,origin in enumerate(origins):
        y = gene_y - 1.0 - r*0.62
        # assignment dot (green=assigned A, orange=assigned B) — same in both panels
        ax.scatter([-9.5],[y], s=46, color=("#1a7f37" if origin=="A" else "#e8730c"), zorder=5)
        # read bar
        ax.add_patch(FancyBboxPatch((EXONS[0][0], y-0.20), EXONS[-1][1]-EXONS[0][0], 0.40,
                     boxstyle="round,pad=0,rounding_size=0.15", facecolor=READ, edgecolor="none", zorder=3))
        # mismatch ticks at PSVs where the read's base != THIS locus' reference
        rb = refA if origin=="A" else refB
        for px in PSV:
            base = rb[px]
            if base != ref[px]:
                ax.add_patch(Rectangle((px-0.7, y-0.20), 1.4, 0.40, facecolor=IGV[base], edgecolor="none", zorder=4))
        # one sequencing error (random, non-PSV) in a single read — shows in both panels
        if r == ERR_READ and ERR_X != ref.get(ERR_X):
            ax.add_patch(Rectangle((ERR_X-0.7, y-0.20), 1.4, 0.40, facecolor=IGV[ERR_BASE], edgecolor="none", zorder=4))

fig, (axA, axB) = plt.subplots(2, 1, figsize=(12.0, 8.4))
fig.suptitle("IGV-style pileup of both paralog copies — a read is 'clean' at the copy it belongs to, "
             "and mismatches the sister at every PSV", fontsize=12, fontweight="bold", y=0.985)
GY = len(origins)*0.62 + 1.0
draw_locus(axA, "Copy A   NC_073235.2:32,712,981–32,744,133", refA, refB, GY)
draw_locus(axB, "Copy B   NC_073235.2:31,684,679–31,709,055  (the starved copy)", refB, refA, GY)

# legends / caption
from matplotlib.lines import Line2D
base_leg = [Line2D([0],[0],marker="s",color="w",markerfacecolor=IGV[b],markersize=10,label=f"{b}") for b in "ACGT"]
asn_leg  = [Line2D([0],[0],marker="o",color="w",markerfacecolor="#1a7f37",markersize=9,label="assigned → Copy A"),
            Line2D([0],[0],marker="o",color="w",markerfacecolor="#e8730c",markersize=9,label="assigned → Copy B")]
fig.legend(handles=base_leg, loc="lower left", ncol=4, frameon=False, fontsize=8.5,
           bbox_to_anchor=(0.07,0.005), title="read base at mismatch (IGV colors)", title_fontsize=8.5)
fig.legend(handles=asn_leg, loc="lower right", ncol=2, frameon=False, fontsize=8.5, bbox_to_anchor=(0.95,0.01))
fig.text(0.5, 0.085, "Each dot = the read's assignment — corroborate by eye: a read is clean at its assigned copy,\n"
                     "and shows a column of PSV mismatches at the sister.   A lone off-column tick (read 5) = sequencing error.",
         ha="center", fontsize=9.0, style="italic", color="#444")
fig.subplots_adjust(top=0.94, bottom=0.13, left=0.02, right=0.99, hspace=0.18)
out="/mnt/c/Users/jfris/Desktop/Rustle/bench/psv_igv_figure.png"
fig.savefig(out, dpi=160); print(f"wrote {out}")
