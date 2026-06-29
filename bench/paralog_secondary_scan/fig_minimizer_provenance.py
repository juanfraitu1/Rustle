#!/usr/bin/env python3
"""Slide: ONE minimizer definition (minimap2's), used two ways — to ALIGN the reads,
and (same definition, reference-only) to PREDICT identifiability. Nothing custom.
"""
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch

BLUE, GREEN, INK = "#1f6fb2", "#2e9e5b", "#222"
fig = plt.figure(figsize=(13.6, 8.4)); ax = fig.add_axes([0,0,1,1]); ax.axis("off")
ax.set_xlim(0,100); ax.set_ylim(0,100)

ax.text(50, 96, "One minimizer definition — minimap2's — used two ways",
        ha="center", fontsize=16, fontweight="bold")
ax.text(50, 91.6, "We did NOT invent a hash. Both the alignment and the predictor use minimap2's own minimizers (mm_sketch: canonical k-mer + hash64, k=15, w=5).",
        ha="center", fontsize=10.3, style="italic", color="#444")

# shared definition banner
ax.add_patch(FancyBboxPatch((18,80),64,6.2,boxstyle="round,pad=0.3,rounding_size=2",fc="#f3f0e6",ec="#cbb96a",lw=1.4))
ax.text(50, 83.1, "minimap2 minimizer:  canonical k-mer (strand-symmetric)  →  hash64 bit-mixer  →  smallest per window  ·  k = 15, w = 5",
        ha="center", fontsize=10, family="monospace", fontweight="bold", color="#7a5c12")

# two arrows down
ax.add_patch(FancyArrowPatch((38,79.5),(27,69),arrowstyle="-|>",mutation_scale=16,color="#999",lw=1.6))
ax.add_patch(FancyArrowPatch((62,79.5),(73,69),arrowstyle="-|>",mutation_scale=16,color="#999",lw=1.6))

def box(x,y,w,h,fc,ec):
    ax.add_patch(FancyBboxPatch((x,y),w,h,boxstyle="round,pad=0.3,rounding_size=2.2",fc=fc,ec=ec,lw=1.6,zorder=2))

# LEFT: used to ALIGN
box(3, 30, 44, 39, "#eef5fb", BLUE)
ax.text(25, 65.0, "① used to ALIGN the reads", ha="center", fontsize=12.5, fontweight="bold", color=BLUE)
ax.text(25, 61.7, "(this already made GGO.bam)", ha="center", fontsize=9, color="#555", style="italic")
ax.add_patch(FancyBboxPatch((5,54.5),40,4.6,boxstyle="round,pad=0.2,rounding_size=1",fc="white",ec="#bcd",lw=1))
ax.text(25, 56.4, "minimap2 -ax splice:hq -uf  …   (@PG in the BAM)", ha="center",
        fontsize=8.2, family="monospace", color=INK)
for i,t in enumerate([
    "minimap2 seeds on these minimizers, chains,",
    "and places each read — marking it primary",
    "(best hit) or secondary (an equally-good hit",
    "on another copy).",
]):
    ax.text(6, 50.2-i*3.1, t, fontsize=9.5, color=INK)
ax.text(25, 35.0, "WHAT WE SEE IN THE BAM:", ha="center", fontsize=9.6, fontweight="bold", color=BLUE)
ax.text(25, 32.0, "the read pileup + primary/secondary colours.", ha="center", fontsize=9.3, color="#444")

# RIGHT: used to PREDICT
box(53, 30, 44, 39, "#eef8f1", GREEN)
ax.text(75, 65.0, "② used to PREDICT identifiability", ha="center", fontsize=12.5, fontweight="bold", color=GREEN)
ax.text(75, 61.7, "(same definition — but on the reference, no reads)", ha="center", fontsize=9, color="#555", style="italic")
for i,t in enumerate([
    "Take the SAME minimap2 minimizers of each",
    "copy from the reference genome, and ask:",
    "how many are NOT shared with the sibling copy?",
]):
    ax.text(6+50, 56.0-i*3.1, t, fontsize=9.5, color=INK)
ax.text(56, 45.4, "disc_frac high → tie-breakers exist (RABL2)", fontsize=9.5, color=GREEN, fontweight="bold")
ax.text(56, 42.3, "disc_frac ~0 → none exist; reads tie (DAZ)", fontsize=9.5, color="#c0392b", fontweight="bold")
ax.text(75, 35.0, "WHAT WE DRAW:", ha="center", fontsize=9.6, fontweight="bold", color=GREEN)
ax.text(75, 32.0, "the green/grey track — computed from the genome,\nbefore any read is assembled.", ha="center", fontsize=9.3, color="#444")

# bottom takeaway
box(3, 6, 94, 18, "#f7f7f8", "#cfd5da")
ax.text(50, 20.5, "The two sides use the IDENTICAL minimizer definition. The only difference is the INPUT:",
        ha="center", fontsize=11, fontweight="bold", color=INK)
ax.text(50, 16.0, "left consumes the reads (and reports primary/secondary);  right uses the reference sequence alone (and predicts whether a read COULD be placed).",
        ha="center", fontsize=9.8, color="#444")
ax.text(50, 11.7, "So the predictor is simply minimap2's seeding logic asked one step earlier — it agrees with the alignment because it IS the same minimizers.",
        ha="center", fontsize=9.8, color="#444")
ax.text(50, 8.0, "Nothing here is a new or unapproved method.",
        ha="center", fontsize=10, fontweight="bold", color=GREEN)

out="/mnt/c/Users/jfris/Desktop/Rustle/bench/paralog_secondary_scan/fig0_minimizer_provenance.png"
fig.savefig(out, dpi=160, facecolor="white"); print("wrote", out)
