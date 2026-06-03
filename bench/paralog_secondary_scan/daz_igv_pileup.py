#!/usr/bin/env python3
"""IGV-style pileup of the DAZ region from GGO.bam (gorilla IsoSeq).

Each read is drawn as IGV draws it: thick blocks = aligned exon segments,
thin connector line = intron (CIGAR N). Reads are greedily packed into rows
(a pileup). Colour = alignment class: PRIMARY (blue) vs SECONDARY (red).

Factual claim only: DAZ1 is piled with primary reads; DAZ3 is piled with
secondary reads. Read length is visible as bar extent.
"""
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
from matplotlib.lines import Line2D
import pysam

BAM = "/mnt/c/Users/jfris/Desktop/GGO.bam"
CHR = "NC_073248.2"
DAZ1 = (42_783_133, 42_859_657)   # - strand
DAZ3 = (42_879_918, 42_945_552)   # + strand (LOC129530216)
SPAN = (42_780_000, 42_948_000)

PRIM = "#1f6fb2"   # blue
SEC  = "#c0392b"   # red
mb = 1e6

def fetch():
    """Return list of (start, end, blocks, is_secondary) for primary+secondary alns."""
    reads = []
    with pysam.AlignmentFile(BAM, "rb") as bam:
        for a in bam.fetch(CHR, SPAN[0], SPAN[1]):
            if a.is_unmapped or a.is_supplementary:
                continue
            blocks = a.get_blocks()  # list of (ref_start, ref_end) for M/=/X
            if not blocks:
                continue
            reads.append((a.reference_start, a.reference_end, blocks, a.is_secondary))
    return reads

def pack(reads, pad=400):
    """Greedy IGV-style row packing. Returns list of (row, start, end, blocks, sec)."""
    reads = sorted(reads, key=lambda r: r[0])
    row_end = []          # last ref_end placed in each row
    placed = []
    for start, end, blocks, sec in reads:
        row = None
        for i, e in enumerate(row_end):
            if start > e + pad:
                row = i; row_end[i] = end; break
        if row is None:
            row = len(row_end); row_end.append(end)
        placed.append((row, start, end, blocks, sec))
    return placed, len(row_end)

reads = fetch()
n_prim = sum(1 for r in reads if not r[3])
n_sec  = sum(1 for r in reads if r[3])
print("fetched %d alns: %d primary, %d secondary" % (len(reads), n_prim, n_sec))

# Pack primary and secondary into SEPARATE stacks so the two classes read cleanly,
# primary on top (rows grow up), secondary mirrored below the axis.
prim = [r for r in reads if not r[3]]
sec  = [r for r in reads if r[3]]
pprim, nrp = pack(prim)
psec,  nrs = pack(sec)

fig, ax = plt.subplots(figsize=(13, 8.5))
H = 0.8   # bar thickness fraction of a row

def draw(packed, sign, color):
    for row, start, end, blocks, _sec in packed:
        y = sign * (row + 1)
        # intron connector (thin line spanning the read)
        ax.plot([start/mb, end/mb], [y, y], color=color, lw=0.5, alpha=0.6, zorder=1)
        # exon blocks (thick)
        for bs, be in blocks:
            ax.add_patch(Rectangle((bs/mb, y - H/2), (be - bs)/mb, H,
                                   facecolor=color, edgecolor="none", alpha=0.9, zorder=2))

draw(pprim, +1, PRIM)
draw(psec,  -1, SEC)

# zero line + locus shading
ax.axhline(0, color="0.3", lw=0.8)
ax.axvspan(DAZ1[0]/mb, DAZ1[1]/mb, color=PRIM, alpha=0.05, zorder=0)
ax.axvspan(DAZ3[0]/mb, DAZ3[1]/mb, color=SEC,  alpha=0.05, zorder=0)

ax.set_xlim(SPAN[0]/mb, SPAN[1]/mb)
ax.set_ylim(-(nrs + 4), nrp + 6)
ax.set_xlabel("%s position (Mb)" % CHR, fontsize=12)
ax.set_yticks([])
for s in ("top", "right", "left"):
    ax.spines[s].set_visible(False)

# side labels for the two stacks
ax.text(SPAN[0]/mb - 0.001, nrp/2, "PRIMARY\n%d reads" % n_prim, ha="right", va="center",
        fontsize=12, fontweight="bold", color=PRIM, rotation=0)
ax.text(SPAN[0]/mb - 0.001, -nrs/2, "SECONDARY\n%d reads" % n_sec, ha="right", va="center",
        fontsize=12, fontweight="bold", color=SEC, rotation=0)

# locus brackets/labels at top
c1 = (DAZ1[0]+DAZ1[1])/2/mb
c3 = (DAZ3[0]+DAZ3[1])/2/mb
ax.annotate("DAZ1  (−strand)", xy=(c1, 1), xycoords=("data", "axes fraction"),
            xytext=(0, 4), textcoords="offset points", ha="center", va="bottom",
            fontsize=12, fontweight="bold", color=PRIM)
ax.annotate("DAZ3 / LOC129530216  (+strand)", xy=(c3, 1), xycoords=("data", "axes fraction"),
            xytext=(0, 4), textcoords="offset points", ha="center", va="bottom",
            fontsize=12, fontweight="bold", color=SEC)
for lo, hi in (DAZ1, DAZ3):
    ax.annotate("", xy=(lo/mb, 1.005), xytext=(hi/mb, 1.005),
                xycoords=("data", "axes fraction"),
                arrowprops=dict(arrowstyle="-", color="0.5", lw=1))

leg = [Line2D([0], [0], color=PRIM, lw=6, label="Primary alignment"),
       Line2D([0], [0], color=SEC,  lw=6, label="Secondary alignment")]
ax.legend(handles=leg, loc="lower right", fontsize=11, frameon=False)

fig.suptitle("DAZ region pileup (IGV-style) — gorilla IsoSeq, GGO.bam",
             fontsize=15, fontweight="bold", y=0.98)
fig.text(0.5, 0.015,
         "Thick = aligned exon blocks, thin = intron (full-length spliced reads). "
         "DAZ1 is piled with PRIMARY alignments; DAZ3 is covered almost entirely by SECONDARY alignments.",
         ha="center", fontsize=10, style="italic", color="0.35")

fig.subplots_adjust(left=0.10, right=0.985, top=0.92, bottom=0.075)
out = "/mnt/c/Users/jfris/Desktop/Rustle/bench/paralog_secondary_scan/daz_igv_pileup.png"
fig.savefig(out, dpi=160)
print("wrote", out)
