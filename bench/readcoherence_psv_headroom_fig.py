#!/usr/bin/env python3
"""Figure: read-coherence recall vs PSV copy-resolution are DISJOINT (the headroom no-go).

Left  — funnel: 23,791 read-coherence recall isoforms collapse to 0 at real collapsed-paralog loci.
Right — the two levers act on disjoint locus populations (single-copy alt-isoforms vs collapsed paralogs).
Numbers come straight from bench/readcoherence_psv_headroom.py.
"""
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, Circle

NAVY = "#22313f"; SLATE = "#8898aa"; ORANGE = "#e8590c"
GREEN = "#1e8449"; RED = "#c0392b"; LIGHT = "#eceff1"
plt.rcParams["font.family"] = "DejaVu Sans"

fig = plt.figure(figsize=(13.2, 6.4))
axL = fig.add_axes([0.045, 0.07, 0.46, 0.82])
axR = fig.add_axes([0.560, 0.07, 0.42, 0.82])
for ax in (axL, axR):
    ax.set_xticks([]); ax.set_yticks([])
    for s in ax.spines.values():
        s.set_visible(False)

# ---------------- LEFT: funnel ----------------
axL.set_xlim(0, 10); axL.set_ylim(0, 10)
axL.text(5, 9.7, "read-coherence recall → PSV-resolvable?",
         fontsize=13.5, fontweight="bold", color=NAVY, ha="center")

# (label, value, width-fraction, color)
stages = [
    ("read-coherence recall isoforms", "23,791", 1.00, NAVY),
    ("at a multi-copy family locus", "88", 0.62, SLATE),
    ("at a COLLAPSED locus (PSVs could split)", "35", 0.40, ORANGE),
    ("", "0", 0.16, RED),
]
y = 8.4
dy = 1.95
cx = 5.0
maxw = 8.6
for i, (lab, val, frac, col) in enumerate(stages):
    w = maxw * (0.16 + 0.84 * frac)   # keep last bar visible
    h = 1.15
    axL.add_patch(FancyBboxPatch((cx - w / 2, y - h / 2), w, h,
                                 boxstyle="round,pad=0.02,rounding_size=0.08",
                                 linewidth=0, facecolor=col, alpha=0.92, zorder=3))
    axL.text(cx, y + 0.14, val, fontsize=19, fontweight="bold", color="white",
             ha="center", va="center", zorder=4)
    lab_fs = 9.3 if frac > 0.3 else 8.0
    axL.text(cx, y - 0.36, lab, fontsize=lab_fs, color="white",
             ha="center", va="center", zorder=4)
    if i < len(stages) - 1:
        axL.annotate("", xy=(cx, y - h / 2 - dy + h / 2 + 0.06),
                     xytext=(cx, y - h / 2 - 0.04),
                     arrowprops=dict(arrowstyle="-|>", color=SLATE, lw=2))
    y -= dy

# the last (narrowest) bar is centred at y = 8.4 - 3*1.95 = 2.55; label it to the side
axL.text(cx + 1.55, 2.55, "real paralog\n(not domain-sharer)",
         fontsize=8.2, color=RED, ha="left", va="center", fontweight="bold")
axL.text(5, 0.35, "0.147% reach a collapsed locus · 0 reach a real collapsed paralog",
         fontsize=9.4, color=RED, ha="center", style="italic", fontweight="bold")

# ---------------- RIGHT: disjoint channels ----------------
axR.set_xlim(0, 10); axR.set_ylim(0, 10)
axR.text(5, 9.7, "the two levers act on disjoint loci",
         fontsize=13.5, fontweight="bold", color=NAVY, ha="center")

# two barely-touching circles
c1 = Circle((3.55, 5.4), 2.55, facecolor=GREEN, alpha=0.16, edgecolor=GREEN, lw=2.2, zorder=2)
c2 = Circle((6.65, 5.4), 2.05, facecolor=ORANGE, alpha=0.16, edgecolor=ORANGE, lw=2.2, zorder=2)
axR.add_patch(c1); axR.add_patch(c2)

axR.text(3.0, 7.35, "read-coherence", fontsize=11.5, fontweight="bold", color=GREEN, ha="center")
axR.text(3.0, 5.85, "23,791", fontsize=17, fontweight="bold", color=GREEN, ha="center")
axR.text(3.0, 5.25, "recall isoforms", fontsize=9.2, color=NAVY, ha="center")
axR.text(3.0, 4.55, "single-copy genes:\nalt-splice, novel junctions",
         fontsize=8.4, color=SLATE, ha="center", va="center")

axR.text(7.05, 7.2, "PSV copy-split", fontsize=11.0, fontweight="bold", color=ORANGE, ha="center")
axR.text(7.05, 5.75, "collapsed", fontsize=10.5, fontweight="bold", color=ORANGE, ha="center")
axR.text(7.05, 5.2, "paralogs", fontsize=10.5, fontweight="bold", color=ORANGE, ha="center")
axR.text(7.05, 4.5, "needs shared\ncoordinate frame", fontsize=8.2, color=SLATE, ha="center", va="center")

# the intersection ~ empty
axR.text(5.15, 5.4, "∩ ≈ 0", fontsize=12.5, fontweight="bold", color=RED, ha="center", va="center",
         zorder=5)

axR.text(5, 1.85,
         "NO-GO for a unified molecule+PSV graph\njustified by combining recall + copy-resolution",
         fontsize=10.8, color=RED, ha="center", va="center", fontweight="bold",
         bbox=dict(boxstyle="round,pad=0.5", facecolor="#fdecea", edgecolor=RED, lw=1.3))
axR.text(5, 0.5,
         "read-coherence ships its recall additively today (--read-coherence);\nthe rebuild would buy architecture, not new recall or copy-resolution",
         fontsize=8.0, color=SLATE, ha="center", va="center", style="italic")

fig.suptitle("Headroom probe: does threading PSVs through the molecule graph pay off?",
             fontsize=14.5, fontweight="bold", color=NAVY, y=0.995)

out = "/mnt/c/Users/jfris/Desktop/Rustle/bench/readcoherence_psv_headroom.png"
fig.savefig(out, dpi=185, bbox_inches="tight", facecolor="white")
print("saved", out)
