#!/usr/bin/env python3
"""Explanatory figure: why POA contiguous-core coverage separates true
multi-copy paralogs (one long co-aligned block) from domain-sharers
(scattered chance matches whose all-column coverage looks misleadingly high).

Deck style: navy/orange/green/slate, DejaVu Sans, minimal text.
"""
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, Rectangle
import numpy as np

# ----- palette -----
NAVY = "#22313f"
ORANGE = "#e8590c"
GREEN = "#1e8449"
SLATE = "#8898aa"
GREEN_LT = "#a9dfbf"
SLATE_LT = "#dfe5ec"
WHITE = "#ffffff"

plt.rcParams.update({
    "font.family": "DejaVu Sans",
    "font.size": 10.5,
    "text.color": NAVY,
    "axes.edgecolor": NAVY,
})

fig, ax = plt.subplots(figsize=(10.2, 8.0))
ax.set_xlim(0, 100)
ax.set_ylim(0, 100)
ax.axis("off")
fig.patch.set_facecolor(WHITE)

rng = np.random.default_rng(7)


def exon_box(x, y, w, h, color, ec=NAVY, lw=1.2, z=3):
    ax.add_patch(FancyBboxPatch(
        (x, y), w, h, boxstyle="round,pad=0.0,rounding_size=0.6",
        facecolor=color, edgecolor=ec, linewidth=lw, zorder=z))


def gene_track(x0, y, exons, color, lw=1.2):
    """exons: list of (start_offset, width). Draw thin intron line + boxes."""
    xs = x0 + exons[0][0]
    xe = x0 + exons[-1][0] + exons[-1][1]
    ax.plot([xs, xe], [y + 1.6, y + 1.6], color=NAVY, lw=1.0, zorder=2)
    for (off, w) in exons:
        exon_box(x0 + off, y, w, 3.2, color)


# ============================================================
# TITLE
# ============================================================
ax.text(50, 97.5, "Why contiguous-core separates copies from domain-sharing",
        ha="center", va="top", fontsize=16.5, fontweight="bold", color=NAVY)

# ============================================================
# ROW 1 : TRUE COPIES (paralogs) -> MERGE
# ============================================================
r1_top = 90
# row label chip
ax.add_patch(FancyBboxPatch((2, r1_top - 4.0), 47, 4.0,
             boxstyle="round,pad=0.2,rounding_size=0.8",
             facecolor=GREEN, edgecolor="none", zorder=4))
ax.text(3.2, r1_top - 2.0, "TRUE COPIES (paralogs)", ha="left", va="center",
        fontsize=11.5, fontweight="bold", color=WHITE, zorder=5)
ax.text(60.5, r1_top - 2.0, "MERGE", ha="left", va="center",
        fontsize=13, fontweight="bold", color=GREEN, zorder=5)
ax.annotate("", xy=(60, r1_top - 2.0), xytext=(50, r1_top - 2.0),
            arrowprops=dict(arrowstyle="-|>", color=GREEN, lw=2.2))

# gene geometry (shared x axis 6..94)
gx0 = 6
# copy A and copy B: many homologous exons spanning most of the gene
exA = [(0, 9), (12, 7), (22, 8), (33, 9), (45, 7), (55, 9), (67, 7), (77, 9)]
exB = [(0, 9), (12, 7), (22, 8), (33, 9), (45, 7), (55, 9), (67, 7), (77, 9)]
yA1 = r1_top - 11
yB1 = r1_top - 18
gene_track(gx0, yA1, exA, NAVY)
gene_track(gx0, yB1, exB, NAVY)
ax.text(gx0 - 2.5, yA1 + 1.6, "copy A", ha="right", va="center", fontsize=9.5, color=NAVY)
ax.text(gx0 - 2.5, yB1 + 1.6, "copy B", ha="right", va="center", fontsize=9.5, color=NAVY)

# POA co-alignment track beneath: ONE long green block (contiguous core)
yTrack1 = r1_top - 25.5
core_x0 = gx0 + 0
core_w = 86  # spans most of both genes
# faint baseline of the alignment lane
ax.add_patch(Rectangle((gx0 - 0.5, yTrack1 - 0.4), 87.5, 3.6,
             facecolor=SLATE_LT, edgecolor="none", zorder=1))
exon_box(core_x0, yTrack1, core_w, 2.8, GREEN, ec=GREEN, lw=1.0, z=3)
ax.text(gx0 - 2.5, yTrack1 + 1.4, "POA", ha="right", va="center",
        fontsize=9.0, color=SLATE)
ax.text(core_x0 + core_w / 2, yTrack1 + 1.4, "contiguous core",
        ha="center", va="center", fontsize=10, fontweight="bold", color=WHITE, zorder=5)

# connector ticks (homology) between the two genes -> all consistent
for (off, w) in exA:
    cx = gx0 + off + w / 2
    ax.plot([cx, cx], [yB1, yA1 + 3.2], color=GREEN, lw=0.7, alpha=0.45, zorder=1)

# row-1 caption
ax.text(50, yTrack1 - 3.6,
        "longest co-aligned block = most of the gene  →  core coverage HIGH ( ≥ T )",
        ha="center", va="top", fontsize=11, color=GREEN, fontweight="bold")

# divider
ax.plot([4, 96], [49.0, 49.0], color=SLATE_LT, lw=1.4, zorder=0)

# ============================================================
# ROW 2 : DOMAIN-SHARERS (false family) -> REJECT
# ============================================================
r2_top = 46.5
ax.add_patch(FancyBboxPatch((2, r2_top - 4.0), 52, 4.0,
             boxstyle="round,pad=0.2,rounding_size=0.8",
             facecolor=ORANGE, edgecolor="none", zorder=4))
ax.text(3.2, r2_top - 2.0, "DOMAIN-SHARERS (false family)", ha="left", va="center",
        fontsize=11.5, fontweight="bold", color=WHITE, zorder=5)
ax.text(65.5, r2_top - 2.0, "REJECT", ha="left", va="center",
        fontsize=13, fontweight="bold", color=ORANGE, zorder=5)
ax.annotate("", xy=(65, r2_top - 2.0), xytext=(55, r2_top - 2.0),
            arrowprops=dict(arrowstyle="-|>", color=ORANGE, lw=2.2))

# two genes; share ONLY one small exon/domain
yA2 = r2_top - 11
yB2 = r2_top - 18
exA2 = [(0, 8), (11, 7), (21, 7), (31, 9), (43, 7), (53, 8), (64, 7), (74, 9)]
exB2 = [(0, 9), (12, 8), (23, 7), (33, 9), (45, 7), (55, 8), (66, 8), (77, 7)]
gene_track(gx0, yA2, exA2, NAVY)
gene_track(gx0, yB2, exB2, NAVY)
ax.text(gx0 - 2.5, yA2 + 1.6, "gene A", ha="right", va="center", fontsize=9.5, color=NAVY)
ax.text(gx0 - 2.5, yB2 + 1.6, "gene B", ha="right", va="center", fontsize=9.5, color=NAVY)

# the ONE shared small domain (highlight both copies green)
sh_off = 31
exon_box(gx0 + sh_off, yA2, 9, 3.2, GREEN)
exon_box(gx0 + 33, yB2, 9, 3.2, GREEN)
ax.plot([gx0 + sh_off + 4.5, gx0 + 33 + 4.5], [yB2 + 3.2, yA2],
        color=GREEN, lw=1.6, zorder=2)

# POA co-alignment track: ONE short green block + SCATTERED orange ticks
yTrack2 = r2_top - 25.5
ax.add_patch(Rectangle((gx0 - 0.5, yTrack2 - 0.4), 87.5, 3.6,
             facecolor=SLATE_LT, edgecolor="none", zorder=1))
ax.text(gx0 - 2.5, yTrack2 + 1.4, "POA", ha="right", va="center",
        fontsize=9.0, color=SLATE)
# the single real short core block
short_core_x = gx0 + 34
short_core_w = 6
exon_box(short_core_x, yTrack2, short_core_w, 2.8, GREEN, ec=GREEN, lw=1.0, z=4)
# scattered chance-match ticks across the rest of the lane
tick_positions = []
lane_x0, lane_x1 = gx0 + 1, gx0 + 85
xs = rng.uniform(lane_x0, lane_x1, 46)
for x in xs:
    if short_core_x - 1.5 < x < short_core_x + short_core_w + 1.5:
        continue
    h = rng.uniform(1.0, 2.6)
    ax.add_patch(Rectangle((x, yTrack2 + (2.8 - h) / 2 + 0.3), 0.55, h,
                 facecolor=ORANGE, edgecolor="none", zorder=3))
    tick_positions.append(x)

# bracket marking the short real core
ax.annotate("", xy=(short_core_x, yTrack2 - 0.9),
            xytext=(short_core_x + short_core_w, yTrack2 - 0.9),
            arrowprops=dict(arrowstyle="-", color=GREEN, lw=1.4))

# row-2 captions: misleading all-column vs true core
ax.text(50, yTrack2 - 3.3,
        "all-column coverage looks high (scattered chance matches sum up)  —  MISLEADING",
        ha="center", va="top", fontsize=10.2, color=ORANGE, fontweight="bold")
ax.text(50, yTrack2 - 6.9,
        "longest co-aligned block = a few %  →  core coverage LOW ( < T )  →  rejected",
        ha="center", va="top", fontsize=11, color=NAVY, fontweight="bold")

# ============================================================
# BOTTOM STRIP : criterion + result
# ============================================================
strip_y = 1.5
ax.add_patch(FancyBboxPatch((3, strip_y), 94, 7.2,
             boxstyle="round,pad=0.2,rounding_size=1.2",
             facecolor=NAVY, edgecolor="none", zorder=4))
ax.text(50, strip_y + 5.0,
        "POA contiguous-core coverage  =  longest matched (identical) co-aligned block / shorter gene   "
        "( T = 0.13 )",
        ha="center", va="center", fontsize=10.8, color=WHITE, fontweight="bold", zorder=5)
ax.text(50, strip_y + 2.0,
        "accepts 5/5 real copies      ✓          rejects 7/7 domain-sharers      ✗",
        ha="center", va="center", fontsize=11.5, color=GREEN_LT, fontweight="bold", zorder=5)

plt.tight_layout(pad=0.6)
fig.savefig("/mnt/c/Users/jfris/Desktop/Rustle/bench/fig_core_coverage.png",
            dpi=190, bbox_inches="tight", pad_inches=0.12, facecolor=WHITE)
print("wrote fig_core_coverage.png")
