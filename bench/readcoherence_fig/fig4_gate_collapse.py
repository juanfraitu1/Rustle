import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch
from matplotlib.lines import Line2D

# ---- palette (shared style guide) ----
NAVY   = "#22313f"   # annotation
SLATE  = "#8898aa"   # StringTie / flow
ORANGE = "#e8590c"   # read-coherence (accent)
LIGHT  = "#cfd8dc"   # read alignments
AMBER  = "#f0932b"   # spanning alignments
GREEN  = "#1e8449"   # gate-pass / canonical
RED    = "#c0392b"   # gate-fail / non-canonical

EH = 0.30            # exon half-height


def exon(ax, x0, x1, yc, color, h=EH):
    r = FancyBboxPatch((x0, yc - h), x1 - x0, 2 * h,
                       boxstyle="round,pad=0.0,rounding_size=0.06",
                       linewidth=0, facecolor=color, zorder=3)
    ax.add_patch(r)


def intron(ax, x0, x1, yc, color, lw=1.3):
    # peaked intron line (IGV-style hat)
    xm = (x0 + x1) / 2.0
    ax.plot([x0, xm, x1], [yc, yc + 0.14, yc], color=color, lw=lw,
            zorder=2, solid_capstyle="round", solid_joinstyle="round")


def chain(ax, exons, yc, exon_color, intron_color, lw=1.3):
    for (a, b) in zip(exons[:-1], exons[1:]):
        intron(ax, a[1], b[0], yc, intron_color, lw=lw)
    for (x0, x1) in exons:
        exon(ax, x0, x1, yc, exon_color)


def junc_label(ax, x0, x1, yc, donor, acceptor, color, dy=0.40):
    xm = (x0 + x1) / 2.0
    ax.text(xm, yc + dy, f"{donor}–{acceptor}", ha="center", va="bottom",
            fontsize=9.6, color=color, fontweight="bold",
            family="DejaVu Sans Mono")


# =================================================================
fig = plt.figure(figsize=(11.2, 8.0))
gs = fig.add_gridspec(2, 1, height_ratios=[1.05, 1.0], hspace=0.28)
axT = fig.add_subplot(gs[0])
axB = fig.add_subplot(gs[1])

for ax in (axT, axB):
    ax.set_xticks([]); ax.set_yticks([])
    for s in ax.spines.values():
        s.set_visible(False)

# =================================================================
# TOP: realness gate (annotation-free)
# =================================================================
axT.set_xlim(0, 12.6)
axT.set_ylim(0, 6.2)

axT.text(0.05, 5.95, "realness gate", fontsize=14, fontweight="bold",
         color=NAVY, ha="left", va="top")
axT.text(0.05, 5.55, "annotation-free: genome sequence + read depth only",
         fontsize=10.2, color=SLATE, ha="left", va="top", style="italic")

VERDICT_X = 11.35   # where check / cross + verdict sit
DEPTH_X   = 9.55    # read-depth column

# ---- candidate 1: PASS (all canonical, depth ok) ----
y1 = 4.55
chain(axT, [(0.5, 2.2), (4.0, 5.8), (7.4, 9.0)], y1, NAVY, "#4a5b6a")
junc_label(axT, 2.2, 4.0, y1, "GT", "AG", GREEN)
junc_label(axT, 5.8, 7.4, y1, "GC", "AG", GREEN)
axT.text(DEPTH_X, y1, "depth 7 ≥ 3", fontsize=9.8, color=GREEN,
         ha="center", va="center", fontweight="bold")
axT.text(VERDICT_X, y1, "✓", fontsize=20, color=GREEN, ha="center",
         va="center", fontweight="bold")
axT.text(VERDICT_X + 0.55, y1, "kept", fontsize=10.2, color=GREEN,
         ha="left", va="center", fontweight="bold")

# ---- candidate 2: DROP (non-canonical junction) ----
y2 = 3.05
chain(axT, [(0.5, 2.4), (4.3, 6.0), (7.4, 9.0)], y2, NAVY, "#4a5b6a")
junc_label(axT, 2.4, 4.3, y2, "GT", "AG", GREEN)
junc_label(axT, 6.0, 7.4, y2, "AC", "TG", RED)      # non-canonical -> red
axT.text(DEPTH_X, y2, "depth 6 ≥ 3", fontsize=9.8, color=SLATE,
         ha="center", va="center")
axT.text(VERDICT_X, y2, "✗", fontsize=18, color=RED, ha="center",
         va="center", fontweight="bold")
axT.text(VERDICT_X + 0.55, y2, "non-canonical", fontsize=9.6, color=RED,
         ha="left", va="center", fontweight="bold")

# ---- candidate 3: DROP (read-depth below min_cov) ----
y3 = 1.55
chain(axT, [(0.5, 2.2), (4.0, 5.8), (7.4, 9.0)], y3, "#b9c3cc", "#aab4bd")
junc_label(axT, 2.2, 4.0, y3, "GT", "AG", SLATE)
junc_label(axT, 5.8, 7.4, y3, "GT", "AG", SLATE)
axT.text(DEPTH_X, y3, "depth 1 < 3", fontsize=9.8, color=RED,
         ha="center", va="center", fontweight="bold")
axT.text(VERDICT_X, y3, "✗", fontsize=18, color=RED, ha="center",
         va="center", fontweight="bold")
axT.text(VERDICT_X + 0.55, y3, "below min_cov", fontsize=9.6, color=RED,
         ha="left", va="center", fontweight="bold")

# column headers (small)
axT.text(DEPTH_X, 5.25, "read depth", fontsize=9.4, color=NAVY,
         ha="center", va="center", style="italic")

# accepted canonical classes chip (small, bottom-left)
axT.text(0.5, 0.55, "accepted:", fontsize=9.6, color=NAVY, ha="left",
         va="center", style="italic")
for i, cls in enumerate(["GT–AG", "GC–AG", "AT–AC"]):
    cx = 1.85 + i * 1.45
    chip = FancyBboxPatch((cx - 0.6, 0.28), 1.2, 0.54,
                          boxstyle="round,pad=0.0,rounding_size=0.10",
                          linewidth=1.0, edgecolor=GREEN,
                          facecolor="#e7f3ec", zorder=2)
    axT.add_patch(chip)
    axT.text(cx, 0.55, cls, fontsize=9.4, color=GREEN, ha="center",
             va="center", fontweight="bold", family="DejaVu Sans Mono")

# =================================================================
# BOTTOM: degradation-aware collapse
# =================================================================
axB.set_xlim(0, 12.6)
axB.set_ylim(0, 6.2)

axB.text(0.05, 5.95, "degradation-aware collapse", fontsize=14,
         fontweight="bold", color=NAVY, ha="left", va="top")
axB.text(0.05, 5.55, "fragments fold into their full-length parent (ISM → FSM)",
         fontsize=10.2, color=SLATE, ha="left", va="top", style="italic")

# full-length parent (FSM) on the right.
# Use the parent exon coordinates as the reference frame for all fragments
# so each fragment reads as an exact sub-path (sharing exon boundaries).
PE = [(8.8, 10.4), (10.9, 11.9), (12.3, 12.4)]   # parent exons (right column)
# widen 3rd exon a touch for legibility
PE = [(8.8, 10.4), (10.9, 11.9), (12.2, 12.45)]
yP = 2.55
chain(axB, PE, yP, ORANGE, "#d98a2b", lw=1.5)
axB.text((PE[0][0] + PE[-1][1]) / 2, yP + 0.62, "full-length parent  (FSM)",
         fontsize=10.2, color=ORANGE, ha="center", va="bottom",
         fontweight="bold")

# fragments on the left: each is a sub-path of the parent, drawn in a
# left column. Boundaries are aligned to the parent's so collapse reads
# as exact containment.
FX = [(3.05, 4.65), (5.15, 6.15), (6.45, 6.7)]   # left-column exon coords
FRAGS = [
    # (label, parent-exon indices present, ycenter)
    ("5ʹ-truncated",   [1, 2],    4.75),   # missing first exon
    ("internal",        [1],       2.55),   # one internal exon only
    ("3ʹ-truncated",   [0, 1],    1.35),   # missing last exon
]
frag_y = []
for lbl, idxs, yc in FRAGS:
    exons = [FX[i] for i in idxs]
    chain(axB, exons, yc, AMBER, "#d98a2b", lw=1.3)
    axB.text(exons[0][0] - 0.2, yc, lbl, fontsize=9.6, color=NAVY,
             ha="right", va="center")
    frag_y.append((yc, exons[-1][1]))

# fold-in arrows: from each fragment's right end into the parent's left edge
target = (PE[0][0] - 0.2, yP)
for (yc, x_right) in frag_y:
    rad = 0.18 if yc > yP else (-0.18 if yc < yP else 0.0)
    a = FancyArrowPatch((x_right + 0.25, yc), target,
                        arrowstyle="-|>", mutation_scale=13,
                        connectionstyle=f"arc3,rad={rad}",
                        color=ORANGE, lw=1.6, zorder=4,
                        shrinkA=3, shrinkB=4, capstyle="round")
    axB.add_patch(a)

# label the fold action in the clear gap between the two columns (top of gap)
axB.text(7.55, 4.55, "fold into\nparent", fontsize=10.2, color=ORANGE,
         ha="center", va="center", fontweight="bold", linespacing=1.15)

# note in genuinely clear space (centered, below the lowest fragment)
axB.text(6.3, 0.35, "fragments  →  not separate isoforms",
         fontsize=9.8, color=SLATE, ha="center", va="center", style="italic")

# =================================================================
# shared legend (bottom, minimal)
# =================================================================
handles = [
    Line2D([0], [0], marker="s", color="none", markerfacecolor=NAVY,
           markersize=11, label="candidate / parent exon"),
    Line2D([0], [0], marker="s", color="none", markerfacecolor=AMBER,
           markersize=11, label="fragment exon"),
    Line2D([0], [0], marker="s", color="none", markerfacecolor=ORANGE,
           markersize=11, label="read-coherence"),
    Line2D([0], [0], marker="o", color="none", markerfacecolor=GREEN,
           markersize=11, label="canonical / pass"),
    Line2D([0], [0], marker="o", color="none", markerfacecolor=RED,
           markersize=11, label="non-canonical / drop"),
]
fig.legend(handles=handles, loc="lower center", ncol=5, frameon=False,
           fontsize=10.0, bbox_to_anchor=(0.5, -0.005), handletextpad=0.35,
           columnspacing=1.5)

fig.suptitle("Two precision mechanisms", fontsize=15,
             fontweight="bold", color=NAVY, y=0.985)

fig.subplots_adjust(left=0.02, right=0.985, top=0.93, bottom=0.07)
out = "/mnt/c/Users/jfris/Desktop/Rustle/bench/readcoherence_fig/fig4_gate_collapse.png"
fig.savefig(out, dpi=190, bbox_inches="tight", facecolor="white")
print("saved", out)
