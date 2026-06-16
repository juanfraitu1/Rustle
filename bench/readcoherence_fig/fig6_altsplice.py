import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch
from matplotlib.path import Path
from matplotlib.patches import PathPatch

# ---- palette ----
NAVY   = "#22313f"   # exon / annotation
ORANGE = "#e8590c"   # highlighted / differing junction
SLATE  = "#8898aa"   # normal junction & secondary
RED    = "#c0392b"   # drop / non-canonical (unused here, reserved)

plt.rcParams["font.family"] = "DejaVu Sans"

EH = 0.26   # exon half-height


def exon_box(ax, x0, x1, yc, color=NAVY, dashed=False, fill=True, lw=1.4):
    """Solid (filled) or dashed-outline exon box."""
    if fill and not dashed:
        r = FancyBboxPatch((x0, yc - EH), x1 - x0, 2 * EH,
                           boxstyle="round,pad=0.0,rounding_size=0.06",
                           linewidth=0, edgecolor="white",
                           facecolor=color, zorder=3)
        ax.add_patch(r)
    else:
        r = FancyBboxPatch((x0, yc - EH), x1 - x0, 2 * EH,
                           boxstyle="round,pad=0.0,rounding_size=0.06",
                           linewidth=lw, edgecolor=color,
                           facecolor="none", linestyle=(0, (4, 3)),
                           zorder=3)
        ax.add_patch(r)


def junction_line(ax, x0, x1, yc, color=SLATE, lw=1.4, z=2):
    ax.plot([x0, x1], [yc, yc], color=color, lw=lw, zorder=z,
            solid_capstyle="round")


def skip_arc(ax, x0, x1, yc, color=SLATE, lw=1.4, h=0.95, z=2):
    """Curved arc above the baseline from x0 to x1."""
    xm = (x0 + x1) / 2.0
    verts = [(x0, yc), (xm, yc + h), (x1, yc)]
    codes = [Path.MOVETO, Path.CURVE3, Path.CURVE3]
    ax.add_patch(PathPatch(Path(verts, codes), facecolor="none",
                           edgecolor=color, lw=lw, zorder=z,
                           capstyle="round"))


# =====================================================================
# Toy integer coordinates (shared across all rows so the differing
# junction is visibly a different number between isoforms).
# Exon spans (start, end):
#   E1: 10-20
#   E2: 30-40 (donor d_a = 40)   d_b = 45 (alt 5' donor, distal)
#   E3: 55-65 (cassette)
#   E4: 80-90 (acceptor a_a = 80) a_b = 75 (alt 3' acceptor, distal)
#   E5: 100-110
# Junctions are (donor_end -> acceptor_start) intron coordinates.
# =====================================================================
E1 = (10, 20)
E2 = (30, 40); D_A = 40; D_B = 45
E3 = (55, 65)
E4 = (80, 90); A_A = 80; A_B = 75
E5 = (100, 110)

fig = plt.figure(figsize=(11.8, 8.2))

# --- two stacked regions: gene model (top), isoforms (bottom) ---
axT = fig.add_axes([0.04, 0.66, 0.92, 0.26])   # top gene model
axB = fig.add_axes([0.04, 0.05, 0.92, 0.55])   # bottom isoforms

for ax in (axT, axB):
    ax.set_xticks([]); ax.set_yticks([])
    for s in ax.spines.values():
        s.set_visible(False)

# left gutter reserved for isoform row labels; exon track shifted right by GUT
GUT = 44

# coordinate -> plot-x mapping (toy coord + gutter offset)
def X(c):
    return c + GUT

XMIN, XMAX = 4, 132 + GUT

# =====================================================================
# TOP: gene model / splice graph
# =====================================================================
axT.set_xlim(XMIN, XMAX)
axT.set_ylim(-1.55, 1.85)
yc = 0.0

# baseline introns (slate) connecting the canonical chain
junction_line(axT, X(E1[1]), X(E2[0]), yc)      # E1 -> E2
junction_line(axT, X(D_A), X(E3[0]), yc)        # E2(d_a) -> E3
junction_line(axT, X(E3[1]), X(A_A), yc)        # E3 -> E4(a_a)
junction_line(axT, X(E4[1]), X(E5[0]), yc)      # E4 -> E5

# exons
exon_box(axT, X(E1[0]), X(E1[1]), yc)                    # E1 solid
exon_box(axT, X(E2[0]), X(D_A), yc)                      # E2 core solid (to d_a)
# E2 dashed extension to d_b (alt 5' donor)
exon_box(axT, X(D_A), X(D_B), yc, color=ORANGE, dashed=True)
exon_box(axT, X(E3[0]), X(E3[1]), yc, color=NAVY, dashed=True)  # E3 cassette
# E4: dashed distal acceptor extension (a_b proximal-to-graph, distal in seq)
exon_box(axT, X(A_B), X(A_A), yc, color=ORANGE, dashed=True)
exon_box(axT, X(A_A), X(E4[1]), yc)                      # E4 core solid (from a_a)
exon_box(axT, X(E5[0]), X(E5[1]), yc)                    # E5 solid

# skip arc E2 -> E4 over E3
skip_arc(axT, X(D_A), X(A_A), yc, color=SLATE, lw=1.4, h=1.15)

# exon labels (inside boxes, navy ones white)
def elabel(ax, x0, x1, yc, txt, white=True):
    ax.text((x0 + x1) / 2, yc, txt, fontsize=9.5,
            color="white" if white else NAVY, ha="center", va="center",
            fontweight="bold", zorder=5)

elabel(axT, X(E1[0]), X(E1[1]), yc, "E1")
elabel(axT, X(E2[0]), X(D_A), yc, "E2")
elabel(axT, X(E3[0]), X(E3[1]), yc, "E3", white=False)
elabel(axT, X(A_A), X(E4[1]), yc, "E4")
elabel(axT, X(E5[0]), X(E5[1]), yc, "E5")

# branch-point annotations (small italic)
axT.text(X((D_A + D_B) / 2), yc + 0.75, "alt 5' donor", fontsize=8.6,
         color=ORANGE, style="italic", ha="center", va="bottom")
axT.annotate("", xy=(X(D_B), yc + 0.40), xytext=(X((D_A + D_B) / 2), yc + 0.72),
             arrowprops=dict(arrowstyle="-", color=ORANGE, lw=0.9))

axT.text(X((D_A + A_A) / 2), yc + 1.30, "cassette (skippable)", fontsize=8.6,
         color=SLATE, style="italic", ha="center", va="bottom")

axT.text(X((A_B + A_A) / 2), yc - 0.95, "alt 3' acceptor", fontsize=8.6,
         color=ORANGE, style="italic", ha="center", va="top")
axT.annotate("", xy=(X(A_B), yc - 0.40), xytext=(X((A_B + A_A) / 2), yc - 0.92),
             arrowprops=dict(arrowstyle="-", color=ORANGE, lw=0.9))

axT.set_title("Alternative splicing = distinct junction chains",
              fontsize=15, fontweight="bold", color=NAVY, pad=12)

# =====================================================================
# BOTTOM: four isoforms, one per row.
# Each row: left label, exon chain w/ the differing junction in orange,
# right: exact ordered junction chain as a monospace tuple list.
# =====================================================================
# share the SAME x-axis as the top panel so exon columns line up vertically;
# add headroom on the right for the monospace tuple lists.
axB.set_xlim(XMIN, XMAX + 58)
axB.set_ylim(-0.6, 4.4)

rows_y = [3.7, 2.5, 1.3, 0.1]

LBL_X = XMIN + 1            # label sits in the left gutter (track starts at X(10))
CHAIN_X = X(E5[1]) + 14     # monospace junction chain, right of the exon track


def draw_isoform(ax, yc, label, segments, hi_idx, chain_str):
    """
    segments: list of exon (x0, x1) boxes in toy coords, in order.
    Between consecutive exons we draw a junction line; hi_idx marks which
    inter-exon junction (0-based among the drawn junctions) is highlighted.
    A skip is just a longer junction line (handled by the exon gap).
    """
    seg = [(X(a), X(b)) for (a, b) in segments]
    # junctions between consecutive boxes
    n_j = len(seg) - 1
    for j in range(n_j):
        x_from = seg[j][1]
        x_to = seg[j + 1][0]
        col = ORANGE if j == hi_idx else SLATE
        lw = 2.4 if j == hi_idx else 1.4
        junction_line(ax, x_from, x_to, yc, color=col, lw=lw)
    # exons
    for (x0, x1) in seg:
        exon_box(ax, x0, x1, yc)
    # left label (in gutter)
    ax.text(LBL_X, yc, label, fontsize=10.0, color=NAVY, ha="left",
            va="center", fontweight="bold")
    # right monospace junction chain
    ax.text(CHAIN_X, yc, chain_str, fontsize=8.7, color=NAVY, ha="left",
            va="center", family="monospace")


# Junction tuples are (donor_end, acceptor_start) intron coords.
# Canonical baseline chain for reference:
#   (20,30) (40,55) (65,80) (90,100)

# 1. inclusion : E1-E2(d_a)-E3-E4(a_a)-E5  -- no highlight (reference)
draw_isoform(
    axB, rows_y[0], "inclusion",
    [(E1[0], E1[1]), (E2[0], D_A), (E3[0], E3[1]), (A_A, E4[1]), (E5[0], E5[1])],
    hi_idx=-1,
    chain_str="(20,30) (40,55) (65,80) (90,100)",
)

# 2. exon skip (E3) : E1-E2(d_a)---skip---E4(a_a)-E5
#    highlight the long E2->E4 junction (index 1 among 3 junctions)
draw_isoform(
    axB, rows_y[1], "exon skip (E3)",
    [(E1[0], E1[1]), (E2[0], D_A), (A_A, E4[1]), (E5[0], E5[1])],
    hi_idx=1,
    chain_str="(20,30) (40,80) (90,100)",
)

# 3. alt 5' donor : E1-E2(d_b)-E3-E4(a_a)-E5
#    E2 extends to d_b; highlight the E2(d_b)->E3 junction (index 1)
draw_isoform(
    axB, rows_y[2], "alt 5' donor",
    [(E1[0], E1[1]), (E2[0], D_B), (E3[0], E3[1]), (A_A, E4[1]), (E5[0], E5[1])],
    hi_idx=1,
    chain_str="(20,30) (45,55) (65,80) (90,100)",
)

# 4. alt 3' acceptor : E1-E2(d_a)-E3-E4(a_b)-E5
#    E4 starts at a_b (distal); highlight the E3->E4(a_b) junction (index 2)
draw_isoform(
    axB, rows_y[3], "alt 3' acceptor",
    [(E1[0], E1[1]), (E2[0], D_A), (E3[0], E3[1]), (A_B, E4[1]), (E5[0], E5[1])],
    hi_idx=2,
    chain_str="(20,30) (40,55) (65,75) (90,100)",
)

# tiny header above the tuple column
axB.text(CHAIN_X, rows_y[0] + 0.55, "junction chain (exact key)", fontsize=8.4,
         color=SLATE, style="italic", ha="left", va="center")

# caption (<=12 words)
fig.text(0.5, 0.012,
         "each event changes a junction → a different exact key → a distinct transcript",
         fontsize=10.5, color=NAVY, ha="center", va="bottom")

out = "/mnt/c/Users/jfris/Desktop/Rustle/bench/readcoherence_fig/fig6_altsplice.png"
fig.savefig(out, dpi=190, bbox_inches="tight", facecolor="white")
print("saved", out)
