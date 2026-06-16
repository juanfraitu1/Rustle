import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patheffects as pe
from matplotlib.patches import FancyBboxPatch

# white halo so colored allele letters stay legible on navy exon boxes
HALO = [pe.withStroke(linewidth=2.4, foreground="white")]

# ---- palette ----
NAVY   = "#22313f"   # exon / annotation
ORANGE = "#e8590c"   # read-coherence accent / Copy B
SLATE  = "#8898aa"   # PSV guide lines / muted
BLUE   = "#2c7fb8"   # Copy A

plt.rcParams["font.family"] = "DejaVu Sans"

# ---------------------------------------------------------------------------
# Shared geometry (toy coords) for the 3-exon read structure, reused L & R.
# ---------------------------------------------------------------------------
EXONS = [(6, 20), (30, 44), (54, 68)]   # (x0, x1) of the three exons
RH = 0.16                                # read box half-height
N_READS = 8
ROW_TOP = 7.4
ROW_DY = 0.62                            # vertical gap between reads
READ_YS = [ROW_TOP - i * ROW_DY for i in range(N_READS)]

# PSV columns (x positions inside the exons) and per-copy alleles
PSV_X = [13, 37, 61]                     # one PSV per exon
PSV_LABELS = ["P1", "P2", "P3"]
ALLELE_A = ["A", "C", "T"]              # Copy A allele vector (blue)
ALLELE_B = ["G", "T", "A"]             # Copy B allele vector (orange)

# transcript box geometry
TH = 0.30                                # transcript box half-height


def junction_segments(exons, psv_x):
    """
    Return the x-intervals where a horizontal connector line is drawn between
    consecutive exons. (Used so PSV guides only ever cross exon gaps, but kept
    here for clarity / symmetry.)
    """
    segs = []
    for i in range(len(exons) - 1):
        segs.append((exons[i][1], exons[i + 1][0]))
    return segs


def read_row(ax, yc, exons):
    """One read: three navy exon boxes joined by two thin junction lines."""
    for i in range(len(exons) - 1):
        x_from = exons[i][1]
        x_to = exons[i + 1][0]
        ax.plot([x_from, x_to], [yc, yc], color=SLATE, lw=1.0,
                zorder=2, solid_capstyle="round")
    for (x0, x1) in exons:
        r = FancyBboxPatch((x0, yc - RH), x1 - x0, 2 * RH,
                           boxstyle="round,pad=0.0,rounding_size=0.05",
                           linewidth=0, facecolor=NAVY, zorder=3)
        ax.add_patch(r)


# ===========================================================================
# Figure scaffold: two panels side by side + a bottom key strip.
# ===========================================================================
fig = plt.figure(figsize=(12.4, 8.4))

axL = fig.add_axes([0.045, 0.165, 0.435, 0.75])
axR = fig.add_axes([0.520, 0.165, 0.435, 0.75])
axK = fig.add_axes([0.045, 0.010, 0.910, 0.135])   # bottom key strip

XMIN, XMAX = 0, 74
YMIN, YMAX = -0.2, 8.7
for ax in (axL, axR):
    ax.set_xlim(XMIN, XMAX)
    ax.set_ylim(YMIN, YMAX)
    ax.set_xticks([]); ax.set_yticks([])
    for s in ax.spines.values():
        s.set_visible(False)

axK.set_xlim(0, 1); axK.set_ylim(0, 1)
axK.set_xticks([]); axK.set_yticks([])
for s in axK.spines.values():
    s.set_visible(False)

# ---------------------------------------------------------------------------
# Layout anchors shared by both panels (so reads & transcripts never overlap).
# ---------------------------------------------------------------------------
READ_STACK_BOT = READ_YS[-1] - RH      # bottom edge of the lowest read box
ARROW_TOP = READ_STACK_BOT - 0.18      # arrow starts just below the read stack
ARROW_BOT = 1.85                       # arrow ends above the transcript region

# ---------------------------------------------------------------------------
# LEFT PANEL: read-coherence (junctions only) -> 1 transcript
# ---------------------------------------------------------------------------
axL.text((XMIN + XMAX) / 2, 8.45, "read-coherence (junctions only)",
         fontsize=13, fontweight="bold", color=NAVY, ha="center", va="center")

for y in READ_YS:
    read_row(axL, y, EXONS)

# arrow down to the single transcript (points DOWN: from above to below)
axL.annotate("", xy=((XMIN + XMAX) / 2, ARROW_BOT),
             xytext=((XMIN + XMAX) / 2, ARROW_TOP),
             arrowprops=dict(arrowstyle="-|>", color=SLATE, lw=2.0,
                             shrinkA=0, shrinkB=0))

# single navy transcript box (3-exon)
TX_Y = 1.30
tx_exons = EXONS
for i in range(len(tx_exons) - 1):
    axL.plot([tx_exons[i][1], tx_exons[i + 1][0]], [TX_Y, TX_Y],
             color=NAVY, lw=1.4, zorder=2, solid_capstyle="round")
for (x0, x1) in tx_exons:
    r = FancyBboxPatch((x0, TX_Y - TH), x1 - x0, 2 * TH,
                       boxstyle="round,pad=0.0,rounding_size=0.08",
                       linewidth=0, facecolor=NAVY, zorder=3)
    axL.add_patch(r)
axL.text((tx_exons[0][0] + tx_exons[-1][1]) / 2, TX_Y - TH - 0.30,
         "1 transcript", fontsize=11, color=NAVY, ha="center", va="top",
         fontweight="bold")

axL.text((XMIN + XMAX) / 2, 0.05,
         "near-identical paralogs share the junction chain\n"
         "→ one bucket (copies merged)",
         fontsize=9.2, color=SLATE, ha="center", va="center")

# ---------------------------------------------------------------------------
# RIGHT PANEL: read-coherence + PSVs -> two copies
# ---------------------------------------------------------------------------
axR.text((XMIN + XMAX) / 2, 8.45, "read-coherence + PSVs",
         fontsize=13, fontweight="bold", color=NAVY, ha="center", va="center")

# PSV vertical dashed guide lines across the whole read stack.
# Drawn BEHIND the exon boxes (low zorder): the opaque navy boxes hide the line
# wherever it crosses a block, so it never strikes through the allele letters.
# It stays visible only in the inter-exon gaps + the label cap above the stack.
psv_line_top = READ_YS[0] + RH + 0.05
psv_line_bot = READ_YS[-1] - RH - 0.05
for px, plab in zip(PSV_X, PSV_LABELS):
    axR.plot([px, px], [psv_line_bot, psv_line_top], color=SLATE, lw=0.9,
             linestyle=(0, (3, 3)), zorder=1)
    axR.text(px, psv_line_top + 0.18, plab, fontsize=8.6, color=SLATE,
             ha="center", va="bottom", zorder=8)

# draw the reads (zorder 2-3) ON TOP of the dashed guides
for idx, y in enumerate(READ_YS):
    read_row(axR, y, EXONS)

# the 8 reads: first 4 = Copy A (blue), last 4 = Copy B (orange)
for idx, y in enumerate(READ_YS):
    is_a = idx < 4
    alleles = ALLELE_A if is_a else ALLELE_B
    col = BLUE if is_a else ORANGE
    for px, base in zip(PSV_X, alleles):
        axR.text(px, y, base, fontsize=9.5, color=col, ha="center",
                 va="center", fontweight="bold", zorder=7,
                 path_effects=HALO)

# arrow down to two transcripts (split) — points DOWN
axR.annotate("", xy=((XMIN + XMAX) / 2, ARROW_BOT),
             xytext=((XMIN + XMAX) / 2, ARROW_TOP),
             arrowprops=dict(arrowstyle="-|>", color=SLATE, lw=2.0,
                             shrinkA=0, shrinkB=0))


# two transcript boxes, stacked: Copy A (blue, top) and Copy B (orange, lower).
# Each is a single solid block per exon with ONE centered white allele letter,
# so Copy A and Copy B are symmetric (contrast is color only).
def copy_transcript(ax, yc, color, label, alleles):
    for i in range(len(EXONS) - 1):
        ax.plot([EXONS[i][1], EXONS[i + 1][0]], [yc, yc],
                color=color, lw=1.4, zorder=2, solid_capstyle="round")
    for (x0, x1) in EXONS:
        r = FancyBboxPatch((x0, yc - 0.24), x1 - x0, 0.48,
                           boxstyle="round,pad=0.0,rounding_size=0.07",
                           linewidth=0, facecolor=color, zorder=3)
        ax.add_patch(r)
    # ONE clean allele letter per exon (white on the colored block)
    for px, base in zip(PSV_X, alleles):
        ax.text(px, yc, base, fontsize=9.0, color="white", ha="center",
                va="center", fontweight="bold", zorder=6)
    ax.text(EXONS[-1][1] + 2.0, yc, label, fontsize=10.5, color=color,
            ha="left", va="center", fontweight="bold")


copy_transcript(axR, 1.32, BLUE, "Copy A", ALLELE_A)
copy_transcript(axR, 0.55, ORANGE, "Copy B", ALLELE_B)

axR.text((XMIN + XMAX) / 2, -0.12,
         "PSV-allele vector splits the bucket → two paralog copies",
         fontsize=9.2, color=SLATE, ha="center", va="center")

# ---------------------------------------------------------------------------
# BOTTOM KEY STRIP (spans both panels)
# ---------------------------------------------------------------------------
# "transcript key  =  junction chain  (+)  PSV-allele vector"
ky = 0.72            # main key line
sy = 0.45            # axis sublabels just under it
ny = 0.10            # tiny note, own lower row
axK.text(0.105, ky, "transcript key  =", fontsize=12.5, color=NAVY,
         ha="left", va="center", fontweight="bold")
axK.text(0.305, ky, "junction chain", fontsize=12.5, color=NAVY,
         ha="left", va="center", fontweight="bold")
axK.text(0.379, sy, "structural axis", fontsize=8.2, color=SLATE,
         ha="center", va="center", style="italic")

# "(+)"
axK.text(0.512, ky, "(+)", fontsize=12.5, color=ORANGE, ha="center",
         va="center", fontweight="bold")

axK.text(0.553, ky, "PSV-allele vector", fontsize=12.5, color=ORANGE,
         ha="left", va="center", fontweight="bold")
axK.text(0.645, sy, "copy axis", fontsize=8.2, color=SLATE,
         ha="center", va="center", style="italic")

# tiny note (own lower row, centered, clear of all sublabels)
axK.text(0.5, ny, "PSV = base difference consistent across reads "
         "(not sequencing error)", fontsize=8.0, color=SLATE,
         ha="center", va="center", style="italic")

# thin separator above the key strip
axK.plot([0.0, 1.0], [1.04, 1.04], color=SLATE, lw=0.7,
         transform=axK.transAxes, clip_on=False)

# ---------------------------------------------------------------------------
# Title
# ---------------------------------------------------------------------------
fig.suptitle("Using PSVs: same junction chain, split by allele into copies",
             fontsize=15.5, fontweight="bold", color=NAVY, y=0.975)

out = "/mnt/c/Users/jfris/Desktop/Rustle/bench/readcoherence_fig/fig7_psv_bridge.png"
fig.savefig(out, dpi=190, bbox_inches="tight", facecolor="white")
print("saved", out)
