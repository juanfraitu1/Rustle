import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patheffects as pe
from matplotlib.patches import FancyBboxPatch, Rectangle

# ---- palette ----
NAVY   = "#22313f"   # bases / text
SLATE  = "#8898aa"   # muted / guides
ORANGE = "#e8590c"   # PSV / read-coherence accent
GREEN  = "#1e8449"   # consistent / real
RED    = "#c0392b"   # error

MONO = "DejaVu Sans Mono"
plt.rcParams["font.family"] = "DejaVu Sans"

# white halo so colored letters stay legible on tinted bands
HALO = [pe.withStroke(linewidth=2.2, foreground="white")]

# ===========================================================================
# Figure scaffold: two panels + a bottom strip.
# ===========================================================================
fig = plt.figure(figsize=(12.8, 6.9))

axL = fig.add_axes([0.040, 0.150, 0.450, 0.730])
axR = fig.add_axes([0.530, 0.150, 0.445, 0.730])
axB = fig.add_axes([0.040, 0.012, 0.935, 0.110])   # bottom strip

for ax in (axL, axR, axB):
    ax.set_xticks([]); ax.set_yticks([])
    for s in ax.spines.values():
        s.set_visible(False)

# ===========================================================================
# LEFT PANEL: align the copies (POA) -> homologous columns
# ===========================================================================
axL.set_xlim(0, 14.0)
axL.set_ylim(0, 10.0)

axL.text(0.2, 9.55, "1. align the copies (POA)  →  homologous columns",
         fontsize=12.5, fontweight="bold", color=NAVY, ha="left", va="center")

# A small 3-row, 12-column MSA. Most columns identical; cols 4 and 9 differ.
COPIES = ["Copy A", "Copy B", "Copy C"]
#         col:      0    1    2    3    4*   5    6    7    8    9*   10   11
seq_A = list("G C A T A C G T G T A C".split())
seq_B = list("G C A T G C G T G T A C".split())
seq_C = list("G C A T A C G T G C A C".split())
PSV_COLS = {4, 9}      # columns where fixed bases differ across copies

NCOL = 12
X0 = 2.05              # left edge of the alignment grid
CW = 0.92             # column width
ROW_Y = [7.55, 6.55, 5.55]   # one y per copy row

# orange PSV bands behind the two differing columns (full grid height)
band_top = ROW_Y[0] + 0.50
band_bot = ROW_Y[-1] - 0.50
for c in PSV_COLS:
    cx = X0 + c * CW
    axL.add_patch(Rectangle((cx - CW / 2, band_bot), CW, band_top - band_bot,
                            facecolor=ORANGE, alpha=0.16, edgecolor="none",
                            zorder=1))
    # PSV tick + label on top of the band
    axL.plot([cx, cx], [band_top + 0.02, band_top + 0.24], color=ORANGE,
             lw=1.8, zorder=4, solid_capstyle="round")
    axL.text(cx, band_top + 0.40, "PSV", fontsize=9.0, color=ORANGE,
             ha="center", va="bottom", fontweight="bold", zorder=4)

# row labels
for name, ry in zip(COPIES, ROW_Y):
    axL.text(X0 - 0.55, ry, name, fontsize=10.0, color=NAVY,
             ha="right", va="center", fontweight="bold")

# the base letters, column-aligned
for ri, (seq, ry) in enumerate(zip([seq_A, seq_B, seq_C], ROW_Y)):
    for c in range(NCOL):
        cx = X0 + c * CW
        is_psv = c in PSV_COLS
        col = ORANGE if is_psv else NAVY
        fw = "bold" if is_psv else "normal"
        eff = HALO if is_psv else None
        axL.text(cx, ry, seq[c], fontsize=13.5, color=col, family=MONO,
                 ha="center", va="center", fontweight=fw, zorder=3,
                 path_effects=eff)

# caption under the panel
axL.text(7.0, 4.55,
         "PSV = column where copies' fixed bases differ",
         fontsize=9.6, color=SLATE, ha="center", va="center", style="italic")

# ===========================================================================
# RIGHT PANEL: confirm across reads (PSV vs sequencing error)
# ===========================================================================
axR.set_xlim(0, 14.0)
axR.set_ylim(0, 10.0)

axR.text(0.2, 9.55, "2. confirm across reads (PSV vs sequencing error)",
         fontsize=12.5, fontweight="bold", color=NAVY, ha="left", va="center")

# read geometry
READ_X0 = 1.4         # left edge of read bars
READ_W = 6.6         # read bar width
RH = 0.155            # read bar half-height
FOCAL_X = 5.6        # x of the focal-position base letter (inside the bars)
N_READS = 6

def pileup(ax, y_top, bases, base_colors, dy=0.46):
    """Stack N reads as light bars, each with one focal base letter."""
    ys = [y_top - i * dy for i in range(len(bases))]
    for y in ys:
        ax.add_patch(FancyBboxPatch((READ_X0, y - RH), READ_W, 2 * RH,
                                    boxstyle="round,pad=0.0,rounding_size=0.05",
                                    linewidth=0, facecolor="#eceff1", zorder=2))
    for y, b, c in zip(ys, bases, base_colors):
        ax.text(FOCAL_X, y, b, fontsize=11.5, color=c, family=MONO,
                ha="center", va="center", fontweight="bold", zorder=4)
    return ys

# faint focal-column guide spanning both pileups
axR.plot([FOCAL_X, FOCAL_X], [0.95, 8.55], color=SLATE, lw=0.8,
         linestyle=(0, (3, 3)), zorder=1)

# --- TOP pileup: a PSV column, all 6 reads identical (GREEN) ---
TOP_Y = 8.30
pileup(axR, TOP_Y, ["A"] * N_READS, [GREEN] * N_READS)
axR.text(READ_X0 + READ_W + 0.45, TOP_Y - 2.5 * 0.46 + 0.0,
         "consistent  ~100%", fontsize=10.0, color=GREEN,
         ha="left", va="center", fontweight="bold")
axR.text(READ_X0 + READ_W + 0.45, TOP_Y - 2.5 * 0.46 - 0.42,
         "→  PSV (real)", fontsize=10.0, color=GREEN,
         ha="left", va="center", fontweight="bold")

# divider between the two pileups
axR.plot([0.6, 13.4], [4.95, 4.95], color=SLATE, lw=0.7, zorder=1)

# --- BOTTOM pileup: an error position, 5 consensus + 1 random (RED) ---
BOT_Y = 4.35
bot_bases  = ["C", "C", "T", "C", "C", "C"]   # one stray 'T'
bot_colors = [NAVY, NAVY, RED, NAVY, NAVY, NAVY]
pileup(axR, BOT_Y, bot_bases, bot_colors)
axR.text(READ_X0 + READ_W + 0.45, BOT_Y - 2.5 * 0.46 + 0.0,
         "sporadic  ~0.1–1%", fontsize=10.0, color=RED,
         ha="left", va="center", fontweight="bold")
axR.text(READ_X0 + READ_W + 0.45, BOT_Y - 2.5 * 0.46 - 0.42,
         "→  sequencing error", fontsize=10.0, color=RED,
         ha="left", va="center", fontweight="bold")

# caption under the panel
axR.text(7.0, 0.55,
         "fixed & consistent across every read  vs  random / rare",
         fontsize=9.6, color=SLATE, ha="center", va="center", style="italic")

# ===========================================================================
# BOTTOM STRIP (spans both panels)
# ===========================================================================
axB.set_xlim(0, 1); axB.set_ylim(0, 1)
axB.plot([0.0, 1.0], [0.98, 0.98], color=SLATE, lw=0.7,
         transform=axB.transAxes, clip_on=False)
axB.text(0.5, 0.66,
         "identifiable iff  ≥ K  PSV columns co-segregate"
         "   (K independent error-agreements are improbable)",
         fontsize=11.5, color=NAVY, ha="center", va="center", fontweight="bold")
axB.text(0.5, 0.20,
         "read → copy by which PSV alleles it carries  (+ NM / AS)",
         fontsize=8.6, color=SLATE, ha="center", va="center", style="italic")

# ===========================================================================
# Title
# ===========================================================================
fig.suptitle("Calling a PSV: a fixed inter-copy column, confirmed across reads",
             fontsize=15.0, fontweight="bold", color=NAVY, y=0.985)

out = "/mnt/c/Users/jfris/Desktop/Rustle/bench/readcoherence_fig/fig8_psv_calc.png"
fig.savefig(out, dpi=190, bbox_inches="tight", facecolor="white")
print("saved", out)
