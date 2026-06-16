import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch

# ---- Palette (shared style guide) ----
NAVY = "#22313f"
GRAY = "#8898aa"
ORANGE = "#e8590c"
GREEN = "#1e8449"

CELL_FILL = "#f6f8fa"
CELL_EDGE = "#c4ccd4"
ORANGE_FILL = "#fdebd9"

plt.rcParams["font.family"] = "DejaVu Sans"

fig, ax = plt.subplots(figsize=(12.6, 7.2))
ax.set_xlim(0, 126)
ax.set_ylim(0, 72)
ax.axis("off")

# ===================== Grid geometry =====================
# 3 columns (identity resolved), 2 rows (task)
COL_X = [38.5, 70.5, 102.5]          # centers of the 3 data columns
COL_LABELS = ["— (none)", "allele / haplotype", "paralog copy"]
ROW_Y = [44.0, 22.0]                 # centers: top = quantify, bottom = de-novo
ROW_LABELS = ["quantify\nknown transcripts", "de-novo\nisoform assembly"]

CW = 29.0   # cell width
CH = 16.5   # cell height
HDR_Y = 60.0          # column header row y
ROWHDR_X = 12.5       # row header x

def cell(cx, cy, w, h, lines, fill, edge, lw=1.6, rounding=1.2):
    p = FancyBboxPatch((cx - w / 2.0, cy - h / 2.0), w, h,
                       boxstyle=f"round,pad=0.02,rounding_size={rounding}",
                       linewidth=lw, edgecolor=edge, facecolor=fill, zorder=3)
    ax.add_patch(p)
    return p

# ===================== Title =====================
ax.text(63, 69.2, "Thesis contribution: the open cell",
        ha="center", va="center", color=NAVY, fontsize=17, fontweight="bold")

# ===================== Axis super-labels =====================
# Top axis caption (columns)
ax.text((COL_X[0] + COL_X[-1]) / 2.0, 65.4, "identity resolved",
        ha="center", va="center", color=GRAY, fontsize=10.5,
        fontweight="bold", style="italic")
# Left axis caption (rows) -- rotated
ax.text(3.4, (ROW_Y[0] + ROW_Y[1]) / 2.0, "task", rotation=90,
        ha="center", va="center", color=GRAY, fontsize=10.5,
        fontweight="bold", style="italic")

# ===================== Column headers =====================
for x, lab in zip(COL_X, COL_LABELS):
    ax.text(x, HDR_Y, lab, ha="center", va="center", color=NAVY,
            fontsize=11.5, fontweight="bold")

# ===================== Row headers =====================
for y, lab in zip(ROW_Y, ROW_LABELS):
    ax.text(ROWHDR_X, y, lab, ha="center", va="center", color=NAVY,
            fontsize=11.5, fontweight="bold", linespacing=1.3)

# ===================== Cell contents =====================
# format: (row_idx, col_idx) -> dict
# row 0 = quantify, row 1 = de-novo

def draw_cell(ri, ci, lines, tcolor=NAVY, fill=CELL_FILL, edge=CELL_EDGE,
              fs=10.0, weight="normal", lw=1.6, sub=None, sub_color=GRAY):
    cx = COL_X[ci]
    cy = ROW_Y[ri]
    cell(cx, cy, CW, CH, lines, fill, edge, lw=lw)
    if sub is not None:
        # main text nudged up, footnote below
        ax.text(cx, cy + 1.6, lines, ha="center", va="center", color=tcolor,
                fontsize=fs, fontweight=weight, linespacing=1.3, zorder=4)
        ax.text(cx, cy - 4.6, sub, ha="center", va="center", color=sub_color,
                fontsize=7.6, style="italic", linespacing=1.2, zorder=4)
    else:
        ax.text(cx, cy, lines, ha="center", va="center", color=tcolor,
                fontsize=fs, fontweight=weight, linespacing=1.35, zorder=4)

# --- top row: quantify ---
# (quantify, none) -- muted slate
draw_cell(0, 0, "Salmon / kallisto\n(short-read)",
          tcolor=GRAY, fs=10.0)
# (quantify, allele)
draw_cell(0, 1, "RPVG · pantas\n(spliced pangenome)",
          tcolor=NAVY, fs=10.0)
# (quantify, paralog copy) -- Paraphase with footnote
draw_cell(0, 2, "Paraphase*",
          tcolor=NAVY, fs=11.5, weight="bold",
          sub="* DNA / gene level,\nnot isoforms")

# --- bottom row: de-novo ---
# (de-novo, none)
draw_cell(1, 0, "IsoQuant · Cupcake\nFLAIR · LRAA",
          tcolor=NAVY, fs=9.6)
# (de-novo, allele)
draw_cell(1, 1, "HapIso\n(diploid)",
          tcolor=NAVY, fs=10.5)
# (de-novo, paralog copy) -- THE highlighted cell
cx = COL_X[2]
cy = ROW_Y[1]
cell(cx, cy, CW + 0.6, CH + 0.6, "", ORANGE_FILL, ORANGE, lw=3.0, rounding=1.4)
ax.text(cx, cy + 2.4, "★ THIS THESIS", ha="center", va="center",
        color=ORANGE, fontsize=13.5, fontweight="bold", zorder=5)
ax.text(cx, cy - 4.0, "the open cell", ha="center", va="center",
        color=ORANGE, fontsize=8.4, style="italic", zorder=5)

# ===================== BOTTOM STRIP =====================
STRIP_TOP = 11.4
# thin separator rule
ax.plot([7, 119], [STRIP_TOP, STRIP_TOP], color=CELL_EDGE, lw=1.1, zorder=1)

# object line
ax.text(7, 7.6,
        "object:",
        ha="left", va="center", color=NAVY, fontsize=9.6, fontweight="bold")
ax.text(20.5, 7.6,
        "spliced family variation graph: a transcript = a path over junctions (+) PSV alleles;  "
        "a copy = an allele-consistent haplotype + its isoforms",
        ha="left", va="center", color=NAVY, fontsize=9.4)

# theorem line
ax.text(7, 3.0,
        "theorem:",
        ha="left", va="center", color=GREEN, fontsize=9.6, fontweight="bold")
ax.text(20.5, 3.0,
        "identifiable iff copies differ at ≥ K PSV columns spanned by enough reads "
        "(else provably merged)",
        ha="left", va="center", color=NAVY, fontsize=9.4)

plt.savefig("/mnt/c/Users/jfris/Desktop/Rustle/bench/readcoherence_fig/fig9_contribution.png",
            dpi=190, bbox_inches="tight", facecolor="white")
print("done")
