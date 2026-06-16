#!/usr/bin/env python3
"""IGV-style locus panel: read-coherence recovers an isoform StringTie skips."""
import json
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle

# ---- palette (shared style guide) ----
C_REF = "#22313f"      # RefSeq / annotation
C_ST = "#8898aa"       # StringTie / flow
C_RCG = "#e8590c"      # read-coherence (accent)
C_READ = "#cfd8dc"     # plain read alignments
C_SPAN = "#f0932b"     # spanning alignments (evidence)

DATA = "/tmp/cre_guided/figdata/igv_locus.json"
OUT = "/mnt/c/Users/jfris/Desktop/Rustle/bench/readcoherence_fig/fig1_igv_locus.png"

d = json.load(open(DATA))
W0, W1 = d["window"]
SPAN = W1 - W0


def clamp(x):
    return max(W0, min(W1, x))


def exon_box(ax, y, s, e, h, color, z=3):
    s, e = clamp(s), clamp(e)
    if e <= s:
        return
    ax.add_patch(Rectangle((s, y - h / 2), e - s, h,
                           linewidth=0, facecolor=color, zorder=z))


def intron_line(ax, y, xs, xe, color, lw=1.2, z=2):
    xs, xe = clamp(xs), clamp(xe)
    if xe > xs:
        ax.plot([xs, xe], [y, y], color=color, lw=lw, zorder=z,
                solid_capstyle="butt")


def draw_transcript(ax, y, exons, color, h, line_lw=1.4):
    xs = min(e[0] for e in exons)
    xe = max(e[1] for e in exons)
    intron_line(ax, y, xs, xe, color, lw=line_lw)
    for (a, b) in exons:
        exon_box(ax, y, a, b, h, color)


def draw_read(ax, y, read, color, h):
    blocks = read["blocks"]
    xs = clamp(min(b[0] for b in blocks))
    xe = clamp(max(b[1] for b in blocks))
    intron_line(ax, y, xs, xe, color, lw=0.7, z=2)
    for (a, b) in blocks:
        a, b = clamp(a), clamp(b)
        if b > a:
            ax.add_patch(Rectangle((a, y - h / 2), b - a, h,
                                   linewidth=0, facecolor=color, zorder=3))


# ---- internal exon region StringTie skips (from RefSeq exons) ----
ref_exons = d["ref"]["exons"]
SKIP_L = min(e[0] for e in ref_exons)
SKIP_R = max(e[1] for e in ref_exons)

# ---- layout: rows from top to bottom ----
fig, ax = plt.subplots(figsize=(10.2, 7.4))

ROW = 1.0
y = 0.0
yticks, ylabels = [], []

TH = 0.62   # transcript exon height
RH = 0.40   # read block height

# highlight band over the skipped internal-exon region (drawn first, behind)
band = ax.axvspan(SKIP_L, SKIP_R, color="#fff3e0", zorder=0)

top_y = None

# (1) RefSeq -- single row
y -= ROW
draw_transcript(ax, y, ref_exons, C_REF, TH, line_lw=1.6)
yticks.append(y); ylabels.append("RefSeq")
top_y = y

# (2) StringTie -- 4 variant rows (plus strand)
st = d["stringtie"]
st_names = sorted(st.keys())
st_rows = []  # (y, intron_start, intron_end) for strand chevrons
for name in st_names:
    y -= ROW
    exons = st[name]["exons"]
    draw_transcript(ax, y, exons, C_ST, TH)
    st_rows.append((y, min(e[0] for e in exons), max(e[1] for e in exons)))
    yticks.append(y); ylabels.append(name)

# (3) read-coherence -- single row
rcg_id = d["rcg_id"]
rcg = d["rcg"][rcg_id]
y -= ROW
draw_transcript(ax, y, rcg["exons"], C_RCG, TH, line_lw=1.8)
yticks.append(y); ylabels.append("read-coherence")
rcg_y = y

# gap before reads
y -= ROW * 0.7

# (4) alignments -- spanning first (grouped, amber), then a few capped others
reads = d["reads"]
spans = [r for r in reads if r.get("spans_ref")]


def overlaps_band(r):
    # block-level overlap with the skipped internal-exon region
    return any(b[1] > SKIP_L and b[0] < SKIP_R for b in r["blocks"])


# keep only non-spanning alignments that touch the highlighted band and stay
# near it (drop the stray far-left singleton and very-long far-right reads that
# distract from the amber-vs-StringTie contrast)
others = [r for r in reads if not r.get("spans_ref")
          and overlaps_band(r)
          and r["start"] >= SKIP_L - 0.05 * SPAN
          and r["end"] <= SKIP_R + 0.10 * SPAN]
# sort spanning by start for tidy stacking
spans.sort(key=lambda r: (r["start"], r["end"]))
others.sort(key=lambda r: (r["start"], r["end"]))

RROW = 0.46
reads_top_y = y

# amber spanning alignments (ALL kept)
for r in spans:
    y -= RROW
    draw_read(ax, y, r, C_SPAN, RH)
span_bottom_y = y

# a small separator, then a few capped light-gray others (thinned)
CAP_OTHERS = 4
shown_others = others[:CAP_OTHERS]
y -= RROW * 0.4
for r in shown_others:
    y -= RROW
    draw_read(ax, y, r, C_READ, RH)
reads_bottom_y = y

# single tick for the alignments block
align_mid = (reads_top_y + reads_bottom_y) / 2
yticks.append(align_mid)
ylabels.append("alignments")

# ---- StringTie group label (one bold) ----
st_ys = yticks[1:1 + len(st_names)]
st_mid = (st_ys[0] + st_ys[-1]) / 2

# ---- y axis ----
ax.set_yticks(yticks)
ax.set_yticklabels(ylabels, fontsize=9)
# color the tick labels by track
label_colors = [C_REF] + [C_ST] * len(st_names) + [C_RCG, "#4a4a4a"]
for tlab, col in zip(ax.get_yticklabels(), label_colors):
    tlab.set_color(col)
ax.get_yticklabels()[0].set_fontweight("bold")
# rcg label index
ax.get_yticklabels()[1 + len(st_names)].set_fontweight("bold")

ax.set_ylim(reads_bottom_y - RROW * 0.6, top_y + ROW * 0.7)

# ---- x axis in kb ----
ax.set_xlim(W0, W1)
import numpy as np
xticks = np.linspace(W0, W1, 6)
ax.set_xticks(xticks)
ax.set_xticklabels([f"{x/1000:.1f}" for x in xticks], fontsize=9)
ax.set_xlabel(f"{d['chrom']}  position (kb)", fontsize=10)

# ---- strand chevrons -------------------------------------------------------
def chevrons(y, x0, x1, color, minus, n=6):
    """Draw small strand chevrons along [x0,x1]; minus=left, plus=right."""
    x0, x1 = clamp(x0), clamp(x1)
    if x1 - x0 < 0.02 * SPAN:
        return
    dx = 0.006 * SPAN
    for frac in np.linspace(0.08, 0.92, n):
        xc = x0 + frac * (x1 - x0)
        if minus:   # points left
            tail, head = xc + dx, xc - dx
        else:       # points right
            tail, head = xc - dx, xc + dx
        ax.annotate("", xy=(head, y), xytext=(tail, y),
                    arrowprops=dict(arrowstyle="-|>", color=color, lw=0.9),
                    zorder=4)


# RefSeq + read-coherence are on the minus strand (chevrons point left)
chevrons(top_y, SKIP_L, SKIP_R, "#9fb0bf", minus=True, n=7)
# StringTie variants are on the plus strand (chevrons point right)
for (sy, sx0, sx1) in st_rows:
    chevrons(sy, sx0, sx1, "#5a6b7a", minus=False, n=6)

# ---- spines ----
for sp in ["top", "right", "left"]:
    ax.spines[sp].set_visible(False)
ax.spines["bottom"].set_color("#888888")
ax.tick_params(axis="y", length=0)
ax.tick_params(axis="x", color="#888888")

# ---- minimal annotations ----
# bracket over the skipped region (above RefSeq)
bar_y = top_y + ROW * 0.42
ax.annotate("", xy=(SKIP_L, bar_y), xytext=(SKIP_R, bar_y),
            arrowprops=dict(arrowstyle="-", color=C_RCG, lw=1.6))
ax.text((SKIP_L + SKIP_R) / 2, bar_y + 0.06, "exons StringTie skips",
        ha="center", va="bottom", fontsize=9, color=C_RCG, fontweight="bold")

# spanning-evidence count label, near the amber block
ax.text(W1 - 0.015 * SPAN, (reads_top_y + span_bottom_y) / 2,
        f"{d['n_reads_span_ref']} spanning\nalignments",
        ha="right", va="center", fontsize=8.5, color=C_SPAN, fontweight="bold")

# small StringTie group marker
ax.text(W0 - 0.135 * SPAN, st_mid, "StringTie", rotation=90,
        ha="center", va="center", fontsize=9.5, color=C_ST,
        fontweight="bold")

# ---- title ----
ax.set_title("Read-coherence recovers a read-supported isoform StringTie skips",
             fontsize=12.5, fontweight="bold", color="#1a1a1a", pad=12)

fig.set_facecolor("white")
ax.set_facecolor("white")
plt.tight_layout()
fig.savefig(OUT, dpi=190, bbox_inches="tight", facecolor="white")
print("saved", OUT)
