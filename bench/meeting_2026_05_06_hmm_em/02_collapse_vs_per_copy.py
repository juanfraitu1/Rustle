#!/usr/bin/env python3
"""Side-by-side: averaged consensus profile (info collapse, BAD) vs per-copy
profiles (each copy keeps its own emission, GOOD). Directly addresses the
advisor's "the VG just averages everything" worry.

Layout uses a vertical 5-band stack on each panel so labels never overlap:
  [1] paralog sequences        (the input MSA)
  [2] arrow + transformation label
  [3] emission-probability bar charts P(base | profile, column)
  [4] A/C/G/T axis labels      (so the "bands" are obvious)
  [5] read-base callout + conclusion
"""
import matplotlib.pyplot as plt
import matplotlib.patches as patches
from matplotlib.patches import FancyBboxPatch
import numpy as np

fig, axes = plt.subplots(1, 2, figsize=(15.5, 8.6))
fig.suptitle("Why a single consensus profile is wrong for paralog assignment",
             fontsize=17, fontweight='bold', y=0.99)

# Toy MSA: 4 columns × 3 paralogs.
# Column 2 is the paralog-distinguishing SNP (A=C, B=G, C=G).
copies = ['A', 'B', 'C']
seqs = [
    ['A', 'C', 'T', 'A'],   # paralog A
    ['A', 'G', 'T', 'A'],   # paralog B
    ['A', 'G', 'T', 'T'],   # paralog C
]
copy_colors = ['#1f77b4', '#ff7f0e', '#2ca02c']
base_to_idx = {'A': 0, 'C': 1, 'G': 2, 'T': 3}
base_colors = {'A': '#1f77b4', 'C': '#ff7f0e', 'G': '#2ca02c', 'T': '#9467bd'}

n_cols = 4
col_x = np.array([1.4, 3.4, 5.4, 7.4])  # x of each profile column's center

def emission_consensus(col):
    counts = np.zeros(4)
    for s in seqs:
        counts[base_to_idx[s[col]]] += 1
    return counts / counts.sum()

def emission_singleton(seq, col):
    eps = 0.02
    e = np.full(4, eps / 3.0)
    e[base_to_idx[seq[col]]] = 1.0 - eps
    return e

# Vertical y positions (one definition, used by both panels for symmetry).
Y = {
    "msa_top":   9.6,   # paralog letters
    "msa_lbl":   9.95,  # "MSA columns" header
    "arrow_top": 8.4,
    "arrow_bot": 7.7,
    "bar_base":  4.6,   # bottom of bar chart
    "bar_top":   7.4,   # top of bar chart
    "axis_lbl":  4.30,  # "A C G T" axis labels under bars
    "col_lbl":   3.80,  # "col 1, col 2..." labels
    "read_box":  2.55,  # read-base callout
    "concl":     1.20,  # conclusion text
    "icon":      0.45,
}

def draw_panel(ax, mode):
    """mode = 'consensus' or 'per_copy'."""
    ax.set_xlim(0.0, 9.5)
    ax.set_ylim(0.0, 10.5)
    ax.axis('off')
    if mode == 'consensus':
        ax.set_title("Collapsed consensus  (one shared profile per node)",
                     fontsize=13, fontweight='bold', color='#a02020', pad=6)
    else:
        ax.set_title("Per-copy profiles  (one profile per paralog, NOT averaged)",
                     fontsize=13, fontweight='bold', color='#207020', pad=6)

    # ── BAND 1: input MSA ────────────────────────────────────────────────────
    ax.text(0.45, Y["msa_lbl"], "MSA columns:",
            fontsize=9.5, color='#666', va='center', fontstyle='italic')
    for i, (s, c) in enumerate(zip(seqs, copy_colors)):
        y = Y["msa_top"] - i * 0.42
        ax.text(0.45, y, f"copy {copies[i]}",
                fontsize=10.5, color=c, fontweight='bold', va='center')
        for j, b in enumerate(s):
            ax.text(col_x[j], y, b, fontsize=13, ha='center', va='center',
                    color=c, fontweight='bold')

    # ── BAND 2: transformation arrow ─────────────────────────────────────────
    if mode == 'consensus':
        msg = "average across copies"
    else:
        msg = "keep each copy's own base"
    ax.annotate(msg, xy=((col_x[0] + col_x[-1]) / 2, Y["arrow_bot"]),
                xytext=((col_x[0] + col_x[-1]) / 2, Y["arrow_top"]),
                fontsize=10, ha='center', color='#444',
                arrowprops=dict(arrowstyle='->', color='#444', lw=1.2))

    # ── BAND 3: emission-probability bar charts ──────────────────────────────
    bar_h_max = Y["bar_top"] - Y["bar_base"]      # max height = probability 1.0
    # y-axis hint
    ax.plot([0.55, 0.55], [Y["bar_base"], Y["bar_top"]], color='#aaa', lw=0.6)
    ax.text(0.20, (Y["bar_base"] + Y["bar_top"]) / 2,
            "P(base | profile, col)",
            fontsize=8.5, color='#666', va='center', rotation=90)
    # tick marks at 0 and 1
    for v, lbl in [(0.0, "0"), (1.0, "1")]:
        yy = Y["bar_base"] + v * bar_h_max
        ax.plot([0.50, 0.55], [yy, yy], color='#aaa', lw=0.6)
        ax.text(0.40, yy, lbl, fontsize=8, color='#666', ha='right', va='center')

    if mode == 'consensus':
        # ONE distribution per column (shared across paralogs).
        for j in range(n_cols):
            e = emission_consensus(j)
            for k, b in enumerate('ACGT'):
                bar_x = col_x[j] - 0.45 + 0.30 * k
                ax.bar(bar_x, e[k] * bar_h_max, width=0.22,
                       bottom=Y["bar_base"], color=base_colors[b],
                       edgecolor='black', lw=0.4)
                # value label above bar
                if e[k] > 0.05:
                    ax.text(bar_x, Y["bar_base"] + e[k] * bar_h_max + 0.10,
                            f"{e[k]:.2f}", fontsize=7.5, ha='center', color='#333')
                # base letter under bar
                ax.text(bar_x, Y["axis_lbl"], b, fontsize=9, ha='center',
                        color=base_colors[b], fontweight='bold')
    else:
        # THREE distributions per column (one per copy), interleaved.
        copy_offsets = [-0.55, 0.0, 0.55]
        for j in range(n_cols):
            for ci, (s, col) in enumerate(zip(seqs, copy_colors)):
                e = emission_singleton(s, j)
                center = col_x[j] + copy_offsets[ci]
                for k, b in enumerate('ACGT'):
                    bar_x = center - 0.21 + 0.14 * k
                    ax.bar(bar_x, e[k] * bar_h_max, width=0.10,
                           bottom=Y["bar_base"], color=base_colors[b],
                           edgecolor='black', lw=0.3)
                # copy label above its mini-distribution
                ax.text(center, Y["bar_top"] + 0.18, copies[ci],
                        fontsize=8.5, ha='center', color=col, fontweight='bold')
                # ACGT axis under each mini-distribution
                for k, b in enumerate('ACGT'):
                    bar_x = center - 0.21 + 0.14 * k
                    ax.text(bar_x, Y["axis_lbl"], b, fontsize=6.5, ha='center',
                            color=base_colors[b])

    # column labels at the bottom of the bar band
    for j in range(n_cols):
        ax.text(col_x[j], Y["col_lbl"], f"col {j+1}",
                fontsize=9.5, ha='center', color='#444', fontweight='bold')

    # Highlight column 2 with a soft band (the SNP column) — drawn behind bars.
    snp_band = patches.Rectangle(
        (col_x[1] - 0.75, Y["bar_base"] - 0.05),
        1.5, bar_h_max + 0.2,
        facecolor='#fff3cd', alpha=0.40, edgecolor='#caa845', lw=0.8, zorder=0,
    )
    ax.add_patch(snp_band)
    if mode == 'consensus':
        # Place the SNP-column callout to the RIGHT, off the bars.
        ax.text(col_x[3] + 0.6, Y["bar_top"] - 0.4,
                "↑ paralog-distinguishing\n   SNP column",
                fontsize=8.5, ha='left', va='top', color='#a04040', fontstyle='italic')
    else:
        ax.text(col_x[3] + 0.6, Y["bar_top"] - 0.4,
                "↑ paralog-distinguishing\n   SNP column",
                fontsize=8.5, ha='left', va='top', color='#a04040', fontstyle='italic')

    # ── BAND 5: read-base callout + conclusion ───────────────────────────────
    ax.text(col_x[1], Y["read_box"],
            "Now: a read aligns with base 'G' at column 2.  What does the profile say?",
            fontsize=10, ha='center', va='center', color='#222',
            bbox=dict(boxstyle='round,pad=0.30', facecolor='#fff', edgecolor='#888'))

    if mode == 'consensus':
        concl = ("P(G | consensus, col 2)  ≈  2/3   for ALL paralogs\n"
                 "→ no paralog-distinguishing signal  →  EM cannot tell A from B or C")
        bg, ec, color = '#fde0e0', '#a02020', '#a02020'
        icon, icon_color = '✗', '#a02020'
    else:
        concl = ("P(G | A, col 2) = 0.02     P(G | B, col 2) = 0.98     P(G | C, col 2) = 0.98\n"
                 "→ A is decisively wrong; B and C are candidates  →  EM has signal to redistribute")
        bg, ec, color = '#e0f5e0', '#207020', '#207020'
        icon, icon_color = '✓', '#207020'

    ax.text((col_x[0] + col_x[-1]) / 2, Y["concl"], concl,
            fontsize=10, ha='center', va='center', color=color,
            bbox=dict(boxstyle='round,pad=0.40', facecolor=bg, edgecolor=ec))
    ax.text(0.65, Y["icon"], icon, fontsize=42, fontweight='bold', color=icon_color)


draw_panel(axes[0], 'consensus')
draw_panel(axes[1], 'per_copy')

# Footer: explain what a "band" / "bar group" is.
fig.text(0.5, -0.005,
         "Each bar group at a profile column is the emission distribution P(base | column) — one bar per nucleotide (A C G T). "
         "Heights sum to 1.\nLeft: a single distribution per column (the consensus profile, paralogs collapsed).  "
         "Right: one distribution per paralog (the per-copy profiles).",
         fontsize=10, ha='center', color='#333')

plt.tight_layout(rect=[0, 0.02, 1, 0.96])
plt.savefig('/scratch/jxi21/Assembler/Rustle/bench/meeting_2026_05_06_hmm_em/02_collapse_vs_per_copy.png',
            dpi=160, bbox_inches='tight')
print("02_collapse_vs_per_copy.png written")
