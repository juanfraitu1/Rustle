#!/usr/bin/env python3
"""Side-by-side: averaged consensus profile (info collapse, BAD) vs per-copy
profiles (each copy keeps its own emission, GOOD). This is the diagram that
directly addresses the advisor's "the VG just averages everything" worry."""
import matplotlib.pyplot as plt
import matplotlib.patches as patches
from matplotlib.patches import FancyBboxPatch
import numpy as np

fig, axes = plt.subplots(1, 2, figsize=(14, 7.0))
fig.suptitle("Why a single consensus profile is wrong for paralog assignment",
             fontsize=17, fontweight='bold', y=0.99)

# A small example: 4 columns, 3 paralogs, with 2 paralog-distinguishing positions.
# Column 1: all paralogs have A (conserved)
# Column 2: A=C, B=G, C=G  ← paralog-distinguishing SNP
# Column 3: all paralogs have T (conserved)
# Column 4: A=A, B=A, C=T  ← paralog-distinguishing SNP
copies = ['A', 'B', 'C']
seqs = [
    ['A', 'C', 'T', 'A'],   # paralog A
    ['A', 'G', 'T', 'A'],   # paralog B
    ['A', 'G', 'T', 'T'],   # paralog C
]
copy_colors = ['#1f77b4', '#ff7f0e', '#2ca02c']

n_cols = 4
col_x = np.array([1.0, 2.4, 3.8, 5.2])

base_to_idx = {'A':0, 'C':1, 'G':2, 'T':3}
def emission_consensus(col):
    """Build a 4-vector emission distribution for one MSA column over ACGT,
    averaging across all copies."""
    counts = np.zeros(4)
    for s in seqs:
        counts[base_to_idx[s[col]]] += 1
    counts /= counts.sum()
    return counts

def emission_singleton(seq, col):
    """Singleton emission: one-hot at the copy's actual base, with a small
    pseudocount so logs are finite."""
    eps = 0.02
    e = np.full(4, eps / 3.0)
    e[base_to_idx[seq[col]]] = 1.0 - eps
    return e

# ── Left panel: collapsed consensus ──────────────────────────────────────────
ax = axes[0]
ax.set_xlim(0, 7.0)
ax.set_ylim(0, 8.4)
ax.axis('off')
ax.set_title("Collapsed consensus (a single shared profile)",
             fontsize=13, fontweight='bold', color='#a02020')

# show the 3 paralog seqs at the top
for i, (s, c) in enumerate(zip(seqs, copy_colors)):
    y = 7.5 - i * 0.45
    ax.text(0.4, y, copies[i], fontsize=11, color=c, fontweight='bold', va='center')
    for j, b in enumerate(s):
        ax.text(col_x[j], y, b, fontsize=12, ha='center', va='center', color=c)

# arrow down
ax.annotate("average", xy=(3.0, 5.7), xytext=(3.0, 6.2),
            fontsize=10, ha='center', color='#444',
            arrowprops=dict(arrowstyle='->', color='#444', lw=1.2))

# Show the consensus emission as bar charts per column
bar_y = 4.0
bar_h = 1.2
for j in range(n_cols):
    e = emission_consensus(j)
    for k, b in enumerate('ACGT'):
        ax.bar(col_x[j] - 0.18 + 0.12 * k, e[k] * bar_h, width=0.10, bottom=bar_y,
               color=['#1f77b4','#ff7f0e','#2ca02c','#9467bd'][k])
    ax.text(col_x[j], bar_y - 0.25, f"col {j+1}", fontsize=9, ha='center', color='#666')

# Draw a "read base = G" at column 2
ax.text(col_x[1], 2.2, "read base 'G'\nat column 2",
        fontsize=10.5, ha='center', va='center',
        bbox=dict(boxstyle='round,pad=0.3', facecolor='#fff', edgecolor='#888'))

# Score: P(G | consensus col 2) is the average — uniform-ish, no discrimination.
ax.text(3.5, 1.0,
        "P(G | consensus, col 2) ≈ 2/3 for ALL paralogs\n"
        "→ no paralog-distinguishing signal\n"
        "→ EM cannot tell A apart from B or C",
        fontsize=10.5, ha='center', va='center', color='#a02020',
        bbox=dict(boxstyle='round,pad=0.4', facecolor='#fde0e0', edgecolor='#a02020'))

# Big red X
ax.text(0.6, 0.35, "✗", fontsize=42, fontweight='bold', color='#a02020')

# ── Right panel: per-copy profiles ───────────────────────────────────────────
ax = axes[1]
ax.set_xlim(0, 7.0)
ax.set_ylim(0, 8.4)
ax.axis('off')
ax.set_title("Per-copy profiles (one emission per paralog, NOT averaged)",
             fontsize=13, fontweight='bold', color='#207020')

for i, (s, c) in enumerate(zip(seqs, copy_colors)):
    y = 7.5 - i * 0.45
    ax.text(0.4, y, copies[i], fontsize=11, color=c, fontweight='bold', va='center')
    for j, b in enumerate(s):
        ax.text(col_x[j], y, b, fontsize=12, ha='center', va='center', color=c)

ax.annotate("each copy keeps its own", xy=(3.0, 5.7), xytext=(3.0, 6.2),
            fontsize=10, ha='center', color='#444',
            arrowprops=dict(arrowstyle='->', color='#444', lw=1.2))

# Three stacks of per-copy emissions per column, color-coded by copy.
# We'll show only column 2 in detail (the SNP column) and indicate the others.
# Actually let's show all 4 columns as small triple-stacks.
copy_offset = [-0.18, 0.0, 0.18]  # x offsets for the 3 copies

bar_y = 4.0
bar_h = 1.0
for j in range(n_cols):
    for ci, (s, col) in enumerate(zip(seqs, copy_colors)):
        e = emission_singleton(s, j)
        # show a small bar for the dominant base, color-coded by copy
        dom = base_to_idx[s[j]]
        ax.bar(col_x[j] + copy_offset[ci], 1.0 * (bar_h - 0.1), width=0.13,
               bottom=bar_y + ci * 0.05, color=col, edgecolor='black', lw=0.4)
        ax.text(col_x[j] + copy_offset[ci], bar_y + 0.55,
                s[j], fontsize=9, ha='center', va='center',
                color='white', fontweight='bold')
    ax.text(col_x[j], bar_y - 0.25, f"col {j+1}", fontsize=9, ha='center', color='#666')

# annotate the SNP column
ax.annotate("paralog-distinguishing\nSNP column", xy=(col_x[1], bar_y + 1.2),
            xytext=(col_x[1] + 1.2, bar_y + 2.0),
            fontsize=9.5, ha='left',
            arrowprops=dict(arrowstyle='->', color='#444', lw=0.9))

# Show "read base = G" at column 2
ax.text(col_x[1], 2.2, "read base 'G'\nat column 2",
        fontsize=10.5, ha='center', va='center',
        bbox=dict(boxstyle='round,pad=0.3', facecolor='#fff', edgecolor='#888'))

# Per-copy scores
ax.text(3.5, 1.0,
        "P(G | col 2):  A → 0.02  ·  B → 0.98  ·  C → 0.98\n"
        "→ A is decisively wrong, B and C are candidates\n"
        "→ EM has signal to redistribute weight",
        fontsize=10.5, ha='center', va='center', color='#207020',
        bbox=dict(boxstyle='round,pad=0.4', facecolor='#e0f5e0', edgecolor='#207020'))

ax.text(0.6, 0.35, "✓", fontsize=42, fontweight='bold', color='#207020')

plt.tight_layout(rect=[0, 0, 1, 0.96])
plt.savefig('/scratch/jxi21/Assembler/Rustle/bench/meeting_2026_05_06_hmm_em/02_collapse_vs_per_copy.png',
            dpi=160, bbox_inches='tight')
print("02_collapse_vs_per_copy.png written")
