#!/usr/bin/env python3
"""Problem setup: an FLNC read multi-maps across paralogs; the BAM aligner
hands us 1/NH weights — we want a way to actually decide which copy this
read came from."""
import matplotlib.pyplot as plt
import matplotlib.patches as patches
from matplotlib.patches import FancyArrowPatch, FancyBboxPatch

fig, ax = plt.subplots(figsize=(13, 6.5))
ax.set_xlim(0, 13)
ax.set_ylim(0, 8)
ax.axis('off')

# Title
ax.text(6.5, 7.5, "The assignment problem", fontsize=20, fontweight='bold', ha='center')
ax.text(6.5, 7.05, "An FLNC multi-mapper has the same alignment score against several paralogs — "
                    "which copy did it come from?", fontsize=11, ha='center', color='#444')

# Three paralog tracks
paralog_y = [5.2, 4.0, 2.8]
paralog_labels = ["Paralog A  (e.g. AMY2A)", "Paralog B  (e.g. AMY2B)", "Paralog C  (e.g. LOC101133335)"]
paralog_colors = ['#1f77b4', '#ff7f0e', '#2ca02c']

# Schematic exon layouts for each paralog (very stylized).
exon_layouts = [
    [(1.5, 2.3), (3.0, 3.6), (4.2, 5.0), (5.6, 6.4), (7.1, 7.9)],
    [(1.5, 2.3), (3.0, 3.6), (4.2, 5.0), (5.6, 6.4), (7.1, 7.9)],
    [(1.5, 2.3), (3.0, 3.6), (4.2, 5.0), (5.6, 6.4), (7.1, 7.9)],
]

for y, lbl, col, exons in zip(paralog_y, paralog_labels, paralog_colors, exon_layouts):
    # backbone line
    ax.plot([1.0, 8.5], [y, y], color='#888', lw=0.8, zorder=1)
    # exons
    for (x0, x1) in exons:
        rect = patches.Rectangle((x0, y - 0.18), x1 - x0, 0.36,
                                 facecolor=col, edgecolor='black', lw=0.6, zorder=2)
        ax.add_patch(rect)
    ax.text(0.85, y, lbl, fontsize=10.5, ha='right', va='center', color=col, fontweight='bold')

# A multi-mapping read shown 3 times (once over each paralog) with thin dashed arrows
read_x0, read_x1 = 4.4, 6.3
for y in paralog_y:
    rect = patches.Rectangle((read_x0, y + 0.45), read_x1 - read_x0, 0.18,
                             facecolor='#d62728', edgecolor='black', lw=0.5, alpha=0.85, zorder=3)
    ax.add_patch(rect)
    ax.text(read_x0 + (read_x1 - read_x0)/2, y + 0.54,
            "read R", fontsize=7.5, ha='center', va='center', color='white', fontweight='bold')

# Annotation: "weight = 1/NH"
ax.annotate("BAM aligner gives this read three placements\nwith equal weight (1/NH ≈ 0.33 each)",
            xy=(6.4, 5.5), xytext=(9.0, 5.5),
            fontsize=10.5, va='center',
            arrowprops=dict(arrowstyle='->', color='#444', lw=1.0))
ax.annotate("", xy=(6.4, 4.3), xytext=(9.0, 5.3),
            arrowprops=dict(arrowstyle='->', color='#888', lw=0.6))
ax.annotate("", xy=(6.4, 3.1), xytext=(9.0, 5.3),
            arrowprops=dict(arrowstyle='->', color='#888', lw=0.6))

# Right-side question box
qbox = FancyBboxPatch((9.0, 1.0), 3.6, 1.65,
                     boxstyle="round,pad=0.12",
                     facecolor='#fff3cd', edgecolor='#caa845', lw=1.5)
ax.add_patch(qbox)
ax.text(10.8, 2.35, "?", fontsize=28, fontweight='bold', ha='center', color='#a07000')
ax.text(10.8, 1.55, "Which paralog\nactually expressed this read?",
        fontsize=10.5, ha='center', va='center', color='#5a4a00')

# Bottom: what we'd lose with naive 1/NH
ax.text(0.5, 1.0,
        "If we keep 1/NH:\n"
        "  • all three copies appear equally expressed\n"
        "  • paralog-specific isoforms get diluted below the assembly threshold\n"
        "  • copies that are actually silent get phantom coverage",
        fontsize=10, va='top', color='#444',
        bbox=dict(boxstyle='round,pad=0.4', facecolor='#f8f8f8', edgecolor='#ccc'))

plt.tight_layout()
plt.savefig('/scratch/jxi21/Assembler/Rustle/bench/meeting_2026_05_06_hmm_em/01_problem_setup.png',
            dpi=160, bbox_inches='tight')
print("01_problem_setup.png written")
