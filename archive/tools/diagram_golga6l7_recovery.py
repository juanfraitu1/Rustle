#!/usr/bin/env python3
"""Visualize the StringTie vs Rustle+VG comparison on GOLGA6L7."""
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle, FancyArrowPatch

fig, axes = plt.subplots(3, 1, figsize=(12, 7), sharex=True)
fig.suptitle("GOLGA6L7 on chr19: three tandem paralogs (~40 kb apart)",
             fontsize=13, y=0.995)

# Real coordinates (simplified to range 0-90 for display)
# L7_1: 104789647-104796276 → 0-6629
# L7_2: 104830536-104837094 → 40889-47447  
# L7_3: 104871356-104877901 → 81709-88254
offset = 104789647
def rel(x): return x - offset

# Reference transcripts (9 exons each)
refs = [
    ('L7_1', [(104789647, 104789918), (104791223, 104791352), (104792167, 104792196),
              (104792466, 104792516), (104792605, 104792698), (104792781, 104792887),
              (104794130, 104794188), (104794518, 104794659), (104794954, 104796276)]),
    ('L7_2', [(104830536, 104830737), (104832042, 104832171), (104832986, 104833015),
              (104833284, 104833335), (104833424, 104833517), (104833600, 104833706),
              (104834953, 104835011), (104835341, 104835482), (104835777, 104837094)]),
    ('L7_3', [(104871356, 104871545), (104872851, 104872979), (104873794, 104873823),
              (104874093, 104874143), (104874232, 104874325), (104874408, 104874514),
              (104875759, 104875817), (104876147, 104876288), (104876583, 104877901)]),
]

# StringTie output (only L7_1 matches)
stringtie = [
    ('L7_1', refs[0][1]),
    ('L7_2', []),
    ('L7_3', []),
]
# Rustle+VG output (all match)
rustle = [
    ('L7_1', refs[0][1]),
    ('L7_2', refs[1][1]),
    ('L7_3', refs[2][1]),
]

def draw_track(ax, transcripts, title, arrow_dir='right'):
    """Draw a set of transcripts as exon boxes."""
    ax.set_xlim(-1000, 90000)
    ax.set_ylim(-0.5, len(transcripts) - 0.5 + 0.5)
    ax.set_title(title, loc='left', fontsize=11, fontweight='bold')
    ax.set_yticks(range(len(transcripts)))
    ax.set_yticklabels([t[0] for t in transcripts])
    for i, (name, exons) in enumerate(transcripts):
        y = i
        if not exons:
            ax.text(5000, y, "— no assembly emitted —", fontsize=9,
                    color='red', va='center', style='italic')
            continue
        # Sort exons, draw
        exs = sorted(exons)
        # connecting line
        ax.plot([rel(exs[0][0]), rel(exs[-1][1])], [y, y], 
                color='gray', linewidth=0.6)
        for es, ee in exs:
            rect = Rectangle((rel(es), y - 0.22), rel(ee) - rel(es), 0.44,
                           facecolor='steelblue', edgecolor='black', linewidth=0.5)
            ax.add_patch(rect)
    ax.set_xticks([])
    for spine in ['top', 'right']:
        ax.spines[spine].set_visible(False)

# Row 1: reference
draw_track(axes[0], refs, "Reference annotation  (ground truth: 3 expressed copies)")

# Row 2: StringTie
draw_track(axes[1], stringtie, "StringTie -L (long-read mode)  →  1/3 recovered")

# Row 3: Rustle+VG
draw_track(axes[2], rustle, "Rustle -L --vg (family rescue + extend)  →  3/3 recovered")

# X axis label on bottom
axes[2].set_xlabel(f"genomic coordinate (relative to {offset:,})", fontsize=10)
axes[2].set_xticks([0, 20000, 40000, 60000, 80000])
axes[2].set_xticklabels(['0', '+20 kb', '+40 kb', '+60 kb', '+80 kb'])

# Legend
from matplotlib.lines import Line2D
legend = [
    Line2D([0],[0], marker='s', color='w', markerfacecolor='steelblue',
           markeredgecolor='black', markersize=10, label='exon (class = exact match)'),
    Line2D([0],[0], marker='_', color='gray', markersize=15, label='intron'),
]
axes[0].legend(handles=legend, loc='upper right', fontsize=9)

plt.tight_layout()
plt.savefig('/storage/group/kdm16/default/jxi21/apes_transcriptome_analysis/clusterer/figures/fig1_golga6l7_recovery.png', dpi=150, bbox_inches='tight')
plt.savefig('/storage/group/kdm16/default/jxi21/apes_transcriptome_analysis/clusterer/figures/fig1_golga6l7_recovery.pdf', bbox_inches='tight')
print("Saved fig1_golga6l7_recovery.{png,pdf}")
