#!/usr/bin/env python3
"""Headline empirical result: per-paralog recovery on the medium-similarity
demo (AMY + NBPF). Shows where HMM-EM uniquely rescues a paralog at moderate
divergence — the regime the advisor was worried about."""
import matplotlib.pyplot as plt
import numpy as np

# Per-paralog recovery: 1 = recovered exact ('=' in gffcompare refmap), 0 = missed.
# Source: bench/medium_similarity_demo/DEMO.md (commit 0ac3475 verified).
rows = [
    # paralog,                jaccard, ST, EM_default, HMM (gap+norm)
    ("AMY · AMY2B",           0.59,    1, 0, 1),
    ("AMY · LOC101133335",    0.52,    0, 0, 1),
    ("NBPF · LOC101133271",   0.50,    1, 0, 1),
    ("NBPF · LOC115933275",   0.44,    0, 1, 1),
]

fig, ax = plt.subplots(figsize=(15, 6.0))

# Reorder rows to put the headline-LOC101133335 second from top so it pops.
n_rows = len(rows)
y_positions = np.arange(n_rows)[::-1]  # top-to-bottom

# Draw labels and method dots — pushed right so jaccard label has room.
methods = ["StringTie", "heuristic EM", "HMM-EM (per-copy + gap rule)"]
method_x = [7.5, 9.5, 11.8]

for y, (label, jacc, st, em, hmm) in zip(y_positions, rows):
    # paralog label and jaccard at the left
    ax.text(0.0, y, label, fontsize=11.5, va='center', fontweight='bold')
    ax.text(3.4, y, f"jaccard {jacc:.2f}", fontsize=10.5, va='center', color='#666')
    # color-coded jaccard band hint (medium-similarity = 30-60%)
    if 0.30 <= jacc <= 0.60:
        ax.text(4.8, y, "(medium-similarity)", fontsize=8.5, va='center',
                color='#a02020', fontstyle='italic')

    # Method status icons
    statuses = [st, em, hmm]
    for x, s, m in zip(method_x, statuses, methods):
        if s:
            ax.scatter(x, y, s=600, color='#207020', edgecolor='#0a4a0a', lw=1.0, zorder=3)
            ax.text(x, y, "✓", fontsize=18, ha='center', va='center', color='white',
                    fontweight='bold', zorder=4)
        else:
            ax.scatter(x, y, s=600, color='#fdecec', edgecolor='#a02020', lw=1.0, zorder=3)
            ax.text(x, y, "✗", fontsize=16, ha='center', va='center', color='#a02020',
                    fontweight='bold', zorder=4)

# Method headers
for x, m in zip(method_x, methods):
    ax.text(x, n_rows - 0.5 + 0.55, m, fontsize=11, ha='center', fontweight='bold')

# Highlight box around LOC101133335 row
loc_y = y_positions[1]
ax.add_patch(plt.Rectangle((3.2, loc_y - 0.42), 9.4, 0.84,
                           facecolor='#fff3cd', alpha=0.55, edgecolor='#caa845',
                           lw=1.4, zorder=0))
ax.annotate("HMM-EM is the ONLY method\nthat recovers this paralog —\nat 0.52 jaccard, the medium-\nsimilarity regime.",
            xy=(11.8, loc_y), xytext=(13.0, loc_y),
            fontsize=10, va='center', color='#7a5500',
            arrowprops=dict(arrowstyle='->', color='#7a5500', lw=1.3))

# Bottom commentary: "what this proves"
ax.text(7.5, -1.1,
        "Per-paralog exact-match recovery (gffcompare class '=') on subset BAMs.  Source: bench/medium_similarity_demo/DEMO.md.",
        fontsize=9, ha='center', color='#666', fontstyle='italic')

ax.set_xlim(-0.2, 16.0)
ax.set_ylim(-1.5, n_rows + 0.5)
ax.axis('off')
ax.set_title("Per-paralog recovery — HMM-EM rescues a medium-similarity paralog (jaccard 0.52)",
             fontsize=15, fontweight='bold', pad=20)

plt.tight_layout()
plt.savefig('/scratch/jxi21/Assembler/Rustle/bench/meeting_2026_05_06_hmm_em/07_amy_result.png',
            dpi=160, bbox_inches='tight')
print("07_amy_result.png written")
