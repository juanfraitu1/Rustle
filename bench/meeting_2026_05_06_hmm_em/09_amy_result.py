#!/usr/bin/env python3
"""Headline empirical result: per-paralog recovery on the medium-similarity
demo (AMY + NBPF). Adds an explicit jaccard-band scale at the top so the
"medium-similarity band" is no longer jargon — it's a named region on a
1-D scale of paralog k-mer similarity, with a mechanistic explanation of
what each band means biologically and why each is hard or easy."""
import matplotlib.pyplot as plt
import matplotlib.patches as patches
from matplotlib.patches import FancyBboxPatch
import numpy as np

# Per-paralog recovery: 1 = recovered exact ('=' in gffcompare refmap), 0 = missed.
# Source: bench/medium_similarity_demo/DEMO.md (commit 0ac3475 verified).
rows = [
    # paralog,                 jaccard, ST, EM_default, HMM (gap+norm)
    ("AMY · AMY2B",            0.59, 1, 0, 1),
    ("AMY · LOC101133335",     0.52, 0, 0, 1),
    ("NBPF · LOC101133271",    0.50, 1, 0, 1),
    ("NBPF · LOC115933275",    0.44, 0, 1, 1),
]

fig, ax = plt.subplots(figsize=(16, 8.6))
ax.set_xlim(-0.2, 16.5)
ax.set_ylim(-1.8, 9.0)
ax.axis('off')
ax.set_title("Per-paralog recovery — HMM-EM rescues a medium-similarity paralog (jaccard 0.52)",
             fontsize=15, fontweight='bold', pad=14)

# ── Top: jaccard-band scale ─────────────────────────────────────────────────
# Scale spans x from 0.5 to 15.5, mapping jaccard 0.0 → 1.0.
SCALE_X0, SCALE_X1 = 0.8, 15.5
SCALE_Y_TOP = 8.4
SCALE_Y_BOT = 7.7
def jacc_to_x(j):
    return SCALE_X0 + (SCALE_X1 - SCALE_X0) * j

bands = [
    (0.00, 0.30, "low",            "#dde2eb", "different gene families\n(very few shared k-mers)"),
    (0.30, 0.60, "medium",         "#fff3cd", "paralogs with substantial divergence\n(many SNPs, indels, exon turnover)"),
    (0.60, 0.90, "high",           "#d4edda", "paralogs with mostly conserved sequence\n(SNP density much lower)"),
    (0.90, 1.00, "near-identical", "#cce4f7", "recent duplicates\n(SNPs sparse — VG primitive shines here)"),
]
for j0, j1, name, color, mech in bands:
    x0, x1 = jacc_to_x(j0), jacc_to_x(j1)
    rect = patches.Rectangle((x0, SCALE_Y_BOT), x1 - x0, SCALE_Y_TOP - SCALE_Y_BOT,
                             facecolor=color, edgecolor='#888', lw=0.8, zorder=1)
    ax.add_patch(rect)
    cx = (x0 + x1) / 2
    ax.text(cx, SCALE_Y_TOP - 0.20, f"{name}-similarity",
            fontsize=10, fontweight='bold', ha='center', color='#222')
    ax.text(cx, SCALE_Y_BOT + 0.18, mech,
            fontsize=8, ha='center', va='bottom', color='#444', fontstyle='italic')

# X-axis ticks at 0.0, 0.3, 0.6, 0.9, 1.0
for j in [0.0, 0.3, 0.6, 0.9, 1.0]:
    x = jacc_to_x(j)
    ax.plot([x, x], [SCALE_Y_BOT - 0.10, SCALE_Y_BOT - 0.02], color='#666', lw=0.8)
    ax.text(x, SCALE_Y_BOT - 0.30, f"{j:.1f}", fontsize=9, ha='center', color='#444')
ax.text((SCALE_X0 + SCALE_X1) / 2, SCALE_Y_BOT - 0.65,
        "k-mer Jaccard between a paralog and its nearest sibling   "
        "(higher = more sequence shared at the k-mer level)",
        fontsize=9.5, ha='center', color='#444', fontstyle='italic')

# Tick marks for each paralog on the scale
mark_y = SCALE_Y_TOP + 0.18
for label, jacc, *_ in rows:
    x = jacc_to_x(jacc)
    ax.plot([x, x], [SCALE_Y_TOP, SCALE_Y_TOP + 0.20], color='#a02060', lw=1.6)
    ax.plot(x, SCALE_Y_TOP + 0.30, marker='v', color='#a02060', markersize=8)
ax.text(jacc_to_x(0.51), SCALE_Y_TOP + 0.55,
        "↓ all four paralogs in this demo land in the medium-similarity band",
        fontsize=9, ha='center', color='#a02060', fontweight='bold')

# ── Recovery table ──────────────────────────────────────────────────────────
n_rows = len(rows)
y_positions = np.arange(n_rows)[::-1] * 1.05  # top-to-bottom, with row spacing
methods = ["StringTie", "heuristic EM", "HMM-EM (per-copy + gap rule)"]
method_x = [7.8, 10.0, 12.4]
table_y0 = 1.6  # base y for the lowest row in the table area

for yi, (label, jacc, st, em, hmm) in zip(y_positions, rows):
    y = table_y0 + yi
    ax.text(0.2, y, label, fontsize=11.5, va='center', fontweight='bold')
    ax.text(3.6, y, f"jaccard {jacc:.2f}", fontsize=10.5, va='center', color='#666')
    statuses = [st, em, hmm]
    for x, s in zip(method_x, statuses):
        if s:
            ax.scatter(x, y, s=600, color='#207020', edgecolor='#0a4a0a',
                       lw=1.0, zorder=3)
            ax.text(x, y, "✓", fontsize=18, ha='center', va='center',
                    color='white', fontweight='bold', zorder=4)
        else:
            ax.scatter(x, y, s=600, color='#fdecec', edgecolor='#a02020',
                       lw=1.0, zorder=3)
            ax.text(x, y, "✗", fontsize=16, ha='center', va='center',
                    color='#a02020', fontweight='bold', zorder=4)

# Method headers
for x, m in zip(method_x, methods):
    ax.text(x, table_y0 + (n_rows - 1) * 1.05 + 0.85, m,
            fontsize=11, ha='center', fontweight='bold')

# Highlight box on LOC101133335 row (second from top)
loc_y = table_y0 + y_positions[1]
ax.add_patch(plt.Rectangle((3.4, loc_y - 0.45), 9.4, 0.90,
                           facecolor='#fff3cd', alpha=0.55, edgecolor='#caa845',
                           lw=1.4, zorder=0))
ax.annotate("HMM-EM is the ONLY method\nthat recovers this paralog —\n"
            "at jaccard 0.52, the regime\nadvisor #2 was worried about.",
            xy=(12.4, loc_y), xytext=(13.4, loc_y),
            fontsize=10, va='center', color='#7a5500',
            arrowprops=dict(arrowstyle='->', color='#7a5500', lw=1.3))

# Bottom commentary
ax.text(8.0, -1.4,
        "Per-paralog exact-match recovery (gffcompare class '=') on subset BAMs.  "
        "Source: bench/medium_similarity_demo/DEMO.md.",
        fontsize=9, ha='center', color='#666', fontstyle='italic')

plt.tight_layout()
plt.savefig('/scratch/jxi21/Assembler/Rustle/bench/meeting_2026_05_06_hmm_em/09_amy_result.png',
            dpi=160, bbox_inches='tight')
print("09_amy_result.png written")
