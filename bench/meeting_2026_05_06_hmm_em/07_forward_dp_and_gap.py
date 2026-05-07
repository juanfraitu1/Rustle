#!/usr/bin/env python3
"""Forward DP per paralog path + score-gap rule.

Shows: a single read R is scored against each paralog's path through the
family graph using the per-copy profile chain. The result is a log-likelihood
PER paralog. EM redistributes weight by the softmax — but the score-gap rule
abstains when the top-vs-runner-up gap is too small ("don't gamble")."""
import matplotlib.pyplot as plt
import matplotlib.patches as patches
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch
import numpy as np

fig, ax = plt.subplots(figsize=(14, 7.6))
ax.set_xlim(0, 14)
ax.set_ylim(0, 8.6)
ax.axis('off')

ax.text(7, 8.2, "E-step: score the read against each paralog path",
        fontsize=17, fontweight='bold', ha='center')
ax.text(7, 7.78,
        "Per-copy profiles + path topology give one log-likelihood per paralog. EM uses these as the posterior; "
        "the score-gap rule abstains when uncertain.",
        fontsize=10.5, ha='center', color='#444')

# Paralog paths (3 shown). Each path is a chain of nodes labeled E1..E5
# rendered as colored circles.
paralogs = ['A', 'B', 'C']
path_y = [6.2, 4.4, 2.6]
copy_colors = ['#1f77b4', '#ff7f0e', '#2ca02c']
paths = [
    [1, 2, 3, 5],         # A: skips E4
    [1, 2, 3, 4, 5],      # B: full
    [1, 3, 4, 5],         # C: skips E2
]
node_x_grid = {1: 2.0, 2: 3.5, 3: 5.0, 4: 6.5, 5: 8.0}
exon_label = {1: 'E1', 2: 'E2', 3: 'E3', 4: 'E4', 5: 'E5'}

for paralog, y, col, pth in zip(paralogs, path_y, copy_colors, paths):
    # path label
    ax.text(0.55, y, f"paralog {paralog}", fontsize=11.5, fontweight='bold',
            color=col, va='center')
    # backbone line
    ax.plot([1.4, 9.0], [y, y], color='#bbb', lw=0.8, zorder=0)
    # nodes
    for k, nid in enumerate(pth):
        cx = node_x_grid[nid]
        circ = patches.Circle((cx, y), 0.32, facecolor=col, alpha=0.85,
                              edgecolor='black', lw=0.7, zorder=2)
        ax.add_patch(circ)
        ax.text(cx, y, exon_label[nid], fontsize=8.5, ha='center', va='center',
                color='white', fontweight='bold')
    # arrow to score
    ax.annotate("", xy=(11.0, y), xytext=(9.0, y),
                arrowprops=dict(arrowstyle='->', color=col, lw=1.4))

# Read R representation between path A and path B
ax.text(2.0, 7.3, "read R", fontsize=10, color='#d62728', fontweight='bold')
read_rect = patches.Rectangle((2.4, 7.20), 4.5, 0.20,
                              facecolor='#d62728', edgecolor='black', lw=0.5)
ax.add_patch(read_rect)
# little annotation arrows from read to each path
for y in path_y:
    arr = FancyArrowPatch((4.5, 7.18), (4.5, y + 0.4),
                          arrowstyle='-', linestyle=':', color='#888', lw=0.6)
    ax.add_patch(arr)

# Forward DP labels at the right side: log-likelihoods
ll = [-118.3, -149.7, -120.9]   # toy values
posteriors = np.exp(np.array(ll) - max(ll))
posteriors /= posteriors.sum()
gap_with_runner_up = max(ll) - sorted(ll)[-2]   # 2.6 nats — small!

for paralog, y, col, val, post in zip(paralogs, path_y, copy_colors, ll, posteriors):
    box = FancyBboxPatch((11.05, y - 0.35), 2.55, 0.7,
                         boxstyle="round,pad=0.05",
                         facecolor='white', edgecolor=col, lw=1.4)
    ax.add_patch(box)
    ax.text(11.18, y + 0.06, f"log P(R | {paralog}) = {val:.1f}",
            fontsize=9.5, color='#333')
    ax.text(11.18, y - 0.20, f"posterior  ≈  {post:.2f}",
            fontsize=9.5, color=col, fontweight='bold')

# Score-gap rule box at bottom
gap_box = FancyBboxPatch((1.0, 0.35), 12.0, 1.55,
                         boxstyle="round,pad=0.10",
                         facecolor='#fff7e0', edgecolor='#caa845', lw=1.4)
ax.add_patch(gap_box)
ax.text(7, 1.65, "Score-gap rule:  redistribute weight ONLY if  best − runner-up  ≥  Δ",
        fontsize=12, fontweight='bold', ha='center', color='#7a5500')
ax.text(7, 1.20,
        f"This read:  best = paralog A ({ll[0]:.1f}),  runner-up = C ({ll[2]:.1f}),  gap = {gap_with_runner_up:.1f} nats",
        fontsize=10.5, ha='center', color='#444')
ax.text(7, 0.75,
        "if gap < Δ → keep BAM aligner's 1/NH (don't gamble);  if gap ≥ Δ → use posterior to redistribute weight.",
        fontsize=10, ha='center', color='#444', fontstyle='italic')
ax.text(7, 0.40,
        "Default Δ = 10 nats (≈ 4 unambiguous SNPs of evidence).  Mechanically prevents EM from regressing on cases ST already gets.",
        fontsize=9.5, ha='center', color='#7a5500')

plt.tight_layout()
plt.savefig('/scratch/jxi21/Assembler/Rustle/bench/meeting_2026_05_06_hmm_em/07_forward_dp_and_gap.png',
            dpi=160, bbox_inches='tight')
print("07_forward_dp_and_gap.png written")
