#!/usr/bin/env python3
"""HMM forward trellis + boundary threading.

Two panels.

LEFT — the inner profile-HMM trellis. Each profile column j has three states
(M_j match, I_j insert, D_j delete). The forward algorithm fills a 2-D grid
(profile column × read position) using log-space transitions. This is the
SAME forward algorithm as a textbook profile HMM (Krogh '94 / HMMER), unchanged.

RIGHT — boundary threading. The graph extension does NOT modify the trellis.
It just chains per-node trellises: node k's exit-boundary (a vector of log-prob
over read positions) becomes node k+1's entry-boundary. Each node's trellis is
still the standard forward DP. The graph adds nothing magic to the math —
only the chain rule for log probabilities along a path.
"""
import matplotlib.pyplot as plt
import matplotlib.patches as patches
from matplotlib.patches import FancyArrowPatch, FancyBboxPatch
import numpy as np

fig, axes = plt.subplots(1, 2, figsize=(15.5, 7.6))
fig.suptitle("Yes — every per-copy profile is a textbook profile HMM (M/I/D forward trellis). "
             "The family graph chains them via boundary threading.",
             fontsize=14, fontweight='bold', y=1.00)

# ── LEFT panel: profile HMM trellis ─────────────────────────────────────────
ax = axes[0]
ax.set_xlim(0, 9)
ax.set_ylim(-0.5, 7.5)
ax.axis('off')
ax.set_title("Profile HMM forward trellis  (one node)", fontsize=12.5, fontweight='bold')

# Grid: profile columns j = 1..4 (x-axis), read positions i = 0..3 (y-axis).
# At each (j, i) we have three states stacked: M, I, D.
n_cols = 4
n_pos = 4
col_x = np.linspace(1.4, 7.5, n_cols + 1)   # extra slot for entry column j=0
pos_y = np.linspace(2.0, 6.2, n_pos)

# Draw axis labels
ax.text(0.0, 6.6, "read pos i", fontsize=10, color='#444')
ax.text(7.8, 0.4, "profile column j", fontsize=10, color='#444', ha='right')

# Column headers: "j=0 (entry)", "j=1", ...
for jx, x in enumerate(col_x):
    if jx == 0:
        ax.text(x, 6.85, "entry", fontsize=9.5, ha='center', color='#666', fontstyle='italic')
    else:
        ax.text(x, 6.85, f"j={jx}", fontsize=10, ha='center', fontweight='bold')

# Row labels at left
for ix, y in enumerate(pos_y):
    ax.text(0.6, y, f"i={ix}", fontsize=10, ha='right', va='center', color='#444')

# Draw each cell as 3 small dots (M, I, D) stacked vertically inside the cell
state_colors = {'M': '#1f77b4', 'I': '#ff7f0e', 'D': '#888888'}
state_dy = 0.18
for jx, x in enumerate(col_x):
    for ix, y in enumerate(pos_y):
        # column 0 (entry) only has the M boundary — show as single dot
        if jx == 0:
            ax.plot(x, y, 'o', color=state_colors['M'], markersize=10, zorder=3)
        else:
            for k, s in enumerate(['M', 'I', 'D']):
                ax.plot(x, y + state_dy * (1 - k), 'o',
                        color=state_colors[s], markersize=8, zorder=3)

# Annotate one cell with the state labels
ann_x = col_x[2]; ann_y = pos_y[2]
ax.text(ann_x + 0.18, ann_y + state_dy + 0.05, " M_j (match)",
        fontsize=8.5, color=state_colors['M'], va='center')
ax.text(ann_x + 0.18, ann_y + 0.05, " I_j (insert)",
        fontsize=8.5, color=state_colors['I'], va='center')
ax.text(ann_x + 0.18, ann_y - state_dy + 0.05, " D_j (delete)",
        fontsize=8.5, color=state_colors['D'], va='center')

# Show three transitions from cell (j-1, i-1) into cell (j, i): M→M, I→M, D→M.
# These are the three terms in the M_j(i) recurrence.
src_x = col_x[1]; src_y = pos_y[1]
dst_x = col_x[2]; dst_y = pos_y[2]
for k, (s, lbl) in enumerate(zip(['M','I','D'], ['a_MM', 'a_IM', 'a_DM'])):
    sy = src_y + state_dy * (1 - k)
    arrow = FancyArrowPatch((src_x + 0.10, sy), (dst_x - 0.10, dst_y + state_dy),
                            arrowstyle='->', color=state_colors[s], lw=1.0,
                            connectionstyle="arc3,rad=0.10", alpha=0.85)
    ax.add_patch(arrow)
    # transition label near mid-arrow
    mx, my = (src_x + dst_x)/2, (sy + dst_y + state_dy)/2 + (0.10 - 0.05*k)
    ax.text(mx, my, lbl, fontsize=7.5, color=state_colors[s], ha='center')

# Recurrence equation at bottom (plain ASCII to avoid matplotlib mathtext quirks).
recurrence = (
    "log M_j(i)  =  e_M(read[i] | j)  +  logsumexp(\n"
    "                  M_{j-1}(i-1) + log a_MM,    "
    "I_{j-1}(i-1) + log a_IM,    "
    "D_{j-1}(i-1) + log a_DM )"
)
ax.text(4.5, 1.05, recurrence, fontsize=9, ha='center', color='#222',
        family='monospace',
        bbox=dict(boxstyle='round,pad=0.4', facecolor='#f4f6fa', edgecolor='#bbb'))

# What is e_M ? — point to slide 4 (POA tilt)
ax.text(4.5, 0.05,
        "e_M(b | j) is the per-copy emission distribution from slide 4\n"
        "(α-blend of singleton and POA-MSA column).",
        fontsize=8.5, ha='center', va='center', color='#555', fontstyle='italic')

# Top-right: log-likelihood of the read against THIS node
ax.text(8.2, 5.5, "Final answer\n(this node):",
        fontsize=9.5, fontweight='bold', color='#222', ha='left')
ax.text(8.2, 4.2,
        "α_exit[i] =\n  logsumexp(\n    M_M(i),\n    I_M(i),\n    D_M(i)\n  )",
        fontsize=8.5, ha='left', va='top', color='#222', family='monospace')

# ── RIGHT panel: boundary threading ─────────────────────────────────────────
ax = axes[1]
ax.set_xlim(0, 10)
ax.set_ylim(0, 7.5)
ax.axis('off')
ax.set_title("Boundary threading  (graph extension — no new math)",
             fontsize=12.5, fontweight='bold')

# Three node trellises in a row, each rendered as a labeled box.
node_x = [0.7, 4.0, 7.3]
node_w = 2.4
node_h = 3.6
node_y = 2.0

for k, x in enumerate(node_x):
    box = FancyBboxPatch((x, node_y), node_w, node_h,
                         boxstyle="round,pad=0.05",
                         facecolor='#eef6ff', edgecolor='#5a82bc', lw=1.2)
    ax.add_patch(box)
    ax.text(x + node_w/2, node_y + node_h - 0.35,
            f"node {k+1}\nprofile HMM trellis",
            fontsize=10.5, ha='center', va='center', fontweight='bold', color='#1f3a6a')
    # mini-trellis dots inside
    for jj in range(4):
        for ii in range(3):
            ax.plot(x + 0.4 + jj*0.50, node_y + 0.6 + ii*0.50,
                    'o', color='#5a82bc', markersize=4, alpha=0.55)

# Boundary vectors between nodes
def boundary_arrow(x_start, x_end, label, color='#a02060'):
    arr = FancyArrowPatch((x_start, node_y + node_h/2),
                          (x_end, node_y + node_h/2),
                          arrowstyle='->', color=color, lw=1.8,
                          mutation_scale=15)
    ax.add_patch(arr)
    ax.text((x_start + x_end)/2, node_y + node_h/2 + 0.45,
            label, fontsize=9, ha='center', color=color, fontweight='bold')

boundary_arrow(node_x[0] + node_w, node_x[1], "α¹_exit[r]")
boundary_arrow(node_x[1] + node_w, node_x[2], "α²_exit[r]")

# Entry to node 1
ax.text(node_x[0] - 0.05, node_y + node_h/2,
        "α¹_entry[r]\n= [0, -∞, …, -∞]",
        fontsize=9, ha='right', va='center', color='#a02060', family='monospace')

# Exit from last node = read log-likelihood
ax.annotate("log P(read | path)\n  = αᴷ_exit[L]",
            xy=(node_x[2] + node_w + 0.05, node_y + node_h/2),
            xytext=(node_x[2] + node_w + 0.3, node_y + node_h/2),
            fontsize=9.5, ha='left', va='center', color='#a02060',
            fontweight='bold', family='monospace')

# The key statement of the slide
ax.text(5, 6.9, "Each node runs the SAME forward trellis (left).",
        fontsize=11, ha='center', fontweight='bold', color='#1f3a6a')
ax.text(5, 6.45,
        "The graph just chains them: node k's per-read-position exit log-prob\n"
        "becomes node k+1's entry log-prob. No new equations.",
        fontsize=10, ha='center', color='#444')

# Per-paralog version
ax.text(5, 1.4,
        "For paralog c: walk c's path through the family graph,\n"
        "use c's per-copy profile at each node (slide 4), get  log P(read | paralog c).",
        fontsize=10, ha='center', color='#225522',
        bbox=dict(boxstyle='round,pad=0.35', facecolor='#e8f5e8', edgecolor='#207020', lw=1.0))

# Bottom: code reference
ax.text(5, 0.4,
        "Code: profile_forward_with_boundary_banded  (per-node trellis)  +  "
        "forward_against_path_for_copy  (boundary threading along c's path).\n"
        "src/rustle/vg_hmm/scorer.rs",
        fontsize=8.5, ha='center', color='#666', fontstyle='italic')

plt.tight_layout(rect=[0, 0, 1, 0.96])
plt.savefig('/scratch/jxi21/Assembler/Rustle/bench/meeting_2026_05_06_hmm_em/05_hmm_trellis.png',
            dpi=160, bbox_inches='tight')
print("05_hmm_trellis.png written")
