#!/usr/bin/env python3
"""Family graph with exon-class nodes, where each node carries a per-copy
profile bundle. Emphasizes: SHARED structure (which exons are orthologous,
what junctions exist) but UN-SHARED sequence (each copy has its own
profile inside the node)."""
import matplotlib.pyplot as plt
import matplotlib.patches as patches
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch

fig, ax = plt.subplots(figsize=(14, 7.5))
ax.set_xlim(0, 14)
ax.set_ylim(0, 8.5)
ax.axis('off')

ax.text(7, 8.1, "Family graph: structure shared, per-copy sequence preserved",
        fontsize=18, fontweight='bold', ha='center')
ax.text(7, 7.65,
        "Each exon-class node groups orthologous exons across paralogs — but stores one profile PER copy, not a consensus.",
        fontsize=11, ha='center', color='#444')

# Five exon-class nodes laid out left-to-right
node_x = [1.5, 4.0, 6.5, 9.0, 11.5]
node_labels = ["E1", "E2", "E3", "E4", "E5"]
copies = ['A', 'B', 'C']
copy_colors = ['#1f77b4', '#ff7f0e', '#2ca02c']
# which copies contribute to each node (E2 is paralog-A-specific, E4 paralog-C-specific)
node_copies = [
    [0, 1, 2],   # E1 — all
    [0, 1],      # E2 — A,B only
    [0, 1, 2],   # E3 — all
    [1, 2],      # E4 — B,C only
    [0, 1, 2],   # E5 — all
]

# Draw each node as a stack: header (exon class) + per-copy profile lozenges
node_w = 1.6
node_h_header = 0.5
profile_h = 0.42

for x, lbl, contrib in zip(node_x, node_labels, node_copies):
    n_copies_here = len(contrib)
    total_h = node_h_header + n_copies_here * profile_h + 0.18
    y0 = 4.0 - total_h / 2
    # outer dashed box
    box = FancyBboxPatch((x - node_w/2, y0), node_w, total_h,
                         boxstyle="round,pad=0.05", linestyle='--',
                         facecolor='#f5f5fa', edgecolor='#5a5a8a', lw=1.0)
    ax.add_patch(box)
    # header
    hdr = patches.Rectangle((x - node_w/2 + 0.05, y0 + total_h - node_h_header - 0.05),
                            node_w - 0.10, node_h_header,
                            facecolor='#dcdcf0', edgecolor='#5a5a8a', lw=0.8)
    ax.add_patch(hdr)
    ax.text(x, y0 + total_h - node_h_header/2 - 0.05, f"exon-class {lbl}",
            fontsize=10, ha='center', va='center', fontweight='bold', color='#222')
    # per-copy profiles
    for i, ci in enumerate(contrib):
        py = y0 + 0.05 + i * profile_h
        rect = patches.Rectangle((x - node_w/2 + 0.10, py), node_w - 0.20, profile_h - 0.06,
                                 facecolor=copy_colors[ci], alpha=0.85,
                                 edgecolor='black', lw=0.5)
        ax.add_patch(rect)
        ax.text(x, py + (profile_h - 0.06)/2, f"profile · copy {copies[ci]}",
                fontsize=8.5, ha='center', va='center', color='white', fontweight='bold')

# Edges between nodes (junctions)
def edge(ax, x0, x1, y, color='#888', lw=1.4, dashed=False):
    style = ':' if dashed else '-'
    ax.plot([x0, x1], [y, y], color=color, lw=lw, linestyle=style, zorder=0)

# Draw simplified edges: E1-E2, E1-E3 (paralog C skips E2), E2-E3, E3-E4, E3-E5 (A skips E4),
# E4-E5
edge(ax, 1.5 + node_w/2, 4.0 - node_w/2, 4.0)
edge(ax, 4.0 + node_w/2, 6.5 - node_w/2, 4.0)
edge(ax, 6.5 + node_w/2, 9.0 - node_w/2, 4.0)
edge(ax, 9.0 + node_w/2, 11.5 - node_w/2, 4.0)
# A bypass arc for the paralog-skipping junctions
arc1 = FancyArrowPatch((1.5 + node_w/2, 4.6), (6.5 - node_w/2, 4.6),
                       connectionstyle="arc3,rad=0.25", arrowstyle='-',
                       color='#2ca02c', lw=1.6)
ax.add_patch(arc1)
ax.text((1.5 + 6.5)/2, 5.55, "C skips E2", fontsize=8.5, color='#2ca02c',
        ha='center', fontstyle='italic')

arc2 = FancyArrowPatch((6.5 + node_w/2, 3.4), (11.5 - node_w/2, 3.4),
                       connectionstyle="arc3,rad=-0.25", arrowstyle='-',
                       color='#1f77b4', lw=1.6)
ax.add_patch(arc2)
ax.text((6.5 + 11.5)/2, 2.4, "A skips E4", fontsize=8.5, color='#1f77b4',
        ha='center', fontstyle='italic')

# Legend / what is shared vs not
legend_y = 0.5
ax.text(0.3, legend_y + 0.95, "What's shared:", fontsize=11, fontweight='bold', color='#207020')
ax.text(0.3, legend_y + 0.55,
        "• exon-class identity (which exons are orthologous)\n"
        "• junction topology (which paths exist in the family)\n"
        "• POA-MSA per node (used as a smoothing prior, not a consensus)",
        fontsize=9.5, va='top')

ax.text(7.5, legend_y + 0.95, "What stays per-copy:", fontsize=11, fontweight='bold', color='#a02020')
ax.text(7.5, legend_y + 0.55,
        "• match emissions at every column → preserves SNPs\n"
        "• donor / acceptor microshifts (per-copy column count)\n"
        "• small per-copy indels (each copy's own M/I/D states)",
        fontsize=9.5, va='top')

plt.tight_layout()
plt.savefig('/scratch/jxi21/Assembler/Rustle/bench/meeting_2026_05_06_hmm_em/03_family_graph_per_copy.png',
            dpi=160, bbox_inches='tight')
print("03_family_graph_per_copy.png written")
