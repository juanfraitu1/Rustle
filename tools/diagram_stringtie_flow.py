#!/usr/bin/env python3
"""Diagram: how StringTie decomposes a splice graph into isoforms via max flow."""
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle, Circle

fig, axes = plt.subplots(3, 1, figsize=(14, 12), gridspec_kw={'hspace': 0.35})

# ---- Panel A ----
ax1 = axes[0]
ax1.set_xlim(0, 115); ax1.set_ylim(-0.3, 6.2)
ax1.set_title("Step 1.  Reads aligned to the reference → build splice graph",
              loc='left', fontsize=12, fontweight='bold')

exons = [(5, 18), (28, 38), (48, 57), (66, 75), (85, 96)]
for i, (s, e) in enumerate(exons):
    ax1.add_patch(Rectangle((s, 4.5), e-s, 0.9, facecolor='steelblue',
                           edgecolor='black', linewidth=0.8))
    ax1.text((s+e)/2, 4.95, f"E{i+1}", ha='center', va='center',
             color='white', fontsize=10, fontweight='bold')
ax1.plot([exons[0][1], exons[-1][0]], [4.95, 4.95], color='gray', lw=0.4, linestyle=':')

read_patterns = [
    ([(5, 18), (28, 38), (48, 57), (66, 75), (85, 96)], "read 1 — full transcript"),
    ([(28, 38), (48, 57), (66, 75), (85, 96)],          "read 2 — skips E1"),
    ([(5, 18), (28, 38), (48, 57), (85, 96)],           "read 3 — skips E4"),
    ([(5, 18), (48, 57), (66, 75), (85, 96)],           "read 4 — skips E2"),
    ([(5, 18), (28, 38), (48, 57), (66, 75)],           "read 5 — no E5"),
]
for idx, (exons_r, lbl) in enumerate(read_patterns):
    y = 3.4 - idx * 0.55
    for es, ee in exons_r:
        ax1.add_patch(Rectangle((es, y-0.12), ee-es, 0.25,
                               facecolor='#ccac42', edgecolor='black', linewidth=0.3))
    xs = sorted(exons_r)
    ax1.plot([xs[0][0], xs[-1][1]], [y, y], color='gray', lw=0.4, linestyle=':')
    ax1.text(99, y, lbl, fontsize=9, va='center')
ax1.axis('off')

# ---- Panel B: splice graph DAG ----
ax2 = axes[1]
ax2.set_xlim(-5, 115); ax2.set_ylim(-4.5, 8)
ax2.set_title("Step 2.  Splice graph — nodes = exons, edges = junctions, edge capacity = read support",
              loc='left', fontsize=12, fontweight='bold')

node_positions = {
    'S':  (0, 1), 'E1': (15, 1), 'E2': (32, 1), 'E3': (50, 1),
    'E4': (68, 1), 'E5': (86, 1), 'T': (102, 1),
}
edges = [
    ('S', 'E1', 4, False), ('S', 'E2', 1, True),
    ('E1', 'E2', 3, False), ('E1', 'E3', 1, True),
    ('E2', 'E3', 4, False),
    ('E3', 'E4', 4, False), ('E3', 'E5', 1, True),
    ('E4', 'E5', 3, False), ('E4', 'T', 1, True),
    ('E5', 'T', 4, False),
]

# Nodes
for name, (x, y) in node_positions.items():
    if name in ('S', 'T'):
        ax2.add_patch(Circle((x, y), 2.5, facecolor='lightgray', edgecolor='black', lw=1.2))
    else:
        ax2.add_patch(Circle((x, y), 3, facecolor='steelblue', edgecolor='black', lw=1.2))
    col = 'black' if name in ('S','T') else 'white'
    ax2.text(x, y, name, ha='center', va='center', fontsize=11, fontweight='bold', color=col)

# Edges
for u, v, cap, is_skip in edges:
    ux, uy = node_positions[u]; vx, vy = node_positions[v]
    ru = 2.5 if u in ('S','T') else 3
    rv = 2.5 if v in ('S','T') else 3
    dx = vx - ux; dy = vy - uy
    L = (dx**2 + dy**2)**0.5
    ux2 = ux + ru*dx/L; uy2 = uy + ru*dy/L
    vx2 = vx - rv*dx/L; vy2 = vy - rv*dy/L
    if is_skip:
        ax2.annotate("", xy=(vx2, vy2), xytext=(ux2, uy2),
                    arrowprops=dict(arrowstyle='->', color='#d62728', lw=1.4,
                                    connectionstyle="arc3,rad=-0.35"))
        midx = (ux + vx)/2
        # Fixed label height per arc-span so labels stay inside the panel
        arc_span = abs(vx - ux)
        if arc_span < 25: midy = uy + 4.0
        elif arc_span < 45: midy = uy + 5.5
        else: midy = uy + 6.5
        ax2.text(midx, midy, str(cap), ha='center', va='center', fontsize=10,
                 color='#d62728', fontweight='bold',
                 bbox=dict(facecolor='white', edgecolor='#d62728', pad=2, boxstyle='round,pad=0.3'))
    else:
        ax2.annotate("", xy=(vx2, vy2), xytext=(ux2, uy2),
                    arrowprops=dict(arrowstyle='->', color='black', lw=1.4))
        midx = (ux + vx)/2
        ax2.text(midx, uy - 2.2, str(cap), ha='center', va='center', fontsize=10,
                 color='black', fontweight='bold',
                 bbox=dict(facecolor='white', edgecolor='none', pad=1))

# Legend
ax2.text(4, -3.5, "●", fontsize=14, color='black'); ax2.text(7, -3.5, "consecutive junction", fontsize=9, va='center')
ax2.text(45, -3.5, "●", fontsize=14, color='#d62728'); ax2.text(48, -3.5, "exon-skipping junction (alternative splice)", fontsize=9, va='center')
ax2.axis('off')

# ---- Panel C ----
ax3 = axes[2]
ax3.set_xlim(0, 115); ax3.set_ylim(-0.5, 3.5)
ax3.set_title("Step 3.  Max-flow decomposition → extract heaviest path, subtract flow, repeat",
              loc='left', fontsize=12, fontweight='bold')

isoforms = [
    (3, [(5,18),(28,38),(48,57),(66,75),(85,96)], "iter 1:  heaviest path (flow = 3)   →   ISOFORM 1   (full transcript)"),
    (1, [(5,18),(28,38),(48,57),(85,96)],         "iter 2:  next heaviest (flow = 1)   →   ISOFORM 2   (skips E4)"),
    (1, [(5,18),(48,57),(66,75),(85,96)],         "iter 3:  next heaviest (flow = 1)   →   ISOFORM 3   (skips E2)"),
]
colors = ['#1b5e20', '#388e3c', '#66bb6a']
for i, (flow, exons_i, lbl) in enumerate(isoforms):
    y = 2.6 - i * 0.9
    exs = sorted(exons_i)
    ax3.plot([exs[0][0], exs[-1][1]], [y, y], color='gray', lw=0.6)
    for s, e in exs:
        ax3.add_patch(Rectangle((s, y-0.2), e-s, 0.4,
                               facecolor=colors[i], edgecolor='black', lw=0.5))
    ax3.text(99, y, lbl, fontsize=9, va='center')

ax3.text(50, -0.1,
         "Flow conservation:  Σ (isoform flow at edge)  =  edge capacity",
         ha='center', fontsize=9.5, style='italic', color='#555555')
ax3.axis('off')

fig.suptitle("How StringTie assembles transcripts — network-flow decomposition",
             fontsize=15, y=0.985, fontweight='bold')

plt.savefig('/storage/group/kdm16/default/jxi21/apes_transcriptome_analysis/clusterer/figures/fig2_stringtie_flow.png', dpi=150, bbox_inches='tight')
plt.savefig('/storage/group/kdm16/default/jxi21/apes_transcriptome_analysis/clusterer/figures/fig2_stringtie_flow.pdf', bbox_inches='tight')
print("Saved fig2_stringtie_flow.{png,pdf}")
