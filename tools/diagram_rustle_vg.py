#!/usr/bin/env python3
"""Diagram: how Rustle's VG mode extends StringTie's flow approach to
multi-copy gene families."""
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle, Circle, FancyBboxPatch, FancyArrowPatch

fig, axes = plt.subplots(4, 1, figsize=(14, 13), gridspec_kw={'hspace': 0.45,
    'height_ratios': [1, 1.4, 1.2, 1.2]})

# ===========================================================================
# Panel A: The starting point — three paralog splice graphs (same topology)
# ===========================================================================
ax1 = axes[0]
ax1.set_xlim(-2, 115); ax1.set_ylim(-0.6, 3.5)
ax1.set_title("Step 0.  Three tandem paralog copies on the genome (same topology, shifted coords)",
              loc='left', fontsize=12, fontweight='bold')

def draw_gene(ax, y, start_x, offset_lbl, color='steelblue', label=''):
    """Draw a gene as exons at (start_x) + pattern."""
    exon_pattern = [(0, 4), (8, 10), (15, 17), (21, 24), (28, 34)]  # 5-exon
    for s, e in exon_pattern:
        ax.add_patch(Rectangle((start_x + s, y - 0.2), e - s, 0.4,
                              facecolor=color, edgecolor='black', lw=0.5))
    ax.plot([start_x + exon_pattern[0][0], start_x + exon_pattern[-1][1]],
            [y, y], color='gray', lw=0.4, linestyle=':')
    ax.text(start_x - 2, y, label, fontsize=10, va='center', ha='right', fontweight='bold')
    ax.text(start_x + 17, y + 0.5, offset_lbl, fontsize=8, color='#555', ha='center', style='italic')

draw_gene(ax1, 2.5, 10, "ref_start",  '#2e7d32', "L7_1")
draw_gene(ax1, 1.5, 45, "Δ = +40 kb", '#2e7d32', "L7_2")
draw_gene(ax1, 0.5, 78, "Δ = +80 kb", '#2e7d32', "L7_3")

# Add note
ax1.text(7, 3.15, "Genomic organization:",
         fontsize=10, fontweight='bold', color='#333')
ax1.text(110, 2.5, "26 primary reads", fontsize=9, va='center', color='#333')
ax1.text(110, 1.5, "6 primary reads",  fontsize=9, va='center', color='#c00')
ax1.text(110, 0.5, "5 primary reads",  fontsize=9, va='center', color='#c00')
ax1.axis('off')

# ===========================================================================
# Panel B: StringTie's view — only primary alignments
# ===========================================================================
ax2 = axes[1]
ax2.set_xlim(-2, 115); ax2.set_ylim(-0.5, 4)
ax2.set_title("Step 1.  StringTie view — drops secondary alignments → copies 2 & 3 starve",
              loc='left', fontsize=12, fontweight='bold')

def draw_copy_graph(ax, y, start_x, n_reads, label, read_color='#ccac42'):
    """Draw a mini splice graph with reads."""
    exon_pattern = [(0, 4), (8, 10), (15, 17), (21, 24), (28, 34)]
    # Exons
    for s, e in exon_pattern:
        ax.add_patch(Rectangle((start_x + s, y - 0.2), e - s, 0.4,
                              facecolor='steelblue', edgecolor='black', lw=0.5))
    ax.plot([start_x + exon_pattern[0][0], start_x + exon_pattern[-1][1]],
            [y, y], color='gray', lw=0.4, linestyle=':')
    # Reads — stack vertically under the gene
    for r in range(min(n_reads, 4)):
        ry = y - 0.4 - r * 0.15
        for s, e in exon_pattern:
            ax.add_patch(Rectangle((start_x + s, ry - 0.04), e - s, 0.1,
                                   facecolor=read_color, edgecolor='none'))
    if n_reads > 4:
        ax.text(start_x + 16, y - 1.15, f"+ {n_reads-4} more", fontsize=7, ha='center',
                style='italic', color='#666')
    ax.text(start_x - 2, y, label, fontsize=10, va='center', ha='right', fontweight='bold')

# L7_1 gets all reads: assembles
draw_copy_graph(ax2, 3, 10, 26, "L7_1")
ax2.text(46, 3, "→  ISOFORM ASSEMBLED  ✓", fontsize=10, color='#1b5e20', fontweight='bold', va='center')

# L7_2 gets 6 reads, too few to reach TSS
draw_copy_graph(ax2, 1.5, 10, 6, "L7_2")
ax2.text(46, 1.5, "→  NO ASSEMBLY — reads don't reach TSS exon  ✗",
         fontsize=10, color='#b71c1c', fontweight='bold', va='center')

# L7_3 gets 5 reads
draw_copy_graph(ax2, 0, 10, 5, "L7_3")
ax2.text(46, 0, "→  NO ASSEMBLY — same problem  ✗",
         fontsize=10, color='#b71c1c', fontweight='bold', va='center')
ax2.axis('off')

# ===========================================================================
# Panel C: Rustle's VG view — all 3 copies + multi-mapping reads
# ===========================================================================
ax3 = axes[2]
ax3.set_xlim(-2, 115); ax3.set_ylim(-1.5, 4.5)
ax3.set_title("Step 2.  Rustle VG view — add (a) cross-copy multi-mapping reads and (b) chain-homology detection",
              loc='left', fontsize=12, fontweight='bold')

# Three copies' splice graphs
for i, (y, lbl, x0) in enumerate([(3, "L7_1", 10), (1.8, "L7_2", 10), (0.6, "L7_3", 10)]):
    exon_pattern = [(0, 4), (8, 10), (15, 17), (21, 24), (28, 34)]
    for s, e in exon_pattern:
        ax3.add_patch(Rectangle((x0 + s, y - 0.2), e - s, 0.4,
                               facecolor='steelblue', edgecolor='black', lw=0.5))
    ax3.plot([x0 + exon_pattern[0][0], x0 + exon_pattern[-1][1]],
             [y, y], color='gray', lw=0.4, linestyle=':')
    ax3.text(x0 - 2, y, lbl, fontsize=10, va='center', ha='right', fontweight='bold')

# Arrows showing multi-mapping reads between copies (with weight 1/N)
for arrow_y_from, arrow_y_to, caption in [
    (3, 1.8, "1/3"),
    (3, 0.6, "1/3"),
    (1.8, 0.6, "1/3"),
]:
    arrow = FancyArrowPatch((48, arrow_y_from), (48, arrow_y_to),
                            arrowstyle='<->',
                            mutation_scale=12, color='#ff6f00',
                            linewidth=1.3, linestyle='--',
                            connectionstyle="arc3,rad=-0.3")
    ax3.add_patch(arrow)

ax3.text(52, 2.1, "multi-mapper share", fontsize=9, color='#ff6f00', va='center', style='italic')
ax3.text(52, 1.2, "(weight = 1/NH)",   fontsize=9, color='#ff6f00', va='center', style='italic')

# Chain-homology annotation
ax3.annotate("", xy=(10, -0.8), xytext=(44, -0.8),
             arrowprops=dict(arrowstyle='<->', color='#00796b', lw=1.5, linestyle='-'))
ax3.text(27, -0.5, "chain homology — same 8-intron topology, ±20bp drift",
         fontsize=9.5, color='#00796b', ha='center', style='italic', fontweight='bold')

ax3.text(68, 2.4, "Each copy now gets evidence from:", fontsize=9.5, color='#333', fontweight='bold')
ax3.text(68, 1.9, "• its own primary reads", fontsize=9)
ax3.text(68, 1.4, "• fractional (1/NH) secondary reads from sister copies", fontsize=9)
ax3.text(68, 0.9, "• chain-homology bootstrap from a completed sister assembly", fontsize=9)

ax3.axis('off')

# ===========================================================================
# Panel D: Final result — SNPs help copy assignment; all 3 assembled
# ===========================================================================
ax4 = axes[3]
ax4.set_xlim(-2, 115); ax4.set_ylim(-2, 4.5)
ax4.set_title("Step 3.  Resolve copy identity via SNPs + flow redistribution → all 3 assembled",
              loc='left', fontsize=12, fontweight='bold')

# SNP panel on left
exon_pattern = [(0, 4), (8, 10), (15, 17), (21, 24), (28, 34)]
for idx, (y, lbl, col, snps) in enumerate([
    (3.3, "L7_1", '#388e3c', {2: 'A', 9: 'C', 22: 'G'}),
    (2.0, "L7_2", '#388e3c', {2: 'T', 9: 'C', 22: 'G'}),  # differs at pos 2
    (0.7, "L7_3", '#388e3c', {2: 'A', 9: 'G', 22: 'G'}),  # differs at pos 9
]):
    for s, e in exon_pattern:
        ax4.add_patch(Rectangle((10 + s, y - 0.2), e - s, 0.4,
                               facecolor=col, edgecolor='black', lw=0.5))
    ax4.plot([10 + exon_pattern[0][0], 10 + exon_pattern[-1][1]],
             [y, y], color='gray', lw=0.4, linestyle=':')
    ax4.text(8, y, lbl, fontsize=10, va='center', ha='right', fontweight='bold')
    ax4.text(45, y, "✓ EXACT MATCH (= class)",
             fontsize=10, color='#1b5e20', fontweight='bold', va='center')
    # SNP ticks
    for pos, base in snps.items():
        ax4.plot([10 + pos, 10 + pos], [y + 0.25, y + 0.4], color='crimson', lw=1.5)
        ax4.text(10 + pos, y + 0.55, base, ha='center', va='center',
                 fontsize=8, color='crimson', fontweight='bold')

# Summary
ax4.text(68, 3.5, "For each multi-mapping read:", fontsize=10, color='#333', fontweight='bold')
ax4.text(68, 3.0, "• compare its bases at SNP positions vs each copy", fontsize=9)
ax4.text(68, 2.5, "• assign weight proportional to alignment score (NM, AS)", fontsize=9)
ax4.text(68, 2.0, "• EM iterates: reads → weights → assembly → weights", fontsize=9)
ax4.text(68, 1.5, "• extend partial reconstructions via family-graph projection", fontsize=9)

ax4.text(50, -1.3,
         "Result on GOLGA6L7:  StringTie 1/3 exact matches  →  Rustle VG  3/3 exact matches",
         ha='center', fontsize=12, fontweight='bold', color='#1b5e20',
         bbox=dict(boxstyle='round,pad=0.6', facecolor='#e8f5e9', edgecolor='#1b5e20', lw=1.5))

ax4.axis('off')

fig.suptitle("How Rustle extends StringTie's flow model for multi-copy gene families",
             fontsize=14.5, y=0.995, fontweight='bold')

plt.savefig('/storage/group/kdm16/default/jxi21/apes_transcriptome_analysis/clusterer/figures/fig3_rustle_vg_extension.png', dpi=150, bbox_inches='tight')
plt.savefig('/storage/group/kdm16/default/jxi21/apes_transcriptome_analysis/clusterer/figures/fig3_rustle_vg_extension.pdf', bbox_inches='tight')
print("Saved fig3_rustle_vg_extension.{png,pdf}")
