#!/usr/bin/env python3
"""FLNC ≠ post-cluster2 isoform.

Advisor clarification: FLNC reads are Full-Length Non-Concatemer reads from
the IsoSeq pipeline — they have polyA + 5' cap signals intact, but they have
NOT been collapsed to per-isoform consensus by `isoseq cluster2`. Two
implications:

1. An FLNC read may NOT span the full transcript — it can be a partial
   sub-sequence. Truncations and 5'/3' biases are common.
2. A multi-mapping FLNC read may only span a CONSERVED DOMAIN shared
   across paralogs (the only region where they look alike).

Why per-copy profile HMM-EM still works in case (2): scoring is local along
the path the read covers. SNPs WITHIN that conserved domain still tip the
forward log-likelihood — even if the read never sees the easy paralog-
distinguishing 5' or 3' UTRs."""
import matplotlib.pyplot as plt
import matplotlib.patches as patches
from matplotlib.patches import FancyBboxPatch

fig, ax = plt.subplots(figsize=(15, 8.5))
ax.set_xlim(0, 15)
ax.set_ylim(0, 9.0)
ax.axis('off')

ax.text(7.5, 8.55, "FLNC ≠ post-cluster2 isoform — and why per-copy profiles still work",
        fontsize=16, fontweight='bold', ha='center')
ax.text(7.5, 8.15,
        "FLNC reads have polyA+cap intact but no per-isoform consensus has been built yet. "
        "Many are partial. Multi-mappers often only span conserved domains.",
        fontsize=10.5, ha='center', color='#444')

# Three rows of illustrations — exons start at x=4.5 to leave room for row labels.
EX_LEFT = 4.5
exon_xs = [(EX_LEFT + 0.0, EX_LEFT + 1.0),
           (EX_LEFT + 1.3, EX_LEFT + 2.4),
           (EX_LEFT + 2.7, EX_LEFT + 3.9),
           (EX_LEFT + 4.2, EX_LEFT + 5.1),
           (EX_LEFT + 5.4, EX_LEFT + 6.3),
           (EX_LEFT + 6.6, EX_LEFT + 7.8),
           (EX_LEFT + 8.1, EX_LEFT + 9.0)]
TRACK_LEFT = EX_LEFT - 0.4
TRACK_RIGHT = EX_LEFT + 9.1
exon_color_paralog_A = '#1f77b4'
exon_color_paralog_B = '#ff7f0e'

# ── Row 1: cluster2 isoform = full-length consensus ─────────────────────────
y1 = 6.7
ax.text(0.3, y1, "(1) post-cluster2 isoform\n(ideal but downstream)",
        fontsize=10, ha='left', va='center', color='#333')
# transcript backbone
ax.plot([TRACK_LEFT, TRACK_RIGHT], [y1, y1], color='#888', lw=0.8)
for x0, x1 in exon_xs:
    ax.add_patch(patches.Rectangle((x0, y1 - 0.18), x1 - x0, 0.36,
                                   facecolor='#5dade2', edgecolor='black', lw=0.5))
ax.text(TRACK_RIGHT + 0.2, y1, "spans the entire\ntranscript",
        fontsize=9.5, color='#666', va='center', fontstyle='italic')

# ── Row 2: full-length FLNC read (still possible) ───────────────────────────
y2 = 5.0
ax.text(0.3, y2, "(2) FLNC, full-length\n(some reads — best case)",
        fontsize=10, ha='left', va='center', color='#333')
ax.plot([TRACK_LEFT, TRACK_RIGHT], [y2, y2], color='#888', lw=0.8)
for x0, x1 in exon_xs:
    ax.add_patch(patches.Rectangle((x0, y2 - 0.18), x1 - x0, 0.36,
                                   facecolor='#d62728', alpha=0.85,
                                   edgecolor='black', lw=0.5))
ax.text(TRACK_RIGHT + 0.2, y2, "polyA + cap intact,\nspans full transcript",
        fontsize=9.5, color='#666', va='center', fontstyle='italic')

# ── Row 3: PARTIAL FLNC read ────────────────────────────────────────────────
y3 = 3.3
ax.text(0.3, y3, "(3) FLNC, partial\n(common — internal\ndegradation, truncation)",
        fontsize=10, ha='left', va='center', color='#333')
ax.plot([TRACK_LEFT, TRACK_RIGHT], [y3, y3], color='#bbb', lw=0.5, linestyle=':')
# only middle 3 exons present
for x0, x1 in exon_xs[2:5]:
    ax.add_patch(patches.Rectangle((x0, y3 - 0.18), x1 - x0, 0.36,
                                   facecolor='#d62728', alpha=0.85,
                                   edgecolor='black', lw=0.5))
ax.text(TRACK_RIGHT + 0.2, y3, "spans only a sub-region;\nedges may be missing",
        fontsize=9.5, color='#666', va='center', fontstyle='italic')

# ── Row 4: multi-mapper at a conserved domain ───────────────────────────────
y4 = 1.6
ax.text(0.3, y4, "(4) FLNC multi-mapper\nat a conserved domain",
        fontsize=10, ha='left', va='center', color='#a02020', fontweight='bold')
# Two paralog tracks at this row, one slightly above + below
y4a = y4 + 0.45
y4b = y4 - 0.45
ax.text(TRACK_RIGHT + 0.2, y4a, "paralog A locus",
        fontsize=8.5, color=exon_color_paralog_A, va='center')
ax.text(TRACK_RIGHT + 0.2, y4b, "paralog B locus",
        fontsize=8.5, color=exon_color_paralog_B, va='center')
# paralog A backbone + exons
ax.plot([TRACK_LEFT, TRACK_RIGHT], [y4a, y4a], color='#aaa', lw=0.7)
for x0, x1 in exon_xs:
    ax.add_patch(patches.Rectangle((x0, y4a - 0.13), x1 - x0, 0.26,
                                   facecolor=exon_color_paralog_A, alpha=0.5,
                                   edgecolor='black', lw=0.4))
# paralog B backbone + exons
ax.plot([TRACK_LEFT, TRACK_RIGHT], [y4b, y4b], color='#aaa', lw=0.7)
for x0, x1 in exon_xs:
    ax.add_patch(patches.Rectangle((x0, y4b - 0.13), x1 - x0, 0.26,
                                   facecolor=exon_color_paralog_B, alpha=0.5,
                                   edgecolor='black', lw=0.4))

# Highlight conserved-domain region (e.g. exons 4 + 5)
cd_x0 = exon_xs[3][0] - 0.15
cd_x1 = exon_xs[4][1] + 0.15
ax.add_patch(patches.Rectangle((cd_x0, y4b - 0.7), cd_x1 - cd_x0, 1.55,
                               facecolor='#fff3cd', alpha=0.55,
                               edgecolor='#caa845', lw=1.0,
                               zorder=0))
ax.text((cd_x0 + cd_x1)/2, y4b - 0.95, "conserved domain\n(only region shared at high identity)",
        fontsize=9, ha='center', color='#7a5500')

# Show the read mapping to BOTH paralogs in the conserved domain
for y in [y4a, y4b]:
    rd = patches.Rectangle((cd_x0 + 0.10, y - 0.08), cd_x1 - cd_x0 - 0.20, 0.16,
                            facecolor='#d62728', edgecolor='black', lw=0.5,
                            zorder=3)
    ax.add_patch(rd)

# Right-side commentary
commentary = FancyBboxPatch((0.3, 0.05), 14.4, 0.85,
                            boxstyle="round,pad=0.10",
                            facecolor='#e8f5e8', edgecolor='#207020', lw=1.2)
ax.add_patch(commentary)
ax.text(7.5, 0.65,
        "Even when an FLNC multi-mapper only sees a conserved domain, per-copy profile HMM-EM still scores it locally:",
        fontsize=10.5, ha='center', fontweight='bold', color='#207020')
ax.text(7.5, 0.30,
        "the forward DP integrates only over the read's path slice — every per-copy SNP within that slice contributes; "
        "if the gap is too small, the score-gap rule abstains and falls back to 1/NH.",
        fontsize=10, ha='center', color='#225522', fontstyle='italic')

plt.tight_layout()
plt.savefig('/scratch/jxi21/Assembler/Rustle/bench/meeting_2026_05_06_hmm_em/08_flnc_partial_coverage.png',
            dpi=160, bbox_inches='tight')
print("08_flnc_partial_coverage.png written")
