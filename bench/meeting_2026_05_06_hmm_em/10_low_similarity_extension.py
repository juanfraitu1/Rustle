#!/usr/bin/env python3
"""Slide 10 — additive extension concept for low-similarity paralogs (jaccard < 0.30).

Two-pronged extension that REUSES the existing HMM machinery, doesn't override it:

  (a) Family-discovery side — add a complementary linking signal that survives
      sequence drift (splice-site / exon-boundary Jaccard, conserved-domain
      anchors). The k-mer linker stays at jaccard ≥ 0.30 default.

  (b) Scoring side — for low-sim families, swap forward_against_path_for_copy
      (rigid: read MUST fit one paralog's exact path) for assign_via_graph_viterbi
      (relaxed: read picks its own best path through the graph; assign by
      node-overlap with each paralog's known path). Same M/I/D trellis inside
      every node.

Dispatch by similarity band — auto-solver routes high/medium → existing path
forward, low → graph-Viterbi. Nothing existing changes."""
import matplotlib.pyplot as plt
import matplotlib.patches as patches
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch

fig, ax = plt.subplots(figsize=(16, 10))
ax.set_xlim(0, 16)
ax.set_ylim(0, 10.5)
ax.axis('off')

ax.text(8, 10.05, "Slide 10  —  Extending to low-similarity paralogs (jaccard < 0.30)",
        fontsize=17, fontweight='bold', ha='center')
ax.text(8, 9.65,
        "Additive — nothing existing changes.  HMM trellis stays the same.  Auto-solver routes per family by similarity band.",
        fontsize=11, ha='center', color='#444', fontstyle='italic')

# ── Top band: jaccard scale (mirrors slide 9) ───────────────────────────────
SX0, SX1 = 1.0, 15.0
SYT, SYB = 9.0, 8.55

def jx(j):
    return SX0 + (SX1 - SX0) * j

bands = [
    (0.00, 0.30, "low",            "#dde2eb", True),
    (0.30, 0.60, "medium",         "#fff3cd", False),
    (0.60, 0.90, "high",           "#d4edda", False),
    (0.90, 1.00, "near-identical", "#cce4f7", False),
]
for j0, j1, name, color, highlight in bands:
    x0, x1 = jx(j0), jx(j1)
    rect = patches.Rectangle(
        (x0, SYB), x1 - x0, SYT - SYB,
        facecolor=color,
        edgecolor='#a02060' if highlight else '#888',
        lw=2.0 if highlight else 0.8,
        zorder=1,
    )
    ax.add_patch(rect)
    cx = (x0 + x1) / 2
    ax.text(cx, (SYT + SYB) / 2 + 0.10, f"{name}-similarity",
            fontsize=9.5, ha='center', va='center', fontweight='bold', color='#222')
    ax.text(cx, (SYT + SYB) / 2 - 0.13,
            f"jaccard {j0:.1f} – {j1:.1f}",
            fontsize=8.5, ha='center', va='center', color='#444')
# Tick label
ax.text((SX0 + SX1)/2, SYB - 0.18,
        "k-mer Jaccard (paralog vs nearest sibling)",
        fontsize=9, ha='center', color='#666', fontstyle='italic')

# Pointer to the new target band — placed BELOW the scale so it doesn't
# collide with the title.
ax.annotate("↑ this band — today's POC target",
            xy=(jx(0.15), SYB - 0.05), xytext=(jx(0.15), SYB - 0.55),
            fontsize=10, ha='center', color='#a02060', fontweight='bold',
            arrowprops=dict(arrowstyle='->', color='#a02060', lw=1.4))

# ── Two-pronged extension diagram ───────────────────────────────────────────
# LEFT column: family-discovery side
# RIGHT column: scoring side
COL_W = 7.2
LX0, RX0 = 0.5, 8.3

# Box style helper
def proposal_box(ax, x, y, w, h, title, color, body_lines, footer):
    box = FancyBboxPatch((x, y), w, h,
                         boxstyle="round,pad=0.10",
                         facecolor='#fafafa', edgecolor=color, lw=1.6)
    ax.add_patch(box)
    ax.text(x + 0.20, y + h - 0.30, title,
            fontsize=12, fontweight='bold', color=color, va='top')
    cur_y = y + h - 0.85
    for line in body_lines:
        ax.text(x + 0.30, cur_y, line, fontsize=9.8, va='top', color='#222')
        cur_y -= 0.40
    ax.text(x + w/2, y + 0.25, footer,
            fontsize=9, ha='center', va='center', color='#666', fontstyle='italic')

# (a) family-discovery side
left_lines = [
    "Existing k-mer linker (jaccard ≥ 0.30) stays the default.",
    "ADD a complementary signal that survives sequence drift:",
    "  •  splice-site Jaccard  —  donor / acceptor positions",
    "      are often conserved when the sequence has diverged",
    "  •  conserved-domain anchors  —  short Pfam / InterPro motifs",
    "      identifiable directly from genomic k-mers",
    "Two bundles linked into the same family if EITHER signal fires.",
    "Gated by --vg-low-sim-link  (off by default).",
]
proposal_box(ax, LX0, 4.0, COL_W, 4.10,
             "(a)  Family-discovery side",
             '#1f77b4',
             left_lines,
             "Otherwise: the two paralogs never enter the same FamilyGraph.")

# (b) scoring side
right_lines = [
    "Existing per-copy path forward DP stays the default.",
    "ADD: graph-relaxed assignment for low-sim families:",
    "  •  run viterbi_path(fg, read)  —  same M/I/D trellis,",
    "      no per-paralog path constraint",
    "  •  for each known paralog c, score by node-overlap:",
    "      recall(c)    =  | viterbi ∩ path_c |  /  | path_c |",
    "      precision(c) =  | viterbi ∩ path_c |  /  | viterbi |",
    "Sort by recall desc.  Top paralog = the assignment.",
]
proposal_box(ax, RX0, 4.0, COL_W, 4.10,
             "(b)  Scoring side",
             '#207020',
             right_lines,
             "Function added: assign_via_graph_viterbi() in scorer.rs (POC, 3 tests passing).")

# ── Dispatch logic at the bottom ────────────────────────────────────────────
disp_box = FancyBboxPatch((0.5, 0.50), 15.0, 3.20,
                          boxstyle="round,pad=0.10",
                          facecolor='#f4f6fa', edgecolor='#5a5a8a', lw=1.4)
ax.add_patch(disp_box)
ax.text(8, 3.45, "Dispatch — per family, automatic by similarity band",
        fontsize=12, fontweight='bold', ha='center', color='#222')

# ASCII-style code box centered
code_lines = [
    "match family_similarity_band(fg) {",
    "    Band::High  | Band::Medium  =>  forward_against_path_for_copy(fg, read, path_c, c)   // existing — slides 5–7",
    "    Band::Low                   =>  assign_via_graph_viterbi(fg, read)                   // new — slide 10",
    "    Band::NearIdentical         =>  forward_against_path_for_copy + SNP-aware E-step    // existing --vg-snp",
    "}",
]
y = 2.85
for line in code_lines:
    ax.text(1.0, y, line, fontsize=10, family='monospace', color='#222')
    y -= 0.35

ax.text(8, 0.95,
        "M-step (priors) and gap rule are unchanged — only the per-(read, copy) score in the E-step differs.\n"
        "Gap rule generalizes: for graph-Viterbi, abstain if  top.recall − second.recall < δ_recall  (proposed default δ_recall = 0.20).",
        fontsize=10, ha='center', va='center', color='#333')

# ── Right margin: POC results callout ───────────────────────────────────────
poc_box = FancyBboxPatch((11.0, 0.55), 4.5, 2.6,
                         boxstyle="round,pad=0.08",
                         facecolor='#ffffff', edgecolor='#a02060', lw=1.2,
                         zorder=5)
# Put it at top-right of the dispatch area instead — overlapping might collide.
# Actually let me put a small POC result note overlaying the bottom-right empty
# area of the figure — but only if there's space. Let me skip it; the dispatch
# block already mentions "POC, 3 tests passing".

plt.tight_layout()
plt.savefig('/scratch/jxi21/Assembler/Rustle/bench/meeting_2026_05_06_hmm_em/10_low_similarity_extension.png',
            dpi=160, bbox_inches='tight')
print("10_low_similarity_extension.png written")
