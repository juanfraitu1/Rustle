#!/usr/bin/env python3
"""Slide 13 — Where the per-paralog likelihoods log P(r | c) actually
come from.

The advisor's deepest skepticism: "deriving the likelihoods looks like
an unsolvable problem." Answer: we are NOT deriving them. Every
probability is an observed frequency from the multi-copy MSA. The
likelihood is just a chain of log-lookups summed along the read's
alignment path. This slide walks from raw data → emission tables →
per-base log additions → the concrete log P(r | c) numbers that feed
slides 6 and 12."""
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch

fig, ax = plt.subplots(figsize=(16, 19))
ax.set_xlim(0, 16)
ax.set_ylim(0, 19)
ax.axis('off')

# ── Title ───────────────────────────────────────────────────────────────
ax.text(8, 18.55, "Slide 13  —  Where do the per-paralog likelihoods come from?",
        fontsize=17, fontweight='bold', ha='center')
ax.text(8, 18.05,
        "log P(r | c) is NOT derived analytically.  It is a sum of observed-frequency log-emissions counted directly from the multi-copy MSA.",
        fontsize=10.8, ha='center', color='#444', fontstyle='italic')
ax.text(8, 17.70,
        "No optimization, no MCMC, no closed-form derivation — just lookups + addition.",
        fontsize=10.8, ha='center', color='#444', fontstyle='italic')


def block(x, y, w, h, title, lines, color='#5a82bc', bg='#eef6ff',
          font=9.6, line_h=0.30, title_font=12, top_pad=0.85):
    box = FancyBboxPatch((x, y), w, h, boxstyle="round,pad=0.10",
                         facecolor=bg, edgecolor=color, lw=1.2)
    ax.add_patch(box)
    ax.text(x + 0.20, y + h - 0.30, title,
            fontsize=title_font, fontweight='bold', color=color, va='top')
    cur_y = y + h - top_pad
    for line in lines:
        if isinstance(line, tuple):
            text, fam = line
        else:
            text, fam = line, 'monospace'
        ax.text(x + 0.30, cur_y, text, fontsize=font, family=fam,
                va='top', color='#222')
        cur_y -= line_h


# ── Step A: Observed MSA → per-copy emission distributions ─────────────
stepA_lines = [
    ("INPUT (observed): the multi-copy MSA at one exon-class node.  "
     "Each row is one observed copy of a paralog, aligned to common columns.",
     'sans-serif'),
    "",
    "                          col 1   col 2   col 3   col 4   col 5   col 6   col 7",
    "  paralog A, copy 1:        A       T       G       C       A       G       T",
    "  paralog A, copy 2:        A       T       G       C       A       G       T",
    "  paralog A, copy 3:        A       T       G       C       A       G       T",
    "  paralog B, copy 1:        A       T       G       C       G       G       T     ← col 5 differs",
    "  paralog B, copy 2:        A       T       G       C       G       G       T",
    "",
    ("Per-copy emission distribution = α · paralog-specific column counts  +  (1−α) · POA-MSA shared column.",
     'sans-serif'),
    ("With α = 0.7 and ε = 0.02 pseudocount (slide 4):",
     'sans-serif'),
    "",
    "  e_M(b | col 5, paralog A) = (A=0.95,  C=0.02,  G=0.02,  T=0.01)    ← from 3/3 'A'",
    "  e_M(b | col 5, paralog B) = (A=0.02,  C=0.02,  G=0.95,  T=0.01)    ← from 2/2 'G'",
    "",
    ("Cols 1-4 and 6-7 are consensus columns (all 5 copies agree) — same emission distribution for both paralogs.",
     'sans-serif'),
]
block(0.40, 12.65, 15.20, 4.90, "A.  Where the emission distributions come from  (observed MSA)",
      stepA_lines, color='#1f77b4', bg='#eef6ff', line_h=0.275)

# ── Step B: Read + alignment → per-base column lookups ─────────────────
stepB_lines = [
    ("INPUT (observed): the read's BAM placement gives genomic coordinates.  "
     "The family-graph maps each genomic position → (exon-class node, column index).",
     'sans-serif'),
    ("This is data, not computation — the aligner did this work already.",
     'sans-serif'),
    "",
    "  Read r:    A   T   G   C   A   G   T          (7 bases, no indels in this example)",
    "             │   │   │   │   │   │   │",
    "             ▼   ▼   ▼   ▼   ▼   ▼   ▼",
    "  Column:    1   2   3   4   5   6   7          (in node E2 of paralog A's path)",
    "                                 ↑",
    "                               discriminating column",
]
block(0.40, 9.60, 15.20, 2.95, "B.  Where the per-base column assignments come from  (BAM + family graph)",
      stepB_lines, color='#207020', bg='#e8f5e8', line_h=0.275)

# ── Step C: Sum log-emissions → log P(r | c) ───────────────────────────
# Two side-by-side columns for paralog A and paralog B
parA_lines = [
    "log P(r | A)  =  Σ_j  log e_M(r[j] | col j, A)",
    "",
    "  log e_M(A | col 1, A) = log(0.95) = −0.051",
    "  log e_M(T | col 2, A) = log(0.94) = −0.062",
    "  log e_M(G | col 3, A) = log(0.94) = −0.062",
    "  log e_M(C | col 4, A) = log(0.94) = −0.062",
    "  log e_M(A | col 5, A) = log(0.95) = −0.051   ← match",
    "  log e_M(G | col 6, A) = log(0.94) = −0.062",
    "  log e_M(T | col 7, A) = log(0.94) = −0.062",
    "                                     ─────────",
    "  log P(r | A)                    =   −0.412",
]
block(0.40, 4.45, 7.50, 5.05, "C.  Sum log-emissions  —  paralog A",
      parA_lines, color='#a02060', bg='#fff0f5', line_h=0.295)

parB_lines = [
    "log P(r | B)  =  Σ_j  log e_M(r[j] | col j, B)",
    "",
    "  log e_M(A | col 1, B) = log(0.95) = −0.051",
    "  log e_M(T | col 2, B) = log(0.94) = −0.062",
    "  log e_M(G | col 3, B) = log(0.94) = −0.062",
    "  log e_M(C | col 4, B) = log(0.94) = −0.062",
    "  log e_M(A | col 5, B) = log(0.02) = −3.912   ← MISMATCH (B has G here)",
    "  log e_M(G | col 6, B) = log(0.94) = −0.062",
    "  log e_M(T | col 7, B) = log(0.94) = −0.062",
    "                                     ─────────",
    "  log P(r | B)                    =   −4.273",
]
block(8.10, 4.45, 7.50, 5.05, "C.  Sum log-emissions  —  paralog B",
      parB_lines, color='#a02060', bg='#fff0f5', line_h=0.295)

# ── Result: gap from one SNP ───────────────────────────────────────────
result_lines = [
    "log P(r | A) − log P(r | B)  =  −0.412 − (−4.273)  =  +3.861 nats",
    "",
    ("One paralog-distinguishing SNP at column 5 alone delivers ~4 nats of evidence in favor of paralog A.",
     'sans-serif'),
    ("Each additional discriminating column adds independently — gaps reach Δ ≥ 10 (the score-gap threshold) with as few as 3-4 SNPs.",
     'sans-serif'),
]
block(0.40, 2.50, 15.20, 1.85, "Result  —  the gap that drives EM",
      result_lines, color='#7a5500', bg='#fff7e0', line_h=0.31)

# ── Step D: Generalizing to indels (transitions also observed) ─────────
stepD_lines = [
    ("Same idea, slightly more bookkeeping.  Inserted bases get an I-state log emission + transition costs.  "
     "Deleted columns get D-state transitions only (silent — 0 bases).",
     'sans-serif'),
    "",
    "  Insertion of k bases between cols j-1 and j:   add   log a_MI + Σ log e_I(b)  +  (k−1)·log a_II  +  log a_IM",
    "  Deletion of k consecutive cols:                add   log a_MD  +  (k−1)·log a_DD  +  log a_DM",
    "",
    ("Transition probabilities a_MM, a_MI, a_MD, … are also OBSERVED — counted from indel events in the MSA.",
     'sans-serif'),
    ("e_I(b) typically uniform 0.25 (insertions are noise-driven).  No new inference; same lookup table strategy.",
     'sans-serif'),
]
block(0.40, 0.20, 15.20, 2.20, "D.  Same arithmetic for indels  (transition probabilities also come from the MSA)",
      stepD_lines, color='#5a5a82', bg='#f0f0fa', line_h=0.275)

plt.tight_layout()
plt.savefig('/scratch/jxi21/Assembler/Rustle/bench/meeting_2026_05_06_hmm_em/13_likelihood_origin.png',
            dpi=160, bbox_inches='tight')
print("13_likelihood_origin.png written")
