#!/usr/bin/env python3
"""Slide 12 — Worked EM example with concrete numbers.

Three paralogs (A, B, C), four reads, Δ = 10 nats, ε = 1e-3. Every
arithmetic step is shown so the advisor can verify by hand:
initialization, M-step priors, E-step softmax, score-gap abstention,
and convergence in two iterations."""
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch

fig, ax = plt.subplots(figsize=(16, 22))
ax.set_xlim(0, 16)
ax.set_ylim(0, 22)
ax.axis('off')

# ── Title ───────────────────────────────────────────────────────────────
ax.text(8, 21.55, "Slide 12  —  Worked EM example: numbers all the way down",
        fontsize=17, fontweight='bold', ha='center')
ax.text(8, 21.10,
        "3 paralogs (A, B, C),  4 reads,  Δ = 10 nats,  ε = 1e-3.   "
        "Every arithmetic step shown — verify by hand if you like.",
        fontsize=10.5, ha='center', color='#444', fontstyle='italic')


def block(x, y, w, h, title, lines, color='#5a82bc', bg='#eef6ff',
          font=9.5, line_h=0.30, title_font=12, top_pad=0.85):
    """Render a labeled block. `lines` is a list of strings (monospace
    by default) or (text, family) tuples for non-monospace lines."""
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


# ── 1. Setup  (the data going in) ──────────────────────────────────────
setup_lines = [
    ("Per-paralog likelihoods log P(r | c) come from forward_against_path_for_copy (slide 5).",
     'sans-serif'),
    ("Computed once up front from the HMM trellis; constant across all EM iterations.",
     'sans-serif'),
    "",
    "  read    placements    log P(r|A)   log P(r|B)   log P(r|C)   gap = best − 2nd",
    "  ────    ───────────   ──────────   ──────────   ──────────   ────────────────",
    "  r1      {A, B, C}        −100         −103         −105         3   < Δ → ABSTAIN",
    "  r2      {A, B}            −95         −120          —          25",
    "  r3      {A, C}           −110          —            −95        15",
    "  r4      {A, B, C}        −118         −100         −125        18",
]
block(0.40, 17.85, 15.20, 3.00, "Setup  —  the data going in",
      setup_lines, color='#444', bg='#fafafa', line_h=0.275)

# ── 2. Initialization (t = 0) ──────────────────────────────────────────
init_lines = [
    "w_{r,c}^(0) = 1 / |P_r|   for c ∈ P_r,   else 0       ← BAM aligner's 1/NH",
    "",
    "  r1:  w_{1,A}^(0) = w_{1,B}^(0) = w_{1,C}^(0) = 1/3 ≈ 0.333",
    "  r2:  w_{2,A}^(0) = w_{2,B}^(0)               = 1/2 = 0.500",
    "  r3:  w_{3,A}^(0) = w_{3,C}^(0)               = 1/2 = 0.500",
    "  r4:  w_{4,A}^(0) = w_{4,B}^(0) = w_{4,C}^(0) = 1/3 ≈ 0.333",
]
block(0.40, 14.85, 15.20, 2.90, "1.  Initialization  (t = 0)",
      init_lines, color='#1f77b4', bg='#eef6ff', line_h=0.30)

# ── 3a. M-step at t = 1  (left half) ───────────────────────────────────
mstep_lines = [
    "n_c = Σ_r  w_{r,c}^(0)      π_c^(1) = n_c / Σ n_{c'} + ε",
    "",
    "n_A = 0.333 + 0.500 + 0.500 + 0.333 = 1.667",
    "n_B = 0.333 + 0.500 +   0   + 0.333 = 1.167",
    "n_C = 0.333 +   0   + 0.500 + 0.333 = 1.167",
    "Σ n_c = 4.000",
    "",
    "π_A^(1) = 1.667 / 4.000 + ε = 0.417",
    "π_B^(1) = 1.167 / 4.000 + ε = 0.292",
    "π_C^(1) = 1.167 / 4.000 + ε = 0.292",
    "",
    "log π_A^(1) = ln(0.417) = −0.875",
    "log π_B^(1) = ln(0.292) = −1.232",
    "log π_C^(1) = ln(0.292) = −1.232",
]
block(0.40, 9.90, 7.50, 4.85, "2a.  M-step at t = 1",
      mstep_lines, color='#207020', bg='#e8f5e8', line_h=0.27)

# ── 3b. E-step at t = 1  (right half) ──────────────────────────────────
estep_lines = [
    "log Q_{r,c} = log P(r|c) + log π_c^(1)        softmax over c ∈ P_r",
    "",
    "r1:  gap = −100 − (−103) = 3 < Δ → ABSTAIN",
    "     w_{1,*}^(1) = (0.333, 0.333, 0.333)   [unchanged]",
    "",
    "r2:  gap = 25 ≥ Δ → update",
    "     log Q_{2,A} = −95  + (−0.875) = −95.875",
    "     log Q_{2,B} = −120 + (−1.232) = −121.232",
    "     subtract max:  0.000, −25.357",
    "     softmax:  w_{2,A}^(1) ≈ 1.000,  w_{2,B}^(1) ≈ 0.000",
    "",
    "r3:  gap = 15 ≥ Δ →  log Q = (−110.875, —, −96.232)",
    "     softmax:  w_{3,A}^(1) ≈ 0.000,  w_{3,C}^(1) ≈ 1.000",
    "",
    "r4:  gap = 18 ≥ Δ →  log Q = (−118.875, −101.232, −126.232)",
    "     softmax:  w_{4,A}^(1) ≈ 0,  w_{4,B}^(1) ≈ 1,  w_{4,C}^(1) ≈ 0",
]
block(8.10, 9.90, 7.50, 4.85, "2b.  E-step at t = 1  (gap rule + softmax)",
      estep_lines, color='#7a5500', bg='#fff7e0', line_h=0.245)

# ── 3c. Convergence check at t = 1 ─────────────────────────────────────
conv1_lines = [
    "Δ_1 = max_{r,c} |w_{r,c}^(1) − w_{r,c}^(0)|",
    "",
    "    = max(  0.000,                       // r1 abstained",
    "            |1.000 − 0.500| = 0.500,     // r2",
    "            |1.000 − 0.500| = 0.500,     // r3",
    "            |1.000 − 0.333| = 0.667 )    // r4",
    "",
    "    = 0.667   >>  0.001     →   continue to t = 2",
]
block(0.40, 6.85, 15.20, 2.95, "2c.  Convergence check at t = 1",
      conv1_lines, color='#a04040', bg='#fde0e0', line_h=0.27)

# ── 4. Iteration t = 2  (M-step + E-step + convergence) ────────────────
t2_lines = [
    "M-step at t = 2:    π_c^(2) = n_c / 4.000 + ε   from w^(1)",
    "  n_A = 0.333 + 1.000 + 0.000 + 0.000 = 1.333    →   π_A^(2) = 0.334",
    "  n_B = 0.333 + 0.000 +   0   + 1.000 = 1.333    →   π_B^(2) = 0.334     (uniform)",
    "  n_C = 0.333 +   0   + 1.000 + 0.000 = 1.333    →   π_C^(2) = 0.334",
    "",
    "E-step at t = 2:    forward log P(r | c) is constant  →  gaps unchanged.",
    "  Equal log π_c cancels in the softmax difference  →  same w^(2) as w^(1).",
    "",
    "Convergence:   Δ_2 = max |w^(2) − w^(1)|  ≈  0   <  0.001    →   STOP",
]
block(0.40, 3.30, 15.20, 3.45, "3.  Iteration t = 2  (M-step + E-step + convergence)",
      t2_lines, color='#5a5a82', bg='#f0f0fa', line_h=0.275)

# ── 5. Final result ────────────────────────────────────────────────────
final_lines = [
    "  read    placements   w_A      w_B      w_C       outcome",
    "  ────    ──────────   ──────   ──────   ──────    ──────────────────────────────────",
    "  r1      {A,B,C}       0.333    0.333    0.333    abstained (gap = 3 < Δ) — keep BAM 1/NH",
    "  r2      {A,B}         1.000    0.000     —       → A   (gap = 25 nats)",
    "  r3      {A,C}         0.000     —       1.000    → C   (gap = 15 nats)",
    "  r4      {A,B,C}       0.000    1.000    0.000    → B   (gap = 18 nats)",
    "",
    ("Final priors  π = (0.334, 0.334, 0.334) — uniform, because each non-abstained read voted for a different copy.",
     'sans-serif'),
    ("With biased data the priors shift accordingly; slide 9 (AMY) is a real-data run with non-uniform converged priors.",
     'sans-serif'),
]
block(0.40, 0.20, 15.20, 3.10, "Result  —  EM converged in 2 iterations",
      final_lines, color='#a02060', bg='#fff0f5', line_h=0.265)

plt.tight_layout()
plt.savefig('/scratch/jxi21/Assembler/Rustle/bench/meeting_2026_05_06_hmm_em/12_em_worked_example.png',
            dpi=160, bbox_inches='tight')
print("12_em_worked_example.png written")
