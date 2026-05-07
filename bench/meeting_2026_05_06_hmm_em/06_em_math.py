#!/usr/bin/env python3
"""EM math — exact equations for priors and posteriors, tied directly to the
code path in src/rustle/vg.rs (run_pre_assembly_em_hmm, lines 2587–2647).

This is the slide that answers "where do the priors come from?" precisely.
No hand-waving — every term is named, every constant is exposed, and the
score-gap rule abstention is shown as part of the same loop."""
import matplotlib.pyplot as plt
import matplotlib.patches as patches
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch
import numpy as np

fig, ax = plt.subplots(figsize=(15.5, 9.8))
ax.set_xlim(0, 15.5)
ax.set_ylim(0, 11)
ax.axis('off')

ax.text(7.75, 10.55, "Where do the priors come from?  (standard EM, exactly)",
        fontsize=17, fontweight='bold', ha='center')
ax.text(7.75, 10.1,
        "Notation: read r has placements P_r ⊂ {1…C}; w_{r,c} is read r's weight on copy c "
        "at the current EM iteration; π_c is copy c's prior probability.",
        fontsize=10.5, ha='center', color='#444')

# Helper to draw a labeled equation block
def block(ax, x, y, w, h, title, body, color='#5a82bc', bg='#eef6ff'):
    box = FancyBboxPatch((x, y), w, h,
                         boxstyle="round,pad=0.10",
                         facecolor=bg, edgecolor=color, lw=1.2)
    ax.add_patch(box)
    ax.text(x + 0.20, y + h - 0.30, title,
            fontsize=11.5, fontweight='bold', color=color, va='top')
    ax.text(x + 0.20, y + h - 0.80, body,
            fontsize=9.5, va='top', family='monospace', color='#222')

# ── 1. Initialization ────────────────────────────────────────────────────────
init_body = (
    "for each read r and each placement c ∈ P_r:\n"
    "    w_{r,c}^(0)  =  1 / |P_r|        ← BAM aligner's 1/NH\n"
    "    w_{r,c'}     =  0   for c' ∉ P_r\n\n"
    "If a read maps to 3 paralogs, each gets weight 1/3 to start."
)
block(ax, 0.4, 7.55, 7.0, 2.10,
      "1.  Initialization  (t = 0)",
      init_body, color='#1f77b4', bg='#eef6ff')

# ── 2. M-step ────────────────────────────────────────────────────────────────
mstep_body = (
    "for each copy c:\n"
    "    n_c     =  Σ_r  w_{r,c}^(t-1)        ← sum of weights at copy c\n"
    "    π_c^(t) =  n_c / Σ_{c'} n_{c'}  + ε   ← normalize, add floor ε = 1e-3\n"
    "    log π_c^(t) = ln( π_c^(t) )\n\n"
    "Floor keeps log-prior finite when a copy has no current support.\n"
    "Code: vg.rs:2589–2598."
)
block(ax, 8.0, 7.55, 7.1, 2.10,
      "2.  M-step  (re-estimate copy priors)",
      mstep_body, color='#207020', bg='#e8f5e8')

# ── 3. Forward DP recap ──────────────────────────────────────────────────────
fwd_body = (
    "for each read r and each placement c ∈ P_r:\n"
    "    log P(r | c)  =  forward_against_path_for_copy(\n"
    "                          family_graph, read_r, path_c, copy_id = c )\n"
    "                    ← chain of profile HMM trellises along c's path  (slide 5)\n\n"
    "Computed ONCE up-front in PHASE 2 (parallel via rayon), cached for all EM iterations."
)
block(ax, 0.4, 5.3, 14.7, 2.05,
      "3.  Per-paralog likelihood from the HMM  (constant across EM iterations)",
      fwd_body, color='#a02060', bg='#fff0f5')

# ── 4. E-step (posterior + softmax) ─────────────────────────────────────────
estep_body = (
    "for each read r and each placement c ∈ P_r:\n"
    "    log Q_{r,c}  =  log P(r | c)  +  log π_c^(t)       ← un-normalized log-posterior\n"
    "                                                        (slide 5)        (M-step)\n\n"
    "Softmax over c ∈ P_r (other placements ignored — paralog-restricted EM):\n"
    "    w_{r,c}^(t)  =  exp( log Q_{r,c} − maxᶜ log Q_{r,c} )  /  Σ_{c'} exp( · )\n\n"
    "Code: vg.rs:2603–2640."
)
block(ax, 0.4, 2.65, 9.1, 2.40,
      "4.  E-step  (posterior weight per placement)",
      estep_body, color='#7a5500', bg='#fff7e0')

# ── 5. Score-gap rule (abstention) ───────────────────────────────────────────
gap_body = (
    "Before applying the softmax, check confidence:\n\n"
    "  best   = maxᶜ  log P(r | c)\n"
    "  second = 2nd-max\n\n"
    "  if (best − second) < Δ:   skip update\n"
    "                           (keep w_{r,c}^(t-1))\n\n"
    "Default Δ = 10 nats (env: RUSTLE_VG_EM_SCORE_GAP).\n"
    "Code: vg.rs:2614–2624."
)
block(ax, 9.7, 2.65, 5.4, 2.40,
      "5.  Score-gap rule  (don't gamble)",
      gap_body, color='#a04040', bg='#fde0e0')

# ── 6. Convergence ───────────────────────────────────────────────────────────
conv_body = (
    "After each iteration t:\n"
    "    Δ_t  =  max_{r, c ∈ P_r}  | w_{r,c}^(t)  −  w_{r,c}^(t-1) |\n\n"
    "Stop when  Δ_t < 0.001  or  t = max_iter (default 5).\n"
    "Most families converge in 2–3 iterations on the AMY/NBPF demo."
)
block(ax, 0.4, 0.20, 14.7, 2.05,
      "6.  Convergence",
      conv_body, color='#5a5a82', bg='#f0f0fa')

# Outer "EM loop" arrow — drawn in the empty margin, NOT through any block.
loop_arrow = FancyArrowPatch((15.20, 4.0), (15.20, 8.5),
                             arrowstyle='->', color='#7a5500', lw=1.6,
                             connectionstyle='arc3,rad=0.0')
ax.add_patch(loop_arrow)
ax.text(15.35, 6.3, "loop t → t+1",
        fontsize=9, color='#7a5500', fontweight='bold',
        ha='left', va='center', rotation=90)

plt.tight_layout()
plt.savefig('/scratch/jxi21/Assembler/Rustle/bench/meeting_2026_05_06_hmm_em/06_em_math.png',
            dpi=160, bbox_inches='tight')
print("06_em_math.png written")
