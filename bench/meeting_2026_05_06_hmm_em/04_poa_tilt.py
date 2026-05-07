#!/usr/bin/env python3
"""POA-MSA tilt: how the per-copy profile blends the singleton (one-hot at the
copy's own base) with the column-wise MSA observation. tilt α = 1 → pure
singleton (sharp, harsh penalty per error). α = 0 → pure shared MSA
(information-collapse). α = 0.7 (default) → keeps paralog-distinguishing
signal AND softens the singleton's 6.9 log-unit per-base penalty."""
import matplotlib.pyplot as plt
import numpy as np

fig, axes = plt.subplots(1, 3, figsize=(15, 5.5))
fig.suptitle("Per-copy profile = α·(this copy's base) + (1−α)·(MSA column)",
             fontsize=16, fontweight='bold', y=1.02)

# Same toy column as in slide 2: paralogs A,B,C with bases C,G,G respectively.
# Target = paralog A → singleton emits C with prob 1-eps, others eps/3.
# MSA column observation = {C: 1/3, G: 2/3, A=T=0}.
bases = ['A', 'C', 'G', 'T']
copy_colors = ['#1f77b4', '#ff7f0e', '#2ca02c']

eps = 0.02
singleton = np.array([eps/3, 1.0 - eps, eps/3, eps/3])  # 'C' = idx 1
msa_col   = np.array([0.00, 1/3, 2/3, 0.00])

# log-unit conversion (natural log) for axis #2
def to_log(p):
    p = np.clip(p, 1e-10, 1.0)
    return np.log(p)

alphas = [1.0, 0.7, 0.0]
titles = ["α = 1.0  (pure singleton)",
          "α = 0.7  (default — POA-smoothed)",
          "α = 0.0  (pure consensus — INFO COLLAPSE)"]
title_colors = ['#1a4f8a', '#207020', '#a02020']
descs = [
    "harsh: P(G | A) = 0.007  →  log-penalty ≈ 5.0 nats per error",
    "balanced: P(G | A) = 0.227  →  log-penalty ≈ 1.5 nats per error\n"
    "still concentrates probability on 'C' (paralog-A's base)",
    "averaged: P(C | A) = P(G | A) = identical for ALL paralogs\n"
    "→ NO discrimination — exactly what we want to avoid",
]

for ax, a, t, tc, desc in zip(axes, alphas, titles, title_colors, descs):
    blended = a * singleton + (1 - a) * msa_col
    blended = blended / blended.sum()  # renormalize
    bars = ax.bar(bases, blended, color=['#888'] * 4, edgecolor='black', lw=0.7)
    # paralog A's actual base: highlight col index 1 ('C')
    bars[1].set_color('#1f77b4')
    bars[1].set_label("paralog A's actual base")
    # MSA-supported alt: 'G'
    bars[2].set_color('#ff7f0e')
    bars[2].set_label("MSA-supported alt")

    ax.set_ylim(0, 1.0)
    ax.set_ylabel("emission probability  P(base | profile, col 2)")
    ax.set_title(t, fontsize=12, fontweight='bold', color=tc)
    ax.set_xlabel(desc, fontsize=9.5, color='#444')
    ax.grid(axis='y', alpha=0.3)
    # value labels
    for b, v in zip(bars, blended):
        if v > 0.01:
            ax.text(b.get_x() + b.get_width()/2, v + 0.02, f"{v:.2f}",
                    ha='center', fontsize=9)
    ax.legend(loc='upper right', fontsize=8)

# bottom note
fig.text(0.5, -0.04,
         "α controls the trade-off between per-copy detail (singleton) and family-wide regularization (POA-MSA).\n"
         "rustle's HMM-EM uses α = 0.7 by default — measured to recover medium-similarity paralogs (jaccard 30–60%) "
         "without losing the per-copy SNP signal that pure consensus discards.",
         ha='center', fontsize=10.5, color='#333')

plt.tight_layout(rect=[0, 0.02, 1, 0.95])
plt.savefig('/scratch/jxi21/Assembler/Rustle/bench/meeting_2026_05_06_hmm_em/04_poa_tilt.png',
            dpi=160, bbox_inches='tight')
print("04_poa_tilt.png written")
