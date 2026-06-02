#!/usr/bin/env python3
"""Worked-example derivation of the copy-assignment (fingerprint) EM.
Latent z_r = copy a read came from; parameters pi_c = per-copy abundance;
per-read likelihood = match/mismatch Bernoulli over distinguishing positions.
Renders a single self-contained PNG suitable for slides / thesis.
(mathtext-safe: no \\begin{cases}, \\text, or \\dfrac.)
"""
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch

fig = plt.figure(figsize=(12.5, 15))
ax = fig.add_axes([0, 0, 1, 1]); ax.axis('off'); ax.set_xlim(0, 1); ax.set_ylim(0, 1)

def box(x, y, w, h, fc, ec='#666'):
    ax.add_patch(FancyBboxPatch((x, y), w, h, boxstyle="round,pad=0.006",
                                fc=fc, ec=ec, lw=1.2, zorder=1))
def txt(x, y, s, size=12, ha='left', va='top', color='black', weight='normal', mono=False, style='normal'):
    ax.text(x, y, s, size=size, ha=ha, va=va, color=color, weight=weight, style=style,
            family='monospace' if mono else 'sans-serif', zorder=3)

txt(0.5, 0.985, "Deriving the copy-assignment EM (fingerprint EM)", 19, ha='center', weight='bold')
txt(0.5, 0.963, r"latent = which copy a read came from   ·   parameters = per-copy abundances $\pi_c$   ·   "
                "likelihood = match/mismatch at distinguishing positions",
    11.5, ha='center', color='#444', style='italic')

# ---- 1. The picture ----
box(0.04, 0.78, 0.92, 0.155, '#eef3f7')
txt(0.06, 0.928, "1.  The picture — two paralog copies differ at a few positions; a read is a string of bases", 13, weight='bold')
cols = [0.30, 0.37, 0.44, 0.51, 0.58]
posn = ['p1', 'p2', 'p3', 'p4', 'p5']
copyA = ['A', 'C', 'G', 'T', 'A']
copyB = ['G', 'C', 'G', 'A', 'A']     # differ at p1, p4
read  = ['A', 'C', 'G', 'T', 'A']     # matches A; 2 mismatches vs B
for c, p in zip(cols, posn): txt(c, 0.905, p, 10, ha='center', color='#888', mono=True)
txt(0.27, 0.878, "copy A :", 12, ha='right', mono=True, color='#2c7fb8', weight='bold')
txt(0.27, 0.855, "copy B :", 12, ha='right', mono=True, color='#2ca25f', weight='bold')
txt(0.27, 0.828, "read  r :", 12, ha='right', mono=True, weight='bold')
for i, c in enumerate(cols):
    dist = copyA[i] != copyB[i]
    txt(c, 0.878, copyA[i], 13, ha='center', mono=True, color='#2c7fb8', weight='bold')
    txt(c, 0.855, copyB[i], 13, ha='center', mono=True, color='#2ca25f', weight='bold')
    txt(c, 0.828, read[i], 13, ha='center', mono=True, weight='bold')
    if dist:
        ax.add_patch(plt.Rectangle((c-0.022, 0.815), 0.044, 0.105, fill=False, ec='#c0392b', lw=1.4, ls='--', zorder=2))
txt(0.62, 0.880, "distinguishing positions (red): p1, p4.", 11, color='#a8480a')
txt(0.62, 0.856, "read matches A at all 5 -> 0 mismatches to A,", 11)
txt(0.62, 0.832, r"2 mismatches to B -> $x_A=0,\ x_B=2$.", 11)

# ---- 2. Nucleotides -> probability ----
box(0.04, 0.585, 0.92, 0.18, '#f5f0fa')
txt(0.06, 0.758, "2.  Nucleotides $\\rightarrow$ probability  (this is the whole 'fingerprint')", 13, weight='bold')
txt(0.06, 0.730, r"Generative model: pick a copy $c$ with prob $\pi_c$, then copy its bases with error rate "
                 r"$\varepsilon$ (match w.p. $1-\varepsilon$, mismatch w.p. $\varepsilon$).", 11.5)
txt(0.06, 0.700, r"$P(r \mid c)=\prod_{pos} p_{pos}$,   $p_{pos}=1-\varepsilon$ if read matches $c$ there, else "
                 r"$\varepsilon$   $=(1-\varepsilon)^{m_c}\,\varepsilon^{x_c}$", 13)
txt(0.06, 0.665, r"$\log P(r\mid c)=\mathrm{const}-x_c\cdot\ln\frac{1-\varepsilon}{\varepsilon}$"
                 r"     $\Rightarrow$  more mismatches $x_c$ $\rightarrow$ exponentially less likely.", 13)
txt(0.06, 0.624, r"Worked ($\varepsilon=0.05$):   $P(r\mid A)=0.95^5=0.774$    ·    "
                 r"$P(r\mid B)=0.95^3\,0.05^2=0.00214$", 12.5, color='#333')
txt(0.06, 0.601, r"$x_c$ is the edit distance (NM) at the distinguishing positions — that is why raw $\Delta$NM is the signal.",
    11.5, color='#a8480a', style='italic')

# ---- 3. E-step ----
box(0.04, 0.40, 0.45, 0.165, '#e3eef5')
txt(0.06, 0.558, "3.  E-step — soft-assign each read (Bayes)", 12.5, weight='bold')
txt(0.06, 0.520, r"$\gamma_{rc}=P(c\mid r)=\frac{\pi_c\,P(r\mid c)}{\sum_{c'}\pi_{c'}P(r\mid c')}$", 15)
txt(0.06, 0.468, r"with $\pi_A=\pi_B=0.5$:", 11.5)
txt(0.06, 0.440, r"$\gamma_{rA}=\frac{0.5\cdot0.774}{0.5\cdot0.774+0.5\cdot0.00214}=0.997$", 12.5, color='#2c7fb8')
txt(0.06, 0.410, "read r is 99.7% assigned to copy A.", 11.5)

# ---- 4. M-step ----
box(0.51, 0.40, 0.45, 0.165, '#d7efe0')
txt(0.53, 0.558, "4.  M-step — re-estimate abundances", 12.5, weight='bold')
txt(0.53, 0.515, r"$\pi_c=\frac{1}{N}\sum_{r}\gamma_{rc}$", 15)
txt(0.53, 0.466, "(soft read counts / total). Each copy's", 11.5)
txt(0.53, 0.443, "new abundance = the responsibility", 11.5)
txt(0.53, 0.420, "mass it collected in the E-step.", 11.5)

ax.add_patch(FancyArrowPatch((0.49, 0.50), (0.51, 0.50), arrowstyle='-|>', mutation_scale=18, lw=1.6, color='#444'))
ax.add_patch(FancyArrowPatch((0.51, 0.46), (0.49, 0.46), arrowstyle='-|>', mutation_scale=18, lw=1.6, color='#444'))
txt(0.50, 0.585, "iterate", 10.5, ha='center', color='#444')

# ---- 5. Convergence + the deltaNM connection ----
box(0.04, 0.215, 0.92, 0.165, '#fff4e6')
txt(0.06, 0.372, r"5.  Iterate E $\leftrightarrow$ M until $\pi$ stops changing  —  the tie / decisive boundary", 13, weight='bold')
txt(0.06, 0.338, r"The likelihood ratio between two copies depends ONLY on the mismatch difference $\Delta=x_B-x_A$:", 11.5)
txt(0.06, 0.300, r"$\frac{P(r\mid A)}{P(r\mid B)}=\left(\frac{1-\varepsilon}{\varepsilon}\right)^{x_B-x_A}"
                 r"=\left(\frac{0.95}{0.05}\right)^{\Delta}=19^{\,\Delta}$", 15)
txt(0.06, 0.256, r"$\Delta$NM $=0$  (no distinguishing position covered) $\rightarrow$ ratio $=1$ $\rightarrow$ "
                 "50/50 = NON-IDENTIFIABLE.", 11.5, color='#7570b3')
txt(0.06, 0.232, r"$\Delta$NM $=2$ $\rightarrow$ ratio $=361$ $\rightarrow$ decisive   "
                 r"(exactly our anchor threshold $T=2$).", 11.5, color='#2c8255')

# ---- footer ----
box(0.04, 0.045, 0.92, 0.15, '#f2f4f6')
txt(0.06, 0.185, "The algorithm, in full:", 12.5, weight='bold')
steps = [
    r"  0.  initialise abundances $\pi_c$ (uniform, or from pileup depth)",
    r"  1.  E-step:  for every read, $\gamma_{rc}\propto\pi_c\,(1-\varepsilon)^{m_c}\varepsilon^{x_c}$   (soft copy assignment)",
    r"  2.  M-step:  $\pi_c\leftarrow\frac{1}{N}\sum_r\gamma_{rc}$   (re-estimate abundances)",
    r"  3.  repeat 1-2 until $\pi$ converges; reads with $\Delta$NM$=0$ stay split by the prior $\pi$ (honestly unresolved)",
]
for i, s in enumerate(steps):
    txt(0.06, 0.158 - i*0.028, s, 11.5)

fig.savefig('bench/paralog_secondary_scan/em_derivation.png', dpi=140, bbox_inches='tight')
print("wrote bench/paralog_secondary_scan/em_derivation.png")
