#!/usr/bin/env python3
"""Companion to em_derivation: (A) EM iterations converging on DAZ-like numbers,
(B) the structural-linkage channel — copy-specific junctions as a second source
of distinguishing mismatches (dNM). mathtext-safe.
"""
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch, Arc

fig = plt.figure(figsize=(12.5, 14))
ax = fig.add_axes([0, 0, 1, 1]); ax.axis('off'); ax.set_xlim(0, 1); ax.set_ylim(0, 1)

def box(x, y, w, h, fc, ec='#666'):
    ax.add_patch(FancyBboxPatch((x, y), w, h, boxstyle="round,pad=0.006", fc=fc, ec=ec, lw=1.2, zorder=1))
def txt(x, y, s, size=12, ha='left', va='top', color='black', weight='normal', mono=False, style='normal'):
    ax.text(x, y, s, size=size, ha=ha, va=va, color=color, weight=weight, style=style,
            family='monospace' if mono else 'sans-serif', zorder=3)

txt(0.5, 0.985, "EM in action, and the structural-linkage channel", 19, ha='center', weight='bold')

# ================= Panel A: iterations converging (DAZ-like) =================
box(0.04, 0.52, 0.92, 0.43, '#eef3f7')
txt(0.06, 0.94, "A.  Running the EM to convergence — a DAZ-like family (copy A = DAZ1 dominant, copy B = DAZ3 weak)",
    13, weight='bold')

# read groups table
txt(0.06, 0.905, "10 reads, three kinds:", 12, weight='bold')
hdr = ["", "$x_A$", "$x_B$", "count", "$\\gamma_{rA}$ (E-step)"]
xc = [0.06, 0.27, 0.34, 0.42, 0.52]
for x, h in zip(xc, hdr): txt(x, 0.878, h, 11.5, weight='bold', color='#333')
rows = [
    ("reads that fit DAZ1", "0", "3", "6", r"$\approx 1$"),
    ("reads that fit DAZ3", "3", "0", "1", r"$\approx 0$"),
    ("tied reads (no distinguishing pos)", "1", "1", "3", r"$=\pi_A$  (prior!)"),
]
for i, r in enumerate(rows):
    y = 0.852 - i*0.026
    for x, v in zip(xc, r):
        txt(x, y, v, 11, color=('#7570b3' if 'prior' in v else 'black'))

txt(0.06, 0.760, r"M-step: $\pi_A=\frac{1}{10}\sum_r\gamma_{rA}=\frac{6\cdot1+1\cdot0+3\cdot\pi_A}{10}$"
                 r"   $\Rightarrow$ fixed point $\pi_A=\frac{6}{7}=0.857$.", 12.5)

# iteration table
txt(0.06, 0.715, "iteration", 11.5, weight='bold', color='#333')
txt(0.20, 0.715, r"$\pi_A$ (DAZ1)", 11.5, weight='bold', color='#2c7fb8')
txt(0.35, 0.715, r"$\pi_B$ (DAZ3)", 11.5, weight='bold', color='#2ca25f')
iters = [(0, .500, .500), (1, .750, .250), (2, .825, .175), (3, .848, .152), (4, .854, .146)]
for i, (t, a, b) in enumerate(iters):
    y = 0.690 - i*0.024
    txt(0.09, y, str(t), 11, mono=True)
    txt(0.22, y, f"{a:.3f}", 11, mono=True, color='#2c7fb8')
    txt(0.37, y, f"{b:.3f}", 11, mono=True, color='#2ca25f')
txt(0.09, 0.690 - 5*0.024, "$\\infty$", 11)
txt(0.22, 0.690 - 5*0.024, "0.857", 11, mono=True, color='#2c7fb8', weight='bold')
txt(0.37, 0.690 - 5*0.024, "0.143", 11, mono=True, color='#2ca25f', weight='bold')

# convergence line plot (inset)
axp = fig.add_axes([0.60, 0.585, 0.33, 0.27])
its = list(range(7))
pa = [0.5]
for _ in range(6): pa.append((6 + 3*pa[-1]) / 10)
pb = [1 - x for x in pa]
axp.plot(its, pa, 'o-', color='#2c7fb8', label='$\\pi_A$ DAZ1', ms=4)
axp.plot(its, pb, 'o-', color='#2ca25f', label='$\\pi_B$ DAZ3', ms=4)
axp.axhline(6/7, ls='--', color='#2c7fb8', lw=0.8); axp.axhline(1/7, ls='--', color='#2ca25f', lw=0.8)
axp.set_xlabel('EM iteration', size=9); axp.set_ylabel('abundance $\\pi$', size=9)
axp.set_ylim(0, 1); axp.tick_params(labelsize=8); axp.legend(fontsize=8, loc='center right')
axp.set_title('converges in a few iterations', size=9)

txt(0.06, 0.545, "Honest behavior: DAZ3 settles at ~14% — its 1 genuine read plus the prior-share of the 3 tied reads. "
                 "The ties are APPORTIONED by $\\pi$, never fabricated into a confident split.",
    11.5, color='#a8480a', style='italic')

# ================= Panel B: structural-linkage channel =================
box(0.04, 0.05, 0.92, 0.43, '#f5f0fa')
txt(0.06, 0.47, "B.  The structural-linkage channel — a copy-specific JUNCTION is also a distinguishing position",
    13, weight='bold')
txt(0.06, 0.442, "When copies share their SNPs but differ in STRUCTURE, the splice junctions carry the signal. "
                 "Treat a junction exactly like a base:", 11.5)

# mini splice graphs
def exon(x, y, lab, col):
    ax.add_patch(plt.Rectangle((x, y), 0.05, 0.03, fc='#dfe6ec', ec=col, lw=1.5, zorder=2))
    txt(x+0.025, y+0.022, lab, 10, ha='center', mono=True)
# copy A: E1-E2-E4 (skips E3)
txt(0.15, 0.405, "copy A :", 11, ha='right', mono=True, color='#2c7fb8', weight='bold')
for i, (x, lab) in enumerate([(0.17, 'E1'), (0.30, 'E2'), (0.56, 'E4')]):
    exon(x, 0.378, lab, '#2c7fb8')
ax.add_patch(FancyArrowPatch((0.22, 0.393), (0.30, 0.393), arrowstyle='-', lw=1.5, color='#2c7fb8'))
ax.add_patch(Arc((0.43, 0.408), 0.26, 0.06, theta1=0, theta2=180, color='#2c7fb8', lw=1.8))
txt(0.43, 0.452, "E2$\\rightarrow$E4 junction (skips E3)", 9.5, ha='center', color='#2c7fb8')
# copy B: E1-E2-E3-E4
txt(0.15, 0.345, "copy B :", 11, ha='right', mono=True, color='#2ca25f', weight='bold')
for i, (x, lab) in enumerate([(0.17, 'E1'), (0.30, 'E2'), (0.43, 'E3'), (0.56, 'E4')]):
    exon(x, 0.318, lab, '#2ca25f')
for (x1, x2) in [(0.22, 0.30), (0.35, 0.43), (0.48, 0.56)]:
    ax.add_patch(FancyArrowPatch((x1, 0.333), (x2, 0.333), arrowstyle='-', lw=1.5, color='#2ca25f'))
txt(0.64, 0.345, "includes E3 (no E2$\\rightarrow$E4 junction)", 9.5, color='#2ca25f', va='center')

# a read using the E2->E4 junction
txt(0.06, 0.285, "A read whose alignment uses the E2$\\rightarrow$E4 junction MATCHES copy A's structure and "
                 "MISMATCHES copy B (which has no such junction).", 11)

txt(0.06, 0.245, r"$P(r\mid c)=(1-\varepsilon)^{m_c}\,\varepsilon^{x_c}\times\prod_{j} q_j$"
                 "    (sequence term  x  one factor per junction)", 13.5)
txt(0.06, 0.213, r"$q_j=1-\varepsilon_J$  if junction $j$ is in copy $c$'s structure,  else  $\varepsilon_J$", 12)
txt(0.06, 0.178, "Each copy-specific junction adds to ΔNM exactly like a SNP:   "
                 "ΔNM(total) = ΔNM(sequence) + ΔNM(structure).", 12.5)
txt(0.06, 0.145, "A read TIED on sequence (ΔNM_seq = 0) but using a copy-specific junction "
                 "becomes DECISIVE — same EM, more signal.", 11.5, color='#2c8255')
txt(0.06, 0.108, "This is the thesis's novel channel: structure as input. It pushes the identifiability boundary "
                 "— resolving reads that SNPs alone cannot", 11.5, color='#a8480a', style='italic')
txt(0.06, 0.084, "(and it is exactly why DAZ1$-$/DAZ3$+$, near-identical in sequence but oppositely structured, "
                 "are resolvable at all).", 11.5, color='#a8480a', style='italic')

fig.savefig('bench/paralog_secondary_scan/em_derivation_part2.png', dpi=140, bbox_inches='tight')
print("wrote bench/paralog_secondary_scan/em_derivation_part2.png")
