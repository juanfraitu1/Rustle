#!/usr/bin/env python3
"""End-to-end mechanism on real DAZ reads: BAM tags -> de events -> EM
responsibility gamma -> bundle reweight -> max-flow assembly. mathtext-safe."""
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch

fig = plt.figure(figsize=(13, 14))
ax = fig.add_axes([0, 0, 1, 1]); ax.axis('off'); ax.set_xlim(0, 1); ax.set_ylim(0, 1)
def box(x, y, w, h, fc, ec='#666'):
    ax.add_patch(FancyBboxPatch((x, y), w, h, boxstyle="round,pad=0.006", fc=fc, ec=ec, lw=1.2, zorder=1))
def txt(x, y, s, size=12, ha='left', va='top', color='black', weight='normal', mono=False, style='normal'):
    ax.text(x, y, s, size=size, ha=ha, va=va, color=color, weight=weight, style=style,
            family='monospace' if mono else 'sans-serif', zorder=3)
def arrow(x1, y1, x2, y2):
    ax.add_patch(FancyArrowPatch((x1, y1), (x2, y2), arrowstyle='-|>', mutation_scale=20, lw=1.8, color='#444'))

txt(0.5, 0.985, "Read $\\rightarrow$ copy $\\rightarrow$ bundle $\\rightarrow$ max-flow  (real DAZ reads)",
    19, ha='center', weight='bold')
txt(0.5, 0.962, "every number is from  samtools view GGO.bam NC_073248.2:42700000-43000000", 11,
    ha='center', color='#666', style='italic', mono=True)

# Stage 1: BAM tags + de events
box(0.04, 0.74, 0.92, 0.20, '#eef3f7')
txt(0.06, 0.928, "1.  BAM tags  $\\rightarrow$  events = de $\\times$ aligned-length  (gap-compressed: each indel = 1)", 13, weight='bold')
hdr = ["read", "de(DAZ1)", "de(DAZ3)", "NM(D1/D3)", "events D1/D3", "$\\Delta$Events", "call"]
xs = [0.06, 0.24, 0.36, 0.48, 0.62, 0.78, 0.89]
for x, h in zip(xs, hdr): txt(x, 0.900, h, 10.5, weight='bold', color='#333')
rows = [
    ("…46596327", "0.0004", "0.0063", "1 / 14", "1 / 14", "+13", "DAZ1", '#2c7fb8'),
    ("…53870980", "0.0278", "0.0229", "88 / 565", "80 / 54", "-26", "DAZ3", '#2ca25f'),
    ("…24709554", "0.0000", "0.0000", "0 / 0", "0 / 0", "0", "tie", '#7570b3'),
]
for i, r in enumerate(rows):
    y = 0.872 - i*0.030
    for x, v in zip(xs, r[:7]):
        txt(x, y, v, 10.5, mono=True, color=(r[7] if x in (xs[5], xs[6]) else 'black'),
            weight=('bold' if x in (xs[5], xs[6]) else 'normal'))
txt(0.06, 0.762, "Row 2 is the proof: NM says DAZ1 (88<565) but the 565 is indel slippage; de says DAZ3 (54<80 events).",
    11, color='#a8480a', style='italic')

arrow(0.50, 0.738, 0.50, 0.716)

# Stage 2: EM responsibility
box(0.04, 0.55, 0.92, 0.165, '#f5f0fa')
txt(0.06, 0.703, "2.  EM E-step:  $\\gamma_{rc}=\\mathrm{softmax}(\\log\\,\\mathrm{score}_c+\\log\\,\\mathrm{prior}_c)$    (vg.rs:2436) — soft, never hard", 13, weight='bold')
txt(0.06, 0.672, "responsibility from $\\Delta$Events (each event $\\approx$ log-odds $\\ln\\frac{1-\\varepsilon}{\\varepsilon}\\approx2.94$):", 11.5)
gam = [("…46596327", "1.000", "0.000", '#2c7fb8'), ("…53870980", "0.000", "1.000", '#2ca25f'), ("…24709554", "0.500", "0.500", '#7570b3')]
txt(0.10, 0.640, "read", 11, weight='bold'); txt(0.34, 0.640, "$\\gamma$(DAZ1)", 11, weight='bold', color='#2c7fb8'); txt(0.52, 0.640, "$\\gamma$(DAZ3)", 11, weight='bold', color='#2ca25f')
for i, (q, a, b, c) in enumerate(gam):
    y = 0.616 - i*0.024
    txt(0.10, y, q, 10.5, mono=True); txt(0.36, y, a, 10.5, mono=True, color='#2c7fb8'); txt(0.54, y, b, 10.5, mono=True, color='#2ca25f')
txt(0.66, 0.616, "tie (row 3) splits 50/50", 11, color='#7570b3')
txt(0.66, 0.594, "by the prior — apportioned,", 11, color='#7570b3')
txt(0.66, 0.572, "not invented.", 11, color='#7570b3')

arrow(0.50, 0.548, 0.50, 0.526)

# Stage 3: bundles + max-flow
box(0.04, 0.05, 0.92, 0.47, '#eef7f0')
txt(0.06, 0.508, "3.  $\\gamma$ becomes the read's weight in EACH copy's bundle (vg.rs:2480) $\\rightarrow$ splice graph $\\rightarrow$ MAX-FLOW", 13, weight='bold')
txt(0.06, 0.481, "A multimapper stays in BOTH bundles; its weight there = $\\gamma$. Coverage = $\\Sigma$ weights. Max-flow pushes the heaviest", 11)
txt(0.06, 0.461, "source$\\rightarrow$sink path = an isoform (path_extract.rs:11). Spillover (weight$\\approx$0 in sibling) adds $\\approx$0 flow.", 11)

# two copy panels
def copy_panel(x0, title, col, reads, flow, made):
    box(x0, 0.10, 0.42, 0.33, '#ffffff', ec=col)
    txt(x0+0.21, 0.415, title, 13, ha='center', weight='bold', color=col)
    txt(x0+0.02, 0.388, "reads (weight = $\\gamma$):", 10.5, color='#555')
    for i, (lab, w, ww) in enumerate(reads):
        yy = 0.366 - i*0.024
        ax.add_patch(plt.Rectangle((x0+0.03, yy-0.008), 0.16*ww+0.002, 0.013, fc=col, ec=col, alpha=0.25+0.7*ww, zorder=2))
        txt(x0+0.21, yy, f"{lab}  w={w}", 9.5, mono=True, color=('#333' if ww>0.3 else '#aaa'))
    # splice graph (3 exon nodes + flow path)
    gy = 0.20
    for j, ex in enumerate(['E1','E2','E3']):
        ax.add_patch(plt.Rectangle((x0+0.05+j*0.12, gy), 0.06, 0.03, fc='#dfe6ec', ec=col, lw=1.3, zorder=2))
        txt(x0+0.08+j*0.12, gy+0.022, ex, 9, ha='center', mono=True)
    for j in range(2):
        ax.add_patch(FancyArrowPatch((x0+0.11+j*0.12, gy+0.015), (x0+0.17+j*0.12, gy+0.015),
                                     arrowstyle='-|>', mutation_scale=11, lw=1.5+2*flow, color=col, alpha=0.4+0.6*flow))
    txt(x0+0.21, gy-0.018, "max-flow  $\\Rightarrow$  " + made, 10.5, ha='center', color=col, weight='bold')

copy_panel(0.05, "DAZ1 bundle", '#2c7fb8',
           [("…46596327", "1.0", 1.0), ("…53870980", "0.0", 0.0), ("…24709554", "0.5", 0.5), ("+ unique DAZ1", "1.0", 1.0)],
           1.0, "DAZ1 isoforms (full abundance)")
copy_panel(0.53, "DAZ3 bundle", '#2ca25f',
           [("…46596327", "0.0", 0.0), ("…53870980", "1.0", 1.0), ("…24709554", "0.5", 0.5), ("(150 spillover)", "~0", 0.05)],
           0.35, "DAZ3 at HONEST low abundance")

txt(0.5, 0.075, "Spillover reads that ALIGN to DAZ3 but BELONG to DAZ1 carry weight ~0 here -> no phantom DAZ3 isoform.",
    11, ha='center', color='#a8480a', style='italic')

fig.savefig('bench/paralog_secondary_scan/mechanism_pipeline.png', dpi=135, bbox_inches='tight')
print("wrote bench/paralog_secondary_scan/mechanism_pipeline.png")
