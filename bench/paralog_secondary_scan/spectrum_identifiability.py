#!/usr/bin/env python3
"""Spectrum-of-identifiability figure. Rebuts "one different SNP makes it trivial,
existing tools already solve it" by showing the trivial case is a SINGLE row at
the top while real expressed paralog families live in the rows below it, where
existing tools fail. mathtext-safe. Output: spectrum_identifiability.png

Layout is budgeted: every cell has at most 3 body lines on a fixed 0.020 grid,
plus an optional motif box, so nothing overlaps.
"""
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, Rectangle, FancyArrowPatch

fig = plt.figure(figsize=(16, 13))
ax = fig.add_axes([0, 0, 1, 1]); ax.axis('off'); ax.set_xlim(0, 1); ax.set_ylim(0, 1)

GREEN, RED, BLUE, GREY = '#1a9850', '#d73027', '#2c7fb8', '#555'

def box(x, y, w, h, fc, ec='#888', lw=1.0, z=1):
    ax.add_patch(FancyBboxPatch((x, y), w, h, boxstyle="round,pad=0.003",
                                fc=fc, ec=ec, lw=lw, zorder=z))
def rect(x, y, w, h, fc, ec='none', lw=0, z=1):
    ax.add_patch(Rectangle((x, y), w, h, fc=fc, ec=ec, lw=lw, zorder=z))
def txt(x, y, s, size=12, ha='left', va='top', color='black', weight='normal',
        mono=False, style='normal', rot=0):
    ax.text(x, y, s, size=size, ha=ha, va=va, color=color, weight=weight, style=style,
            family='monospace' if mono else 'sans-serif', rotation=rot, zorder=4)

LH = 0.020   # body line height (figure fraction)

# ---------------------------------------------------------------- title
txt(0.5, 0.987, "Read-to-copy assignment is a spectrum, not a single-SNP lookup",
    23, ha='center', weight='bold')
txt(0.5, 0.963,
    "\"One different SNP and it's trivial\" is TRUE -- for the top row only. "
    "Real expressed paralog families live in the rows below it.",
    13.5, ha='center', color=GREY, style='italic')

# ---------------------------------------------------------------- columns
c1x, w1 = 0.035, 0.140
c2x, w2 = 0.185, 0.260
c3x, w3 = 0.455, 0.180
c4x, w4 = 0.645, 0.200
c5x, w5 = 0.855, 0.130
HDR_Y = 0.930
for x, lab in [(c1x, "Regime"), (c2x, "What the reads actually look like"),
               (c3x, "Existing tools"), (c4x, "This method"), (c5x, "Output")]:
    txt(x, HDR_Y, lab, 13, weight='bold', color='#222')
rect(0.03, 0.918, 0.955, 0.002, '#999')

# left "increasing difficulty" arrow
ax.add_patch(FancyArrowPatch((0.016, 0.89), (0.016, 0.16), arrowstyle='-|>',
             mutation_scale=15, lw=2.0, color='#aaa', zorder=2))
txt(0.009, 0.52, "increasing difficulty", 11, ha='center', va='center',
    color='#888', rot=90, style='italic')

# ---------------------------------------------------------------- rows
rows = [
    dict(diff='TRIVIAL', dc=GREEN, letter='A', name=["One clean", "fixed SNP"],
         sit=["A fixed difference between the copies,",
              "and the read covers it."],
         motif="C1 ..A..   read = A\nC2 ..G..   -> copy 1",
         ex_ok=True,
         ex=["Any variant-aware assigner",
             "or phaser does this.",
             "We claim nothing here."],
         me=["Same answer --", "not what the tool is for."],
         out="hard call,\nhigh confidence"),

    dict(diff='HARD', dc='#fdae61', letter='B', name=["Sparse &", "fragile markers"],
         sit=["~1 marker per kb, each a single base",
              "= fragile at ~1% error; a read may",
              "cover few markers, or an erroneous one."],
         motif="..|.....|......|..\n(few, weak, error-prone)",
         ex_ok=False,
         ex=["Single-base call is un-",
             "calibrated; short-read tools",
             "lose linkage across markers."],
         me=["Phase MANY markers per long",
             "read -> calibrated confidence",
             "+ self-detected bad reads."],
         out="calibrated call\n+ confidence"),

    dict(diff='HARDER', dc='#f16913', letter='C', name=["The 'SNP'", "assumption", "breaks"],
         sit=["Markers are INDELS in repeat arrays",
              "(DAZ) where the sequencer errs; or a",
              "'marker' is a within-copy het (PSV vs SNP)."],
         motif="AAAA / AAA  (indel)\nbase X: PSV or het?",
         ex_ok=False,
         ex=["Identity over-penalizes indels;",
             "trusting annotated SNPs assigns",
             "on NON-markers. Calling PSVs",
             "is itself the open problem."],
         me=["Error-model-weighted sites +",
             "array length profile; discover",
             "& VALIDATE true markers by",
             "cross-read linkage."],
         out="validated-\nmarker call"),

    dict(diff='HARDER', dc='#d73027', letter='D', name=["Partial", "separability", "(>2 copies)"],
         sit=["A marker splits copy 1 from {2,3}",
              "but cannot separate 2 from 3."],
         motif="{ C1 | C2 , C3 }\n(class-resolvable only)",
         ex_ok=False,
         ex=["Force one best copy ->",
             "CONFIDENT misassignment",
             "among the indistinguishable."],
         me=["Identifiability classes:",
             "assign the class,",
             "apportion within it."],
         out="class call +\nwithin-class\nfraction"),

    dict(diff='IMPOSSIBLE\nfrom read', dc='#555', letter='E', name=["Non-identifiable", "or absent copy"],
         sit=["Read spans 0 distinguishing positions,",
              "OR the copy is absent from the reference",
              "(no marker can even be defined)."],
         motif="======= (no markers)\ncopy NOT in reference",
         ex_ok=False,
         ex=["Drop the read or dump it on",
             "the primary -> coverage on the",
             "WRONG copy; ref-only tools",
             "cannot see an absent copy."],
         me=["EM apportionment by resolvable-",
             "read expression (no fabrication);",
             "de novo discovery from read",
             "variant-linkage."],
         out="fractional weight\n/ discovered copy"),
]

top = 0.905
row_h = 0.150
gap = 0.012
for i, r in enumerate(rows):
    y1 = top - i * row_h          # row top
    h = row_h - gap
    y0 = y1 - h                   # row bottom
    rect(0.03, y0, 0.955, h, '#fbfbfb' if i % 2 == 0 else '#eef2f5', z=0.5)
    body_top = y1 - 0.014

    # ---- col1: difficulty pill + letter + wrapped name
    box(c1x, y0 + 0.010, w1, h - 0.020, '#ffffff', ec=r['dc'], lw=1.8)
    rect(c1x, y1 - 0.030, w1, 0.030, r['dc'], z=1.2)
    txt(c1x + w1/2, y1 - 0.015, r['diff'], 10, ha='center', va='center',
        color='white', weight='bold')
    txt(c1x + 0.012, y1 - 0.052, r['letter'], 22, color=r['dc'], weight='bold')
    for k, nl in enumerate(r['name']):
        txt(c1x + 0.050, y1 - 0.050 - k*0.020, nl, 11.5, weight='bold', color='#222')

    # ---- col2: situation lines + motif box
    for k, line in enumerate(r['sit']):
        txt(c2x, body_top - k*LH, line, 10.5, color='#333')
    motif_top = body_top - len(r['sit'])*LH - 0.008
    box(c2x, motif_top - 0.038, w2 - 0.012, 0.040, '#eef3f7', ec='#bcd', lw=1.0)
    txt(c2x + 0.010, motif_top - 0.006, r['motif'], 9.5, mono=True, color='#225')

    # ---- col3: existing tools (color-coded = the rhetorical wall)
    cell_fc = '#e6f4ea' if r['ex_ok'] else '#fdecea'
    cell_ec = GREEN if r['ex_ok'] else RED
    box(c3x, y0 + 0.010, w3 - 0.010, h - 0.020, cell_fc, ec=cell_ec, lw=1.5)
    txt(c3x + 0.012, body_top, ("✓ SOLVED" if r['ex_ok'] else "✗ FAILS"),
        12, weight='bold', color=cell_ec)
    for k, line in enumerate(r['ex']):
        txt(c3x + 0.012, body_top - 0.026 - k*LH, line, 9.5, color='#333')

    # ---- col4: this method
    for k, line in enumerate(r['me']):
        txt(c4x, body_top - k*LH, line, 10, color=BLUE,
            weight=('bold' if i == 0 else 'normal'))

    # ---- col5: output
    txt(c5x, body_top, r['out'], 10.5, color='#222', weight='bold')

# ---------------------------------------------------------------- bottom banner
by0 = top - len(rows)*row_h - 0.006
box(0.03, by0 - 0.090, 0.955, 0.092, '#eef0f7', ec=BLUE, lw=1.8)
txt(0.048, by0 - 0.012,
    "The \"existing tools\" column is one green row, then a wall of red.",
    14.5, weight='bold', color='#222')
txt(0.048, by0 - 0.040,
    "Real expressed primate paralog families -- DAZ repeat arrays, NBPF, segmental dups at 99-99.9% identity -- live in rows B-E, NOT in A.",
    12.5, color='#333')
txt(0.048, by0 - 0.066,
    "The contribution is ONE principled, calibrated treatment across the whole spectrum: confident where it can be, fractional where it must be, "
    "and it never fabricates a copy.",
    12.5, color=BLUE, weight='bold')

fig.savefig("spectrum_identifiability.png", dpi=150, bbox_inches='tight')
print("wrote spectrum_identifiability.png")
