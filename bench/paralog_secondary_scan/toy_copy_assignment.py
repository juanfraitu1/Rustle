#!/usr/bin/env python3
"""Toy worked example: assigning a read to one of two near-identical paralog
copies when the ALIGNMENT SCORE TIES. Shows (1) why a whole-read scalar (AS/NM/de)
ties, (2) the reframe to copy-distinguishing positions only, (3) how a family
graph supplies those positions and reads the alleles in one pass, (4) the
likelihood-ratio math, (5) the honest non-identifiable caveat.

All numbers are computed here, not hand-typed.  mathtext-safe.
Output: toy_copy_assignment.png
"""
import math
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch, Rectangle

# ----------------------------------------------------------------------------
# MODEL NUMBERS (the realistic regime: long read, near-identical copies)
# ----------------------------------------------------------------------------
L        = 6000      # read length (bp)
eps      = 0.01      # per-base sequencing error rate (1%)
copy_div = 0.001     # divergence between the two copies (99.9% identical seg-dup)
n_diag   = round(L * copy_div)          # diagnostic sites the read spans  -> 6
n_err    = round(L * eps)               # sequencing errors per alignment  -> 60
nm_A     = n_err                        # read is from copy A: errors only
nm_B     = n_err + n_diag               # + every diagnostic site mismatches B
id_A     = 100 * (1 - nm_A / L)
id_B     = 100 * (1 - nm_B / L)
w_site   = math.log((1 - eps) / (eps / 3))      # nats per decisive diagnostic site
logLR    = n_diag * w_site
# posterior P(A) = 1/(1+e^-logLR); report the complement for readability
p_not_A  = 1.0 / (1.0 + math.exp(logLR))

# ----------------------------------------------------------------------------
# ZOOM WINDOW (a 35-bp slice of the read; shows 2 of the 6 diagnostic sites
# plus 2 sequencing errors so the contrast is visible)
# ----------------------------------------------------------------------------
COPY_A = list("AGCTTACGATGCAACGTTGCAAGTCCATGGCATAC")
W      = len(COPY_A)
DIAG   = {8: ('A', 'T'), 23: ('G', 'C')}     # col -> (A-allele, B-allele)
ERR    = {3: 'T', 16: 'A'}                    # col -> erroneous read base
COPY_B = COPY_A[:]
for c, (a, b) in DIAG.items():
    COPY_A[c] = a
    COPY_B[c] = b
READ = COPY_A[:]                              # read truly comes from copy A
for c, base in ERR.items():
    READ[c] = base                            # inject sequencing errors

# ----------------------------------------------------------------------------
# CANVAS
# ----------------------------------------------------------------------------
fig = plt.figure(figsize=(14, 17))
ax = fig.add_axes([0, 0, 1, 1]); ax.axis('off'); ax.set_xlim(0, 1); ax.set_ylim(0, 1)

BLUE, GREEN, RED, PURP, GREY = '#2c7fb8', '#2ca25f', '#d7301f', '#7570b3', '#666'

def box(x, y, w, h, fc, ec=GREY, lw=1.2):
    ax.add_patch(FancyBboxPatch((x, y), w, h, boxstyle="round,pad=0.004",
                                fc=fc, ec=ec, lw=lw, zorder=1))
def txt(x, y, s, size=12, ha='left', va='top', color='black', weight='normal',
        mono=False, style='normal'):
    ax.text(x, y, s, size=size, ha=ha, va=va, color=color, weight=weight, style=style,
            family='monospace' if mono else 'sans-serif', zorder=4)
def arrow(x1, y1, x2, y2, color='#444', lw=1.8):
    ax.add_patch(FancyArrowPatch((x1, y1), (x2, y2), arrowstyle='-|>',
                                 mutation_scale=18, lw=lw, color=color, zorder=3))

# sequence-row drawing -------------------------------------------------------
SEQ_X0, SEQ_X1 = 0.085, 0.945
dx = (SEQ_X1 - SEQ_X0) / W
def colx(i): return SEQ_X0 + dx * (i + 0.5)

def draw_seq_row(y, seq, label, kind):
    """kind: 'A','B','read'. Colors diagnostic + error cells."""
    txt(0.075, y, label, 10.5, ha='right', va='center', weight='bold', color='#222')
    for i, ch in enumerate(seq):
        cellc, ec, tc, wt = 'white', 'none', 'black', 'normal'
        if i in DIAG:
            cellc = '#eef3f7'
            if kind == 'A': tc, wt = BLUE, 'bold'
            elif kind == 'B': tc, wt = GREEN, 'bold'
            else:  # read: color by which copy it matches
                tc = BLUE if ch == DIAG[i][0] else (GREEN if ch == DIAG[i][1] else RED)
                wt = 'bold'
        if kind == 'read' and i in ERR:
            tc, wt = RED, 'bold'
        if cellc != 'white':
            ax.add_patch(Rectangle((colx(i) - dx*0.5, y - 0.011), dx, 0.022,
                                   fc=cellc, ec='none', zorder=1.5))
        txt(colx(i), y, ch, 11.5, ha='center', va='center', color=tc, weight=wt, mono=True)

# =====================================================================
# TITLE
# =====================================================================
txt(0.5, 0.985, "Assigning a read to a paralog copy when the alignment score ties",
    21, ha='center', weight='bold')
txt(0.5, 0.963,
    "two copies 99.9% identical  $\\cdot$  read 6 kb  $\\cdot$  error rate 1%  $\\cdot$  read truly from COPY A",
    12.5, ha='center', color=GREY, style='italic')

# =====================================================================
# PANEL 1 — THE SETUP
# =====================================================================
box(0.03, 0.815, 0.94, 0.135, '#f7f9fb')
txt(0.05, 0.940, "1.  The setup — copies are identical EXCEPT at a few \"diagnostic\" sites",
    14, weight='bold')
txt(0.05, 0.918,
    "Shaded columns are the only positions where the copies differ. A 35-bp window of the 6 kb read (2 of its 6 diagnostic sites shown):",
    11, color='#333')
draw_seq_row(0.892, COPY_A, "COPY A", 'A')
draw_seq_row(0.866, COPY_B, "COPY B", 'B')
draw_seq_row(0.838, READ,  "READ",  'read')
# markers under the window
for i in DIAG:
    txt(colx(i), 0.823, "diag", 8, ha='center', va='center', color=PURP, weight='bold')
for i in ERR:
    txt(colx(i), 0.823, "err", 8, ha='center', va='center', color=RED, weight='bold')
txt(0.05, 0.806, "", 1)

# =====================================================================
# PANEL 2 — WHY THE SCALAR TIES
# =====================================================================
box(0.03, 0.615, 0.94, 0.185, '#fdf2f0')
txt(0.05, 0.790, "2.  Why a whole-read scalar (AS / NM / de) ties",
    14, weight='bold', color=RED)
txt(0.05, 0.768,
    "Over 6 kb the read collects ~60 sequencing errors against EITHER copy. The 6 diagnostic differences are a rounding error on top:",
    11, color='#333')
# the tally
tx = 0.07
txt(tx, 0.744, f"mismatches vs COPY A  =  {n_err} errors  +  0 diagnostic        =  {nm_A}     -> identity {id_A:.2f}%",
    11.5, mono=True, color=BLUE)
txt(tx, 0.722, f"mismatches vs COPY B  =  {n_err} errors  +  {n_diag} diagnostic        =  {nm_B}     -> identity {id_B:.2f}%",
    11.5, mono=True, color=GREEN)
txt(tx, 0.700, f"difference            =  {nm_B-nm_A} bp  =  {100*(nm_B-nm_A)/L:.2f}% of the read   <<  the 1% error noise",
    11.5, mono=True, color='#333', weight='bold')
# verdict bar
box(0.07, 0.628, 0.86, 0.055, '#ffffff', ec=RED, lw=1.4)
txt(0.085, 0.675, "minimap2 verdict:", 11.5, weight='bold')
txt(0.085, 0.655,
    "AS(A) $\\approx$ AS(B)  $\\rightarrow$  reported as a TIE, MAPQ 0  (aligner picks a \"primary\" essentially at random)",
    11, color='#222')
txt(0.085, 0.637,
    "and the secondary alignment stores  SEQ = *  $\\rightarrow$  you cannot even recompute per-base divergence to break it",
    11, color=RED, weight='bold')

# =====================================================================
# PANEL 3 — THE REFRAME
# =====================================================================
box(0.03, 0.445, 0.94, 0.155, '#eef7f0')
txt(0.05, 0.590, "3.  The reframe — score ONLY the copy-distinguishing positions",
    14, weight='bold', color=GREEN)
txt(0.05, 0.568,
    "Errors land everywhere and are uninformative; throw them away. Keep the read's base at each diagnostic site = one VOTE for whichever copy's allele it matches:",
    11, color='#333')
# vote table for the 6 sites (2 shown explicitly, all 6 vote A)
votes_y = 0.540
hdr = ["site", "COPY A", "COPY B", "READ base", "vote"]
hx  = [0.10, 0.27, 0.41, 0.57, 0.74]
for x, h in zip(hx, hdr): txt(x, votes_y, h, 10.5, weight='bold', color='#333')
sample = [("exon2 +112", "A", "T", "A", "A"),
          ("exon3 +ornt", "G", "C", "G", "A"),
          ("... 4 more ...", "-", "-", "-", "A")]
for k, r in enumerate(sample):
    y = votes_y - 0.020 - k*0.020
    cols = [r[0], r[1], r[2], r[3], r[4]]
    colors = ['black', BLUE, GREEN,
              (BLUE if r[3]==r[1] else 'black'), BLUE]
    for x, v, c in zip(hx, cols, colors):
        txt(x, y, v, 10.5, mono=True, color=c,
            weight=('bold' if x in (hx[3], hx[4]) else 'normal'))
txt(0.10, votes_y - 0.020 - 3*0.020 - 0.004,
    f"all {n_diag} diagnostic sites the read spans agree -> COPY A   (the {nm_B-nm_A}-vote signal that was invisible to AS)",
    11, color=GREEN, weight='bold')

# =====================================================================
# PANEL 4 — THE GRAPH
# =====================================================================
box(0.03, 0.255, 0.94, 0.175, '#eef0f7')
txt(0.05, 0.420, "4.  Where the graph comes in", 14, weight='bold', color=PURP)
txt(0.05, 0.399,
    "Align the copies into ONE family graph: shared sequence = a single chain; each diagnostic site = a BUBBLE labelled with each copy's allele.",
    11, color='#333')
# draw a small graph: backbone nodes + 3 bubbles, read path on top alleles
gy = 0.330
def node(x, y, w, s, fc, tc='black'):
    box(x, y-0.012, w, 0.024, fc, ec='#888', lw=1.1)
    txt(x + w/2, y, s, 10, ha='center', va='center', color=tc, weight='bold', mono=True)
# backbone segments and bubbles alternating
xb = 0.07
seg_w, gap = 0.075, 0.045
positions = []
node(xb, gy, seg_w, "shared", '#e8e8e8'); positions.append(xb+seg_w)
xb += seg_w + gap
for bi in range(3):
    # bubble: top = A allele, bottom = B allele
    bx = xb
    node(bx, gy+0.030, 0.045, ["A","G","A"][bi], '#d6e6f2', BLUE)
    node(bx, gy-0.030, 0.045, ["T","C","G"][bi], '#d8efe0', GREEN)
    txt(bx+0.0225, gy+0.058, f"site {bi+1}", 8, ha='center', color=PURP)
    # connectors from previous backbone to both alleles and back
    arrow(positions[-1], gy, bx, gy+0.030, color='#aaa', lw=1.0)
    arrow(positions[-1], gy, bx, gy-0.030, color='#aaa', lw=1.0)
    xb += 0.045 + gap*0.6
    node(xb, gy, seg_w, "shared", '#e8e8e8'); positions.append(xb+seg_w)
    arrow(bx+0.045, gy+0.030, xb, gy, color='#aaa', lw=1.0)
    arrow(bx+0.045, gy-0.030, xb, gy, color='#aaa', lw=1.0)
    xb += seg_w + gap
# read path: thick translucent blue line threading the backbone + top (A) alleles
path_pts = [
    (0.070, gy), (0.145, gy),
    (0.190, gy+0.030), (0.235, gy+0.030),
    (0.262, gy), (0.337, gy),
    (0.382, gy+0.030), (0.427, gy+0.030),
    (0.454, gy), (0.529, gy),
    (0.574, gy+0.030), (0.619, gy+0.030),
    (0.646, gy), (0.721, gy),
]
ax.plot([p[0] for p in path_pts], [p[1] for p in path_pts],
        color=BLUE, lw=6, alpha=0.30, solid_capstyle='round', zorder=3.5)
txt(0.745, gy+0.030, "<- READ", 9.5, va='center', color=BLUE, weight='bold')
txt(0.07, gy-0.075,
    "READ's path threads the COPY-A allele at every bubble  =  it spells out copy A's haplotype",
    11, color=BLUE, weight='bold')
# three roles
roles = [
    "(a) the graph DEFINES the diagnostic sites (the bubbles) and each copy's allele - no guessing which mismatch matters",
    "(b) one read, one alignment, ALL copies at once - no arbitrary primary, no reference bias, NO  SEQ=*  loss on secondaries",
    "(c) the alleles are read JOINTLY: linkage across bubbles = haplotype phasing  ->  decisive even when no single site is",
]
for k, r in enumerate(roles):
    txt(0.06, 0.290 - k*0.0145, r, 10, color='#222')

# =====================================================================
# PANEL 5 — THE MATH
# =====================================================================
box(0.03, 0.095, 0.94, 0.145, '#fbf7ee')
txt(0.05, 0.230, "5.  The likelihood ratio — why a few clean votes crush a whole-read tie",
    14, weight='bold')
txt(0.06, 0.207,
    "log LR(A:B)  =  $\\sum_{sites}$ log [ P(base | A) / P(base | B) ]   "
    "$\\approx$   (A-votes  -  B-votes)  $\\times$  log[ (1-$\\varepsilon$) / ($\\varepsilon$/3) ]",
    13)
txt(0.06, 0.170,
    f"per decisive site:  log((1-{eps})/({eps}/3))  =  {w_site:.2f} nats        "
    f"6 sites, all A:  6 x {w_site:.2f}  =  {logLR:.1f} nats",
    11.5, mono=True, color='#222')
txt(0.06, 0.150,
    f"posterior  P(COPY A | read)  =  1 - {p_not_A:.1e}   (essentially certain)   "
    "<-> the scalar said 'tie'",
    11.5, mono=True, color=GREEN, weight='bold')
txt(0.06, 0.122,
    "Key: sequencing errors are SPREAD over thousands of positions (each negligible); the few diagnostic votes are CONCENTRATED and unanimous.",
    10.5, color='#333', style='italic')

# =====================================================================
# PANEL 6 — THE HONEST CAVEAT
# =====================================================================
box(0.03, 0.018, 0.94, 0.068, '#f2eef7', ec=PURP, lw=1.3)
txt(0.05, 0.077, "6.  The honest caveat — when it is genuinely non-identifiable",
    13, weight='bold', color=PURP)
txt(0.05, 0.055,
    "A read covering a stretch where the copies are IDENTICAL spans 0 bubbles -> 0 votes. No method (graph or BAM) can assign it from sequence.",
    10.5, color='#222')
txt(0.05, 0.036,
    "Correct response: do NOT force a guess. APPORTION it fractionally by the copies' expression measured from the reads that ARE resolvable (EM).",
    10.5, color=PURP, weight='bold')

fig.savefig("toy_copy_assignment.png", dpi=150, bbox_inches='tight')
print("wrote toy_copy_assignment.png")
print(f"  n_diag={n_diag} n_err={n_err} nm_A={nm_A} nm_B={nm_B} "
      f"id_A={id_A:.2f} id_B={id_B:.2f} w_site={w_site:.3f} logLR={logLR:.2f} p_not_A={p_not_A:.2e}")
