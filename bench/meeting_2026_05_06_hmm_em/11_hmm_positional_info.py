#!/usr/bin/env python3
"""Slide 11 — How the profile HMM preserves positional information across variant types.

Each variant — SNP, indel, alt-donor, alt-acceptor, UTR, exon-skipping — maps
to a specific feature of the M/I/D trellis or graph topology. The position
anchor (state index) is never lost: column j in the exon-class profile is
column j whether or not the read has a base there."""
import matplotlib.pyplot as plt
import matplotlib.patches as patches
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch

fig, ax = plt.subplots(figsize=(16, 13))
ax.set_xlim(0, 16)
ax.set_ylim(0, 13)
ax.axis('off')

# ── Title ───────────────────────────────────────────────────────────────
ax.text(8, 12.55,
        "Slide 11  —  How profile HMMs preserve positional information across variant types",
        fontsize=16.5, fontweight='bold', ha='center')
ax.text(8, 12.10,
        "Every variant maps to a specific feature of the M/I/D trellis or graph topology.  "
        "The state index acts as a ruler — column j is column j of the exon-class profile, regardless of read length.",
        fontsize=10.5, ha='center', color='#444', fontstyle='italic')

# ── Reference panel ─────────────────────────────────────────────────────
ref_box = FancyBboxPatch((0.45, 9.75), 15.10, 1.85,
                         boxstyle="round,pad=0.08",
                         facecolor='#f0f7fb', edgecolor='#1f77b4', lw=1.2)
ax.add_patch(ref_box)
ax.text(0.65, 11.45, "Reference profile  —  one exon-class node, 10 match columns",
        fontsize=11.5, fontweight='bold', color='#1f77b4', va='top')

ref_y = 10.65
col_x_ref = lambda j: 2.5 + j * 1.20
for j in range(10):
    cx = col_x_ref(j)
    c = patches.Circle((cx, ref_y), 0.28, facecolor='#cce4f7', edgecolor='#1f77b4', lw=1.3)
    ax.add_patch(c)
    ax.text(cx, ref_y, f"M{j+1}", ha='center', va='center', fontsize=7.5, fontweight='bold')
    ax.text(cx, ref_y - 0.55, "ATGCAGTACG"[j], ha='center', va='center',
            fontsize=11, family='monospace')
for j in range(9):
    a = FancyArrowPatch((col_x_ref(j) + 0.28, ref_y), (col_x_ref(j+1) - 0.28, ref_y),
                        arrowstyle='->', color='#888', lw=0.8)
    ax.add_patch(a)
ax.text(1.10, ref_y, "5'", fontsize=10, fontweight='bold', color='#666', va='center')
ax.text(1.10, ref_y - 0.55, "bases:", fontsize=8.5, color='#666', va='center')

ax.text(8, 9.95,
        "Each match-state owns a P(base | column) emission distribution.  "
        "The six panels below show what each variant type does to this ruler — and which feature absorbs the change.",
        fontsize=9.5, ha='center', va='center', color='#1a1a1a', fontstyle='italic')

# ── Panel layout ────────────────────────────────────────────────────────
PANEL_W, PANEL_H = 5.00, 4.30
ROW1_Y, ROW2_Y = 5.00, 0.50
COL_X = [0.45, 5.55, 10.65]


def frame(x, y, w, h, idx, title, color):
    box = FancyBboxPatch((x, y), w, h, boxstyle="round,pad=0.08",
                         facecolor='#ffffff', edgecolor=color, lw=1.4)
    ax.add_patch(box)
    ax.text(x + 0.18, y + h - 0.30, f"({idx})  {title}",
            fontsize=11.5, fontweight='bold', color=color, va='top')


def takeaway(x, y, w, text):
    ax.text(x + w/2, y + 0.20, "▸ " + text,
            fontsize=8.8, ha='center', va='center', color='#1a1a1a', fontweight='bold')


def m_circle(cx, cy, lbl, r=0.16, fc='#cce4f7', ec='#1f77b4'):
    c = patches.Circle((cx, cy), r, facecolor=fc, edgecolor=ec, lw=1.0)
    ax.add_patch(c)
    ax.text(cx, cy, lbl, ha='center', va='center', fontsize=5.8, fontweight='bold')


def i_square(cx, cy, lbl='I', r=0.16, fc='#fff3cd', ec='#d4a04a'):
    rect = patches.Rectangle((cx - r, cy - r), 2*r, 2*r, facecolor=fc, edgecolor=ec, lw=1.0)
    ax.add_patch(rect)
    ax.text(cx, cy, lbl, ha='center', va='center', fontsize=6, fontweight='bold', color='#7a5500')


def d_diamond(cx, cy, lbl, r=0.18, fc='#dddddd', ec='#666'):
    rho = patches.RegularPolygon((cx, cy), 4, radius=r, orientation=0.785,
                                 facecolor=fc, edgecolor=ec, lw=1.0)
    ax.add_patch(rho)
    ax.text(cx, cy, lbl, ha='center', va='center', fontsize=5.5, fontweight='bold', color='#444')


# ─── Panel 1: SNP ───────────────────────────────────────────────────────
x, y = COL_X[0], ROW1_Y
frame(x, y, PANEL_W, PANEL_H, 1, "SNP at one column", '#a02060')

mx = x + 1.05
my = y + PANEL_H - 0.85
# State row
ax.text(mx - 0.55, my, "state:", fontsize=8, ha='right', va='center', color='#666')
for j in range(8):
    cx = mx + j * 0.45
    fc = '#fcd0d6' if j == 4 else '#cce4f7'
    ec = '#a02060' if j == 4 else '#1f77b4'
    m_circle(cx, my, str(j+1), fc=fc, ec=ec)
# Par A bases
ax.text(mx - 0.55, my - 0.50, "par A:", fontsize=8, ha='right', va='center', color='#666')
for j, b in enumerate("ATGCAGTA"):
    cx = mx + j * 0.45
    if j == 4:
        rect = patches.Rectangle((cx - 0.18, my - 0.50 - 0.16), 0.36, 0.32,
                                 facecolor='#fdebef', edgecolor='#a02060', lw=0.8)
        ax.add_patch(rect)
    ax.text(cx, my - 0.50, b, ha='center', va='center',
            fontsize=9.5, family='monospace',
            color='#a02060' if j == 4 else '#222',
            fontweight='bold' if j == 4 else 'normal')
# Par B bases
ax.text(mx - 0.55, my - 0.95, "par B:", fontsize=8, ha='right', va='center', color='#666')
for j, b in enumerate("ATGCGGTA"):
    cx = mx + j * 0.45
    if j == 4:
        rect = patches.Rectangle((cx - 0.18, my - 0.95 - 0.16), 0.36, 0.32,
                                 facecolor='#fdebef', edgecolor='#a02060', lw=0.8)
        ax.add_patch(rect)
    ax.text(cx, my - 0.95, b, ha='center', va='center',
            fontsize=9.5, family='monospace',
            color='#a02060' if j == 4 else '#222',
            fontweight='bold' if j == 4 else 'normal')

exp_y = y + 1.75
ax.text(x + 0.30, exp_y, "HMM mapping:",
        fontsize=9.8, fontweight='bold', va='top', color='#a02060')
exp_y -= 0.35
for line in [
    "Same M_5 state for both paralogs.",
    "Per-copy emission differs (slide 4):",
    "    e_M(b|5)  for A:  P(A) = 0.95",
    "    e_M(b|5)  for B:  P(G) = 0.95",
]:
    ax.text(x + 0.30, exp_y, line, fontsize=8.7, va='top',
            family='monospace' if line.startswith("    ") else 'sans-serif',
            color='#222')
    exp_y -= 0.30

takeaway(x, y, PANEL_W,
         "Column index preserved; SNP signal lives in per-copy emission.")

# ─── Panel 2: Indel (insertion + deletion stacked) ──────────────────────
x, y = COL_X[1], ROW1_Y
frame(x, y, PANEL_W, PANEL_H, 2, "Indel  (insertion / deletion)", '#7a5500')

# Insertion (top half)
mx = x + 1.05
my = y + PANEL_H - 0.85
ax.text(x + 0.30, my + 0.32, "Insertion — par A: +3 nt between cols 4 and 5",
        fontsize=8.7, fontweight='bold', color='#7a5500')

# state row: M1..M4, I, I, I, M5..M7
ins_states = [('M', '1'), ('M', '2'), ('M', '3'), ('M', '4'),
              ('I', 'I'), ('I', 'I'), ('I', 'I'),
              ('M', '5'), ('M', '6'), ('M', '7')]
for j, (kind, lbl) in enumerate(ins_states):
    cx = mx + j * 0.40
    if kind == 'I':
        i_square(cx, my, lbl, r=0.14)
    else:
        m_circle(cx, my, lbl, r=0.14)
# par A bases
ax.text(mx - 0.55, my - 0.40, "par A:", fontsize=8, ha='right', va='center', color='#666')
ins_bases = "ATGCgatAGT"  # lowercase = inserted bases
for j, b in enumerate(ins_bases):
    cx = mx + j * 0.40
    color = '#7a5500' if b.islower() else '#222'
    fw = 'bold' if b.islower() else 'normal'
    ax.text(cx, my - 0.40, b.upper(), ha='center', va='center',
            fontsize=8.5, family='monospace', color=color, fontweight=fw)

# Deletion (lower half)
my2 = y + PANEL_H - 2.15
ax.text(x + 0.30, my2 + 0.32, "Deletion — par B: cols 5–7 missing (no bases consumed)",
        fontsize=8.7, fontweight='bold', color='#7a5500')

del_states = [('M', '1'), ('M', '2'), ('M', '3'), ('M', '4'),
              ('D', '5'), ('D', '6'), ('D', '7'),
              ('M', '8'), ('M', '9'), ('M', '10')]
for j, (kind, lbl) in enumerate(del_states):
    cx = mx + j * 0.40
    if kind == 'D':
        d_diamond(cx, my2, lbl, r=0.16)
    else:
        m_circle(cx, my2, lbl, r=0.14)
ax.text(mx - 0.55, my2 - 0.40, "par B:", fontsize=8, ha='right', va='center', color='#666')
del_bases = "ATGC---ACG"
for j, b in enumerate(del_bases):
    cx = mx + j * 0.40
    if b == '-':
        ax.text(cx, my2 - 0.40, '–', ha='center', va='center',
                fontsize=10, family='monospace', color='#aaa')
    else:
        ax.text(cx, my2 - 0.40, b, ha='center', va='center',
                fontsize=8.5, family='monospace', color='#222')

exp_y = y + 1.10
ax.text(x + 0.30, exp_y, "HMM mapping:",
        fontsize=9.8, fontweight='bold', va='top', color='#7a5500')
exp_y -= 0.32
for line in [
    "Insertion → I-states emit between M-cols",
    "Deletion  → D-states (silent, 0 bases)",
]:
    ax.text(x + 0.30, exp_y, line, fontsize=8.7, va='top', color='#222')
    exp_y -= 0.28

takeaway(x, y, PANEL_W,
         "Column anchors stay put; indels live in I (between) and D (skip) states.")

# ─── Panel 3: Alt donor (5' splice) ─────────────────────────────────────
x, y = COL_X[2], ROW1_Y
frame(x, y, PANEL_W, PANEL_H, 3, "Alternate 5' splice donor", '#207020')

mx = x + 0.40
# par A: 6 columns then GT donor
parA_y = y + PANEL_H - 0.85
ax.text(x + 0.30, parA_y + 0.32, "par A — exon ends col 6:",
        fontsize=8.7, fontweight='bold', color='#207020')
for j in range(6):
    cx = mx + 0.45 + j * 0.34
    m_circle(cx, parA_y, str(j+1), r=0.12)
# arrow + GT
arr_x_end = mx + 0.45 + 6*0.34 + 0.05
ax.annotate("", xy=(arr_x_end + 0.18, parA_y), xytext=(arr_x_end, parA_y),
            arrowprops=dict(arrowstyle='->', color='#207020', lw=1.4))
ax.text(arr_x_end + 0.40, parA_y, "GT", fontsize=9, color='#207020',
        fontweight='bold', va='center')
ax.text(arr_x_end + 0.78, parA_y - 0.25, "→ next exon",
        fontsize=7.5, color='#207020', va='center', fontstyle='italic')

# par B: 9 columns (donor +3) then GT
parB_y = parA_y - 1.20
ax.text(x + 0.30, parB_y + 0.32, "par B — exon ends col 9 (donor +3 nt):",
        fontsize=8.7, fontweight='bold', color='#207020')
for j in range(9):
    cx = mx + 0.45 + j * 0.34
    fc = '#cce4f7' if j < 6 else '#d4edda'
    m_circle(cx, parB_y, str(j+1), r=0.12, fc=fc, ec='#207020' if j >= 6 else '#1f77b4')
arr_x_end_b = mx + 0.45 + 9*0.34 + 0.05
ax.annotate("", xy=(arr_x_end_b + 0.18, parB_y), xytext=(arr_x_end_b, parB_y),
            arrowprops=dict(arrowstyle='->', color='#207020', lw=1.4))
ax.text(arr_x_end_b + 0.40, parB_y, "GT", fontsize=9, color='#207020',
        fontweight='bold', va='center')

exp_y = y + 1.55
ax.text(x + 0.30, exp_y, "HMM mapping:",
        fontsize=9.8, fontweight='bold', va='top', color='#207020')
exp_y -= 0.32
for line in [
    "Each paralog's exon-class node has",
    "its own LENGTH (per-copy profile).",
    "Boundary threading exits at",
    "col 6 vs col 9 — different junction.",
]:
    ax.text(x + 0.30, exp_y, line, fontsize=8.6, va='top', color='#222')
    exp_y -= 0.27

takeaway(x, y, PANEL_W,
         "Per-copy node length encodes alt-donor positions.")

# ─── Panel 4: Alt acceptor (3' splice) ──────────────────────────────────
x, y = COL_X[0], ROW2_Y
frame(x, y, PANEL_W, PANEL_H, 4, "Alternate 3' splice acceptor", '#207020')

mx = x + 0.45
# par A: 8 columns starting with AG acceptor
parA_y = y + PANEL_H - 0.85
ax.text(x + 0.30, parA_y + 0.32, "par A — exon starts col 1:",
        fontsize=8.7, fontweight='bold', color='#207020')
ax.text(mx + 0.10, parA_y, "AG", fontsize=9, color='#207020',
        fontweight='bold', va='center')
ax.annotate("", xy=(mx + 0.55 - 0.05, parA_y), xytext=(mx + 0.40, parA_y),
            arrowprops=dict(arrowstyle='->', color='#207020', lw=1.4))
for j in range(8):
    cx = mx + 0.70 + j * 0.42
    m_circle(cx, parA_y, str(j+1), r=0.13)
ax.text(mx + 0.10, parA_y - 0.25, "prev exon →", fontsize=7.5,
        color='#207020', va='center', fontstyle='italic')

# par B: shifted +3 — first 3 columns are paralog-specific, share cols 4..8
parB_y = parA_y - 1.20
ax.text(x + 0.30, parB_y + 0.32, "par B — exon starts col 4 (acceptor +3 nt):",
        fontsize=8.7, fontweight='bold', color='#207020')
ax.text(mx + 0.10, parB_y, "AG", fontsize=9, color='#207020',
        fontweight='bold', va='center')
ax.annotate("", xy=(mx + 0.55 - 0.05, parB_y), xytext=(mx + 0.40, parB_y),
            arrowprops=dict(arrowstyle='->', color='#207020', lw=1.4))
# show par B starting at col 4 (skip first 3)
for j in range(8):
    cx = mx + 0.70 + j * 0.42
    if j < 3:
        # ghost columns — par B doesn't use cols 1..3
        c = patches.Circle((cx, parB_y), 0.13, facecolor='#f5f5f5',
                           edgecolor='#bbb', lw=0.8, linestyle='--')
        ax.add_patch(c)
        ax.text(cx, parB_y, str(j+1), ha='center', va='center',
                fontsize=5.8, color='#999')
    else:
        m_circle(cx, parB_y, str(j+1), r=0.13)

exp_y = y + 1.55
ax.text(x + 0.30, exp_y, "HMM mapping:",
        fontsize=9.8, fontweight='bold', va='top', color='#207020')
exp_y -= 0.32
for line in [
    "Per-copy profile starts at a different",
    "effective column.  Equivalent: graph",
    "carries paralog-specific entry edges",
    "into the same exon-class node.",
]:
    ax.text(x + 0.30, exp_y, line, fontsize=8.6, va='top', color='#222')
    exp_y -= 0.27

takeaway(x, y, PANEL_W,
         "Paralog-specific entry column / entry edge per copy.")

# ─── Panel 5: UTR variation ─────────────────────────────────────────────
x, y = COL_X[1], ROW2_Y
frame(x, y, PANEL_W, PANEL_H, 5, "UTR length / sequence variation", '#6a3d9a')

# Show two paralogs with different 5'UTR lengths feeding into shared CDS
mx = x + 0.30
y_top = y + PANEL_H - 0.85

# par A: short UTR → CDS
ax.text(x + 0.30, y_top + 0.32, "par A — 5'UTR short (4 cols), then CDS exon:",
        fontsize=8.6, fontweight='bold', color='#6a3d9a')
# UTR-A node (purple)
utr_a_x0 = mx + 0.30
for j in range(4):
    cx = utr_a_x0 + j * 0.36
    m_circle(cx, y_top, str(j+1), r=0.13, fc='#e6d6ed', ec='#6a3d9a')
# arrow into CDS
cds_x0 = utr_a_x0 + 4*0.36 + 0.20
ax.annotate("", xy=(cds_x0 - 0.04, y_top), xytext=(utr_a_x0 + 4*0.36 + 0.04, y_top),
            arrowprops=dict(arrowstyle='->', color='#444', lw=1.0))
for j in range(6):
    cx = cds_x0 + 0.10 + j * 0.36
    m_circle(cx, y_top, str(j+1), r=0.13)
ax.text(cds_x0 + 0.10 + 6*0.36 + 0.10, y_top, "CDS", fontsize=8,
        color='#1f77b4', fontweight='bold', va='center')

# par B: long UTR → CDS
y_bot = y_top - 1.20
ax.text(x + 0.30, y_bot + 0.32, "par B — 5'UTR long (8 cols), then CDS exon:",
        fontsize=8.6, fontweight='bold', color='#6a3d9a')
utr_b_x0 = mx + 0.30
for j in range(8):
    cx = utr_b_x0 + j * 0.30
    m_circle(cx, y_bot, str(j+1), r=0.11, fc='#e6d6ed', ec='#6a3d9a')
cds_x0_b = utr_b_x0 + 8*0.30 + 0.20
ax.annotate("", xy=(cds_x0_b - 0.04, y_bot), xytext=(utr_b_x0 + 8*0.30 + 0.04, y_bot),
            arrowprops=dict(arrowstyle='->', color='#444', lw=1.0))
for j in range(6):
    cx = cds_x0_b + 0.10 + j * 0.30
    m_circle(cx, y_bot, str(j+1), r=0.11)

exp_y = y + 1.50
ax.text(x + 0.30, exp_y, "HMM mapping:",
        fontsize=9.8, fontweight='bold', va='top', color='#6a3d9a')
exp_y -= 0.32
for line in [
    "UTR is its own exon-class node.",
    "Each paralog has its own per-copy",
    "profile of paralog-specific length.",
    "Read need not span UTR (partial OK).",
]:
    ax.text(x + 0.30, exp_y, line, fontsize=8.6, va='top', color='#222')
    exp_y -= 0.27

takeaway(x, y, PANEL_W,
         "UTR = separate node, paralog-specific length, partial coverage allowed.")

# ─── Panel 6: Cassette / exon skipping ──────────────────────────────────
x, y = COL_X[2], ROW2_Y
frame(x, y, PANEL_W, PANEL_H, 6, "Cassette exon  /  skipping", '#a05a00')

# Show graph topology with E1, E2, E3 nodes and bypass edge
node_y_top = y + PANEL_H - 1.15
node_y_bot = node_y_top - 1.00

# E1 node
e1_x = x + 0.55
e1 = FancyBboxPatch((e1_x - 0.40, node_y_top - 0.30), 0.80, 0.60,
                    boxstyle="round,pad=0.02", facecolor='#cce4f7', edgecolor='#1f77b4', lw=1.2)
ax.add_patch(e1)
ax.text(e1_x, node_y_top, "E1", ha='center', va='center', fontsize=10, fontweight='bold')

# E2 node (cassette — top path; orange-tinted to mark it as the skip-able one)
e2_x = e1_x + 1.55
e2 = FancyBboxPatch((e2_x - 0.40, node_y_top - 0.30), 0.80, 0.60,
                    boxstyle="round,pad=0.02", facecolor='#fce4b3', edgecolor='#a05a00', lw=1.4)
ax.add_patch(e2)
ax.text(e2_x, node_y_top, "E2", ha='center', va='center', fontsize=10, fontweight='bold')

# E3 node
e3_x = e2_x + 1.55
e3 = FancyBboxPatch((e3_x - 0.40, node_y_top - 0.30), 0.80, 0.60,
                    boxstyle="round,pad=0.02", facecolor='#cce4f7', edgecolor='#1f77b4', lw=1.2)
ax.add_patch(e3)
ax.text(e3_x, node_y_top, "E3", ha='center', va='center', fontsize=10, fontweight='bold')

# Top path: E1 → E2 → E3 (par A includes the cassette)
ax.annotate("", xy=(e2_x - 0.40, node_y_top), xytext=(e1_x + 0.40, node_y_top),
            arrowprops=dict(arrowstyle='->', color='#1f77b4', lw=1.6))
ax.annotate("", xy=(e3_x - 0.40, node_y_top), xytext=(e2_x + 0.40, node_y_top),
            arrowprops=dict(arrowstyle='->', color='#1f77b4', lw=1.6))
ax.text((e1_x + e3_x) / 2, node_y_top + 0.45, "par A path",
        fontsize=8, color='#1f77b4', fontweight='bold', ha='center')

# Bottom path: E1 → E3 (par B skips cassette)
arc_arrow = FancyArrowPatch((e1_x + 0.40, node_y_top - 0.05),
                            (e3_x - 0.40, node_y_top - 0.05),
                            connectionstyle="arc3,rad=-0.45",
                            arrowstyle='->', color='#a05a00', lw=1.6, linestyle='--')
ax.add_patch(arc_arrow)
ax.text((e1_x + e3_x) / 2, node_y_bot - 0.10, "par B path  (skips E2)",
        fontsize=8, color='#a05a00', fontweight='bold', ha='center', fontstyle='italic')

exp_y = y + 1.45
ax.text(x + 0.30, exp_y, "HMM mapping:",
        fontsize=9.8, fontweight='bold', va='top', color='#a05a00')
exp_y -= 0.32
for line in [
    "Same M/I/D trellis inside every node.",
    "Different paralogs traverse different",
    "EDGES through the family graph.",
    "forward_against_path_for_copy chains",
    "the trellises along each copy's path.",
]:
    ax.text(x + 0.30, exp_y, line, fontsize=8.4, va='top', color='#222')
    exp_y -= 0.25

takeaway(x, y, PANEL_W,
         "Skipping = graph topology; trellises themselves don't change.")

# ── Bottom unifying takeaway ───────────────────────────────────────────
# (none — leave room for natural figure margins; the per-panel takeaways carry the message)

plt.tight_layout()
plt.savefig('/scratch/jxi21/Assembler/Rustle/bench/meeting_2026_05_06_hmm_em/11_hmm_positional_info.png',
            dpi=160, bbox_inches='tight')
print("11_hmm_positional_info.png written")
