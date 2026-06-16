import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch
from matplotlib.lines import Line2D

# ---- palette ----
NAVY   = "#22313f"   # annotation
SLATE  = "#8898aa"   # StringTie / flow
ORANGE = "#e8590c"   # read-coherence (accent)
LIGHT  = "#cfd8dc"   # read alignments
AMBER  = "#f0932b"   # spanning alignments
GREEN  = "#1e8449"
RED    = "#c0392b"

# Toy locus: shared 5' (A) + alt middle (B/C/D) + shared 3' (E).
E = {
    "A": (0.4, 1.5),
    "B": (2.2, 2.9),
    "C": (3.2, 3.9),
    "D": (4.2, 4.9),
    "E": (5.6, 6.7),
}
EH = 0.30  # exon half-height

MOLS = [
    ["A", "B", "E"],
    ["A", "C", "E"],
    ["A", "D", "E"],
]


def draw_molecule(ax, exons, ycenter, exon_color, intron_color):
    xs = [E[e] for e in exons]
    for (x0a, x1a), (x0b, x1b) in zip(xs[:-1], xs[1:]):
        ax.plot([x1a, x0b], [ycenter, ycenter], color=intron_color,
                lw=1.3, zorder=2, solid_capstyle="round")
    for (x0, x1) in xs:
        r = FancyBboxPatch((x0, ycenter - EH), x1 - x0, 2 * EH,
                           boxstyle="round,pad=0.0,rounding_size=0.10",
                           linewidth=0, edgecolor="white",
                           facecolor=exon_color, zorder=3)
        ax.add_patch(r)


def style_axis(ax, title):
    ax.set_xlim(-0.3, 8.0)
    ax.set_ylim(-0.2, 9.7)
    ax.set_xticks([]); ax.set_yticks([])
    for s in ax.spines.values():
        s.set_visible(False)
    ax.set_title(title, fontsize=13, fontweight="bold", color=NAVY, pad=10)


fig, (axL, axR) = plt.subplots(1, 2, figsize=(11.6, 5.6))

# ============ y layout ============
y_mol   = [8.7, 7.8, 6.9]   # three input molecules (top)
y_graph = 4.4               # splice-graph band (middle)
y_out_L = [2.1, 1.1]        # left outputs (collapsed -> 2)
y_out_R = [2.6, 1.6, 0.6]   # right outputs (one per molecule)

LBL_X = -0.25

# ---------------- LEFT: flow assembly (StringTie) ----------------
style_axis(axL, "flow assembly (StringTie)")

axL.text(LBL_X, y_mol[0] + 0.55, "molecules", fontsize=8.5, color=NAVY,
         style="italic", ha="left")
for mol, yc in zip(MOLS, y_mol):
    draw_molecule(axL, mol, yc, LIGHT, "#b0bec5")

# converging funnel arrows: molecules -> graph
funnel_apex = (3.55, y_graph + 1.25)
for yc in y_mol:
    a = FancyArrowPatch((3.55, yc - 0.45), funnel_apex,
                        arrowstyle="-", color="#b8c2cc", lw=1.1, zorder=1)
    axL.add_patch(a)
a = FancyArrowPatch(funnel_apex, (3.55, y_graph + 0.65),
                    arrowstyle="-|>", mutation_scale=12,
                    color="#b8c2cc", lw=1.1, zorder=1)
axL.add_patch(a)

# shared splice graph
axL.text(LBL_X, y_graph + 1.0, "splice graph", fontsize=8.5, color=NAVY,
         style="italic", ha="left")
node_pos = {
    "A": (1.0, y_graph),
    "B": (3.1, y_graph + 0.8),
    "C": (3.6, y_graph),
    "D": (4.1, y_graph - 0.8),
    "E": (6.1, y_graph),
}
for mid in ["B", "C", "D"]:
    for (na, nb) in [("A", mid), (mid, "E")]:
        x0, y0 = node_pos[na]; x1, y1 = node_pos[nb]
        axL.plot([x0, x1], [y0, y1], color=SLATE, lw=1.4, zorder=2,
                 solid_capstyle="round")
for n, (x, y) in node_pos.items():
    axL.add_patch(plt.Circle((x, y), 0.26, facecolor=SLATE,
                             edgecolor="white", lw=1.2, zorder=3))
    axL.text(x, y, n, fontsize=8.5, color="white", ha="center",
             va="center", fontweight="bold", zorder=4)

# arrow graph -> output (with parsimony label beside it)
a = FancyArrowPatch((3.55, y_graph - 1.15), (3.55, y_out_L[0] + 0.55),
                    arrowstyle="-|>", mutation_scale=14, color=SLATE,
                    lw=1.9, zorder=2)
axL.add_patch(a)
axL.text(4.0, (y_graph - 1.15 + y_out_L[0]) / 2, "parsimony:\nmin #transcripts",
         fontsize=8.3, color=SLATE, ha="left", va="center", fontweight="bold")

# outputs: collapsed -> 2 transcripts (B-form and D-form); C lost
axL.text(LBL_X, y_out_L[0] + 0.55, "output", fontsize=8.5, color=NAVY,
         style="italic", ha="left")
draw_molecule(axL, ["A", "B", "E"], y_out_L[0], SLATE, SLATE)
draw_molecule(axL, ["A", "D", "E"], y_out_L[1], SLATE, SLATE)
axL.text(6.95, (y_out_L[0] + y_out_L[1]) / 2, "2 tx", fontsize=11.5,
         color=SLATE, ha="right", va="center", fontweight="bold")
axL.text(3.4, y_out_L[1] - 0.85, "C merged / dropped", fontsize=8.3,
         color=SLATE, ha="center", va="center")

# ---------------- RIGHT: read-coherent extraction ----------------
style_axis(axR, "read-coherent extraction")

axR.text(LBL_X, y_mol[0] + 0.55, "molecules", fontsize=8.5, color=NAVY,
         style="italic", ha="left")
for mol, yc in zip(MOLS, y_mol):
    draw_molecule(axR, mol, yc, AMBER, "#d98a2b")

# one short straight-down arrow per molecule -> its own output
# (mirrors the left funnel: one molecule -> one transcript, read instantly)
RAIL = 7.15
for yc_in, yc_out in zip(y_mol, y_out_R):
    a = FancyArrowPatch((RAIL, yc_in - 0.45), (RAIL, yc_out + 0.45),
                        arrowstyle="-|>", mutation_scale=12,
                        color=ORANGE, lw=1.6, zorder=1)
    axR.add_patch(a)

axR.text(3.5, (y_mol[-1] + y_out_R[0]) / 2, "emit each path",
         fontsize=9.0, color=ORANGE, ha="center", va="center",
         fontweight="bold")

# outputs: one per molecule, preserved
axR.text(LBL_X, y_out_R[0] + 0.55, "output", fontsize=8.5, color=NAVY,
         style="italic", ha="left")
for mol, yc in zip(MOLS, y_out_R):
    draw_molecule(axR, mol, yc, ORANGE, ORANGE)
# "3 tx": match the weight/placement of left's "2 tx" (same size, to the
# right of the output stack) so the 2-vs-3 count contrast is symmetric
axR.text(7.95, y_out_R[1], "3 tx", fontsize=11.5, color=ORANGE,
         ha="right", va="center", fontweight="bold")

# ---- shared legend (bottom, minimal) ----
handles = [
    Line2D([0], [0], marker="s", color="none", markerfacecolor=LIGHT,
           markersize=11, label="molecule (alignments)"),
    Line2D([0], [0], marker="s", color="none", markerfacecolor=SLATE,
           markersize=11, label="flow / StringTie"),
    Line2D([0], [0], marker="s", color="none", markerfacecolor=ORANGE,
           markersize=11, label="read-coherence"),
]
fig.legend(handles=handles, loc="lower center", ncol=3, frameon=False,
           fontsize=9, bbox_to_anchor=(0.5, 0.0), handletextpad=0.4,
           columnspacing=1.8)

fig.suptitle("Why the two regimes differ", fontsize=14.5,
             fontweight="bold", color=NAVY, y=1.01)

fig.tight_layout(rect=[0, 0.05, 1, 0.95])
out = "/mnt/c/Users/jfris/Desktop/Rustle/bench/readcoherence_fig/fig2_mechanism.png"
fig.savefig(out, dpi=190, bbox_inches="tight", facecolor="white")
print("saved", out)
