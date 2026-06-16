import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch

# ---- Palette (shared style guide) ----
NAVY = "#22313f"
GRAY = "#8898aa"
ORANGE = "#e8590c"
GREEN = "#1e8449"

SRC_FILL = "#dde3e8"
FLOW_FILL = "#eef1f4"
RC_FILL = "#fdece1"
GATE_FILL = "#e7f3ec"

fig, ax = plt.subplots(figsize=(12.2, 6.8))
ax.set_xlim(0, 120)
ax.set_ylim(0, 66)
ax.axis("off")


def box(cx, cy, w, h, label, fill, edge, tcolor=NAVY, lw=1.7, fs=10.5, weight="normal"):
    p = FancyBboxPatch((cx - w / 2.0, cy - h / 2.0), w, h,
                       boxstyle="round,pad=0.02,rounding_size=1.4",
                       linewidth=lw, edgecolor=edge, facecolor=fill, zorder=3)
    ax.add_patch(p)
    ax.text(cx, cy, label, ha="center", va="center", color=tcolor,
            fontsize=fs, fontweight=weight, zorder=4, linespacing=1.25)
    return dict(cx=cx, cy=cy, w=w, h=h, l=cx - w / 2, r=cx + w / 2,
                t=cy + h / 2, b=cy - h / 2)


def arrow(p0, p1, color, lw=2.2, rad=0.0, mut=17):
    a = FancyArrowPatch(p0, p1, arrowstyle="-|>",
                        connectionstyle=f"arc3,rad={rad}",
                        mutation_scale=mut, linewidth=lw, color=color,
                        zorder=2, shrinkA=1, shrinkB=1, capstyle="round")
    ax.add_patch(a)


def elbow(pts, color, lw=2.2):
    xs = [p[0] for p in pts]
    ys = [p[1] for p in pts]
    ax.plot(xs[:-1], ys[:-1], color=color, lw=lw, solid_capstyle="round",
            solid_joinstyle="round", zorder=2)
    arrow(pts[-2], pts[-1], color, lw=lw)


# ===== Lanes (y) =====
Y_FLOW = 50
Y_RC = 12
Y_OUT = 31

# ===== Columns (x) =====
X_SRC = 12
X_EXT = 33
X_FLR = 53
X_FILT = 53
X_TERM = 80    # flow+guide output / union back column
X_OUT = 106    # OUTPUT

# ---- Source ----
src = box(X_SRC, Y_OUT, 18, 8, "transfrags\n(bundle)", SRC_FILL, NAVY, fs=11, weight="bold")

# ---- FLOW lane (gray) ----
fext = box(X_EXT, Y_FLOW, 18, 8, "flow\nextractor", FLOW_FILL, GRAY, fs=10.5)
ffloor = box(X_FLR, Y_FLOW, 16, 7, "flow floor", FLOW_FILL, GRAY, fs=10)
filt = box(X_FILT, Y_OUT, 26, 9.5,
           "per-bundle + global filters\npairwise / isofrac / predcluster",
           FLOW_FILL, GRAY, fs=9.4)
fgo = box(X_TERM, Y_FLOW, 18, 8, "flow + guide\noutput", FLOW_FILL, GRAY, fs=10)

# ---- READ-CHAIN lane (orange) ----
rcext = box(X_EXT, Y_RC, 19, 8.5, "read-chain\nextractor", RC_FILL, ORANGE,
            tcolor=ORANGE, fs=10.5, weight="bold")
rcset = box(X_FLR, Y_RC, 16, 7, "read-chain set", RC_FILL, ORANGE,
            tcolor=ORANGE, fs=10)
gate = box(X_TERM, Y_RC, 25, 9,
           "realness gate\ncanonical + RT-switch + read-depth",
           GATE_FILL, GREEN, tcolor=GREEN, fs=9.0)
union = box(X_TERM, Y_RC + 13.5, 22, 7.5, "union back\nNOVEL chains only",
            RC_FILL, ORANGE, tcolor=ORANGE, fs=9.2)

# ---- OUTPUT ----
out = box(X_OUT, Y_OUT, 13, 11, "OUTPUT", NAVY, NAVY, tcolor="white",
          fs=13.5, weight="bold")

# ================= ARROWS =================
# Source splits into the two extractors
elbow([(src["r"], src["cy"] + 2), (X_EXT, src["cy"] + 2),
       (X_EXT, fext["b"])], GRAY)
elbow([(src["r"], src["cy"] - 2), (X_EXT, src["cy"] - 2),
       (X_EXT, rcext["t"])], ORANGE)

# FLOW lane: extractor -> flow floor -> (down) filters -> (right then up) flow+guide out -> OUTPUT
arrow((fext["r"], fext["cy"]), (ffloor["l"], ffloor["cy"]), GRAY)
elbow([(ffloor["cx"], ffloor["b"]), (filt["cx"], filt["t"])], GRAY)
elbow([(filt["r"], filt["cy"]), (X_TERM, filt["cy"]),
       (X_TERM, fgo["b"])], GRAY)
# flow + guide output -> OUTPUT (clean diagonal into top-left of OUTPUT)
arrow((fgo["r"], fgo["cy"]), (out["cx"], out["t"]), GRAY, rad=-0.16)

# READ-CHAIN lane: extractor -> set -> BYPASS -> gate -> union -> OUTPUT
arrow((rcext["r"], rcext["cy"]), (rcset["l"], rcset["cy"]), ORANGE)
arrow((rcset["r"], rcset["cy"]), (gate["l"], gate["cy"]), ORANGE)
# gate -> union (straight up)
arrow((gate["cx"], gate["t"]), (union["cx"], union["b"]), ORANGE)
# union back -> OUTPUT (clean diagonal into bottom-left of OUTPUT)
arrow((union["r"], union["cy"]), (out["cx"], out["b"]), ORANGE, rad=0.16)

# ================= ANNOTATIONS =================
# placed in clear space below OUTPUT
ax.text(out["cx"], Y_OUT - 9.5, r"output $\supseteq$ flow",
        ha="center", va="center", color=NAVY, fontsize=12.5, fontweight="bold")
ax.text(out["cx"], Y_OUT - 13.0, "never displaces\na flow find",
        ha="center", va="center", color=ORANGE, fontsize=8.6, style="italic",
        linespacing=1.2)

# small "bypass" tag on the read-chain route (between set and gate)
ax.text((rcset["r"] + gate["l"]) / 2.0, Y_RC + 3.4, "bypass",
        ha="center", va="center", color=ORANGE, fontsize=8.2, style="italic")

# Title
ax.text(60, 64.0, "Additive architecture: read-chain only adds",
        ha="center", va="center", color=NAVY, fontsize=14.5, fontweight="bold")

# Tiny lane legend (top-left clear area)
ax.plot([3.5, 8.3], [60.5, 60.5], color=GRAY, lw=3.2, solid_capstyle="round")
ax.text(9.3, 60.5, "flow path", ha="left", va="center", fontsize=9, color=NAVY)
ax.plot([3.5, 8.3], [57.3, 57.3], color=ORANGE, lw=3.2, solid_capstyle="round")
ax.text(9.3, 57.3, "read-chain path", ha="left", va="center", fontsize=9, color=ORANGE)

plt.tight_layout()
plt.savefig("/mnt/c/Users/jfris/Desktop/Rustle/bench/readcoherence_fig/fig3_architecture.png",
            dpi=190, bbox_inches="tight", facecolor="white")
print("done")
