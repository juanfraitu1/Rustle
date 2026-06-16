import matplotlib
matplotlib.use("Agg")
import json
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch

# ---- Palette (shared style guide) ----
NAVY = "#22313f"
GRAY = "#8898aa"
ORANGE = "#e8590c"
GREEN = "#1e8449"
RED = "#c0392b"

with open("/tmp/cre_guided/figdata/quant.json") as fh:
    Q = json.load(fh)

sweep = sorted(Q["sweep"], key=lambda d: d["min_cov"])
shipped_mc = Q["shipped_min_cov"]
st_tx = Q["stringtie_tx"]
flow_tx = Q["flow_tx"]
ung_rc = Q["tools"]["ungated_readchain"]["beat_st"]

fig, (axL, axR) = plt.subplots(
    1, 2, figsize=(13.6, 6.2),
    gridspec_kw={"width_ratios": [1.25, 1.0]})

# =====================================================================
# LEFT: recall vs transcript count (min_cov sweep)
# =====================================================================
xs = [d["total_tx"] / 1000.0 for d in sweep]
ys = [d["beat_st"] for d in sweep]

# reference lines
axL.axvline(st_tx / 1000.0, color=GRAY, ls="--", lw=1.4, zorder=1)
axL.axvline(flow_tx / 1000.0, color=GRAY, ls="--", lw=1.4, zorder=1)
axL.axhline(ung_rc, color=ORANGE, ls="--", lw=1.4, alpha=0.7, zorder=1)

# sweep curve
axL.plot(xs, ys, "-", color=ORANGE, lw=2.6, zorder=4, solid_capstyle="round")

for d in sweep:
    x = d["total_tx"] / 1000.0
    y = d["beat_st"]
    if d["min_cov"] == shipped_mc:
        axL.plot(x, y, marker="*", ms=22, color=ORANGE,
                 markeredgecolor=NAVY, markeredgewidth=1.2, zorder=6)
        axL.annotate("default\nmin_cov 3", (x, y),
                     textcoords="offset points", xytext=(14, 4),
                     ha="left", va="center", fontsize=9.5,
                     fontweight="bold", color=NAVY, linespacing=1.1)
    else:
        axL.plot(x, y, marker="o", ms=8.5, color="white",
                 markeredgecolor=ORANGE, markeredgewidth=2.2, zorder=5)
        axL.annotate(f"{d['min_cov']}", (x, y),
                     textcoords="offset points", xytext=(8, -13),
                     ha="left", va="top", fontsize=9, color=ORANGE)

# reference labels
ymax = max(ys)
axL.text(st_tx / 1000.0, ymax * 1.02, "StringTie", rotation=90,
         ha="right", va="top", fontsize=9, color=GRAY)
axL.text(flow_tx / 1000.0, ymax * 1.02, "flow", rotation=90,
         ha="left", va="top", fontsize=9, color=GRAY)
axL.text(xs[0], ung_rc, "ungated read-chain", ha="right", va="bottom",
         fontsize=8.8, color=ORANGE, style="italic")

# clarify what the two vertical reference lines mean + what labels are
axL.text(0.985, 0.10,
         "vertical lines = each tool's emitted count",
         transform=axL.transAxes, ha="right", va="bottom",
         fontsize=8.5, color=GRAY, style="italic")
axL.text(0.985, 0.05, "labels = min_cov", transform=axL.transAxes,
         ha="right", va="bottom", fontsize=8.5, color=GRAY, style="italic")

axL.set_xlabel("transcripts emitted  (thousands)", fontsize=11, color=NAVY)
axL.set_ylabel("StringTie-missed isoforms recovered", fontsize=11, color=NAVY)
axL.set_title("Recall vs transcript count", fontsize=13,
              fontweight="bold", color=NAVY, pad=10)
axL.set_xlim(65, 97)
axL.set_ylim(0, ymax * 1.18)
axL.spines["top"].set_visible(False)
axL.spines["right"].set_visible(False)
axL.spines["left"].set_color(GRAY)
axL.spines["bottom"].set_color(GRAY)
axL.tick_params(colors=NAVY, labelsize=9.5)

# =====================================================================
# RIGHT: what the read-chain extras are (SQANTI3 stacked bar)
# =====================================================================
ex = Q["sqanti_extras"]
n = sum(ex.values())

real = ex["full-splice_match"] + ex["novel_in_catalog"] + ex["novel_not_in_catalog"]
ism = ex["incomplete-splice_match"]
noise = (ex["intergenic"] + ex["antisense"] + ex["fusion"]
         + ex["genic_intron"] + ex["genic"])

# ordered segments: REAL (green shades), ISM (gray), NOISE (red shades)
segs = [
    ("FSM", ex["full-splice_match"], "#1e8449", "white"),
    ("NIC", ex["novel_in_catalog"], "#52be80", "white"),
    ("NNC", ex["novel_not_in_catalog"], "#a9dfbf", NAVY),
    ("ISM", ism, GRAY, "white"),
    ("", ex["intergenic"], "#c0392b", "white"),
    ("", ex["fusion"], "#d98880", NAVY),
    ("", ex["antisense"], "#e6b0aa", NAVY),
    ("", ex["genic_intron"] + ex["genic"], "#f2d7d5", NAVY),
]

y0 = 0.40
h = 0.30
left = 0.0
for label, val, col, tcol in segs:
    axR.barh(y0, val, left=left, height=h, color=col,
             edgecolor="white", linewidth=1.0, zorder=3)
    if label and val / n > 0.05:
        axR.text(left + val / 2.0, y0, label,
                 ha="center", va="center", fontsize=8.6,
                 color=tcol, fontweight="bold")
    left += val

# group brackets / labels above the bar (alternating heights to avoid overlap)
def group_label(x0, x1, text, color, lift):
    ytop = y0 + h / 2.0
    yb = ytop + 0.03
    axR.plot([x0, x1], [yb, yb], color=color, lw=2.0,
             solid_capstyle="butt", zorder=4)
    xc = (x0 + x1) / 2.0
    axR.plot([xc, xc], [yb, yb + lift - 0.02], color=color, lw=1.0, zorder=4)
    axR.text(xc, yb + lift, text, ha="center", va="bottom",
             fontsize=10, fontweight="bold", color=color)

group_label(0, real, f"REAL  {real:,}", GREEN, 0.035)
group_label(real, real + ism, f"ISM  {ism:,}", GRAY, 0.135)
group_label(real + ism, n, f"NOISE  {noise:,}", RED, 0.035)

# three short facts (one clean row below the bar)
pct = Q["sqanti_pct"]
facts = [
    (f"{pct['canonical']}% canonical", NAVY),
    (f"{pct['real_FSM_NIC_NNC']}% real", GREEN),
    (f"{pct['noise']}% noise", RED),
]
fx = [0.16, 0.5, 0.84]
for (txt, col), fpos in zip(facts, fx):
    axR.text(fpos * n, y0 - h / 2.0 - 0.13, txt, ha="center", va="center",
             fontsize=12, fontweight="bold", color=col,
             transform=axR.transData)

# ---- per-segment category legend (decodes every coloured sub-bar) ----
# two compact rows of swatch + label, so the lighter red / green shades and
# the inner FSM/NIC/NNC/ISM segments are unambiguous.
leg = [
    ("FSM", "#1e8449"), ("NIC", "#52be80"), ("NNC", "#a9dfbf"),
    ("ISM", GRAY),
    ("intergenic", "#c0392b"), ("fusion", "#d98880"),
    ("antisense", "#e6b0aa"), ("genic", "#f2d7d5"),
]
leg_row_y = [-0.085, -0.20]
sw = 0.022 * n          # swatch width (data units)
sh = 0.045              # swatch height (axes-y units)
col_x = [0.0, 0.265, 0.53, 0.80]   # 4 columns across the panel (fractions of n)
for k, (lbl, col) in enumerate(leg):
    row = k // 4
    cx = col_x[k % 4] * n
    ly = leg_row_y[row]
    axR.add_patch(plt.Rectangle((cx, ly - sh / 2), sw, sh,
                                facecolor=col, edgecolor="white",
                                linewidth=0.6, zorder=3, clip_on=False))
    axR.text(cx + sw * 1.3, ly, lbl, ha="left", va="center",
             fontsize=8.2, color=NAVY)

axR.set_xlim(-0.01 * n, 1.02 * n)
axR.set_ylim(-0.27, 0.95)
axR.set_yticks([])
axR.set_xticks([0, n / 2.0, n])
axR.set_xticklabels(["0", f"{int(n/2):,}", f"{n:,}"], fontsize=9, color=NAVY)
axR.set_xlabel("read-chain extras", fontsize=11, color=NAVY)
axR.set_title(f"What the extras are  (SQANTI3, n={n:,})",
              fontsize=13, fontweight="bold", color=NAVY, pad=10)
for sp in ("top", "right", "left"):
    axR.spines[sp].set_visible(False)
axR.spines["bottom"].set_color(GRAY)
axR.tick_params(colors=NAVY)

plt.tight_layout(w_pad=3.0)
plt.savefig("/mnt/c/Users/jfris/Desktop/Rustle/bench/readcoherence_fig/fig5_results.png",
            dpi=190, bbox_inches="tight", facecolor="white")
print("done", "real", real, "ism", ism, "noise", noise, "n", n)
