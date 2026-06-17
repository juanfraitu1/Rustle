#!/usr/bin/env python3
"""Figure for task (c): the direct-BAM hidden-collapse scan, after confound control, finds 0.

Left  — the decisive discriminator: 300/306 COLLAPSED_LIKE loci are UNIQUELY mapped (frac_mq0<0.10).
        A real second genomic copy MUST multimap, so 98% of the raw hits cannot be collapse.
Right — confound waterfall: 306 raw COLLAPSED_LIKE -> remove diploid het / segdup-spillover / RNA
        editing -> 0 confound-controlled hidden collapses (matching the annotation-based probe's 0).
Numbers from bench/hidden_collapse_aggregate.py + the deterministic joint gate.
"""
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch

NAVY = "#22313f"; SLATE = "#8898aa"; ORANGE = "#e8590c"
GREEN = "#1e8449"; RED = "#c0392b"
plt.rcParams["font.family"] = "DejaVu Sans"

fig = plt.figure(figsize=(13.4, 6.2))
axL = fig.add_axes([0.055, 0.12, 0.38, 0.74])
axR = fig.add_axes([0.55, 0.10, 0.42, 0.78])
for ax in (axL, axR):
    ax.set_xticks([]); ax.set_yticks([])
    for s in ax.spines.values():
        s.set_visible(False)

# ---------------- LEFT: frac_mq0 distribution ----------------
axL.set_xlim(0, 10); axL.set_ylim(0, 10)
axL.text(5, 9.6, "do the reads actually multimap?", fontsize=13, fontweight="bold",
         color=NAVY, ha="center")
axL.text(5, 8.9, "a real hidden copy MUST map to ≥2 places", fontsize=9.5,
         color=SLATE, ha="center", style="italic")

bars = [("uniquely mapped\n(frac_mq0 < 0.10)", 300, RED),
        ("weak\n(0.10–0.30)", 4, ORANGE),
        ("genuine local\nmultimap (≥0.30)", 2, GREEN)]
y = 7.2
for lab, val, col in bars:
    w = 8.2 * (val / 306) ** 0.42        # sqrt-ish so the 2 and 4 are visible
    w = max(w, 0.9)
    axL.add_patch(FancyBboxPatch((0.7, y - 0.5), w, 1.0,
                                 boxstyle="round,pad=0.02,rounding_size=0.06",
                                 facecolor=col, edgecolor="none", alpha=0.9))
    axL.text(0.85, y, f"{val}", fontsize=16, fontweight="bold", color="white", va="center")
    axL.text(0.7 + w + 0.25, y, lab, fontsize=9.2, color=NAVY, va="center")
    y -= 2.0
axL.text(5, 0.7, "300 / 306 (98%) cannot be a collapsed copy",
         fontsize=10.5, color=RED, ha="center", fontweight="bold")

# ---------------- RIGHT: confound waterfall ----------------
axR.set_xlim(0, 10); axR.set_ylim(0, 10)
axR.text(5, 9.6, "confound control: 306 → 0", fontsize=13, fontweight="bold",
         color=NAVY, ha="center")

steps = [
    ("306", "raw COLLAPSED_LIKE (895 isoforms, 59 FSM)", NAVY, 1.00),
    ("−167", "diploid heterozygosity (uniquely mapped)", SLATE, 0.74),
    ("−~120", "segdup spillover (no-SEQ secondaries)", SLATE, 0.48),
    ("−editing", "RNA editing (hits TIER A too)", SLATE, 0.30),
    ("0", "multimap ∧ ≥8 cols ∧ non-edit", RED, 0.13),
]
y = 8.2
maxw = 8.4
cx = 5.0
for i, (val, lab, col, frac) in enumerate(steps):
    w = maxw * (0.14 + 0.86 * frac)
    h = 1.0
    axR.add_patch(FancyBboxPatch((cx - w / 2, y - h / 2), w, h,
                                 boxstyle="round,pad=0.02,rounding_size=0.07",
                                 facecolor=col, edgecolor="none", alpha=0.92))
    axR.text(cx, y + 0.13, val, fontsize=15, fontweight="bold", color="white",
             ha="center", va="center")
    fs = 8.0 if frac > 0.25 else 7.3
    axR.text(cx, y - 0.30, lab, fontsize=fs, color="white", ha="center", va="center")
    if i < len(steps) - 1:
        axR.annotate("", xy=(cx, y - 0.5 - 0.62), xytext=(cx, y - 0.5 - 0.02),
                     arrowprops=dict(arrowstyle="-|>", color=SLATE, lw=1.8))
    y -= 1.62

axR.text(5, 0.25, "confound-controlled hidden-collapse headroom = 0   (= the annotation-based probe)",
         fontsize=9.8, color=RED, ha="center", fontweight="bold")

fig.suptitle("Task (c): scanning the BAM for HIDDEN collapsed copies the annotation misses — none survive",
             fontsize=14, fontweight="bold", color=NAVY, y=0.99)

out = "/mnt/c/Users/jfris/Desktop/Rustle/bench/hidden_collapse_headroom.png"
fig.savefig(out, dpi=185, bbox_inches="tight", facecolor="white")
print("saved", out)
