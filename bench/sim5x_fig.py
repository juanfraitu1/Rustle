#!/usr/bin/env python3
"""Figure: the 5-equally-good-places copy-assignment identifiability benchmark."""
import json
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

OUT = "/home/juanfra/winloci_scratch/sim5x"
d = json.load(open(f"{OUT}/sim5x_summary.json"))
S = d["summary"]; EC = d["err_curve"]
NAVY = "#22313f"; ORANGE = "#e8590c"; GREEN = "#1e8449"; RED = "#c0392b"; SLATE = "#8898aa"
plt.rcParams["font.family"] = "DejaVu Sans"

fig, (axA, axB) = plt.subplots(1, 2, figsize=(13.0, 5.2))

K = [s["K"] for s in S]
mq0 = [s["pct_MAPQ0"] if "pct_MAPQ0" in s else s["pct_mq0"] for s in S]
acc = [s["psv_acc"] * 100 for s in S]
ident = [s["identifiable"] / 5 * 100 for s in S]
x = range(len(K))
axA.plot(x, mq0, "-o", color=SLATE, lw=2, label="reads minimap2 leaves MAPQ=0\n(coordinates can't assign)")
axA.plot(x, acc, "-o", color=GREEN, lw=2.4, label="PSV-based assignment accuracy")
axA.plot(x, ident, "--s", color=ORANGE, lw=1.6, ms=5, label="copies identifiable (of 5)")
axA.axhline(20, color=RED, ls=":", lw=1, alpha=0.6)
axA.text(len(K)-1, 22, "random (1/5)", color=RED, fontsize=7.5, ha="right")
axA.set_xticks(list(x)); axA.set_xticklabels(K)
axA.set_xlabel("PSVs per copy  (K)"); axA.set_ylabel("%")
axA.set_title("5 identical copies (K=0) = unassignable;\n≥2 PSVs make all 5 separable", fontsize=11, color=NAVY)
axA.set_ylim(-5, 108); axA.legend(fontsize=7.6, loc="center right"); axA.grid(alpha=0.2)
axA.annotate("K=0: 5 EQUALLY good places\n(MAPQ 0, no PSV info)", (0, 100), fontsize=8, color=RED,
             xytext=(0.4, 70), arrowprops=dict(arrowstyle="->", color=RED, lw=1))

ee = [e for e, _ in EC]; aa = [a * 100 for _, a in EC]
axB.plot(range(len(ee)), aa, "-o", color=GREEN, lw=2.4)
axB.axvline(0, color=ORANGE, ls="--", lw=1.2); axB.text(0.1, 84, "HiFi (e=0.003)", color=ORANGE, fontsize=8)
axB.set_xticks(range(len(ee))); axB.set_xticklabels([f"{e:.0%}" for e in ee])
axB.set_xlabel("per-base error rate  (at K=4 PSVs)"); axB.set_ylabel("PSV assignment accuracy (%)")
axB.set_title("the error floor: PSVs must clear the\nper-base error to assign copies", fontsize=11, color=NAVY)
axB.set_ylim(75, 102); axB.grid(alpha=0.2)

fig.suptitle(f"Benchmark: 5 near-identical copies, '5 equally good places' (base gene {d['base_gene']})\n"
             "coordinates cannot assign reads to a copy — PSVs can, iff ≥K columns clear the error floor",
             fontsize=12.5, fontweight="bold", color=NAVY, y=1.02)
fig.tight_layout(rect=[0, 0, 1, 0.9])
out = "/mnt/c/Users/jfris/Desktop/Rustle/bench/sim5x_benchmark.png"
fig.savefig(out, dpi=160, facecolor="white", bbox_inches="tight")
print("saved", out)
