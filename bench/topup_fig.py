#!/usr/bin/env python3
"""Figure: Dataset 1 — ideal-coverage top-up. Real per-gene coverage distribution with the
under-covered tail (topped up to 40x) highlighted."""
import os
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

NAVY = "#22313f"; ORANGE = "#e8590c"; GREEN = "#1e8449"; SLATE = "#8898aa"
plt.rcParams["font.family"] = "DejaVu Sans"
nr = []
for line in open("/home/juanfra/winloci_scratch/gene_cov.tsv"):
    if line.startswith("gene\t"):
        continue
    nr.append(int(line.rstrip("\n").split("\t")[5]))
nr = np.array(nr)
TARGET = 40
under = int((nr < TARGET).sum()); zero = int((nr == 0).sum())

fig, ax = plt.subplots(figsize=(9.2, 5.0))
bins = np.linspace(0, 200, 41)
ax.hist(np.clip(nr, 0, 200), bins=bins, color=SLATE, alpha=0.85, edgecolor="white", linewidth=0.4)
ax.axvline(TARGET, color=ORANGE, lw=2, ls="--")
ax.axvspan(0, TARGET, color=ORANGE, alpha=0.10)
ax.text(TARGET + 3, ax.get_ylim()[1]*0.92, f"top-up target {TARGET}x", color=ORANGE, fontsize=10, fontweight="bold")
ax.text(2, ax.get_ylim()[1]*0.70,
        f"{under:,} genes under {TARGET}x\n({zero:,} with ZERO reads)\n→ simulated to {TARGET}x\n(96% achieved, 100% map)",
        color=NAVY, fontsize=9.5, va="top",
        bbox=dict(boxstyle="round,pad=0.4", facecolor="#fff4ec", edgecolor=ORANGE))
ax.set_xlabel("real IsoSeq reads per gene (clipped at 200)"); ax.set_ylabel("genes")
ax.set_title("Dataset 1: 'ideal coverage' GGO — top up only what the real data lacks\n"
             f"22,983 genes (median 43 reads); GGO_ideal.bam = real + 344,868 synthetic reads (1.8 GB)",
             fontsize=11.5, color=NAVY, fontweight="bold")
for s in ax.spines.values():
    s.set_color("#cccccc")
out = "/mnt/c/Users/jfris/Desktop/Rustle/bench/topup_coverage.png"
fig.tight_layout(); fig.savefig(out, dpi=160, facecolor="white")
print("saved", out)
