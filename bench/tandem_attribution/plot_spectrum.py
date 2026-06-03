#!/usr/bin/env python3
"""Plot the per-copy attribution identifiability spectrum from the sweep.

Data are the mean-over-seeds results from the identity gradient (reads/copy=150,
3 seeds). x-axis: copy-distinguishing EXONIC SNPs per copy (the quantity that
actually reaches a spliced read), which is round((1 - identity)/2 * exon_bp).
"""
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt

# identity, diag_snp/copy, decisive_acc, decisive_frac, abstain_frac, attr_acc
ROWS = [
    (0.9995, 0,  None, 0.000, 1.000, 0.554),
    (0.999,  1,  0.624, 0.811, 0.165, 0.621),
    (0.995,  4,  0.954, 0.979, 0.000, 0.955),
    (0.99,   9,  0.970, 0.971, 0.029, 0.942),
    (0.98,   18, 0.926, 0.885, 0.014, 0.858),
    (0.97,   27, 0.930, 0.791, 0.209, 0.780),
    (0.95,   45, 0.749, 0.871, 0.123, 0.741),
]
ROWS.sort(key=lambda r: r[1])  # ascending diagnostic SNPs

snps = [r[1] for r in ROWS]
dec_acc = [r[2] for r in ROWS]
dec_frac = [r[3] for r in ROWS]
abstain = [r[4] for r in ROWS]

x = list(range(len(snps)))
xl = [str(s) for s in snps]

fig, ax = plt.subplots(figsize=(9, 5.5))
ax.plot(x, abstain, "o-", color="#d73027", lw=2.2, ms=7, label="abstain fraction (EM refuses to call)")
ax.plot(x, dec_frac, "s-", color="#4575b4", lw=2.0, ms=6, label="decisive fraction (EM commits)")
# decisive accuracy: skip the None (0 SNP has no decisive calls)
xa = [xi for xi, d in zip(x, dec_acc) if d is not None]
da = [d for d in dec_acc if d is not None]
ax.plot(xa, da, "^--", color="#1a9850", lw=2.0, ms=7, label="decisive accuracy (when it commits)")
ax.axhline(1/3, color="grey", ls=":", lw=1, label="chance (3 copies)")

ax.set_xticks(x); ax.set_xticklabels(xl)
ax.set_xlabel("copy-distinguishing EXONIC SNPs per copy  (← more identical)", fontsize=11)
ax.set_ylabel("fraction", fontsize=11)
ax.set_ylim(-0.03, 1.05)
ax.set_title("Per-copy read attribution: the identifiability spectrum\n"
             "fingerprint-EM on ambiguous (multimapping) synthetic reads, 3 near-identical copies",
             fontsize=12, weight="bold")
ax.legend(loc="center left", fontsize=9.5, framealpha=0.95)
ax.grid(alpha=0.25)

# annotate the two ends
ax.annotate("0 exonic diagnostics = the DAZ case\n→ 100% ABSTAIN, never fabricates",
            xy=(0, 1.0), xytext=(0.35, 0.86), fontsize=9, color="#d73027",
            arrowprops=dict(arrowstyle="->", color="#d73027"))
ax.annotate("≥4 diagnostics → 95–97% decisive accuracy\n(identifiable copies resolved)",
            xy=(2, 0.954), xytext=(2.4, 0.45), fontsize=9, color="#1a9850",
            arrowprops=dict(arrowstyle="->", color="#1a9850"))

fig.tight_layout()
out = "bench/tandem_attribution/attribution_spectrum.png"
fig.savefig(out, dpi=150, bbox_inches="tight")
print("wrote", out)
