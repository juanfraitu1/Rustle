#!/usr/bin/env python
"""Concordance with Soto et al. 2025 (Cell, 'Human-specific gene expansions contribute to brain evolution').

A two-sided, cross-lab, cross-modal (their DNA/T2T vs our RNA) validation on THEIR gene set:
  SPECIFICITY  — their human-specific expansions -> our RNA method calls single-copy in gorilla (we AGREE they
                 are not gorilla families; we do not fabricate them). All loci are expressed (read counts shown).
  RECOVERY     — ancestral / shared gorilla families -> recovered at concordant copy number; 3 are genuine SD98
                 (>=98% identity) segdups, the exact objects Soto's SD98+famCN pipeline clusters.

All numbers measured this session on GGO_mm.bam / GGO.fasta / GGO_genomic.gff / SEDEF final.bed (see SOTO_CONCORDANCE.md).
Output: bench/slides/soto_concordance.png
Run: /home/juanfra/miniforge3/bin/python bench/make_soto_concordance.py
"""
import os
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch

OUT = os.path.join(os.path.dirname(__file__), "slides")
os.makedirs(OUT, exist_ok=True)
plt.rcParams.update({"font.family": "DejaVu Sans"})
GREEN, BLUE, GREY, GOLD, NAVY = "#4daf4a", "#377eb8", "#999999", "#d4a017", "#14285a"

# SPECIFICITY: Soto human-specific families -> gorilla single-copy (our RNA), with testis read coverage.
# (family, human paralog Soto studies, gorilla ortholog, our RNA call, reads)  reads None => absent in gorilla
SPEC = [
    ("GPR89",   "GPR89B*",   "GPR89A",     "single-copy", 70),
    ("FRMPD2",  "FRMPD2B*",  "FRMPD2",     "single-copy", 363),
    ("CD8B",    "CD8B2*",    "CD8B",       "single-copy", 4),
    ("SRGAP2",  "SRGAP2B/C/D", "SRGAP2",   "single-copy", 248),
    ("ARHGAP11","ARHGAP11B", "ARHGAP11A",  "single-copy", 310),
    ("HYDIN",   "HYDIN2",    "HYDIN",      "single-copy", 187),
    ("ROCK1",   "ROCK1P1",   "ROCK1",      "single-copy", 579),
    ("FAM72",   "FAM72B/C/D","FAM72A",     "single-copy", None),   # covered within SRGAP2 locus
    ("GPRIN2",  "GPRIN2B",   "GPRIN2",     "single-copy", 8),
    ("DUSP22",  "DUSP22B",   "DUSP22",     "single-copy", 277),
    ("NPY4R",   "NPY4R2",    "—",          "absent",      None),
    ("CFC1",    "CFC1B",     "—",          "absent",      None),
    ("NOTCH2NL","NOTCH2NL",  "—",          "absent",      None),
]
# RECOVERY: ancestral gorilla families -> our chi(H); SD98 = max SEDEF identity >= 98% (Soto clusters these)
REC = [
    ("RBMY",  6, 99.74, True),
    ("TSPY",  5, 99.45, True),
    ("DAZ",   2, 99.63, True),
    ("GSTM",  3, 95.29, False),
    ("MAGEA", 2, 94.85, False),
    ("PCDHB", 5, 88.86, False),
]


def build():
    fig = plt.figure(figsize=(13.2, 6.7))
    gs = fig.add_gridspec(1, 2, width_ratios=[1.28, 1.0], wspace=0.12)
    axL = fig.add_subplot(gs[0]); axR = fig.add_subplot(gs[1])
    for ax in (axL, axR):
        ax.axis("off"); ax.set_xlim(0, 1); ax.set_ylim(0, 1)

    # ---- LEFT: SPECIFICITY
    axL.add_patch(FancyBboxPatch((0.0, 0.90), 1.0, 0.085, boxstyle="round,pad=0.005,rounding_size=0.02",
                                 fc=GREEN, ec=GREEN))
    axL.text(0.5, 0.943, "SPECIFICITY:  human-specific  →  single-copy in gorilla",
             ha="center", va="center", color="white", fontsize=11.5, weight="bold")
    n = len(SPEC); y0, y1 = 0.855, 0.055
    dy = (y0 - y1) / n
    axL.text(0.02, 0.878, "gene family", fontsize=8.3, color="#888", weight="bold")
    axL.text(0.41, 0.878, "human paralog (Soto)", fontsize=8.3, color="#888", weight="bold", ha="center")
    axL.text(0.72, 0.878, "our RNA — gorilla", fontsize=8.3, color="#888", weight="bold", ha="center")
    for i, (fam, hpar, gorth, call, reads) in enumerate(SPEC):
        y = y0 - (i + 0.5) * dy
        absent = call == "absent"
        col = GREY if absent else GREEN
        axL.text(0.02, y, fam, fontsize=10.5, weight="bold", va="center", color=NAVY)
        axL.text(0.41, y, hpar, fontsize=9, va="center", ha="center", color="#333")
        badge = "absent (0)" if absent else "✓  single-copy"
        axL.add_patch(FancyBboxPatch((0.58, y - 0.5 * dy + 0.006), 0.28, dy - 0.012,
                      boxstyle="round,pad=0.003,rounding_size=0.02", fc=col, ec=col, alpha=0.92))
        axL.text(0.72, y, badge, fontsize=9.3, va="center", ha="center", color="white", weight="bold")
        rtxt = "" if reads is None else f"{reads} reads"
        axL.text(0.995, y, rtxt, fontsize=8.2, va="center", ha="right", color="#999")
    axL.text(0.02, 0.028, "13/13 concordant — not gorilla families (we agree)",
             fontsize=9.5, color=GREEN, weight="bold")
    axL.text(0.02, -0.002, "* GPR89B, FRMPD2B, CD8B2 = Soto's functionally-validated genes (brain size · synapse · T-cell selection)",
             fontsize=8, style="italic", color="#888")

    # ---- RIGHT: RECOVERY
    axR.add_patch(FancyBboxPatch((0.0, 0.90), 1.0, 0.085, boxstyle="round,pad=0.005,rounding_size=0.02",
                                 fc=BLUE, ec=BLUE))
    axR.text(0.5, 0.943, "RECOVERY:  ancestral gorilla families, resolved",
             ha="center", va="center", color="white", fontsize=11.5, weight="bold")
    axR.text(0.03, 0.862, "family", fontsize=8.3, color="#888", weight="bold")
    axR.text(0.40, 0.862, "our RNA χ(H)", fontsize=8.3, color="#888", weight="bold", ha="center")
    axR.text(0.985, 0.862, "gorilla segdup identity", fontsize=8.3, color="#888", weight="bold", ha="right")
    m = len(REC); ry0 = 0.815
    rdy = (ry0 - 0.10) / m
    for i, (fam, cn, ident, sd98) in enumerate(REC):
        y = ry0 - (i + 0.5) * rdy
        axR.text(0.03, y, fam, fontsize=13, weight="bold", va="center", color=NAVY)
        axR.text(0.40, y, f"{cn}", fontsize=17, weight="bold", va="center", ha="center", color=BLUE)
        axR.text(0.475, y, "copies", fontsize=10, va="center", color=BLUE)
        tag = f"SD98  {ident:.1f}%" if sd98 else f"paralog  {ident:.1f}%"
        tcol = GOLD if sd98 else GREY
        axR.add_patch(FancyBboxPatch((0.66, y - 0.5 * rdy + 0.01), 0.325, rdy - 0.02,
                      boxstyle="round,pad=0.003,rounding_size=0.02", fc="none", ec=tcol, lw=1.6))
        axR.text(0.8225, y, tag, fontsize=9.5, va="center", ha="center", color=tcol, weight="bold")
    axR.text(0.5, 0.03,
             "DAZ · RBMY · TSPY are genuine SD98 segdups (the objects Soto clusters);\nGSTM · MAGEA · PCDHB are more-divergent paralog families our PSV method still resolves.",
             ha="center", fontsize=8.4, style="italic", color="#888")

    fig.suptitle("Concordance with Soto et al. 2025 (Cell) — their DNA/T2T human-specific catalog vs our gorilla RNA method",
                 fontsize=13.5, weight="bold", color=NAVY, y=1.02)
    p = os.path.join(OUT, "soto_concordance.png")
    fig.savefig(p, bbox_inches="tight", facecolor="white", dpi=150)
    plt.close()
    return p


if __name__ == "__main__":
    print("wrote", build())
