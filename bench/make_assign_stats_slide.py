#!/usr/bin/env python
"""Slide: what it means to assign a read to a copy, and the statistical test that decides it.

Left  — the mechanism: a MAPQ-0 read matches one copy's allele fingerprint at the distinguishing PSVs.
Right — the IsoCon significance certificate (copy_assign.rs): p = Poisson-binomial upper tail of matches
        under H0 "read came from the competitor"; assign iff p < alpha/(K-1); else certify TIED, never 1/k.
Run: /home/juanfra/miniforge3/bin/python bench/make_assign_stats_slide.py
"""
import os
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, Rectangle, FancyArrowPatch

OUT = os.path.join(os.path.dirname(__file__), "slides")
os.makedirs(OUT, exist_ok=True)
plt.rcParams.update({"font.family": "DejaVu Sans"})
NAVY = "#14285a"; BLUE = "#377eb8"; RED = "#e41a1c"; GREEN = "#1e8449"; GREY = "#8a8f98"; DARK = "#222"
BASE = {"A": "#3cb44b", "C": "#4363d8", "G": "#f58231", "T": "#e6194b"}


def build():
    fig = plt.figure(figsize=(15.5, 8.6))
    gs = fig.add_gridspec(1, 2, width_ratios=[1.03, 1.0], wspace=0.13)
    axL = fig.add_subplot(gs[0]); axR = fig.add_subplot(gs[1])
    for ax in (axL, axR):
        ax.axis("off"); ax.set_xlim(0, 10); ax.set_ylim(0, 10)

    fig.suptitle("Assigning a read to a copy — and the test that decides it",
                 fontsize=17, weight="bold", color=NAVY, y=0.99)

    # ================= LEFT: the mechanism =================
    axL.text(0.2, 9.5, "The question", fontsize=12.5, weight="bold", color=NAVY)
    axL.text(0.2, 9.0, "A read (one transcript molecule) fits several near-identical", fontsize=10.2, color=DARK)
    axL.text(0.2, 8.65, "copies equally well → MAPQ 0.  Which paralog produced it?", fontsize=10.2, color=DARK)
    axL.text(0.2, 8.2, "Not where it aligns — which COPY is its source.", fontsize=9.6, style="italic", color="#555")

    # worked matrix: copy A / copy B fingerprints + the READ, over 6 distinguishing PSVs
    cols = ["A", "C", "G", "T", "A", "C"]      # copy A alleles at 6 distinguishing PSVs
    copB = ["G", "T", "A", "C", "G", "T"]      # copy B alleles (all differ)
    read = ["A", "C", "G", "T", "A", "C"]      # observed read bases (match A)
    x0, y0, cw, ch = 1.6, 4.6, 0.92, 0.62
    labels = [("copy A", cols, BLUE), ("copy B", copB, RED), ("the read", read, "#333")]
    axL.text(x0 + 3 * cw, 7.65, "distinguishing PSV bubbles (where the copies differ)", ha="center", fontsize=8.8, color="#555")
    for r, (name, seq, col) in enumerate(labels):
        yy = y0 + (2 - r) * (ch + 0.12)
        axL.text(x0 - 0.2, yy + ch / 2, name, ha="right", va="center", fontsize=10.5, weight="bold", color=col)
        for c, b in enumerate(seq):
            axL.add_patch(Rectangle((x0 + c * cw, yy), cw - 0.08, ch, fc=BASE[b], ec="white", lw=1.5))
            axL.text(x0 + c * cw + (cw - 0.08) / 2, yy + ch / 2, b, ha="center", va="center",
                     fontsize=11, weight="bold", color="white")
        if name == "the read":
            for c in range(len(seq)):
                axL.text(x0 + c * cw + (cw - 0.08) / 2, yy - 0.32, "✓", ha="center", fontsize=11, weight="bold", color=GREEN)
    axL.text(x0 + 3 * cw, y0 - 0.95, "the read agrees with copy A at every distinguishing PSV",
             ha="center", fontsize=9.6, weight="bold", color=GREEN)
    axL.text(0.2, 2.7, "Each copy is a PATH through the variation graph (a fixed allele", fontsize=10, color=DARK)
    axL.text(0.2, 2.35, "at each bubble); the read is a PARTIAL path (its observed", fontsize=10, color=DARK)
    axL.text(0.2, 2.0, "alleles). Assignment = which copy-path the read matches.", fontsize=10, color=DARK)
    axL.text(0.2, 1.35, "k = # distinguishing features the read matches the best copy", fontsize=9.6, color=NAVY, weight="bold")
    axL.text(0.2, 1.0, "D = # distinguishing features the read spans   (here k = D = 6)", fontsize=9.6, color=NAVY, weight="bold")

    # ================= RIGHT: the statistics =================
    axR.add_patch(FancyBboxPatch((0.1, 0.2), 9.8, 9.35, boxstyle="round,pad=0.02,rounding_size=0.06",
                                 fc="#f4f7fb", ec=NAVY, lw=1.6))
    axR.text(0.45, 9.15, "The significance test  (IsoCon certificate)", fontsize=12.5, weight="bold", color=NAVY)

    axR.text(0.45, 8.5, "H₀ :  the read came from the competitor copy B.", fontsize=10.6, weight="bold", color=DARK)
    axR.text(0.45, 8.02, "Under H₀, each match to copy A at a distinguishing position", fontsize=9.8, color="#333")
    axR.text(0.45, 7.66, "must be a sequencing error:", fontsize=9.8, color="#333")
    axR.text(0.9, 7.2, "εⱼ = e/3 ≈ 10⁻³   (HiFi base;  from the read's phred Q)", fontsize=9.8, color=RED, weight="bold")
    axR.text(0.9, 6.82, "εⱼ = 10⁻⁴          (a copy-specific junction)", fontsize=9.8, color=RED, weight="bold")

    axR.text(0.45, 6.2, "Test statistic (per read, best copy A vs its hardest rival):", fontsize=10, weight="bold", color=NAVY)
    axR.add_patch(FancyBboxPatch((0.6, 5.0), 8.9, 1.0, boxstyle="round,pad=0.02,rounding_size=0.04",
                                 fc="white", ec=GREY, lw=1.2))
    axR.text(5.05, 5.72, "p  =  P( ≥ k of the D features match A by error )",
             ha="center", fontsize=10.6, weight="bold", color=DARK)
    axR.text(5.05, 5.24, "=  Poisson-binomial upper tail (k ; ε₁ … ε_D)   ≈  εᵏ",
             ha="center", fontsize=10.0, color=DARK)

    axR.text(0.45, 4.5, "Decision:", fontsize=10.5, weight="bold", color=NAVY)
    axR.text(0.9, 4.06, "ASSIGN to A   ⟺   p  <  α / (K−1),   α = 10⁻³", fontsize=10.6, weight="bold", color=GREEN)
    axR.text(0.95, 3.7, "(Bonferroni over the K−1 rival copies; log-LR margin > 0)", fontsize=8.9, color="#555")
    axR.text(0.9, 3.2, "ABSTAIN → certify TIED   if the read spans 0 distinguishing", fontsize=10.2, weight="bold", color=GREY)
    axR.text(0.95, 2.85, "features, or cannot reach α even if perfect  (Πεⱼ ≥ α):", fontsize=9.6, color="#555")
    axR.text(0.95, 2.5, "the identifiability floor (K = 0).   Never split 1/k.", fontsize=9.6, color="#555", style="italic")

    # real example
    axR.add_patch(FancyBboxPatch((0.5, 0.5), 9.0, 1.65, boxstyle="round,pad=0.02,rounding_size=0.04",
                                 fc="#eef6ee", ec=GREEN, lw=1.3))
    axR.text(0.75, 1.86, "Real GSTM read → copy 0:  k = 414 → p ≈ 10⁻¹²⁰⁰ ≪ α  →  ASSIGNED",
             fontsize=9.6, weight="bold", color=GREEN)
    axR.text(0.75, 1.42, "Read spanning 0 distinguishing features  →  TIED (abstain)", fontsize=9.6, color="#555")
    axR.text(0.75, 0.92, "2654 assigned · 16 tied · unique-mapper agreement 100% (1341/1341)",
             fontsize=9.4, weight="bold", color=DARK)

    p = os.path.join(OUT, "assign_stats.png")
    fig.savefig(p, bbox_inches="tight", facecolor="white", dpi=150)
    plt.close()
    print("wrote", p)


if __name__ == "__main__":
    build()
