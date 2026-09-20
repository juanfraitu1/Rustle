#!/usr/bin/env python
"""Cross-species generalization — the IDENTICAL method + code on HUMAN testis Iso-Seq.

Rebuts "overfit to gorilla / to this sample": we ran the same copy_assign on a human testis Iso-Seq (ERR13885926,
GENCODE, Sequel II HiFi — a different lab, individual, and species) aligned to T2T-CHM13v2.0 with the same
-N 50 --secondary=yes recipe, at the human orthologs of the families validated in gorilla. The method reports the
REAL, species-specific human copy numbers (matching the human annotation) — it tracks the biology, not the sample.

Measured this session (see HUMAN_CROSSSPECIES.md). Output: bench/slides/human_crossspecies.png
Run: /home/juanfra/miniforge3/bin/python bench/make_human_crossspecies.py
"""
import os
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

OUT = os.path.join(os.path.dirname(__file__), "slides")
os.makedirs(OUT, exist_ok=True)
plt.rcParams.update({"font.family": "DejaVu Sans"})
GORILLA, HUMAN, ANNOT, GREY = "#8c6d1f", "#377eb8", "#4daf4a", "#bbbbbb"
NAVY = "#14285a"

# (family, gorilla chi(H), human chi(H), human annot copies, human reads, note)  human=None => coverage-limited
DATA = [
    ("RBMY",  6, 6,  6,  147, ""),
    ("TSPY",  5, 33, 35, 683, ""),
    ("MAGEA", 2, 11, 11, 3241, ""),
    ("GSTM",  3, 2,  5,  2337, "partial (GSTM2/5 expressed)"),
    ("PCDHB", 5, None, 16, 77, "77 reads — under-expressed here"),
    ("DAZ",   2, None, 4,  16, "16 reads — under-expressed here"),
]


def build():
    fig, ax = plt.subplots(figsize=(12.6, 6.0))
    fams = [d[0] for d in DATA]
    x = np.arange(len(fams))
    w = 0.27
    g = [d[1] for d in DATA]
    h = [d[2] if d[2] is not None else 0 for d in DATA]
    a = [d[3] for d in DATA]
    b1 = ax.bar(x - w, g, w, color=GORILLA, ec="black", lw=.6, label="gorilla  χ(H)  (our RNA)")
    b2 = ax.bar(x,     h, w, color=HUMAN,   ec="black", lw=.6, label="HUMAN  χ(H)  (our RNA, same code)")
    b3 = ax.bar(x + w, a, w, color="none",  ec=ANNOT, lw=2, ls="--",
                label="human annotation (T2T)")
    ax.bar_label(b1, fontsize=10, weight="bold", color=GORILLA, padding=1)
    for i, d in enumerate(DATA):
        if d[2] is not None:
            ax.text(x[i], h[i] + 0.6, str(h[i]), ha="center", fontsize=11.5, weight="bold", color=HUMAN)
        else:
            ax.text(x[i], 1.2, "n/a", ha="center", fontsize=10, weight="bold", color=GREY, rotation=0)
    ax.bar_label(b3, fontsize=9.5, color="#2e7d32", padding=1)
    ax.set_xticks(x); ax.set_xticklabels(fams, fontsize=13, weight="bold")
    ax.set_ylabel("copy number", fontsize=12)
    ax.set_ylim(0, 40)
    ax.legend(loc="upper right", fontsize=11, framealpha=.95)
    ax.set_title("Same method, same code — run on HUMAN testis Iso-Seq (T2T-CHM13). It tracks the species, not the sample.",
                 fontsize=13.5, weight="bold", color=NAVY, pad=12)
    ax.spines["top"].set_visible(False); ax.spines["right"].set_visible(False)
    # key message in clear space (below the legend, above the low bars)
    ax.text(2.35, 26.5, "MAGEA 2→11  ·  TSPY 5→33:\nthe same code reports the REAL human\nexpansions (≈ the human annotation),\nnot the gorilla numbers",
            fontsize=11, color=HUMAN, weight="bold", va="top")
    ax.text(0.0, -5.6, "Human = ERR13885926 (GENCODE testis, Sequel II HiFi) — different lab, individual, species. "
            "PCDHB/DAZ n/a = coverage-limited in this single library (77 / 16 reads), an honest depth limit, not a method failure.",
            fontsize=9, style="italic", color="#777", transform=ax.transData)
    p = os.path.join(OUT, "human_crossspecies.png")
    fig.savefig(p, bbox_inches="tight", facecolor="white", dpi=150)
    plt.close()
    return p


if __name__ == "__main__":
    print("wrote", build())
