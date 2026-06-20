#!/usr/bin/env python3
"""Compare two family-definition criteria on a Compara-labelled set: sequence SIMILARITY (minimizer-Jaccard,
the proxy the current 0.13 definition thresholds) vs the OPERATIONAL read-conflict (mutual-mappability) metric.

Reads bench/family_def_readconflict.tsv. Renders:
  A  minimizer-Jaccard by truth class -> the two classes OVERLAP (a domain-sharer scores highest); no threshold
     separates true paralogs from domain-sharers.
  B  read-conflict reads by truth class -> ZERO on every domain-sharer; fires only on the paralogs whose reads
     are actually confused (RABL2, APOBEC3). It is silent on RFPL because RFPL reads place uniquely (0 secondary
     / 0 MAPQ-0) -- a true paralog that poses no copy-assignment problem, correctly excluded from the
     resolution-needed set.

Run with /home/juanfra/miniforge3/bin/python."""
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

ROWS = []
for ln in open("/mnt/c/Users/jfris/Desktop/Rustle/bench/family_def_readconflict.tsv").read().splitlines()[1:]:
    f = ln.split("\t")
    ROWS.append(dict(label=f[0], a=f[1], b=f[2], n=int(f[3]), mut=float(f[4]),
                     ra=int(f[7]), rb=int(f[8]), jac=float(f[9])))

CLASSES = ["paralog", "domain-sharer"]
COL = {"paralog": "#2ca02c", "domain-sharer": "#d62728"}


def strip(ax, yfn, ylabel, title, logy=False, annot=None):
    rng = np.random.default_rng(3)
    for ci, cls in enumerate(CLASSES):
        pts = [r for r in ROWS if r["label"] == cls]
        ys = [yfn(r) for r in pts]
        xs = ci + (rng.random(len(pts)) - 0.5) * 0.28
        ax.scatter(xs, ys, c=COL[cls], s=90, edgecolor="black", linewidth=0.7, zorder=3, alpha=0.9)
        for x, y, r in zip(xs, ys, pts):
            if annot:
                lab = annot(r)
                if lab:
                    ax.annotate(lab, (x, y), fontsize=6.8, ha="center", va="bottom",
                                xytext=(0, 7), textcoords="offset points", color="#333")
    ax.set_xticks([0, 1]); ax.set_xticklabels([f"true paralogs\n(n={sum(r['label']=='paralog' for r in ROWS)})",
                                               f"domain-sharers\n(n={sum(r['label']=='domain-sharer' for r in ROWS)})"])
    ax.set_ylabel(ylabel); ax.set_title(title, fontsize=10, loc="left")
    if logy:
        ax.set_yscale("symlog", linthresh=1)
    ax.grid(axis="y", alpha=0.25)
    ax.set_xlim(-0.5, 1.5)


def main():
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(13.5, 6))
    # A: similarity overlaps
    jp = [r["jac"] for r in ROWS if r["label"] == "paralog"]
    jd = [r["jac"] for r in ROWS if r["label"] == "domain-sharer"]
    lo, hi = max(min(jp), min(jd)), min(max(jp), max(jd))
    ax1.axhspan(lo, hi, color="#999", alpha=0.13, zorder=0)
    strip(ax1, lambda r: r["jac"], "minimizer-Jaccard (sequence similarity)",
          "A.  Similarity does NOT separate the classes",
          annot=lambda r: f"{r['a']}~{r['b']}" if (r["jac"] > 0.29 or r["jac"] < 0.12) else None)
    ax1.text(0.5, hi + 0.01, "overlap band — no threshold separates\n(domain-sharer CREB1~METTL21A is the highest)",
             ha="center", fontsize=7.5, color="#555")
    # B: read-conflict separates
    strip(ax2, lambda r: r["n"], "read-conflict reads (AS-tied cross-mapping, symlog)",
          "B.  Read-conflict fires only where reads are actually confused", logy=True,
          annot=lambda r: (f"{r['a']}~{r['b']}" if r["n"] > 0 else
                           ("RFPL (reads resolve)" if r["a"].startswith("RFPL") and r["b"].startswith("RFPL")
                            and r["a"] == "RFPL1" and r["b"] == "RFPL2" else None)))
    ax2.set_ylim(-0.5, 2000)
    ax2.axhline(0.5, color="#d62728", ls=":", lw=1)
    ax2.text(1.0, 0.62, "every domain-sharer = 0", ha="center", fontsize=8, color="#d62728")
    fig.suptitle("Defining a multi-copy family at the RNA level: sequence similarity vs read-conflict "
                 "(Compara-labelled)", fontsize=12.5, fontweight="bold")
    fig.text(0.5, 0.005,
             "Read-conflict = reads with a primary in one gene and an AS-tied SECONDARY in the other (a real "
             "alternative placement, not one read overlapping two nested genes). It has ZERO false positives on "
             "domain-sharers, and flags exactly the families where copy-assignment is needed (RABL2 190, "
             "APOBEC3 6). It is silent on RFPL — a true paralog whose reads place uniquely (0 MAPQ-0) — i.e. it "
             "measures identifiability-RELEVANCE, not evolutionary paralogy.", ha="center", fontsize=8,
             wrap=True, style="italic")
    fig.tight_layout(rect=[0, 0.05, 1, 0.96])
    out = "/mnt/c/Users/jfris/Desktop/Rustle/bench/family_def_comparison.png"
    fig.savefig(out, dpi=140)
    print(f"[wrote {out}]")


if __name__ == "__main__":
    main()
