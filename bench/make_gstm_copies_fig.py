#!/usr/bin/env python
"""Real-data illustrative example: the gorilla GSTM multi-copy family.

Shows, from ACTUAL copy_assign output on GGO_mm.bam (region NC_073224.2:129.16-129.23M,
--homology-primary --dump-psv), that the 3 called copies are the annotated GST-Mu paralogs
GSTM3/GSTM5/GSTM1, and that the 2673 reads sort cleanly into 3 copy-blocks by their PSV alleles.

Panel A  read x PSV genotype matrix (SDA/Vollger-style): rows = reads (grouped by assigned copy),
         columns = informative PSV positions, cells coloured by base. Three clean blocks = three copies.
Panel B  corroboration: copy -> annotated gene, per-copy private SUN markers, unique-mapper agreement,
         and the measured pairwise similarity of the clean pair.

Inputs (winloci_scratch): gstm_vg.psv_{reads,copies,cols}.tsv, gstm_vg.families.tsv
Run: /home/juanfra/miniforge3/bin/python bench/make_gstm_copies_fig.py
"""
import os
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Patch, FancyBboxPatch
import numpy as np

SC = "/home/juanfra/winloci_scratch"
OUT = os.path.join(os.path.dirname(__file__), "slides")
os.makedirs(OUT, exist_ok=True)
plt.rcParams.update({"font.family": "DejaVu Sans"})

GENES = ["GSTM3", "GSTM5", "GSTM1"]           # copy 0/1/2 by genomic order
CCOL = ["#e41a1c", "#377eb8", "#4daf4a"]      # per-copy colour
BASE = {"A": "#3cb44b", "C": "#4363d8", "G": "#f58231", "T": "#e6194b", ".": "#f0f0f0", "-": "#f0f0f0"}


def load():
    copies = []
    for l in open(f"{SC}/gstm_vg.psv_copies.tsv"):
        if l.startswith("family_id"):
            continue
        _, ci, tid, al = l.rstrip("\n").split("\t")
        copies.append((int(ci), al))
    copies.sort()
    cop_al = [c[1] for c in copies]
    reads = []
    for l in open(f"{SC}/gstm_vg.psv_reads.tsv"):
        if l.startswith("read_name"):
            continue
        p = l.rstrip("\n").split("\t")
        reads.append((int(p[2]) if p[2].isdigit() else -1, p[3], p[6]))  # assigned_copy, status, alleles
    return cop_al, reads


def main():
    cop_al, reads = load()
    ncol = len(cop_al[0])

    # informative columns: all 3 copies defined AND pairwise differ (true PSVs distinguishing copies)
    info = [k for k in range(ncol)
            if all(cop_al[i][k] != "." for i in range(3)) and len({cop_al[i][k] for i in range(3)}) >= 2]
    # spread ~70 across the transcript
    sel = info[:: max(1, len(info) // 70)][:70]

    # reads with enough coverage over selected cols, assigned to a copy
    def cov(al):
        return sum(al[k] != "." for k in sel)
    good = [(c, st, al) for (c, st, al) in reads if c in (0, 1, 2) and cov(al) >= 6]
    good.sort(key=lambda r: (r[0], -cov(r[2])))
    # subsample per copy for rendering (keep proportions readable)
    per = {0: [], 1: [], 2: []}
    for c, st, al in good:
        per[c].append(al)
    MAXP = 260
    # top strip: the 3 COPY CONSENSUS rows (each copy's own allele), then a gap, then the reads
    rows, rowcopy = [], []
    for c in (0, 1, 2):
        rows.append([cop_al[c][k] for k in sel]); rowcopy.append(("cons", c))
    GAP = 4
    for _ in range(GAP):
        rows.append(["." for _ in sel]); rowcopy.append(("gap", -1))
    n_consensus = 3 + GAP
    for c in (0, 1, 2):
        chosen = per[c][:MAXP]
        for al in chosen:
            rows.append([al[k] for k in sel])
            rowcopy.append(("read", c))

    # base -> int grid for imshow
    palette = ["A", "C", "G", "T", "."]
    idx = {b: i for i, b in enumerate(palette)}
    grid = np.array([[idx.get(b, 4) for b in r] for r in rows])
    from matplotlib.colors import ListedColormap
    cmap = ListedColormap([BASE[b] for b in palette])

    # SUNs per copy + stats
    suns = [0, 0, 0]
    for k in range(ncol):
        bases = [cop_al[i][k] for i in range(3) if cop_al[i][k] != "."]
        for i in range(3):
            b = cop_al[i][k]
            if b != "." and len(bases) >= 2 and sum(1 for x in bases if x == b) == 1:
                suns[i] += 1
    fam = [l.rstrip().split("\t") for l in open(f"{SC}/gstm_vg.families.tsv")]
    h = fam[0]; d = dict(zip(h, fam[1]))
    n_assigned = int(d["assigned_j"]); n_reads = int(d["n_reads"])
    uniq = int(d["uniq"]); uniq_agree = int(d["uniq_agree"])

    fig = plt.figure(figsize=(15, 8.2))
    gs = fig.add_gridspec(1, 2, width_ratios=[2.35, 1.0], wspace=0.16)
    axH = fig.add_subplot(gs[0]); axR = fig.add_subplot(gs[1]); axR.axis("off")

    axH.imshow(grid, aspect="auto", cmap=cmap, interpolation="nearest")
    nrow = len(rows)
    # top: 3 copy-consensus rows (colour tick each) + a "copies" bracket
    for c in (0, 1, 2):
        axH.add_patch(plt.Rectangle((-2.6, c), 1.8, 1, color=CCOL[c], clip_on=False))
    axH.text(-4.4, 1.0, "copies", rotation=90, va="center", ha="center", fontsize=8.5, weight="bold", color="#333")
    axH.axhline(n_consensus - 0.5, color="#222", lw=1.3)  # separator: consensus above, reads below
    # read blocks, coloured by assigned copy
    b0 = n_consensus
    for c in (0, 1, 2):
        n = sum(1 for rc in rowcopy if rc == ("read", c))
        axH.add_patch(plt.Rectangle((-2.6, b0), 1.8, n, color=CCOL[c], clip_on=False))
        axH.text(-4.4, b0 + n / 2, GENES[c], rotation=90, va="center", ha="center",
                 fontsize=11, weight="bold", color=CCOL[c])
        b0 += n
    axH.set_xticks([]); axH.set_yticks([])
    axH.set_xlabel(f"{len(sel)} of {len(info)} PSV columns (copy-distinguishing positions)", fontsize=10)
    axH.set_ylabel(f"top: 3 copy consensus  ·  below: {nrow - n_consensus} reads (of {n_reads}) by assigned copy", fontsize=9.5)
    axH.set_title("Real gorilla reads sort into 3 copies by their PSV alleles",
                  fontsize=13.5, weight="bold", pad=10)
    axH.legend(handles=[Patch(fc=BASE[b], ec="#999", label=b) for b in "ACGT"],
               loc="upper right", ncol=4, fontsize=8, title="read base", title_fontsize=8,
               framealpha=0.9)

    # ---- right: corroboration ----
    def line(y, s, **kw):
        axR.text(0.0, y, s, transform=axR.transAxes, va="top", fontsize=kw.pop("fs", 10.5), **kw)
    axR.add_patch(FancyBboxPatch((-0.02, 0.02), 1.02, 0.97, transform=axR.transAxes,
                                 boxstyle="round,pad=0.01,rounding_size=0.02", fc="#f4f7fb", ec="#14285a", lw=1.4))
    line(0.955, "GSTM — a real multi-copy family", weight="bold", fs=13.5, color="#14285a")
    line(0.895, "gorilla GGO_mm.bam · NC_073224.2:129.17–129.22 Mb", fs=8.8, color="#666", style="italic")

    line(0.855, "The 3 called copies ARE annotated paralogs:", weight="bold", fs=10.5)
    for i, g in enumerate(GENES):
        axR.add_patch(plt.Rectangle((0.02, 0.788 - i * 0.05), 0.028, 0.028, transform=axR.transAxes, color=CCOL[i]))
        line(0.807 - i * 0.05, f"    copy {i} → {g}   ({suns[i]} private SUN markers)", fs=9.8)

    line(0.63, "Grouped into ONE family by homology", weight="bold", fs=10.5)
    line(0.585, "   (exon-sum identity ≥ 80% — the membership gate;", fs=9.5, color="#333")
    line(0.545, "    GSTM5↔GSTM1 measured 83%; GSTM3 the divergent", fs=9.5, color="#333")
    line(0.505, "    member — 15% PSV-site agreement).", fs=9.5, color="#333")

    line(0.40, "Why believe the assignment:", weight="bold", fs=10.5, color="#14285a")
    line(0.355, f"   • unique-mapper agreement  {uniq_agree}/{uniq}  = 100%", fs=10, color="#1a7d1a", weight="bold")
    line(0.315, "     (where minimap2 is confident, we agree)", fs=8.8, color="#666")
    line(0.265, f"   • {n_assigned}/{n_reads} reads assigned; 16 certified TIED", fs=10)
    line(0.225, "     (abstain, never 1/k)", fs=8.8, color="#666")
    line(0.175, f"   • every copy carries 100s of private SUNs →", fs=10)
    line(0.135, "     one read can deterministically pin its copy", fs=8.8, color="#666")

    line(0.06, "A PSV column = a bubble; a copy = a path; a read", fs=8.6, color="#555", style="italic")
    line(0.028, "follows its alleles to a copy or abstains.", fs=8.6, color="#555", style="italic")

    p = os.path.join(OUT, "gstm_real_copies.png")
    fig.savefig(p, bbox_inches="tight", facecolor="white", dpi=150)
    plt.close()
    print("wrote", p, "|", len(rows), "reads x", len(sel), "cols; SUNs", dict(zip(GENES, suns)),
          "; uniq-agree", f"{uniq_agree}/{uniq}")


if __name__ == "__main__":
    main()
