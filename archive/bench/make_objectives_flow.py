#!/usr/bin/env python
"""One figure, four objectives: reads -> variation graph -> result, each row a real gorilla example.

O1  reads at homologous loci        -> homology graph -> family + copy number chi_H   (GSTM, real)
O2  reads thread the PSV bubbles     -> variation graph -> assigned copy / tied         (GSTM, real)
O3  reads split by an allele         -> a junction      -> allele-specific splicing      (ABCC4, real)
O4  reads pile up beyond the ref     -> variation graph -> a copy NOT in the assembly    (DAZ, real)

Numbers are from copy_assign / asj on GGO_mm.bam (see GSTM_REAL_COPIES.md, OBJECTIVES_FLOW.md).
Run: /home/juanfra/miniforge3/bin/python bench/make_objectives_flow.py
"""
import json, os
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch, Circle, Rectangle

SC = "/home/juanfra/winloci_scratch"
OUT = os.path.join(os.path.dirname(__file__), "slides")
os.makedirs(OUT, exist_ok=True)
plt.rcParams.update({"font.family": "DejaVu Sans"})
NAVY = "#14285a"; BLUE = "#377eb8"; RED = "#e41a1c"; GREEN = "#4daf4a"; GREY = "#8a8f98"
ORANGE = "#e08a00"; DARK = "#222"; SPINE = "#cccccc"; CC = [RED, BLUE, GREEN]


def o4_numbers():
    d = {"n_copies": None, "chi_H": None, "depth_cn": None, "dna_needs": 0, "n_reads": None}
    fp = f"{SC}/daz_o4.famcn_readonly.tsv"
    if os.path.exists(fp):
        rows = [l.rstrip("\n").split("\t") for l in open(fp)]
        h = rows[0]
        if len(rows) > 1:
            r = dict(zip(h, rows[1]))
            d.update(n_copies=r.get("n_copies"), chi_H=r.get("chi_H"),
                     depth_cn=r.get("depth_cn"), n_reads=r.get("n_reads"))
    dn = f"{SC}/daz_o4.dna_needs.tsv"
    if os.path.exists(dn):
        d["dna_needs"] = max(0, sum(1 for _ in open(dn)) - 1)
    return d


def stage(ax, x, y, w, h, title, tcolor):
    ax.add_patch(FancyBboxPatch((x, y), w, h, boxstyle="round,pad=0.01,rounding_size=0.03",
                                fc="white", ec=tcolor, lw=1.8, zorder=2))
    ax.text(x + w / 2, y + h - 0.16, title, ha="center", va="top", fontsize=9.5, weight="bold", color=tcolor)


def arrow(ax, x0, x1, y, label=""):
    ax.add_patch(FancyArrowPatch((x0, y), (x1, y), arrowstyle="-|>", mutation_scale=16, lw=2, color=GREY, zorder=3))
    if label:
        ax.text((x0 + x1) / 2, y + 0.12, label, ha="center", va="bottom", fontsize=7.8, style="italic", color="#555")


def pileup(ax, cx, cy, n, color, w=0.5):
    for i in range(min(n, 7)):
        ax.plot([cx - w, cx + w], [cy + i * 0.14, cy + i * 0.14], color=color, lw=2.4, solid_capstyle="round")


def build():
    o4 = o4_numbers()
    fig, axes = plt.subplots(4, 1, figsize=(15.5, 15.2))
    for ax in axes:
        ax.axis("off"); ax.set_xlim(0, 16); ax.set_ylim(0, 3.2)
    XA, XB, XC = 0.6, 6.4, 11.4      # stage x-origins
    WA, WB, WC = 4.2, 4.0, 4.2

    # ================= O1 — reads -> homology -> family + chi_H =================
    ax = axes[0]
    ax.text(0.1, 3.05, "O1", fontsize=20, weight="bold", color=NAVY)
    ax.text(0.85, 3.06, "define the family + count copies", fontsize=11.5, color=DARK, va="center")
    stage(ax, XA, 0.35, WA, 2.35, "READS at 3 loci", NAVY)
    for i, (g, nm) in enumerate([("GSTM3", 187), ("GSTM5", "…"), ("GSTM1", "…")]):
        pileup(ax, XA + 0.9 + i * 1.15, 0.8, 6, CC[i])
        ax.text(XA + 0.9 + i * 1.15, 0.55, g, ha="center", fontsize=8.5, weight="bold", color=CC[i])
    ax.text(XA + WA / 2, 2.18, "2673 reads · 3 homologous loci", ha="center", fontsize=8.2, color="#555")
    arrow(ax, XA + WA, XB, 1.5, "homology ≥ 80%")
    stage(ax, XB, 0.35, WB, 2.35, "HOMOLOGY GRAPH", GREEN)
    nodes = [(XB + 1.1, 1.9), (XB + 2.9, 2.0), (XB + 2.0, 0.95)]
    for a in range(3):
        for b in range(a + 1, 3):
            ax.plot([nodes[a][0], nodes[b][0]], [nodes[a][1], nodes[b][1]], color=GREEN, lw=1.8, zorder=2)
    for i, (x, y) in enumerate(nodes):
        ax.add_patch(Circle((x, y), 0.22, fc=CC[i], ec="black", lw=1, zorder=4))
    ax.add_patch(FancyBboxPatch((XB + 0.6, 0.6), 2.8, 1.8, boxstyle="round,pad=0.02", fc="none",
                                ec=GREEN, lw=1.6, ls=(0, (4, 2)), zorder=1))
    ax.text(XB + WB / 2, 0.5, "one family (γ-quasi-clique)", ha="center", fontsize=8.2, color="#555", style="italic")
    arrow(ax, XB + WB, XC, 1.5, "colour conflict graph")
    stage(ax, XC, 0.35, WC, 2.35, "COPY NUMBER", RED)
    for i in range(3):
        ax.add_patch(Rectangle((XC + 0.6, 1.9 - i * 0.55), 3.0, 0.34, fc=CC[i], ec="black", lw=0.8))
        ax.text(XC + 0.4, 2.07 - i * 0.55, f"c{i}", ha="right", va="center", fontsize=9, weight="bold", color=CC[i])
    ax.text(XC + WC / 2, 0.62, "χ_H = 3   (= GSTM3, GSTM5, GSTM1)", ha="center", fontsize=10.5, weight="bold", color=RED)

    # ================= O2 — reads -> variation graph -> assigned copy =================
    ax = axes[1]
    ax.text(0.1, 3.05, "O2", fontsize=20, weight="bold", color=NAVY)
    ax.text(0.85, 3.06, "assign each read to a copy (or abstain)", fontsize=11.5, color=DARK, va="center")
    stage(ax, XA, 0.35, WA, 2.35, "READS", NAVY)
    for i in range(4):
        ax.plot([XA + 0.7, XA + 3.5], [2.1 - i * 0.42, 2.1 - i * 0.42], color="#666", lw=3, solid_capstyle="round")
    ax.text(XA + WA / 2, 0.55, "multi-mapping reads (MAPQ 0)", ha="center", fontsize=8.2, color="#555")
    arrow(ax, XA + WA, XB, 1.5, "thread PSV bubbles")
    # variation graph: spine + bubbles + 3 copy-paths
    stage(ax, XB, 0.35, WB, 2.35, "VARIATION GRAPH", ORANGE)
    sy = 1.5; sx0 = XB + 0.5; sx1 = XB + 3.6
    ax.add_patch(Rectangle((sx0, sy - 0.05), sx1 - sx0, 0.1, fc=SPINE, ec="none"))
    bx = [XB + 1.2, XB + 2.0, XB + 2.8]
    for x in bx:
        ax.add_patch(Circle((x, sy + 0.5), 0.13, fc=RED, ec="k", lw=0.6))
        ax.add_patch(Circle((x, sy), 0.13, fc=BLUE, ec="k", lw=0.6))
        ax.add_patch(Circle((x, sy - 0.5), 0.13, fc=GREEN, ec="k", lw=0.6))
    for off, c in [(0.5, RED), (0.0, BLUE), (-0.5, GREEN)]:
        ax.plot([sx0] + bx + [sx1], [sy] + [sy + off] * 3 + [sy], color=c, lw=2, alpha=.85)
    ax.text(XB + WB / 2, 2.15, "PSV bubbles = where copies differ", ha="center", fontsize=8, color="#555")
    ax.text(XB + WB / 2, 0.5, "a copy = a path", ha="center", fontsize=8.2, style="italic", color="#555")
    arrow(ax, XB + WB, XC, 1.5, "significance gate")
    stage(ax, XC, 0.35, WC, 2.35, "ASSIGNED  (certificate)", GREEN)
    rows = [("read → c0", "414 decisive PSVs, p≈0", GREEN),
            ("read → c1", "347 decisive PSVs, p≈0", GREEN),
            ("read → c2", "328 decisive PSVs, p≈0", GREEN),
            ("read → TIED", "0 distinguishing bubbles → abstain", GREY)]
    for i, (a, b, c) in enumerate(rows):
        ax.text(XC + 0.35, 2.15 - i * 0.5, a, fontsize=9.5, weight="bold", color=c)
        ax.text(XC + 1.7, 2.15 - i * 0.5, b, fontsize=8.0, color="#444")
    ax.text(XC + WC / 2, 0.45, "2654 assigned · 16 tied · never 1/k", ha="center", fontsize=8.4, weight="bold", color=DARK)

    # ================= O3 — reads by allele -> junction -> allele-specific splicing =================
    ax = axes[2]
    j = json.load(open(f"{SC}/abcc4_o3.json")) if os.path.exists(f"{SC}/abcc4_o3.json") else \
        dict(nT=173, nC=117, T_use=0, T_skip=173, C_use=45, C_skip=72)
    psiT = j["T_use"] / j["nT"]; psiC = j["C_use"] / j["nC"]
    ax.text(0.1, 3.05, "O3", fontsize=20, weight="bold", color=NAVY)
    ax.text(0.85, 3.06, "allele-specific junctions (ABCC4)", fontsize=11.5, color=DARK, va="center")
    stage(ax, XA, 0.35, WA, 2.35, "READS split by allele", NAVY)
    ax.add_patch(Circle((XA + 0.8, 2.0), 0.16, fc="#d62728", ec="k", lw=0.6)); ax.text(XA + 1.15, 2.0, f"allele T  ({j['nT']} reads)", va="center", fontsize=9.2, weight="bold")
    ax.add_patch(Circle((XA + 0.8, 1.1), 0.16, fc="#1f77b4", ec="k", lw=0.6)); ax.text(XA + 1.15, 1.1, f"allele C  ({j['nC']} reads)", va="center", fontsize=9.2, weight="bold")
    ax.text(XA + WA / 2, 0.55, "het SNP at NC_073238.2:109,943,482", ha="center", fontsize=7.8, color="#555")
    arrow(ax, XA + WA, XB, 1.5, "same molecule")
    stage(ax, XB, 0.35, WB, 2.35, "A JUNCTION", ORANGE)
    ax.add_patch(Rectangle((XB + 0.5, 1.42), 1.0, 0.16, fc="#8fbf8f", ec="k", lw=0.6))
    ax.add_patch(Rectangle((XB + 2.6, 1.42), 1.0, 0.16, fc="#8fbf8f", ec="k", lw=0.6))
    ax.add_patch(FancyArrowPatch((XB + 1.5, 1.5), (XB + 2.6, 1.5), connectionstyle="arc3,rad=-0.7",
                                 arrowstyle="-", lw=2, color=ORANGE, ls=(0, (4, 2))))
    ax.text(XB + WB / 2, 2.05, "donor 109,947,543", ha="center", fontsize=8.4, color=ORANGE, weight="bold")
    ax.text(XB + WB / 2, 0.85, "use it, or splice it out?", ha="center", fontsize=8.2, style="italic", color="#555")
    arrow(ax, XB + WB, XC, 1.5, "per molecule")
    stage(ax, XC, 0.35, WC, 2.35, "ALLELE-SPECIFIC", GREEN)
    ax.text(XC + 0.35, 2.15, f"T allele:  {j['T_use']}/{j['nT']} use it   (PSI {psiT:.2f})", fontsize=10, weight="bold", color="#d62728")
    ax.text(XC + 0.35, 1.65, f"C allele:  {j['C_use']}/{j['nC']} use it   (PSI {psiC:.2f})", fontsize=10, weight="bold", color="#1f77b4")
    ax.text(XC + 0.35, 1.1, f"ΔPSI = {psiC - psiT:.2f}     q = 3×10⁻¹⁸", fontsize=10.5, weight="bold", color=DARK)
    ax.text(XC + WC / 2, 0.55, "allele → junction on single reads, no phasing", ha="center", fontsize=8.4, style="italic", color="#555")

    # ================= O4 — reads -> variation graph -> copy NOT in the genome =================
    ax = axes[3]
    ax.text(0.1, 3.05, "O4", fontsize=20, weight="bold", color=NAVY)
    ax.text(0.85, 3.06, "a copy NOT in the reference assembly", fontsize=11.5, color=DARK, va="center")
    stage(ax, XA, 0.35, WA, 2.35, "READS map confidently, yet split", NAVY)
    pileup(ax, XA + 2.1, 0.8, 7, "#666", w=1.4)
    ax.text(XA + 2.1, 0.55, "all MAPQ 60 · depth > 1 copy", ha="center", fontsize=8.2, color="#555")
    ax.text(XA + 2.1, 2.25, "siblings absent → reads pile on one locus", ha="center", fontsize=7.9, color="#555", style="italic")
    arrow(ax, XA + WA, XB, 1.5, "discover PSVs")
    stage(ax, XB, 0.35, WB, 2.35, "VARIATION GRAPH", ORANGE)
    sy = 1.5; sx0 = XB + 0.5; sx1 = XB + 3.6
    ax.add_patch(Rectangle((sx0, sy - 0.05), sx1 - sx0, 0.1, fc=SPINE, ec="none"))
    for x in [XB + 1.3, XB + 2.1, XB + 2.9]:
        ax.add_patch(Circle((x, sy + 0.45), 0.13, fc=RED, ec="k", lw=0.6))
        ax.add_patch(Circle((x, sy - 0.45), 0.13, fc=BLUE, ec="k", lw=0.6))
    ax.plot([sx0, XB + 1.3, XB + 2.1, XB + 2.9, sx1], [sy, sy + .45, sy + .45, sy + .45, sy], color=RED, lw=2)
    ax.plot([sx0, XB + 1.3, XB + 2.1, XB + 2.9, sx1], [sy, sy - .45, sy - .45, sy - .45, sy], color=BLUE, lw=2)
    ax.text(XB + WB / 2, 2.15, "co-segregating PSVs → ≥ 2 copies", ha="center", fontsize=8, color="#555")
    ax.text(XB + WB / 2, 0.5, "at ONE reference locus", ha="center", fontsize=8.2, style="italic", color="#555")
    arrow(ax, XB + WB, XC, 1.5, "count vs reference")
    stage(ax, XC, 0.35, WC, 2.35, "REFERENCE-ABSENT COPY", RED)
    ax.text(XC + 0.35, 2.22, "gorilla LOC115932956 (NC_073236.2):", fontsize=9.6, weight="bold", color=DARK)
    ax.text(XC + 0.35, 1.90, "95 reads, ALL MAPQ 60 → 3 co-segregating", fontsize=9.2, color="#333")
    ax.text(XC + 0.35, 1.63, "copy-haplotypes (16 / 14 / 6 reads)", fontsize=9.2, color="#333")
    ax.text(XC + 0.35, 1.28, "paralog divergence, NOT RNA editing", fontsize=9.0, color="#555", style="italic")
    ax.text(XC + 0.35, 0.94, "→ 2 copies collapsed / absent from", fontsize=9.6, weight="bold", color=RED)
    ax.text(XC + 0.35, 0.70, "   the reference assembly", fontsize=9.6, weight="bold", color=RED)
    ax.text(XC + WC / 2, 0.47, "high MAPQ + depth = collapse (SDA), not ambiguity", ha="center", fontsize=7.6, style="italic", color="#777")

    fig.suptitle("From reads to results — the four objectives on real gorilla data",
                 fontsize=16, weight="bold", color=NAVY, y=0.997)
    fig.tight_layout(rect=[0, 0, 1, 0.985], h_pad=1.6)
    p = os.path.join(OUT, "objectives_flow.png")
    fig.savefig(p, bbox_inches="tight", facecolor="white", dpi=150)
    plt.close()
    print("wrote", p, "| O4:", o4)


if __name__ == "__main__":
    build()
