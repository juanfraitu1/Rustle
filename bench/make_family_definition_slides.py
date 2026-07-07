#!/usr/bin/env python3
"""
make_family_definition_slides.py — illustrate what a multi-copy gene family is,
using a real example from the refined catalog:
  1. definition slide
  2. birdseye quasi-clique view of one family
  3. zoom on one copy showing isoforms

Run: /home/juanfra/miniforge3/bin/python bench/make_family_definition_slides.py
Output: bench/family_definition_slides.pptx (+ 3 PNGs)
"""
import os
import csv
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch, Rectangle
from pptx import Presentation
from pptx.util import Inches, Pt
from pptx.dml.color import RGBColor
from pptx.enum.text import PP_ALIGN

HERE = os.path.dirname(os.path.abspath(__file__))
REFINE_TSV = os.path.join(HERE, "family_rna_refine.tsv")
META_TSV = os.path.join("/home/juanfra/winloci_scratch", "denovo_transcripts.meta.tsv")
SKEL_TSV = os.path.join("/home/juanfra/winloci_scratch", "denovo_skeletons.tsv")
EDGES_TSV = os.path.join(HERE, "denovo_family_edges.tsv")

FIG_DEF = os.path.join(HERE, "famdef_1_definition.png")
FIG_BIRD = os.path.join(HERE, "famdef_2_birdseye.png")
FIG_ISO = os.path.join(HERE, "famdef_3_isoforms.png")
OUT = os.path.join(HERE, "family_definition_slides.pptx")

NAVY = RGBColor(0x1F, 0x2D, 0x5A)
TEAL = RGBColor(0x12, 0x7A, 0x6E)
GREY = RGBColor(0x44, 0x44, 0x44)
ORANGE = RGBColor(0xD9, 0x6B, 0x27)
PURPLE = RGBColor(0x8E, 0x44, 0xAD)

CN, CT, CG, CO, CP = "#1F2D5A", "#127A6E", "#444444", "#D96B27", "#8E44AD"


def load_meta():
    d = {}
    with open(META_TSV) as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            d[r["id"]] = dict(
                chrom=r["chrom"],
                start=int(r["start"]),
                end=int(r["end"]),
                strand=r["strand"],
                n_exon=int(r["n_exon"]),
                n_reads=int(r["n_reads"]),
            )
    return d


def load_skeletons():
    d = {}
    with open(SKEL_TSV) as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            c, s, e = r["chrom"], int(r["start"]), int(r["end"])
            introns = []
            if r["introns"].strip():
                for tok in r["introns"].split(";"):
                    a, b = tok.split("-")
                    introns.append((int(a), int(b)))
            d[(c, s, e)] = introns
    return d


def exon_intervals(meta, skel, tid):
    m = meta[tid]
    c, s, e = m["chrom"], m["start"], m["end"]
    introns = skel.get((c, s, e), [])
    exons = []
    cur = s
    for d, a in introns:
        exons.append((cur, d))
        cur = a
    exons.append((cur, e + 1))
    return exons


def load_family(fid="45"):
    members = []
    with open(REFINE_TSV) as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            if r["family_id"] == fid:
                members.append(r)
    return members


def load_edges(members):
    ids = {m["member_dn"] for m in members}
    edges = []
    with open(EDGES_TSV) as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            if r["a"] in ids and r["b"] in ids:
                edges.append((r["a"], r["b"], float(r["core_recip"])))
    return edges


# ----------------------------------------------------------------- figure 1: definition
def make_definition():
    fig, ax = plt.subplots(figsize=(12.0, 6.0))
    ax.set_xlim(0, 12); ax.set_ylim(0, 6); ax.axis("off")

    # chromosome strip
    ax.add_patch(Rectangle((0.8, 4.5), 10.4, 0.25, fc="#dfe7f2", ec=CN, lw=1.5, zorder=1))
    ax.text(6.0, 5.05, "chromosome", ha="center", fontsize=10, color=CN, style="italic")

    # copy A locus
    ax.add_patch(FancyBboxPatch((1.5, 4.15), 2.2, 0.95, boxstyle="round,pad=0.05,rounding_size=0.1",
                 fc=CT, ec=CN, lw=2.0, zorder=2))
    ax.text(2.6, 4.85, "copy A", ha="center", va="center", fontsize=13, color="white", fontweight="bold")
    ax.text(2.6, 4.45, "locus 1", ha="center", va="center", fontsize=9, color="white")

    # copy B locus
    ax.add_patch(FancyBboxPatch((7.0, 4.15), 2.2, 0.95, boxstyle="round,pad=0.05,rounding_size=0.1",
                 fc=CO, ec=CN, lw=2.0, zorder=2))
    ax.text(8.1, 4.85, "copy B", ha="center", va="center", fontsize=13, color="white", fontweight="bold")
    ax.text(8.1, 4.45, "locus 2", ha="center", va="center", fontsize=9, color="white")

    # homology edge
    ax.annotate("", xy=(7.0, 4.62), xytext=(3.7, 4.62),
                arrowprops=dict(arrowstyle="<->", color=CN, lw=2.5))
    ax.text(5.35, 4.9, "homology edge", ha="center", fontsize=11, color=CN, fontweight="bold")

    # definition box
    ax.add_patch(FancyBboxPatch((1.2, 1.0), 9.6, 2.6, boxstyle="round,pad=0.1,rounding_size=0.15",
                 fc="#f6f8fb", ec=CN, lw=2.0, zorder=1))
    ax.text(6.0, 3.15, "Multi-copy gene family", ha="center", fontsize=18, color=CN, fontweight="bold")
    ax.text(6.0, 2.4,
            "= two or more distinct genomic loci\n"
            "  that carry homologous copies of the same gene\n"
            "  and are expressed as RNA",
            ha="center", fontsize=14, color=CG, linespacing=1.4)

    # bullet distinctions
    ax.text(2.0, 1.6, "• copies = different loci", ha="left", fontsize=11, color=CG)
    ax.text(2.0, 1.2, "• isoforms = same locus, different splicing", ha="left", fontsize=11, color=CG)

    fig.text(0.5, 0.02,
             "A family is a connected, cohesive subgraph (γ-quasi-clique) of the homology graph; "
             "each node is a distinct locus, each edge is a sequence-homology link.",
             ha="center", fontsize=9, style="italic", color=CG)
    fig.savefig(FIG_DEF, dpi=150, bbox_inches="tight")
    plt.close(fig)


# ----------------------------------------------------------------- figure 2: birdseye quasi-clique
def make_birdseye():
    members = load_family("45")
    meta = load_meta()
    edges = load_edges(members)

    fig, (axL, axR) = plt.subplots(1, 2, figsize=(13.3, 6.0),
                                   gridspec_kw={"width_ratios": [1.3, 1]})

    # ---- left: genomic view ----
    axL.set_xlim(8220000, 8250000)
    axL.set_ylim(0, 4)
    axL.axis("off")

    # chromosome bar
    ychrom = 2.0
    axL.add_patch(Rectangle((8220000, ychrom - 0.08), 30000, 0.16, fc="#dfe7f2", ec=CN, lw=1.5))
    axL.text(8235000, ychrom + 0.45, "NC_073224.2  (gorilla)", ha="center", fontsize=10, color=CN)

    # copy 1 is one locus; copy 2 has two isoforms that map to the same locus
    # draw copy 1
    copy1 = [m for m in members if m["member_dn"] == "DN_NC_073224.2_8222521_2"][0]
    s1, e1 = int(copy1["start"]), int(copy1["end"])
    axL.add_patch(Rectangle((s1, ychrom + 0.55), e1 - s1, 0.32, fc=CT, ec=CN, lw=1.5))
    axL.text((s1 + e1) / 2, ychrom + 0.71, "copy 1", ha="center", va="center",
             fontsize=8, color="white", fontweight="bold")

    # draw copy 2 locus as union, with inner shading for each isoform
    copy2_ids = ["DN_NC_073224.2_8222533_2", "DN_NC_073224.2_8222533_5"]
    copy2_members = [m for m in members if m["member_dn"] in copy2_ids]
    s2 = min(int(m["start"]) for m in copy2_members)
    e2 = max(int(m["end"]) for m in copy2_members)
    axL.add_patch(Rectangle((s2, ychrom - 0.87), e2 - s2, 0.32, fc=CO, ec=CN, lw=1.5))
    axL.text((s2 + e2) / 2, ychrom - 0.71, "copy 2  (2 isoforms)", ha="center", va="center",
             fontsize=8, color="white", fontweight="bold")

    # axis ticks
    for x in range(8221000, 8249001, 5000):
        axL.plot([x, x], [ychrom - 0.18, ychrom - 0.28], color=CG, lw=1.0)
        axL.text(x, ychrom - 0.45, f"{x/1e6:.3f}M", ha="center", fontsize=8, color=CG)

    axL.set_title("Family 45 — LOC101124778  (real refined family)",
                  fontsize=14, fontweight="bold", color=CN, pad=10)
    axL.text(8235000, 0.4,
             "2 distinct loci  ·  copy 2 has 2 isoforms  ·  3 member transcripts",
             ha="center", fontsize=10, color=CG, style="italic")

    # ---- right: homology graph ----
    axR.set_xlim(0, 6); axR.set_ylim(0, 6); axR.axis("off")
    axR.set_title("the homology graph (γ-quasi-clique)", fontsize=14, fontweight="bold", color=CN)

    pos = {
        "DN_NC_073224.2_8222521_2": (1.5, 4.5),
        "DN_NC_073224.2_8222533_2": (4.5, 5.5),
        "DN_NC_073224.2_8222533_5": (4.5, 3.5),
    }
    labels = {
        "DN_NC_073224.2_8222521_2": "copy 1\n2 exons",
        "DN_NC_073224.2_8222533_2": "copy 2\nisoform a",
        "DN_NC_073224.2_8222533_5": "copy 2\nisoform b",
    }
    colors = {"DN_NC_073224.2_8222521_2": CT,
              "DN_NC_073224.2_8222533_2": CO,
              "DN_NC_073224.2_8222533_5": CO}

    # missing edges shown as faint grey
    possible = [("DN_NC_073224.2_8222521_2", "DN_NC_073224.2_8222533_5"),
                ("DN_NC_073224.2_8222533_2", "DN_NC_073224.2_8222533_5")]
    for a, b in possible:
        if not any((a == e[0] and b == e[1]) or (a == e[1] and b == e[0]) for e in edges):
            axR.plot([pos[a][0], pos[b][0]], [pos[a][1], pos[b][1]],
                     color="#cccccc", lw=1.5, ls="--", zorder=1)

    # real edges
    edge_seen = set()
    for a, b, cr in edges:
        if (b, a) in edge_seen:
            continue
        edge_seen.add((a, b))
        axR.annotate("", xy=pos[b], xytext=pos[a],
                     arrowprops=dict(arrowstyle="-", color=CN, lw=2.5))
        mx, my = (pos[a][0] + pos[b][0]) / 2, (pos[a][1] + pos[b][1]) / 2
        axR.text(mx + 0.15, my, f"{cr:.2f}", fontsize=9, color=CN, fontweight="bold")

    # nodes
    for tid, (x, y) in pos.items():
        fc = CT if colors[tid] == CT else CO
        axR.add_patch(FancyBboxPatch((x - 0.65, y - 0.45), 1.3, 0.9,
                     boxstyle="round,pad=0.05,rounding_size=0.1",
                     fc=fc, ec=CN, lw=2.0, zorder=3))
        axR.text(x, y, labels[tid], ha="center", va="center",
                 fontsize=9, color="white", fontweight="bold", zorder=4)

    # legend / note
    axR.text(3.0, 0.6,
             "solid = real homology edge  ·  dashed = no edge (not required in a quasi-clique)\n"
             "density = 1/3 edges  ≥  γ=0.20  →  valid family",
             ha="center", fontsize=9, color=CG, style="italic", linespacing=1.3)

    fig.savefig(FIG_BIRD, dpi=150, bbox_inches="tight")
    plt.close(fig)


# ----------------------------------------------------------------- figure 3: isoform zoom
def make_isoforms():
    meta = load_meta()
    skel = load_skeletons()

    fig, ax = plt.subplots(figsize=(12.0, 5.5))
    ax.set_xlim(8221500, 8247500)
    ax.set_ylim(0, 4.5)
    ax.axis("off")

    # chromosome bar
    ax.add_patch(Rectangle((8222000, 0.3), 24500, 0.12, fc="#dfe7f2", ec=CN, lw=1.5))
    ax.text(8234250, 0.15, "NC_073224.2", ha="center", fontsize=9, color=CN)

    isoforms = [
        ("DN_NC_073224.2_8222533_2", "isoform a  (2 exons, 3 reads)", 3.6, CT),
        ("DN_NC_073224.2_8222533_5", "isoform b  (5 exons, 26 reads)", 2.0, CO),
    ]

    for tid, label, y, col in isoforms:
        m = meta[tid]
        exons = exon_intervals(meta, skel, tid)
        # transcript backbone
        ax.plot([m["start"], m["end"]], [y, y], color=col, lw=2.0, zorder=1)
        # exons
        for s, e in exons:
            ax.add_patch(Rectangle((s, y - 0.18), e - s, 0.36, fc=col, ec=CN, lw=1.0, zorder=2))
        # intron tick marks
        for d, a in skel[(m["chrom"], m["start"], m["end"])]:
            ax.plot([d, a], [y, y], color=col, lw=2.0, zorder=1)
        ax.text(8221800, y, label, ha="right", va="center", fontsize=11, color=col, fontweight="bold")

    # coordinate ticks
    for x in range(8223000, 8247001, 5000):
        ax.plot([x, x], [0.3, 0.22], color=CG, lw=1.0)
        ax.text(x, 0.08, f"{x/1e6:.3f}M", ha="center", fontsize=8, color=CG)

    ax.set_title("Zoom: copy 2 of Family 45 has two different isoforms",
                 fontsize=14, fontweight="bold", color=CN, pad=10)
    ax.text(8234250, 4.2,
            "Same locus (start 8,222,533)  ·  different splicing / 3' ends  ·  shared first exon",
            ha="center", fontsize=10, color=CG, style="italic")

    # legend
    ax.add_patch(Rectangle((8238000, 3.85), 120, 0.25, fc=CT, ec=CN))
    ax.text(8238200, 3.97, "exon", ha="left", va="center", fontsize=9, color=CG)
    ax.plot([8238000, 8238120], [3.55, 3.55], color=CO, lw=2.0)
    ax.text(8238200, 3.55, "intron", ha="left", va="center", fontsize=9, color=CG)

    fig.savefig(FIG_ISO, dpi=150, bbox_inches="tight")
    plt.close(fig)


# ----------------------------------------------------------------- pptx assembly
def add_figure_slide(prs, title, img, caption):
    from PIL import Image
    s = prs.slides.add_slide(prs.slide_layouts[6])
    tb = s.shapes.add_textbox(Inches(0.5), Inches(0.22), Inches(12.3), Inches(0.8))
    tp = tb.text_frame.paragraphs[0]
    tp.text = title
    tp.font.size = Pt(24)
    tp.font.bold = True
    tp.font.color.rgb = NAVY

    w, h = Image.open(img).size
    scale = min(Inches(12.4) / w, Inches(5.4) / h)
    iw, ih = int(w * scale), int(h * scale)
    s.shapes.add_picture(img, int((Inches(13.33) - iw) / 2), Inches(1.15), width=iw, height=ih)

    cb = s.shapes.add_textbox(Inches(0.6), Inches(6.75), Inches(12.1), Inches(0.6))
    cp = cb.text_frame.paragraphs[0]
    cp.text = caption
    cp.font.size = Pt(12)
    cp.font.color.rgb = GREY
    cp.alignment = PP_ALIGN.CENTER
    cb.text_frame.word_wrap = True
    return s


def build():
    make_definition()
    make_birdseye()
    make_isoforms()

    prs = Presentation()
    prs.slide_width = Inches(13.33)
    prs.slide_height = Inches(7.5)

    add_figure_slide(prs,
                     "What is a multi-copy gene family?",
                     FIG_DEF,
                     "A family is ≥2 distinct genomic loci carrying homologous, expressed copies of one gene. "
                     "Copies are different loci; isoforms are different transcripts from the same locus.")

    add_figure_slide(prs,
                     "Birdseye view: a real family as a quasi-clique",
                     FIG_BIRD,
                     "Family 45 (LOC101124778) has two copies on gorilla NC_073224.2. The homology graph is a "
                     "γ-quasi-clique: not every pair needs an edge, but the group is cohesive enough to be one family.")

    add_figure_slide(prs,
                     "Zoom: one copy, multiple isoforms",
                     FIG_ISO,
                     "Copy 2 at 8.22 Mb produces two isoforms. They share the first exon but differ in splicing "
                     "and 3' extent. Isoforms are collapsed to one locus before families are called.")

    prs.save(OUT)
    print(f"[+] wrote {OUT} ({len(prs.slides._sldIdLst)} slides)")


if __name__ == "__main__":
    build()
