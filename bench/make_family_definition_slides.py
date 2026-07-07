#!/usr/bin/env python3
"""
make_family_definition_slides.py — illustrate what a multi-copy gene family is,
using a real example from the refined catalog (Family 69, RABL2):
  1. definition slide
  2. birdseye quasi-clique view of one family
  3. zoom on one copy showing isoforms

Run: /home/juanfra/miniforge3/bin/python bench/make_family_definition_slides.py
Output: bench/family_definition_slides.pptx (+ 3 PNGs)
"""
import os
import csv
import math
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, Rectangle
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


def load_family(fid="69"):
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
    fig, ax = plt.subplots(figsize=(12.0, 6.2))
    ax.set_xlim(0, 12); ax.set_ylim(0, 6.2); ax.axis("off")

    # chromosome strip
    ax.add_patch(Rectangle((0.8, 4.6), 10.4, 0.25, fc="#dfe7f2", ec=CN, lw=1.5, zorder=1))
    ax.text(6.0, 5.15, "chromosome", ha="center", fontsize=10, color=CN, style="italic")

    # copy A locus
    ax.add_patch(FancyBboxPatch((1.5, 4.25), 2.2, 0.95, boxstyle="round,pad=0.05,rounding_size=0.1",
                 fc=CT, ec=CN, lw=2.0, zorder=2))
    ax.text(2.6, 4.95, "copy A", ha="center", va="center", fontsize=13, color="white", fontweight="bold")
    ax.text(2.6, 4.55, "locus 1", ha="center", va="center", fontsize=9, color="white")

    # copy B locus
    ax.add_patch(FancyBboxPatch((7.0, 4.25), 2.2, 0.95, boxstyle="round,pad=0.05,rounding_size=0.1",
                 fc=CO, ec=CN, lw=2.0, zorder=2))
    ax.text(8.1, 4.95, "copy B", ha="center", va="center", fontsize=13, color="white", fontweight="bold")
    ax.text(8.1, 4.55, "locus 2", ha="center", va="center", fontsize=9, color="white")

    # homology edge
    ax.annotate("", xy=(7.0, 4.72), xytext=(3.7, 4.72),
                arrowprops=dict(arrowstyle="<->", color=CN, lw=2.5))
    ax.text(5.35, 5.0, "homology edge", ha="center", fontsize=11, color=CN, fontweight="bold")

    # definition box
    ax.add_patch(FancyBboxPatch((1.2, 0.75), 9.6, 3.2, boxstyle="round,pad=0.1,rounding_size=0.15",
                 fc="#f6f8fb", ec=CN, lw=2.0, zorder=1))

    # title at the top of the box
    ax.text(6.0, 3.75, "Multi-copy gene family", ha="center", fontsize=20, color=CN, fontweight="bold")

    # definition equation below title with clear separation
    ax.text(6.0, 2.55,
            "= two or more distinct genomic loci\n"
            "  that carry homologous copies of the same gene\n"
            "  and are expressed as RNA",
            ha="center", fontsize=14, color=CG, linespacing=1.5)

    # bullet distinctions
    ax.text(2.0, 1.55, "• copies = different loci (often different chromosomes)", ha="left", fontsize=11, color=CG)
    ax.text(2.0, 1.15, "• isoforms = same locus, different splicing", ha="left", fontsize=11, color=CG)

    fig.text(0.5, 0.02,
             "A family is a connected, cohesive subgraph (γ-quasi-clique) of the homology graph; "
             "each node is a distinct locus, each edge is a sequence-homology link.",
             ha="center", fontsize=9, style="italic", color=CG)
    fig.savefig(FIG_DEF, dpi=150, bbox_inches="tight")
    plt.close(fig)


# ----------------------------------------------------------------- figure 2: birdseye quasi-clique
def make_birdseye():
    members = load_family("69")
    meta = load_meta()
    edges = load_edges(members)

    fig, (axL, axR) = plt.subplots(1, 2, figsize=(13.3, 6.2),
                                   gridspec_kw={"width_ratios": [1.25, 1]})

    # ---- left: genomic view (one strip per chromosome) ----
    axL.set_xlim(0, 10)
    axL.set_ylim(0, 6.2)
    axL.axis("off")

    # sort by chromosome name for stable layout
    members_sorted = sorted(members, key=lambda m: (m["chrom"], int(m["start"])))
    colors = [CT, CO, CP, "#2874A6", "#117864"]
    y0 = 5.1
    dy = 1.05

    for i, m in enumerate(members_sorted):
        y = y0 - i * dy
        chrom = m["chrom"]
        s, e = int(m["start"]), int(m["end"])
        gene = m["member_gene"] if m["member_gene"] != "NA" else f"copy {i+1}"
        tid = m["member_dn"]
        n_exon = meta[tid]["n_exon"]

        # chromosome bar
        axL.add_patch(Rectangle((2.2, y - 0.08), 7.0, 0.16, fc="#dfe7f2", ec=CN, lw=1.2))
        # locus rectangle (width proportional to log10 length, capped)
        length = e - s
        width = min(2.8, max(0.6, math.log10(length) * 0.55))
        color = colors[i % len(colors)]
        axL.add_patch(FancyBboxPatch((3.0, y - 0.18), width, 0.36,
                     boxstyle="round,pad=0.02,rounding_size=0.05",
                     fc=color, ec=CN, lw=1.5, zorder=2))
        # label inside the locus rectangle (coordinates omitted to avoid clutter)
        label = f"{gene}"
        axL.text(3.0 + width / 2, y, label, ha="center", va="center",
                 fontsize=8, color="white", fontweight="bold")
        # chromosome label on the left
        axL.text(2.05, y, chrom, ha="right", va="center", fontsize=9, color=CN, fontweight="bold")

    axL.set_title("Family 69 — RABL2  (real refined family)",
                  fontsize=14, fontweight="bold", color=CN, pad=12)
    axL.text(5.0, 0.35,
             "5 distinct loci on 5 different chromosomes  ·  RABL2A and RABL2B are highlighted paralogs",
             ha="center", fontsize=10, color=CG, style="italic")

    # ---- right: homology graph ----
    axR.set_xlim(0, 6); axR.set_ylim(0, 6); axR.axis("off")
    axR.set_title("the homology graph (perfect clique)", fontsize=14, fontweight="bold", color=CN)

    # pentagon layout, large enough to avoid label collisions
    ids = [m["member_dn"] for m in members_sorted]
    labels = {}
    for i, m in enumerate(members_sorted):
        gene = m["member_gene"] if m["member_gene"] != "NA" else f"copy {i+1}"
        chrom = m["chrom"].replace("NC_", "").replace(".1", "").replace(".2", "")
        labels[m["member_dn"]] = f"{gene}\nchr {chrom}"

    center = (3.0, 3.0)
    radius = 2.1
    pos = {}
    for idx, tid in enumerate(ids):
        angle = math.pi / 2 - idx * (2 * math.pi / len(ids))
        pos[tid] = (center[0] + radius * math.cos(angle),
                    center[1] + radius * math.sin(angle))

    node_colors = {tid: colors[i % len(colors)] for i, tid in enumerate(ids)}

    # edges (complete graph)
    edge_seen = set()
    for a, b, cr in edges:
        if (b, a) in edge_seen:
            continue
        edge_seen.add((a, b))
        axR.plot([pos[a][0], pos[b][0]], [pos[a][1], pos[b][1]],
                 color=CN, lw=1.8, alpha=0.75, zorder=1)

    # nodes
    for tid, (x, y) in pos.items():
        fc = node_colors[tid]
        axR.add_patch(FancyBboxPatch((x - 0.62, y - 0.38), 1.24, 0.76,
                     boxstyle="round,pad=0.04,rounding_size=0.08",
                     fc=fc, ec=CN, lw=2.0, zorder=3))
        axR.text(x, y, labels[tid], ha="center", va="center",
                 fontsize=7.5, color="white", fontweight="bold", zorder=4, linespacing=1.05)

    # legend / note
    axR.text(3.0, 0.55,
             "all 10 possible edges present  →  perfect clique\n"
             "every locus is homologous to every other locus",
             ha="center", fontsize=9, color=CG, style="italic", linespacing=1.3)

    fig.savefig(FIG_BIRD, dpi=150, bbox_inches="tight")
    plt.close(fig)


# ----------------------------------------------------------------- figure 3: isoform zoom
def make_isoforms():
    meta = load_meta()
    skel = load_skeletons()

    fig, ax = plt.subplots(figsize=(12.0, 5.5))
    ax.set_xlim(48817500, 48832500)
    ax.set_ylim(0, 5.0)
    ax.axis("off")

    # chromosome bar
    ax.add_patch(Rectangle((48818000, 0.35), 14000, 0.12, fc="#dfe7f2", ec=CN, lw=1.5))
    ax.text(48825000, 0.15, "NC_086018.1  (RABL2B locus)", ha="center", fontsize=9, color=CN)

    isoforms = [
        ("DN_NC_086018.1_48818439_9", "isoform a  (9 exons, 4 reads)", 3.8, CT),
        ("DN_NC_086018.1_48818439_10", "isoform b  (10 exons, 3 reads)", 2.2, CO),
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
        ax.text(48817300, y, label, ha="right", va="center", fontsize=11, color=col, fontweight="bold")

    # coordinate ticks
    for x in range(48819000, 48832001, 5000):
        ax.plot([x, x], [0.35, 0.26], color=CG, lw=1.0)
        ax.text(x, 0.08, f"{x/1e6:.3f}M", ha="center", fontsize=8, color=CG)

    ax.set_title("Zoom: RABL2B produces multiple isoforms from one locus",
                 fontsize=14, fontweight="bold", color=CN, pad=12)
    ax.text(48825000, 4.6,
            "Same locus (start 48,818,439)  ·  different splicing / 3' ends  ·  shared first exon",
            ha="center", fontsize=10, color=CG, style="italic")

    # legend
    ax.add_patch(Rectangle((48830500, 4.15), 120, 0.25, fc=CT, ec=CN))
    ax.text(48830700, 4.27, "exon", ha="left", va="center", fontsize=9, color=CG)
    ax.plot([48830500, 48830620], [3.85, 3.85], color=CO, lw=2.0)
    ax.text(48830700, 3.85, "intron", ha="left", va="center", fontsize=9, color=CG)

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
    scale = min(Inches(12.4) / w, Inches(5.5) / h)
    iw, ih = int(w * scale), int(h * scale)
    s.shapes.add_picture(img, int((Inches(13.33) - iw) / 2), Inches(1.05), width=iw, height=ih)

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
                     "Copies are different loci (often on different chromosomes); isoforms are different transcripts from the same locus.")

    add_figure_slide(prs,
                     "Birdseye view: a real family as a clique",
                     FIG_BIRD,
                     "Family 69 (RABL2) has five copies on five gorilla chromosomes. The homology graph is a perfect clique: "
                     "every locus is linked to every other locus, so the group is one cohesive family.")

    add_figure_slide(prs,
                     "Zoom: one copy, multiple isoforms",
                     FIG_ISO,
                     "RABL2B at 48.82 Mb on NC_086018.1 produces two representative isoforms. They share the first exon but differ "
                     "in splicing and 3' extent. Isoforms are collapsed to one locus before families are called.")

    prs.save(OUT)
    print(f"[+] wrote {OUT} ({len(prs.slides._sldIdLst)} slides)")


if __name__ == "__main__":
    build()
