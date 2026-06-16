#!/usr/bin/env python3
"""
Build the 7-slide read-coherence deck.

Renders slide1.png .. slide7.png (16:9, figsize=(13.33,7.5), dpi=150, white bg)
and read_coherence_deck.pdf (pages 1..7) via matplotlib's PdfPages.

Run with: python3 build_slides.py
"""
import os
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.image as mpimg
from matplotlib.backends.backend_pdf import PdfPages

# ----------------------------------------------------------------------------
# Paths
# ----------------------------------------------------------------------------
HERE = os.path.dirname(os.path.abspath(__file__))
SRC = os.path.dirname(HERE)  # the readcoherence_fig directory (read-only inputs)
OUT = HERE

FIG_DIR = SRC

# ----------------------------------------------------------------------------
# Global style
# ----------------------------------------------------------------------------
NAVY = "#22313f"
ORANGE = "#e8590c"
SLATE = "#8898aa"
WHITE = "#ffffff"

FIGSIZE = (13.33, 7.5)
DPI = 150

plt.rcParams["font.family"] = "DejaVu Sans"  # clean sans
plt.rcParams["pdf.fonttype"] = 42

TOTAL = 7


def new_slide():
    """Create a blank 16:9 white slide; return (fig, ax) spanning [0,1]x[0,1]."""
    fig = plt.figure(figsize=FIGSIZE, dpi=DPI)
    fig.patch.set_facecolor(WHITE)
    ax = fig.add_axes([0, 0, 1, 1])
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.axis("off")
    ax.set_facecolor(WHITE)
    return fig, ax


def add_footer(ax, n):
    """Thin footer: 'rustle . read-coherence' left, 'n / 7' right, small slate."""
    ax.text(0.045, 0.035, "rustle · read-coherence",
            color=SLATE, fontsize=9.5, ha="left", va="center")
    ax.text(0.955, 0.035, f"{n} / {TOTAL}",
            color=SLATE, fontsize=9.5, ha="right", va="center")
    # thin separator line above footer
    ax.plot([0.045, 0.955], [0.072, 0.072], color=SLATE, lw=0.6, alpha=0.45)


def add_fitted_image(ax, path):
    """
    Place an image filling the content area while preserving aspect ratio
    (centered, no distortion, no cropping). Content area sits above the footer.
    """
    img = mpimg.imread(path)
    ih, iw = img.shape[0], img.shape[1]
    img_aspect = iw / ih  # width / height

    # Content box in figure-fraction coords (x0, y0, x1, y1).
    # Leave generous margins; reserve bottom strip for footer.
    cx0, cy0, cx1, cy1 = 0.035, 0.10, 0.965, 0.965
    box_w_frac = cx1 - cx0
    box_h_frac = cy1 - cy0

    # Convert box to physical inches to compare aspect ratios correctly.
    fig_w_in, fig_h_in = FIGSIZE
    box_w_in = box_w_frac * fig_w_in
    box_h_in = box_h_frac * fig_h_in
    box_aspect = box_w_in / box_h_in

    if img_aspect > box_aspect:
        # image is relatively wider -> width-constrained
        draw_w_in = box_w_in
        draw_h_in = draw_w_in / img_aspect
    else:
        # image is relatively taller -> height-constrained
        draw_h_in = box_h_in
        draw_w_in = draw_h_in * img_aspect

    draw_w_frac = draw_w_in / fig_w_in
    draw_h_frac = draw_h_in / fig_h_in

    # center inside the content box
    x0 = cx0 + (box_w_frac - draw_w_frac) / 2.0
    y0 = cy0 + (box_h_frac - draw_h_frac) / 2.0

    fig = ax.figure
    iax = fig.add_axes([x0, y0, draw_w_frac, draw_h_frac])
    iax.imshow(img, interpolation="bilinear", aspect="auto")
    iax.axis("off")
    return iax


# ----------------------------------------------------------------------------
# Slide 1 : Title (custom, no figure)
# ----------------------------------------------------------------------------
def slide1():
    fig, ax = new_slide()

    # Big bold navy title
    ax.text(0.5, 0.665, "Read-coherence",
            color=NAVY, fontsize=58, fontweight="bold",
            ha="center", va="center")

    # Medium subtitle
    ax.text(0.5, 0.535, "Per-molecule isoform recovery, additive over flow assembly",
            color=NAVY, fontsize=21, ha="center", va="center")

    # thin accent rule
    ax.plot([0.34, 0.66], [0.46, 0.46], color=ORANGE, lw=2.0, alpha=0.9)

    # One orange thesis line
    ax.text(0.5, 0.385,
            "Recovers real isoforms StringTie's parsimony merges away —",
            color=ORANGE, fontsize=16.5, ha="center", va="center")
    ax.text(0.5, 0.335,
            "as a provable superset of the flow output",
            color=ORANGE, fontsize=16.5, ha="center", va="center")

    # Small footer line (subtitle-level, not the page footer)
    ax.text(0.5, 0.165, "rustle · long-read transcript assembly",
            color=SLATE, fontsize=12.5, ha="center", va="center")

    add_footer(ax, 1)
    return fig


# ----------------------------------------------------------------------------
# Slide 2 : figure (fig2_mechanism) - no slide title
# ----------------------------------------------------------------------------
def slide2():
    fig, ax = new_slide()
    add_fitted_image(ax, os.path.join(FIG_DIR, "fig2_mechanism.png"))
    add_footer(ax, 2)
    return fig


# ----------------------------------------------------------------------------
# Slide 3 : figure (fig1_igv_locus)
# ----------------------------------------------------------------------------
def slide3():
    fig, ax = new_slide()
    add_fitted_image(ax, os.path.join(FIG_DIR, "fig1_igv_locus.png"))
    add_footer(ax, 3)
    return fig


# ----------------------------------------------------------------------------
# Slide 4 : Rigor slide (custom text, two columns)
# ----------------------------------------------------------------------------
def slide4():
    fig, ax = new_slide()

    # Title
    ax.text(0.5, 0.88, "It's not clustering — and where transcripts begin & end",
            color=NAVY, fontsize=23, fontweight="bold",
            ha="center", va="center")

    left_x = 0.075
    right_x = 0.545
    bullet_indent = 0.018

    header_y = 0.745
    first_bullet_y = 0.635
    line_step = 0.105

    left_header = "Not clustering"
    left_bullets = [
        "Identity = exact intron chain (ordered splice junctions)",
        "Same junctions → same transcript; different → different",
        "No similarity metric, no threshold, no merge radius",
        "Same identity as gffcompare FSM; only knob is a\nread-depth floor (evidence, not a radius)",
    ]

    right_header = "Transcript boundaries"
    right_bullets = [
        "Junctions are base-precise\n(confirmed GT-AG / GC-AG / AT-AC)",
        "Ends = outer envelope of the supporting alignments",
        "Truncated fragments fold into the full-length parent",
        "Single-exon reads excluded (no chain → would\nneed span-clustering)",
    ]

    def render_col(x, header, bullets):
        ax.text(x, header_y, header, color=ORANGE,
                fontsize=18.5, fontweight="bold", ha="left", va="center")
        # underline accent under header
        ax.plot([x, x + 0.30], [header_y - 0.045, header_y - 0.045],
                color=ORANGE, lw=1.4, alpha=0.55)
        y = first_bullet_y
        for b in bullets:
            ax.text(x, y, "•", color=ORANGE, fontsize=15,
                    ha="left", va="top")
            ax.text(x + bullet_indent, y, b, color=NAVY, fontsize=13.5,
                    ha="left", va="top", linespacing=1.25)
            y -= line_step

    render_col(left_x, left_header, left_bullets)
    render_col(right_x, right_header, right_bullets)

    # subtle vertical divider between columns
    ax.plot([0.505, 0.505], [0.18, 0.70], color=SLATE, lw=0.7, alpha=0.35)

    add_footer(ax, 4)
    return fig


# ----------------------------------------------------------------------------
# Slide 5 : figure (fig4_gate_collapse)
# ----------------------------------------------------------------------------
def slide5():
    fig, ax = new_slide()
    add_fitted_image(ax, os.path.join(FIG_DIR, "fig4_gate_collapse.png"))
    add_footer(ax, 5)
    return fig


# ----------------------------------------------------------------------------
# Slide 6 : figure (fig3_architecture)
# ----------------------------------------------------------------------------
def slide6():
    fig, ax = new_slide()
    add_fitted_image(ax, os.path.join(FIG_DIR, "fig3_architecture.png"))
    add_footer(ax, 6)
    return fig


# ----------------------------------------------------------------------------
# Slide 7 : figure (fig5_results)
# ----------------------------------------------------------------------------
def slide7():
    fig, ax = new_slide()
    add_fitted_image(ax, os.path.join(FIG_DIR, "fig5_results.png"))
    add_footer(ax, 7)
    return fig


def main():
    builders = [slide1, slide2, slide3, slide4, slide5, slide6, slide7]

    figs = []
    for i, b in enumerate(builders, start=1):
        fig = b()
        png_path = os.path.join(OUT, f"slide{i}.png")
        fig.savefig(png_path, dpi=DPI, facecolor=WHITE)
        print("wrote", png_path)
        figs.append(fig)

    pdf_path = os.path.join(OUT, "read_coherence_deck.pdf")
    with PdfPages(pdf_path) as pdf:
        for fig in figs:
            pdf.savefig(fig, facecolor=WHITE)
    print("wrote", pdf_path)

    for fig in figs:
        plt.close(fig)


if __name__ == "__main__":
    main()
