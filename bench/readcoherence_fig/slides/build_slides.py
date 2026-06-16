#!/usr/bin/env python3
"""
Build the 13-slide read-coherence deck.

Renders slide1.png .. slide13.png (16:9, figsize=(13.33,7.5), dpi=150, white bg)
and read_coherence_deck.pdf (pages 1..13) via matplotlib's PdfPages.

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

TOTAL = 13


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
    """Thin footer: 'rustle . read-coherence' left, 'n / 13' right, small slate."""
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
# Slide 4 : Identity = GROUP BY (custom text, two columns) — airy
# ----------------------------------------------------------------------------
def slide4():
    fig, ax = new_slide()

    # Title (navy, bold)
    ax.text(0.5, 0.885,
            "Identity = GROUP BY the exact junction string (not clustering)",
            color=NAVY, fontsize=22, fontweight="bold",
            ha="center", va="center")

    left_x = 0.075
    right_x = 0.545
    bullet_indent = 0.022

    header_y = 0.735
    first_bullet_y = 0.605
    line_step = 0.125

    left_header = "How a transcript is defined"
    left_bullets = [
        "Read → its ordered list of splice junctions (introns)",
        "Same exact list → same transcript\n(hash key / GROUP BY)",
        "An equivalence relation: one unique,\ndeterministic partition",
        "Truncated fragments → folded into their\ncontaining chain (sub-path)",
    ]

    right_header = "Why it isn't clustering"
    right_bullets = [
        "No distance, no similarity threshold, no merge radius",
        "No objective to optimize → result can't drift\nwith a parameter",
        "Distinct junctions never merge (an exon-skip is not\nits inclusion isoform)",
        "Junctions matched exactly; small tolerance only at\nthe 5'/3' ends",
    ]

    def render_col(x, header, bullets):
        ax.text(x, header_y, header, color=ORANGE,
                fontsize=18, fontweight="bold", ha="left", va="center")
        # underline accent under header
        ax.plot([x, x + 0.345], [header_y - 0.05, header_y - 0.05],
                color=ORANGE, lw=1.4, alpha=0.55)
        y = first_bullet_y
        for b in bullets:
            ax.text(x, y, "•", color=ORANGE, fontsize=15,
                    ha="left", va="top")
            ax.text(x + bullet_indent, y, b, color=NAVY, fontsize=14.5,
                    ha="left", va="top", linespacing=1.3)
            y -= line_step

    render_col(left_x, left_header, left_bullets)
    render_col(right_x, right_header, right_bullets)

    # subtle vertical divider between columns
    ax.plot([0.505, 0.505], [0.16, 0.685], color=SLATE, lw=0.7, alpha=0.35)

    add_footer(ax, 4)
    return fig


# ----------------------------------------------------------------------------
# Slide 5 : figure (fig6_altsplice)
# ----------------------------------------------------------------------------
def slide5():
    fig, ax = new_slide()
    add_fitted_image(ax, os.path.join(FIG_DIR, "fig6_altsplice.png"))
    add_footer(ax, 5)
    return fig


# ----------------------------------------------------------------------------
# Slide 6 : figure (fig4_gate_collapse)
# ----------------------------------------------------------------------------
def slide6():
    fig, ax = new_slide()
    add_fitted_image(ax, os.path.join(FIG_DIR, "fig4_gate_collapse.png"))
    add_footer(ax, 6)
    return fig


# ----------------------------------------------------------------------------
# Slide 7 : figure (fig3_architecture)
# ----------------------------------------------------------------------------
def slide7():
    fig, ax = new_slide()
    add_fitted_image(ax, os.path.join(FIG_DIR, "fig3_architecture.png"))
    add_footer(ax, 7)
    return fig


# ----------------------------------------------------------------------------
# Slide 8 : Pseudocode (custom code slide) — alignments -> transcripts
# ----------------------------------------------------------------------------
# Verbatim code block. Preserve indentation and line order exactly.
CODE_LINES = [
    "# 1. group alignments by their EXACT splice-junction chain",
    "groups = {}                              # key = ordered list of (donor, acceptor)",
    "for aln in long_read_alignments:",
    "    chain = introns(aln)                 # splice junctions, in order",
    "    if chain == []: continue             # single-exon -> skip",
    "    g = groups[chain]                    # GROUP BY exact chain (hash key, not similarity)",
    "    g.depth += aln.weight",
    "    g.start  = min(g.start, aln.start)   # 5'/3' ends = outer envelope",
    "    g.end    = max(g.end,   aln.end)",
    "",
    "# 2. fold truncated fragments into their full-length parent",
    "for B in groups:",
    "    A = longest chain that CONTAINS B as a sub-path   # terminus within tol",
    "    if A: A.depth += B.depth ; drop B                 # ISM -> FSM",
    "",
    "# 3. emit one transcript per surviving chain",
    "read_chains = [ Transcript(exons(chain, g.start, g.end), g.depth)",
    "                for chain, g in groups ]",
    "",
    "# 4. add over the flow floor, gated   (output >= flow)",
    "output = flow_assembly(reads)            # StringTie-style baseline",
    "seen   = { introns(t) for t in output }",
    "for t in read_chains:                    # held out of the flow filters",
    "    if introns(t) in seen:            continue   # already there -> no-op",
    "    if any junction not canonical(t): continue   # GT-AG / GC-AG / AT-AC",
    "    if rt_switch(t) or t.depth < K:   continue   # realness gate",
    "    output.add(t)                        # novel + real -> add (never removes a flow tx)",
]

# The four step-header lines are wholly orange + bold.
STEP_HEADER_PREFIXES = ("# 1.", "# 2.", "# 3.", "# 4.")


def slide8():
    fig, ax = new_slide()

    # Title (navy, bold)
    ax.text(0.5, 0.915,
            "Alignments -> transcripts  (pseudocode)",
            color=NAVY, fontsize=22, fontweight="bold",
            ha="center", va="center")

    # Code block geometry. Left-aligned single column, fixed line spacing.
    code_x = 0.075          # left margin for the code
    top_y = 0.815           # y of the first code line (top), below the title
    bottom_limit = 0.155    # code must not go below this (above caption/footer)
    n_lines = len(CODE_LINES)
    # fixed line spacing so all ~26 lines fit between top_y and bottom_limit
    line_step = (top_y - bottom_limit) / (n_lines - 1)
    code_fontsize = 10.5

    # Use a monospace font so columns/indentation align.
    mono = {"family": "DejaVu Sans Mono"}

    y = top_y
    for line in CODE_LINES:
        if line == "":
            y -= line_step
            continue

        if line.startswith(STEP_HEADER_PREFIXES):
            # whole line orange + bold
            ax.text(code_x, y, line, color=ORANGE, fontsize=code_fontsize,
                    fontweight="bold", ha="left", va="center",
                    fontfamily="DejaVu Sans Mono")
            y -= line_step
            continue

        # Split at the FIRST '#': code part navy, comment part slate.
        hash_idx = line.find("#")
        if hash_idx == -1:
            # pure code line
            ax.text(code_x, y, line, color=NAVY, fontsize=code_fontsize,
                    ha="left", va="center", fontfamily="DejaVu Sans Mono")
        else:
            code_part = line[:hash_idx]
            comment_part = line[hash_idx:]
            # Render the full line transparently first to get exact glyph metrics
            # for where the comment starts; simpler: render code in navy at code_x
            # then render comment in slate at the same character column by padding
            # the code part with spaces (monospace -> column-accurate).
            # Navy code (keep trailing spaces so monospace columns hold).
            ax.text(code_x, y, code_part, color=NAVY, fontsize=code_fontsize,
                    ha="left", va="center", fontfamily="DejaVu Sans Mono")
            # Slate comment: re-emit the line but blank out the code chars so the
            # comment lands in its true column (monospace alignment preserved).
            padded_comment = (" " * hash_idx) + comment_part
            ax.text(code_x, y, padded_comment, color=SLATE,
                    fontsize=code_fontsize, ha="left", va="center",
                    fontfamily="DejaVu Sans Mono")
        y -= line_step

    # Small slate caption under the code block.
    ax.text(code_x, 0.108,
            "steps 1-3 = the read-chain path in the previous slide;  "
            "step 4 = hold-out + gate + union-back",
            color=SLATE, fontsize=12.5, ha="left", va="center")

    add_footer(ax, 8)
    return fig


# ----------------------------------------------------------------------------
# Slide 9 : figure (fig5_results)
# ----------------------------------------------------------------------------
def slide9():
    fig, ax = new_slide()
    add_fitted_image(ax, os.path.join(FIG_DIR, "fig5_results.png"))
    add_footer(ax, 9)
    return fig


# ----------------------------------------------------------------------------
# Slide 10 : Positioning vs prior art (custom text, two columns) — airy
# ----------------------------------------------------------------------------
def slide10():
    fig, ax = new_slide()

    # Title (navy, bold)
    ax.text(0.5, 0.895,
            "Positioning: read-coherence is the established collapse paradigm",
            color=NAVY, fontsize=21, fontweight="bold",
            ha="center", va="center")

    left_x = 0.075
    right_x = 0.555
    bullet_indent = 0.022

    header_y = 0.745
    first_bullet_y = 0.625
    line_step = 0.135

    left_header = "Prior art (per-molecule collapse)"
    left_bullets = [
        "IsoSeq collapse / cDNA_Cupcake —\ncollapse full-length reads by junction chains",
        "IsoQuant — per-read intron graph;\nalready beats StringTie on novel isoforms (lower FP)",
        "FLAIR / TALON / Mandalorion —\nper-read isoforms by splice structure",
        "So \"per-molecule > flow assembly on recall\"\nis already known",
    ]

    right_header = "What this adds"
    right_bullets = [
        "A provably-additive layer over a StringTie\nflow port (output ⊇ flow), in one tool",
        "Annotation-free realness gate (canonical incl.\nGC-AG/AT-AC + RT-switch + depth)",
        "The substrate for the novel part:\nparalog-copy resolution via PSVs →",
    ]

    def render_col(x, header, bullets):
        ax.text(x, header_y, header, color=ORANGE,
                fontsize=17, fontweight="bold", ha="left", va="center")
        # underline accent under header
        ax.plot([x, x + 0.355], [header_y - 0.05, header_y - 0.05],
                color=ORANGE, lw=1.4, alpha=0.55)
        y = first_bullet_y
        for b in bullets:
            ax.text(x, y, "•", color=ORANGE, fontsize=15,
                    ha="left", va="top")
            ax.text(x + bullet_indent, y, b, color=NAVY, fontsize=13.5,
                    ha="left", va="top", linespacing=1.3)
            y -= line_step

    render_col(left_x, left_header, left_bullets)
    render_col(right_x, right_header, right_bullets)

    # subtle vertical divider between columns
    ax.plot([0.515, 0.515], [0.155, 0.695], color=SLATE, lw=0.7, alpha=0.35)

    # tiny citation line (small slate, above the page footer)
    ax.text(0.5, 0.112,
            "IsoQuant, Nat Biotech 2023  ·  FLAIR2, Genome Biol 2024  ·  "
            "Paraphase, Nat Commun 2025",
            color=SLATE, fontsize=10.5, ha="center", va="center", style="italic")

    add_footer(ax, 10)
    return fig


# ----------------------------------------------------------------------------
# Slide 11 : figure (fig8_psv_calc) - "Calling a PSV"; figure carries its own
# title, so no separate slide title. Inserted right before the PSV-bridge.
# ----------------------------------------------------------------------------
def slide11():
    fig, ax = new_slide()
    add_fitted_image(ax, os.path.join(FIG_DIR, "fig8_psv_calc.png"))
    add_footer(ax, 11)
    return fig


# ----------------------------------------------------------------------------
# Slide 12 : figure (fig7_psv_bridge) - has its own title, no slide title
# ----------------------------------------------------------------------------
def slide12():
    fig, ax = new_slide()
    add_fitted_image(ax, os.path.join(FIG_DIR, "fig7_psv_bridge.png"))
    add_footer(ax, 12)
    return fig


# ----------------------------------------------------------------------------
# Slide 13 : figure (fig9_contribution) - the contribution-in-one-slide grid.
# The figure carries its own title, so no separate slide title.
# ----------------------------------------------------------------------------
def slide13():
    fig, ax = new_slide()
    add_fitted_image(ax, os.path.join(FIG_DIR, "fig9_contribution.png"))
    add_footer(ax, 13)
    return fig


def main():
    builders = [slide1, slide2, slide3, slide4,
                slide5, slide6, slide7, slide8,
                slide9, slide10, slide11, slide12,
                slide13]

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
