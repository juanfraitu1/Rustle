#!/usr/bin/env python3
"""Render THESIS_CONTRIBUTION.md into a clean, readable portrait document PDF.

Layout philosophy: a real document (not slides). Generous margins, comfortable
line spacing, justified/wrapped body via textwrap. Navy title, navy/orange bold
section headers, a boxed identifiability callout, and a tight positioning list.

The content does not fit legibly on one page at >= 9pt, so we lay it out across
TWO clean pages with automatic page breaks (matplotlib PdfPages). A flowing
cursor model keeps margins and line spacing consistent across both pages and
never lets a section header dangle at the foot of a page.
"""

import os
import textwrap

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.patches import FancyBboxPatch

# ---------------------------------------------------------------- palette ----
NAVY = "#22313f"
SLATE = "#8898aa"
ORANGE = "#e8590c"
GREEN = "#1e8449"
WHITE = "#ffffff"
CALLOUT_BG = "#fdf2ec"   # very faint orange wash for the theorem box
CALLOUT_EDGE = ORANGE

HERE = os.path.dirname(os.path.abspath(__file__))
PDF_PATH = os.path.join(HERE, "THESIS_CONTRIBUTION.pdf")

# DejaVu Sans (matplotlib default sans) has the unicode glyphs we use.
plt.rcParams["font.family"] = "DejaVu Sans"
plt.rcParams["pdf.fonttype"] = 42  # embed TrueType so text stays selectable

# ----------------------------------------------------------- page geometry ---
PAGE_W, PAGE_H = 8.5, 11.0           # inches (portrait)
DPI = 200
M_LEFT = 0.090                       # margins as figure fraction
M_RIGHT = 0.090
M_TOP = 0.945
M_BOTTOM = 0.060
CONTENT_W = 1.0 - M_LEFT - M_RIGHT

# Font sizes (points)
FS_TITLE = 17
FS_SUB = 11
FS_HEAD = 12.5
FS_BODY = 10.0
FS_BULLET = 10.0
FS_REF = 9.0

# Vertical rhythm (figure-fraction per line at the body font size)
LINE = 0.0168          # body line height
HEAD_GAP_BEFORE = 0.0150
HEAD_GAP_AFTER = 0.0050
PARA_GAP = 0.0110
BULLET_INDENT = 0.024
HANG_INDENT = 0.024

# Characters-per-line for wrapping. Calibrated against measured DejaVu Sans
# glyph widths at each font size (content width ~1394px @200dpi); values include
# a safety margin so no line touches the right margin.
WRAP_BODY = 84       # full content width @ 10pt
WRAP_BULLET = 78     # content width minus bullet hanging indent @ 10pt
WRAP_CALLOUT = 82    # inside the callout box @ 10pt
WRAP_REF = 96        # full content width @ 9pt

PNG_PREVIEW = False  # set True via --png to also emit per-page PNG previews


class Doc:
    """A simple top-down flowing document across one or more pages."""

    def __init__(self):
        self.pages = []
        self.fig = None
        self.ax = None
        self.y = 0.0
        self._new_page()

    def _new_page(self):
        fig = plt.figure(figsize=(PAGE_W, PAGE_H), dpi=DPI)
        fig.patch.set_facecolor(WHITE)
        ax = fig.add_axes([0, 0, 1, 1])
        ax.set_facecolor(WHITE)
        ax.set_xlim(0, 1)
        ax.set_ylim(0, 1)
        ax.axis("off")
        self.pages.append(fig)
        self.fig = fig
        self.ax = ax
        self.y = M_TOP

    def need(self, height):
        """Ensure `height` fig-fraction is available; else start a new page."""
        if self.y - height < M_BOTTOM:
            self._new_page()

    # -- primitives --------------------------------------------------------
    def text(self, x, s, *, fontsize, color, weight="normal", style="normal",
             zorder=2):
        self.ax.text(x, self.y, s, fontsize=fontsize, color=color,
                     weight=weight, style=style, ha="left", va="top",
                     transform=self.ax.transAxes, zorder=zorder)

    def wrapped(self, s, *, x=M_LEFT, width=WRAP_BODY, fontsize=FS_BODY,
                color=NAVY, weight="normal", style="normal", line=LINE,
                hang=0.0):
        lines = textwrap.wrap(s, width=width, break_long_words=False,
                              break_on_hyphens=False)
        for i, ln in enumerate(lines):
            self.need(line)
            xx = x if i == 0 else x + hang
            self.text(xx, ln, fontsize=fontsize, color=color, weight=weight,
                      style=style)
            self.y -= line

    def header(self, s):
        # keep header with at least the first wrapped body line
        self.y -= HEAD_GAP_BEFORE
        self.need(LINE * 3)
        self.text(M_LEFT, s, fontsize=FS_HEAD, color=NAVY, weight="bold")
        self.y -= LINE
        self.y -= HEAD_GAP_AFTER

    def bullet(self, s, *, emphasize=False):
        bx = M_LEFT + BULLET_INDENT
        tx = M_LEFT + BULLET_INDENT + HANG_INDENT
        if emphasize:
            body_color, weight, mk, mk_color = ORANGE, "bold", "▶", ORANGE
        else:
            body_color, weight, mk, mk_color = NAVY, "normal", "•", ORANGE
        lines = textwrap.wrap(s, width=WRAP_BULLET, break_long_words=False,
                              break_on_hyphens=False)
        if not lines:
            return
        self.need(LINE)
        self.text(bx, mk, fontsize=FS_BULLET, color=mk_color, weight="bold")
        for i, ln in enumerate(lines):
            if i > 0:
                self.need(LINE)
            self.text(tx, ln, fontsize=FS_BULLET, color=body_color,
                      weight=weight)
            self.y -= LINE

    def gap(self, h=PARA_GAP):
        self.y -= h

    def callout(self, title, body):
        body_wrap = textwrap.wrap(body, width=WRAP_CALLOUT,
                                  break_long_words=False, break_on_hyphens=False)
        pad_top, pad_bottom, pad_x = 0.012, 0.013, 0.016
        title_h = LINE + 0.004
        body_h = len(body_wrap) * LINE
        box_h = pad_top + title_h + body_h + pad_bottom
        # keep the whole callout on a single page
        self.y -= HEAD_GAP_BEFORE
        self.need(box_h + 0.004)
        box_top = self.y
        box_bottom = box_top - box_h
        box_left = M_LEFT - 0.008
        box_right = 1 - M_RIGHT + 0.008
        box = FancyBboxPatch(
            (box_left, box_bottom), box_right - box_left, box_h,
            boxstyle="round,pad=0.004,rounding_size=0.009",
            linewidth=1.5, edgecolor=CALLOUT_EDGE, facecolor=CALLOUT_BG,
            transform=self.ax.transAxes, clip_on=False, zorder=1)
        self.ax.add_patch(box)
        self.y = box_top - pad_top
        self.text(box_left + pad_x, title, fontsize=FS_HEAD, color=ORANGE,
                  weight="bold")
        self.y -= title_h
        for ln in body_wrap:
            self.text(box_left + pad_x, ln, fontsize=FS_BODY, color=NAVY)
            self.y -= LINE
        self.y = box_bottom

    def save(self, png_preview=False):
        with PdfPages(PDF_PATH) as pp:
            for f in self.pages:
                pp.savefig(f, facecolor=WHITE)
        if png_preview:
            for i, f in enumerate(self.pages):
                f.savefig(os.path.join(HERE, "_verify_p%d.png" % (i + 1)),
                          facecolor=WHITE, dpi=150)
        for f in self.pages:
            plt.close(f)


def render():
    d = Doc()

    # ---- title block -----------------------------------------------------
    # Title is wider than the content area at a legible size, so set it on two
    # lines rather than shrinking it.
    d.text(M_LEFT, "Copy-Resolved Isoform Assembly on a",
           fontsize=FS_TITLE, color=NAVY, weight="bold")
    d.y -= 0.029
    d.text(M_LEFT, "Spliced Family Variation Graph",
           fontsize=FS_TITLE, color=NAVY, weight="bold")
    d.y -= 0.030
    d.text(M_LEFT,
           "Proposed thesis contribution — rustle (long-read transcript "
           "assembly)",
           fontsize=FS_SUB, color=SLATE, style="italic")
    d.y -= 0.022
    d.ax.plot([M_LEFT, 1 - M_RIGHT], [d.y, d.y], color=NAVY, lw=1.2,
              transform=d.ax.transAxes, clip_on=False)
    d.y -= 0.014

    # ---- Gap -------------------------------------------------------------
    d.header("Gap")
    d.wrapped(
        "Per-molecule long-read isoform assembly is mature: \"collapse\" by "
        "intron chain (cDNA_Cupcake, IsoSeq, IsoQuant, FLAIR, LRAA). "
        "Graph/allele-aware methods are HAPLOTYPIC: per-molecule allelic "
        "isoforms (HapIso), or haplotype-aware quantification of KNOWN "
        "transcripts on spliced pangenome graphs (RPVG, pantas). PSV-based "
        "PARALOG resolution exists only at the DNA/gene level (Paraphase). "
        "The unoccupied cell: de-novo, per-molecule isoform assembly resolved "
        "to paralog COPIES.")
    d.gap()

    # ---- Setting ---------------------------------------------------------
    d.header("Setting")
    d.wrapped(
        "A gene family of unknown copy number; copies near-identical, "
        "separated by paralog-sequence variants (PSVs = fixed inter-copy base "
        "differences). Long reads multimap across copies; each read is a "
        "SPLICED PATH (intron chain) carrying an ALLELE VECTOR at PSV columns. "
        "Build a spliced family variation graph G: nodes = exonic segments + "
        "PSV columns (allele variants); a TRANSCRIPT = a path specifying "
        "junctions and allele choices; a COPY = a maximal allele-consistent "
        "haplotype with its isoform set.")
    d.gap()

    # ---- Problem ---------------------------------------------------------
    d.header("Problem (joint decomposition)")
    d.wrapped(
        "Given the multiset of read paths on G, find the minimum set of "
        "(copy, isoform) paths and a read assignment such that:")
    d.gap(0.002)
    d.bullet(
        "(read-coherence) every emitted isoform is spanned by at least one "
        "molecule whose junctions match;")
    d.bullet(
        "(allele-linkage) each read's alleles are consistent along its single "
        "path;")
    d.bullet(
        "(copy-consistency) all isoforms of a copy share one "
        "allele-haplotype;")
    d.bullet(
        "(shared evidence) reads without decisive PSV evidence are apportioned "
        "across copies, not forced to one.")
    d.gap(0.002)
    d.wrapped(
        "This is a constrained path-cover / flow-decomposition on G under "
        "linkage constraints. It FLIPS the resolve-one (facility-location) "
        "frame: multimappers are pooled shared evidence, disambiguated "
        "per-molecule only where PSVs are decisive.")
    d.gap()

    # ---- Identifiability theorem (boxed callout) -------------------------
    d.callout(
        "Identifiability theorem (the guarantee)",
        "Copies c_i, c_j are DISTINGUISHABLE iff they differ at >= K PSV "
        "columns that are jointly spanned by >= m reads above the per-base "
        "error floor e. Under this condition the (copy, isoform) "
        "decomposition, restricted to spanned isoforms, is unique up to the "
        "unidentifiable residue; below it the copies are provably merged. This "
        "characterizes the recoverable regime (K independent error-agreements "
        "are improbable at long-read error rates) and turns the empirical "
        "\"indistinguishability wall\" into a theorem.")
    d.gap()

    # ---- What exists vs what is new -------------------------------------
    d.header("What exists vs. what is new")
    d.wrapped(
        "rustle already provides the ingredients: family variation graph, PSV "
        "columns + linkage, EM read-reweighting, read-chain (per-molecule) "
        "extraction, copy certificates. NEW: the joint formulation — a "
        "single per-molecule key (intron chain (+) PSV-haplotype) decomposed "
        "on G — the identifiability theorem, and the regime characterization.")
    d.gap()

    # ---- Evaluation ------------------------------------------------------
    d.header("Evaluation")
    d.wrapped(
        "Restrict to identifiable families (>= K PSVs); report copy-resolved "
        "exact-isoform (FSM) recovery vs StringTie / IsoQuant / per-copy "
        "truth; present the non-identifiable regime as a limit-of-the-data "
        "result with honest identifiability accounting.")
    d.gap()

    # ---- Positioning -----------------------------------------------------
    d.header("Positioning")
    d.bullet(
        "IsoQuant / Cupcake / FLAIR / LRAA : de-novo per-molecule isoform "
        "— but single-locus; paralogs out of scope.")
    d.bullet(
        "HapIso : de-novo per-molecule isoform split by variant — ALLELIC "
        "(2 haplotypes of one gene), not paralog copies.")
    d.bullet(
        "RPVG / pantas (spliced pangenome) : graph + variants — QUANTIFY "
        "known transcripts, allelic/population, not de-novo copy assembly.")
    d.bullet(
        "Paraphase : PSV-based paralog resolution — DNA/gene level, not "
        "isoforms.")
    d.bullet(
        "THIS THESIS : de-novo per-molecule isoform assembly x paralog-copy "
        "resolution, on one spliced family variation graph, with provable "
        "identifiability.", emphasize=True)
    d.gap()

    # ---- Refs ------------------------------------------------------------
    d.wrapped(
        "Refs: cDNA_Cupcake; IsoSeq; IsoQuant (Nat Biotech 2023); FLAIR2 "
        "(Genome Biol 2024); HapIso; RPVG / spliced pangenome (Nat Methods "
        "2022); pantas (PLOS CB 2024); Paraphase (Nat Commun 2025).",
        fontsize=FS_REF, color=SLATE, width=WRAP_REF, line=LINE * 0.96)

    d.save(png_preview=PNG_PREVIEW)
    return len(d.pages), d.y


if __name__ == "__main__":
    import sys
    PNG_PREVIEW = "--png" in sys.argv
    n_pages, final_y = render()
    print("rendered:", PDF_PATH)
    print("pages:", n_pages)
    print("final_y on last page (fig fraction): {:.4f}  (bottom margin {:.3f})"
          .format(final_y, M_BOTTOM))
    if final_y < M_BOTTOM - 1e-9:
        print("WARNING: last block overran bottom margin by",
              round(M_BOTTOM - final_y, 4))
