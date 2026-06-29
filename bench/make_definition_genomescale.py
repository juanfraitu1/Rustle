#!/usr/bin/env python3
"""make_definition_genomescale.py — "Part 2" figures: the READ-ONLY definition at genome scale,
plus the ORTHOGONAL (independent) validation. Matches make_definition_visuals.py style/palette.

IMPORTANT framing: the family definition uses ONLY the IsoSeq reads (de-novo loci → ~R → ~B → 196
families). The cDNA all-vs-all + protein homology NEVER see the reads — they are an ORTHOGONAL
witness used to validate, NOT part of the definition.

Generates: fig_readdef.png, fig_subB.png, fig_ortho.png, fig_gallery.png
Run: /home/juanfra/miniforge3/bin/python bench/make_definition_genomescale.py
"""
import os
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, Rectangle, FancyArrow
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
NAVY, TEAL, ORANGE, GREY = "#1F2D5A", "#127A6E", "#D96B27", "#9aa0a6"
RED = "#b2362f"


def _box(ax, x, y, w, h, fc, ec, lw=1.6):
    ax.add_patch(FancyBboxPatch((x, y), w, h, boxstyle="round,pad=0.04", facecolor=fc, edgecolor=ec, linewidth=lw))


def _arrow(ax, x, y, dx):
    ax.add_patch(FancyArrow(x, y, dx, 0, width=0.02, head_width=0.16, head_length=0.16,
                            length_includes_head=True, color=GREY))


def fig_readdef():
    """The definition scales — built ONLY from the reads."""
    fig, ax = plt.subplots(figsize=(13, 6.4)); ax.axis("off"); ax.set_xlim(0, 13); ax.set_ylim(0, 7)
    ax.text(6.5, 6.6, "Does it scale?  —  the whole definition is built ONLY from the IsoSeq reads",
            ha="center", fontsize=15, weight="bold", color=NAVY)
    steps = [
        ("IsoSeq\nreads", "full-length\nmolecules", "#fdeee6", ORANGE),
        ("~101k DE-NOVO loci", "assembled from reads\n(NO annotation)", "#eef5f4", TEAL),
        ("~R   read-conflict", "reads TIE loci that\ncan't be told apart", "#eef5f4", TEAL),
        ("~B   shared backbone", "copy MODELS (read\nexon-unions) align", "#eef5f4", TEAL),
        ("196 multi-copy\nFAMILIES", "genome-wide,\nfrom reads alone", "#e8f0ea", "#2e7d32"),
    ]
    x = 0.35
    w = [1.9, 2.35, 2.25, 2.25, 2.2]
    for i, (t, sub, fc, ec) in enumerate(steps):
        if i:
            _arrow(ax, x - 0.28, 4.95, 0.25)
        _box(ax, x, 4.2, w[i], 1.55, fc, ec, 2.0 if i in (0, 4) else 1.5)
        ax.text(x + w[i] / 2, 5.35, t, ha="center", fontsize=9.8, weight="bold",
                color=NAVY if i not in (4,) else "#1b5e20")
        ax.text(x + w[i] / 2, 4.62, sub, ha="center", fontsize=8.2, color="#444")
        x += w[i] + 0.28
    ax.text(6.5, 3.2, "No cDNA all-vs-all, no protein, no DNA homology enters the definition — only what the reads say.",
            ha="center", fontsize=12, weight="bold", color=TEAL)
    # mini stats
    ax.text(6.5, 2.2, "1,136 de-novo loci fall into 196 families  ·  largest 47 loci  ·  73 are 2-copy, 24 are >10-copy",
            ha="center", fontsize=11, color="#333")
    ax.text(6.5, 1.4, "validated on a 17-candidate panel with hand-checked truth:  precision = recall = 1.00",
            ha="center", fontsize=11.5, weight="bold", color=NAVY)
    fig.savefig(os.path.join(HERE, "fig_readdef.png"), dpi=145, bbox_inches="tight")
    plt.close(fig)


def fig_subB():
    """~B sharpened: copy MODELS (from reads) share the whole backbone, not a domain."""
    fig, ax = plt.subplots(1, 2, figsize=(13, 5.6))
    rng = np.random.default_rng(3)
    a = ax[0]
    n = 800
    d = rng.uniform(0, 1, n)
    a.scatter(d, d + rng.normal(0, 0.012, n), s=3, color=TEAL)
    a.plot([0, 1], [0, 1], color=TEAL, lw=0.6, alpha=0.4)
    a.set_title("REAL family: copy models share the WHOLE backbone\nreciprocal coverage high  →  ~B ✓",
                fontsize=12, weight="bold", color=NAVY)
    a.set_xlabel("copy A model (read exon-union)"); a.set_ylabel("copy B model")
    b = ax[1]
    seg = rng.uniform(0.35, 0.6, 220)
    b.scatter(seg, seg + 0.05 + rng.normal(0, 0.01, 220), s=3, color=RED)
    b.plot([0, 1], [0, 1], color=GREY, lw=0.5, ls=":", alpha=0.5)
    b.axvspan(0.35, 0.6, color=RED, alpha=0.06)
    b.set_title("BRIDGE: loci share only a sub-region\npartial coverage  →  ~B ✗  (cut)",
                fontsize=12, weight="bold", color=RED)
    b.set_xlabel("locus A model"); b.set_ylabel("locus B model")
    for x in ax:
        x.set_xlim(0, 1); x.set_ylim(0, 1); x.set_xticks([]); x.set_yticks([])
    fig.suptitle("~B is a SECOND minimap2 run — copy-model vs copy-model (both built from reads): "
                 "a diagonal over the WHOLE of each = shared backbone; a partial diagonal = a bridge, cut.",
                 fontsize=12, weight="bold", color=NAVY, y=1.02)
    fig.savefig(os.path.join(HERE, "fig_subB.png"), dpi=145, bbox_inches="tight")
    plt.close(fig)


def fig_complete():
    """The complete, property-based definition: ~B DEFINES the family; expression & identifiability
    are LABELS, not gates. Identifiability (resolving copies) is next meeting (PSVs)."""
    fig, ax = plt.subplots(figsize=(13, 6.9)); ax.axis("off"); ax.set_xlim(0, 13); ax.set_ylim(0, 7)
    ax.text(6.5, 6.65, "A family is a STRUCTURAL fact — it holds whether or not it is expressed",
            ha="center", fontsize=15, weight="bold", color=NAVY)
    # center: ~B defines the family
    _box(ax, 4.2, 4.5, 4.6, 1.7, "#e8f0ea", "#2e7d32", 2.4)
    ax.text(6.5, 5.78, "~B  —  copies share a", ha="center", fontsize=12, weight="bold", color="#1b5e20")
    ax.text(6.5, 5.42, "WHOLE-GENE backbone", ha="center", fontsize=12, weight="bold", color="#1b5e20")
    ax.text(6.5, 4.95, "DEFINES the family   →   ~1,000 structural families", ha="center",
            fontsize=10.5, weight="bold", color=NAVY)
    ax.text(6.5, 4.62, "(tested by copy-vs-copy alignment — today)", ha="center", fontsize=8.5, color="#444", style="italic")
    # left: expression = a label (observation)
    ax.add_patch(FancyArrow(4.6, 4.45, -1.4, -1.0, width=0.015, head_width=0.16, head_length=0.16,
                            length_includes_head=True, color=GREY))
    _box(ax, 0.4, 1.5, 4.5, 2.0, "#eef5f4", TEAL, 1.8)
    ax.text(2.65, 3.18, "EXPRESSION  (observation)", ha="center", fontsize=11, weight="bold", color=TEAL)
    ax.text(2.65, 2.68, "EXPRESSED  495   (seen in reads)", ha="center", fontsize=10, color="#222")
    ax.text(2.65, 2.28, "SILENT  518   (not in this sample)", ha="center", fontsize=10, color="#222")
    ax.text(2.65, 1.78, "a LABEL, not a gate — silent\nfamilies are families too", ha="center",
            fontsize=8.5, weight="bold", color=NAVY)
    # right: identifiability = next meeting
    ax.add_patch(FancyArrow(8.4, 4.45, 1.4, -1.0, width=0.015, head_width=0.16, head_length=0.16,
                            length_includes_head=True, color=GREY))
    _box(ax, 8.1, 1.5, 4.5, 2.0, "#fdeee6", ORANGE, 2.0)
    ax.text(10.35, 3.18, "IDENTIFIABILITY  (~R)", ha="center", fontsize=11, weight="bold", color=ORANGE)
    ax.text(10.35, 2.62, "can we tell the copies apart\nand assign reads to them?", ha="center", fontsize=9.5, color="#222")
    ax.text(10.35, 1.82, "▶  NEXT MEETING: resolve with\nPSVs & copy-assignment", ha="center",
            fontsize=9.5, weight="bold", color=RED)
    ax.text(6.5, 0.55, "Today = the FAMILY (structure, ~B).   Next = resolving the COPIES (identifiability, PSVs).",
            ha="center", fontsize=11.5, weight="bold", color=NAVY)
    fig.savefig(os.path.join(HERE, "fig_complete.png"), dpi=145, bbox_inches="tight")
    plt.close(fig)


def fig_alncrit():
    """What counts as a GOOD vs BAD copy-vs-copy (~B) alignment — REAL dot-plots + the criteria."""
    from make_definition_visuals import load_seqs, dotplot
    fig = plt.figure(figsize=(13, 7.0))
    gs = fig.add_gridspec(2, 3, height_ratios=[1.05, 0.95], hspace=0.95, wspace=0.34,
                          left=0.05, right=0.96, top=0.86, bottom=0.05)
    fig.suptitle("What makes a GOOD vs BAD copy-vs-copy alignment?   (real dot-plots — every dot is a shared 12-mer)",
                 fontsize=14, weight="bold", color=NAVY, y=0.97)
    # REAL examples: GOOD = two array copies; BRIDGE = CD200R1~CD200R1L (share one Ig domain only);
    # NO HOMOLOGY = OCLN~SEPTIN7.
    panels = [
        ("GOOD   ~B ✓", "LOC129528712", "LOC129528713", TEAL,
         "two array copies — one continuous\ndiagonal over the WHOLE of both", "copy A", "copy B"),
        ("BRIDGE   ~B ✗", "CD200R1", "CD200R1L", RED,
         "CD200R1 ~ CD200R1L share one DOMAIN:\nmost of CD200R1 unmatched → cut", "CD200R1 (bp)", "CD200R1L (bp)"),
        ("NO HOMOLOGY   ~B ✗", "OCLN", "SEPTIN7", "#777",
         "OCLN vs SEPTIN7 — no diagonal,\nscattered hits (~0% covered)", "OCLN (bp)", "SEPTIN7 (bp)"),
    ]
    seqs = load_seqs([g for _, a, b, *_ in panels for g in (a, b)])
    for i, (ttl, a, b, col, sub, xl, yl) in enumerate(panels):
        ax = fig.add_subplot(gs[0, i])
        A, B = seqs.get(a, ""), seqs.get(b, "")
        if A and B:
            dotplot(ax, A, B, color=col)
            ax.set_xlim(0, len(A)); ax.set_ylim(0, len(B))
        ax.set_title(ttl, color=col if col != "#777" else "#555", weight="bold", fontsize=12)
        ax.text(0.5, -0.32, sub, transform=ax.transAxes, ha="center", fontsize=8.3, color="#333")
        ax.set_xlabel(xl, fontsize=8); ax.set_ylabel(yl, fontsize=8)
        ax.tick_params(labelsize=6.5); ax.grid(alpha=0.18)
    # criteria table
    axt = fig.add_subplot(gs[1, :]); axt.axis("off"); axt.set_xlim(0, 13); axt.set_ylim(0, 5)
    rows = [
        ("what minimap2 reports", "GOOD alignment  (~B ✓ → same family)", "BAD alignment  (~B ✗ → bridge, cut)"),
        ("the diagonal", "one continuous, full-length", "a short fragment, scattered, or none"),
        ("reciprocal coverage  (aligned span / length)", "≥ 0.30 of BOTH copies (often whole-gene)", "< 0.30 — only a domain, or 0"),
        ("identity  (matches / block length)", "≥ ~85%  (our copies are ≥ 90%)", "< 80% → asm20 cannot even seed"),
        ("means", "the copies share a BACKBONE", "shared domain only, or no homology"),
    ]
    colx = [0.15, 4.7, 8.85]; colw = [4.4, 4.0, 4.0]
    yh = 0.92; y = 4.55
    for ri, (a, b, c) in enumerate(rows):
        if ri == 0:
            axt.add_patch(Rectangle((0.1, y - yh + 0.06), 12.8, yh, facecolor=NAVY))
            cols = [("white", a), ("#bfe3dc", b), ("#f3c9c0", c)]
            fw = "bold"
        else:
            sh = "#eef5f4" if ri % 2 else "#ffffff"
            axt.add_patch(Rectangle((0.1, y - yh + 0.06), 12.8, yh, facecolor=sh, edgecolor="#dde"))
            cols = [(NAVY, a), (TEAL, b), (RED, c)]
            fw = "normal"
        for (cc, txt), x, w in zip(cols, colx, colw):
            axt.text(x, y - yh / 2 + 0.06, txt, fontsize=9 if ri else 9.5, color=cc,
                     va="center", weight="bold" if (ri == 0 or fw == "bold") else "normal")
        y -= yh
    axt.text(6.5, 0.18, "~B is a SEPARATE minimap2 run (copy-model vs copy-model) — it is the test that turns a read-tie (~R) into a real family.",
             ha="center", fontsize=10, weight="bold", color=NAVY)
    fig.savefig(os.path.join(HERE, "fig_alncrit.png"), dpi=145, bbox_inches="tight")
    plt.close(fig)


def fig_clique():
    """From pairwise ~B to the whole family: a CLIQUE (all members mutually share a backbone),
    not a connected chain (which over-merges via a bridge)."""
    fig, ax = plt.subplots(1, 2, figsize=(13, 5.9))
    # LEFT: a 5-clique (every pair an edge)
    th = np.linspace(0, 2 * np.pi, 6)[:-1] + np.pi / 2
    px, py = np.cos(th), np.sin(th)
    a = ax[0]
    for i in range(5):
        for j in range(i + 1, 5):
            a.plot([px[i], px[j]], [py[i], py[j]], color=TEAL, lw=2.2, zorder=1)
    a.scatter(px, py, s=1100, color="#e8f0ea", edgecolors=TEAL, linewidths=2.6, zorder=2)
    for i in range(5):
        a.text(px[i], py[i], str(i + 1), ha="center", va="center", fontsize=13, weight="bold", color=NAVY)
    a.set_title("FAMILY = a CLIQUE\nevery pair shares a backbone (~B ✓)", fontsize=13, weight="bold", color=TEAL)
    a.text(0, -1.55, "all 5 copies mutually align → ONE family\n(81% of families are full cliques)",
           ha="center", fontsize=9.5, color="#333")
    a.set_xlim(-1.5, 1.5); a.set_ylim(-1.9, 1.5); a.axis("off")
    # RIGHT: two cliques joined by a single bridge -> connected, but NOT a clique
    b = ax[1]
    Ax, Ay = [-1.65, -1.05, -1.75, -1.0], [0.55, 0.75, -0.35, -0.55]
    Bx, By = [1.6, 1.05, 1.55], [0.6, -0.25, -0.65]
    for i in range(4):
        for j in range(i + 1, 4):
            b.plot([Ax[i], Ax[j]], [Ay[i], Ay[j]], color=TEAL, lw=2.0, zorder=1)
    for i in range(3):
        for j in range(i + 1, 3):
            b.plot([Bx[i], Bx[j]], [By[i], By[j]], color=NAVY, lw=2.0, zorder=1)
    b.plot([Ax[1], Bx[1]], [Ay[1], By[1]], color=RED, lw=2.6, ls="--", zorder=1)
    b.scatter(Ax, Ay, s=780, color="#eef5f4", edgecolors=TEAL, linewidths=2.3, zorder=2)
    b.scatter(Bx, By, s=780, color="#eaecf4", edgecolors=NAVY, linewidths=2.3, zorder=2)
    b.text((Ax[1] + Bx[1]) / 2, (Ay[1] + By[1]) / 2 + 0.18, "bridge\n(shared domain/repeat)",
           color=RED, fontsize=9, weight="bold", ha="center")
    b.text(-1.35, 1.15, "family A", color=TEAL, fontsize=10.5, weight="bold", ha="center")
    b.text(1.35, 1.1, "family B", color=NAVY, fontsize=10.5, weight="bold", ha="center")
    b.set_title("BRIDGE = connected, but NOT a clique\nshares a backbone with only PART → cut",
                fontsize=13, weight="bold", color=RED)
    b.text(0, -1.55, "connected-components would MERGE A+B via the bridge;\nthe clique / dense-community rule cuts it → two families",
           ha="center", fontsize=9.5, color="#333")
    b.set_xlim(-2.0, 2.0); b.set_ylim(-1.9, 1.5); b.axis("off")
    fig.suptitle("From pairwise ~B to the whole family: include a member only if it shares a backbone with the REST (a clique), not just one (a chain)",
                 fontsize=12.5, weight="bold", color=NAVY, y=1.02)
    from matplotlib.lines import Line2D
    leg = [Line2D([0], [0], marker="o", color="none", markerfacecolor="#e8f0ea", markeredgecolor=TEAL,
                  markersize=14, markeredgewidth=2.4, label="  a COPY  (its exon-union copy model)"),
           Line2D([0], [0], color=TEAL, lw=2.8, label="  a ~B edge:  the two copies share a backbone")]
    fig.legend(handles=leg, loc="lower center", ncol=2, fontsize=11, frameon=True, edgecolor="#ccc",
               bbox_to_anchor=(0.5, -0.06), handletextpad=0.3, columnspacing=2.5)
    fig.subplots_adjust(bottom=0.1)
    fig.savefig(os.path.join(HERE, "fig_clique.png"), dpi=145, bbox_inches="tight")
    plt.close(fig)


def fig_ortho():
    """Orthogonal validation: independent DNA/protein homology that never sees the reads."""
    fig, ax = plt.subplots(figsize=(13, 6.6)); ax.axis("off"); ax.set_xlim(0, 13); ax.set_ylim(0, 7)
    ax.text(6.5, 6.65, "Orthogonal validation — an independent witness that never sees the reads",
            ha="center", fontsize=15, weight="bold", color=NAVY)
    # left: the definition (reads)
    _box(ax, 0.4, 3.7, 5.4, 2.4, "#eef5f4", TEAL, 2.2)
    ax.text(3.1, 5.75, "THE DEFINITION  (reads only)", ha="center", fontsize=12, weight="bold", color=TEAL)
    ax.text(3.1, 5.25, "IsoSeq reads → de-novo loci → ~R ∩ ~B", ha="center", fontsize=10.5, color="#222")
    ax.text(3.1, 4.7, "196 multi-copy families", ha="center", fontsize=12, weight="bold", color=NAVY)
    ax.text(3.1, 4.2, "RABL2 · APOBEC3 · MAGEA · RFPL4A · …", ha="center", fontsize=9, color="#333")
    # right: independent homology
    _box(ax, 7.2, 3.7, 5.4, 2.4, "#eaecf4", NAVY, 2.2)
    ax.text(9.9, 5.75, "INDEPENDENT homology  (annotation)", ha="center", fontsize=12, weight="bold", color=NAVY)
    ax.text(9.9, 5.25, "cDNA all-vs-all + protein (mmseqs)", ha="center", fontsize=10.5, color="#222")
    ax.text(9.9, 4.7, "paralog families — never uses the reads", ha="center", fontsize=10, color="#333")
    ax.text(9.9, 4.2, "a separate ground truth", ha="center", fontsize=9.5, style="italic", color=GREY)
    # converge
    ax.add_patch(FancyArrow(5.85, 4.9, 1.3, 0, width=0.02, head_width=0.0, color=GREY))
    ax.text(6.5, 3.25, "→  they AGREE on the real families  ←", ha="center", fontsize=13, weight="bold", color="#1b5e20")
    # caution box
    _box(ax, 1.2, 0.5, 10.6, 2.2, "#fdeee6", ORANGE, 1.6)
    ax.text(6.5, 2.35, "Caution — and why the read-based ~B is the right test:", ha="center", fontsize=11.5, weight="bold", color=ORANGE)
    ax.text(6.5, 1.78, "Sequence homology ALONE over-merges (a 238-gene blob of unrelated / lncRNA / retro genes bridged by a shared repeat).",
            ha="center", fontsize=10, color="#333")
    ax.text(6.5, 1.28, "To clean it you must add: whole-gene reciprocal coverage · drop intronless retrocopies · drop antisense pairs —",
            ha="center", fontsize=10, color="#333")
    ax.text(6.5, 0.82, "exactly the conditions the read-based ~B (copy-model backbone) enforces FOR FREE.",
            ha="center", fontsize=10.5, weight="bold", color=NAVY)
    fig.savefig(os.path.join(HERE, "fig_ortho.png"), dpi=145, bbox_inches="tight")
    plt.close(fig)


def fig_gallery():
    """Real families the READ definition finds (orthogonally confirmed)."""
    fig, ax = plt.subplots(figsize=(13, 7.6)); ax.axis("off"); ax.set_xlim(0, 13); ax.set_ylim(0, 14)
    ax.text(6.5, 13.5, "Real families the definition covers  —  validated, named examples",
            ha="center", fontsize=15, weight="bold", color=NAVY)
    groups = [
        ("Tandem arrays (co-located copies)", TEAL, [
            ("RFPL4A array", "5 tandem copies · NC_073244.2"),
            ("MAGEA9 / A2 / A3 clusters", "cancer-testis antigen arrays · chrX"),
            ("TMEM185A + LOC copy", "83 tied reads")]),
        ("Dispersed paralogs (different loci / chromosomes)", ORANGE, [
            ("RABL2A / RABL2B", "47 tied reads · two chromosomes"),
            ("APOBEC3D / APOBEC3F", "APOBEC3 cytidine-deaminase cluster"),
            ("AK6 · CCDC196  (LOC-named copies)", "found despite Gnomon LOC naming")]),
        ("Scale & validation", NAVY, [
            ("1,013 structural families (~B)", "495 expressed · 518 silent — all are families"),
            ("17-candidate panel", "precision = recall = 1.00 (hand-checked)"),
            ("each a real paralog set", "confirmed at the protein level")]),
        ("Correctly EXCLUDED (~B ✗ — no shared backbone)", RED, [
            ("OCLN ~ SEPTIN7", "share no whole-gene backbone → not a family"),
            ("CD200R1 ~ CD200R1L (share one domain)", "no whole-gene diagonal → cut"),
            ("processed retrocopies (intronless)", "no shared intron-exon backbone → cut")]),
    ]
    y = 12.6
    for title, col, rows in groups:
        ax.add_patch(Rectangle((0.2, y - 0.05), 12.6, 0.55, facecolor=col))
        ax.text(0.4, y + 0.22, title, color="white", fontsize=11, weight="bold", va="center")
        y -= 0.75
        for name, desc in rows:
            ax.text(0.55, y, "•  " + name, fontsize=10, color=NAVY, weight="bold", va="center")
            ax.text(6.3, y, desc, fontsize=9, color="#333", va="center")
            y -= 0.62
        y -= 0.35
    ax.text(6.5, 0.35, "The family is structural (~B): ~1,000 families, expressed or silent — tandem arrays, dispersed paralogs, immune loci. "
            "Today we DEFINE them; next we RESOLVE the copies.", ha="center", fontsize=10.5, weight="bold", color=TEAL)
    fig.savefig(os.path.join(HERE, "fig_gallery.png"), dpi=145, bbox_inches="tight")
    plt.close(fig)


if __name__ == "__main__":
    fig_complete(); fig_alncrit(); fig_ortho(); fig_gallery()
    print("[+] wrote fig_complete.png, fig_alncrit.png, fig_ortho.png, fig_gallery.png")
