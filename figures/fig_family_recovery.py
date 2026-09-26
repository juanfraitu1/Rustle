"""Figure 7 — family recovery by Rustle's two modes, de novo and guided.

A Rustle-internal comparison: the two MODES of Rustle, both grouped by the same family rule (MCL with `--min-exonic-bp
1 --min-shared-exon-frac 0.60`). De novo IS Rustle's one default de novo family definition (user decision 2026-09-25:
reads -> seeded assembly loci -> one representative per locus, its "positional exon sum" -> families; the families
copy assignment consumes); guided = the same rule on the annotated gene and pseudogene bodies. It is not a comparison
with other tools. The thesis goal is to reduce the difference between the modes.

References are EXTERNAL (docs/PREREG_genome_wide_families_2026-09-25.md, Amendment 1): Ensembl Compara families
(duplications within primates; human; the headline), Soto et al. 2025 (human; not independent of the exon threshold),
the NPIP reference set (human chr16, an inset) and Liftoff copy pairs (the Fig. 8 self-lift; every species; the de
novo mode only, since Liftoff's extra copies are unannotated by construction). Protein-homology families are a
secondary reference of the supplementary figure fig7s_protein_homology only.

Genome-wide version (default): one genome-wide run per mode (de novo per sample, six samples; guided per species, it
reads no RNA), scored with `family_score --chrom ALL --per-family --pairwise`.
(a) Pairwise sensitivity and precision; (b) one-to-one bipartite sensitivity, precision and F (pooled), genome minus
    development contigs (open ring: minus every contig used for a family decision); (c) reference families exact /
    partial / missed, per mode (Liftoff rows: copy pairs in one family).
(d) The same genome-wide run broken down by contig (Compara families, human), development and reused contigs marked:
    no contig was picked.
(e-g) Per-family F, de novo vs guided: human A119b and human testis against Compara, human A119b against Soto 2025.

Development version (`fig7_scope dev`, or `make.py data fig7 --recorded`): per-chromosome runs on human chr16, chr2,
chr6, chr8, chr10 and gorilla NC_073244.2, NC_073234.2 (tables fig7_summary, fig7_per_family), drawn by the same
module when the genome-wide tables are absent: the main figure from the Compara / Soto / NPIP rows (human), the
supplement from the protein-homology rows (secondary; human and gorilla).
"""
from __future__ import annotations

import re
import sys
from pathlib import Path

import figlib
import _o1_recovery as R

GW_TABLES = ["fig7_gw_summary", "fig7_gw_contigs", "fig7_gw_per_family"]
GW_SUPPLEMENTARY = ["fig7_gw_compara_levels", "fig7_gw_bridge"]
MAIN_REFS = ["compara", "soto", "liftoff", "npip_u2"]      # genome-wide main figure (external references)
SUPP_REFS = ["protein_homology"]                            # supplement fig7s_protein_homology (secondary)
DEV_TABLES = ["fig7_summary", "fig7_per_family"]


def _gw_ready(data_dir: Path) -> bool:
    return all((Path(data_dir) / f"{t}.tsv").exists() for t in GW_TABLES)


DEV_CLAIM = (
    "Development tables (per-chromosome runs of both modes on human chr16, chr2, chr6, chr8, chr10; a Rustle-internal "
    "comparison of its two modes, not a tool comparison), scored against Ensembl Compara families (duplications "
    "within primates, restricted to each chromosome), Soto 2025 (not independent) and the NPIP reference set; the "
    "numbers are printed by `python3 figures/fig_family_recovery.py summary` and quoted in captions/fig7.md. The "
    "protein-homology rows (human and gorilla) are a secondary reference in the supplementary figure. The "
    "genome-wide version (every sample; Liftoff copy pairs for every species) replaces these tables when it is "
    "built.")
GW_CLAIM = (
    "Genome-wide, every sample; a Rustle-internal comparison of its two modes, not a tool comparison (pre-registered "
    "claims F7.1-F7.4, docs/PREREG_genome_wide_families_2026-09-25.md, Amendment 1; the numbers are printed by "
    "`python3 figures/fig_family_recovery.py summary` and quoted in captions/fig7.md): Rustle's default de novo "
    "families per sample vs guided families per species against Ensembl Compara families (primates), Soto 2025 and "
    "the NPIP reference set (human), and Liftoff copy pairs (every species; de novo only), on the genome minus "
    "development contigs, with a per-contig breakdown.")

META = {
    "id": "fig7",
    "title": "Family recovery by Rustle's two modes (de novo vs guided; Rustle-internal, not a tool comparison)",
    "claim": GW_CLAIM if _gw_ready(figlib.DATA_DIR) else DEV_CLAIM,
    "tables": GW_TABLES if _gw_ready(figlib.DATA_DIR) else DEV_TABLES,
    "supplementary_tables": GW_SUPPLEMENTARY,
}
GEN = "figures/fig_family_recovery.py"
PROVISIONAL = ("provisional: make.py data fig7 --recorded, scored from the recorded runs of registers 988-1019 and "
               "1100/1101 (2026-09-20/23), which include no de novo arm on a held-out human chromosome and no guided "
               "arm on a gorilla contig; make.py data fig7 builds every arm with the current binaries and driver "
               "defaults")

MODE_COLOR = {"denovo": figlib.BLUE[650], "guided": figlib.BLUE[300]}
CLASS_COLOR = {"exact": figlib.INK_2, "partial": "#aeada7", "missed": figlib.GRID}
# substrate status in d-f: marker SHAPE and grey level (all filled: open markers mean "primaries only" elsewhere)
STATUS_STYLE = {R.UNTOUCHED: dict(marker="o", color=figlib.INK, s=9),
                R.REUSED: dict(marker="D", color=figlib.INK_2, s=6),
                R.DEVELOPMENT: dict(marker="^", color=figlib.INK_3, s=10)}
LEGACY_STATUS = {"held out": R.UNTOUCHED}   # tables written before the three-way status (drawn as held out)
TRUTH_SHORT = {"compara": "Compara, primates", "referee": "protein homology", "soto": "Soto 2025",
               "npip_u2": "NPIP reference set"}
TRUTH_ORDER = ["compara", "soto", "npip_u2", "referee"]
# display labels of the dev tables' status values (GLOSSARY: name the decision)
STATUS_DISPLAY = {R.UNTOUCHED: "held out, never used", R.REUSED: "held out, reused verdict set",
                  R.DEVELOPMENT: "development"}
DODGE = 0.19


# ================================================================ data
def build(cfg: dict, data_dir: Path, force: bool = False, recorded: bool = False):
    """Genome scope (default): the genome-wide tables (build_genome; HEAVY units resumable, exit 75 = call again).
    Dev scope (`--set fig7_scope=dev`) or recorded=True: the per-chromosome tables (build_dev)."""
    if recorded or R.scope(cfg) == "dev":
        return build_dev(cfg, data_dir, force, recorded)
    return build_genome(cfg, data_dir, force)


def build_dev(cfg: dict, data_dir: Path, force: bool = False, recorded: bool = False):
    """Write fig7_summary and fig7_per_family (per-chromosome development tables).

    recorded=False (make.py data fig7): regenerate both modes on every substrate (HEAVY; _o1_recovery.ensure_sources,
    cached under ${work}/fig7/). recorded=True (make.py data fig7 --recorded): score the recorded runs (seconds)."""
    w = figlib.work_dir(cfg, R.FIG) / ("recorded" if recorded else "current")
    w.mkdir(parents=True, exist_ok=True)
    cached = str(cfg.get("fig7_dev_cached", "")).strip().lower() in ("1", "true", "yes")
    src = (R.recorded_sources(cfg, w) if recorded else R.cached_sources(cfg, w) if cached
           else R.ensure_sources(cfg, w, force=force))
    status = {(sp, c): st for sp, c, st in R.substrates(cfg)}
    notes = ([PROVISIONAL] if recorded else []) + [src["source"]]
    inputs, summary, per_family = dict(src.get("extra_inputs", {})), [], []
    for (sp, c, mode), (clusters, how) in sorted(src["arms"].items(), key=lambda kv: _order(cfg, kv[0])):
        genes = src["genes"][(sp, c)]
        inputs[f"clusters_{sp}_{c}_{mode}"] = clusters
        inputs[f"genes_{sp}_{c}"] = genes
        for (sp2, c2, truth), tpath in sorted(src["truths"].items()):
            if (sp2, c2) != (sp, c):
                continue
            inputs[f"truth_{sp}_{c}_{truth}"] = tpath
            label = f"{sp}_{c}_{mode}_{truth}"
            fs = R.family_score(cfg, clusters, genes, tpath, c, label, w / f"{label}.family_score.txt")
            mine = R.score_arm(clusters, str(genes), tpath, c)
            R.check_against_family_score(mine, fs, label)
            summary.append([sp, c, status[(sp, c)], truth, mode, how, mine["truth_families"], mine["truth_genes"],
                            mine["clusters_total"], mine["loci_total"], mine["clusters_scored"], mine["pair_tp"],
                            mine["truth_pairs"], mine["pred_pairs"], mine["pair_sens"], mine["pair_prec"],
                            mine["matched"], mine["truth_genes"], mine["pred_members"], mine["sens"], mine["prec"],
                            mine["f"], mine["exact"], mine["touched"], fs.get("collapsed"), fs.get("no_locus")])
            for r in mine["per_family"]:
                per_family.append([sp, c, status[(sp, c)], truth, mode, r["family"], r["label"], r["flagship"],
                                   r["n_truth"], r["cluster"], r["hit"], r["n_pred"], r["sens"], r["prec"], r["f"],
                                   r["jaccard"]])
    built = {(r[0], r[1]) for r in summary}
    pending = [f"{sp} {c}" for sp, c, _ in R.substrates(cfg) if (sp, c) not in built]
    missing = sorted({f"{sp} {c} {m}" for (sp, c) in built for m in R.MODES
                      if (sp, c, m) not in src["arms"]})
    status_notes = [f"{sp} {c}: {R.STATUS_NOTE.get((sp, c), st)}" for sp, c, st in R.substrates(cfg)]
    common = notes + status_notes + [
        "not built (no arm): " + (", ".join(pending) if pending else "none"),
        "mode missing on a built substrate: " + (", ".join(missing) if missing else "none"),
        "denovo = loci assembled from the IsoSeq reads (driver default: also built from secondary alignments "
        "within 2% of the read's best score) -> "
        "rustle_pipeline.sh families (mcl_families --from-gtf --min-exonic-bp 1 --min-shared-exon-frac 0.60); "
        "guided = annotated gene + pseudogene bodies -> minimap2 -x asm20 -c --eqx -P all-vs-all -> mcl_families "
        "--paf --gff --min-exonic-bp 1 --min-shared-exon-frac 0.60",
        "reference (column truth): compara = Ensembl Compara release 116 families at Primates (connected components "
        "of protein-coding paralogue pairs whose duplication node is at or below Primates; built genome-wide, "
        "restricted to the contig; the main figure's headline); soto = bench/soto/soto_famCN_S1C.tsv (first Family ID "
        "per gene; not independent of the 0.60 exon threshold); npip_u2 = the NPIP reference set of register 990 "
        "(built by bench/union_truth_npip.py, retired in 8db314c7); referee = the contig's protein-homology families "
        "(bench/truth.py protein-homology; SECONDARY, supplementary figure only); only reference families with >= 2 "
        "members on the contig are scored",
        "scored with family_score semantics: locus -> ONE gene by max overlap; predicted families intersected with "
        "the reference genes (an unlabelled member is not scored: precision is an upper bound, register 770/991)",
    ]
    if not recorded:
        common.append("guided all-vs-all: " + "; ".join(f"{sp} {c}: {v}" for (sp, c), v in src["paf_source"].items()))
    if "compara" not in {r[3] for r in summary}:
        common.append("compara: not scored (the genome-wide Compara table or the human annotation cache is absent)")
    figlib.write_table(
        "fig7_summary",
        ["species", "contig", "status", "truth", "mode", "source", "truth_families", "truth_genes", "clusters_total",
         "loci_in_clusters", "clusters_scored", "pair_tp", "truth_pairs", "pred_pairs", "pair_sens", "pair_prec",
         "bip_matched", "bip_truth_members", "bip_pred_members", "bip_recall", "bip_precision", "bip_f", "exact",
         "touched", "collapsed", "no_locus"],
        summary, generator=GEN, inputs=inputs, data_dir=data_dir,
        notes=common + [
            "pair_*: reference pairs (column truth_pairs) = within-family gene pairs of the reference; predicted pairs = "
            "within-family gene pairs over the reference genes; pair_sens = tp / truth_pairs, pair_prec = tp / pred_pairs",
            "bip_*: one-to-one assignment maximising total member overlap (scipy linear_sum_assignment = family_score);"
            " recall (sensitivity) = matched / reference members; precision = matched / members of the matched predicted "
            "families; asserted equal "
            "to family_score's printed sens / prec / F for every row",
            "exact = reference families with F = 1; touched = reference families whose matched family holds >= 1 member; "
            "collapsed / no_locus = family_score's columns (reference genes sharing a best-covering locus / with no locus "
            "in any cluster)",
        ])
    figlib.write_table(
        "fig7_per_family",
        ["species", "contig", "status", "truth", "mode", "family", "label", "flagship", "n_truth", "cluster", "hit",
         "n_pred", "sens", "prec", "f", "jaccard"],
        per_family, generator=GEN, inputs=inputs, data_dir=data_dir,
        notes=common + [
            "one row per reference family (with >= 2 members on the contig) per mode: its one-to-one matched predicted "
            "family, hit = members shared, n_pred = matched family size over the reference genes; sens = hit / n_truth, prec = "
            "hit / n_pred, f = harmonic mean, jaccard = hit / (n_truth + n_pred - hit); an unmatched family scores 0",
            f"label = longest common prefix of the named members (display only); flagship = first of {R.FLAGSHIP} "
            "that prefixes a member",
        ])


def _order(cfg, key):
    subs = [(sp, c) for sp, c, _ in R.substrates(cfg)]
    sp, c, mode = key
    return (subs.index((sp, c)) if (sp, c) in subs else 99, R.MODES.index(mode))


# ================================================================ plot
def _f(x):
    return None if x in ("", None) else float(x)


def _status_order(st):
    st = LEGACY_STATUS.get(st, st)
    return R.STATUSES.index(st) if st in R.STATUSES else len(R.STATUSES)


def _layout_rows(summary, truths=R.MAIN_DEV_TRUTHS):
    """[(kind, payload)] top to bottom: ('species', title), ('status', title) or ('row', (species, contig, status,
    truth)); statuses in the order untouched, reused, development (any other label after them); only the references
    in `truths` (main figure: Compara, Soto, NPIP; supplement: the protein-homology families)."""
    summary = [r for r in summary if r["truth"] in truths]
    seen, subs = set(), []
    for r in summary:
        k = (r["species"], r["contig"], r["status"])
        if k not in seen:
            seen.add(k)
            subs.append(k)
    truths = {}
    for r in summary:
        truths.setdefault((r["species"], r["contig"]), [])
        if r["truth"] not in truths[(r["species"], r["contig"])]:
            truths[(r["species"], r["contig"])].append(r["truth"])
    out = []
    for sp in ("human", "gorilla"):
        mine = [s for s in subs if s[0] == sp]
        if not mine:
            continue
        out.append(("species", {"human": "Human A119b (T2T-CHM13 v2.0)",
                                "gorilla": "Gorilla OR6737, testis (mGorGor1)"}.get(sp, sp)))
        statuses = sorted({s[2] for s in mine}, key=lambda st: (_status_order(st), [s[2] for s in mine].index(st)))
        for st in statuses:
            out.append(("status", st))
            for _, c, _ in (s for s in mine if s[2] == st):
                for t in sorted(truths[(sp, c)], key=TRUTH_ORDER.index):
                    out.append(("row", (sp, c, st, t)))
    return out


def _dot_column(ax, ys, by, key, title):
    for (y, k) in ys:
        vals = {m: _f(by[k + (m,)][key]) if k + (m,) in by else None for m in R.MODES}
        dn, g = vals["denovo"], vals["guided"]
        if dn is not None and g is not None:
            ax.plot([dn, g], [y - DODGE, y + DODGE], color=figlib.INK_3, linewidth=0.6, zorder=1)
        for m, v in vals.items():
            if v is not None:
                ax.plot([v], [y + (-DODGE if m == "denovo" else DODGE)], marker="o", markersize=3.4, linestyle="none",
                        markerfacecolor=MODE_COLOR[m], markeredgecolor=MODE_COLOR[m], markeredgewidth=0.4, zorder=3)
    ax.set_xlim(-0.04, 1.04)
    ax.set_xticks([0, 0.5, 1])
    ax.set_xticklabels(["0", "0.5", "1"])
    ax.grid(axis="x", color=figlib.GRID, linewidth=0.5)
    ax.grid(axis="y", visible=False)
    ax.tick_params(axis="y", left=False, labelleft=False)
    ax.tick_params(axis="x", labelsize=6, pad=1.5)
    ax.spines["left"].set_visible(False)
    ax.set_title(title, fontsize=6.5, pad=3)


def _class_column(ax, ys, by, perfam):
    for (y, k) in ys:
        for m in R.MODES:
            yy = y + (-DODGE if m == "denovo" else DODGE)
            if k + (m,) not in by:
                ax.text(0.02, yy, f"{R.MODE_LABEL[m]}: not built yet", fontsize=5.5, color=figlib.INK_3,
                        va="center", ha="left", style="italic")
                continue
            fams = perfam.get(k + (m,), [])
            n = len(fams)
            cls = {"exact": sum(1 for f in fams if f == 1.0), "partial": sum(1 for f in fams if 0 < f < 1),
                   "missed": sum(1 for f in fams if f == 0)}
            left = 0.0
            for c in ("exact", "partial", "missed"):
                wdt = cls[c] / n if n else 0
                if wdt:
                    ax.barh(yy, wdt, left=left, height=2 * DODGE * 0.86, color=CLASS_COLOR[c],
                            edgecolor=figlib.SURFACE, linewidth=0.4)
                    if c != "missed" and wdt >= 0.075:
                        ax.text(left + wdt / 2, yy, str(cls[c]), ha="center", va="center", fontsize=5.5,
                                color=figlib.SURFACE if c == "exact" else figlib.INK)
                left += wdt
            ax.text(1.02, yy, f"{cls['exact']}·{cls['partial']}·{cls['missed']}", ha="left", va="center",
                    fontsize=5.5, color=figlib.INK_2)
    ax.set_xlim(0, 1)
    ax.set_xticks([0, 0.5, 1])
    ax.set_xticklabels(["0", "0.5", "1"])
    ax.tick_params(axis="y", left=False, labelleft=False)
    ax.tick_params(axis="x", labelsize=6, pad=1.5)
    ax.spines["left"].set_visible(False)
    ax.grid(False)


def _natural(c):
    m = re.match(r"([A-Za-z_]*)(\d+)(.*)", c)
    return (m.group(1), int(m.group(2)), m.group(3)) if m else (c, 0, "")


def _seg_rect(p, q, x0, y0, x1, y1):
    """True if segment p-q meets the rectangle [x0, x1] x [y0, y1] (Liang-Barsky clipping)."""
    dx, dy = q[0] - p[0], q[1] - p[1]
    t0, t1 = 0.0, 1.0
    for pk, qk in ((-dx, p[0] - x0), (dx, x1 - p[0]), (-dy, p[1] - y0), (dy, y1 - p[1])):
        if pk == 0:
            if qk < 0:
                return False
            continue
        t = qk / pk
        if pk < 0:
            t0 = max(t0, t)
        else:
            t1 = min(t1, t)
        if t0 > t1:
            return False
    return True


def _seg_point(p, q, m):
    dx, dy = q[0] - p[0], q[1] - p[1]
    L2 = dx * dx + dy * dy
    t = 0.0 if L2 == 0 else max(0.0, min(1.0, ((m[0] - p[0]) * dx + (m[1] - p[1]) * dy) / L2))
    return ((p[0] + t * dx - m[0]) ** 2 + (p[1] + t * dy - m[1]) ** 2) ** 0.5


def _seg_seg(a, b, c, d):
    def orient(p, q, r):
        v = float((q[0] - p[0]) * (r[1] - p[1]) - (q[1] - p[1]) * (r[0] - p[0]))
        return 1 if v > 1e-9 else (-1 if v < -1e-9 else 0)
    return (orient(a, b, c) * orient(a, b, d) < 0) and (orient(c, d, a) * orient(c, d, b) < 0)


# label positions around a point, in points: (side, horizontal gap to the label's near edge, vertical offset of its
# centre); tried nearest first, left and right alike
LABEL_CANDIDATES = sorted(((side, gx, dy) for side in ("L", "R") for gx in (10, 17, 26, 38, 52)
                           for dy in (4, -5, 11, -12, 19, -20, 28, -29, 37, -38, 47, -48)),
                          key=lambda c: (c[1] ** 2 + c[2] ** 2, c[0] == "R"))


def _scatter(ax, rows, species_truth, fig):
    """Per-family F, de novo (x) vs guided (y), for the families at least one mode recovers; each family is one point,
    jittered by at most 0.015 (a fixed hash of its id) so that families with identical values stay visible; marker
    shape and grey level = substrate status."""
    import collections
    import hashlib

    pair = collections.defaultdict(dict)
    meta = {}
    for r in rows:
        key = (r["contig"], r["family"])
        pair[key][r["mode"]] = _f(r["f"])
        meta[key] = r
    both = {k: v for k, v in pair.items() if "denovo" in v and "guided" in v}
    contigs = sorted({k[0] for k in both}, key=_natural)
    ax.set_xlim(-0.06, 1.06)
    ax.set_ylim(-0.06, 1.06)
    for axis in (ax.xaxis, ax.yaxis):
        axis.set_ticks([0, 0.5, 1])
        axis.set_ticklabels(["0", "0.5", "1"])
    ax.set_aspect("equal")
    ax.grid(axis="both", color=figlib.GRID, linewidth=0.5)
    ax.set_xlabel("de novo, per-family F", labelpad=2)
    ax.set_ylabel("guided, per-family F", labelpad=2)
    if not both:
        ax.set_title(f"{species_truth}\nnot built: no substrate has both modes\n\n", loc="left", fontsize=6, pad=3,
                     linespacing=1.25)
        ax.text(0.5, 0.5, "not built yet", ha="center", va="center", fontsize=6, color=figlib.INK_3,
                style="italic")
        return
    ax.plot([0, 1], [0, 1], color=figlib.INK_3, linewidth=0.6, linestyle=(0, (2, 2)), zorder=1)
    drawn = {k: v for k, v in both.items() if v["denovo"] > 0 or v["guided"] > 0}
    none = len(both) - len(drawn)
    marks = []
    for k, v in sorted(drawn.items()):
        h = hashlib.md5(f"{k[0]}|{k[1]}".encode()).digest()
        jx, jy = (h[0] / 255 - 0.5) * 0.03, (h[1] / 255 - 0.5) * 0.03
        x, y = v["denovo"] + jx, v["guided"] + jy
        st = STATUS_STYLE.get(LEGACY_STATUS.get(meta[k]["status"], meta[k]["status"]), STATUS_STYLE[R.UNTOUCHED])
        ax.scatter([x], [y], s=st["s"], marker=st["marker"], color=st["color"], alpha=0.6, linewidths=0, zorder=3,
                   clip_on=False)
        marks.append((k, x, y))
    g_hi = sum(1 for v in drawn.values() if v["guided"] > v["denovo"])
    d_hi = sum(1 for v in drawn.values() if v["denovo"] > v["guided"])
    eq = len(drawn) - g_hi - d_hi
    ax.set_title(f"{species_truth}\n{', '.join(contigs)}\nguided higher {g_hi} · de novo higher {d_hi}\n"
                 f"equal {eq} · both 0 (not drawn) {none}", loc="left", fontsize=6, pad=3, linespacing=1.25)
    _label_flagships(ax, fig, marks, meta)


def _label_flagships(ax, fig, marks, meta, suffix=True):
    """Name the flagship families with leader lines from the middle of the label's near side. A position collides
    when its label covers a point or another label, or when its leader crosses (or runs within 3 pt of) another label,
    crosses another leader or passes over another point. Labels are placed greedily, in every order until one places
    them all without a collision; if none does, the order with the fewest collisions is kept (crossings weigh more
    than a covered point; a flagship label is never dropped silently)."""
    import itertools

    fig.canvas.draw()
    renderer = fig.canvas.get_renderer()
    px = fig.dpi / 72.0
    pad = 0.4 * 5.5 * px                     # the label's white box pad
    mark_px = [tuple(float(v) for v in ax.transData.transform((x, y))) for _, x, y in marks]
    rad = 1.8 * px + 1.0                     # marker radius (largest marker, s = 10 pt^2) + a little air
    air = 3 * px                             # a leader keeps this far from every other label
    cover = 1.0 * px                         # a point whose centre is this close to a label box is covered
    axbb = ax.get_window_extent(renderer)
    # a label may overhang the axes' right edge into the gap before the next panel (never past the figure)
    x_max = min(axbb.x1 + 14 * px, fig.get_window_extent(renderer).x1 - 2 * px)
    labels = []                              # (i, text, [(static cost, box, left, cy, start)])
    for i, (k, x, y) in enumerate(marks):
        fl = meta[k]["flagship"]
        if not fl:
            continue
        text = f"{meta[k]['label'] if meta[k]['label'].startswith(fl) else fl}" + (f" ({k[0]})" if suffix else "")
        probe = ax.text(0, 0, text, fontsize=5.5)
        tb = probe.get_window_extent(renderer)
        probe.remove()
        w, h = tb.width, tb.height
        mx, my = mark_px[i]
        others = [m for j, m in enumerate(mark_px) if j != i and (abs(m[0] - mx) > 0.5 or abs(m[1] - my) > 0.5)]
        cands = []
        fbb = fig.get_window_extent(renderer)
        for side, gx, dy in LABEL_CANDIDATES:
            left = mx + gx * px if side == "R" else mx - gx * px - w
            cy = my + dy * px
            box = (left - pad, cy - h / 2 - pad, left + w + pad, cy + h / 2 + pad)
            inside = box[0] >= axbb.x0 - 2 and box[2] <= x_max and box[1] >= axbb.y0 - 2 and box[3] <= axbb.y1 + 2
            if not inside:
                continue
            start = (box[2] if side == "L" else box[0], cy)          # middle of the side facing the point
            cost = sum(1 for m in others
                       if box[0] - cover <= m[0] <= box[2] + cover and box[1] - cover <= m[1] <= box[3] + cover)
            # a leader ends on its own point, so the points that touch it (a jittered pile) are not counted
            cost += sum(1 for m in others if _seg_point(start, (mx, my), m) < rad
                        and ((m[0] - mx) ** 2 + (m[1] - my) ** 2) ** 0.5 > 2 * rad)
            cands.append((cost, box, left, cy, start))
        if not cands:   # nothing inside the axes: allow the label anywhere inside the figure (said, never silent)
            for side, gx, dy in LABEL_CANDIDATES:
                left = mx + gx * px if side == "R" else mx - gx * px - w
                cy = my + dy * px
                box = (left - pad, cy - h / 2 - pad, left + w + pad, cy + h / 2 + pad)
                if box[0] >= fbb.x0 + 2 and box[2] <= fbb.x1 - 2 and box[1] >= fbb.y0 + 2 and box[3] <= fbb.y1 - 2:
                    start = (box[2] if side == "L" else box[0], cy)
                    cands.append((5, box, left, cy, start))
            if cands:
                print(f"[fig7] WARNING: the {text} label is placed outside its axes (no free position inside)",
                      file=sys.stderr)
        if not cands:
            raise RuntimeError(f"no position inside the figure for the {text} label")
        labels.append((i, text, cands))

    def place(order):
        total, chosen, placed, leaders = 0, {}, [], []
        for n in order:
            i, _, cands = labels[n]
            mxy = mark_px[i]
            best = None
            for c0, box, left, cy, start in cands:
                seg = (start, mxy)
                cost = c0
                cost += 3 * sum(1 for b in placed if box[0] < b[2] and b[0] < box[2] and box[1] < b[3] and b[1] < box[3])
                cost += 3 * sum(1 for b in placed if _seg_rect(*seg, b[0] - air, b[1] - air, b[2] + air, b[3] + air))
                cost += 3 * sum(1 for s2 in leaders if _seg_seg(*seg, *s2)
                                or _seg_rect(*s2, box[0] - air, box[1] - air, box[2] + air, box[3] + air))
                if best is None or cost < best[0]:
                    best = (cost, box, left, cy, start)
                if cost == 0:
                    break
            total += best[0]
            chosen[n] = best
            placed.append(best[1])
            leaders.append((best[4], mxy))
        return total, chosen

    best = None
    for order in itertools.islice(itertools.permutations(range(len(labels))), 720):
        total, chosen = place(order)
        if best is None or total < best[0]:
            best = (total, chosen)
        if total == 0:
            break
    inv = ax.transData.inverted()
    for n, (i, text, _) in enumerate(labels):
        _, box, left, cy, start = best[1][n]
        _, x, y = marks[i]
        ax.annotate(text, xy=(x, y), xytext=tuple(inv.transform((left, cy))), textcoords="data", fontsize=5.5,
                    color=figlib.INK, ha="left", va="center", zorder=5, annotation_clip=False,
                    bbox=dict(facecolor=figlib.SURFACE, edgecolor="none", pad=0.4, alpha=0.9))
        ax.annotate("", xy=(x, y), xytext=tuple(inv.transform(start)), textcoords="data", annotation_clip=False,
                    arrowprops=dict(arrowstyle="-", color=figlib.INK_2, linewidth=0.5, shrinkA=0, shrinkB=2))
    for n, (i, text, _) in enumerate(labels):   # never silent: a label that had to collide says so
        if best[1][n][0]:
            print(f"[fig7] WARNING: the {text} label collides (weight {best[1][n][0]}): no free position",
                  file=sys.stderr)
    return best[0]


# layout, in inches (the figure grows with the number of substrate x truth rows)
IN_PER_UNIT = 0.215      # one row of the dot plot
TOP_IN = 0.58            # key, panel headers, column titles
BOTTOM_IN = 2.78         # scatter block: headline, 4-line titles, 1.45 in axes, tick labels, x label, stamp


DEV_MAIN_SCATTERS = [("human", "compara", "d", "Human A119b · Compara families (primates)"),
                     ("human", "soto", "e", "Human A119b · Soto 2025 (not independent)"),
                     ("human", "npip_u2", "f", "Human A119b · NPIP set (chr16)")]
DEV_SUPP_SCATTERS = [("human", "referee", "d", "Human A119b · protein-homology families"),
                     ("gorilla", "referee", "e", "Gorilla OR6737 · protein-homology families")]


def _plot_dev(data_dir: Path, out_dir: Path, truths=R.MAIN_DEV_TRUTHS, name="fig7_family_recovery",
              scatters=DEV_MAIN_SCATTERS, banner=None):
    import matplotlib.pyplot as plt
    import matplotlib.transforms as mtransforms
    from matplotlib.patches import Patch

    summary = [r for r in figlib.read_table("fig7_summary", data_dir) if r["truth"] in truths]
    perfam_rows = [r for r in figlib.read_table("fig7_per_family", data_dir) if r["truth"] in truths]
    if not summary:
        return []
    by = {(r["species"], r["contig"], r["status"], r["truth"], r["mode"]): r for r in summary}
    perfam = {}
    for r in perfam_rows:
        perfam.setdefault((r["species"], r["contig"], r["status"], r["truth"], r["mode"]), []).append(float(r["f"]))
    ys, y, labels, groups = [], 0.0, [], []   # groups: (y, kind, title, separator y or None)
    prev = None
    for kind, payload in _layout_rows(summary, truths):
        if kind == "species":
            y += 0.35 if groups else 0.0
            groups.append((y, kind, payload, y - 0.5 if groups else None))
            y += 0.8
        elif kind == "status":
            if prev != "species":
                y += 0.2
            groups.append((y, kind, payload, y - 0.5 if prev not in ("species", None) else None))
            y += 0.8
        else:
            sp, c, st, t = payload
            r0 = next(by[payload + (m,)] for m in R.MODES if payload + (m,) in by)
            labels.append((y, f"{c} · {TRUTH_SHORT[t]}",
                           f"{int(r0['truth_families'])} families · {int(r0['truth_genes'])} genes"))
            ys.append((y, payload))
            y += 1.0
        prev = kind
    y_lo, y_hi = -0.75, y - 0.4
    dot_in = (y_hi - y_lo) * IN_PER_UNIT
    H = TOP_IN + dot_in + 0.30 + BOTTOM_IN
    W = figlib.WIDTH_DOUBLE
    fy = lambda inch_from_top: 1 - inch_from_top / H           # noqa: E731

    fig = plt.figure(figsize=(W, H))
    top = fig.add_gridspec(1, 7, left=0.215, right=0.935, top=fy(TOP_IN), bottom=fy(TOP_IN + dot_in), wspace=0.13,
                           width_ratios=[1, 1, 0.35, 1, 1, 1, 1.9])
    axes = [fig.add_subplot(top[0, i]) for i in (0, 1, 3, 4, 5, 6)]
    for ax in axes[1:]:
        ax.sharey(axes[0])
    metrics = [("pair_sens", "Sensitivity"), ("pair_prec", "Precision"), ("bip_recall", "Sensitivity"),
               ("bip_precision", "Precision"), ("bip_f", "F")]
    for ax, (key, title) in zip(axes[:5], metrics):
        _dot_column(ax, ys, by, key, title)
    _class_column(axes[5], ys, by, perfam)
    axes[0].set_ylim(y_hi, y_lo)
    # row labels (two lines: what, and its n) and group headers in the left margin
    tr = mtransforms.blended_transform_factory(fig.transFigure, axes[0].transData)
    for yy, text, n in labels:
        fig.text(0.207, yy - 0.2, text, transform=tr, ha="right", va="center", fontsize=6, color=figlib.INK)
        fig.text(0.207, yy + 0.25, n, transform=tr, ha="right", va="center", fontsize=5.5, color=figlib.INK_2)
    for yy, kind, text, sep in groups:
        if kind == "species":
            fig.text(0.012, yy, text, transform=tr, ha="left", va="center", fontsize=6.5, fontweight="bold",
                     color=figlib.INK)
        else:
            fig.text(0.024, yy, STATUS_DISPLAY.get(LEGACY_STATUS.get(text, text), text), transform=tr, ha="left",
                     va="center", fontsize=6.2, color=figlib.INK_2)
        if sep is not None:
            for ax in axes:
                ax.axhline(sep, color=figlib.INK_3 if kind == "species" else figlib.GRID, linewidth=0.5, zorder=0)

    def header(ax0, ax1, letter, text):
        x0, x1 = ax0.get_position().x0, ax1.get_position().x1
        fig.text(x0, fy(0.24), letter, fontsize=9, fontweight="bold", ha="left", va="bottom")
        fig.text(x0 + 0.018, fy(0.235), text, fontsize=6.5, ha="left", va="bottom", color=figlib.INK)
        fig.add_artist(plt.Line2D([x0, x1], [fy(0.27), fy(0.27)], transform=fig.transFigure, color=figlib.INK_3,
                                  linewidth=0.5))
    header(axes[0], axes[1], "a", "Pairwise (gene pairs)")
    header(axes[2], axes[4], "b", "One-to-one bipartite (family members)")
    header(axes[5], axes[5], "c", "Reference families, per mode")
    # mode key: colour AND position (upper = de novo, lower = guided), once, in the label column
    for i, m in enumerate(R.MODES):
        yk = fy(0.1 + 0.15 * i)
        fig.add_artist(plt.Line2D([0.016], [yk], transform=fig.transFigure, marker="o", markersize=3.4,
                                  markerfacecolor=MODE_COLOR[m], markeredgecolor=MODE_COLOR[m], linestyle="none"))
        fig.text(0.024, yk, f"{R.MODE_LABEL[m]} ({'upper' if m == 'denovo' else 'lower'} dot and bar)",
                 fontsize=6, va="center", ha="left", color=figlib.INK)
    fig.text(0.012, fy(0.4), "same family rule in both modes; Rustle-internal, not a tool comparison", fontsize=5.5,
             va="center", ha="left", color=figlib.INK_3)
    axes[5].legend(handles=[Patch(facecolor=CLASS_COLOR[c], edgecolor=figlib.INK_3 if c == "missed" else
                                  CLASS_COLOR[c], linewidth=0.4, label=c) for c in ("exact", "partial", "missed")],
                   loc="lower center", bbox_to_anchor=(0.5, 1.0), ncol=3, fontsize=6, handlelength=0.9,
                   handleheight=0.8, columnspacing=0.9, handletextpad=0.35, borderaxespad=0.2)

    # bottom block: per-family scatters
    b_top = TOP_IN + dot_in + 0.30
    head_y = fy(b_top + 0.12)
    sc_top, sc_bot = b_top + 0.2 + 0.62, b_top + 0.2 + 0.62 + 1.45
    bot = fig.add_gridspec(1, 3, left=0.075, right=0.975, top=fy(sc_top), bottom=fy(sc_bot), wspace=0.55)
    sc = [fig.add_subplot(bot[0, i]) for i in range(len(scatters))]
    for ax, (sp, t, letter, title) in zip(sc, scatters):
        _scatter(ax, [r for r in perfam_rows if r["species"] == sp and r["truth"] == t], title, fig)
        fig.text(ax.get_position().x0 - 0.075, fy(sc_top - 0.43), letter, fontsize=9, fontweight="bold", ha="left",
                 va="bottom")
    head = fig.text(0.075, head_y, "Per-family F: one point per reference family that at least one mode recovers "
                    "(jitter ≤ 0.015)", fontsize=6, va="center", ha="left", color=figlib.INK)
    # status key, right-aligned on the same line: shape AND grey level
    renderer = fig.canvas.get_renderer()
    fw = fig.get_window_extent(renderer).width
    xr = 0.975
    for st in reversed(R.STATUSES):
        t = fig.text(xr, head_y, STATUS_DISPLAY.get(st, st), fontsize=5.8, va="center", ha="right")
        xm = xr - t.get_window_extent(renderer).width / fw - 0.008
        sty = STATUS_STYLE[st]
        fig.add_artist(plt.Line2D([xm], [head_y], transform=fig.transFigure, marker=sty["marker"],
                                  markersize=sty["s"] ** 0.5 * 1.15, markerfacecolor=sty["color"],
                                  markeredgecolor=sty["color"], linestyle="none", alpha=0.75))
        xr = xm - 0.022
    if head.get_window_extent(renderer).x1 / fw > xr + 0.01:
        raise RuntimeError("fig7: the status key runs into the scatter heading")
    if banner:
        fig.text(0.012, fy(H - 0.05), banner + " Development tables (per-chromosome runs).", fontsize=5.5,
                 color=figlib.INK_2, ha="left", va="bottom")
    else:
        fig.text(0.975, fy(H - 0.05), "Development tables (per-chromosome runs); the genome-wide version replaces them",
                 fontsize=5.5, color=figlib.INK_3, ha="right", va="bottom")
    figlib.stamp_provisional(fig, DEV_TABLES, data_dir)
    paths = figlib.save(fig, name, out_dir)
    plt.close(fig)
    return paths


# ================================================================ genome-wide build
GW_UNIT_VERSION = "1"   # bump when a cached scoring unit's rows change
SUMMARY_COLS = ["arm", "mode", "sample", "sample_label", "species", "reference", "reference_label", "substrate",
                "truth_families", "truth_genes", "clusters_total", "loci_in_clusters", "clusters_scored", "bip_matched",
                "bip_pred_members", "bip_sens", "bip_prec", "bip_f", "pair_truth", "pair_pred", "pair_tp", "pair_sens",
                "pair_prec", "exact", "partial", "missed", "collapsed", "no_locus", "pair_covered", "status"]
CONTIG_COLS = ["arm", "mode", "sample", "sample_label", "species", "reference", "contig", "chromosome", "exposure",
               "truth_families", "truth_genes", "bip_sens", "bip_prec", "bip_f", "pair_sens", "pair_prec", "exact",
               "partial", "missed"]
PERFAM_COLS = ["arm", "mode", "sample", "species", "reference", "substrate", "family_id", "label", "flagship",
               "n_truth", "cluster", "hit", "n_pred", "sens", "prec", "f", "n_contigs"]


def _unit(w: Path, name: str, key, fn, budget, est_s: float):
    """A cached scoring unit (rows), recomputed only when `key` (input fingerprints) changes."""
    import json
    p = w / f"{name}.json"
    k = json.dumps([GW_UNIT_VERSION, key], sort_keys=True)
    if p.exists():
        j = json.loads(p.read_text())
        if j.get("key") == k:
            return j["value"]
    R._o1.need(budget, est_s, name)
    v = fn()
    p.write_text(json.dumps({"key": k, "value": v}))
    return v


def _names(members: str) -> list:
    return [m.split(":", 1)[1] if ":" in m else m for m in members.split(",") if m and m != "-"]


def _arm_rows(cfg, arm, ref, tpath, w) -> dict:
    """Every scoring of one arm (de novo sample or guided species) against one reference: S0 / S1 / S2 (S1 with
    per-family rows), the per-contig breakdown, or the chr16 inset for the NPIP set."""
    o1 = R._o1
    sp = arm["species"]
    genes_gff = o1.annotation_cache(cfg, sp)["genes_gff"]
    n_cl, n_loci = R.clusters_profile(Path(arm["clusters"]))
    wd = w / arm["arm"].replace(":", "_") / ref
    base = [arm["arm"], arm["mode"], arm["sample"], arm["label"], sp, ref, R.GW_TRUTH_LABEL[ref]]

    def srow(sub, r):
        return base + [sub, r["truth_families"], r["truth_genes"], n_cl, n_loci, r["clusters_scored"], r["matched"],
                       r["pred_members"], r["sens"], r["prec"], r["f"], r["pair_truth"], r["pair_pred"], r["pair_tp"],
                       r["pair_sens"], r["pair_prec"], r["exact"], r["partial"], r["missed"], r["collapsed"],
                       r["no_locus"], "", "ok"]
    out = {"summary": [], "contigs": [], "per_family": []}
    if ref == "npip_u2":
        r = R.fs_score(cfg, arm["clusters"], genes_gff, tpath, f"{arm['arm']}_{ref}", wd, keep={"chr16"})
        out["summary"].append(srow("chr16_inset", r))
        subs = [("chr16_inset", r)]
    else:
        subs = []
        for sub in ("S0", "S1", "S2"):
            r = R.fs_score(cfg, arm["clusters"], genes_gff, tpath, f"{arm['arm']}_{ref}_{sub}", wd,
                           drop=o1.substrate_drop(sp, sub))
            out["summary"].append(srow(sub, r))
            if sub == "S1":
                subs.append((sub, r))
        names = o1.contig_names(cfg, sp)
        by_contig = _split_truth_contigs(tpath)
        for c in R.annotation_contigs(cfg, sp):
            if len(by_contig.get(c, ())) < 2:
                continue   # no reference family can have >= 2 members here
            r = R.fs_score(cfg, arm["clusters"], genes_gff, tpath, f"{arm['arm']}_{ref}_{c}", wd / "contigs",
                           keep={c})
            if not r["truth_families"]:
                continue
            out["contigs"].append([arm["arm"], arm["mode"], arm["sample"], arm["label"], sp, ref, c, names.get(c, c),
                                   o1.exposure(sp, c), r["truth_families"], r["truth_genes"], r["sens"], r["prec"],
                                   r["f"], r["pair_sens"], r["pair_prec"], r["exact"], r["partial"], r["missed"]])
    for sub, r in subs:
        for p in r["per_family"]:
            nm = _names(p["members"])
            out["per_family"].append([arm["arm"], arm["mode"], arm["sample"], sp, ref, sub, p["family_id"],
                                      R.family_label(nm), R.flagship_of(nm), int(p["n_truth"]),
                                      "" if p["cluster"] == "-" else p["cluster"], int(p["hit"]), int(p["n_pred"]),
                                      float(p["sens"]), float(p["prec"]), float(p["f"]), int(p["n_contigs"])])
    return out


_TRUTH_CONTIGS: dict = {}


def _split_truth_contigs(tpath) -> dict:
    """{contig: set of gene names} of a contig-keyed reference table (memoised per file)."""
    key = str(tpath)
    if key not in _TRUTH_CONTIGS:
        d: dict = {}
        with open(tpath) as fh:
            head = fh.readline().rstrip("\n").split("\t")
            ci = head.index("Contig")
            for ln in fh:
                f = ln.rstrip("\n").split("\t")
                d.setdefault(f[ci], set()).add(f[0])
        _TRUTH_CONTIGS[key] = d
    return _TRUTH_CONTIGS[key]


def _liftoff_rows(cfg, arm, n_cl, n_loci) -> list:
    """Liftoff copy pairs (prereg Amendment 1, item 4): the Fig. 8 self-lift's (record, extra copy) pairs, extra copy
    sequence_ID >= 0.95, both exon unions >= 200 bp, both loci read-supported in the sample; recovered when a copy of
    one DE NOVO family covers >= 50% of each locus's exon union (Liftoff's -a). Pairwise sensitivity only (the relation
    certifies pairs; it is not a partition). S0 / S1 / S2 of the de novo arm (the guided arm is not scored: its loci
    are the annotation, and an extra copy is unannotated by construction)."""
    import _liftoff as L
    o1 = R._o1
    sp, sid = arm["species"], arm["sample"]
    base = [arm["arm"], arm["mode"], sid, arm["label"], sp, "liftoff", R.GW_TRUTH_LABEL["liftoff"]]
    loci, sup, why = L.support_if_ready(cfg, sid, sp)
    if loci is None or not arm.get("copies"):
        why = why or arm.get("copies_missing", "families copy table absent")
        return [base + ["S1"] + [""] * 21 + [why]]
    fam_loci = L.catalog_loci(Path(arm["copies"]))
    out = []
    for sub in ("S0", "S1", "S2"):
        drop = o1.substrate_drop(sp, sub)
        res = L.pair_families(L.copy_pairs(loci, 0.95, sup, (lambda c: c not in drop) if drop else None), fam_loci)
        k, n = sum(1 for _, _, shared in res if shared), len(res)
        out.append(base + [sub, "", "", n_cl, n_loci, "", "", "", "", "", "", n, "", k, k / n if n else None, "",
                           "", "", "", "", "", sum(1 for _, both, _ in res if both), "ok"])
    return out


BRIDGE_CONTIGS = (("human", "chr16"), ("gorilla", "NC_073244.2"))   # the two development contigs (PREREG 2.2)


def _bridge_rows(cfg, budget, missing) -> list:
    """The pre-registered bridging measurement: per-chromosome guided families on the development contigs with the
    recorded flags (-x asm20 -c --eqx -P) vs the de novo flags, scored against the contig's own protein-homology
    families and (human) Soto 2025, per-chromosome family_score semantics. Reported; decides nothing."""
    wdir = figlib.work_dir(cfg, R.FIG) / "current"          # shares the dev tables' caches (gff slices, references)
    gdir = wdir / "gff"
    gdir.mkdir(parents=True, exist_ok=True)
    rows = []
    for sp, contig in BRIDGE_CONTIGS:
        try:
            full = R._o1.annotation_gff(cfg, sp)
            fasta_key = f"{sp}_fasta"
            if fasta_key not in cfg:
                raise R._o1.NotBuilt(f"inputs key {fasta_key} absent")
        except R._o1.NotBuilt as e:
            missing.append(f"bridge {sp}: {e}")
            continue
        genes = R.gff_slices(full, {contig: gdir / f"{sp}_{contig}.gff"})[contig]
        R._o1.need(budget, 200, f"bridge {sp} {contig}")
        truths = {"referee": R.referee_build(cfg, sp, contig, genes, wdir)}
        if sp == "human":
            truths["soto"] = R.bench(cfg) / "soto" / "soto_famCN_S1C.tsv"
        for flags in ("recipe", "denovo"):
            R._o1.need(budget, 200, f"bridge {sp} {contig} {flags}")
            clusters, how = R.guided_families(cfg, sp, contig, genes, full, wdir, flags=flags)
            n_cl, n_loci = R.clusters_profile(clusters)
            for ref, tpath in truths.items():
                lab = f"bridge_{sp}_{contig}_{flags}_{ref}"
                fs = R.family_score(cfg, clusters, genes, tpath, contig, lab, wdir / f"{lab}.family_score.txt")
                mine = R.score_arm(clusters, str(genes), tpath, contig)
                R.check_against_family_score(mine, fs, lab)
                rows.append([sp, contig, "protein_homology" if ref == "referee" else ref, flags,
                             R.GUIDED_MM2[flags], how, mine["truth_families"], mine["truth_genes"], mine["sens"],
                             mine["prec"], mine["f"], mine["pair_sens"], mine["pair_prec"], n_cl, n_loci])
    return rows


def build_genome(cfg: dict, data_dir: Path, force: bool = False):
    o1 = R._o1
    budget = o1.gw_budget(cfg, "fig7")
    w = figlib.work_dir(cfg, R.FIG) / "gw"
    w.mkdir(parents=True, exist_ok=True)
    by_sp = o1.samples_by_species(cfg)
    missing: list = []
    fs_bin = Path(cfg["bin"]) / "family_score"

    def miss(e):
        if str(e) not in missing:
            missing.append(str(e))
            print(f"[fig7] not built: {e}", file=sys.stderr)

    # ---- arms: de novo per sample (run cache only), guided per species (HEAVY, resumable)
    arms = []
    for sp, sids in by_sp.items():
        for sid in sids:
            try:
                arm = {"arm": f"denovo:{sid}", "mode": "denovo", "sample": sid, "label": o1.sample_label(cfg, sid),
                       "species": sp, "clusters": str(o1.stage_product(cfg, sid, "families", "clusters"))}
            except o1.NotBuilt as e:
                miss(e)
                continue
            try:
                arm["copies"] = str(o1.families_copies(cfg, sid, "families"))   # the default families' copy table
            except o1.NotBuilt as e:
                arm["copies_missing"] = str(e)
                miss(e)
            arms.append(arm)
        try:
            arms.append({"arm": f"guided:{sp}", "mode": "guided", "sample": "", "label": f"{sp} annotation (guided)",
                         "species": sp, "clusters": str(R.guided_gw(cfg, sp, budget))})
        except o1.NotBuilt as e:
            miss(e)
    # ---- references
    truths = {}
    ph = str(cfg.get("fig7_protein_homology", "")).strip().lower() in ("1", "true", "yes")
    for sp in by_sp:
        for t in R.gw_truths(sp, ph) + (["npip_u2"] if sp == "human" else []):
            try:
                p = R.gw_truth_table(cfg, sp, t)
                if p:
                    truths[(sp, t)] = p
            except o1.NotBuilt as e:
                miss(e)
    bridge_key = [figlib.file_fingerprint(p) for p in (fs_bin, Path(cfg["bin"]) / "mcl_families")]
    for sp in ("human", "gorilla"):
        try:
            bridge_key.append(figlib.file_fingerprint(o1.annotation_gff(cfg, sp)))
        except o1.NotBuilt:
            pass
    bridge = _unit(w, "bridge", bridge_key,
                   lambda: _bridge_rows(cfg, budget, missing), budget, 400)
    if not arms:
        raise o1.NotBuilt("fig7 genome scope: nothing to score yet:\n  " + "\n  ".join(missing))

    # ---- scoring units (cached per arm x reference)
    summary, contigs, perfam, levels, inputs = [], [], [], [], {"family_score": fs_bin}
    for arm in arms:
        inputs[f"clusters_{arm['arm']}"] = arm["clusters"]
        if arm["mode"] == "denovo":
            import _liftoff as L
            n_cl, n_loci = R.clusters_profile(Path(arm["clusters"]))
            lp = L.loci_path(cfg, arm["species"])
            key = [figlib.file_fingerprint(p) for p in (arm.get("copies"), lp, L.support_path(cfg, arm["sample"]))
                   if p and Path(p).exists()]
            lrows = _unit(w, f"liftoff_{arm['arm'].replace(':', '_')}", key,
                          lambda arm=arm, n_cl=n_cl, n_loci=n_loci: _liftoff_rows(cfg, arm, n_cl, n_loci), budget, 60)
            if any(r[-1] != "ok" for r in lrows):
                miss(next(r[-1] for r in lrows if r[-1] != "ok"))
            else:
                inputs[f"copies_{arm['arm']}"] = arm["copies"]
                inputs[f"liftoff_loci_{arm['species']}"] = lp
                inputs[f"liftoff_support_{arm['sample']}"] = L.support_path(cfg, arm["sample"])
            summary += lrows
        for (sp, ref), tpath in truths.items():
            if sp != arm["species"]:
                continue
            inputs[f"reference_{sp}_{ref}"] = tpath
            genes_gff = o1.annotation_cache(cfg, sp)["genes_gff"]
            rows = _unit(w, f"score_{arm['arm'].replace(':', '_')}_{ref}",
                         [figlib.file_fingerprint(p) for p in (arm["clusters"], tpath, genes_gff, fs_bin)],
                         lambda arm=arm, ref=ref, tpath=tpath: _arm_rows(cfg, arm, ref, tpath, w), budget, 150)
            summary += rows["summary"]
            contigs += rows["contigs"]
            perfam += rows["per_family"]
        if arm["species"] == "human" and ("human", "compara") in truths:
            def sweep(arm=arm):
                out = []
                for level in R.COMPARA_SWEEP:
                    tp, st = R.compara_families(cfg, level)
                    r = R.fs_score(cfg, arm["clusters"], o1.annotation_cache(cfg, "human")["genes_gff"], tp,
                                   f"{arm['arm']}_compara_{level}", w / arm["arm"].replace(":", "_") / "levels",
                                   drop=o1.substrate_drop("human", "S1"))
                    out.append([arm["arm"], arm["mode"], arm["sample"], arm["label"], level, "S1", st["families"],
                                st["genes"], r["truth_families"], r["truth_genes"], r["sens"], r["prec"], r["f"],
                                r["pair_sens"], r["pair_prec"], r["exact"], r["partial"], r["missed"]])
                return out
            levels += _unit(w, f"levels_{arm['arm'].replace(':', '_')}",
                            [figlib.file_fingerprint(p) for p in (arm["clusters"], o1.compara_gw(cfg), fs_bin)],
                            sweep, budget, 120)
    # per-family rows: keep the families that at least one arm of their species recovers (F > 0) — the scatters;
    # the exact / partial / missed counts of every family are in fig7_gw_summary
    hit = {(r[3], r[4], r[6]) for r in perfam if r[15] > 0}
    if not summary:
        raise o1.NotBuilt("fig7 genome scope: nothing scored yet:\n  " + "\n  ".join(missing))
    perfam = [r for r in perfam if (r[3], r[4], r[6]) in hit]

    notes = ["genome-wide, pre-registered in docs/PREREG_genome_wide_families_2026-09-25.md; a Rustle-internal "
             "comparison of its two modes (de novo, guided), not a comparison with other tools"]
    if missing:
        notes.append("provisional: not every sample, mode or reference is built yet — " + " | ".join(missing))
    common = notes + [
        "de novo (arm denovo:<sample>) = the run-cache families stage of the sample (tools/rustle_pipeline.sh families "
        "= mcl_families --from-gtf --min-exonic-bp 1 --min-shared-exon-frac 0.60 on the genome-wide default assembly; "
        "minimap2 -x asm20 -c -X -N 50 -p 0.1 --secondary=yes all-vs-all through tools/mm2_shard.sh); guided (arm "
        "guided:<species>) = every annotated gene and pseudogene body of the species -> the same all-vs-all flags "
        f"(fig7_guided_flags = {cfg.get('fig7_guided_flags', 'denovo')}) -> mcl_families --paf --gff "
        "--min-exonic-bp 1 --min-shared-exon-frac 0.60; guided reads no RNA (one run per species)",
        "references (external; prereg Amendment 1): Compara = Ensembl Compara release 116, connected components of "
        "protein-coding paralogue pairs whose duplication node is at or below Primates (the headline); Soto = Soto et "
        "al. 2025 S1C, first family per gene (NOT independent: the 0.60 exon threshold was chosen against it); npip_u2 "
        "= the NPIP reference set (register 990), human chr16 only (an inset); liftoff = Liftoff's self-lift (Fig. 8) "
        "(record, extra copy) pairs, extra copy sequence_ID >= 0.95, both exon unions >= 200 bp, both loci supported "
        "by >= 2 reads of the sample (primary alignment, aligned block on the exon union), recovered when one de novo "
        "family's copies (the families stage copy table, representative exons) cover >= 50% of each locus's exon "
        "union; de novo only (an extra copy is unannotated by construction, so the guided node set cannot hold it); "
        "pair_truth = pairs, pair_tp = recovered, pair_covered = pairs whose two loci are both covered; "
        "protein_homology = SECONDARY (supplementary figure only; bench/truth.py protein-homology --chrom ALL; only "
        "with fig7_protein_homology 1)",
        "substrates: S0 whole genome; S1 genome minus development contigs (human chr16, gorilla NC_073244.2; the "
        "headline); S2 genome minus every contig used for a family decision (human chr16, chr5, chr7, chr21, chr2, "
        "chr8, chr10; gorilla NC_073244.2); a restriction keeps each reference family's members and each predicted "
        "family's loci on the kept contigs (a family with < 2 members left is not scored)",
        "scored by family_score --chrom ALL --per-family --pairwise: a locus is labelled with ONE gene (max overlap); "
        "predicted families are intersected with the reference genes, so precision is an upper bound for both modes",
    ]
    figlib.write_table("fig7_gw_summary", SUMMARY_COLS, summary, generator=GEN, inputs=inputs, data_dir=data_dir,
                       notes=common + [
                           "bip_*: one-to-one assignment maximising shared members; bip_sens = matched / reference "
                           "members; bip_prec = matched / members of the matched predicted families; bip_f = harmonic "
                           "mean (pooled, not a mean of per-family F); recomputed exactly from the per-family rows and "
                           "checked against family_score's printed line",
                           "pair_*: within-family gene pairs of the reference (pair_truth), within-family gene pairs "
                           "of the predicted families over the reference genes (pair_pred), shared (pair_tp)",
                           "exact / partial / missed = reference families with F = 1 / 0 < F < 1 / no member in the "
                           "matched predicted family; collapsed / no_locus = family_score's columns"])
    figlib.write_table("fig7_gw_contigs", CONTIG_COLS, contigs, generator=GEN, inputs=inputs, data_dir=data_dir,
                       notes=common + [
                           "the GENOME-WIDE run restricted to one contig (not a per-chromosome run: a family built "
                           "genome-wide can have lost members on other contigs); every contig with >= 2 reference genes",
                           "exposure: development / threshold selection (chr5, chr7, chr21: the 0.60 exon threshold, "
                           "register 903) / reused verdict set (chr2, chr8, chr10: about 30 guided-mode tests) / scored "
                           "once, no decision (human chr6, gorilla NC_073234.2: the per-chromosome Fig. 7) / never used "
                           "for a family decision"])
    figlib.write_table("fig7_gw_per_family", PERFAM_COLS, perfam, generator=GEN, inputs=inputs, data_dir=data_dir,
                       notes=common + [
                           "substrate S1 (and the chr16 inset of the NPIP set); one row per reference family per arm, "
                           "kept only for families that at least one arm of their species recovers (F > 0): the "
                           "scatters; every family's class is counted in fig7_gw_summary",
                           "label = longest common prefix of the named members (display only); flagship = first of "
                           f"{R.FLAGSHIP} that prefixes a member"])
    if levels:
        figlib.write_table(
            "fig7_gw_compara_levels",
            ["arm", "mode", "sample", "sample_label", "level", "substrate", "families_genome", "genes_genome",
             "truth_families", "truth_genes", "bip_sens", "bip_prec", "bip_f", "pair_sens", "pair_prec", "exact",
             "partial", "missed"], levels, generator=GEN, inputs=inputs, data_dir=data_dir,
            notes=common + ["Compara families at cumulative duplication levels (a family = connected component of the "
                            "protein-coding pairs whose duplication node is at or below the level); Primates is the "
                            "pre-registered headline, the others are a sweep (table only)"])
    if bridge:
        figlib.write_table(
            "fig7_gw_bridge",
            ["species", "contig", "reference", "flags", "minimap2_flags", "paf_source", "truth_families", "truth_genes",
             "bip_sens", "bip_prec", "bip_f", "pair_sens", "pair_prec", "clusters_total", "loci_in_clusters"], bridge,
            generator=GEN, inputs=inputs, data_dir=data_dir,
            notes=notes + ["the pre-registered bridging measurement (reported, decides nothing): per-chromosome guided "
                           "families on the two development contigs with the recorded flags (recipe) and with the de "
                           "novo mode's flags, scored against the contig's own protein-homology families (bench/truth.py "
                           "protein-homology --chrom C) and Soto 2025 with per-chromosome family_score semantics"])


# ================================================================ genome-wide plot
EXPO_STYLE = {R._o1.DEVELOPMENT: dict(marker="^", s=12), R._o1.THRESHOLD: dict(marker="D", s=7),
              R._o1.REUSED: dict(marker="D", s=7), R._o1.SCORED_ONCE: dict(marker="o", s=9, open=True),
              R._o1.NEVER: dict(marker="o", s=8)}
EXPO_LEGEND = [(R._o1.NEVER, "never used for a family decision"), (R._o1.SCORED_ONCE, "scored once, no decision"),
               (R._o1.REUSED, "threshold selection or reused verdict set"), (R._o1.DEVELOPMENT, "development")]


GENOME_LABEL = {"human": "Human (T2T-CHM13 v2.0)", "gorilla": "Gorilla (mGorGor1)", "chimpanzee": "Chimpanzee (mPanTro3)",
                "orangutan": "Orangutan (mPonPyg2)"}


def _short(sid: str) -> str:
    """'human_A119b' -> 'A119b', 'gorilla_KB3781' -> 'KB3781' (the species is the row group's title)."""
    return sid.split("_", 1)[1] if "_" in sid else sid


def _gw_rows(summary, keep_refs=MAIN_REFS):
    """Dot-table rows top to bottom: [('species', title) | ('row', (species, sample, reference))], references in
    `keep_refs` only."""
    samples_by_sp, refs = {}, {}
    for r in summary:
        if r["reference"] not in keep_refs:
            continue
        if r["mode"] == "denovo":
            samples_by_sp.setdefault(r["species"], [])
            if r["sample"] not in samples_by_sp[r["species"]]:
                samples_by_sp[r["species"]].append(r["sample"])
            refs.setdefault(r["species"], [])
            if r["reference"] not in refs[r["species"]]:
                refs[r["species"]].append(r["reference"])
    order = ["compara", "soto", "liftoff", "protein_homology", "npip_u2"]
    out = []
    for sp in [s for s in R._o1.SPECIES_ORDER if s in samples_by_sp]:
        out.append(("species", GENOME_LABEL.get(sp, sp.capitalize())))
        rs = sorted(refs[sp], key=order.index)
        for sid in samples_by_sp[sp]:
            for ref in rs:
                if ref != "npip_u2":
                    out.append(("row", (sp, sid, ref)))
        if "npip_u2" in rs:
            out.append(("status", "Inset: NPIP reference set, chr16 (development chromosome)"))
            for sid in samples_by_sp[sp]:
                out.append(("row", (sp, sid, "npip_u2")))
    return out


def _gw_dot(ax, ys, val, key, title, ring=None):
    for y, k in ys:
        dn, g = val(k, "denovo", key), val(k, "guided", key)
        if dn is not None and g is not None:
            ax.plot([dn, g], [y - DODGE, y + DODGE], color=figlib.INK_3, linewidth=0.6, zorder=1)
        for m, v in (("denovo", dn), ("guided", g)):
            if v is not None:
                ax.plot([v], [y + (-DODGE if m == "denovo" else DODGE)], marker="o", markersize=3.4, linestyle="none",
                        markerfacecolor=MODE_COLOR[m], markeredgecolor=MODE_COLOR[m], zorder=3)
            if ring:
                v2 = ring(k, m)
                if v2 is not None:
                    ax.plot([v2], [y + (-DODGE if m == "denovo" else DODGE)], marker="o", markersize=5.2,
                            linestyle="none", markerfacecolor="none", markeredgecolor=MODE_COLOR[m],
                            markeredgewidth=0.6, zorder=2)
    ax.set_xlim(-0.04, 1.04)
    ax.set_xticks([0, 0.5, 1])
    ax.set_xticklabels(["0", "0.5", "1"])
    ax.grid(axis="x", color=figlib.GRID, linewidth=0.5)
    ax.grid(axis="y", visible=False)
    ax.tick_params(axis="y", left=False, labelleft=False)
    ax.tick_params(axis="x", labelsize=6, pad=1.5)
    ax.spines["left"].set_visible(False)
    ax.set_title(title, fontsize=6.5, pad=3)


def _gw_class(ax, ys, row_of):
    for y, k in ys:
        for m in R.MODES:
            yy = y + (-DODGE if m == "denovo" else DODGE)
            r = row_of(k, m)
            if k[2] == "liftoff":
                if m == "guided":
                    txt = "guided: not scored (extra copies are unannotated)"
                elif r is None or r.get("status", "ok") != "ok":
                    txt = "de novo: not built yet"
                else:
                    txt = f"{int(r['pair_tp']):,} of {int(r['pair_truth']):,} pairs in one family"
                ax.text(0.02, yy, txt, fontsize=5.3, color=figlib.INK_2 if r is not None and m == "denovo" else
                        figlib.INK_3, va="center", ha="left", style="normal" if m == "denovo" else "italic")
                continue
            if r is None or r.get("status", "ok") != "ok":
                ax.text(0.02, yy, f"{R.MODE_LABEL[m]}: not built yet", fontsize=5.5, color=figlib.INK_3,
                        va="center", ha="left", style="italic")
                continue
            cls = {c: int(r[c]) for c in ("exact", "partial", "missed")}
            n = sum(cls.values())
            left = 0.0
            for c in ("exact", "partial", "missed"):
                wdt = cls[c] / n if n else 0
                if wdt:
                    ax.barh(yy, wdt, left=left, height=2 * DODGE * 0.86, color=CLASS_COLOR[c],
                            edgecolor=figlib.SURFACE, linewidth=0.4)
                left += wdt
            ax.text(1.02, yy, f"{cls['exact']}·{cls['partial']}·{cls['missed']}", ha="left", va="center",
                    fontsize=5.5, color=figlib.INK_2)
    ax.set_xlim(0, 1)
    ax.set_xticks([0, 0.5, 1])
    ax.set_xticklabels(["0", "0.5", "1"])
    ax.tick_params(axis="y", left=False, labelleft=False)
    ax.tick_params(axis="x", labelsize=6, pad=1.5)
    ax.spines["left"].set_visible(False)
    ax.grid(False)


def _gw_strip(ax, contigs, summary, ref="compara"):
    """Per-contig bipartite F against `ref` (main: Compara families, human): one row per de novo sample and one per
    guided species."""
    rows = []
    arms_seen = []
    for r in contigs:
        if r["reference"] == ref and r["arm"] not in arms_seen:
            arms_seen.append(r["arm"])
    species = {r["arm"]: r["species"] for r in contigs}
    s1 = {r["arm"]: _f(r["bip_f"]) for r in summary if r["reference"] == ref and r["substrate"] == "S1"}
    for sp in R._o1.SPECIES_ORDER:
        rows += [a for a in arms_seen if species[a] == sp and a.startswith("denovo:")]
        rows += [a for a in arms_seen if species[a] == sp and a.startswith("guided:")]
    import hashlib
    for i, a in enumerate(rows):
        mode = "denovo" if a.startswith("denovo:") else "guided"
        for r in (x for x in contigs if x["arm"] == a and x["reference"] == ref):
            st = EXPO_STYLE.get(r["exposure"], EXPO_STYLE[R._o1.NEVER])
            h = hashlib.md5(f"{a}|{r['contig']}".encode()).digest()
            jy = (h[0] / 255 - 0.5) * 0.36
            ax.scatter([_f(r["bip_f"])], [i + jy], s=st["s"], marker=st["marker"], linewidths=0.6,
                       facecolors="none" if st.get("open") else MODE_COLOR[mode], edgecolors=MODE_COLOR[mode],
                       alpha=0.85, zorder=3)
        if s1.get(a) is not None:
            ax.plot([s1[a], s1[a]], [i - 0.38, i + 0.38], color=figlib.INK, lw=1.0, zorder=4)
    if not rows:
        ax.text(0.5, 0.5, "not built yet", transform=ax.transAxes, ha="center", va="center", fontsize=6,
                color=figlib.INK_3, style="italic")
    ax.set_yticks(range(len(rows)))
    ax.set_yticklabels([f"de novo, {species[a]} {_short(a.split(':', 1)[1])}" if a.startswith("denovo:")
                        else f"guided, {species[a]} annotation" for a in rows], fontsize=5.8)
    ax.set_ylim(len(rows) - 0.5, -0.5)
    ax.set_xlim(-0.02, 1.02)
    ax.set_xlabel("Bipartite F per contig (the genome-wide run restricted to one contig); black tick = genome minus "
                  "development contigs", labelpad=2, fontsize=6)
    ax.grid(axis="y", visible=False)
    ax.grid(axis="x", color=figlib.GRID, linewidth=0.5)


def _gw_scatter(ax, fig, perfam, sid, species, ref, title):
    """Per-family F, de novo (sample) vs guided (species), families at least one mode recovers."""
    import hashlib
    dn = {r["family_id"]: r for r in perfam if r["arm"] == f"denovo:{sid}" and r["reference"] == ref
          and r["substrate"] in ("S1", "chr16_inset")}
    gd = {r["family_id"]: r for r in perfam if r["arm"] == f"guided:{species}" and r["reference"] == ref
          and r["substrate"] in ("S1", "chr16_inset")}
    ax.set_xlim(-0.06, 1.06)
    ax.set_ylim(-0.06, 1.06)
    for axis in (ax.xaxis, ax.yaxis):
        axis.set_ticks([0, 0.5, 1])
        axis.set_ticklabels(["0", "0.5", "1"])
    ax.set_aspect("equal")
    ax.grid(axis="both", color=figlib.GRID, linewidth=0.5)
    ax.set_xlabel("de novo, per-family F", labelpad=2)
    ax.set_ylabel("guided, per-family F", labelpad=2)
    keys = sorted(set(dn) | set(gd))
    have_dn = any(r["arm"] == f"denovo:{sid}" and r["reference"] == ref for r in perfam)
    have_g = any(r["arm"] == f"guided:{species}" and r["reference"] == ref for r in perfam)
    if not keys or not (have_dn and have_g):
        miss = [m for m, h in (("de novo", have_dn), ("guided", have_g)) if not h]
        ax.set_title(f"{title}\n{' and '.join(miss) or 'both modes'} not built yet", loc="left", fontsize=6, pad=3)
        return
    ax.plot([0, 1], [0, 1], color=figlib.INK_3, linewidth=0.6, linestyle=(0, (2, 2)), zorder=1)
    marks, meta = [], {}
    g_hi = d_hi = eq = 0
    for k in keys:
        x = _f(dn[k]["f"]) if k in dn else 0.0
        y = _f(gd[k]["f"]) if k in gd else 0.0
        if x == 0 and y == 0:
            continue
        g_hi += y > x
        d_hi += x > y
        eq += x == y
        h = hashlib.md5(f"{ref}|{k}".encode()).digest()
        xj, yj = x + (h[0] / 255 - 0.5) * 0.03, y + (h[1] / 255 - 0.5) * 0.03
        ax.scatter([xj], [yj], s=7, marker="o", color=figlib.INK_2, alpha=0.55, linewidths=0, zorder=3, clip_on=False)
        r = dn.get(k) or gd.get(k)
        kk = ("", k)
        meta[kk] = {"flagship": r["flagship"], "label": r["label"], "status": R.UNTOUCHED}
        marks.append((kk, xj, yj))
    ax.set_title(f"{title}\nguided higher {g_hi} · de novo higher {d_hi} · equal {eq}", loc="left", fontsize=6,
                 pad=3, linespacing=1.25)
    if any(meta[k]["flagship"] for k, _, _ in marks):
        # flagship labels without the contig suffix of the per-chromosome figure
        for k in meta:
            meta[k]["label"] = meta[k]["label"] or meta[k]["flagship"]
        _label_flagships_gw(ax, fig, marks, meta)


def _label_flagships_gw(ax, fig, marks, meta):
    """_label_flagships (collision-free placement) for genome-wide families: its label reads '<name> (<key[0]>)', so
    the key's contig slot carries a placeholder that is stripped from the placed labels."""
    tag = "genome"
    _label_flagships(ax, fig, [((tag, k[1]), x, y) for k, x, y in marks], {(tag, k[1]): v for k, v in meta.items()},
                     suffix=False)


GW_MAIN_SCATTERS = [("human_A119b", "human", "compara", "e", "Compara families (primates)"),
                    ("human_testis", "human", "compara", "f", "Compara families (primates)"),
                    ("human_A119b", "human", "soto", "g", "Soto 2025 (not independent)")]
GW_SUPP_SCATTERS = [("human_A119b", "human", "protein_homology", "e", "protein-homology families"),
                    ("gorilla_OR6737", "gorilla", "protein_homology", "f", "protein-homology families"),
                    ("chimp_PTR", "chimpanzee", "protein_homology", "g", "protein-homology families")]


def _plot_genome(data_dir: Path, out_dir: Path, refs=MAIN_REFS, name="fig7_family_recovery", strip_ref="compara",
                 scatters=GW_MAIN_SCATTERS, headline=None):
    import matplotlib.pyplot as plt
    import matplotlib.transforms as mtransforms
    from matplotlib.patches import Patch
    from matplotlib.lines import Line2D

    summary = [r for r in figlib.read_table("fig7_gw_summary", data_dir) if r["reference"] in refs]
    contigs = [r for r in figlib.read_table("fig7_gw_contigs", data_dir) if r["reference"] in refs]
    perfam = [r for r in figlib.read_table("fig7_gw_per_family", data_dir) if r["reference"] in refs]
    if not summary:
        return []
    idx = {}
    for r in summary:
        idx[(r["arm"], r["reference"], r["substrate"])] = r
    labels = {r["sample"]: r["sample_label"] for r in summary if r["mode"] == "denovo"}

    def row_of(k, mode, sub=None):
        sp, sid, ref = k
        arm = f"denovo:{sid}" if mode == "denovo" else f"guided:{sp}"
        return idx.get((arm, ref, sub or ("chr16_inset" if ref == "npip_u2" else "S1")))

    def val(k, mode, key):
        r = row_of(k, mode)
        return _f(r[key]) if r is not None and r.get("status", "ok") == "ok" else None

    def ring(k, mode):
        if k[2] == "npip_u2":
            return None
        r = row_of(k, mode, "S2")
        return _f(r["bip_f"]) if r is not None else None

    ys, y, rowlabels, groups = [], 0.0, [], []
    for kind, payload in _gw_rows(summary, refs):
        if kind == "species":
            y += 0.35 if groups else 0.0
            groups.append((y, kind, payload, y - 0.5 if groups else None))
            y += 0.8
        elif kind == "status":
            y += 0.2
            groups.append((y, kind, payload, y - 0.5))
            y += 0.8
        else:
            sp, sid, ref = payload
            r0 = row_of(payload, "denovo") or row_of(payload, "guided")
            if ref == "liftoff":
                n = (f"{int(r0['pair_truth']):,} copy pairs, both loci read-supported"
                     if r0 and r0.get("status", "ok") == "ok" else "not built yet")
            else:
                n = f"{int(r0['truth_families']):,} families · {int(r0['truth_genes']):,} genes" if r0 else ""
            rowlabels.append((y, f"{_short(sid)} · {R.GW_TRUTH_SHORT[ref]}", n))
            ys.append((y, payload))
            y += 1.0
    y_lo, y_hi = -0.75, y - 0.4
    dot_in = (y_hi - y_lo) * IN_PER_UNIT
    strip_rows = len({r["arm"] for r in contigs if r["reference"] == strip_ref})
    strip_in = max(1.0, 0.16 * strip_rows + 0.45)
    TOP = 0.78
    H = TOP + dot_in + 0.45 + strip_in + 0.55 + BOTTOM_IN - 0.15
    W = figlib.WIDTH_DOUBLE
    fy = lambda inch_from_top: 1 - inch_from_top / H           # noqa: E731

    fig = plt.figure(figsize=(W, H))
    fig.text(0.012, fy(0.10), headline or ("Rustle-internal comparison of its two family modes (de novo = the default "
             "definition, from the reads; guided = the same rule on the annotation). Not a tool comparison."),
             fontsize=6.4, fontweight="bold", ha="left", va="center", color=figlib.INK)
    top = fig.add_gridspec(1, 7, left=0.255, right=0.935, top=fy(TOP), bottom=fy(TOP + dot_in), wspace=0.13,
                           width_ratios=[1, 1, 0.35, 1, 1, 1, 1.9])
    axes = [fig.add_subplot(top[0, i]) for i in (0, 1, 3, 4, 5, 6)]
    for ax in axes[1:]:
        ax.sharey(axes[0])
    cols = [("pair_sens", "Sensitivity"), ("pair_prec", "Precision"), ("bip_sens", "Sensitivity"),
            ("bip_prec", "Precision"), ("bip_f", "F (pooled)")]
    for i, (ax, (key, title)) in enumerate(zip(axes[:5], cols)):
        _gw_dot(ax, ys, val, key, title, ring=ring if key == "bip_f" else None)
    _gw_class(axes[5], ys, row_of)
    axes[0].set_ylim(y_hi, y_lo)
    tr = mtransforms.blended_transform_factory(fig.transFigure, axes[0].transData)
    for yy, text, n in rowlabels:
        fig.text(0.247, yy - 0.2, text, transform=tr, ha="right", va="center", fontsize=5.8, color=figlib.INK)
        fig.text(0.247, yy + 0.25, n, transform=tr, ha="right", va="center", fontsize=5.3, color=figlib.INK_2)
    for yy, kind, text, sep in groups:
        fig.text(0.012 if kind == "species" else 0.024, yy, text, transform=tr, ha="left", va="center",
                 fontsize=6.5 if kind == "species" else 6.0, fontweight="bold" if kind == "species" else "normal",
                 color=figlib.INK if kind == "species" else figlib.INK_2)
        if sep is not None:
            for ax in axes:
                ax.axhline(sep, color=figlib.INK_3 if kind == "species" else figlib.GRID, linewidth=0.5, zorder=0)

    def header(ax0, ax1, letter, text):
        x0, x1 = ax0.get_position().x0, ax1.get_position().x1
        fig.text(x0, fy(TOP - 0.34), letter, fontsize=9, fontweight="bold", ha="left", va="bottom")
        fig.text(x0 + 0.018, fy(TOP - 0.335), text, fontsize=6.3, ha="left", va="bottom", color=figlib.INK)
        fig.add_artist(plt.Line2D([x0, x1], [fy(TOP - 0.30), fy(TOP - 0.30)], transform=fig.transFigure,
                                  color=figlib.INK_3, linewidth=0.5))
    header(axes[0], axes[1], "a", "Pairwise (gene pairs)")
    header(axes[2], axes[4], "b", "One-to-one bipartite (family members)")
    header(axes[5], axes[5], "c", "Reference families")
    for i, m in enumerate(R.MODES):
        yk = fy(0.28 + 0.14 * i)
        fig.add_artist(plt.Line2D([0.016], [yk], transform=fig.transFigure, marker="o", markersize=3.4,
                                  markerfacecolor=MODE_COLOR[m], markeredgecolor=MODE_COLOR[m], linestyle="none"))
        fig.text(0.024, yk, f"{R.MODE_LABEL[m]} ({'upper: per sample' if m == 'denovo' else 'lower: per species'})",
                 fontsize=5.8, va="center", ha="left", color=figlib.INK)
    fig.text(0.012, fy(0.62), "dots: genome minus development contigs (human chr16,\ngorilla chr20); ring in F: minus "
             "every contig used\nfor a family decision; same family rule in both modes", fontsize=5.3, va="center",
             ha="left", color=figlib.INK_3, linespacing=1.15)
    axes[5].legend(handles=[Patch(facecolor=CLASS_COLOR[c], edgecolor=figlib.INK_3 if c == "missed" else
                                  CLASS_COLOR[c], linewidth=0.4, label=c) for c in ("exact", "partial", "missed")],
                   loc="lower center", bbox_to_anchor=(0.5, 1.0), ncol=3, fontsize=5.8, handlelength=0.9,
                   handleheight=0.8, columnspacing=0.8, handletextpad=0.35, borderaxespad=0.2)

    # ---- d: per-contig strip
    s_top = TOP + dot_in + 0.45
    ax_s = fig.add_axes([0.255, fy(s_top + strip_in), 0.68, strip_in / H])
    _gw_strip(ax_s, contigs, summary, strip_ref)
    fig.text(0.012, fy(s_top - 0.12), "d", fontsize=9, fontweight="bold", ha="left", va="bottom")
    fig.text(0.030, fy(s_top - 0.115), f"{R.GW_TRUTH_SHORT.get(strip_ref, strip_ref).capitalize()} families, per "
             "contig: no contig carries the result", fontsize=6.3, ha="left", va="bottom")
    handles = []
    for cls, text in EXPO_LEGEND:
        st = EXPO_STYLE[cls]
        handles.append(Line2D([], [], linestyle="none", marker=st["marker"], markersize=st["s"] ** 0.5 * 1.1,
                              markerfacecolor="none" if st.get("open") else figlib.INK_2,
                              markeredgecolor=figlib.INK_2, label=text))
    ax_s.legend(handles=handles, loc="lower right", bbox_to_anchor=(1.0, 1.0), ncol=4, fontsize=5.5, frameon=False,
                handletextpad=0.3, columnspacing=0.9, borderaxespad=0.1)

    # ---- e-g: per-family scatters
    b_top = s_top + strip_in + 0.55
    sc_top, sc_bot = b_top + 0.2 + 0.52, b_top + 0.2 + 0.52 + 1.45
    bot = fig.add_gridspec(1, 3, left=0.075, right=0.975, top=fy(sc_top), bottom=fy(sc_bot), wspace=0.55)
    sc = [fig.add_subplot(bot[0, i]) for i in range(len(scatters))]
    for ax, (sid, sp, ref, letter, t) in zip(sc, scatters):
        _gw_scatter(ax, fig, perfam, sid, sp, ref, f"{sp.capitalize()} {_short(sid)} · {t}")
        fig.text(ax.get_position().x0 - 0.075, fy(sc_top - 0.34), letter, fontsize=9, fontweight="bold", ha="left",
                 va="bottom")
    fig.text(0.075, fy(b_top + 0.12), "Per-family F, genome minus development contigs: one point per reference family "
             "that at least one mode recovers (jitter ≤ 0.015)", fontsize=6, va="center", ha="left", color=figlib.INK)
    figlib.stamp_provisional(fig, GW_TABLES, data_dir)
    paths = figlib.save(fig, name, out_dir)
    plt.close(fig)
    return paths


SUPP_HEADLINE = ("Supplementary · SECONDARY reference: the annotation's protein-homology families (fold level). Protein "
                 "is not part of Rustle's family definition.")


def plot(data_dir: Path, out_dir: Path):
    figlib.use_style()
    data_dir, out_dir = Path(data_dir), Path(out_dir)
    if _gw_ready(data_dir):
        return (_plot_genome(data_dir, out_dir)
                + _plot_genome(data_dir, out_dir, refs=SUPP_REFS, name="fig7s_protein_homology",
                               strip_ref="protein_homology", scatters=GW_SUPP_SCATTERS, headline=SUPP_HEADLINE))
    return (_plot_dev(data_dir, out_dir)
            + _plot_dev(data_dir, out_dir, truths=("referee",), name="fig7s_protein_homology",
                        scatters=DEV_SUPP_SCATTERS, banner=SUPP_HEADLINE))


# ================================================================ caption numbers
def summary(data_dir: Path = figlib.DATA_DIR):
    """Print every number captions/fig7.md quotes (genome-wide tables when present, else the development tables),
    with the pre-registered claims."""
    data_dir = Path(data_dir)
    if not _gw_ready(data_dir):
        if all((data_dir / f"{t}.tsv").exists() for t in DEV_TABLES):
            _summary_dev(data_dir)
        else:
            print("genome-wide tables absent: build with `make.py data fig7` (repeat while it exits 75)")
        return
    rows = figlib.read_table("fig7_gw_summary", data_dir)
    for n in figlib.table_meta("fig7_gw_summary", data_dir).get("note", []):
        if n.lower().startswith("provisional"):
            print("PROVISIONAL:", n)
    idx = {(r["arm"], r["reference"], r["substrate"]): r for r in rows}
    print(f"{'sample':15s} {'reference':17s} {'sub':11s} {'fams·genes':>13s}  {'dn sens/prec/F':>20s}  "
          f"{'g sens/prec/F':>20s}  gap(g-dn)  dn e·p·m / g e·p·m")
    for r in rows:
        if r["mode"] != "denovo" or r["reference"] == "liftoff":
            continue
        g = idx.get((f"guided:{r['species']}", r["reference"], r["substrate"]))

        def t(x):
            return f"{_f(x['bip_sens']):.3f}/{_f(x['bip_prec']) or 0:.3f}/{_f(x['bip_f']):.3f}" if x else "-"
        gap = f"{_f(g['bip_f']) - _f(r['bip_f']):+.3f}" if g else "-"
        print(f"{r['sample']:15s} {r['reference']:17s} {r['substrate']:11s} {r['truth_families']:>5s}·"
              f"{r['truth_genes']:>6s}  {t(r):>20s}  {t(g):>20s}  {gap:>8s}  {r['exact']}·{r['partial']}·{r['missed']}"
              f" / {g['exact'] + '·' + g['partial'] + '·' + g['missed'] if g else '-'}")
    print("== Liftoff copy pairs (de novo only; F7.3, descriptive)")
    for r in rows:
        if r["reference"] == "liftoff":
            if r["status"] != "ok":
                print(f"  {r['sample']:15s} {r['status']}")
            else:
                print(f"  {r['sample']:15s} {r['substrate']} {r['pair_tp']}/{r['pair_truth']} = {r['pair_sens']} "
                      f"(both loci covered: {r['pair_covered']})")
    for t in ("fig7_gw_compara_levels", "fig7_gw_bridge"):
        if (Path(data_dir) / f"{t}.tsv").exists():
            print(f"== {t}")
            for r in figlib.read_table(t, data_dir):
                print("  " + " ".join(f"{k}={v}" for k, v in r.items() if k not in ("sample_label", "minimap2_flags")))


def _summary_dev(data_dir: Path):
    rows = figlib.read_table("fig7_summary", data_dir)
    idx = {(r["species"], r["contig"], r["truth"], r["mode"]): r for r in rows}
    print("== development tables (per-chromosome runs); main: compara / soto / npip_u2; supplement: referee")
    for (sp, c, t, m), r in sorted(idx.items(), key=lambda kv: (kv[0][0], kv[0][1], TRUTH_ORDER.index(kv[0][2])
                                                               if kv[0][2] in TRUTH_ORDER else 9, kv[0][3])):
        if m != "denovo":
            continue
        g = idx.get((sp, c, t, "guided"))
        gap = f"{float(g['bip_f']) - float(r['bip_f']):+.3f}" if g else "-"
        print(f"  {sp:7s} {c:12s} {r['status'][:28]:28s} {t:8s} {r['truth_families']:>3s} fam {r['truth_genes']:>4s} "
              f"genes | dn sens/prec/F {float(r['bip_recall']):.3f}/{float(r['bip_precision']):.3f}/"
              f"{float(r['bip_f']):.3f} | g "
              + (f"{float(g['bip_recall']):.3f}/{float(g['bip_precision']):.3f}/{float(g['bip_f']):.3f}" if g else "-")
              + f" | gap {gap} | pairs dn {r['pair_tp']}/{r['truth_pairs']} g "
              + (f"{g['pair_tp']}/{g['truth_pairs']}" if g else "-")
              + f" | exact dn {r['exact']} g {g['exact'] if g else '-'}")


if __name__ == "__main__":
    if sys.argv[1:2] == ["summary"]:
        summary(Path(sys.argv[2]) if len(sys.argv) > 2 else figlib.DATA_DIR)
    else:
        sys.exit("usage: python3 figures/fig_family_recovery.py summary [DATA_DIR]")
