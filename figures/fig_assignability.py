"""Figure 4 — what happens to the simulated reads the aligner cannot place (MAPQ 0), one bar per sample; the full
UpSet of every sample is the supplementary Figure 4s.

Reads are simulated from every catalog copy with a spliced sequence of >= 300 bp in a family of >= 2 copies (the copy's
own spliced sequence; the read name records the source copy), mapped to the whole genome with the shipped minimap2
settings, and assigned by `copy_assign --families` (default output: one result per read and family) and
`--union-certificate` (the union test). Genome-wide on every sample (docs/archive/2026-09/PREREG_genome_wide_copy_assignment_2026-09-25.md;
cfg `o2_scope genome`), or the development tables (`o2_scope dev`: human A119b chr16 catalog, gorilla OR6737 chr20
(NC_073244.2) catalog). Samples and species are never pooled.

Main figure (fig4_assignability): per sample, the MAPQ-0 reads split into five fates (_o2.FATES, fixed in the
pre-registration), scored within the read's source family, which only the simulation knows; beside the bar, the
counts a user can apply (default output, union test) and the NM-identical twins.

Supplement (fig4s_assignability_upset): one UpSet per sample over ALL simulated reads; each read sits in exactly one
column (its exact combination of sets); bars are stacked by the identity band of the read's source copy to its most
similar directly aligned copy (figlib.IDENTITY_BANDS; the top band, "100%*", is identity over the aligned segment,
which covers >= 50% of the shorter copy, not full-length identity). Sets (_o2.SET_LABEL; verdicts from
`bench/score.py reads --per-read`):
  unique      the aligner's primary has MAPQ > 0 (copy_assign's gate leaves these to the aligner)
  mapq0       the primary has MAPQ 0: the aligner could not choose a placement
  not_scored  MAPQ 0 but no result for the read in the default output
  identical   0 decisive sites over every placement the read has (the union run's result for its source family)
  family_psv  >= 1 decisive site (a PSV or splice junction where the candidate copies differ) in its source family
  assigned    assigned in its source family's result (not origin-rejected); the simulation selects that result
  correct     ... to the source copy (or a catalog copy at the same locus: score.py reads' locus rule)
"""
from __future__ import annotations

import collections
from pathlib import Path

import figlib
import _o2

META = {
    "id": "fig4",
    "title": "Fate of the simulated reads the aligner cannot place, per sample",
    "claim": ("Development tables (human A119b chr16 catalog; gorilla OR6737 chr20 catalog): scored within the "
              "read's source family (known only in simulation), the copy-assignment test made no wrong call where it "
              "assigned: 163 of the 1,263 human MAPQ-0 reads, all to the source copy (21 source copies; 153 through a "
              "decisive site, 10 with a single candidate copy). 858 of the 1,263 come from copies identical to a "
              "sibling over the aligned segment, and 1 of them is assigned. The outputs a user can apply do not "
              "reproduce this: the default output is correct for 60 of the 648 MAPQ-0 reads it assigns, and the "
              "union test assigns none, because 1,255 of the 1,263 have an NM-identical twin at another locus. "
              "Gorilla chr20: 30 MAPQ-0 reads from three copies, none assigned. The genome-wide tables of all six "
              "samples replace these numbers (docs/archive/2026-09/PREREG_genome_wide_copy_assignment_2026-09-25.md)."),
    "tables": ["fig4_assignability_upset"],
}
TABLE = "fig4_assignability_upset"
PROVISIONAL = ("provisional: from the recorded 2026-09-23/24 run (o2sim h16*/g44*, hash-seeded simulation); "
               "make.py data regenerates with the stable seed and current defaults")
BAND_LABEL = {"identical": "100%* over the aligned segment", "99.5-100%": "99.5–100%", "99-99.5%": "99–99.5%",
              "98-99%": "98–99%", "<98%": "< 98%", "unknown": "unknown"}
FATE_COLOR = {"correct": figlib.BLUE[650], "wrong": "#b5452f", "psv_unassigned": figlib.BLUE[250],
              "no_site": "#a9a8a2", "no_result": "#dddcd6"}
ROW_LABEL = {"human_A119b": "A119b (CHM13 v2.0)", "human_testis": "testis (CHM13 v2.0)",
             "chimp_PTR": "PTR (mPanTro3)", "gorilla_OR6737": "OR6737, testis (mGorGor1)",
             "gorilla_KB3781": "KB3781, fibroblast (mGorGor1)", "orangutan_PPY": "PPY (mPonPyg2)"}
FOOTNOTE = ("* Scored within the read's source family, which is known only in simulation; a user of the output cannot "
            "choose that result. Decisive site: a PSV or splice junction, covered by the read, at which the candidate "
            "copies differ. NM-identical twin: an alignment at another locus with the best alignment score and the "
            "same number of mismatches and indels as the source-copy alignment.")

# saved with bbox_inches="tight"; the provisional stamp reaches x = 0.995, so the canvas stays under 183 mm
FIG_W = figlib.WIDTH_DOUBLE - 6 * figlib.MM


# ================================================================ data
def build(cfg: dict, data_dir: Path, force: bool = False, recorded: bool = False):
    """Regenerate `fig4_assignability_upset` (HEAVY unless cached: _o2.collect; genome scope = every sample's
    genome-wide simulation and sharded copy_assign runs, bounded per call — re-run until it finishes).

    recorded=True tabulates the recorded runs in cfg['o2sim_dir'] instead (development / provisional tables)."""
    rows, notes, inputs = [], [], {}
    for runs, logs, reads in _o2.collect(cfg, force=force, recorded=recorded):
        sid, scope = runs["sample"], runs["catalog_scope"]
        notes += _o2.run_notes(runs, reads, logs)
        notes += _summary_notes(runs, reads)
        rows += _o2.upset_rows(reads, scope)
        inputs.update({f"{sid}_sim_bam": f"{runs['sim']}.bam", f"{sid}_catalog": runs["catalog"],
                       f"{sid}_assignments_default": f"{runs['o2']}.assignments.tsv",
                       f"{sid}_assignments_union": f"{runs['u2']}.assignments.tsv",
                       f"{sid}_score_reads_o2": logs["o2"], f"{sid}_score_reads_u2": logs["u2"],
                       f"{sid}_per_read_o2": logs["o2_per_read"], f"{sid}_per_read_u2": logs["u2_per_read"]})
    notes += [
        "one row per (sample, exact set combination, outcomes, twin state, aligner outcome, identity band); n_reads "
        "counts simulated reads; catalog_scope = genome (the sample's genome-wide copy catalog) or the one contig of "
        "a development catalog",
        "sets: unique = primary MAPQ > 0; mapq0 = primary MAPQ 0; not_scored = MAPQ 0 and no result for the read in "
        "the default output; identical = the union run's result for the read's source family has 0 decisive sites; "
        "family_psv = the default run's result for the source family has >= 1 decisive site (a PSV or splice junction "
        "covered by the read at which the candidate copies differ); assigned = that result is 'assigned' and not "
        "origin-rejected (scored within the source family, which only the simulation knows); correct = assigned to "
        "the source copy or a catalog copy at the same locus (overlap >= 50% of the shorter span, score.py reads' "
        "rule); own_row = the default output has a result for the read in its source family; the test sets are 0 "
        "for MAPQ > 0 reads (the gate does not score them); all from score.py reads --per-read",
        "foreign_claim = an assigned result of the default output names a locus other than the source's (score.py "
        "reads --per-read wrong_locus_rows > 0; the union test removes these)",
        "any_outcome / union_outcome (MAPQ-0 reads; 'na' otherwise) = the verdict a user can apply: any assigned "
        "result of the default output / of the union run (score.py reads ANY: correct, wrong, conflict = two loci "
        "claimed, abstain, not_scored = no result)",
        "twin (MAPQ-0 reads; 'na' otherwise) = _o2.twin_state over every alignment of the read (-F 2052): nm_twin = an "
        "alignment not overlapping the source copy has the best alignment score and the NM of the source-copy "
        "alignment; as_tie_nm_differs = alignments elsewhere reach the best score, none with that NM; "
        "no_equal_as_elsewhere = the source copy holds the best score alone; true_below_best = the best score is "
        "elsewhere; no_true_placement = no alignment overlaps the source copy",
        "aligner_primary = the catalog copy with the largest raw overlap of the primary alignment's reference span: "
        "true_copy (the source copy's locus) / other_copy / outside_catalog / unmapped",
        "identity_band = figlib.identity_band(max_family_identity of the source copy): its identity to the most "
        "similar directly aligned copy of its family; the top band 'identical' = 100% over the aligned segment, "
        "which covers >= 50% of the shorter copy, not full-length identity",
    ]
    if recorded:
        notes.insert(0, PROVISIONAL)
    return figlib.write_table(TABLE, _o2.UPSET_HEADER, rows, generator="figures/fig_assignability.py build()",
                              inputs=inputs, notes=notes, data_dir=data_dir)


def _summary_notes(runs: dict, reads: list[dict]) -> list[str]:
    """Counts the caption quotes that are not a sum of table rows: copies simulated / passed to copy_assign, the
    source copies behind the assigned reads, and the twin_state tally."""
    def rows(path):
        with open(path) as fh:
            return [l.rstrip("\n").split("\t") for l in list(fh)[1:] if l.strip()]
    sid = runs["sample"]
    used, cat = rows(f"{runs['sim']}.copies_used.tsv"), rows(runs["catalog"])
    s = (f"{sid}: simulated {len(used)} copies ({sum(int(r[6]) for r in used)} reads; spliced sequence >= 300 bp in "
         f"a family of >= 2 copies); catalog passed to copy_assign {len(cat)} copies in {len({r[0] for r in cat})} "
         f"families")
    if runs.get("source_catalog"):
        src = rows(runs["source_catalog"])
        s += (f"; source catalog {runs['source_catalog']}: {len(src)} copies in {len({r[0] for r in src})} families, "
              f"{len(src) - len(cat)} dropped (_o2.derive_catalog: not simulated and no simulated primary over the "
              f"span)")
    z = [r for r in reads if r["mapq"] == 0]
    asg = [r for r in z if r["assigned"]]
    sole = [r for r in asg if not r["family_psv"]]
    sole_copies = sorted({f"{r['family']}|{r['copy']}" for r in sole})
    tw = collections.Counter(r["twin"] for r in z)
    return [s,
            f"{sid}: assigned within the source family: {len(asg)} MAPQ-0 reads from "
            f"{len({(r['family'], r['copy']) for r in asg})} source copies; of them without a decisive site (a single "
            f"candidate copy): {len(sole)} reads from {len(sole_copies)} source "
            f"{'copy' if len(sole_copies) == 1 else 'copies'} ({', '.join(sole_copies[:20]) or 'none'}"
            f"{' ...' if len(sole_copies) > 20 else ''})",
            f"{sid}: twin state of the {len(z)} MAPQ-0 reads: " + "; ".join(f"{k} {tw[k]}" for k in _o2.TWIN_STATES)]


# ================================================================ per-sample tallies (both figures)
def sample_summary(rows: list[dict]) -> dict:
    """Headline numbers of one sample's table rows."""
    n = lambda pred: sum(int(r["n_reads"]) for r in rows if pred(r))  # noqa: E731
    z = lambda r: r["mapq0"] == "1"  # noqa: E731
    made = ("correct", "wrong", "conflict")
    fate = collections.Counter()
    for r in rows:
        f = _o2.fate_of(r)
        if f:
            fate[f] += int(r["n_reads"])
    return {
        "total": n(lambda r: True),
        "unique": n(lambda r: r["unique"] == "1"),
        "mapq0": n(z),
        "unmapped": n(lambda r: r["unique"] == "0" and r["mapq0"] == "0"),
        "fate": fate,
        "assigned": n(lambda r: z(r) and r["assigned"] == "1"),
        "wrong": n(lambda r: z(r) and r["assigned"] == "1" and r["correct"] == "0"),
        "any_correct": n(lambda r: z(r) and r["any_outcome"] == "correct"),
        "any_assigned": n(lambda r: z(r) and r["any_outcome"] in made),
        "union_assigned": n(lambda r: z(r) and r["union_outcome"] in made),
        "twins": n(lambda r: z(r) and r["twin"] == "nm_twin"),
        "aligner_unique_ok": n(lambda r: r["unique"] == "1" and r["aligner_primary"] == "true_copy"),
        "aligner_mapq0_ok": n(lambda r: z(r) and r["aligner_primary"] == "true_copy"),
        "legacy": "own_row" not in rows[0] if rows else False,
    }


def by_sample(rows: list[dict]) -> dict:
    out = collections.defaultdict(list)
    for r in rows:
        out[_o2.table_sample(r)].append(r)
    return out


# ================================================================ main figure: fates
def _fate_panel(fig, per: dict):
    import matplotlib.patches as mpatches

    samples = _o2.SAMPLE_ORDER
    # rows top to bottom, a gap between species
    ys, y, prev = {}, 0.0, None
    for s in samples:
        sp = _o2.SAMPLE_SPECIES[s]
        if prev is not None and sp != prev:
            y += 0.55
        ys[s] = y
        y += 1.0
        prev = sp
    y_max = y - 1.0
    x_bar, w_bar = 0.19, 0.30
    top, bottom = 0.80, 0.215
    ax = fig.add_axes([x_bar, bottom, w_bar, top - bottom])
    ax.set_xlim(0, 1)
    ax.set_ylim(y_max + 0.6, -0.6)
    ax.grid(axis="y", visible=False)
    ax.grid(axis="x", color=figlib.GRID, linewidth=0.5)
    ax.set_xticks([0, 0.25, 0.5, 0.75, 1.0])
    ax.set_xticklabels(["0", "25", "50", "75", "100"])
    ax.set_xlabel("MAPQ-0 reads (%)", labelpad=1.5)
    ax.set_yticks([])
    ax.spines["left"].set_visible(False)
    import matplotlib.transforms as mtrans
    tr = mtrans.blended_transform_factory(fig.transFigure, ax.transData)
    # table columns (figure x), header above the bars
    # (header, figure x of the column centre)
    cols = [("Simulated\nreads", 0.545), ("MAPQ 0\n(share)", 0.615), ("Wrong* /\nassigned*", 0.69),
            ("Default output:\ncorrect / assigned", 0.79), ("Union test:\nassigned", 0.883),
            ("NM-identical\ntwin", 0.955)]
    y_head = -0.95
    for label, cx in cols:
        ax.text(cx, y_head, label, transform=tr, ha="center", va="bottom", fontsize=5.9,
                color=figlib.INK, linespacing=1.05)
    ax.text(x_bar - 0.185, y_head, "Sample (copy catalog)", transform=tr, ha="left", va="bottom", fontsize=5.9)
    ax.text(0.5, y_head, "Fate of the MAPQ-0 reads", ha="center", va="bottom", fontsize=5.9)
    prev = None
    for s in samples:
        yy = ys[s]
        sp = _o2.SAMPLE_SPECIES[s]
        if sp != prev:
            ax.text(x_bar - 0.185, yy - 0.62, _o2.SPECIES_TITLE[sp], transform=tr, ha="left", va="center",
                    fontsize=6.2, fontweight="bold")
        prev = sp
        rows = per.get(s)
        scope = _o2.table_scope(rows[0]) if rows else ""
        label = ROW_LABEL[s]
        cat = ("genome-wide catalog" if scope == "genome" else f"{scope} catalog only" if scope else "")
        ax.text(x_bar - 0.185, yy, f"{label}" + (f"\n{cat}" if cat else ""), transform=tr, ha="left",
                va="center", fontsize=5.9, linespacing=1.0, color=figlib.INK)
        if not rows:
            ax.text(0.01, yy, "not run yet (genome-wide catalog and simulation pending)", ha="left", va="center",
                    fontsize=5.6, color=figlib.INK_3, style="italic")
            for _, cx in cols:
                ax.text(cx, yy, "–", transform=tr, ha="center", va="center", fontsize=5.9,
                        color=figlib.INK_3)
            continue
        m = sample_summary(rows)
        n0 = m["mapq0"]
        left = 0.0
        for f in _o2.FATES:
            v = m["fate"].get(f, 0)
            if not n0 or not v:
                continue
            w = v / n0
            ax.barh(yy, w, left=left, height=0.72, color=FATE_COLOR[f], edgecolor=figlib.SURFACE, linewidth=0.5)
            if w >= 0.075:
                dark = f in ("correct", "wrong")
                ax.text(left + w / 2, yy, f"{v:,}", ha="center", va="center", fontsize=5.4,
                        color="white" if dark else figlib.INK)
            left += w
        vals = [f"{m['total']:,}", f"{n0:,}\n({100 * n0 / m['total']:.1f}%)" if m["total"] else "0",
                f"{m['wrong']:,} / {m['assigned']:,}", f"{m['any_correct']:,} / {m['any_assigned']:,}",
                f"{m['union_assigned']:,}", f"{m['twins']:,}"]
        for (label_, cx), v in zip(cols, vals):
            bold = label_.startswith("Wrong")
            ax.text(cx, yy, v, transform=tr, ha="center", va="center", fontsize=5.9,
                    fontweight="bold" if bold else "normal", color=figlib.INK if bold else figlib.INK_2,
                    linespacing=1.0)
    handles = [mpatches.Patch(color=FATE_COLOR[f], label=_o2.FATE_LABEL[f] + ("*" if f in ("correct", "wrong")
                                                                             else "")) for f in _o2.FATES]
    fig.legend(handles=handles, loc="upper left", bbox_to_anchor=(x_bar - 0.19, 0.975), ncol=3, fontsize=6.0,
               handlelength=1.0, handleheight=0.8, handletextpad=0.4, columnspacing=1.2, frameon=False,
               title="MAPQ-0 reads, scored within the source family*", title_fontsize=6.2, alignment="left")
    return ax


def plot(data_dir: Path, out_dir: Path):
    import matplotlib.pyplot as plt
    import textwrap

    rows = figlib.read_table(TABLE, data_dir)
    per = by_sample(rows)
    fig = plt.figure(figsize=(FIG_W, 3.55))
    _fate_panel(fig, per)
    note = [FOOTNOTE]
    if any(_o2.table_scope(v[0]) != "genome" for v in per.values()):
        note.append("Development tables (one contig's catalog per sample); the genome-wide catalogs of all six "
                    "samples replace them.")
    fig.text(0.008, 0.02, "\n".join(textwrap.wrap(" ".join(note), 175)), fontsize=5.4, color=figlib.INK_2,
             va="bottom", ha="left", linespacing=1.1)
    figlib.stamp_provisional(fig, META["tables"], data_dir)
    paths = figlib.save(fig, "fig4_assignability", out_dir)
    plt.close(fig)
    return paths + plot_supplement(rows, out_dir, data_dir)


# ================================================================ supplement: one UpSet per sample
def _upset_panel(fig, spec, rows: list[dict], min_cols: int):
    import matplotlib.ticker as mticker

    sets = _o2.SETS
    counts: dict = collections.Counter()
    stack: dict = collections.defaultdict(collections.Counter)
    for r in rows:
        k = frozenset(s for s in sets if r[s] == "1")
        n = int(r["n_reads"])
        counts[k] += n
        stack[k][r["identity_band"]] += n
    bands = list(figlib.IDENTITY_BANDS) + (["unknown"] if any(r["identity_band"] == "unknown" for r in rows) else [])
    colors = dict(figlib.IDENTITY_COLOR, unknown="#bdbcb6")
    sizes = sorted((v for v in counts.values() if v > 0), reverse=True)
    uniq = frozenset({"unique"})
    break_at = None
    if len(sizes) > 1 and counts[uniq] > 3 * sizes[1]:
        low_top = 1.22 * sizes[1]
        break_at = (low_top, counts[uniq] - low_top / (2.6 * 1.35))
    ax = figlib.upset(fig, dict(counts), sets, stack=stack, stack_order=bands, stack_colors=colors,
                      set_labels=_o2.SET_LABEL, gridspec=spec, bar_ylabel="Simulated reads", set_xlabel="Reads in set",
                      break_at=break_at, max_columns=12)
    cols = ax["columns"]
    thousands = mticker.FuncFormatter(lambda v, _: f"{v:,.0f}")
    for k in ("bar", "bar_top"):
        if ax.get(k) is not None:
            ax[k].set_xlim(-0.5, max(len(cols), min_cols) - 0.5)
            ax[k].yaxis.set_major_formatter(thousands)
            ax[k].tick_params(axis="y", labelsize=5.2)
            ax[k].yaxis.label.set_size(5.8)
            for t in ax[k].texts:
                t.set_fontsize(4.6)
    if ax.get("bar_top") is not None:
        ax["bar_top"].yaxis.set_major_locator(mticker.MaxNLocator(2, integer=True))
    ax["matrix"].tick_params(axis="y", labelsize=5.2)
    for c in ax["matrix"].collections:
        c.set_sizes([6])
    for t in ax["matrix"].texts:
        t.set_fontsize(4.6)
    set_size = {s: sum(v for k, v in counts.items() if s in k) for s in sets}
    ax_set = ax["sets"]
    n_sets = len(sets)
    small = [v for s, v in set_size.items() if s != "unique"]
    scap = 1.9 * max(small) if small and max(small) > 0 else None
    for s in sets:
        y = n_sets - 1 - sets.index(s)
        v = set_size[s]
        shown = min(v, 0.9 * scap) if scap else v
        if scap and v > scap:
            for p in ax_set.patches:
                if abs(p.get_y() + p.get_height() / 2 - y) < 1e-6:
                    p.set_width(0.9 * scap)
        ax_set.text(shown, y, f"{v:,} ", ha="right", va="center", fontsize=4.6, color=figlib.INK_2)
    if scap:
        ax_set.set_xlim(scap * 1.9, 0)
    ax_set.set_xticks([])
    ax_set.spines["bottom"].set_visible(False)
    ax_set.set_xlabel("Reads in set", labelpad=2, fontsize=5.2)
    return ax, bands, colors


def _relayout(ax: dict, x0: float, x1: float, set_w: float = 0.055, gap: float = 0.175):
    """Set-size bars at the panel's left edge, a gap for the set names, then the matrix and the bars."""
    mat, sets = ax["matrix"], ax["sets"]
    pm, ps = mat.get_position(), sets.get_position()
    x_mat = x0 + set_w + gap
    sets.set_position([x0, ps.y0, set_w, ps.height])
    mat.set_position([x_mat, pm.y0, x1 - x_mat, pm.height])
    for k in ("bar", "bar_top"):
        if ax.get(k) is not None:
            p = ax[k].get_position()
            ax[k].set_position([x_mat, p.y0, x1 - x_mat, p.height])


def plot_supplement(rows: list[dict], out_dir: Path, data_dir: Path):
    import matplotlib.gridspec as gs
    import matplotlib.patches as mpatches
    import matplotlib.pyplot as plt

    per = by_sample(rows)
    fig = plt.figure(figsize=(FIG_W, 8.9))
    grid = _o2.GRID_3x2
    g = gs.GridSpec(len(grid), 2, figure=fig, hspace=0.42, wspace=0.08, left=0.01, right=0.99, top=0.9,
                    bottom=0.03)
    n_cols = {s: len({tuple(r[k] for k in _o2.SETS) for r in v}) for s, v in per.items()}
    min_cols = min(13, max(n_cols.values(), default=1))
    legend_done = False
    for i, row in enumerate(grid):
        for j, s in enumerate(row):
            letter = "abcdef"[2 * i + j]
            sub = per.get(s)
            spec = g[i, j]
            bb = spec.get_position(fig)
            if not sub:
                fig.text(bb.x0, bb.y1 + 0.022, f"{letter}  {_o2.catalog_title(s, 'genome', short=True)}",
                         fontsize=6.2, va="bottom", ha="left")
                fig.text(bb.x0 + 0.5 * bb.width, bb.y0 + 0.5 * bb.height, "not run yet", ha="center", va="center",
                         fontsize=6.5, color=figlib.INK_3, style="italic")
                continue
            ax, bands, colors = _upset_panel(fig, spec, sub, min_cols)
            _relayout(ax, bb.x0 + 0.01, bb.x1 - 0.005)
            total = sum(int(r["n_reads"]) for r in sub)
            title = _o2.catalog_title(s, _o2.table_scope(sub[0]), short=True).replace("), ", "),\n", 1)
            fig.text(bb.x0, bb.y1 + 0.022, f"{letter}  {title}: {total:,} simulated reads".replace("\n", "\n    "),
                     fontsize=6.2, va="bottom", ha="left", linespacing=1.1)
            if not legend_done:
                handles = [mpatches.Patch(color=colors[b], label=BAND_LABEL[b]) for b in reversed(bands)]
                fig.legend(handles=handles, loc="upper left", bbox_to_anchor=(0.01, 0.995), ncol=6, fontsize=5.8,
                           handlelength=1.0, handleheight=0.8, handletextpad=0.4, columnspacing=1.0, frameon=False,
                           title="Bar segments: identity of the source copy to its most similar directly aligned "
                                 "copy (* over the aligned segment, which covers ≥ 50% of the shorter copy)",
                           title_fontsize=6.0, alignment="left")
                legend_done = True
    fig.text(0.01, 0.004, "* Assigned within the read's source family, known only in simulation. Columns beyond the "
             "12 largest are pooled as 'other'.", fontsize=5.4, color=figlib.INK_2, va="bottom", ha="left")
    figlib.stamp_provisional(fig, META["tables"], data_dir)
    paths = figlib.save(fig, "fig4s_assignability_upset", out_dir)
    plt.close(fig)
    return paths
