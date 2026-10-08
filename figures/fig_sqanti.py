"""fig2 — are the assembled transcripts real and complete? SQANTI3 structural categories and rules filter per method.

Tables (figures/data/):

  fig2_sqanti_categories  species, scope, tool, structural_category, category, n, n_total, frac
                          SQANTI3 QC structural_category counts per method; `category` folds them into
                          figlib.SQANTI_ORDER (FSM / ISM / NIC / NNC / Other = antisense, intergenic, genic,
                          genic_intron, fusion, moreJunctions); frac = n / n_total (all isoforms of the method).
  fig2_sqanti_filter      species, scope, tool, n_total, n_pass, pass_frac, n_fsm, fsm_frac, n_fsm_pass,
                          n_multiexon, n_fsm_multiexon, fsm_multiexon_frac
                          SQANTI3 rules filter (default JSON): PASS = filter_result "Isoform"; the multi-exon
                          columns count isoforms with > 1 exon (FSM share without mono-exon emission).
  fig2_sqanti_subcategories  species, scope, tool, structural_category, subcategory, n, n_multiexon, n_total,
                          n_total_multiexon, frac
                          SQANTI3 subcategories (completeness: FSM reference_match = both ends within 50 bp;
                          ISM 3prime_fragment = missing the reference's 5' exons); cited by the caption, not drawn.
  fig2_samples_categories, fig2_samples_filter   (supplementary; drawn as fig2s_samples) the same two tables with
                          `sample` (registry id), `label` and `tissue` first, for Rustle's two configurations on EVERY
                          sample of the registry, each against its own annotation (chimpanzee and orangutan: the RefSeq
                          GFF3 converted by gff_to_gtf, assembly.annotation_gtf).

  fig2_guided_categories, fig2_guided_filter   (supplementary; ONLY when a sample has an annotation-guided StringTie /
                          FLAIR GTF registered, samples.tsv stringtie_guided_gtf / flair_guided_gtf; drawn as
                          fig2g_guided) the samples tables' columns plus `mode` (annotation-guided), guided tools only.
                          Never in a panel with the annotation-free methods; no Rustle row (Rustle has no
                          annotation-guided transcript assembly; docs/archive/2026-09/PREREG_guided_transcript_comparison_2026-09-25.md).

Mode: every method of the main and samples tables is ANNOTATION-FREE (de novo): Rustle assemble (reads + genome),
StringTie -L without -G, FLAIR collapse without annotation (flair correct skipped), IsoSeq collapse; the figure says so
and prints "guided comparison: not available (guided StringTie/FLAIR GTFs not supplied)" until fig2_guided_* exist.

Methods: the SAME transcript GTFs as fig1 (assembly.ensure_arm), restricted to the same GENOME-WIDE scope (every
contig the sample's annotation covers; _sqanti.sqanti_contigs). The main tables hold every sample with the three lab
baselines (assembly.benchmark_samples: gorilla_OR6737, human_A119b), five methods each (cfg['fig2_arms'] narrows);
`species` there is the panel key (gorilla, human).

Heavy steps build() runs (foreground, serial, cached per contig under ${work}/fig2/<sample>/<contig>/; the genome and
annotation cuts are shared by samples on the same genome and annotation, _sqanti.reference_home):
  * the methods' GTFs (assembly.ensure_arm; never re-assembles: a stale assembly raises, run `make.py runs` first)
  * SQANTI3 QC + rules filter per sample x method x contig (x part of <= fig2_chunk_tx transcripts), in this order:
    Rustle's two configurations on every sample, then StringTie, FLAIR and IsoSeq collapse on the benchmark samples;
    an interrupted build resumes at the first uncached call. Every cached product carries a `.key` sidecar (source
    path, size, mtime; the SQANTI3 install and version) and is redone when it changes.
Options: `fig2_part` = all (default) | samples (only the supplementary tables: Rustle on every sample) | main |
guided (only the annotation-guided tables; `all` includes them when a guided GTF is registered);
`figs_budget_s` / `fig2_budget_s` bound one call (exit 75 = run again); `figs_plan=1` prints every SQANTI3 call still
to run with an estimated wall time, and runs nothing.
"""
from __future__ import annotations

import re
import sys
from pathlib import Path

import figlib
import _sqanti as S

META = {
    "id": "fig2",
    "title": "SQANTI3 structural categories and rules filter per method (the same GTFs as Figure 1, all "
             "annotation-free)",
    "claim": "Annotation-free (de novo) comparison: Rustle assemble (reads + genome), StringTie -L without -G, FLAIR "
             "collapse without annotation (flair correct skipped), IsoSeq collapse; guided comparison: not available "
             "(guided StringTie/FLAIR GTFs not supplied). "
             "SQANTI3 5.5.4 against the RefSeq annotation, on the Figure 1 GTFs (current tables: gorilla chr20, chr22 "
             "and chrY, i.e. NC_073244.2, NC_073246.2 and NC_073248.2, and human chr20-22; genome-wide after the "
             "rebuild): Rustle's two configurations have the highest full-splice-match (FSM) share and rules-filter "
             "PASS share of the five methods in both species (point estimates; gorilla FSM 39.8% default / 39.9% "
             "primary alignments only vs StringTie 36.6%, FLAIR 21.8%, IsoSeq collapse 12.1%; PASS 79.3 / 79.5% vs "
             "75.0, 59.8, 49.4%; human FSM 17.0 / 17.5% vs 13.5, 4.7, 3.2%; PASS 59.9 / 63.5% vs 48.3, 29.0, 27.8%), "
             "and the FSM order holds among multi-exon transcripts only (gorilla 39.9 / 39.9% vs 36.8, 22.6, 13.5%; "
             "human 17.3 / 18.7% vs 14.8, 6.4, 3.8%). IsoSeq collapse has the most FSM transcripts in both species "
             "(3,414 and 5,376 vs Rustle's 1,931 and 2,531) and FLAIR slightly more than Rustle in human (2,572). "
             "Completeness is mixed: more of Rustle's FSMs than StringTie's match both reference ends within 50 bp "
             "(gorilla 42.1 vs 35.0%, human 26.8 vs 22.5%; FLAIR 37.9% in human), but Rustle has more 3'-fragment "
             "ISMs, transcripts missing the reference's 5' exons (7.2 vs 3.7% of transcripts in gorilla, 5.8 vs 1.5% "
             "in human). Human FSM shares are low for every method; on multi-exon transcripts FSM counts the same "
             "thing as Figure 1's intron-chain precision, so the independent evidence is the intergenic share "
             "(9.6-20.9% of every human method vs 1.2-5.7% in gorilla): the gap lies in the library-annotation pair, "
             "not in one method. Supplementary (fig2s_samples, after the rebuild): Rustle's two configurations on all "
             "six samples.",
    "tables": ["fig2_sqanti_categories", "fig2_sqanti_filter", "fig2_sqanti_subcategories"],
    # drawn as fig2s_samples / fig2g_guided when present (built by the same `make.py data fig2`; the guided tables
    # only when an annotation-guided GTF is registered)
    "supplementary_tables": ["fig2_samples_categories", "fig2_samples_filter", "fig2_guided_categories",
                             "fig2_guided_filter"],
}
SPECIES = ["gorilla", "human"]


# SQANTI3 cost model for the plan (fit to the 36 recorded per-contig QC calls of 2026-09-25: 170-21,462 isoforms,
# 5.7-274 s; "SQANTI3 complete in X sec" of the qc logs): QC low/high = 5 s + 6 / 13 ms per isoform; the rules
# filter 3 s + 0.7 ms per isoform
QC_S = ((5.0, 0.006), (8.0, 0.013))
FILTER_S = (3.0, 0.0007)


class _Plan:
    def __init__(self):
        self.rows, self.lo, self.hi, self.n = [], 0.0, 0.0, 0

    def add(self, key, tool, contig, part, n_iso, what):
        lo = QC_S[0][0] + QC_S[0][1] * n_iso + FILTER_S[0] + FILTER_S[1] * n_iso
        hi = QC_S[1][0] + QC_S[1][1] * n_iso + FILTER_S[0] + FILTER_S[1] * n_iso
        self.rows.append((key, tool, contig, part, n_iso, what, lo, hi))
        self.lo, self.hi, self.n = self.lo + lo, self.hi + hi, self.n + n_iso


_TID_LINE = re.compile(r'transcript_id "([^"]*)"')


def count_by_contig(gtf, contigs) -> dict:
    """{contig: distinct transcript ids} of a GTF on `contigs` (plan mode; reads the file once)."""
    import gzip
    want = set(contigs)
    ids: dict = {c: set() for c in contigs}
    opener = gzip.open if str(gtf).endswith(".gz") else open
    with opener(gtf, "rt") as fh:
        for line in fh:
            c = line.split("\t", 1)[0]
            if c in want:
                m = _TID_LINE.search(line)
                if m:
                    ids[c].add(m.group(1))
    return {c: len(v) for c, v in ids.items()}


def build_sample(cfg: dict, key: str, tools: list[str], *, force=False, budget=None, plan=None):
    """SQANTI3 QC + rules filter of `tools` for one sample, one contig at a time -> {tool: merged tally}, inputs.
    With `plan` (a _Plan) nothing runs: the calls still to run are recorded with their isoform counts."""
    import assembly
    import samples

    contigs = S.sqanti_contigs(cfg, key)
    home = S.reference_home(cfg, key)
    root = figlib.work_dir(cfg, "fig2") / key
    href = figlib.work_dir(cfg, "fig2") / home
    cdir = {c: root / c for c in contigs}  # per-contig cache: valid for any contig set that contains the contig
    hdir = {c: href / c for c in contigs}
    row = samples.get(cfg, key)
    inputs, tallies_of = {}, {}
    if plan is None:
        for d in list(cdir.values()) + list(hdir.values()):
            d.mkdir(parents=True, exist_ok=True)
        genome = {c: S.subset_genome(row["fasta"], [c], hdir[c] / "genome.fa", force=force) for c in contigs}
        ref = S.sqanti_references(assembly.annotation_gtf(cfg, key), {c: hdir[c] / "ref.gtf" for c in contigs},
                                  force=force)
    else:
        genome = {c: hdir[c] / "genome.fa" for c in contigs}
        ref = {c: hdir[c] / "ref.gtf" for c in contigs}
    inputs[f"{key} genome"] = row["fasta"]
    inputs[f"{key} annotation"] = assembly.annotation_gtf(cfg, key, plan_only=plan is not None) or row["annotation_gff"]
    chunk = int(cfg.get("fig2_chunk_tx", "25000"))
    for tool in tools:
        if plan is not None:
            src = assembly.arm_source(cfg, key, tool, check=False)
            if not Path(src).exists():
                plan.rows.append((key, tool, "-", "-", 0, f"no GTF yet: {src}", 0.0, 0.0))
                continue
            counts = None
            for c in contigs:
                d = cdir[c]
                iso = d / f"{tool}.gtf"
                cached = iso.exists() and all(
                    S.qc_cached(cfg, part, ref[c], genome[c], d / f"qc_{tool}{sfx}", tool)
                    and S.filter_cached(cfg, d / f"qc_{tool}{sfx}" / f"{tool}_classification.txt",
                                        d / f"filt_{tool}{sfx}", tool)
                    for part, sfx in _parts_of(iso, chunk))
                if cached:
                    continue
                if counts is None:
                    counts = count_by_contig(src, contigs)
                n = counts.get(c, 0)
                if n:
                    k = -(-n // chunk)
                    per = -(-n // k)   # chunk_gtf's part size
                    for i in range(k):
                        plan.add(key, tool, c, i if k > 1 else "", min(per, n - i * per), "SQANTI3 QC + filter")
            continue
        # the methods are shared provisioning (fig1/fig3): --force recomputes SQANTI3, never the assemblies
        src = assembly.ensure_arm(cfg, key, tool, force=False)
        iso = S.split_gtf(src, {c: cdir[c] / f"{tool}.gtf" for c in contigs}, force=force)
        inputs[f"{key} {tool} gtf"] = src
        tallies = []
        for c in contigs:
            if S.n_transcripts(iso[c]) == 0:
                tallies.append(dict(S.EMPTY_TALLY))
                continue
            # large methods run in parts of <= fig2_chunk_tx transcripts (per-isoform results; parts sum exactly)
            parts = S.chunk_gtf(iso[c], chunk, force=force)
            for k, part in enumerate(parts):
                sfx = "" if len(parts) == 1 else f".part{k}"
                qdir, fdir = cdir[c] / f"qc_{tool}{sfx}", cdir[c] / f"filt_{tool}{sfx}"
                if budget is not None and (force or not S.qc_cached(cfg, part, ref[c], genome[c], qdir, tool)):
                    budget.check(f"SQANTI3 QC of {key} {tool} {c}{sfx}")
                qc = S.run_qc(cfg, part, ref[c], genome[c], qdir, tool, force=force)  # heavy
                if budget is not None and (force or not S.filter_cached(cfg, qc["classification"], fdir, tool)):
                    budget.check(f"SQANTI3 rules filter of {key} {tool} {c}{sfx}")
                filt = S.run_filter(cfg, qc, fdir, tool, force=force)
                tallies.append(S.tally(qc["classification"], filt))
                inputs[f"{key} {c} {tool}{sfx} classification"] = qc["classification"]
                inputs[f"{key} {c} {tool}{sfx} rules filter"] = filt
        tallies_of[tool] = S.merge_tallies(tallies)
    return tallies_of, inputs


def _parts_of(iso: Path, chunk: int):
    """[(part path, suffix)] of a per-contig GTF as chunk_gtf would cut it (plan mode: nothing is written)."""
    n = S.n_transcripts(iso)
    if n <= chunk:
        return [(iso, "")]
    k = -(-n // chunk)
    return [(iso.with_name(f"{iso.stem}.part{i}.gtf"), f".part{i}") for i in range(k)]


def scope_text(cfg: dict, key: str) -> str:
    return "genome-wide" if S.is_genome_wide(cfg, key) else ",".join(S.sqanti_contigs(cfg, key))


def work_order(cfg: dict, part: str) -> list[tuple[str, str]]:
    """(sample key, method) in the order the SQANTI3 calls run: Rustle's two configurations on every sample first
    (cheap; the supplementary figure), then the lab baselines of the benchmark samples."""
    import assembly
    main, everyone = assembly.benchmark_samples(cfg), assembly.all_samples(cfg)
    if part == "guided":
        return guided_order(cfg)
    first = [(k, t) for t in S.RUSTLE_ARMS for k in (everyone if part in ("all", "samples") else main)
             if t in S.arms(cfg, k)]
    rest = [(k, t) for t in S.DEFAULT_ARMS if t not in S.RUSTLE_ARMS for k in main if t in S.arms(cfg, k)] \
        if part in ("all", "main") else []
    return first + rest + (guided_order(cfg) if part == "all" else [])


def guided_order(cfg: dict) -> list[tuple[str, str]]:
    """(sample key, guided tool) of the annotation-guided path; [] until the user registers a guided GTF."""
    import assembly
    return [(k, t) for k in assembly.guided_samples(cfg) for t in assembly.guided_tools(cfg, k)]


def build(cfg: dict, data_dir: Path, force: bool = False):
    import assembly
    import samples

    part = (cfg.get("fig2_part") or "all").strip().lower()
    if part not in ("all", "samples", "main", "guided"):
        raise ValueError(f"fig2_part must be all, samples, main or guided, not {part!r}")
    print(f"[fig2] {assembly.guided_status(cfg)}", file=sys.stderr)
    order = work_order(cfg, part)
    if assembly.plan_only(cfg):
        plan = _Plan()
        for key in dict.fromkeys(k for k, _ in order):
            build_sample(cfg, key, [t for k, t in order if k == key], plan=plan)
        print("sample\tmethod\tcontig\tpart\tisoforms\tunit\test_s_low\test_s_high")
        for r in plan.rows:
            print("\t".join(f"{x:.0f}" if isinstance(x, float) else str(x) for x in r))
        print(f"TOTAL\t{len(plan.rows)} call(s)\t\t\t{plan.n}\t\t{plan.lo:.0f}\t{plan.hi:.0f}  "
              f"(= {plan.lo / 3600:.1f}-{plan.hi / 3600:.1f} h)")
        return
    budget = assembly.Budget(cfg, "fig2")
    tallies, inputs = {}, {}
    for key, tool in order:
        t, i = build_sample(cfg, key, [tool], force=force, budget=budget)
        tallies[(key, tool)] = t[tool]
        inputs.update(i)
    if part in ("all", "main"):
        cat_rows, filt_rows, sub_rows = [], [], []
        for key in assembly.benchmark_samples(cfg):
            scope = scope_text(cfg, key)
            for tool in S.arms(cfg, key):
                t = tallies[(key, tool)]
                cat_rows += S.category_rows(key, scope, tool, t)
                filt_rows.append(S.filter_row(key, scope, tool, t))
                sub_rows += S.subcategory_rows(key, scope, tool, t)
        write_tables(cfg, data_dir, cat_rows, filt_rows, sub_rows,
                     {k: v for k, v in inputs.items() if k.split(" ")[0] in assembly.benchmark_samples(cfg)})
    if part in ("all", "samples"):
        scat, sfilt = [], []
        for key in assembly.all_samples(cfg):
            row = samples.get(cfg, key)
            lead = [row["id"], assembly.sample_label(cfg, key), row["tissue"]]
            for tool in S.RUSTLE_ARMS:
                t = tallies[(key, tool)]
                scat += [lead + r for r in S.category_rows(row["species"], scope_text(cfg, key), tool, t)]
                sfilt.append(lead + S.filter_row(row["species"], scope_text(cfg, key), tool, t))
        write_samples_tables(cfg, data_dir, scat, sfilt, inputs)
    gorder = guided_order(cfg) if part in ("all", "guided") else []
    if gorder:
        gcat, gfilt = [], []
        for key, tool in gorder:
            row = samples.get(cfg, key)
            lead = [row["id"], assembly.sample_label(cfg, key), row["tissue"], assembly.MODE_GUIDED]
            t = tallies[(key, tool)]
            gcat += [lead + r for r in S.category_rows(row["species"], scope_text(cfg, key), tool, t)]
            gfilt.append(lead + S.filter_row(row["species"], scope_text(cfg, key), tool, t))
        write_guided_tables(cfg, data_dir, gcat, gfilt,
                            {k: v for k, v in inputs.items()
                             if any(re.search(rf"(^| ){t}(\.part\d+)?( |$)", k) for _, t in gorder)})


SAMPLES_LEAD = ["sample", "label", "tissue"]


def write_samples_tables(cfg, data_dir, cat_rows, filt_rows, inputs):
    import assembly
    notes = [_sqanti_note(cfg), "supplementary: Rustle's two configurations on every registry sample, each against its "
             "own annotation (chimpanzee and orangutan: RefSeq GFF3 converted with gff_to_gtf); mode: "
             + assembly.MODE_DENOVO + " (Rustle assemble, reads + genome); scope genome-wide = "
             "every contig the annotation covers (human: all but chrM); never pooled across samples or species"] \
        + [f"{k}: {S.reference_home(cfg, k)}'s per-contig genome and annotation cuts" for k in assembly.all_samples(cfg)]
    gen = "figures/fig_sqanti.py build"
    figlib.write_table("fig2_samples_categories", SAMPLES_LEAD + S.CATEGORY_HEADER, cat_rows, generator=gen,
                       inputs=inputs, notes=notes, data_dir=data_dir)
    figlib.write_table("fig2_samples_filter", SAMPLES_LEAD + S.FILTER_HEADER, filt_rows, generator=gen, inputs=inputs,
                       notes=notes + ["n_multiexon / n_fsm_multiexon: isoforms (FSM isoforms) with more than one exon"],
                       data_dir=data_dir)


GUIDED_LEAD = ["sample", "label", "tissue", "mode"]


def write_guided_tables(cfg, data_dir, cat_rows, filt_rows, inputs):
    import assembly
    notes = [_sqanti_note(cfg), "mode: " + assembly.MODE_GUIDED + " (StringTie -G / FLAIR with the annotation, as "
             "supplied in samples.tsv); " + assembly.GUIDED_CAVEAT, assembly.RUSTLE_NO_GUIDED,
             "pre-registered: docs/archive/2026-09/PREREG_guided_transcript_comparison_2026-09-25.md; each sample against its own "
             "annotation, genome-wide; never pooled across samples or species"]
    gen = "figures/fig_sqanti.py build"
    figlib.write_table("fig2_guided_categories", GUIDED_LEAD + S.CATEGORY_HEADER, cat_rows, generator=gen,
                       inputs=inputs, notes=notes, data_dir=data_dir)
    figlib.write_table("fig2_guided_filter", GUIDED_LEAD + S.FILTER_HEADER, filt_rows, generator=gen, inputs=inputs,
                       notes=notes + ["n_multiexon / n_fsm_multiexon: isoforms (FSM isoforms) with more than one exon"],
                       data_dir=data_dir)


def _sqanti_note(cfg) -> str:
    return ("SQANTI3 " + S.sqanti_version(cfg) + ": per contig, sqanti3_qc.py --isoforms METHOD --refGTF REF "
            "--refFasta GENOME --report skip -t N; sqanti3_filter.py rules --sqanti_class ... --filter_gtf "
            "..._corrected.gtf --skip_report (default rules JSON); per-contig (and per-part) tallies summed")


def write_tables(cfg, data_dir, cat_rows, filt_rows, sub_rows, inputs, extra_notes=()):
    import assembly
    notes = [_sqanti_note(cfg),
             "reference: empty gene_id replaced by the record's transcript_id (gene-level RefSeq records)",
             "scope: " + "; ".join(f"{k} = {scope_text(cfg, k)}" + (f" (every annotated contig; left out for every "
                                                                    f"method: {', '.join(assembly.unannotated_contigs(cfg, k))})"
                                                                    if assembly.unannotated_contigs(cfg, k) else "")
                                   for k in assembly.benchmark_samples(cfg)),
             "methods: " + ", ".join(S.arms(cfg)) + " (the fig1 GTFs, assembly.ensure_arm)",
             "mode: " + assembly.MODE_DENOVO_METHODS + "; " + assembly.guided_status(cfg)] + list(extra_notes)
    gen = "figures/fig_sqanti.py build"
    figlib.write_table("fig2_sqanti_categories", S.CATEGORY_HEADER, cat_rows, generator=gen, inputs=inputs,
                       notes=notes, data_dir=data_dir)
    figlib.write_table("fig2_sqanti_filter", S.FILTER_HEADER, filt_rows, generator=gen, inputs=inputs,
                       notes=notes + ["n_multiexon / n_fsm_multiexon: isoforms (FSM isoforms) with more than one exon "
                                      "(classification column `exons`); fsm_multiexon_frac = their ratio"],
                       data_dir=data_dir)
    figlib.write_table("fig2_sqanti_subcategories", S.SUBCATEGORY_HEADER, sub_rows, generator=gen, inputs=inputs,
                       notes=notes + ["SQANTI3 `subcategory` per structural category (FSM: reference_match = both "
                                      "ends within 50 bp of the reference transcript's; ISM: 3prime_fragment = "
                                      "missing the reference's 5' exons); n_multiexon = those with > 1 exon; "
                                      "frac = n / n_total (all isoforms of the arm)"],
                       data_dir=data_dir)


# ------------------------------------------------------------------------------------------------ plotting
_TEXT_ON = {"FSM": "#ffffff", "ISM": "#ffffff", "NIC": figlib.INK, "NNC": figlib.INK, "Other": figlib.INK}
BAR_H = 0.60        # bar height (row pitch 1.0): the gap between bars holds the call-outs of small segments
ROW_IN = 0.33       # inches per arm
SEG_FS = 5.6        # segment label size (>= 5.5 pt)


# the gorilla SQANTI3 contigs by chromosome (region lines of GGO_genomic.gff, RefSeq GCF_029281585.2): the species
# rows are not matched chromosomes, and gorilla includes chrY (ampliconic multi-copy genes, testis library)
CHROM_NAME = {"NC_073244.2": "chr20", "NC_073246.2": "chr22", "NC_073248.2": "chrY"}


def _scope_text(scope: str) -> str:
    if scope in ("genome-wide", ""):
        return scope
    contigs = [c for c in scope.split(",") if c]
    if contigs and all(c in CHROM_NAME for c in contigs):
        return ", ".join(CHROM_NAME[c] for c in contigs) + " (" + ", ".join(contigs) + ")"
    return ", ".join(contigs)


def _pct(w: float) -> str:
    return f"{100 * w:.1f}%" if 0 < w < 0.01 else f"{100 * w:.0f}%"


def _text_size(ax, text, fontsize, linespacing=1.0):
    r = ax.figure.canvas.get_renderer()
    probe = ax.text(0, 0, text, fontsize=fontsize, linespacing=linespacing)
    bb = probe.get_window_extent(r)
    probe.remove()
    return bb.width, bb.height


def label_segments(ax, bars):
    """Name every segment of the 100% stacked bars ("NNC 23%"): one line inside when it fits, two lines ("NNC" over
    "23%") when only that fits, else a call-out above the bar joined by a leader. `bars` = [(y, [(cat, left, w)])]
    in data units (x 0..1). NIC and NNC are never identified by colour alone."""
    fig = ax.figure
    pt = fig.dpi / 72.0
    to_disp = ax.transData.transform
    bar_px = abs(to_disp((0, BAR_H / 2))[1] - to_disp((0, -BAR_H / 2))[1])
    unit_px = abs(to_disp((1, 0))[0] - to_disp((0, 0))[0])
    for y, segs in bars:
        callouts = []
        for cat, left, w in segs:
            if w <= 0:
                continue
            seg_px = w * unit_px
            one = f"{cat} {_pct(w)}"
            w1, h1 = _text_size(ax, one, SEG_FS)
            if w1 + 3 * pt <= seg_px and h1 <= bar_px:
                ax.text(left + w / 2, y, one, ha="center", va="center", fontsize=SEG_FS, color=_TEXT_ON[cat])
                continue
            two = f"{cat}\n{_pct(w)}"
            w2, h2 = _text_size(ax, two, SEG_FS, linespacing=0.9)
            if w2 + 2 * pt <= seg_px and h2 <= bar_px + 1.0 * pt:
                ax.text(left + w / 2, y, two, ha="center", va="center", fontsize=SEG_FS, color=_TEXT_ON[cat],
                        linespacing=0.9)
                continue
            callouts.append((left + w / 2, one, w1 / unit_px))
        # call-outs above the bar, spread left-to-right so none overlap, kept within [0, 1.02]
        pad = 3 * pt / unit_px
        xs = []
        for xc, text, tw in sorted(callouts):
            lo = xc - tw / 2
            if xs:
                lo = max(lo, xs[-1][1] + pad)
            xs.append((lo, lo + tw))
        if xs and xs[-1][1] > 1.02:
            shift = xs[-1][1] - 1.02
            xs = [(a - shift, b - shift) for a, b in xs]
        ty = y + BAR_H / 2 + 3.2 * pt / abs(to_disp((0, 1))[1] - to_disp((0, 0))[1])
        for (xc, text, tw), (a, b) in zip(sorted(callouts), xs):
            ax.annotate(text, xy=(xc, y + BAR_H / 2), xytext=((a + b) / 2, ty), textcoords="data",
                        ha="center", va="bottom", fontsize=SEG_FS, color=figlib.INK, annotation_clip=False,
                        arrowprops=dict(arrowstyle="-", color=figlib.INK_3, lw=0.4, shrinkA=0, shrinkB=0))


def _not_built(fig, spec, species: str, letters: str):
    """A clearly marked 'not yet built' band in place of a species' panels (no empty axes, no stray letters)."""
    from matplotlib.patches import FancyBboxPatch

    ax = fig.add_subplot(spec)
    ax.axis("off")
    ax.add_patch(FancyBboxPatch((0.0, 0.08), 1.0, 0.84, boxstyle="round,pad=0,rounding_size=0.015",
                                transform=ax.transAxes, facecolor="#f7f6f2", edgecolor=figlib.INK_3, lw=0.6,
                                linestyle=(0, (3, 2)), clip_on=False))
    ax.text(0.015, 0.62, f"{letters}   {figlib.SPECIES_LABEL[species]}: SQANTI3 not yet built", fontsize=7.0,
            fontweight="bold", color=figlib.INK_2, ha="left", va="center", transform=ax.transAxes)
    ax.text(0.015, 0.30, "python3 figures/make.py data fig2 fills this row: the same five methods as Figure 1, "
            "genome-wide (see caption)", fontsize=6.0, color=figlib.INK_2, ha="left", va="center",
            transform=ax.transAxes)


def _species_row(fig, spec, species, cats, filt, letters):
    fr = {r["tool"]: r for r in filt if r["species"] == species}
    tools = [t for t in figlib.TOOL_ORDER if t in fr]
    entries = [(t, f"{figlib.TOOL_LABEL[t]}\nn = {int(fr[t]['n_total']):,}", fr[t],
                [r for r in cats if r["species"] == species and r["tool"] == t]) for t in tools]
    return _bar_rows(fig, spec, entries, letters)


def _bar_rows(fig, spec, entries, letters, label_fs=6.2):
    """One row of three panels: (a) 100% stacked structural categories, (b) rules-filter PASS share, (c) FSM
    transcripts; one bar per entry = (method key for the bar style, tick label, filter row, category rows)."""
    import matplotlib.pyplot as plt
    from matplotlib.ticker import PercentFormatter

    gs = spec.subgridspec(1, 3, width_ratios=[2.7, 1.0, 1.0], wspace=0.20)
    ax_c = fig.add_subplot(gs[0])
    ax_p = fig.add_subplot(gs[1], sharey=ax_c)
    ax_f = fig.add_subplot(gs[2], sharey=ax_c)
    for ax, letter in zip((ax_c, ax_p, ax_f), letters):
        figlib.panel_label(ax, letter, x=-0.24 if ax is ax_c else -0.10, y=1.02)
    n = len(entries)
    yy = list(range(n))[::-1]
    # (a/d) 100% stacked structural categories, every segment named
    bars = []
    for (t, _, fr, crs), y in zip(entries, yy):
        n_total = int(fr["n_total"])
        fold = {c: 0 for c in figlib.SQANTI_ORDER}
        for r in crs:
            fold[r["category"]] += int(r["n"])
        left, segs = 0.0, []
        for c in figlib.SQANTI_ORDER:
            w = fold[c] / n_total if n_total else 0.0
            ax_c.barh(y, w, left=left, height=BAR_H, color=figlib.SQANTI_COLOR[c], edgecolor=figlib.SURFACE,
                      linewidth=0.6)
            segs.append((c, left, w))
            left += w
        bars.append((y, segs))
    ax_c.set_xlim(0, 1)
    ax_c.set_ylim(-0.55, n - 1 + 0.95)   # headroom above the top bar for its call-outs
    ax_c.xaxis.set_major_formatter(PercentFormatter(1.0, decimals=0))
    ax_c.set_yticks(yy)
    ax_c.set_yticklabels([lab for _, lab, _, _ in entries], fontsize=label_fs)
    ax_c.tick_params(axis="y", length=0)
    ax_c.grid(False)
    ax_c.set_xlabel("Share of transcripts (SQANTI3 structural category)")
    for sp in ("left",):
        ax_c.spines[sp].set_visible(False)
    label_segments(ax_c, bars)
    # (b/e) rules-filter PASS share
    for (t, _, fr, _), y in zip(entries, yy):
        if fr["pass_frac"] == "":
            continue
        v = float(fr["pass_frac"])
        ax_p.barh(y, v, height=BAR_H, **figlib.tool_bar_kwargs(t))
        ax_p.text(v + 0.03, y, f"{100 * v:.1f}%", ha="left", va="center", fontsize=5.8, color=figlib.INK_2)
    ax_p.set_xlim(0, 1.38)
    ax_p.set_xticks([0, 0.5, 1.0])
    ax_p.xaxis.set_major_formatter(PercentFormatter(1.0, decimals=0))
    ax_p.set_xlabel("Passes SQANTI3's rules filter")
    # (c/f) full-splice-match transcripts
    fsm = [int(fr["n_fsm"]) for _, _, fr, _ in entries]
    top = max(fsm) or 1
    for (t, _, fr, _), y, v in zip(entries, yy, fsm):
        ax_f.barh(y, v, height=BAR_H, **figlib.tool_bar_kwargs(t))
        ax_f.text(v + top * 0.04, y, f"{v:,}", ha="left", va="center", fontsize=5.8, color=figlib.INK_2)
    ax_f.set_xlim(0, top * 1.5)
    ax_f.xaxis.set_major_locator(plt.MaxNLocator(3, integer=True))
    ax_f.xaxis.set_major_formatter(lambda v, _: (f"{v / 1000:g}k" if v >= 1000 else f"{v:g}") if v else "0")
    ax_f.set_xlabel("FSM transcripts")
    for ax in (ax_p, ax_f):
        ax.tick_params(axis="y", left=False, labelleft=False)
        ax.grid(axis="y", visible=False)
        ax.grid(axis="x")
    return ax_c


def plot(data_dir: Path, out_dir: Path):
    import matplotlib.pyplot as plt
    from matplotlib.patches import Patch

    figlib.use_style()
    cats = figlib.read_table("fig2_sqanti_categories", data_dir)
    filt = figlib.read_table("fig2_sqanti_filter", data_dir)
    n_arms = {sp: sum(1 for r in filt if r["species"] == sp) for sp in SPECIES}
    heights = [ROW_IN * n_arms[sp] + 0.62 if n_arms[sp] else 0.42 for sp in SPECIES]
    top_in = 0.74 if n_arms[SPECIES[0]] else 0.40   # a built first row needs room for its species header
    bottom_in = 0.46                                  # keeps the PROVISIONAL line clear of the x labels
    gap_in = 0.88 if n_arms[SPECIES[0]] else 0.48   # x ticks + x label of the row above, header, letters
    fig_h = top_in + bottom_in + sum(heights) + gap_in * (len(SPECIES) - 1)
    fig = plt.figure(figsize=(figlib.WIDTH_DOUBLE - 0.16, fig_h))
    gs = fig.add_gridspec(len(SPECIES), 1, height_ratios=heights, hspace=gap_in / (sum(heights) / len(heights)),
                          left=0.175, right=0.975, top=1 - top_in / fig_h, bottom=bottom_in / fig_h)
    letters = ["abc", "def"]
    for row, species in enumerate(SPECIES):
        scope = next((r["scope"] for r in filt if r["species"] == species), "")
        if not n_arms[species]:
            _not_built(fig, gs[row], species, "a–c" if species == "gorilla" else "d–f")
            continue
        ax_c = _species_row(fig, gs[row], species, cats, filt, letters[row])
        y_top = ax_c.get_position().y1
        fig.text(0.01, y_top + 0.20 / fig_h, figlib.SPECIES_LABEL[species] + f"  ·  {_scope_text(scope)}",
                 fontsize=7.5, fontweight="bold", ha="left", va="bottom")
    handles = [Patch(facecolor=figlib.SQANTI_COLOR[c], edgecolor="none",
                     label=f"{c}: {figlib.SQANTI_LABEL[c].lower()}") for c in figlib.SQANTI_ORDER]
    fig.legend(handles=handles, loc="upper left", bbox_to_anchor=(0.175, 1 - 0.04 / fig_h), ncol=3, fontsize=6.0,
               handlelength=1.1, handleheight=0.9, columnspacing=1.2, borderaxespad=0.0)
    # the comparison mode, stated on the figure (every row above compares annotation-free methods only)
    fig.text(0.01, 1 + 0.10 / fig_h, figlib.mode_lines("fig2g_guided", (Path(data_dir) / "fig2_guided_filter.tsv")
                                                        .exists()),
             fontsize=6.0, color=figlib.INK_2, ha="left", va="bottom", linespacing=1.3)
    figlib.stamp_provisional(fig, META["tables"], data_dir)
    paths = figlib.save(fig, "fig2_sqanti", out_dir)
    plt.close(fig)
    return paths + plot_samples(data_dir, out_dir) + plot_guided(data_dir, out_dir)


def plot_guided(data_dir: Path, out_dir: Path) -> list:
    """fig2g_guided: SQANTI3 categories, PASS and FSM of the annotation-guided StringTie/FLAIR runs (one bar per
    sample x guided tool; never with the annotation-free methods; no Rustle row). Only when the tables exist."""
    import matplotlib.pyplot as plt
    from matplotlib.patches import Patch

    try:
        cats = figlib.read_table("fig2_guided_categories", data_dir)
        filt = figlib.read_table("fig2_guided_filter", data_dir)
    except FileNotFoundError:
        return []
    entries = [(r["tool"], f"{r['label']}\n{figlib.TOOL_LABEL.get(r['tool'], r['tool'])}, n = {int(r['n_total']):,}",
                r, [c for c in cats if c["sample"] == r["sample"] and c["tool"] == r["tool"]]) for r in filt]
    n = len(entries)
    top_in, bottom_in = 0.95, 0.62
    fig_h = top_in + 0.34 * n + 0.2 + bottom_in
    fig = plt.figure(figsize=(figlib.WIDTH_DOUBLE - 0.16, fig_h))
    gs = fig.add_gridspec(1, 1, left=0.26, right=0.975, top=1 - top_in / fig_h, bottom=bottom_in / fig_h)
    _bar_rows(fig, gs[0], entries, "abc", label_fs=5.9)
    handles = [Patch(facecolor=figlib.SQANTI_COLOR[c], edgecolor="none",
                     label=f"{c}: {figlib.SQANTI_LABEL[c].lower()}") for c in figlib.SQANTI_ORDER]
    fig.legend(handles=handles, loc="upper left", bbox_to_anchor=(0.26, 1 - 0.36 / fig_h), ncol=3, fontsize=6.0,
               handlelength=1.1, handleheight=0.9, columnspacing=1.2, borderaxespad=0.0)
    fig.text(0.01, 1 - 0.04 / fig_h, "Annotation-guided runs only (never compared with the annotation-free methods of "
             "Fig. 2). The guided tools were given the annotation they are scored against.\nRustle has no "
             "annotation-guided transcript assembly, so it has no bar here.", fontsize=5.8, color=figlib.INK_2,
             ha="left", va="top")
    scope = sorted({r["scope"] for r in filt})
    fig.text(0.01, 0.16 / fig_h, "Each sample against its own RefSeq annotation (" + ", ".join(scope) + "). Never "
             "pooled across samples.", fontsize=5.8, color=figlib.INK_2, ha="left", va="bottom")
    figlib.stamp_provisional(fig, ["fig2_guided_categories", "fig2_guided_filter"], data_dir)
    paths = figlib.save(fig, "fig2g_guided", out_dir)
    plt.close(fig)
    return paths


def plot_samples(data_dir: Path, out_dir: Path) -> list:
    """fig2s_samples: SQANTI3 categories, rules-filter PASS and FSM transcripts of Rustle's two configurations on
    every sample (two bars per sample, grouped; never pooled). Drawn only when the fig2_samples tables exist."""
    import matplotlib.pyplot as plt
    from matplotlib.patches import Patch

    try:
        cats = figlib.read_table("fig2_samples_categories", data_dir)
        filt = figlib.read_table("fig2_samples_filter", data_dir)
    except FileNotFoundError:
        return []
    fr = {(r["sample"], r["tool"]): r for r in filt}
    order = list(dict.fromkeys(r["sample"] for r in filt))
    entries = []
    for sid in order:
        for t in ("rustle", "rustle_primary"):
            r = fr.get((sid, t))
            if not r:
                continue
            lab = (f"{r['label']}\nRustle, n = {int(r['n_total']):,}" if t == "rustle" else
                   f"Rustle, primary alignments only\nn = {int(r['n_total']):,}")
            entries.append((t, lab, r, [c for c in cats if c["sample"] == sid and c["tool"] == t]))
    n = len(entries)
    top_in, bottom_in = 0.62, 0.62
    body = 0.34 * n + 0.2
    fig_h = top_in + body + bottom_in
    fig = plt.figure(figsize=(figlib.WIDTH_DOUBLE - 0.16, fig_h))
    gs = fig.add_gridspec(1, 1, left=0.26, right=0.975, top=1 - top_in / fig_h, bottom=bottom_in / fig_h)
    ax_c = _bar_rows(fig, gs[0], entries, "abc", label_fs=5.9)
    for k in range(2, n, 2):   # a rule between samples, under the tick labels only (clear of the call-outs)
        ax_c.plot([-0.45, -0.01], [n - k - 0.5, n - k - 0.5], color=figlib.GRID, lw=0.6, clip_on=False,
                  transform=ax_c.get_yaxis_transform())
    handles = [Patch(facecolor=figlib.SQANTI_COLOR[c], edgecolor="none",
                     label=f"{c}: {figlib.SQANTI_LABEL[c].lower()}") for c in figlib.SQANTI_ORDER]
    fig.legend(handles=handles, loc="upper left", bbox_to_anchor=(0.26, 1 - 0.04 / fig_h), ncol=3, fontsize=6.0,
               handlelength=1.1, handleheight=0.9, columnspacing=1.2, borderaxespad=0.0)
    scope = sorted({r["scope"] for r in filt})
    fig.text(0.01, 0.16 / fig_h, "Each sample against its own RefSeq annotation (" + ", ".join(scope)
             + "); hatched: Rustle, primary alignments only. Never pooled across samples.", fontsize=5.8,
             color=figlib.INK_2, ha="left", va="bottom")
    figlib.stamp_provisional(fig, META["supplementary_tables"], data_dir)
    paths = figlib.save(fig, "fig2s_samples", out_dir)
    plt.close(fig)
    return paths


# ------------------------------------------------------------------------------------------------ caption numbers
def caption_numbers(data_dir=None) -> str:
    """Every number the caption quotes, read from the tables (`python3 figures/fig_sqanti.py summary`; light)."""
    data_dir = Path(data_dir or figlib.DATA_DIR)
    out = []
    cats = figlib.read_table("fig2_sqanti_categories", data_dir)
    filt = figlib.read_table("fig2_sqanti_filter", data_dir)
    try:
        subc = figlib.read_table("fig2_sqanti_subcategories", data_dir)
    except FileNotFoundError:
        subc = []

    def pct(a, b):
        return f"{100 * a / b:.1f}%" if b else "-"

    def block(tag, frows, crows, srows, key):
        for r in frows:
            k = r[key]
            n = int(r["n_total"])
            fold = {c: 0 for c in figlib.SQANTI_ORDER}
            fine = {}
            for c in crows:
                if c[key] == k and c["tool"] == r["tool"]:
                    fold[c["category"]] += int(c["n"])
                    fine[c["structural_category"]] = int(c["n"])
            sub = {(s["structural_category"], s["subcategory"]): int(s["n"]) for s in srows
                   if s[key] == k and s["tool"] == r["tool"]}
            n_fsm = int(r["n_fsm"])
            out.append(f"{tag} {k} ({r['scope']}) {r['tool']:>15}: n {n:,}  FSM {pct(n_fsm, n)} ({n_fsm:,})  PASS "
                       + (f"{100 * float(r['pass_frac']):.1f}%" if r["pass_frac"] else "-")
                       + f"  FSM multi-exon {100 * float(r['fsm_multiexon_frac']):.1f}%  "
                       + "  ".join(f"{c} {pct(v, n)}" for c, v in fold.items())
                       + f"  intergenic {pct(fine.get('intergenic', 0), n)}  genic_intron "
                         f"{pct(fine.get('genic_intron', 0), n)}"
                       + (f"  FSM reference_match {pct(sub.get(('full-splice_match', 'reference_match'), 0), n_fsm)}"
                          f" of FSM  ISM 3prime_fragment "
                          f"{pct(sub.get(('incomplete-splice_match', '3prime_fragment'), 0), n)}" if sub else ""))

    block("main", filt, cats, subc, "species")
    try:
        sf = figlib.read_table("fig2_samples_filter", data_dir)
        sc = figlib.read_table("fig2_samples_categories", data_dir)
        block("sample", sf, sc, [], "sample")
    except FileNotFoundError:
        pass
    try:
        gf = figlib.read_table("fig2_guided_filter", data_dir)
        gc = figlib.read_table("fig2_guided_categories", data_dir)
        block("guided (annotation-guided; given the annotation it is scored against; no Rustle counterpart)", gf, gc,
              [], "sample")
    except FileNotFoundError:
        out.append("guided: not available (guided StringTie/FLAIR GTFs not supplied)")
    return "\n".join(out)


def main(argv=None):
    import argparse
    ap = argparse.ArgumentParser(description="fig2: caption numbers from the tables (light)")
    sub = ap.add_subparsers(dest="cmd", required=True)
    sm = sub.add_parser("summary", help="print every number the caption quotes, from figures/data")
    sm.add_argument("--data", default=str(figlib.DATA_DIR))
    a = ap.parse_args(argv)
    if a.cmd == "summary":
        print(caption_numbers(Path(a.data)))


if __name__ == "__main__":
    main()
