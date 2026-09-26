"""Figure 6 — the default de novo family definition across the paralogue identity spectrum.

The families are Rustle's ONE default de novo family definition (user decision 2026-09-25; docs/
PREREG_genome_wide_families_2026-09-25.md, Amendment 1): reads -> seeded assembly loci -> one representative per locus
(its "positional exon sum") -> families (the driver's `families` stage: mcl_families --from-gtf, exon-sum >= 0.60,
MCL 2.8); its copy table (<id>.fam.copies.tsv) is what copy assignment consumes. The reference is external: Ensembl
Compara. No protein enters the main figure.

Genome-wide version (default), human A119b and human testis, whole genome:
(a) Ensembl Compara (release 116) paralogue pairs, cross-chromosome pairs kept, whose two genes both have a locus in
    the sample's genome-wide assembly; sensitivity by the pair's Compara protein-identity band, for a DIRECT nucleotide
    alignment between the two genes' locus representatives (score.py spectrum --skip-t3: minimap2 asm20 and asm20
    -k11 -w5) and for the default families (both genes in one family; split by whether their loci also align
    directly); dashed outline = both genes are family members (the most the families can recover).
(b) Per human sample: 60-90% pooled, the largest paralogue group vs the others; >= 90%, pairs across chromosomes vs
    pairs on one chromosome.
(c) Per human sample: precision against Compara over judgeable pairs (both genes have Compara data): direct
    alignments, and the families' within-family gene pairs (all copies, multi-exon copies), with the value without the
    largest family.
Development version (`fig6_scope dev`): the same panels on human A119b chr16 (tables fig6_chr16_*), drawn when the
genome-wide tables are absent.

Supplementary figures (same module; none is needed by the main figure):
  fig6s_seeding   seeding loci with secondary alignments vs primary alignments only, each through the families stage,
                  on external references: Compara pairs at >= 90% (human samples; genes with >= 2 primary reads on
                  their exons) and Liftoff copy pairs (every sample; both loci read-supported); protein-homology
                  families only as a secondary reference when built (`fig6s_protein_homology 1`). Development version:
                  gorilla OR6737 NC_073244.2 against that contig's protein-homology families (tables fig6_gorilla_*).
  fig6s_protein   the translated protein search (T3) as a comparator of the direct alignment, human chr16
                  (development spectrum): not part of Rustle's rule, never run genome-wide (table fig6s_protein_tiers).

Scorers: `bench/score.py spectrum` and `bench/score.py pairs` (commands in figures/_o1.py). Species are never pooled.
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import figlib
import _o1

GW_TABLES = ["fig6_gw_recall", "fig6_gw_groups", "fig6_gw_precision"]
DEV_TABLES = ["fig6_chr16_recall", "fig6_chr16_groups", "fig6_chr16_precision"]
SUPP_TABLES = ["fig6s_seeding", "fig6s_protein_tiers", "fig6_gorilla_recall", "fig6_gorilla_clusters"]


def _gw_ready(data_dir: Path) -> bool:
    return all((Path(data_dir) / f"{t}.tsv").exists() for t in GW_TABLES)


def _have(data_dir: Path, *tables) -> bool:
    return all((Path(data_dir) / f"{t}.tsv").exists() for t in tables)


DEV_CLAIM = (
    "Development tables (human A119b chr16, the development chromosome): Ensembl Compara pairs whose two genes both "
    "have a locus in Rustle's chr16 assembly, against Rustle's default de novo families (reads -> seeded loci -> one "
    "representative per locus -> families; the families stage on the chr16 restriction of the genome-wide assembly). "
    "The numbers are printed by `python3 figures/fig_family_spectrum.py summary` and quoted in captions/fig6.md. The "
    "genome-wide version (human A119b and human testis) replaces these tables when it is built.")
GW_CLAIM = (
    "Genome-wide, human A119b and human testis (pre-registered claims F6.1-F6.5, "
    "docs/PREREG_genome_wide_families_2026-09-25.md, Amendment 1; the numbers are printed by `python3 "
    "figures/fig_family_spectrum.py summary` and quoted in captions/fig6.md): Ensembl Compara pairs across "
    "chromosomes, sensitivity by protein-identity band, a direct nucleotide alignment vs Rustle's default de novo "
    "families; precision against Compara.")

META = {
    "id": "fig6",
    "title": "The default de novo family definition across the paralogue identity spectrum",
    "claim": GW_CLAIM if _gw_ready(figlib.DATA_DIR) else DEV_CLAIM,
    "tables": GW_TABLES if _gw_ready(figlib.DATA_DIR) else DEV_TABLES,
    # drawn as fig6s_seeding / fig6s_protein when present (none is needed by the main figure)
    "supplementary_tables": SUPP_TABLES,
}
PROVISIONAL = ("provisional: from the recorded 2026-09-22/24 runs (registers 1096-1101); make.py data regenerates "
               "with current defaults")
GEN = "figures/fig_family_spectrum.py"
GW_UNIT_VERSION = "2"    # bump when a cached genome-wide scoring unit's rows change (2: default families, no catalog)

BAND_LABEL = {">=90": "≥90", "80-90": "80–90", "70-80": "70–80", "60-70": "60–70", "50-60": "50–60",
              "30-50": "30–50", "<30": "<30", "<60": "<60", "none": "none", "60-90": "60–90"}
# colours: the direct alignment is NOT a tool, so it wears ink greys; Rustle's families wear Rustle blue. No light
# blue anywhere in panels a-c (light hatched blue = Rustle primaries only, the seeding supplement).
EDGE_NT = figlib.INK_2                 # a direct nucleotide alignment between the two genes' locus representatives
EDGE_PROT = "#bdbcb6"                  # supplement only: the pairs the translated protein search adds
FAMILY = figlib.TOOL_COLOR["rustle"]   # same default family, and the two genes' loci also align directly
FAMILY_VIA = figlib.BLUE[700]          # same default family, no direct alignment (joined through other loci)
CONFIG_LABEL = {"rustle_primary": "primary alignments only",
                "rustle": "default (+ secondary alignments within 2% of the best score)"}
FAMILIES_WHAT = ("default de novo families: families stage (mcl_families --from-gtf --min-exonic-bp 1 "
                 "--min-shared-exon-frac 0.60, MCL 2.8) on the sample's genome-wide assembly (loci seeded with "
                 "secondary alignments within 2% of the read's best score); one copy per member locus = its "
                 "representative transcript (most reads, tie longer span), the copy table <id>.fam.copies.tsv")


# ================================================================ data
def build(cfg: dict, data_dir: Path, force: bool = False, recorded: bool = False):
    """Genome scope (default): the genome-wide tables (HEAVY units resumable, exit 75 = call again). Dev scope
    (`--set fig6_scope=dev`): the chr16 development tables and the development supplements. recorded=True: the
    supplementary development tables from the recorded runs (no recorded run of the default families exists)."""
    if recorded:
        return build_recorded(cfg, data_dir)
    if _o1.scope(cfg) == "dev":
        return build_dev(cfg, data_dir, force)
    return build_genome(cfg, data_dir, force)


def _flag(cfg: dict, key: str) -> bool:
    return str(cfg.get(key, "")).strip().lower() in ("1", "true", "yes")


# ---------------------------------------------------------------- genome scope
def _unit(w: Path, name: str, key, fn, budget, est_s: float):
    """A cached scoring unit: the rows `fn()` returns, recomputed only when `key` (inputs) changes."""
    p = w / f"{name}.json"
    k = json.dumps([GW_UNIT_VERSION, key], sort_keys=True)
    if p.exists():
        j = json.loads(p.read_text())
        if j.get("key") == k:
            return j["value"]
    _o1.need(budget, est_s, name)
    v = fn()
    p.write_text(json.dumps({"key": k, "value": v}))
    return v


def _fp(*paths) -> list:
    return [figlib.file_fingerprint(p) for p in paths]


def build_genome(cfg: dict, data_dir: Path, force: bool = False):
    budget = _o1.gw_budget(cfg, "fig6")
    w = figlib.work_dir(cfg, _o1.FIG) / "gw_score"
    w.mkdir(parents=True, exist_ok=True)
    by_sp = _o1.samples_by_species(cfg)
    missing: list = []
    score_py, lib_py = Path(cfg["repo"]) / "bench" / "score.py", Path(cfg["repo"]) / "bench" / "lib.py"

    def miss(e):
        if str(e) not in missing:
            missing.append(str(e))
            print(f"[fig6] not built: {e}", file=sys.stderr)

    # ---- heavy units first: a call that runs out of budget stops (exit 75) before any table is written
    compara = None
    try:
        compara = _o1.compara_gw(cfg)
    except _o1.NotBuilt as e:
        miss(e)
    human_src = {}
    for sid in by_sp.get("human", []) if compara else []:
        try:
            copies = _o1.families_copies(cfg, sid, "families")
            human_src[sid] = (_o1.spectrum_gw(cfg, sid, budget), copies)
        except (_o1.NotBuilt, RuntimeError) as e:
            if isinstance(e, RuntimeError) and not isinstance(e, _o1.NotBuilt) and "make.py runs" not in str(e):
                raise
            miss(e)

    # ---- main figure, panels a-c: human samples (cached per sample)
    rec_rows, grp_rows, prec_rows, inputs = [], [], [], {"score_py": score_py, "lib_py": lib_py}
    for sid, (prefix, copies) in human_src.items():
        ann = _o1.annotation_cache(cfg, "human")
        inputs.update({f"spectrum_{sid}": f"{prefix}.spectrum.tsv", f"truth_pairs_{sid}": f"{prefix}.truth_pairs.tsv",
                       f"families_copies_{sid}": copies})
        rows = _unit(w, f"compara_{sid}", _fp(f"{prefix}.spectrum.tsv", f"{prefix}.truth_pairs.tsv", copies, compara,
                                             ann["exons_gtf"], score_py, lib_py),
                     lambda sid=sid, prefix=prefix, copies=copies: _compara_rows(cfg, sid, prefix, copies, compara, w),
                     budget, 240)
        rec_rows += rows["recall"]
        grp_rows += rows["groups"]
        prec_rows += rows["precision"]
    inputs["compara_gw"] = compara or ""

    # ---- supplement fig6s_seeding: every sample, both seeding configurations, external references
    seed_rows, seed_inputs = _seeding_gw(cfg, by_sp, compara, budget, w, miss)

    notes = ["genome-wide, pre-registered in docs/PREREG_genome_wide_families_2026-09-25.md (Amendment 1); "
             "generated by figures/fig_family_spectrum.py build_genome"]
    if missing:
        notes.append("provisional: not every sample is built yet — " + " | ".join(missing))
    uni_note = ("pairs scored: Ensembl Compara (release 116) human paralogue pairs, cross-chromosome pairs kept, whose "
                "two genes both have a locus (representative >= 200 bp) in the sample's genome-wide assembly "
                "(score.py spectrum --chrom ALL --skip-t3); conditioned on the assembler's loci, not on read counts")
    fam_note = "families = " + FAMILIES_WHAT + "; score.py pairs --chrom ALL --universe (copies mapped to genes by span)"
    if rec_rows:
        figlib.write_table(
            "fig6_gw_recall",
            ["sample", "sample_label", "species", "substrate", "reference", "compara_band", "n_pairs", "view", "hits",
             "sensitivity", "ci_lo", "ci_hi"], rec_rows, generator=GEN, inputs=inputs, data_dir=data_dir,
            notes=notes + [uni_note, fam_note,
                           "substrate S0 = whole genome; S1 = pairs with no gene on a development contig (human chr16)",
                           "view edge_asm20 = minimap2 asm20, identity >= 0.80; edge_nt = asm20 or asm20 -k11 -w5 "
                           "(identity >= 0.60); coverage >= 0.50 (query span / shorter length); a DIRECT nucleotide "
                           "alignment between any loci of the two genes (one representative per locus)",
                           "view families = both genes in one default family; families_edge / families_no_edge split it "
                           "by the direct alignment (edge_nt) of the same pair; families_both_members = both genes "
                           "have >= 1 copy in some family (the most the families can recover)",
                           "hits = round(sensitivity x n) from spectrum.tsv, or counted per pair; ci = Wilson 95%, "
                           "which treats pairs as independent (they are not: fig6_gw_groups)"])
        figlib.write_table(
            "fig6_gw_groups",
            ["sample", "sample_label", "band", "stratum", "group_label", "n_pairs", "n_groups", "view", "hits",
             "groups_with_hits", "families_with_hits", "sensitivity", "ci_lo", "ci_hi"], grp_rows, generator=GEN,
            inputs=inputs, data_dir=data_dir,
            notes=notes + [uni_note, fam_note,
                           "group = paralogue group: a connected component of the scored pairs' graph (all bands; "
                           "built here, not an Ensembl gene tree); stratum top_group = the group with the most pairs "
                           "in that band, other_groups = the rest; same_chromosome / cross_chromosome = the two genes "
                           "share a contig or not (a symbol on two contigs counts as shared)",
                           "band 60-90 = the pooled 80-90, 70-80 and 60-70 bands; band all = every band",
                           "edge view: the truth_pairs.tsv `recovered` flag (asm20 or asm20 -k11 -w5; nucleotide)",
                           "families_with_hits = default families holding >= 1 recovered pair; ci = Wilson 95% (pairs "
                           "are not independent)"])
        figlib.write_table(
            "fig6_gw_precision",
            ["sample", "sample_label", "measure", "level", "identity_band", "k", "n", "precision", "ci_lo", "ci_hi",
             "families_with_pairs", "largest_family", "largest_family_pairs", "k_without_largest",
             "n_without_largest", "copies", "single_exon_copies"], prec_rows, generator=GEN, inputs=inputs,
            data_dir=data_dir,
            notes=notes + [fam_note,
                           "precision over JUDGEABLE pairs only: both genes have Compara data (any chromosome); pairs "
                           "with a gene Compara does not list are not judged, so precision is an upper bound",
                           "edge rows: aligned representative pairs above the tier's floor, in bands of the tier's own "
                           "nucleotide identity (matches / block length), 'all' = pooled",
                           "family rows: within-family gene pairs of the default families (score.py pairs --chrom "
                           "ALL); multi_exon = copies with >= 2 exons; largest_family = the family holding the most "
                           "judged pairs; k/n_without_largest = the judged pairs left without it"])
    if seed_rows:
        _write_seeding(seed_rows, {**inputs, **seed_inputs}, data_dir, notes)
    _write_protein_tiers(cfg, data_dir)


def _compara_rows(cfg, sid, prefix, copies, compara, w) -> dict:
    """Main panels a-c rows of one human sample (whole genome, S0; and S1 = pairs without a chr16 gene)."""
    ann = _o1.annotation_cache(cfg, "human")
    label = _o1.sample_label(cfg, sid)
    spec = _o1.parse_spectrum(prefix)
    universe = Path(f"{prefix}.truth_pairs.tsv")
    scored = _o1.score_members(cfg, copies, universe, ann["exons_gtf"], w / f"{sid}.families.score.txt",
                               chrom="ALL", compara=str(compara))
    n_uni = _o1.universe_pairs(prefix)
    if scored["universe_pairs"] != n_uni:
        raise RuntimeError(f"{sid}: universe mismatch, families scorer {scored['universe_pairs']} vs spectrum {n_uni}")
    att = _o1.compara_attribution(cfg, copies, prefix, ann["exons_gtf"], spec, scored, chrom="ALL",
                                  compara=str(compara), dev_contigs=_o1.substrate_drop("human", "S1"))
    n_cp, n_single = _o1.copy_exon_profile(copies, None)
    return _panel_rows([sid, label], spec, att, scored, n_cp, n_single, with_substrates=True)


def _panel_rows(base, spec, att, scored, n_cp, n_single, with_substrates: bool) -> dict:
    """Rows of the recall / groups / precision tables from one attribution (genome scope: `base` = [sample, label]
    and substrates S0/S1; dev scope: `base` = [species, chrom, truth] and one substrate)."""
    rec, grp, prec = [], [], []
    for band in _o1.COMPARA_BANDS:
        if band not in spec["recall"]:
            continue
        n, fr = spec["recall"][band]
        ps = [p for p in att["pairs"] if p["band"] == band]
        subs = (("S0", ps), ("S1", [p for p in ps if not p["dev"]])) if with_substrates else (("", ps),)
        for sub, sp_ in subs:
            nn = len(sp_)
            if not nn:
                continue
            views = []
            if sub in ("S0", ""):
                for tier, view in _o1.SPECTRUM_VIEWS.items():
                    if tier in fr and tier != "T1+T2+T3":
                        views.append((view, round(fr[tier] * n)))
                hits, n2 = scored["recall_universe"].get(band, (0, n))
                if n2 != n:
                    raise RuntimeError(f"band {band}: spectrum n={n} but families universe n={n2}")
                views.append(("families", hits))
            else:
                views.append(("edge_nt", sum(p["edge"] for p in sp_)))
                views.append(("families", sum(p["family"] for p in sp_)))
            views += [("families_edge", sum(1 for p in sp_ if p["family"] and p["edge"])),
                      ("families_no_edge", sum(1 for p in sp_ if p["family"] and not p["edge"])),
                      ("families_both_members", sum(1 for p in sp_ if p["in_families"]))]
            for view, k in views:
                if with_substrates:
                    rec.append([*base, "human", sub, "compara", band, nn, view, k, k / nn, *_o1.wilson(k, nn)])
                else:
                    rec.append([*base, band, nn, view, k, k / nn, *_o1.wilson(k, nn)])
    band_sets = ([(b, [b]) for b in _o1.COMPARA_BANDS if b in spec["recall"]] + [_o1.POOLED]
                 + [("all", list(_o1.COMPARA_BANDS))])
    for blabel, bands in band_sets:
        ps = [p for p in att["pairs"] if p["band"] in bands]
        if not ps:
            continue
        per_group: dict = {}
        for p in ps:
            per_group[p["group"]] = per_group.get(p["group"], 0) + 1
        top = max(per_group, key=lambda g: (per_group[g], g))
        strata = [("all", "all groups", ps),
                  ("top_group", _o1._group_label(att["members"][top]), [p for p in ps if p["group"] == top]),
                  ("other_groups", "", [p for p in ps if p["group"] != top])]
        if with_substrates:
            strata += [("same_chromosome", "both genes on one chromosome", [p for p in ps if not p["cross"]]),
                       ("cross_chromosome", "genes on different chromosomes", [p for p in ps if p["cross"]])]
        for stratum, glabel, sp_ in strata:
            n = len(sp_)
            if not n:
                continue
            groups = {p["group"] for p in sp_}
            if stratum == "other_groups":
                glabel = f"{len(groups)} other groups"
            views = [("edge_nt", [p for p in sp_ if p["edge"]]),
                     ("families_both_members", [p for p in sp_ if p["in_families"]]),
                     ("families", [p for p in sp_ if p["family"]]),
                     ("families_edge", [p for p in sp_ if p["family"] and p["edge"]]),
                     ("families_no_edge", [p for p in sp_ if p["family"] and not p["edge"]]),
                     ("edge_nt_not_families", [p for p in sp_ if p["edge"] and not p["family"]])]
            for view, hit in views:
                fams = {f for p in hit for f in p["families"]} if view.startswith("families_") or view == "families" \
                    else None
                if view == "families_both_members":
                    fams = None
                k = len(hit)
                grp.append([*base, blabel, stratum, glabel, n, len(groups), view, k, len({p["group"] for p in hit}),
                            "" if fams is None else len(fams), k / n, *_o1.wilson(k, n)])
    pooled: dict = {}
    for measure, band, n, frac in spec["precision"]:
        k = round(frac * n)
        pooled.setdefault(measure, [0, 0])
        pooled[measure][0] += k
        pooled[measure][1] += n
        prec.append([*base, measure, "edge", band, k, n, k / n, *_o1.wilson(k, n), "", "", "", "", "", "", ""])
    for measure, (k, n) in pooled.items():
        prec.append([*base, measure, "edge", "all", k, n, k / n, *_o1.wilson(k, n), "", "", "", "", "", "", ""])
    for measure in ("families_all_copies", "families_multi_exon"):
        a = att["precision"][measure]
        k, n = a["k"], a["n"]
        prec.append([*base, measure, "family", "all", k, n, k / n if n else None, *_o1.wilson(k, n), a["families"],
                     f"{a['top']} {a['top_label']}", a["top_pairs"], a["loo_k"], a["loo_n"], n_cp, n_single])
    return {"recall": rec, "groups": grp, "precision": prec}


# ---------------------------------------------------------------- supplement fig6s_seeding (genome scope)
SEED_COLS = ["sample", "sample_label", "species", "substrate", "reference", "config", "metric", "band", "k", "n",
             "value", "ci_lo", "ci_hi", "top_family", "top_family_pairs", "k_without_top", "n_without_top", "covered",
             "universe_genes", "status"]
SEED_CONFIGS = (("rustle_primary", "families_primary"), ("rustle", "families"))


def _seeding_gw(cfg, by_sp, compara, budget, w, miss) -> tuple[list, dict]:
    """Every sample x seeding configuration x substrate (S1 headline, S0): Compara (human), Liftoff copy pairs (every
    sample), protein-homology families (secondary; only with fig6s_protein_homology 1). A missing input is a status
    row, never a number."""
    import _liftoff as L
    rows, inputs = [], {}
    score_py, lib_py = Path(cfg["repo"]) / "bench" / "score.py", Path(cfg["repo"]) / "bench" / "lib.py"
    for sp, sids in by_sp.items():
        for sid in sids:
            label = _o1.sample_label(cfg, sid)
            copies = {}
            for config, stage in SEED_CONFIGS:
                try:
                    copies[config] = _o1.families_copies(cfg, sid, stage)
                    inputs[f"copies_{sid}_{config}"] = copies[config]
                except _o1.NotBuilt as e:
                    miss(e)
                    rows.append([sid, label, sp, "-", "-", config, "-", "-", "", "", "", "", "", "", "", "", "", "",
                                 "", f"{stage}: not built ({str(e).split(':')[0]})"])
            if not copies:
                continue
            # Compara (human samples)
            if sp == "human" and compara:
                try:
                    counts = _o1.primary_counts_gw(cfg, sid, budget)
                    inputs[f"primary_counts_{sid}"] = counts
                    ann = _o1.annotation_cache(cfg, "human")
                    for config, cp in copies.items():
                        for sub in ("S1", "S0"):
                            rows += _unit(w, f"seed_compara_{sid}_{config}_{sub}",
                                          _fp(cp, counts, compara, ann["exons_gtf"], score_py, lib_py),
                                          lambda sid=sid, config=config, cp=cp, sub=sub, counts=counts:
                                          _seed_compara(cfg, sid, label, config, cp, sub, counts, compara, w),
                                          budget, 150)
                except _o1.NotBuilt as e:
                    miss(e)
            # Liftoff copy pairs (every sample)
            loci, sup, why = L.support_if_ready(cfg, sid, sp)
            if loci is None:
                miss(why)
                for config in copies:
                    rows.append([sid, label, sp, "-", "liftoff", config, "-", "-", "", "", "", "", "", "", "", "", "",
                                 "", "", why])
            else:
                inputs[f"liftoff_loci_{sp}"] = L.loci_path(cfg, sp)
                inputs[f"liftoff_support_{sid}"] = L.support_path(cfg, sid)
                for config, cp in copies.items():
                    rows += _unit(w, f"seed_liftoff_{sid}_{config}",
                                  _fp(cp, L.loci_path(cfg, sp), L.support_path(cfg, sid)),
                                  lambda sid=sid, config=config, cp=cp, loci=loci, sup=sup:
                                  _seed_liftoff(sid, label, sp, config, cp, loci, sup), budget, 60)
            # protein-homology families: secondary, supplementary, only on request
            if _flag(cfg, "fig6s_protein_homology"):
                try:
                    rows += _seed_protein_homology(cfg, sid, label, sp, budget, w, inputs, miss)
                except _o1.NotBuilt as e:
                    miss(e)
    return rows, inputs


def _seed_compara(cfg, sid, label, config, copies, sub, counts, compara, w) -> list:
    drop = _o1.substrate_drop("human", sub)
    d = w / sid
    d.mkdir(exist_ok=True)
    uni, n_pairs, n_genes = _o1.compara_universe(cfg, sid, counts, drop)
    cp = _o1.filter_rows(copies, d / f"{config}.{sub}.copies.tsv", "chrom", drop) if drop else copies
    ann = _o1.annotation_cache(cfg, "human")
    scored = _o1.score_members(cfg, cp, uni, ann["exons_gtf"], d / f"seed_{config}.{sub}.compara.score.txt",
                               chrom="ALL", compara=str(compara))
    at = _o1.compara_universe_attribution(cfg, cp, uni, ann["exons_gtf"], str(compara), scored)
    out = []
    for q in (">=90", "80-90", "precision"):
        a = at[q]
        k, n = a["k"], a["n"]
        metric = "precision" if q == "precision" else "sensitivity"
        out.append([sid, label, "human", sub, "compara", config, metric, "all" if q == "precision" else q, k, n,
                    k / n if n else None, *_o1.wilson(k, n), a["top"], a["top_pairs"], a["loo_k"],
                    a["loo_n"] if q == "precision" else n, "", n_genes, "ok"])
    return out


def _seed_liftoff(sid, label, sp, config, copies, loci, sup) -> list:
    """Liftoff (record, extra copy) pairs with both loci read-supported in the sample; recovered = in one family."""
    import collections
    import _liftoff as L
    fam_loci = L.catalog_loci(copies)
    out = []
    for sub in ("S1", "S0"):
        drop = _o1.substrate_drop(sp, sub)
        keep = (lambda c: c not in drop) if drop else None
        for sc in L.SC_ROWS:
            res = L.pair_families(L.copy_pairs(loci, sc, sup, keep), fam_loci)
            hits = [int(bool(shared)) for _, _, shared in res]
            by = collections.Counter(f for _, _, shared in res for f in shared)
            top = max(by, key=lambda f: (by[f], f)) if by else ""
            loo = sum(1 for _, _, shared in res if shared - {top}) if top else sum(hits)
            lo, hi = L.boot_ci([src for src, _, _ in res], hits)
            k, n = sum(hits), len(res)
            out.append([sid, label, sp, sub, "liftoff", config, "sensitivity", f"sc>={sc:.2f}", k, n,
                        k / n if n else None, lo, hi, top, by[top] if top else 0, loo, n,
                        sum(1 for _, both, _ in res if both), "", "ok"])
    return out


def _seed_protein_homology(cfg, sid, label, sp, budget, w, inputs, miss) -> list:
    """Secondary reference (supplement only): the species' protein-homology families at >= 90% annotated-mRNA identity,
    the families-stage CLUSTERS of both configurations (score.py pairs --bands paf:), genes with >= 2 primary reads."""
    ph = _o1.protein_homology_families(cfg, sp)
    ann = _o1.annotation_cache(cfg, sp)
    ppaf = _o1.mrna_pairs_paf(cfg, sp, budget)
    counts = _o1.primary_counts_gw(cfg, sid, budget)
    ph_genes = {(c, n) for c, n, _ in _o1.read_ph(ph)}
    rows = []
    inputs.update({f"protein_homology_{sp}": ph, f"mrna_pairs_paf_{sp}": ppaf})
    for config, stage in SEED_CONFIGS:
        try:
            clusters = _o1.stage_product(cfg, sid, stage, "clusters")
        except _o1.NotBuilt as e:
            miss(e)
            continue
        for sub in ("S1", "S0"):
            drop = _o1.substrate_drop(sp, sub)
            d = w / sid
            d.mkdir(exist_ok=True)
            expr, n_expr, _n_fam = _o1.expressed_gw(cfg, sid, counts, drop, restrict=ph_genes)
            truth = _o1.filter_rows(ph, d / f"ph.{sub}.tsv", "Contig", drop) if drop else ph
            cl = _o1.filter_rows(clusters, d / f"{config}.{sub}.clusters.tsv", "chrom", drop) if drop else clusters
            if _o1.clusters_rows(cl) == 0:
                continue
            kw = dict(genes=ann["genes_gff"], truth=truth, paf=ppaf, chrom="ALL")
            r = _o1.score_referee(cfg, cl, f"{sid}_{config}_{sub}", d / f"{config}.{sub}.ph.score.txt", expr, **kw)
            at = _o1.referee_attribution(cfg, cl, r, expr, **kw)
            for q in (">=90", "precision"):
                if q not in at:
                    continue
                a = at[q]
                k, n = a["k"], a["n"]
                rows.append([sid, label, sp, sub, "protein_homology", config,
                             "precision" if q == "precision" else "sensitivity", "all" if q == "precision" else q, k,
                             n, k / n if n else None, *_o1.wilson(k, n), a["top"] or "", a["top_pairs"], a["loo_k"],
                             a["loo_n"] if q == "precision" else n, "", n_expr, "ok"])
    return rows


def _write_seeding(rows, inputs, data_dir, notes):
    bad = [r for r in rows if len(r) != len(SEED_COLS)]
    if bad:
        raise RuntimeError(f"fig6s_seeding: {len(bad)} rows do not have {len(SEED_COLS)} columns: {bad[0]}")
    figlib.write_table(
        "fig6s_seeding", SEED_COLS, rows, generator=GEN, inputs=inputs, data_dir=data_dir,
        notes=notes + [
            "supplementary: seeding loci with secondary alignments (config rustle = the default: secondary alignments "
            "scoring >= 98% of the read's genome-wide best) vs primary alignments only (config rustle_primary), each "
            "through the same families stage (the default families of each assembly; copy tables)",
            "reference compara (human samples): Ensembl Compara release 116 pairs, cross-chromosome kept, whose two "
            "genes both have a record with >= 2 reads whose PRIMARY alignment (not 0x100/0x800/0x4) has an aligned "
            "block (M/=/X) on its annotated exons, neither gene on a dropped contig (_o1.compara_universe; the same "
            "universe for both configurations); sensitivity = pairs with copies in one family (score.py pairs "
            "--universe), band = Compara protein identity; precision over judgeable within-family pairs (both genes "
            "have Compara data; an upper bound)",
            "reference liftoff (every sample): the Fig. 8 self-lift's (source record's annotated placement, extra "
            "copy) pairs, extra copy sequence_ID >= band's sc, both exon unions >= 200 bp, both loci supported by >= 2 "
            "reads of the sample (primary alignment, -F 2308, aligned block on the exon union); recovered = a copy "
            "covering >= 50% of each locus's exon union and the two covering copies share a family (Liftoff's -a); "
            "covered = pairs whose two loci are both covered (the conditional denominator of Fig. 8's F1); ci = "
            "2,000 bootstrap resamples of source records",
            "reference protein_homology: SECONDARY, supplementary only, scored only with fig6s_protein_homology 1 "
            "(the species' genome-wide protein-homology families; >= 90% annotated-mRNA identity band)",
            "substrate S1 = genome minus development contigs (human chr16, gorilla NC_073244.2; the headline), S0 = "
            "whole genome; top_family = the family holding the most recovered (sensitivity) or judged (precision) "
            "pairs, k/n_without_top = the value without it; status != ok rows: not built (the message says what is "
            "missing)",
            "ci = Wilson 95% for compara / protein_homology rows (treats pairs as independent; they are not)"])


def _write_protein_tiers(cfg, data_dir, prefix: Path | None = None, notes=()) -> bool:
    """Supplement fig6s_protein: the chr16 development spectrum WITH the translated tier (T3), tiers T1 / T1+T2 /
    T1+T2+T3 by Compara band, and each tier's precision. Written whenever that spectrum exists."""
    prefix = prefix or _o1.chr16_t3_example(cfg)
    if prefix is None or not Path(f"{prefix}.spectrum.tsv").exists():
        return False
    spec = _o1.parse_spectrum(prefix)
    if _o1.widest_tier(spec) != "T1+T2+T3":
        return False
    rows = []
    for band in _o1.COMPARA_BANDS:
        if band not in spec["recall"]:
            continue
        n, fr = spec["recall"][band]
        for tier, view in _o1.SPECTRUM_VIEWS.items():
            if tier in fr:
                k = round(fr[tier] * n)
                rows.append(["human", _o1.HUMAN_CHROM, "recall", view, band, k, n, k / n, *_o1.wilson(k, n)])
    pooled: dict = {}
    for measure, band, n, frac in spec["precision"]:
        k = round(frac * n)
        pooled.setdefault(measure, [0, 0])
        pooled[measure][0] += k
        pooled[measure][1] += n
        rows.append(["human", _o1.HUMAN_CHROM, "precision", measure, band, k, n, k / n, *_o1.wilson(k, n)])
    for measure, (k, n) in pooled.items():
        rows.append(["human", _o1.HUMAN_CHROM, "precision", measure, "all", k, n, k / n, *_o1.wilson(k, n)])
    figlib.write_table(
        "fig6s_protein_tiers", ["species", "chrom", "metric", "view", "band", "k", "n", "value", "ci_lo", "ci_hi"],
        rows, generator=GEN, inputs={"spectrum_tsv": f"{prefix}.spectrum.tsv",
                                     "score_py": Path(cfg["repo"]) / "bench" / "score.py"}, data_dir=data_dir,
        notes=list(notes) + [
            "SUPPLEMENTARY comparator, not part of Rustle's rule: the translated protein search (T3: mmseqs "
            "easy-search --search-type 2, protein identity >= 0.30, e <= 1e-5, max(query, target coverage) >= 0.50) "
            "added to the direct nucleotide alignment (T1 asm20 identity >= 0.80; T2 asm20 -k11 -w5 identity >= 0.60; "
            "coverage >= 0.50), human A119b chr16 (the development chromosome; score.py spectrum); never run "
            "genome-wide (hours, >= 13 GB per sample)",
            "recall rows: Compara pairs whose two genes both have a locus in the chr16 assembly, by Compara protein "
            "identity band; view edge_asm20 = T1, edge_nt = T1 or T2, edge_nt_protein = T1, T2 or T3",
            "precision rows: aligned locus pairs above each tier's floor, judgeable (both genes have Compara data), in "
            "bands of the tier's own identity (nucleotide for T1/T2, protein for T3); 'all' = pooled",
            "ci = Wilson 95% (treats pairs as independent; they are not)"])
    return True


# ---------------------------------------------------------------- development scope
def build_recorded(cfg: dict, data_dir: Path):
    """The recorded development runs: the supplementary tables only (T3 comparator on the recorded chr16 spectrum;
    gorilla seeding configurations against the contig's protein-homology families). No recorded run of the default
    families' copy table exists: the main development tables need `make.py data fig6 --set fig6_scope=dev`."""
    src = _o1.recorded_sources(cfg)
    _write_protein_tiers(cfg, data_dir, src["spectrum_t3"], notes=[PROVISIONAL, src["source"]])
    _write_gorilla(cfg, data_dir, src, [PROVISIONAL, src["source"]],
                   figlib.work_dir(cfg, _o1.FIG) / "recorded", force=False)
    print("[fig6] --recorded: supplementary development tables only (fig6s_protein_tiers, fig6_gorilla_*); the main "
          "development tables need `make.py data fig6 --set fig6_scope=dev`", file=sys.stderr)


def build_dev(cfg: dict, data_dir: Path, force: bool = False):
    """Regenerate the development tables (HEAVY unless cached: see _o1.ensure_sources): main (human chr16, the default
    families) and the development supplements (chr16 T3 comparator; gorilla seeding vs protein-homology families)."""
    src = _o1.ensure_sources(cfg, force=force)
    w = figlib.work_dir(cfg, _o1.FIG) / "score"
    w.mkdir(parents=True, exist_ok=True)
    notes = [src["source"]]
    ref16 = _o1.chr16_ref_gtf(cfg)

    spec = _o1.parse_spectrum(src["spectrum"])
    universe = Path(f"{src['spectrum']}.truth_pairs.tsv")
    fams = src["families"]
    scored = _o1.score_members(cfg, fams, universe, ref16, w / "chr16_families.score.txt")
    n_uni = _o1.universe_pairs(src["spectrum"])
    if scored["universe_pairs"] != n_uni:
        raise RuntimeError(f"universe mismatch: families scorer {scored['universe_pairs']} vs spectrum {n_uni} pairs")
    att = _o1.compara_attribution(cfg, fams, src["spectrum"], ref16, spec, scored)
    n_cp, n_single = _o1.copy_exon_profile(fams, _o1.HUMAN_CHROM)
    rows = _panel_rows(["human", _o1.HUMAN_CHROM, "compara"], spec, att, scored, n_cp, n_single,
                       with_substrates=False)
    inputs = {"spectrum_tsv": f"{src['spectrum']}.spectrum.tsv", "spectrum_truth_pairs": universe,
              "chr16_families_copies": fams, "mcl_families": src["families_bin"],
              "compara": cfg["compara_chr16"], "chr16_ref_gtf": ref16, "human_ref_gtf": cfg["human_ref_gtf"],
              "score_pairs_log": w / "chr16_families.score.txt",
              "score_py": Path(cfg["repo"]) / "bench" / "score.py", "lib_py": Path(cfg["repo"]) / "bench" / "lib.py"}
    uni_note = (f"pairs scored: {n_uni} Compara paralogue pairs on chr16 whose two genes both have a locus "
                f"(representative >= 200 bp) in Rustle's chr16 assembly (spectrum truth_pairs, score.py spectrum "
                f"--skip-t3; conditioned on the assembler's loci, not on read counts)")
    fam_note = ("families = the DEFAULT de novo families: the driver's families stage (mcl_families --from-gtf "
                "--min-exonic-bp 1 --min-shared-exon-frac 0.60 --emit-units, MCL 2.8) on the chr16 restriction of the "
                "genome-wide A119b assembly (loci seeded with secondary alignments within 2% of the read's best "
                "score); one copy per member locus = its representative transcript; copies mapped to genes by span "
                f"(score.py pairs): {scored['copies']} copies, {scored['families']} families, largest "
                f"{scored['largest']} genes; single-exon copies {n_single}/{n_cp}")
    figlib.write_table(
        "fig6_chr16_recall",
        ["species", "chrom", "truth", "compara_band", "n_pairs", "view", "hits", "recall", "ci_lo", "ci_hi"],
        rows["recall"], generator=GEN, inputs=inputs, data_dir=data_dir,
        notes=notes + [
            uni_note, fam_note,
            "view edge_asm20 = minimap2 asm20, identity >= 0.80, coverage >= 0.50 (query span / shorter length); "
            "edge_nt = asm20 or asm20 -k11 -w5 (identity >= 0.60, coverage >= 0.50); a DIRECT nucleotide alignment "
            "between any loci of the two genes (one representative per locus); the translated protein tier is in the "
            "supplement (fig6s_protein_tiers)",
            "view families = both genes in one default family (score.py pairs --universe); families_edge / "
            "families_no_edge split it by the direct alignment of the same pair (edge_nt): with one, or joined "
            "without one (through other loci of the family); families_both_members = both genes have >= 1 copy in "
            "some family (the most the families can recover)",
            "hits = round(sensitivity x n) from spectrum.tsv (exact: the scorer writes hits/n); column recall = "
            "sensitivity; ci = Wilson 95%, which treats pairs as independent (they are not: see fig6_chr16_groups)",
        ])
    figlib.write_table(
        "fig6_chr16_groups",
        ["species", "chrom", "truth", "band", "stratum", "group_label", "n_pairs", "n_groups", "view", "hits",
         "groups_with_hits", "families_with_hits", "recall", "ci_lo", "ci_hi"], rows["groups"],
        generator=GEN, inputs=inputs, data_dir=data_dir,
        notes=notes + [
            uni_note, fam_note,
            "group = a paralogue group: a connected component of the scored pairs' graph (all bands; built here, not "
            "an Ensembl gene tree); stratum top_group = the group with the most pairs in that band (or pooled "
            "range), other_groups = the rest",
            "band 60-90 = the pooled 80-90, 70-80 and 60-70 bands; band all = every band",
            "views as fig6_chr16_recall; edge_nt_not_families = pairs with a direct alignment that no family holds",
            "groups_with_hits = paralogue groups with >= 1 recovered pair; families_with_hits = default families "
            "holding >= 1 recovered pair; column recall = sensitivity; ci = Wilson 95% (pairs are not independent)",
        ])
    figlib.write_table(
        "fig6_chr16_precision",
        ["species", "chrom", "truth", "measure", "level", "identity_band", "k", "n", "precision", "ci_lo", "ci_hi",
         "families_with_pairs", "largest_family", "largest_family_pairs", "k_without_largest", "n_without_largest",
         "copies", "single_exon_copies"],
        rows["precision"], generator=GEN, inputs=inputs, data_dir=data_dir,
        notes=notes + [
            fam_note,
            "precision over JUDGEABLE pairs only: pairs whose two genes both have Compara data (an upper bound)",
            "edge rows: aligned representative pairs above the tier's floor, in bands of the tier's own nucleotide "
            "identity (T1 asm20, T2 asm20 -k11 -w5), 'all' = pooled; k = round(precision x n) from spectrum.tsv",
            "family rows: within-family gene pairs of the default families (score.py pairs); multi_exon = pairs of "
            "copies with n_exon >= 2; families_with_pairs = families holding >= 1 judged pair; largest_family = the "
            "one holding the most judged pairs (family id + its gene-name group); k/n_without_largest = the judged "
            "pairs left when that family is removed",
            "ci = Wilson 95% (treats pairs as independent; they are not)",
        ])
    _write_protein_tiers(cfg, data_dir, src["spectrum_t3"], notes=notes)
    _write_gorilla(cfg, data_dir, src, notes, w, force)


def _write_gorilla(cfg, data_dir, src, notes, w, force):
    """Development supplement fig6s_seeding: gorilla OR6737 NC_073244.2 (the seeding decision's contig), both seeding
    configurations through the families stage, against that contig's protein-homology families (a SECONDARY
    reference, supplementary only), by annotated-mRNA identity band. Genes scored: >= 2 primary reads on their exons."""
    w.mkdir(parents=True, exist_ok=True)
    notes = list(notes)
    ref = _o1.referee_paths(cfg)
    contig = _o1.gorilla_contig(cfg)
    counts = _o1.primary_counts(cfg, force=force)
    expr = _o1.expressed_list(cfg, "exon", force=force)
    expr_span = _o1.expressed_list(cfg, "span", force=force)
    pc = _o1.read_primary_counts(counts)
    n_exon = sum(1 for v in pc.values() if v["exon"] >= _o1.EXPRESSED_MIN_PRIMARY)
    n_span = sum(1 for v in pc.values() if v["span"] >= _o1.EXPRESSED_MIN_PRIMARY)
    span_set = {g for g, v in pc.items() if v["span"] >= _o1.EXPRESSED_MIN_PRIMARY}
    rec_path = ref["expressed_recorded"]
    if Path(rec_path).exists():
        recorded_set = {ln.split("\t")[0] for ln in open(rec_path) if not ln.startswith("Gene")}
        rec_note = (f"the span rule reproduces the recorded gene list of registers 1096-1101 exactly "
                    f"({len(recorded_set)} genes)" if recorded_set == span_set else
                    f"the span rule does NOT reproduce the recorded gene list ({len(span_set)} vs {len(recorded_set)} "
                    f"genes; {len(span_set ^ recorded_set)} differ)")
    else:
        rec_note = "recorded gene list absent; not compared"
    g_rows, c_rows = [], []
    g_inputs = {"referee": ref["referee"], "genes": ref["genes"], "mrna_paf": ref["mrna_paf"],
                "gorilla_bam": cfg["gorilla_bam"], "score_py": Path(cfg["repo"]) / "bench" / "score.py",
                "primary_counts": counts, "expressed_exon": expr,
                "expressed_span": expr_span, "expressed_recorded": rec_path}
    g_notes = notes + [
        "SUPPLEMENTARY, development: the reference is a SECONDARY one (the contig's protein-homology families; the "
        "genome-wide supplement uses Ensembl Compara and Liftoff copy pairs; protein is not part of Rustle's default "
        "family definition)",
        f"genes scored (plotted, metric recall = sensitivity): protein-homology genes with >= "
        f"{_o1.EXPRESSED_MIN_PRIMARY} reads whose PRIMARY record has an aligned block (M/=/X; not N, not D) on one of "
        f"the gene's annotated exons (_o1.primary_counts; fig. 3's n_mol_primary rule): {n_exon} of {len(pc)} genes",
        f"metric recall_span_universe: the same scorer over the recorded list's rule, >= {_o1.EXPRESSED_MIN_PRIMARY} "
        f"primary records overlapping the gene SPAN (spliced-over reads included): {n_span} of {len(pc)} genes; "
        + rec_note,
    ]
    for arm in ("rustle_primary", "rustle"):
        log = w / f"gorilla_{arm}.score.txt"
        r = _o1.score_referee(cfg, src["clusters"][arm], f"ggo_{arm}", log, expr)
        r_span = _o1.score_referee(cfg, src["clusters"][arm], f"ggo_{arm}_span", w / f"gorilla_{arm}.span.score.txt",
                                   expr_span)
        g_inputs[f"clusters_{arm}"] = src["clusters"][arm]
        if src.get("arm_gtf", {}).get(arm):
            g_inputs[f"gtf_{arm}"] = src["arm_gtf"][arm]
        g_inputs[f"score_pairs_log_{arm}"] = log
        g_inputs[f"score_pairs_log_{arm}_span"] = w / f"gorilla_{arm}.span.score.txt"
        for metric, res in (("recall", r), ("recall_span_universe", r_span)):
            for band in _o1.MRNA_BANDS:
                if band in res["recall"]:
                    k, n = res["recall"][band]
                    g_rows.append(["gorilla", contig, "protein_referee", arm, metric, band, k, n, k / n,
                                   *_o1.wilson(k, n)])
        if tuple(r_span["precision"]) != tuple(r["precision"]):
            raise RuntimeError(f"gorilla {arm}: precision depends on the genes scored ({r['precision']} vs "
                               f"{r_span['precision']})")
        k, n = r["precision"]
        g_rows.append(["gorilla", contig, "protein_referee", arm, "precision", "all", k, n,
                       k / n if n else None, *_o1.wilson(k, n)])
        comp = ", ".join(f"{f}:{c}" for f, c in r["largest_composition"].items())
        g_rows.append(["gorilla", contig, "protein_referee", arm, "largest_cluster_genes", "all",
                       r["largest"], "", r["largest"], "", ""])
        g_notes.append(f"{arm}: largest predicted family {r['largest']} genes (protein-homology families {comp})")
        at = _o1.gorilla_attribution(cfg, src["clusters"][arm], r, expr)
        for q in [b for b in _o1.MRNA_BANDS if b in at] + ["precision"]:
            a = at[q]
            top = a["top"]
            loo_src = "attribution"
            if top and q in (">=90", "precision"):
                loo_file = _o1.drop_cluster(src["clusters"][arm], top, w / f"gorilla_{arm}.without_{top}.clusters.tsv")
                r2 = _o1.score_referee(cfg, loo_file, f"ggo_{arm}-{top}", w / f"gorilla_{arm}.without_{top}.score.txt",
                                       expr)
                got = r2["recall"][q] if q != "precision" else r2["precision"]
                want = (a["loo_k"], a["n"]) if q != "precision" else (a["loo_k"], a["loo_n"])
                if tuple(got) != want:
                    raise RuntimeError(f"gorilla {arm} {q}: rescored without {top} {got} != attribution {want}")
                loo_src = "score.py pairs"
                g_inputs[f"without_top_{arm}_{'ge90' if q == '>=90' else q}"] = w / f"gorilla_{arm}.without_{top}.score.txt"
            info = at["cluster_info"].get(top) if top else None
            span = f"{contig}:{info[1] + 1}-{info[2]}" if info else ""
            compo = ",".join(f"{f}:{c}" for f, c in sorted(info[3].items())) if info else ""
            metric = "precision" if q == "precision" else "recall"
            c_rows.append(["gorilla", contig, "protein_referee", arm, metric, "all" if q == "precision" else q,
                           a["k"], a["n"], a["clusters"], at["n_clusters"], a["fams"],
                           a.get("fams_in_band", ""), top or "", info[0] if info else "", span, compo, a["top_pairs"],
                           a["loo_k"], a["loo_n"] if q == "precision" else a["n"], loo_src if top else "",
                           a.get("locus_pairs", ""), a.get("gene_pairs", "")])
    figlib.write_table(
        "fig6_gorilla_recall",
        ["species", "contig", "truth", "arm", "metric", "band", "k", "n", "value", "ci_lo", "ci_hi"], g_rows,
        generator=GEN, inputs=g_inputs, data_dir=data_dir,
        notes=g_notes + [
            "configuration (column arm) rustle = loci seeded with primaries + secondary alignments within 2% of the "
            "read's genome-wide best alignment score (driver default); rustle_primary = --no-seed-secondaries; both "
            "through the families stage (the default family rule)",
            "recall (sensitivity): same-family pairs whose two genes are both scored (above), banded by the best "
            "minimap2 (asm20 -k11 -w5) identity between the two genes' annotated mRNAs; 'none' = no alignment record",
            "precision: judgeable within-family gene pairs (both genes in a protein-homology family) against the "
            "COMPLETE families (score.py pairs --bands paf:; independent of the genes scored, checked); see "
            "fig6_gorilla_clusters for the judgeable / gene / locus pair counts",
            "column truth = the reference (protein-homology families of the contig); ci = Wilson 95%, which treats "
            "pairs as independent; they are not (fig6_gorilla_clusters)",
        ])
    figlib.write_table(
        "fig6_gorilla_clusters",
        ["species", "contig", "truth", "arm", "metric", "band", "k", "n", "clusters_with_pairs", "clusters_total",
         "referee_families_with_pairs", "referee_families_in_band", "top_cluster", "top_cluster_loci",
         "top_cluster_span", "top_cluster_referee_families", "top_cluster_pairs", "k_without_top", "n_without_top",
         "without_top_source", "within_cluster_locus_pairs", "within_cluster_gene_pairs"], c_rows,
        generator=GEN, inputs=g_inputs, data_dir=data_dir,
        notes=notes + [
            "SUPPLEMENTARY, development (secondary reference: the contig's protein-homology families): which predicted "
            "families and protein-homology families (columns *referee_families*) the recovered or judged pairs come "
            "from, recomputed with score.py's own loaders and checked band by band against the scorer's totals",
            g_notes[len(notes) + 1],
            "recall rows: k/n = recovered / reference pairs in the band; clusters_with_pairs = predicted families "
            "holding >= 1 recovered pair; referee_families_with_pairs / _in_band = protein-homology families with >= 1 "
            "recovered / >= 1 reference pair in the band",
            "precision rows: k/n = same-family / judgeable within-family gene pairs; within_cluster_locus_pairs = sum "
            "C(loci, 2) over predicted families, within_cluster_gene_pairs = distinct gene pairs after mapping loci to "
            "genes (judgeable = both genes in a protein-homology family)",
            "top_cluster = the predicted family holding the most recovered (recall) or judged (precision) pairs, ties "
            "by id; top_cluster_span is 1-based closed; k/n_without_top = the value with that family removed from the "
            "clusters file: RESCORED by score.py pairs for >=90 and precision (without_top_source), from the "
            "attribution for the other bands (the two agree wherever both are computed)",
        ])


# ================================================================ plot helpers
def _f(x):
    return float(x) if x not in ("", None) else None


def _i(x):
    return int(float(x)) if x not in ("", None) else None


def _count(ax, x, y, text, dy=0.015, **kw):
    kw.setdefault("fontsize", 5.5)
    kw.setdefault("color", figlib.INK_2)
    ax.text(x, y + dy, text, ha="center", va="bottom", linespacing=1.05, **kw)


def _inbar(ax, x, y0, y1, text, color="white"):
    """A rotated label inside a bar segment (only where it fits)."""
    ax.text(x, (y0 + y1) / 2, text, rotation=90, ha="center", va="center", fontsize=5.5, color=color)


def _pretty_group(label: str) -> str:
    """'NPIP* (14 genes)' -> '14 NPIP genes'; other labels unchanged."""
    import re
    m = re.match(r"^(\w+)\* \((\d+) genes\)$", label)
    return f"{m.group(2)} {m.group(1)} genes" if m else label


def _hbar_rows(ax, rows, xmax=1.0):
    """Horizontal bars, each with its label above it, k/n to its right and an optional note below it.

    rows: ("header", text) | None (a divider) | (label, [(value, colour, [in-bar texts, longest first])], "k/n",
    note)."""
    y = 0.0
    for row in rows:
        if row is None:
            ax.axhline(y + 0.30, color=figlib.INK_3, linewidth=0.5, linestyle=(0, (2, 2)))
            y -= 0.35
            continue
        if row[0] == "header":
            ax.text(0, y, row[1], ha="left", va="center", fontsize=6.2, color=figlib.INK, fontweight="bold")
            y -= 0.85
            continue
        label, segs, kn, note = row
        left = 0.0
        for v, c, txts in segs:
            ax.barh(y, v, left=left, height=0.42, color=c, edgecolor=figlib.SURFACE, linewidth=0.6)
            # the longest in-bar text that fits the segment (~0.016 of the axis per character at 5.5 pt)
            fit = [t for t in txts if 0.016 * xmax * (len(t) + 1.5) <= v]
            if fit:
                ax.text(left + 0.012 * xmax, y, fit[0], ha="left", va="center", fontsize=5.5, color="white")
            left += v
        ax.text(left + 0.015 * xmax, y, kn, ha="left", va="center", fontsize=5.8, color=figlib.INK)
        ax.text(0, y + 0.25, label, ha="left", va="bottom", fontsize=6.0, color=figlib.INK)
        if note:
            ax.text(0, y - 0.27, note, ha="left", va="top", fontsize=5.5, color=figlib.INK_2)
            y -= 0.40
        y -= 1
    ax.set_ylim(y + 0.55, 0.40)
    ax.set_yticks([])
    ax.spines["left"].set_visible(False)
    ax.set_xlim(0, xmax)
    ax.grid(axis="y", visible=False)
    ax.grid(axis="x")


def _band_bars(ax, by, bands, labels=True, fam_counts=None):
    """Panel a's bars for one sample / substrate: `by` = {(band, view): row}; left bar = direct nucleotide alignment
    (grey), right bar = the default families (blue: the two genes' loci also align directly; navy: joined without a
    direct alignment) with the dashed ceiling (both genes are family members)."""
    xs = list(range(len(bands)))
    w = 0.36
    for x, b in zip(xs, bands):
        n = int(by[(b, "families")]["n_pairs"])
        hd, hv, hf = (_i(by[(b, v)]["hits"]) for v in ("families_edge", "families_no_edge", "families"))
        if hd + hv != hf:
            raise RuntimeError(f"band {b}: with + without a direct alignment != families")
        nt = _i(by[(b, "edge_nt")]["hits"]) / n
        d, v, fam = hd / n, hv / n, hf / n
        xe, xf = x - w / 2 - 0.01, x + w / 2 + 0.01
        ax.bar(xe, nt, width=w, color=EDGE_NT, edgecolor=figlib.SURFACE, linewidth=0.6)
        ax.bar(xf, d, width=w, color=FAMILY, edgecolor=figlib.SURFACE, linewidth=0.6)
        if v > 0:
            ax.bar(xf, v, bottom=d, width=w, color=FAMILY_VIA, edgecolor=figlib.SURFACE, linewidth=0.6)
        reach_row = by.get((b, "families_both_members"))
        reach = _i(reach_row["hits"]) if reach_row else None
        if reach is not None:
            ax.bar(xf, reach / n, width=w + 0.03, fill=False, edgecolor=figlib.INK_2, linewidth=0.6,
                   linestyle=(0, (2, 1.5)), zorder=3)
        _count(ax, xe, nt, by[(b, "edge_nt")]["hits"])
        kk = f"{hf}/{reach}" if reach is not None else str(hf)
        nfam = (fam_counts or {}).get(b, "")
        _count(ax, xf, max(fam, (reach or 0) / n), f"{kk}\n{nfam} fam." if nfam not in ("", "0") else kk)
        if labels and b == ">=90":
            if nt > 0.3:
                _inbar(ax, xe, 0, nt, "direct")
            if d > 0.3:
                _inbar(ax, xf, 0, d, "same family")
        if labels and v > 0.25 and b in ("60-70", "70-80", "80-90"):
            _inbar(ax, xf, d, d + v, "no direct\nalignment")
            labels = False
    ax.set_xticks(xs)
    ax.tick_params(axis="x", length=0, pad=3)
    ax.set_xlim(-0.6, len(bands) - 0.4)
    ax.set_ylim(0, 1.1)
    ax.set_yticks([0, 0.2, 0.4, 0.6, 0.8, 1.0])


def _dot_table(ax, rows, xlabel, right_head=""):
    """A compact dot table: rows ("header", text) | (label, [(value, marker kwargs)], right text). Top to bottom."""
    y = 0.0
    for row in rows:
        if row[0] == "header":
            ax.text(-0.02, y, row[1], transform=ax.get_yaxis_transform(), ha="right", va="center", fontsize=6.2,
                    fontweight="bold", color=figlib.INK)
            if right_head:
                ax.text(1.02, y, right_head, transform=ax.get_yaxis_transform(), ha="left", va="center",
                        fontsize=5.3, color=figlib.INK_3)
            y += 1.0
            continue
        label, marks, right = row
        ax.text(-0.02, y, label, transform=ax.get_yaxis_transform(), ha="right", va="center", fontsize=5.8,
                color=figlib.INK)
        vals = [v for v, _ in marks if v is not None]
        if len(vals) >= 2:
            ax.plot([min(vals), max(vals)], [y, y], color=figlib.INK_3, linewidth=0.6, zorder=1)
        for v, kw in marks:
            if v is not None:
                ax.plot([v], [y], linestyle="none", zorder=3, **kw)
        ax.text(1.02, y, right, transform=ax.get_yaxis_transform(), ha="left", va="center", fontsize=5.5,
                color=figlib.INK_2)
        y += 1.0
    ax.set_ylim(y - 0.4, -0.6)
    ax.set_xlim(-0.03, 1.03)
    ax.set_xticks([0, 0.5, 1.0])
    ax.set_xticklabels(["0", "0.5", "1"])
    ax.set_yticks([])
    ax.spines["left"].set_visible(False)
    ax.grid(axis="y", visible=False)
    ax.grid(axis="x", color=figlib.GRID, linewidth=0.5)
    ax.set_xlabel(xlabel, labelpad=2, fontsize=6.2)


def _disp(label: str) -> str:
    """A sample label for a panel title (the unrecorded tissue is left to the caption)."""
    return label.replace(" (tissue not recorded)", "")


MK_EDGE = dict(marker="s", markersize=3.6, markerfacecolor=EDGE_NT, markeredgecolor=EDGE_NT)
MK_FAM = dict(marker="o", markersize=3.8, markerfacecolor=FAMILY, markeredgecolor=FAMILY)
MK_FAM_LOO = dict(marker="o", markersize=3.8, markerfacecolor=figlib.SURFACE, markeredgecolor=FAMILY,
                  markeredgewidth=0.8)
MK_PRIM = dict(marker="o", markersize=4.0, markerfacecolor=figlib.SURFACE, markeredgecolor=FAMILY, markeredgewidth=0.9)
MK_DEF = dict(marker="o", markersize=4.0, markerfacecolor=FAMILY, markeredgecolor=FAMILY)
FIG_NOTE = ("Rustle's default de novo families: reads → loci (seeded with secondary alignments within 2% of the best "
            "score) → one representative per locus → families (≥ 70% identity, ≥ 60% of the smaller locus's exonic "
            "sequence shared; Markov clustering). No protein.")


# ================================================================ plot (genome-wide main)
def _plot_genome(data_dir: Path, out_dir: Path):
    import matplotlib.pyplot as plt
    from matplotlib.lines import Line2D

    rec = figlib.read_table("fig6_gw_recall", data_dir)
    grp = figlib.read_table("fig6_gw_groups", data_dir)
    prec = figlib.read_table("fig6_gw_precision", data_dir)
    hs = []
    for r in rec:
        if r["substrate"] == "S0" and r["sample"] not in hs:
            hs.append(r["sample"])
    hlabel = {r["sample"]: _disp(r["sample_label"]) for r in rec}

    W, H = 7.2, 5.45
    fig = plt.figure(figsize=(W, H))

    def ax_at(left, top, w, h):
        return fig.add_axes([left / W, 1 - (top + h) / H, w / W, h / H])

    def title(x_in, top_in, letter, text):
        fig.text(x_in / W, 1 - top_in / H, letter, fontsize=9, fontweight="bold", ha="left", va="baseline")
        fig.text((x_in + 0.17) / W, 1 - top_in / H, text, fontsize=7.0, ha="left", va="baseline")

    # ---- a: per human sample, whole genome (pairs across chromosomes included)
    na = max(1, len(hs))
    wa = (W - 0.62 - 0.08 - 0.12 * (na - 1)) / na
    title(0.02, 0.16, "a", "Ensembl Compara paralogue pairs, both genes with a locus in the sample's genome-wide "
                          "assembly (pairs across chromosomes included)")
    for i in range(na):
        ax = ax_at(0.62 + i * (wa + 0.12), 0.34, wa, 1.62)
        sid = hs[i] if hs else None
        if sid is None:
            ax.text(0.5, 0.5, "not built yet", ha="center", va="center", fontsize=6.5, color=figlib.INK_3,
                    style="italic", transform=ax.transAxes)
            continue
        by = {(r["compara_band"], r["view"]): r for r in rec if r["sample"] == sid and r["substrate"] == "S0"}
        bands = [b for b in _o1.COMPARA_BANDS if (b, "families") in by]
        _band_bars(ax, by, bands, labels=i == 0)
        ax.set_xticklabels([f"{BAND_LABEL[b]}\nn={int(by[(b, 'families')]['n_pairs']):,}" for b in bands],
                           linespacing=1.15, fontsize=5.8)
        n_all = sum(int(by[(b, "families")]["n_pairs"]) for b in bands)
        ax.set_title(f"{hlabel[sid]} · {n_all:,} pairs", loc="left", fontsize=6.6, pad=3)
        if i == 0:
            ax.set_ylabel("Sensitivity (fraction of\nCompara paralogue pairs)", fontsize=6.4)
        else:
            ax.tick_params(axis="y", labelleft=False)
        ax.set_xlabel("Compara protein identity of the pair (%)", labelpad=2, fontsize=6.2)
    fig.text(0.62 / W, 1 - 2.42 / H,
             "Grey: a direct alignment between the two genes' loci (minimap2, ≥ 60% nucleotide identity over ≥ 50% of "
             "the shorter). Blue: both genes in one of Rustle's default families and their loci\nalign directly; "
             "navy: in one family without a direct alignment (joined through other loci). Dashed outline: both genes "
             "are family members (the most the families can recover); k/m over the family bar.",
             ha="left", va="top", fontsize=5.5, color=figlib.INK_2, linespacing=1.2)

    # ---- b: largest paralogue group vs the rest (60-90%), across vs within chromosomes (>= 90%)
    brows = []
    for sid in hs:
        g = {(r["band"], r["stratum"], r["view"]): r for r in grp if r["sample"] == sid}
        brows.append(("header", hlabel[sid]))
        for band, stratum, lab in (("60-90", "top_group", None), ("60-90", "other_groups", None),
                                   (">=90", "cross_chromosome", "≥ 90%, genes on two chromosomes"),
                                   (">=90", "same_chromosome", "≥ 90%, genes on one chromosome")):
            e, f = g.get((band, stratum, "edge_nt")), g.get((band, stratum, "families"))
            if not e or not f:
                continue
            n = int(f["n_pairs"])
            if lab is None:
                gl = f["group_label"]
                lab = (f"60–90%, largest group ({gl.split(' (')[0].rstrip('*')})" if stratum == "top_group"
                       else f"60–90%, {gl.replace('groups', 'other groups').replace('other other', 'other')}")
            brows.append((f"{lab}, n={n:,}", [(_f(e["sensitivity"]), MK_EDGE), (_f(f["sensitivity"]), MK_FAM)],
                          f"{e['hits']} · {f['hits']}"))
    hb = max(1.0, 0.125 * len(brows) + 0.2)
    ax_b = ax_at(1.72, 3.02, 1.0, hb)
    _dot_table(ax_b, brows, "Sensitivity", right_head="alignment · family")
    title(0.02, 2.90, "b", "By paralogue group and chromosome")
    ax_b.legend(handles=[Line2D([], [], linestyle="none", label="direct alignment", **MK_EDGE),
                         Line2D([], [], linestyle="none", label="same default family", **MK_FAM)],
                loc="upper left", bbox_to_anchor=(-1.6, -0.12 * 2.2 / hb), ncol=2, fontsize=5.5, frameon=False,
                handletextpad=0.3, columnspacing=1.0, borderaxespad=0)

    # ---- c: precision against Compara over judgeable pairs
    crows = []
    for sid in hs:
        p = {(r["measure"], r["identity_band"]): r for r in prec if r["sample"] == sid}
        crows.append(("header", hlabel[sid]))
        e = p.get(("edge_nt_sensitive", "all"))
        if e:
            crows.append(("Direct alignment (asm20 -k11 -w5)", [(_f(e["precision"]), MK_EDGE)], f"{e['k']}/{e['n']}"))
        for m, lab in (("families_all_copies", "Same family, all copies"),
                       ("families_multi_exon", "Same family, multi-exon copies")):
            r = p.get((m, "all"))
            if not r:
                continue
            loo = (_i(r["k_without_largest"]) / _i(r["n_without_largest"])
                   if r["n_without_largest"] not in ("", "0", None) else None)
            crows.append((lab, [(_f(r["precision"]), MK_FAM), (loo, MK_FAM_LOO)],
                          f"{r['k']}/{r['n']} ({r['k_without_largest']}/{r['n_without_largest']})"))
    ax_c = ax_at(5.05, 3.02, 1.0, max(1.0, 0.125 * len(crows) + 0.2))
    _dot_table(ax_c, crows, "Precision (judgeable pairs; upper bound)", right_head="k/n (w/o largest)")
    title(3.55, 2.90, "c", "Precision against Compara")
    fig.text(0.02 / W, 0.012, FIG_NOTE, fontsize=5.3, color=figlib.INK_3, ha="left", va="bottom")
    figlib.stamp_provisional(fig, GW_TABLES, data_dir)
    return figlib.save(fig, "fig6_family_spectrum", out_dir)


# ================================================================ plot (development tables)
def _panel_a(ax, rec, groups):
    by = {(r["compara_band"], r["view"]): r for r in rec}
    grp = {(r["band"], r["view"]): r for r in groups if r["stratum"] == "all"}
    bands = [b for b in _o1.COMPARA_BANDS if (b, "families") in by]
    fam_counts = {b: grp.get((b, "families"), {}).get("families_with_hits", "") for b in bands}
    _band_bars(ax, by, bands, labels=True, fam_counts=fam_counts)
    ax.set_xticklabels([f"{BAND_LABEL[b]}\nn={int(by[(b, 'families')]['n_pairs'])}\n"
                        f"{grp[(b, 'families')]['n_groups']} groups" for b in bands], linespacing=1.15)
    ax.set_ylim(0, 1.08)
    ax.set_ylabel("Sensitivity (fraction of\nCompara paralogue pairs)")
    ax.set_xlabel("Compara protein identity of the pair (%)", labelpad=3)
    ax.text(0.5, -0.335, "n = pairs in the band · groups = paralogue groups (connected components of the pairs) · "
            "dashed outline = both genes are family members (the most the families can recover)\nk/m over the family "
            "bar = pairs in one family / pairs whose genes are both family members · fam. = default families holding "
            "the recovered pairs", transform=ax.transAxes, ha="center", va="top", fontsize=5.8,
            color=figlib.INK_2, linespacing=1.2)


def _panel_b(ax, groups):
    pooled = {(r["stratum"], r["view"]): r for r in groups if r["band"] == _o1.POOLED[0]}
    rows = []
    for i, stratum in enumerate(("top_group", "other_groups")):
        e = pooled[(stratum, "edge_nt")]
        f = pooled[(stratum, "families")]
        dd, vv = _i(pooled[(stratum, "families_edge")]["hits"]), _i(pooled[(stratum, "families_no_edge")]["hits"])
        n = int(f["n_pairs"])
        if i:
            rows.append(None)
        head = (f"{_pretty_group(f['group_label'])}: one paralogue group, {n} pairs" if stratum == "top_group"
                else f"{f['group_label'].replace('groups', 'paralogue groups')}: {n} pairs")
        rows.append(("header", head))
        rows.append(("Direct alignment between the genes' loci", [(_i(e["hits"]) / n, EDGE_NT, [])],
                     f"{e['hits']}/{n}", ""))
        segs = [(dd / n, FAMILY, [f"{dd} also aligned directly", f"{dd} aligned", f"{dd}"] if dd else [])]
        if vv:
            segs.append((vv / n, FAMILY_VIA, [f"{vv} without a direct alignment", f"{vv} without", f"{vv}"]))
        rows.append(("Same default family", segs, f"{f['hits']}/{n}", ""))
    _hbar_rows(ax, rows)
    ax.set_xticks([0, 0.25, 0.5, 0.75, 1.0])
    ax.set_xlabel("Sensitivity (fraction of the group's pairs)")


def _panel_c(ax, prec):
    pooled = {r["measure"]: r for r in prec if r["identity_band"] == "all"}
    spec = [("edge_nt_sensitive", "Direct alignment, nucleotide (asm20 -k11 -w5)", EDGE_NT),
            None,
            ("families_all_copies", "Same default family, all copies", FAMILY),
            ("families_multi_exon", "Same default family, multi-exon copies only", FAMILY)]
    rows = []
    for s in spec:
        if s is None:
            rows.append(None)
            continue
        key, label, color = s
        r = pooled.get(key)
        if r is None:
            continue
        note = ""
        if r.get("families_with_pairs"):
            name = r["largest_family"].split(" ", 1)[-1].split(" (")[0].rstrip("*")
            note = (f"from {r['families_with_pairs']} families; {r['largest_family_pairs']} of the {r['n']} pairs in "
                    f"one {name} family; without it {r['k_without_largest']}/{r['n_without_largest']}")
        rows.append((label, [(_f(r["precision"]), color, [])], f"{r['k']}/{r['n']}", note))
    _hbar_rows(ax, rows, xmax=1.13)
    ax.set_xticks([0, 0.25, 0.5, 0.75, 1.0])
    ax.set_xlabel("Precision against Compara (fraction of judgeable pairs)")


def _plot_dev(data_dir: Path, out_dir: Path):
    import matplotlib.pyplot as plt

    rec = figlib.read_table("fig6_chr16_recall", data_dir)
    grp = figlib.read_table("fig6_chr16_groups", data_dir)
    prec = figlib.read_table("fig6_chr16_precision", data_dir)
    n_uni = sum(int(r["n_pairs"]) for r in rec if r["view"] == "families")

    fig = plt.figure(figsize=(7.0, 5.6))
    # rows: a | gap (a's x-axis footnote) | b + c
    outer = fig.add_gridspec(3, 1, height_ratios=[1.0, 0.56, 0.95], hspace=0.0, left=0.085, right=0.985, top=0.955,
                             bottom=0.075)
    ax_a = fig.add_subplot(outer[0])
    mid = outer[2].subgridspec(1, 2, width_ratios=[0.45, 0.55], wspace=0.10)
    ax_b = fig.add_subplot(mid[0])
    ax_c = fig.add_subplot(mid[1])
    _panel_a(ax_a, rec, grp)
    ax_a.set_title(f"Human A119b, chr16 (development chromosome) · {n_uni} Compara paralogue pairs, both genes with a "
                   f"locus in Rustle's chr16 assembly", loc="left", fontsize=7.5, pad=6)
    _panel_b(ax_b, grp)
    ax_b.set_title("Human chr16 · 60–90%, by paralogue group", loc="left", fontsize=7.5, pad=6)
    _panel_c(ax_c, prec)
    ax_c.set_title("Human chr16 · precision over judgeable pairs", loc="left", fontsize=7.5, pad=6)
    figlib.panel_label(ax_a, "a", x=-0.075, y=1.03)
    figlib.panel_label(ax_b, "b", x=-0.05, y=1.03)
    figlib.panel_label(ax_c, "c", x=-0.04, y=1.03)
    fig.text(0.005, 0.008, FIG_NOTE, fontsize=5.3, color=figlib.INK_3, ha="left", va="bottom")
    fig.text(0.995, 0.008, "Development tables: the genome-wide version replaces them", fontsize=5.3,
             color=figlib.INK_3, ha="right", va="bottom")
    figlib.stamp_provisional(fig, DEV_TABLES, data_dir)
    return figlib.save(fig, "fig6_family_spectrum", out_dir)


# ================================================================ supplementary: fig6s_protein (T3 comparator)
def plot_protein(data_dir: Path, out_dir: Path) -> list:
    import matplotlib.pyplot as plt
    if not _have(data_dir, "fig6s_protein_tiers"):
        return []
    rows = figlib.read_table("fig6s_protein_tiers", data_dir)
    rec = {(r["band"], r["view"]): r for r in rows if r["metric"] == "recall"}
    bands = [b for b in _o1.COMPARA_BANDS if (b, "edge_nt") in rec]
    fig = plt.figure(figsize=(figlib.WIDTH_DOUBLE, 2.75))
    gs = fig.add_gridspec(1, 2, width_ratios=[0.64, 0.36], left=0.085, right=0.985, top=0.80, bottom=0.25,
                          wspace=0.42)
    ax = fig.add_subplot(gs[0])
    w = 0.6
    for x, b in enumerate(bands):
        n = int(rec[(b, "edge_nt")]["n"])
        nt = int(rec[(b, "edge_nt")]["k"]) / n
        pr = int(rec[(b, "edge_nt_protein")]["k"]) / n if (b, "edge_nt_protein") in rec else nt
        ax.bar(x, nt, width=w, color=EDGE_NT, edgecolor=figlib.SURFACE, linewidth=0.6)
        if pr > nt:
            ax.bar(x, pr - nt, bottom=nt, width=w, color=EDGE_PROT, edgecolor=figlib.SURFACE, linewidth=0.6)
        _count(ax, x, pr, f"{rec[(b, 'edge_nt')]['k']}→{rec.get((b, 'edge_nt_protein'), rec[(b, 'edge_nt')])['k']}")
    ax.set_xticks(range(len(bands)))
    ax.set_xticklabels([f"{BAND_LABEL[b]}\nn={int(rec[(b, 'edge_nt')]['n'])}" for b in bands], linespacing=1.15)
    ax.tick_params(axis="x", length=0, pad=3)
    ax.set_ylim(0, 1.05)
    ax.set_ylabel("Sensitivity (fraction of\nCompara paralogue pairs)")
    ax.set_xlabel("Compara protein identity of the pair (%)", labelpad=3)
    from matplotlib.patches import Patch
    ax.legend(handles=[Patch(facecolor=EDGE_NT, label="direct nucleotide alignment (T1 or T2)"),
                       Patch(facecolor=EDGE_PROT, label="+ translated protein search (T3)")],
              loc="upper right", fontsize=5.8, frameon=False)
    ax.set_title("Pairs a direct alignment recovers, with and without the translated search", loc="left",
                 fontsize=7.0, pad=5)
    axp = fig.add_subplot(gs[1])
    pooled = {r["view"]: r for r in rows if r["metric"] == "precision" and r["band"] == "all"}
    spec = [("edge_asm20", "T1 asm20 (≥ 80% nt)", EDGE_NT), ("edge_nt_sensitive", "T2 asm20 -k11 -w5 (≥ 60% nt)",
                                                            EDGE_NT),
            ("edge_protein", "T3 translated (≥ 30% protein)", EDGE_PROT)]
    ys = []
    for i, (key, lab, col) in enumerate(spec):
        r = pooled.get(key)
        if not r:
            continue
        axp.barh(i, _f(r["value"]), height=0.55, color=col, edgecolor=figlib.SURFACE)
        axp.text(_f(r["value"]) + 0.02, i, f"{r['k']}/{r['n']}", va="center", fontsize=5.8, color=figlib.INK)
        ys.append((i, lab))
    axp.set_yticks([y for y, _ in ys])
    axp.set_yticklabels([lab for _, lab in ys], fontsize=5.8)
    axp.set_ylim(len(spec) - 0.4, -0.6)
    axp.set_xlim(0, 1.18)
    axp.set_xticks([0, 0.5, 1.0])
    axp.set_xlabel("Precision against Compara (judgeable pairs)")
    axp.grid(axis="y", visible=False)
    axp.set_title("Precision by tier", loc="left", fontsize=7.0, pad=5)
    figlib.panel_label(ax, "a", x=-0.12, y=1.06)
    figlib.panel_label(axp, "b", x=-0.62, y=1.06)
    fig.text(0.01, 0.975, "Supplementary · the translated protein search is a comparator, not part of Rustle's family "
             "rule (human A119b chr16, the development chromosome; not run genome-wide)", fontsize=6.6,
             fontweight="bold", ha="left", va="top")
    figlib.stamp_provisional(fig, ["fig6s_protein_tiers"], data_dir)
    paths = figlib.save(fig, "fig6s_protein", out_dir)
    plt.close(fig)
    return paths


# ================================================================ supplementary: fig6s_seeding
TABLE_ROWS = [("predicted families", "clusters_with_pairs", None),
              ("protein-homology families", "referee_families_with_pairs", None),
              ("in the largest family", "top_cluster_pairs", None),
              ("w/o the largest family", "k_without_top", "n_without_top")]


def _panel_seed_dev(ax_r, ax_p, gor, clus):
    import matplotlib.transforms as mtransforms

    arms = ["rustle_primary", "rustle"]
    rec = {(r["arm"], r["band"]): r for r in gor if r["metric"] == "recall"}
    prec = {r["arm"]: r for r in gor if r["metric"] == "precision"}
    cl = {(r["arm"], r["band"] if r["metric"] == "recall" else "precision"): r for r in clus}
    bands = [b for b in _o1.MRNA_BANDS if all((a, b) in rec for a in arms)]
    for b in bands:
        if any(int(rec[(a, b)]["n"]) != int(rec[(arms[0], b)]["n"]) for a in arms):
            raise RuntimeError(f"gorilla band {b}: the configurations disagree on the reference-pair count")
    ns = [int(rec[(arms[0], b)]["n"]) for b in bands]
    w = 0.36
    off = {"rustle_primary": -w / 2 - 0.01, "rustle": w / 2 + 0.01}
    for a in arms:
        for x, b in enumerate(bands):
            v = _f(rec[(a, b)]["value"])
            ax_r.bar(x + off[a], v, width=w, **figlib.tool_bar_kwargs(a))
            _count(ax_r, x + off[a], v, rec[(a, b)]["k"], dy=0.008)
    b0 = bands.index(">=90") if ">=90" in bands else None
    if b0 is not None:
        for a in arms:
            c = cl[(a, ">=90")]
            yl = int(c["k_without_top"]) / int(c["n_without_top"])
            ax_r.plot([b0 + off[a] - w / 2, b0 + off[a] + w / 2], [yl, yl], color=figlib.INK, lw=1.1,
                      solid_capstyle="butt", zorder=4)
        cp, cg = cl[(arms[0], ">=90")], cl[(arms[1], ">=90")]
        kp, kg = cp["k_without_top"], cg["k_without_top"]

        def _fams(c):
            return {x.split(":")[0] for x in c["top_cluster_referee_families"].split(",") if x}

        one_array = (cp["top_cluster_span"] and cg["top_cluster_span"] and len(_fams(cp)) == 1
                     and _fams(cp) == _fams(cg))
        ax_r.annotate(f"without each configuration's largest predicted family\n(one tandem array, "
                      f"{cp['top_cluster_loci']} vs {cg['top_cluster_loci']} loci): {kp} and {kg}"
                      if one_array else f"without each configuration's largest predicted family: {kp} and {kg}",
                      xy=(b0 + off["rustle"] + w / 2, yl), xytext=(b0 + 0.62, 0.17), fontsize=5.5,
                      color=figlib.INK, ha="left", va="bottom",
                      arrowprops=dict(arrowstyle="-", color=figlib.INK, lw=0.5, shrinkA=0, shrinkB=0))
        vp, vg = _f(rec[(arms[0], ">=90")]["value"]), _f(rec[(arms[1], ">=90")]["value"])
        ax_r.annotate("Rustle, primary alignments only", xy=(b0 + off["rustle_primary"], vp + 0.05),
                      xytext=(b0 + off["rustle_primary"] - 0.08, 0.545), fontsize=6, color=figlib.INK, ha="left",
                      va="bottom", arrowprops=dict(arrowstyle="-", color=figlib.INK_3, lw=0.5, shrinkA=0,
                                                   shrinkB=0, relpos=(0.05, 0)))
        ax_r.annotate("Rustle (default: loci also built from secondary alignments within 2% of the best score)",
                      xy=(b0 + off["rustle"], vg + 0.05), xytext=(b0 + off["rustle"] - 0.08, 0.455), fontsize=6,
                      color=figlib.INK, ha="left", va="bottom",
                      arrowprops=dict(arrowstyle="-", color=figlib.INK_3, lw=0.5, shrinkA=0, shrinkB=0,
                                      relpos=(0.08, 0)))
    ax_r.set_xticks(range(len(bands)))
    ax_r.set_xticklabels([f"{BAND_LABEL[b]}\nn={n:,}" for b, n in zip(bands, ns)], linespacing=1.15)
    ax_r.tick_params(axis="x", length=0, pad=3)
    ax_r.set_xlim(-0.6, len(bands) - 0.4)
    ax_r.set_ylim(0, 0.62)
    ax_r.set_ylabel("Sensitivity (fraction of\nsame-family pairs)")
    ax_r.set_xlabel("Nucleotide identity of the two genes' annotated mRNAs (%)", labelpad=3)

    xs_p = {"rustle_primary": 0.0, "rustle": 1.0}
    for a in arms:
        r = prec.get(a)
        if r is None:
            continue
        v = _f(r["value"]) or 0.0
        ax_p.bar(xs_p[a], v, width=0.62, **figlib.tool_bar_kwargs(a))
        ax_p.text(xs_p[a], v - 0.03, f"{int(r['k'])}/{int(r['n'])}", ha="center", va="top", rotation=90,
                  fontsize=5.8, color=figlib.INK)
    ax_p.set_xticks([xs_p[a] for a in arms])
    ax_p.set_xticklabels(["primary\nonly", "+ second-\nary ≥ 98%"], linespacing=1.15)
    ax_p.tick_params(axis="x", length=0, pad=3)
    ax_p.set_xlim(-0.6, 1.6)
    ax_p.set_ylim(0, 1.05)
    ax_p.set_yticks([0, 0.5, 1.0])
    ax_p.set_ylabel("Pair precision (fraction of\njudgeable predicted pairs)")

    y0, dy = -0.38, 0.075
    tr_r = mtransforms.blended_transform_factory(ax_r.transData, ax_r.transAxes)
    tr_p = mtransforms.blended_transform_factory(ax_p.transData, ax_p.transAxes)
    ax_r.text(-0.02, y0 + dy, "pairs come from", transform=ax_r.transAxes, ha="right", va="center", fontsize=5.5,
              color=figlib.INK_2, style="italic")
    for i, (label, kcol, ncol) in enumerate(TABLE_ROWS):
        y = y0 - i * dy
        ax_r.text(-0.02, y, label, transform=ax_r.transAxes, ha="right", va="center", fontsize=5.5,
                  color=figlib.INK_2)
        for a in arms:
            for x, b in enumerate(bands):
                c = cl.get((a, b))
                if c is not None:
                    ax_r.text(x + off[a], y, c[kcol] if c[kcol] != "" else "–", transform=tr_r, ha="center",
                              va="center", fontsize=5.5, color=figlib.INK)
            c = cl.get((a, "precision"))
            if c is not None:
                t = f"{c[kcol]}/{c[ncol]}" if ncol else (c[kcol] if c[kcol] != "" else "–")
                ax_p.text(xs_p[a], y, t, transform=tr_p, ha="center", va="center", fontsize=5.5, color=figlib.INK)


def _plot_seeding_dev(data_dir: Path, out_dir: Path) -> list:
    import matplotlib.pyplot as plt
    gor = figlib.read_table("fig6_gorilla_recall", data_dir)
    clus = figlib.read_table("fig6_gorilla_clusters", data_dir)
    contig = gor[0]["contig"] if gor else _o1.GORILLA_CONTIG_DEFAULT
    fig = plt.figure(figsize=(7.0, 3.3))
    low = fig.add_gridspec(1, 3, width_ratios=[0.075, 0.655, 0.27], wspace=0.0, left=0.085, right=0.985, top=0.80,
                           bottom=0.36)
    ax_d = fig.add_subplot(low[1])
    ax_dp = fig.add_subplot(low[2])
    pos = ax_dp.get_position()
    ax_dp.set_position([pos.x0 + 0.075, pos.y0, pos.width - 0.075, pos.height])
    _panel_seed_dev(ax_d, ax_dp, gor, clus)
    ax_d.set_title(f"Gorilla OR6737, chr20 ({contig}; used to choose the seeding rule) · same-family pairs of the\n"
                   f"protein-homology families whose two genes both have ≥ 2 primary reads on their exons", loc="left",
                   fontsize=7.0, pad=6, linespacing=1.2)
    fig.text(0.01, 0.975, "Supplementary · seeding loci with secondary alignments vs primary alignments only, each "
             "through the default family rule; development contig, SECONDARY reference (protein-homology families)",
             fontsize=6.4, fontweight="bold", ha="left", va="top")
    fig.text(0.995, 0.008, "Development tables: the genome-wide supplement (Ensembl Compara, Liftoff copy pairs) "
             "replaces them", fontsize=5.3, color=figlib.INK_3, ha="right", va="bottom")
    figlib.stamp_provisional(fig, ["fig6_gorilla_recall", "fig6_gorilla_clusters"], data_dir)
    paths = figlib.save(fig, "fig6s_seeding", out_dir)
    plt.close(fig)
    return paths


def _plot_seeding_genome(data_dir: Path, out_dir: Path) -> list:
    """Genome-wide supplement: per sample, both configurations (open = primaries only, filled = default), black tick =
    the default without its largest family; Compara (human) >= 90% sensitivity and precision; Liftoff copy pairs
    (sequence_ID >= 0.95) sensitivity; protein-homology (secondary) when present. Substrate S1."""
    import matplotlib.pyplot as plt
    from matplotlib.lines import Line2D
    rows = figlib.read_table("fig6s_seeding", data_dir)
    order = [s for s in ["human_A119b", "human_testis", "gorilla_OR6737", "gorilla_KB3781", "chimp_PTR",
                         "orangutan_PPY"] if any(r["sample"] == s for r in rows)]
    order += sorted({r["sample"] for r in rows} - set(order))
    labels = {r["sample"]: _disp(r["sample_label"]) for r in rows}
    ok = {(r["sample"], r["reference"], r["config"], r["metric"], r["band"]): r for r in rows
          if r["status"] == "ok" and r["substrate"] == "S1"}
    status = {}
    for r in rows:
        if r["status"] != "ok":
            status.setdefault((r["sample"], r["reference"]), r["status"])
            status.setdefault((r["sample"], "-"), r["status"])
    panels = [("compara", "sensitivity", ">=90", "Compara pairs ≥ 90% protein identity\n(human; genes with ≥ 2 primary "
               "reads on their exons)", "Sensitivity"),
              ("compara", "precision", "all", "Compara precision (human; judgeable\nwithin-family pairs; upper bound)",
               "Precision"),
              ("liftoff", "sensitivity", "sc>=0.95", "Liftoff copy pairs (record, extra copy ≥ 95%;\nboth loci "
               "read-supported)", "Sensitivity")]
    if any(k[1] == "protein_homology" for k in ok):
        panels.append(("protein_homology", "sensitivity", ">=90", "SECONDARY: protein-homology families,\n≥ 90% "
                       "annotated-mRNA identity", "Sensitivity"))
    fig = plt.figure(figsize=(figlib.WIDTH_DOUBLE, 2.9))
    gs = fig.add_gridspec(1, len(panels), left=0.13, right=0.985, top=0.74, bottom=0.14, wspace=0.18)
    for j, (ref, metric, band, title, xlab) in enumerate(panels):
        ax = fig.add_subplot(gs[j])
        for i, s in enumerate(order):
            vals = {}
            for cfg_ in ("rustle_primary", "rustle"):
                r = ok.get((s, ref, cfg_, metric, band))
                if r is None or r["value"] in ("", None):
                    continue
                v = _f(r["value"])
                vals[cfg_] = v
                ax.plot([v], [i], linestyle="none", zorder=3, **(MK_PRIM if cfg_ == "rustle_primary" else MK_DEF))
                if cfg_ == "rustle" and r["n_without_top"] not in ("", None) and _i(r["n_without_top"]):
                    yl = _i(r["k_without_top"]) / _i(r["n_without_top"])
                    ax.plot([yl, yl], [i - 0.22, i + 0.22], color=figlib.INK, lw=1.0, zorder=4)
            if len(vals) == 2:
                ax.plot([vals["rustle_primary"], vals["rustle"]], [i, i], color=figlib.INK_3, lw=0.6, zorder=2)
            r = ok.get((s, ref, "rustle", metric, band)) or ok.get((s, ref, "rustle_primary", metric, band))
            if r is not None:
                ax.text(1.02, i, f"n={int(r['n']):,}", transform=ax.get_yaxis_transform(), fontsize=5.2,
                        color=figlib.INK_2, va="center")
            elif ref == "compara" and not s.startswith("human"):
                ax.text(0.02, i, "no Compara reference (human only)", fontsize=5.0, color=figlib.INK_3,
                        va="center", style="italic")
            else:
                ax.text(0.02, i, "not built", fontsize=5.0, color=figlib.INK_3, va="center", style="italic")
        ax.set_ylim(len(order) - 0.5, -0.5)
        ax.set_yticks(range(len(order)))
        ax.set_yticklabels([labels.get(s, s) for s in order] if j == 0 else [], fontsize=5.8)
        ax.set_xlim(-0.02, 1.02)
        ax.set_xticks([0, 0.5, 1.0])
        ax.set_xticklabels(["0", "0.5", "1"])
        ax.grid(axis="y", visible=False)
        ax.set_xlabel(xlab, fontsize=6.2)
        ax.set_title(title, loc="left", fontsize=6.0, pad=3, linespacing=1.15)
        figlib.panel_label(ax, "abcd"[j], x=-0.08, y=1.20)
    fig.legend(handles=[Line2D([], [], linestyle="none", label="primary alignments only", **MK_PRIM),
                        Line2D([], [], linestyle="none", label="default: + secondary alignments within 2% of the best "
                                                               "score", **MK_DEF),
                        Line2D([], [], color=figlib.INK, lw=1.0, label="default without its largest family")],
               loc="upper left", bbox_to_anchor=(0.13, 0.965), ncol=3, fontsize=5.6, frameon=False, handlelength=1.3)
    fig.text(0.01, 0.99, "Supplementary · seeding loci with secondary alignments, family level (genome minus development "
             "contigs; both configurations through the default family rule)", fontsize=6.6, fontweight="bold",
             ha="left", va="top")
    figlib.stamp_provisional(fig, ["fig6s_seeding"], data_dir)
    paths = figlib.save(fig, "fig6s_seeding", out_dir)
    plt.close(fig)
    return paths


def plot_seeding(data_dir: Path, out_dir: Path) -> list:
    if _have(data_dir, "fig6s_seeding"):
        return _plot_seeding_genome(Path(data_dir), Path(out_dir))
    if _have(data_dir, "fig6_gorilla_recall", "fig6_gorilla_clusters"):
        return _plot_seeding_dev(Path(data_dir), Path(out_dir))
    return []


def plot(data_dir: Path, out_dir: Path):
    import matplotlib.pyplot as plt
    figlib.use_style()
    main = (_plot_genome if _gw_ready(data_dir) else _plot_dev)(Path(data_dir), Path(out_dir))
    plt.close("all")
    return list(main) + plot_seeding(data_dir, out_dir) + plot_protein(Path(data_dir), Path(out_dir))


# ================================================================ caption numbers
def summary(data_dir: Path = figlib.DATA_DIR):
    """Print every number captions/fig6.md quotes (genome-wide tables when present, else the development tables), the
    pre-registered claims, and the supplements."""
    data_dir = Path(data_dir)
    if _gw_ready(data_dir):
        _summary_genome(data_dir)
    elif _have(data_dir, *DEV_TABLES):
        _summary_dev(data_dir)
    else:
        print("no fig6 tables: build with `make.py data fig6` (repeat while it exits 75)")
    if _have(data_dir, "fig6s_seeding"):
        print("== supplement fig6s_seeding")
        for r in figlib.read_table("fig6s_seeding", data_dir):
            if r["status"] != "ok":
                print(f"  {r['sample']:15s} {r['reference']:16s} {r['config']:15s} {r['status']}")
            elif r["band"] in (">=90", "all", "sc>=0.95"):
                print(f"  {r['sample']:15s} {r['substrate']} {r['reference']:16s} {r['config']:15s} {r['metric']:11s} "
                      f"{r['band']:9s} {r['k']}/{r['n']} = {r['value']}  top {r['top_family']} ({r['top_family_pairs']})"
                      f"  without {r['k_without_top']}/{r['n_without_top']}  covered {r['covered']}")
    elif _have(data_dir, "fig6_gorilla_recall"):
        print("== supplement fig6s_seeding (development, gorilla chr20, protein-homology families: secondary)")
        for r in figlib.read_table("fig6_gorilla_recall", data_dir):
            if r["metric"] in ("recall", "precision"):
                print(f"  {r['arm']:15s} {r['metric']:9s} {r['band']:6s} {r['k']}/{r['n']}")
    if _have(data_dir, "fig6s_protein_tiers"):
        print("== supplement fig6s_protein (chr16 development spectrum; T3 not part of the rule)")
        for r in figlib.read_table("fig6s_protein_tiers", data_dir):
            print(f"  {r['metric']:9s} {r['view']:18s} {r['band']:6s} {r['k']}/{r['n']}")


def _summary_dev(data_dir: Path):
    rec = figlib.read_table("fig6_chr16_recall", data_dir)
    grp = figlib.read_table("fig6_chr16_groups", data_dir)
    prec = figlib.read_table("fig6_chr16_precision", data_dir)
    print("== development: human A119b chr16, default families vs Compara")
    by = {(r["compara_band"], r["view"]): r for r in rec}
    for b in _o1.COMPARA_BANDS:
        f = by.get((b, "families"))
        if f:
            print(f"  {b:6s} n={f['n_pairs']:>4s} edge {by[(b, 'edge_nt')]['hits']:>3s} families {f['hits']:>3s} "
                  f"(with alignment {by[(b, 'families_edge')]['hits']}, without {by[(b, 'families_no_edge')]['hits']})"
                  f" ceiling {by[(b, 'families_both_members')]['hits']}")
    g = {(r["band"], r["stratum"], r["view"]): r for r in grp}
    for band in ("60-90", ">=90", "all"):
        for st in ("all", "top_group", "other_groups"):
            f, e = g.get((band, st, "families")), g.get((band, st, "edge_nt"))
            if f:
                print(f"  {band:6s} {st:13s} {f['group_label'][:36]:36s} n={f['n_pairs']} groups={f['n_groups']} edge "
                      f"{e['hits']} families {f['hits']} (families holding them: {f['families_with_hits']}; "
                      f"edge-not-family {g[(band, st, 'edge_nt_not_families')]['hits']}, family-not-edge "
                      f"{g[(band, st, 'families_no_edge')]['hits']})")
    for r in prec:
        if r["identity_band"] == "all":
            print(f"  precision {r['measure']:22s} {r['k']}/{r['n']} = {r['precision']}  largest {r['largest_family']} "
                  f"{r['largest_family_pairs']}  without {r['k_without_largest']}/{r['n_without_largest']}  "
                  f"copies {r['copies']} single-exon {r['single_exon_copies']}")


def _summary_genome(data_dir: Path):
    rec = figlib.read_table("fig6_gw_recall", data_dir)
    grp = figlib.read_table("fig6_gw_groups", data_dir)
    prec = figlib.read_table("fig6_gw_precision", data_dir)
    for n in figlib.table_meta("fig6_gw_recall", data_dir).get("note", []):
        if n.lower().startswith("provisional"):
            print("PROVISIONAL:", n)
    for sid in sorted({r["sample"] for r in rec if r["substrate"] == "S0"}):
        print(f"== {sid}")
        for sub in ("S0", "S1"):
            by = {(r["compara_band"], r["view"]): r for r in rec if r["sample"] == sid and r["substrate"] == sub}
            for b in _o1.COMPARA_BANDS:
                f = by.get((b, "families"))
                if not f:
                    continue
                e, c = by.get((b, "edge_nt")), by.get((b, "families_both_members"))
                print(f"  {sub} {b:6s} n={f['n_pairs']:>6s} edge {e['hits'] if e else '-':>5s} families {f['hits']:>5s} "
                      f"ceiling {c['hits'] if c else '-':>5s}  edge CI {e['ci_lo'] if e else ''}-{e['ci_hi'] if e else ''}"
                      f"  family CI {f['ci_lo']}-{f['ci_hi']}")
        g = {(r["band"], r["stratum"], r["view"]): r for r in grp if r["sample"] == sid}
        for band, st in (("60-90", "top_group"), ("60-90", "other_groups"), (">=90", "cross_chromosome"),
                         (">=90", "same_chromosome"), ("all", "cross_chromosome"), ("all", "same_chromosome")):
            f, e = g.get((band, st, "families")), g.get((band, st, "edge_nt"))
            if f:
                print(f"  {band} {st:17s} {f['group_label'][:40]:40s} n={f['n_pairs']} groups={f['n_groups']} "
                      f"edge {e['hits'] if e else '-'} families {f['hits']} (families: {f['families_with_hits']})")
        for r in prec:
            if r["sample"] == sid and r["identity_band"] == "all":
                print(f"  precision {r['measure']:22s} {r['k']}/{r['n']} = {r['precision']}  largest "
                      f"{r['largest_family']} {r['largest_family_pairs']}  without {r['k_without_largest']}/"
                      f"{r['n_without_largest']}  copies {r['copies']} single-exon {r['single_exon_copies']}")


if __name__ == "__main__":
    if sys.argv[1:2] == ["summary"]:
        summary(Path(sys.argv[2]) if len(sys.argv) > 2 else figlib.DATA_DIR)
    else:
        sys.exit("usage: python3 figures/fig_family_spectrum.py summary [DATA_DIR]")
