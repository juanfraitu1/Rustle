"""Figure 8 — loci in the Liftoff framework.

Pre-registration: docs/PREREG_liftoff_loci_2026-09-25.md (user decision 2026-09-25: every locus comparison is made in
the Liftoff framework). Data layer: figures/_liftoff.py.

  (a) The baseline, per species: Liftoff v1.6.3 lifts the genome's own RefSeq gene and pseudogene records onto the same
      genome with -copies (-sc 0.95; -a 0.5 -s 0.5 defaults). Records placed back in place, and the extra copies found,
      at sequence_ID >= 0.95 / 0.98 / 0.99 / 1.00.
  (b) Like for like (annotation + genome, no reads), per species: Liftoff's extra copies against the candidate loci of
      Rustle's guided-mode search seeded with every record (bench/guided_pipeline.py finders).
  (c) Rustle de novo loci (reads + genome) scored against Liftoff's loci as a reference, per sample: the fraction of
      read-supported annotated loci (in place) and extra copies that a de novo locus covers (Liftoff's -a 0.5 on exon
      bases), and the fraction of de novo loci that lie at a Liftoff locus.
  (d) The same for Rustle's DEFAULT de novo families (the driver's `families` stage copy table: one copy per member
      locus of a family, its representative's exons; prereg amendment 5). The legacy copy catalog (`catalog` stage)
      is scored by the same code as labelled secondary rows (rustle_set `catalog`) when it exists; it is not drawn.
  (e) Liftoff's (record, extra copy) pairs whose two loci are both covered by family members: the fraction in one
      default family (claim F1).
  (f) Missing-copy flags: how often the flagged record has a Liftoff extra copy (all scanned loci, fired loci,
      reference-absent candidates).

Matching rule everywhere (prereg §3): in one genome Liftoff's sequence_ID of a locus placed on another locus's exons
equals its coverage, so both criteria reduce to cov(G | R) = |exons(G) ∩ exons(R)| / |exons(G)| >= 0.5.

Tables (each listed in META only once it exists; a missing input is a `status` row, never a number):
  fig8_liftoff  per species: Liftoff classes by stratum, extra copies by sequence_ID threshold
  fig8_guided   per species: Liftoff extra copies found by Rustle guided candidates, and the reverse
  fig8_samples  per sample: de novo loci, default families (and the legacy catalog) measures (c-e), genome and genome
                minus development contigs
  fig8_flags    per sample: missing-copy flags against Liftoff's extra copies (f)

Build: `python3 figures/make.py data fig8 --set fig8_budget_s=540` under the heavy-process lock, repeated while it
exits 75 (V-L1 first, then Liftoff per species, read support per sample, optionally the guided search:
`--set fig8_guided_finder=1`). `--set fig8_plan=1` prints what is left.
"""
from __future__ import annotations

import collections
import json
import math
import sys
from pathlib import Path

import figlib

TABLES = ["fig8_liftoff", "fig8_guided", "fig8_samples", "fig8_flags"]
SPECIES_ORDER = ["human", "gorilla", "chimpanzee", "orangutan"]
STRATA = ["all", "protein_coding", "lncRNA", "pseudogene", "other"]
IDENTITY_ROWS = [0.80, 0.95]


def _existing_tables(data_dir: Path = figlib.DATA_DIR) -> list[str]:
    return [t for t in TABLES if (Path(data_dir) / f"{t}.tsv").exists()]


META = {
    "id": "fig8",
    "title": "Loci in the Liftoff framework: Liftoff's self-lift baseline, Rustle's guided search like for like, and "
             "Rustle's read-built loci, default families and missing-copy flags scored against it",
    "claim": ("Liftoff (v1.6.3, -copies -sc 0.95) re-places each genome's own RefSeq gene and pseudogene records and "
              "finds their unannotated extra copies: the annotation-guided locus baseline. Rustle's guided-mode "
              "candidate search (annotation + genome) is compared with it like for like; Rustle's de novo loci, its "
              "default de novo families (reads -> seeded loci -> one representative per locus -> families) and its "
              "missing-copy flags (reads + genome) are scored against it as a reference, never as a competitor. Every "
              "match uses Liftoff's own criterion (>= 50% of the reference locus's exon bases). "
              "Pre-registered claims L1, G1, D1 and F1 (docs/PREREG_liftoff_loci_2026-09-25.md); the numbers are "
              "printed by `python3 figures/fig_loci.py summary`."),
    "tables": _existing_tables(),
}

GEN = "figures/fig_loci.py"


# ================================================================ build
def _budget(cfg):
    import _liftoff as L
    return L.Budget(float(cfg.get("fig8_budget_s") or cfg.get("figs_budget_s") or 0))


def _plan_only(cfg) -> bool:
    return any(str(cfg.get(k, "")).strip().lower() in ("1", "true", "yes") for k in ("fig8_plan", "figs_plan"))


def _finder_on(cfg) -> bool:
    return str(cfg.get("fig8_guided_finder", "")).strip().lower() in ("1", "true", "yes")


def build(cfg: dict, data_dir: Path, force: bool = False):
    import _liftoff as L
    import samples
    if _plan_only(cfg):
        print("\n".join(L.plan_report(cfg)))
        return
    budget = _budget(cfg)
    vl1 = L.liftoff_root(cfg) / "vl1" / "vl1_report.json"
    if not vl1.exists():
        L.validate(cfg, budget)          # prereg §1: V-L1 before the genome-wide runs
    species = [s for s in SPECIES_ORDER if s in L.species_list(cfg)]
    for sp in species:
        L.ensure_species(cfg, sp, budget)
    reg = samples.registry(cfg)
    for sid, row in reg.items():
        L.ensure_support(cfg, sid, row["species"], budget)
    if _finder_on(cfg):
        for sp in species:
            L.ensure_finder(cfg, sp, budget)
    vl1_rep = json.loads(vl1.read_text())
    notes_common = [
        "Liftoff v1.6.3 (conda env liftoff; its minimap2 2.24), self-lift of each genome's RefSeq gene + pseudogene "
        f"records: -copies -sc {L.LIFT_SC} -a {L.LIFT_A} -s {L.LIFT_S} -overlap {L.LIFT_OVERLAP} -p {L.LIFT_THREADS} "
        f"-f pseudogene; default -mm2_options ({L.MM2_OPTIONS})",
        "run as one Liftoff call per contig's records (blocks of <= "
        f"{int(cfg.get('fig8_liftoff_block_records', L.LIFTOFF_BLOCK_RECORDS))} records) against the WHOLE genome; "
        "cross-call overlaps resolved by Liftoff's own rule (prereg §1 M1-M3)",
        "V-L1 (mini-genome chr20-22, one Liftoff run vs per-contig calls and vs blocks of records; prereg amendment "
        "3): " + "; ".join(
            f"{arm}: {r['annotated_A'] - len({x[0] for x in r['annotated_only_A']})} of {r['annotated_A']} annotated "
            f"placements identical, {r['copies_both']} of {len(set(map(tuple, r['copies_only_A'])) | set(map(tuple, r['copies_only_arm']))) + r['copies_both']} "
            f"extra copies identical (Jaccard {r['copies_jaccard']}) -> {'PASS' if r['pass'] else 'FAIL: reported as an approximation of one run'}"
            for arm, r in vl1_rep.items()),
        "matching (prereg §3): cov(G|R) = |exons(G) ∩ exons(R)| / |exons(G)| >= 0.5, same contig, strand ignored; in "
        "one genome Liftoff's sequence_ID of a coordinate match equals its coverage",
    ]
    _write_liftoff(cfg, data_dir, species, notes_common)
    _write_guided(cfg, data_dir, species, notes_common)
    _write_samples(cfg, data_dir, reg, notes_common)
    _write_flags(cfg, data_dir, reg, notes_common)
    META["tables"] = _existing_tables(data_dir)


def _write_liftoff(cfg, data_dir, species, notes):
    import _liftoff as L
    rows, inputs, extra = [], {}, []
    for sp in species:
        t = L.species_target(cfg, sp)
        loci = L.read_loci(t.loci)
        inputs[f"{sp}_liftoff_loci"] = t.loci
        inputs[f"{sp}_gff"] = t.gff
        cnt = collections.Counter()
        for r in loci:
            size = "lt200" if r["short"] else "ge200"
            if r["cls"] == "extra_copy":
                for sc in L.SC_ROWS:
                    if r["seqid_f"] >= sc - 1e-9:
                        for st in ("all", r["stratum"]):
                            cnt[("extra_copy", st, size, f"{sc:.2f}")] += 1
            else:
                for st in ("all", r["stratum"]):
                    cnt[(r["cls"], st, size, "-")] += 1
        for (cls, st, size, sc), n in sorted(cnt.items()):
            rows.append([sp, cls, st, size, sc, n])
        rep = json.loads((t.wdir / "merge_report.json").read_text())
        extra.append(f"{sp}: {rep['records']} records in {rep['shards']} calls; cross-call merge dropped "
                     f"{rep['M1_dropped']} copies overlapping another call's annotated locus (M1) and "
                     f"{rep['M2_dropped']} overlapping copies (M2); {rep['M3_overlaps']} annotated overlaps kept (M3); "
                     f"{rep['n_anomalies']} anomalies")
    figlib.write_table("fig8_liftoff", ["species", "cls", "stratum", "size", "sc_min", "n"], rows, generator=GEN,
                       inputs=inputs, data_dir=data_dir,
                       notes=notes + extra + [
                           "cls: in_place = placed on its own record (same contig, reciprocal span overlap >= 0.5); "
                           "moved = placed elsewhere; partial = Liftoff's partial_mapping / low_identity; unmapped; "
                           "extra_copy = -copies (unannotated by construction); dropped_M1/M2 = removed by the "
                           "cross-call rule",
                           "size: exon union >= 200 bp (ge200, scored for Rustle) or < 200 bp (lt200, reported only); "
                           "sc_min: extra copies with sequence_ID >= sc_min (filters of the -sc 0.95 run)"])


def _write_guided(cfg, data_dir, species, notes):
    import _liftoff as L
    rows, inputs = [], {}
    for sp in species:
        t = L.species_target(cfg, sp)
        cand = L.finder_dir(t) / "candidates.tsv"
        if not (_finder_on(cfg) and cand.exists()):
            rows.append([sp, "-", "-", "-", "", "", "", "", "",
                         "not run: make.py data fig8 --set fig8_guided_finder=1 (heavy; splice index per call)"])
            continue
        inputs[f"{sp}_candidates"] = cand
        loci = L.read_loci(t.loci)
        recs = L.read_records(t)
        spans = L.Index([{"contig": r["contig"], "iv": [(r["start0"], r["end"])]} for r in recs])
        cands = L.read_candidates(cand)
        copies = [r for r in loci if r["cls"] == "extra_copy"]
        for r in copies:
            r["inside"] = any(L.ov(spans.rows[i]["iv"][0][0], spans.rows[i]["iv"][0][1], r["iv"][0][0], r["iv"][-1][1]) > 0
                              for i in spans.near(r["contig"], [(r["iv"][0][0], r["iv"][-1][1])]))
        for idm in IDENTITY_ROWS:
            cs = [c for c in cands if c["identity"] >= idm - 1e-9]
            cidx = L.Index(cs)
            for sc in L.SC_ROWS:
                for item, sel in (("liftoff_extra_found", lambda r: not r["inside"]),
                                  ("liftoff_extra_inside_record", lambda r: r["inside"])):
                    ref = [r for r in copies if not r["short"] and r["seqid_f"] >= sc - 1e-9 and sel(r)]
                    hits = [int(cidx.best_cov(r["contig"], r["iv"])[0] >= L.MATCH_COV) for r in ref]
                    lo, hi = L.boot_ci([r["source_id"] for r in ref], hits)
                    rows.append([sp, item, f"{idm:.2f}", f"{sc:.2f}", len(ref), sum(hits),
                                 sum(hits) / len(ref) if ref else None, lo, hi, "ok"])
            gidx = L.Index(copies)
            hits = [int(gidx.best_cov(c["contig"], c["iv"])[0] >= L.MATCH_COV) for c in cs]
            lo, hi = L.boot_ci([c["source_id"] for c in cs], hits)
            rows.append([sp, "rustle_at_liftoff_extra", f"{idm:.2f}", "0.95", len(cs), sum(hits),
                         sum(hits) / len(cs) if cs else None, lo, hi, "ok"])
    figlib.write_table("fig8_guided", ["species", "item", "identity_min", "sc_min", "n", "k", "frac", "ci_lo", "ci_hi",
                                       "status"], rows, generator=GEN, inputs=inputs, data_dir=data_dir,
                       notes=notes + [
                           "Rustle guided candidates: bench/guided_pipeline.py finders with EVERY gene/pseudogene "
                           "record as a seed (unit -x splice, CDS envelope -x asm20, -c -N 100 -p 0.1, species splice "
                           "index; identity >= 0.80, >= 50% aligned; hits overlapping any record span blocked; "
                           "chain-first loci); a candidate's exons = its transcript hit's blocks, else its span",
                           "liftoff_extra_found: Liftoff extra copies (exon union >= 200 bp, outside every record's "
                           "span) covered >= 50% by a candidate with identity >= identity_min (claim G1 at 0.80); "
                           "liftoff_extra_inside_record: the same for extra copies inside another record's span "
                           "(Liftoff allows the other strand; the finder blocks them)",
                           "rustle_at_liftoff_extra: candidates with identity >= identity_min covered >= 50% by a "
                           "Liftoff extra copy (any size, sequence_ID >= 0.95)",
                           "95% intervals: 2,000 bootstrap resamples of source records (random seed 20260925)"])


def _dev(sp: str, contig: str) -> bool:
    import _liftoff as L
    return contig in L.DEV_CONTIGS.get(sp, set())


def _sample_measures(cfg, sid, sp, rset, rloci, loci, sup) -> list[list]:
    """Rows of fig8_samples for one Rustle locus set of one sample: `denovo` (every assembled locus), `families` (the
    default families' copy table) or `catalog` (the legacy copy catalog, secondary)."""
    import _liftoff as L
    out = []
    ref = [r for r in loci if r["cls"] in ("in_place", "moved", "extra_copy") and not r["short"]
           and sup.get(L.locus_key(r), 0) >= L.SUPPORT_READS]
    ridx = L.Index(rloci)
    cov = {id(r): ridx.best_cov(r["contig"], r["iv"])[0] for r in ref}
    for scope in ("genome", "minus_dev"):
        keep = (lambda c: True) if scope == "genome" else (lambda c: not _dev(sp, c))
        for measure, cls in (("sens_in_place", "in_place"), ("sens_moved", "moved"), ("sens_extra_copy", "extra_copy")):
            for sc in (L.SC_ROWS if cls == "extra_copy" else [None]):
                for st in STRATA:
                    sel = [r for r in ref if r["cls"] == cls and keep(r["contig"]) and (st == "all" or r["stratum"] == st)
                           and (sc is None or r["seqid_f"] >= sc - 1e-9)]
                    hits = [int(cov[id(r)] >= L.MATCH_COV) for r in sel]
                    lo, hi = L.boot_ci([r["source_id"] for r in sel], hits)
                    out.append([sid, sp, rset, scope, measure, st, "-" if sc is None else f"{sc:.2f}", len(sel),
                                sum(hits), sum(hits) / len(sel) if sel else None, lo, hi, "ok"])
        # located: Rustle loci at a Liftoff locus (cov(R | G) >= 0.5 for some placed Liftoff locus)
        gidx = L.Index([r for r in loci if r["cls"] in ("in_place", "moved", "extra_copy")])
        sel = [r for r in rloci if keep(r["contig"]) and r["iv"]]
        hits = [int(gidx.best_cov(r["contig"], r["iv"])[0] >= L.MATCH_COV) for r in sel]
        lo, hi = L.boot_ci([f"{r['contig']}:{r['locus']}" for r in sel], hits)
        out.append([sid, sp, rset, scope, "located", "all", "-", len(sel), sum(hits),
                    sum(hits) / len(sel) if sel else None, lo, hi, "ok"])
        if rset in ("families", "catalog"):
            out += _pairs(sid, sp, rset, rloci, loci, keep, scope)
    return out


def _pairs(sid, sp, rset, fam_loci, loci, keep, scope) -> list[list]:
    """F1: Liftoff (record, extra copy) pairs whose two loci are both covered >= 50% by family members; k = a covering
    member of one shares a family with a covering member of the other (_liftoff.copy_pairs / pair_families; no
    read-support condition here, as pre-registered for F1)."""
    import _liftoff as L
    out = []
    for sc in L.SC_ROWS:
        res = L.pair_families(L.copy_pairs(loci, sc, keep=keep), fam_loci)
        pairs = [(src, int(bool(shared))) for src, both, shared in res if both]
        lo, hi = L.boot_ci([p[0] for p in pairs], [p[1] for p in pairs])
        k = sum(p[1] for p in pairs)
        out.append([sid, sp, rset, scope, "pair_cofamily", "all", f"{sc:.2f}", len(pairs), k,
                    k / len(pairs) if pairs else None, lo, hi, "ok"])
    return out


def _write_samples(cfg, data_dir, reg, notes):
    import _liftoff as L
    rows, inputs = [], {}
    for sid, row in reg.items():
        sp = row["species"]
        lp = L.loci_path(cfg, sp)
        supp = L.support_path(cfg, sid)
        loci = L.read_loci(lp)
        sup = L.read_support(supp)
        inputs[f"{sid}_support"] = supp
        gtf, st = L.product_if_ready(cfg, sid, "assemble", "gtf")
        if gtf is None:
            rows.append([sid, sp, "denovo", "-", "-", "-", "-", "", "", "", "", "", st])
        else:
            inputs[f"{sid}_denovo_gtf"] = gtf
            rl = L.denovo_loci(gtf, L.liftoff_root(cfg) / "rustle_loci")
            rows += _sample_measures(cfg, sid, sp, "denovo", rl, loci, sup)
        # the default de novo families (families stage copy table), then the legacy catalog (secondary, table only)
        for rset, stage in (("families", "families"), ("catalog", "catalog")):
            cp, st = L.product_if_ready(cfg, sid, stage, "copies")
            if cp is None:
                rows.append([sid, sp, rset, "-", "-", "-", "-", "", "", "", "", "", st])
            else:
                inputs[f"{sid}_{rset}_copies"] = cp
                rows += _sample_measures(cfg, sid, sp, rset, L.catalog_loci(cp), loci, sup)
    figlib.write_table("fig8_samples", ["sample", "species", "rustle_set", "scope", "measure", "stratum", "sc_min", "n",
                                        "k", "frac", "ci_lo", "ci_hi", "status"], rows, generator=GEN, inputs=inputs,
                       data_dir=data_dir, notes=notes + [
                           "reference loci: Liftoff in_place / moved / extra_copy with exon union >= 200 bp and >= 2 "
                           "reads of the sample whose primary alignment (-F 2308) has an aligned block on the exon "
                           "union (the rule of Fig. 6d); sens_* = k of n covered >= 50% by one Rustle locus",
                           "denovo = gene_id groups of the sample's genome-wide assembly (run cache stage assemble; exon "
                           "union of every transcript); families = the DEFAULT de novo families' copy table (run cache "
                           "stage families, <id>.fam.copies.tsv: one copy per member locus of a family, its "
                           "representative's exons; prereg amendment 5); catalog = the LEGACY copy catalog (stage "
                           "catalog, its exons column; secondary, no claim, not drawn)",
                           "located = Rustle loci covered >= 50% by one placed Liftoff locus (the rest lie outside "
                           "every Liftoff locus: not an error); pair_cofamily = claim F1 (families; the catalog rows "
                           "are the legacy comparison)",
                           "scope minus_dev: without the development contigs (human chr16; gorilla chr20 = "
                           "NC_073244.2); chimpanzee and orangutan have none",
                           "95% intervals: 2,000 bootstrap resamples of source records (random seed 20260925)"])


def _write_flags(cfg, data_dir, reg, notes):
    import csv
    import _liftoff as L
    rows, inputs, review = [], {}, []
    for sid, row in reg.items():
        sp = row["species"]
        calls, st = L.product_if_ready(cfg, sid, "flag", "calls")
        if calls is None:
            rows.append([sid, sp, "-", "-", "", "", "", st])
            continue
        inputs[f"{sid}_flags"] = calls
        loci = L.read_loci(L.loci_path(cfg, sp))
        with open(calls) as fh:
            fl = list(csv.DictReader(fh, delimiter="\t"))
        ann_ids = {r["source_id"] for r in loci if r["cls"] in ("in_place", "moved", "partial")}
        for sc in L.SC_ROWS:
            cop = collections.defaultdict(list)
            allc = []
            for r in loci:
                if r["cls"] == "extra_copy" and r["seqid_f"] >= sc - 1e-9:
                    cop[r["source_id"]].append(r)
                    allc.append(r)
            aidx = L.Index(allc)
            annl = [r for r in loci if r["cls"] in ("in_place", "moved")]
            anidx = L.Index(annl)

            def other_covers(f, rows_, idx):
                o = f.get("other_locus", "NA")
                if o in ("NA", "", None):
                    return []
                c, _, se = o.rpartition(":")
                s, _, e = se.partition("-")
                s, e = int(s), int(e)
                return [rows_[i] for i in idx.near(c, [(s, e)])
                        if L.ov(s, e, rows_[i]["iv"][0][0], rows_[i]["iv"][-1][1])
                        >= 0.5 * (rows_[i]["iv"][-1][1] - rows_[i]["iv"][0][0])]
            sets = {"scanned": fl, "fired": [f for f in fl if f["status"].startswith("fired")],
                    "candidate": [f for f in fl if f.get("verdict") == "reference_absent_candidate"]}
            for name, sel in sets.items():
                k = sum(1 for f in sel if cop.get(f["locus"]))
                rows.append([sid, sp, f"{name}_with_copy", f"{sc:.2f}", len(sel), k, k / len(sel) if sel else None,
                             "ok"])
            with_copy = [f for f in sets["fired"] if cop.get(f["locus"])]
            home = [f for f in with_copy if any(x["source_id"] == f["locus"] for x in other_covers(f, allc, aidx))]
            rows.append([sid, sp, "fired_with_copy_home_is_copy", f"{sc:.2f}", len(with_copy), len(home),
                         len(home) / len(with_copy) if with_copy else None, "ok"])
            cands = [f for f in sets["candidate"] if cop.get(f["locus"])]
            rev = [f for f in cands if not any(x["source_id"] == f["locus"] for x in other_covers(f, allc, aidx))]
            rows.append([sid, sp, "candidate_with_copy_not_home", f"{sc:.2f}", len(cands), len(rev),
                         len(rev) / len(cands) if cands else None, "ok"])
            if sc == L.SC_ROWS[0]:
                review += [(sid, f["locus"], f["name"], f["chrom"], f["start"], f["end"], f.get("other_locus", "NA"),
                            ",".join(f"{x['contig']}:{x['start0']}-{x['end']}" for x in cop[f["locus"]]))
                           for f in rev]
            par = [f for f in sets["fired"] if f.get("verdict") == "unannotated_paralogue"]
            for item, test in (("paralogue_at_same_record_copy",
                                lambda f: any(x["source_id"] == f["locus"] for x in other_covers(f, allc, aidx))),
                               ("paralogue_at_any_copy", lambda f: bool(other_covers(f, allc, aidx))),
                               ("paralogue_at_other_annotated",
                                lambda f: any(x["source_id"] != f["locus"] for x in other_covers(f, annl, anidx)))):
                k = sum(1 for f in par if test(f))
                rows.append([sid, sp, item, f"{sc:.2f}", len(par), k, k / len(par) if par else None, "ok"])
    rp = None
    if review:
        import _liftoff as L
        rp = L.liftoff_root(cfg) / "flags_review.tsv"
        with open(rp, "w") as fo:
            fo.write("sample\tlocus\tname\tchrom\tstart\tend\tother_locus\tliftoff_extra_copies\n")
            for r in review:
                fo.write("\t".join(map(str, r)) + "\n")
    figlib.write_table("fig8_flags", ["sample", "species", "item", "sc_min", "n", "k", "frac", "status"], rows,
                       generator=GEN, inputs=inputs, data_dir=data_dir, notes=notes + [
                           "flag table: run cache stage flag (missing_copy.tsv; one row per annotated gene/pseudogene "
                           "record with >= 10 reads; `locus` = the record ID = Liftoff's source record)",
                           "<set>_with_copy: k of n loci whose record has >= 1 Liftoff extra copy with sequence_ID >= "
                           "sc_min (scanned = every row; fired = status fired*; candidate = verdict "
                           "reference_absent_candidate)",
                           "fired_with_copy_home_is_copy: the consensus's best other hit (other_locus) spans >= 50% of "
                           "one of that record's extra copies; candidate_with_copy_not_home: candidates whose record "
                           "has an extra copy that other_locus does not span (listed for review"
                           + (f": {rp})" if rp else ")"),
                           "paralogue_*: fired loci with verdict unannotated_paralogue whose other_locus spans >= 50% "
                           "of an extra copy of the same record / any extra copy / another record's placed locus "
                           "(prereg amendment 1)"])


# ================================================================ plot
def _f(x):
    try:
        v = float(x)
        return None if math.isnan(v) else v
    except (TypeError, ValueError):
        return None


def _na(ax, text):
    ax.text(0.5, 0.5, text, transform=ax.transAxes, ha="center", va="center", fontsize=6, color=figlib.INK_3,
            wrap=True)
    ax.set_xticks([])
    ax.set_yticks([])
    for s in ax.spines.values():
        s.set_visible(False)
    ax.grid(False)


SPECIES_LABEL = {"human": "Human (T2T-CHM13 v2.0)", "gorilla": "Gorilla (mGorGor1)",
                 "chimpanzee": "Chimpanzee (mPanTro3)", "orangutan": "Orangutan (mPonPyg2)"}


def _legend_below(ax, ncol=3):
    h, l = ax.get_legend_handles_labels()
    if h:
        ax.legend(h, l, loc="upper center", bbox_to_anchor=(0.5, -0.36), ncol=ncol, fontsize=5.4, handlelength=1.1,
                  columnspacing=1.0, borderaxespad=0.0)
GREY_RAMP = {"0.95": "#d6d5cf", "0.98": "#bdbcb6", "0.99": "#8a8983", "1.00": "#52514e"}


def _vl1_footnote(data_dir: Path) -> str | None:
    """The approximation statement for panel a when V-L1 missed its bar (prereg amendment 3)."""
    try:
        notes = figlib.table_meta("fig8_liftoff", data_dir).get("note", [])
    except FileNotFoundError:
        return None
    n = next((x for x in notes if x.startswith("V-L1")), None)
    if not n or "FAIL" not in n:
        return None
    import re
    m = re.search(r"C_blocks: (\d+) of (\d+) annotated placements identical, (\d+) of (\d+) extra copies", n)
    if not m:
        return "Liftoff ran once per block of records: an approximation of one run (see notes)"
    a, b, c, d = m.groups()
    return (f"Liftoff ran once per block of records, an approximation of one run: on chr20–22 it agrees with one run on "
            f"{int(a):,} of {int(b):,} placements and {int(c):,} of {int(d):,} extra copies (the rest: identical-copy "
            f"arrays of the acrocentric short arms)")


def _panel_liftoff(ax, rows):
    import numpy as np
    sp = [s for s in SPECIES_ORDER if any(r["species"] == s for r in rows)]
    y = np.arange(len(sp))[::-1]
    mx = 1
    labels = []
    for yi, s in zip(y, sp):
        rs = [r for r in rows if r["species"] == s]
        ex = {r["sc_min"]: int(r["n"]) for r in rs if r["cls"] == "extra_copy" and r["stratum"] == "all"
              and r["size"] == "ge200"}
        exs = {r["sc_min"]: int(r["n"]) for r in rs if r["cls"] == "extra_copy" and r["stratum"] == "all"
               and r["size"] == "lt200"}
        for sc in ("0.95", "0.98", "0.99", "1.00"):
            ax.barh(yi, ex.get(sc, 0), height=0.56, color=GREY_RAMP[sc], edgecolor=figlib.SURFACE, linewidth=0.6,
                    label=f"sequence_ID ≥ {sc}" if yi == y[0] else None)
        mx = max(mx, ex.get("0.95", 0))
        ann = collections.Counter()
        for r in rs:
            if r["cls"] in ("in_place", "moved", "partial", "unmapped") and r["stratum"] == "all":
                ann[r["cls"]] += int(r["n"])
        tot = sum(ann.values())
        labels.append(f"{SPECIES_LABEL[s]}\nin place {ann['in_place']:,} of {tot:,}" if tot else SPECIES_LABEL[s])
        ax.text(ex.get("0.95", 0), yi, f" {ex.get('0.95', 0):,} (+{exs.get('0.95', 0):,} < 200 bp)",
                va="center", ha="left", fontsize=5.4, color=figlib.INK_2)
    ax.set_yticks(y)
    ax.set_yticklabels(labels, fontsize=5.8)
    ax.set_xlim(0, mx * 1.6)
    ax.set_xlabel("Extra copies found by Liftoff -copies (exon union ≥ 200 bp)")
    ax.grid(axis="x", color=figlib.GRID, linewidth=0.5)
    ax.grid(axis="y", visible=False)
    _legend_below(ax, ncol=4)


def _panel_guided(ax, rows):
    import numpy as np
    ok = [r for r in rows if r["status"] == "ok"]
    if not ok:
        msg = next((r["status"] for r in rows), "not built")
        return _na(ax, "Rustle guided candidate search:\n" + msg)
    sp = [s for s in SPECIES_ORDER if any(r["species"] == s for r in ok)]
    y = np.arange(len(sp))[::-1] * 1.0
    blue_dark, blue_light, grey = figlib.BLUE[650], figlib.BLUE[250], "#d6d5cf"
    mx = 1
    for yi, s in zip(y, sp):
        f = next((r for r in ok if r["species"] == s and r["item"] == "liftoff_extra_found" and r["identity_min"] == "0.80"
                  and r["sc_min"] == "0.95"), None)
        c = next((r for r in ok if r["species"] == s and r["item"] == "rustle_at_liftoff_extra"
                  and r["identity_min"] == "0.95"), None)
        if f:
            n, k = int(f["n"]), int(f["k"])
            ax.barh(yi + 0.2, k, height=0.34, color=blue_dark, edgecolor=figlib.SURFACE, linewidth=0.6,
                    label="found by both" if yi == y[0] else None)
            ax.barh(yi + 0.2, n - k, left=k, height=0.34, color=grey, edgecolor=figlib.SURFACE, linewidth=0.6,
                    label="Liftoff only" if yi == y[0] else None)
            lo, hi = _f(f["ci_lo"]), _f(f["ci_hi"])
            ax.text(n, yi + 0.2, f" {k:,} of {n:,}" + (f" ({lo:.2f}–{hi:.2f})" if lo is not None else ""),
                    va="center", fontsize=5.4, color=figlib.INK_2)
            mx = max(mx, n)
        if c:
            n, k = int(c["n"]), int(c["k"])
            ax.barh(yi - 0.2, k, height=0.34, color=blue_dark, edgecolor=figlib.SURFACE, linewidth=0.6)
            ax.barh(yi - 0.2, n - k, left=k, height=0.34, color=blue_light, edgecolor=figlib.SURFACE, linewidth=0.6,
                    label="Rustle guided only" if yi == y[0] else None)
            ax.text(n, yi - 0.2, f" {k:,} of {n:,}", va="center", fontsize=5.4, color=figlib.INK_2)
            mx = max(mx, n)
    ax.set_yticks(list(y + 0.2) + list(y - 0.2))
    ax.set_yticklabels([f"{s.capitalize()} · Liftoff" for s in sp]
                       + [f"{s.capitalize()} · Rustle ≥ 95%" for s in sp], fontsize=5.6)
    ax.set_xlim(0, mx * 1.6)
    ax.set_xlabel("New loci (outside every annotated record)")
    ax.grid(axis="x", color=figlib.GRID, linewidth=0.5)
    ax.grid(axis="y", visible=False)
    _legend_below(ax, ncol=3)


def _sample_axis(ax, sids, labels):
    import numpy as np
    y = np.arange(len(sids))[::-1]
    ax.set_yticks(y)
    ax.set_yticklabels([labels.get(s, s) for s in sids], fontsize=5.8)
    ax.set_ylim(-0.6, len(sids) - 0.4)
    ax.grid(axis="x", color=figlib.GRID, linewidth=0.5)
    ax.grid(axis="y", visible=False)
    return dict(zip(sids, y))


def _short(st: str) -> str:
    """A status row's message, shortened for a panel: 'families run: run `python3 ...`' -> 'families: not built yet'."""
    if "`" in st:
        return st.split(":", 1)[0].split()[0] + ": not built yet (see table status)"
    return st[:60]


def _dot(ax, x, y, lo, hi, **kw):
    if x is None:
        return
    if lo is not None and hi is not None:
        ax.plot([lo, hi], [y, y], color=kw.get("color", figlib.INK_2), linewidth=0.9, solid_capstyle="butt", zorder=2)
    ax.plot([x], [y], linestyle="none", zorder=3, **kw)


SAMPLE_ORDER = ["human_A119b", "human_testis", "gorilla_OR6737", "gorilla_KB3781", "chimp_PTR", "orangutan_PPY"]
SAMPLE_LABEL = {"human_A119b": "Human A119b", "human_testis": "Human, testis", "gorilla_OR6737": "Gorilla OR6737, testis",
                "gorilla_KB3781": "Gorilla KB3781, fibroblast", "chimp_PTR": "Chimpanzee PTR",
                "orangutan_PPY": "Orangutan PPY"}


def _panel_sens(ax, rows, rset, title_what):
    sel = [r for r in rows if r["rustle_set"] == rset]
    sids = [s for s in SAMPLE_ORDER if any(r["sample"] == s for r in sel)]
    if not sids or not any(r["status"] == "ok" for r in sel):
        msg = next((r["status"] for r in sel if r["status"] != "ok"), "not built")
        return _na(ax, f"{title_what}:\n{msg}")
    lab = {}
    for s in sids:
        rs = {(r["measure"], r["sc_min"]): r for r in sel if r["sample"] == s and r["scope"] == "genome"
              and r["stratum"] == "all" and r["status"] == "ok"}
        a, e = rs.get(("sens_in_place", "-")), rs.get(("sens_extra_copy", "0.95"))
        lab[s] = SAMPLE_LABEL.get(s, s) + (f"\nn = {int(a['n']):,} · {int(e['n']):,}" if a and e else "")
    ys = _sample_axis(ax, sids, lab)
    blue = figlib.TOOL_COLOR["rustle"]
    first = True
    for s in sids:
        rs = {(r["measure"], r["sc_min"]): r for r in sel if r["sample"] == s and r["scope"] == "genome"
              and r["stratum"] == "all" and r["status"] == "ok"}
        if not rs:
            st = next((r["status"] for r in sel if r["sample"] == s), "")
            ax.text(0.02, ys[s], _short(st), fontsize=5, color=figlib.INK_3, va="center",
                    transform=ax.get_yaxis_transform())
            continue
        a = rs.get(("sens_in_place", "-"))
        e = rs.get(("sens_extra_copy", "0.95"))
        loc = rs.get(("located", "-"))
        if a:
            _dot(ax, _f(a["frac"]), ys[s] + 0.14, _f(a["ci_lo"]), _f(a["ci_hi"]), marker="o", color=blue,
                 markerfacecolor=blue, markeredgecolor=figlib.SURFACE, markersize=5,
                 label="annotated loci, in place" if first else None)
        if e:
            _dot(ax, _f(e["frac"]), ys[s] - 0.14, _f(e["ci_lo"]), _f(e["ci_hi"]), marker="D", color=blue,
                 markerfacecolor=figlib.SURFACE, markeredgecolor=blue, markeredgewidth=1.1, markersize=4.2,
                 label="Liftoff extra copies" if first else None)
        if loc:
            _dot(ax, _f(loc["frac"]), ys[s], None, None, marker="s", color=figlib.INK_3, markerfacecolor=figlib.INK_3,
                 markeredgecolor=figlib.SURFACE, markersize=3.8,
                 label=f"{title_what} at a Liftoff locus" if first else None)
        first = False
    ax.set_xlim(0, 1.0)
    ax.set_xlabel("Fraction covered ≥ 50% (Liftoff's -a on exon bases)")
    _legend_below(ax, ncol=3)


def _panel_pairs(ax, rows, rset="families"):
    sel = [r for r in rows if r["rustle_set"] == rset and r["measure"] == "pair_cofamily" and r["scope"] == "genome"
           and r["sc_min"] == "0.95" and r["status"] == "ok"]
    if not sel:
        msg = next((r["status"] for r in rows if r["rustle_set"] == rset and r["status"] != "ok"), "not built")
        return _na(ax, "Default families vs Liftoff copy pairs:\n" + msg)
    sids = [s for s in SAMPLE_ORDER if any(r["sample"] == s for r in sel)]
    ys = _sample_axis(ax, sids, {r["sample"]: f"{SAMPLE_LABEL.get(r['sample'], r['sample'])}\n{int(r['k']):,} of "
                                              f"{int(r['n']):,} pairs" for r in sel})
    blue = figlib.TOOL_COLOR["rustle"]
    for r in sel:
        _dot(ax, _f(r["frac"]), ys[r["sample"]], _f(r["ci_lo"]), _f(r["ci_hi"]), marker="o", color=blue,
             markerfacecolor=blue, markeredgecolor=figlib.SURFACE, markersize=5)
    ax.axvline(0.90, color=figlib.INK_3, linewidth=0.7, zorder=1)
    ax.text(0.90, 1.0, "pre-registered bar 0.90", transform=ax.get_xaxis_transform(), fontsize=5.2, color=figlib.INK_3,
            ha="right", va="bottom")
    ax.set_xlim(0, 1.0)
    ax.set_xlabel("Liftoff (record, extra copy) pairs in one default family")


def _panel_flags(ax, rows):
    ok = [r for r in rows if r["status"] == "ok" and r["sc_min"] == "0.95"]
    if not ok:
        msg = next((r["status"] for r in rows if r["status"] != "ok"), "not built")
        return _na(ax, "Missing-copy flags vs Liftoff extra copies:\n" + msg)
    sids = [s for s in SAMPLE_ORDER if any(r["sample"] == s for r in rows)]
    lab = {}
    for s in sids:
        rv = next((r for r in ok if r["sample"] == s and r["item"] == "candidate_with_copy_not_home"), None)
        lab[s] = SAMPLE_LABEL.get(s, s) + (f"\n{int(rv['k']):,} candidates to review" if rv else "")
    ys = _sample_axis(ax, sids, lab)
    blue = figlib.TOOL_COLOR["rustle"]
    spec = [("scanned_with_copy", "all scanned loci", dict(marker="o", color=figlib.INK_3, markerfacecolor=figlib.SURFACE,
                                                            markeredgecolor=figlib.INK_3, markeredgewidth=1.0,
                                                            markersize=4.2), 0.18),
            ("fired_with_copy", "fired loci", dict(marker="^", color=blue, markerfacecolor=blue,
                                                   markeredgecolor=figlib.SURFACE, markersize=5), 0.0),
            ("candidate_with_copy", "reference-absent candidates", dict(marker="D", color=blue,
                                                                        markerfacecolor=figlib.SURFACE,
                                                                        markeredgecolor=blue, markeredgewidth=1.1,
                                                                        markersize=4.2), -0.18)]
    first = {k: True for k, *_ in spec}
    for s in sids:
        rs = {r["item"]: r for r in ok if r["sample"] == s}
        if not rs:
            st = next((r["status"] for r in rows if r["sample"] == s), "")
            ax.text(0.02, ys[s], _short(st), fontsize=5, color=figlib.INK_3, va="center",
                    transform=ax.get_yaxis_transform())
            continue
        for item, lab, kw, dy in spec:
            r = rs.get(item)
            if r and _f(r["frac"]) is not None:
                _dot(ax, _f(r["frac"]), ys[s] + dy, None, None, label=lab if first[item] else None, **kw)
                first[item] = False
    ax.set_xlim(0, 1.0)
    ax.set_xlabel("Fraction whose record has a Liftoff extra copy")
    _legend_below(ax, ncol=3)


def plot(data_dir: Path, out_dir: Path):
    import matplotlib.pyplot as plt
    figlib.use_style()
    data_dir = Path(data_dir)
    have = set(_existing_tables(data_dir))
    tab = {t: (figlib.read_table(t, data_dir) if t in have else []) for t in TABLES}
    fig = plt.figure(figsize=(figlib.WIDTH_DOUBLE, 178 * figlib.MM))
    gs = fig.add_gridspec(3, 2, left=0.19, right=0.97, top=0.95, bottom=0.1, hspace=0.95, wspace=0.62)
    axes = {k: fig.add_subplot(gs[i, j]) for k, (i, j) in zip("abcdef", [(0, 0), (0, 1), (1, 0), (1, 1), (2, 0),
                                                                         (2, 1)])}
    if tab["fig8_liftoff"]:
        _panel_liftoff(axes["a"], tab["fig8_liftoff"])
        fn = _vl1_footnote(data_dir)
        if fn:
            import textwrap
            axes["a"].text(-0.02, -0.5, "\n".join(textwrap.wrap(fn, 105)), transform=axes["a"].transAxes,
                           fontsize=5.0, color=figlib.INK_3, ha="left", va="top")
    else:
        _na(axes["a"], "Liftoff self-lift: not built\n(make.py data fig8)")
    _panel_guided(axes["b"], tab["fig8_guided"]) if tab["fig8_guided"] else _na(axes["b"], "not built")
    _panel_sens(axes["c"], tab["fig8_samples"], "denovo", "Rustle de novo loci") if tab["fig8_samples"] \
        else _na(axes["c"], "not built")
    _panel_sens(axes["d"], tab["fig8_samples"], "families", "family members") if tab["fig8_samples"] \
        else _na(axes["d"], "not built")
    _panel_pairs(axes["e"], tab["fig8_samples"]) if tab["fig8_samples"] else _na(axes["e"], "not built")
    _panel_flags(axes["f"], tab["fig8_flags"]) if tab["fig8_flags"] else _na(axes["f"], "not built")
    titles = {"a": "Liftoff self-lift (annotation + genome), per species",
              "b": "Like for like: Liftoff -copies vs Rustle guided search",
              "c": "Rustle de novo loci (reads) vs Liftoff loci, per sample",
              "d": "Rustle default families (their member loci) vs Liftoff loci",
              "e": "Default families contain Liftoff's copy pairs",
              "f": "Missing-copy flags vs Liftoff extra copies"}
    for k, ax in axes.items():
        figlib.panel_label(ax, k, x=-0.02, y=1.08)
        ax.text(0.03, 1.08, titles[k], transform=ax.transAxes, fontsize=6.8, fontweight="bold", va="bottom",
                ha="left")
    figlib.stamp_provisional(fig, TABLES, data_dir)
    paths = figlib.save(fig, "fig8_loci", out_dir)
    plt.close(fig)
    return paths


# ================================================================ caption numbers
def summary(data_dir: Path = figlib.DATA_DIR):
    """Print the pre-registered claims L1, G1, D1, F1 and every caption number from the tables."""
    have = set(_existing_tables(data_dir))
    if "fig8_liftoff" in have:
        rows = figlib.read_table("fig8_liftoff", data_dir)
        for sp in SPECIES_ORDER:
            rs = [r for r in rows if r["species"] == sp and r["stratum"] == "all"]
            if not rs:
                continue
            c = collections.Counter()
            for r in rs:
                if r["cls"] != "extra_copy":
                    c[r["cls"]] += int(r["n"])
            tot = sum(c[k] for k in ("in_place", "moved", "partial", "unmapped"))
            frac = c["in_place"] / tot if tot else float("nan")
            ex = {(r["sc_min"], r["size"]): int(r["n"]) for r in rs if r["cls"] == "extra_copy"}
            print(f"L1 {sp}: in place {c['in_place']:,} of {tot:,} = {frac:.4f} (bar 0.990: "
                  f"{'PASS' if frac >= 0.99 else 'FAIL'}); moved {c['moved']}, partial {c['partial']}, unmapped "
                  f"{c['unmapped']}; extra copies >= 200 bp at sc 0.95/0.98/0.99/1.00: "
                  + "/".join(str(ex.get((s, 'ge200'), 0)) for s in ("0.95", "0.98", "0.99", "1.00"))
                  + f" (+{ex.get(('0.95', 'lt200'), 0)} < 200 bp)")
    if "fig8_guided" in have:
        for r in figlib.read_table("fig8_guided", data_dir):
            if r["status"] != "ok":
                print(f"G1 {r['species']}: {r['status']}")
            elif r["item"] == "liftoff_extra_found" and r["identity_min"] == "0.80" and r["sc_min"] == "0.95":
                fr = _f(r["frac"])
                print(f"G1 {r['species']}: {r['k']} of {r['n']} = {fr:.3f} (95% {r['ci_lo']}-{r['ci_hi']}; bar 0.90: "
                      f"{'PASS' if fr is not None and fr >= 0.90 else 'FAIL'})")
    if "fig8_samples" in have:
        rows = figlib.read_table("fig8_samples", data_dir)
        for sid in SAMPLE_ORDER:
            for scope in ("genome", "minus_dev"):
                idx = {(r["rustle_set"], r["measure"], r["sc_min"]): r for r in rows if r["sample"] == sid
                       and r["scope"] == scope and r["stratum"] == "all" and r["status"] == "ok"}
                a, e = idx.get(("denovo", "sens_in_place", "-")), idx.get(("denovo", "sens_extra_copy", "0.95"))
                if a and e and _f(a["frac"]) is not None and _f(e["frac"]) is not None:
                    ok = _f(e["frac"]) >= _f(a["frac"]) - 0.10
                    print(f"D1 {sid} {scope}: extra copies {_f(e['frac']):.3f} ({e['k']} of {e['n']}) vs in place "
                          f"{_f(a['frac']):.3f} ({a['k']} of {a['n']}): {'PASS' if ok else 'FAIL'}")
                p = idx.get(("families", "pair_cofamily", "0.95"))
                if p and _f(p["frac"]) is not None:
                    print(f"F1 {sid} {scope}: {p['k']} of {p['n']} = {_f(p['frac']):.3f} (bar 0.90: "
                          f"{'PASS' if _f(p['frac']) >= 0.90 else 'FAIL'})")
                p = idx.get(("catalog", "pair_cofamily", "0.95"))
                if p and _f(p["frac"]) is not None:
                    print(f"   legacy catalog (secondary, no claim) {sid} {scope}: {p['k']} of {p['n']} = "
                          f"{_f(p['frac']):.3f}")
            for r in rows:
                if r["sample"] == sid and r["status"] != "ok":
                    print(f"   {sid} {r['rustle_set']}: {r['status']}")
    if "fig8_flags" in have:
        for r in figlib.read_table("fig8_flags", data_dir):
            if r["status"] != "ok":
                print(f"flags {r['sample']}: {r['status']}")
            elif r["sc_min"] == "0.95":
                print(f"flags {r['sample']} {r['item']}: {r['k']} of {r['n']}")


if __name__ == "__main__":
    sys.path.insert(0, str(Path(__file__).resolve().parent))
    if sys.argv[1:2] == ["summary"]:
        summary(Path(sys.argv[2]) if len(sys.argv) > 2 else figlib.DATA_DIR)
    else:
        sys.exit("usage: python3 figures/fig_loci.py summary [DATA_DIR]")
