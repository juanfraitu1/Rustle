#!/usr/bin/env python3
"""family_copy_number.py -- per-family O1 COPY NUMBER: chi(H) counts, Strong-Sep assigns.

Advisor (Canzar) refinement of the E_r RNA-family reframe. The read-CONFLICT signal --
demoted from the O1 family *boundary* by the E_r reframe -- returns to O1 as a derived
family-level *multiplicity*: the copy number K(C) = chi(H_C) = MCC (THEORY.md Lemma 1).
This does NOT reverse the reframe (RNA family membership stays transcript-homology E_r);
it ADDS "how many copies", INCLUDING copies no read can be assigned to.

THE COUNT-vs-ASSIGN DISTINCTION (theorem-level, THEORY.md Lemma 1 / Thm 2 / Thm 4):
  * COUNTING copies = chi(H) -- a functional of the conflict relation ALONE; well-defined
    EVEN WHEN not a single read is assignable. On the de-tied / significance-gated
    (copy-consensus) graph and in the full-span regime, chi(H) is a LOWER BOUND on the true
    copy number (on the RAW allele-disagreement read graph error edges inflate chi to ~3xK,
    so the gate/copy-consensus is a precondition for the bound -- THEORY.md §3 Lemma).
  * ASSIGNING is strictly harder and splits into TWO theorems the shorthand must keep apart:
    (i) SINGLE-READ assignability of a read r (|N(r)|=1) needs a single-position
        Strong-Separation witness = a SUN = Theorem 4(ii) / the Bridge -- the per-read
        certificate, and exactly the Tier-1 property (per-copy, coverage-independent);
    (ii) family-wide UNIQUE-COVER recovery needs Strong Separation (Theorem 2), a
        COVERAGE-dependent condition. So SUN / Tier-1 is NOT the same object as Strong-Sep
        (an all-Tier-1 family with short reads can still fail Strong Separation).
  => a family is COUNTED as K = chi(H) copies while only the Tier-1 subset is single-read
     ASSIGNABLE (Thm 4).

3-TIER copy hierarchy (per copy, from bench/sun_identifiability.*), and its counting semantics:
  Tier1 RESOLVABLE            SUN-identifiable: its own color in chi(H) AND single-read
                              assignable (the O2 nice core).
  Tier2 DISTINGUISHABLE-BUT-  unique hap-vector, NO SUN: its own color in chi(H) but NOT
        UNASSIGNABLE          single-read assignable (needs >=2-PSV co-observation) --
                              the advisor's "counted but not resolvable" copies.
  Tier3 IDENTICAL / collapsed shares a hap-vector, no distinguishing PSV: chi(H) COLLAPSES
                              the group to ONE color -> chi(H) UNDER-counts by (size-1) ->
                              countable only by read-DEPTH / DNA (K=0 frontier, parCN/O4).

Per-family metrics (all machine-checked identities; fails loudly if a source drifts):
  n_loci                 = co-located reference copies fed to the graph = Tier1+Tier2+Tier3.
  chi_H (= K)            = # DISTINCT pairwise-conflicting copy hap-vectors (colors)
                          = copy-consensus copyonly_K = #singleton_groups + #nonsingleton_groups.
  n_resolvable           = Tier1  (counted AND single-read assignable).
  n_counted_unassignable = chi_H - n_resolvable = Tier2 + #nonsingleton-group color-reps
                          (copies chi(H) counts but a single read cannot pin -- advisor's point).
  collapsed_excess       = Tier3 - #nonsingleton_groups = SUM(size-1) over collapse groups
                          = identical copies chi(H) MISSES (need depth/DNA).
  true_copy_lower_bound  = chi_H + collapsed_excess = n_loci (reference-resolved regime).

HONEST BOUND CHAIN (verified genome-wide; REFERENCE-RESOLVED regime -- minimap2 pre-splits
distinct reference loci, so identical copies collapse and chi(H) UNDER-counts vs n_loci):
  n_resolvable <= chi_H <= n_loci = true_copy_lower_bound <= true_copy_number,
  last gap = reference-absent copies (O4). GENERAL FORM unifying both regimes:
  max(n_loci, chi_H) <= true_copy_number; the reference-COLLAPSED (SDA/Vollger) regime has
  the two middle terms swapped (n_loci <= chi_H: one locus hides several hap-vectors, so
  reads reveal MORE copies than the reference shows). Do NOT state "n_loci <= chi(H)" flatly
  -- on THIS reference-resolved GGO substrate it is chi(H)=361 <= n_loci=412.

THIRD PROBE -- read-DEPTH (folded in 2026-07; purely ADDITIVE, new columns only). The O1
cardinality stack is now THREE INDEPENDENT lower-bound probes of the same true_copy_number,
each counting a strictly different copy set:  n_loci  ->  chi_H  ->  depth_cn.
  depth_cn = E_fam / lambda_within, where E_fam = total reads over the family's (deduped)
  reference copies and lambda_within = median per-reference-copy read count WITHIN the family
  (RNA-only self-normalization a la Sudmant-2010 parCN). This is the SAME E_fam/lambda_within
  as the validated rna_copy_number_depth.py (depth_cn_within), reusing
  psv_graph_genomewide.dedup_copies so the per-copy read piles are the identical pipeline path.
  Its unambiguous GGO win is over chi_H, which is BLIND to Tier-3 identical copies: on the
  collapsed_excess>0 families depth_cn > chi_H, recovering exactly the identical copies chi(H)
  collapses away; its genuinely-novel value is the reference-COLLAPSED regime where BOTH n_loci
  and chi_H fail. HONEST CAVEAT (load-bearing): read depth counts EXPRESSED copies, so depth_cn
  is a NOISY LOWER BOUND on genomic CN, NOT an unbiased count -- SILENT copies (expression ~0)
  are INVISIBLE => biased low, and per-copy expression skew is two-sided noise. New columns
  (ADDITIVE; the nine shipped columns above are byte-identical / unchanged):
    depth_cn                        expressed-copy estimate = E_fam / lambda_within.
    regime                          reference_resolved (chi_H <= n_loci) vs reference_collapsed
                                    (n_loci < chi_H: one locus hides several hap-vectors, so
                                    depth is the only multiplicity signal). On this GGO
                                    reference-resolved substrate every family is resolved.
    depth_recovers_collapsed_excess 1 iff collapsed_excess>0 AND depth_cn > chi_H (depth counts
                                    the Tier-3 identical copies chi(H) missed; validated on all
                                    21 collapsed families, e.g. LOC115930538 chi_H=1->depth~11.4).

WHICH GRAPH IS chi(H). chi_H is the COPY-CONSENSUS conflict graph (copyonly_K): copies
distinguishable in principle from the assembled copy sequences, NO read-support gate --
the correct O1 COUNT (it counts a copy even when no read resolved it: the unassignable ones).
The read-level psv_graph K (chi_H_readlevel) is the ASSIGNMENT-realized count (what THIS
read set happened to resolve); it is carried as a cross-check and diverges on 11/154 families
(count-vs-assign made concrete: consensus counts a copy the reads never co-observed).

VERIFY. For a couple families we reconstruct the actual per-copy hap-vectors (reusing the
copyonly asm20 path) and confirm chi_H == #distinct hap-vectors AND every pair of distinct
hap-vectors differs at >=1 PSV column (pairwise conflict => H is complete multipartite =>
chi(H) = #distinct, Lemma C2). So we do NOT recompute chi(H) from scratch; we reuse the
copy-consensus catalog value (copyonly_K) and check it == #distinct-pairwise-conflicting
hap-vectors reconstructed from sequence.

Sources (keyed by integer family id):
  bench/psv_graph_genomewide.json  read-level K, cls   (cross-check + class label)
  bench/sun_identifiability.json   copyonly_K (=chi_H), tier1/2/3, n_singleton_groups
  psv_graph_genomewide.FAM_TSV     per-locus nreads for the depth leg (deduped via
                                   dedup_copies); optional -- if the module/data is
                                   unavailable, depth_cn is left null and the 3 depth
                                   columns degrade gracefully (the shipped 9 are unaffected).
Run from repo root:  /home/juanfra/miniforge3/bin/python bench/family_copy_number.py
Writes bench/family_copy_number.tsv and bench/family_copy_number.json.
Env: FCN_VERIFY_FAMS="0,32,42" families to hap-verify (empty to skip the BAM step).
"""
import collections
import csv
import itertools
import json
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
PG_JSON = os.path.join(HERE, "psv_graph_genomewide.json")
SUN_JSON = os.path.join(HERE, "sun_identifiability.json")
OUT_TSV = os.environ.get("FCN_OUT_TSV", os.path.join(HERE, "family_copy_number.tsv"))
OUT_JSON = os.environ.get("FCN_OUT_JSON", os.path.join(HERE, "family_copy_number.json"))
VERIFY_FAMS = [int(x) for x in os.environ.get("FCN_VERIFY_FAMS", "0,32,42").split(",") if x.strip()]

COLS = ["family", "gene", "n_loci", "chi_H", "n_resolvable", "n_counted_unassignable",
        "collapsed_excess", "true_copy_lower_bound", "cls",
        "chi_H_readlevel", "readlevel_disagrees",
        # --- ADDITIVE third leg (depth); the columns above are byte-identical/unchanged ---
        "depth_cn", "regime", "depth_recovers_collapsed_excess"]


def load():
    pg = {f["family"]: f for f in json.load(open(PG_JSON))["families"]}
    sun = {f["family"]: f for f in json.load(open(SUN_JSON))["families"]}
    return pg, sun


def median(xs):
    """Exact match to rna_copy_number_depth.median (lower-of-two for even n uses the pair mean)."""
    s = sorted(xs)
    n = len(s)
    if n == 0:
        return 0.0
    m = n // 2
    return float(s[m]) if n % 2 else (s[m - 1] + s[m]) / 2.0


def load_depth_cn():
    """Per-family depth_cn = E_fam / lambda_within, REPLICATING rna_copy_number_depth.py exactly
    (E_fam = sum of deduped per-copy read counts; lambda_within = their median) so the value here
    equals that script's `depth_cn_within`. We replicate the tiny arithmetic INLINE rather than
    importing rna_copy_number_depth, because that module reads family_copy_number.tsv (our OWN
    output) as input -- importing it would make this script's output depend on a stale copy of the
    file it is about to write (a fragile self-cycle). We DO reuse psv_graph_genomewide.dedup_copies
    (the shared, deterministic copy-dedup path) via a lazy import guarded like verify_hapvectors():
    if the module/FAM_TSV is unavailable the depth leg degrades to null and the shipped columns are
    untouched. Returns {family: depth_cn or None}."""
    try:
        sys.path.insert(0, HERE)
        import psv_graph_genomewide as pgmod
    except Exception:  # pragma: no cover - environment without pysam/data
        return {}
    if not os.path.exists(pgmod.FAM_TSV):
        return {}
    loci = collections.defaultdict(list)
    with open(pgmod.FAM_TSV) as fh:
        next(fh)
        for line in fh:
            fi, lid, c, s, e, nr = line.rstrip("\n").split("\t")
            loci[int(fi)].append((lid, c, int(s), int(e), int(nr)))
    out = {}
    for fam, rows in loci.items():
        ded = pgmod.dedup_copies(rows)            # [(chrom,start,end,reads_summed)] == n_loci copies
        depths = [r[3] for r in ded]
        e_fam = sum(depths)
        lam_within = median(depths) if depths else 0.0
        dcn = (e_fam / lam_within) if lam_within > 0 else float("nan")
        out[fam] = round(dcn, 2) if dcn == dcn else None   # NaN -> None
    return out


def family_row(fam, pg, sun):
    """Copy-number metric row for one family; every identity asserted (machine-checked)."""
    s, p = sun[fam], pg[fam]
    n_loci = s["n_copies"]                                   # Tier1+Tier2+Tier3
    t1, t2, t3 = s["tier1"], s["tier2"], s["tier3"]
    assert t1 + t2 + t3 == n_loci, (fam, "tier-sum", t1, t2, t3, n_loci)

    chi_H = s["copyonly_K"]                                  # colors = distinct hap-vectors
    nsg = s["n_singleton_groups"]                            # Tier1 + Tier2 (own-color copies)
    assert nsg == t1 + t2, (fam, "singleton", nsg, t1, t2)
    n_nonsingleton_groups = chi_H - nsg                      # collapse-group color reps
    assert n_nonsingleton_groups >= 0, (fam, "nnsg", n_nonsingleton_groups)

    collapsed_excess = t3 - n_nonsingleton_groups           # SUM(size-1); identical copies missed
    assert collapsed_excess >= 0, (fam, "excess", collapsed_excess)

    n_resolvable = t1                                        # counted AND single-read assignable
    n_counted_unassignable = chi_H - n_resolvable           # = t2 + collapse reps
    assert n_counted_unassignable == t2 + n_nonsingleton_groups, (fam, "n_cu")

    true_copy_lower_bound = chi_H + collapsed_excess         # = n_loci (reference-resolved)
    assert true_copy_lower_bound == n_loci, (fam, "tclb", true_copy_lower_bound, n_loci)
    assert n_resolvable <= chi_H <= true_copy_lower_bound, (fam, "chain")

    return {
        "family": fam, "gene": s["gene"], "n_loci": n_loci, "chi_H": chi_H,
        "n_resolvable": n_resolvable, "n_counted_unassignable": n_counted_unassignable,
        "collapsed_excess": collapsed_excess, "true_copy_lower_bound": true_copy_lower_bound,
        "cls": p["cls"], "chi_H_readlevel": p["K"],
        "readlevel_disagrees": int(p["K"] != chi_H),
        "_n_nonsingleton_groups": n_nonsingleton_groups,
        "_t1": t1, "_t2": t2, "_t3": t3,
    }


def verify_hapvectors(fams, expected):
    """Reconstruct actual per-copy hap-vectors for `fams` and confirm, on real sequence, that
    chi_H == #distinct hap-vectors AND those hap-vectors are pairwise-conflicting (Lemma C2).
    Reuses sun_identifiability's copyonly asm20 path. Returns list of per-family check dicts."""
    try:
        sys.path.insert(0, HERE)
        import pysam
        import psv_graph_genomewide as pg
        import sun_identifiability as sun_mod
    except Exception as ex:  # pragma: no cover - environment without BAM/tools
        return [{"skipped": f"deps/data unavailable: {ex}"}]
    if not (os.path.exists(pg.GENOME) and os.path.exists(pg.FAM_TSV)):
        return [{"skipped": f"missing {pg.GENOME} or {pg.FAM_TSV}"}]

    fam_regions = collections.defaultdict(list)
    with open(pg.FAM_TSV) as fh:
        next(fh)
        for line in fh:
            fi, lid, c, s, e, nr = line.rstrip("\n").split("\t")
            fam_regions[int(fi)].append((lid, c, int(s), int(e), int(nr)))

    genome = pysam.FastaFile(pg.GENOME)
    clen = dict(zip(genome.references, genome.lengths))
    out = []
    for fam in fams:
        regions = pg.dedup_copies(fam_regions[fam])
        res = sun_mod.copy_psvs(regions, genome, clen)
        if res is None:
            out.append({"family": fam, "skipped": "backbone<200bp"})
            continue
        names, regions, psvs, _ = res
        hap_rows = {nm: tuple(hap[nm] for _, hap in psvs) for nm in names}
        groups = collections.defaultdict(list)
        for nm in names:
            groups[hap_rows[nm]].append(nm)
        distinct = list(groups)
        pairwise_conflict = all(
            any(a[i] != b[i] for i in range(len(a)))
            for a, b in itertools.combinations(distinct, 2))
        n_distinct = len(distinct)
        exp = expected.get(fam)
        out.append({
            "family": fam,  # gene label lives in the catalog rows, omitted from this seq-check
            "n_copies": len(names), "n_psv": len(psvs),
            "n_distinct_hapvectors": n_distinct,
            "pairwise_conflicting": pairwise_conflict,
            "group_sizes": sorted((len(v) for v in groups.values()), reverse=True),
            "chi_H_catalog": exp,
            "chi_H_equals_distinct_pairwise_conflicting": bool(
                pairwise_conflict and exp is not None and n_distinct == exp),
        })
    genome.close()
    return out


def main():
    pg, sun = load()
    rows = [family_row(fam, pg, sun) for fam in sorted(sun)]

    # --- THIRD LEG (depth): purely ADDITIVE. Annotate each row with depth_cn / regime /
    #     depth_recovers_collapsed_excess. The existing per-family values are untouched. ---
    depth_cn_map = load_depth_cn()
    for r in rows:
        dcn = depth_cn_map.get(r["family"])
        r["depth_cn"] = dcn
        # regime: reference-RESOLVED (chi_H<=n_loci; reference pre-splits distinct loci) vs
        # reference-COLLAPSED (n_loci<chi_H; one locus hides several hap-vectors -> depth-only).
        r["regime"] = "reference_collapsed" if r["n_loci"] < r["chi_H"] else "reference_resolved"
        # depth recovers the Tier-3 identical copies chi(H) collapsed away (need depth/DNA).
        r["depth_recovers_collapsed_excess"] = int(
            r["collapsed_excess"] > 0 and dcn is not None and dcn > r["chi_H"])

    with open(OUT_TSV, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=COLS, delimiter="\t", extrasaction="ignore")
        w.writeheader()
        for r in rows:
            w.writerow(r)

    agg = collections.Counter()
    for r in rows:
        for k in ("n_loci", "chi_H", "n_resolvable", "n_counted_unassignable",
                  "collapsed_excess", "true_copy_lower_bound", "_t1", "_t2", "_t3",
                  "_n_nonsingleton_groups"):
            agg[k] += r[k]

    fam_unassignable = [r for r in rows if r["chi_H"] > r["n_resolvable"]]      # advisor's point
    fam_collapsed = [r for r in rows if r["collapsed_excess"] > 0]              # need depth/DNA
    disagree = [r for r in rows if r["readlevel_disagrees"]]

    # --- depth-leg aggregates (ADDITIVE) ---
    regime_counts = collections.Counter(r["regime"] for r in rows)
    n_depth = sum(1 for r in rows if r["depth_cn"] is not None)
    depth_recovers = [r for r in fam_collapsed if r["depth_recovers_collapsed_excess"]]

    # named examples
    ex_beyond = max((r for r in rows if r["_t2"] > 0), key=lambda r: r["chi_H"], default=None) \
        or max(fam_unassignable, key=lambda r: r["chi_H"])
    ex_collapsed = max(fam_collapsed, key=lambda r: r["collapsed_excess"])

    def _ex(r):
        return {k: r[k] for k in ("family", "gene", "n_loci", "chi_H", "n_resolvable",
                                  "n_counted_unassignable", "collapsed_excess",
                                  "true_copy_lower_bound", "cls")}

    expected_chi = {r["family"]: r["chi_H"] for r in rows}
    verification = verify_hapvectors(VERIFY_FAMS, expected_chi) if VERIFY_FAMS else []
    verify_all_pass = bool(verification) and all(
        v.get("chi_H_equals_distinct_pairwise_conflicting") for v in verification
        if "skipped" not in v)

    summary = {
        "n_families": len(rows),
        "totals": {
            "n_loci": agg["n_loci"],
            "chi_H": agg["chi_H"],
            "n_resolvable_Tier1": agg["n_resolvable"],
            "n_counted_unassignable": agg["n_counted_unassignable"],
            "collapsed_excess": agg["collapsed_excess"],
            "true_copy_lower_bound_eq_n_loci": agg["true_copy_lower_bound"],
            "tier1": agg["_t1"], "tier2": agg["_t2"], "tier3": agg["_t3"],
            "nonsingleton_group_color_reps": agg["_n_nonsingleton_groups"],
        },
        "chi_H_decomposition": (
            f"chi_H({agg['chi_H']}) = Tier1({agg['_t1']}) + Tier2({agg['_t2']}) + "
            f"collapsed_group_reps({agg['_n_nonsingleton_groups']})"),
        "bound_chain": (
            f"n_resolvable({agg['n_resolvable']}) <= chi_H({agg['chi_H']}) <= "
            f"n_loci=true_copy_lower_bound({agg['n_loci']}); upper gap = "
            f"{agg['collapsed_excess']} Tier-3 identical copies chi(H) misses (need depth/DNA)"),
        "families_carry_unassignable_counted_copies": len(fam_unassignable),
        "families_need_depth_or_dna": len(fam_collapsed),
        "readlevel_divergence": {
            "n_families_psvgraphK_ne_copyonlyK": len(disagree),
            "families": [{"family": r["family"], "gene": r["gene"],
                          "chi_H_copyonly": r["chi_H"], "K_readlevel": r["chi_H_readlevel"],
                          "direction": ("read_splits_more" if r["chi_H_readlevel"] > r["chi_H"]
                                        else "consensus_splits_more")} for r in disagree],
        },
        "examples": {
            "chi_H_counts_beyond_assignable_core": _ex(ex_beyond),
            "collapsed_chi_H_undercounts": _ex(ex_collapsed),
        },
        "verification": {
            "families": VERIFY_FAMS,
            "all_pass": verify_all_pass,
            "note": "chi_H == #distinct pairwise-conflicting hap-vectors (reconstructed from "
                    "sequence); confirms K=chi(H) so no from-scratch recompute needed",
            "detail": verification,
        },
        "depth_leg": {
            "note": "THIRD independent lower-bound probe (ADDITIVE; shipped columns unchanged). "
                    "depth_cn = E_fam / lambda_within (family aggregate reads / within-family "
                    "median per-copy read depth; RNA self-normalized parCN, Sudmant 2010). "
                    "Same E_fam/lambda_within as rna_copy_number_depth.py (depth_cn_within), "
                    "reusing psv_graph_genomewide.dedup_copies. HONEST: depth counts EXPRESSED "
                    "copies => a NOISY LOWER BOUND on genomic CN (silent copies invisible => "
                    "biased low; per-copy expression skew = two-sided noise), NOT an unbiased count.",
            "n_families_with_depth_cn": n_depth,
            "regime_counts": dict(regime_counts),
            "n_collapsed_excess_families": len(fam_collapsed),
            "n_depth_recovers_collapsed_excess": len(depth_recovers),
            "collapsed_family_depth": [
                {"family": r["family"], "gene": r["gene"], "chi_H": r["chi_H"],
                 "collapsed_excess": r["collapsed_excess"], "depth_cn": r["depth_cn"],
                 "depth_recovers_collapsed_excess": r["depth_recovers_collapsed_excess"]}
                for r in sorted(fam_collapsed, key=lambda x: -x["collapsed_excess"])],
        },
        "deterministic": True,
    }

    json.dump({"summary": summary, "families":
               [{k: r[k] for k in COLS} for r in rows]},
              open(OUT_JSON, "w"), indent=1)

    print(f"wrote {OUT_TSV} and {OUT_JSON}  ({len(rows)} families)")
    print(summary["bound_chain"])
    print(summary["chi_H_decomposition"])
    print(f"families carrying UNASSIGNABLE-but-counted copies (chi_H>n_resolvable): "
          f"{len(fam_unassignable)}")
    print(f"families needing DEPTH/DNA (collapsed_excess>0): {len(fam_collapsed)}")
    print(f"read-level psv_graph K disagrees with copy-consensus chi_H on {len(disagree)}/154")
    print(f"DEPTH leg: depth_cn computed for {n_depth}/{len(rows)} families; regimes={dict(regime_counts)}; "
          f"depth_cn>chi_H on {len(depth_recovers)}/{len(fam_collapsed)} collapsed_excess>0 families")
    print(f"hap-vector verification {VERIFY_FAMS}: all_pass={verify_all_pass}")
    for v in verification:
        print("  ", v)
    ex = summary["examples"]
    print("EXAMPLE chi_H>core :", ex["chi_H_counts_beyond_assignable_core"])
    print("EXAMPLE collapse   :", ex["collapsed_chi_H_undercounts"])


if __name__ == "__main__":
    main()
