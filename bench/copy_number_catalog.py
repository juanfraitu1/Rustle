#!/usr/bin/env python3
"""
copy_number_catalog.py -- Per-family COPY-NUMBER catalog (O1 enrichment).

Advisor (Canzar) refinement of the E_r RNA-family reframe: the read-CONFLICT
signal returns to O1 as the family's COPY NUMBER -- how many copies a family has,
INCLUDING copies that cannot be assigned. This does NOT reverse the reframe
(RNA family membership = transcript-homology E_r component, family_definition_formal.md);
it ADDS the conflict-derived multiplicity chi(H) as a family-level O1 property, and
stratifies the copies by the SUN 3-tier resolvability ladder.

THE COUNT vs ASSIGN distinction (theorem-level, THEORY.md Lemma 1 / Thm 2 / Thm 4):
  * COUNTING copies = MCC = chi(H) (chromatic number of the copy-conflict graph;
    Lemma 1). Needs ONLY the conflict structure; well-defined EVEN WHEN reads
    cannot be assigned. On the de-tied / significance-gated (copy-consensus) graph and in
    the full-span regime, chi(H) is a LOWER BOUND on the true copy number (the RAW
    allele-disagreement read graph over-counts to ~3xK from error edges, so it is NOT a
    lower bound there -- the gate/copy-consensus is a precondition).
  * ASSIGNING is strictly harder: (i) SINGLE-READ assignability of a read (|N(r)|=1) needs
    a SUN = Theorem 4(ii) / the Bridge (the per-read certificate, the per-copy Tier-1
    property); (ii) family-wide UNIQUE-COVER recovery needs Strong Separation (Theorem 2,
    COVERAGE-dependent). Tier-1 / SUN is per-copy, NOT the coverage-dependent Strong-Sep.
  => a family is COUNTED as chi(H) copies even when only the Tier-1 subset is single-read
     ASSIGNABLE. Bound chain (reference-RESOLVED here): n_resolvable <= chi(H) <= n_loci <=
     true; general form max(n_loci, chi(H)) <= true (reference-COLLAPSED SDA/Vollger swaps
     the two middle terms, n_loci <= chi(H)).

3-tier copy hierarchy (per copy, from bench/sun_identifiability.*):
  Tier1 RESOLVABLE            SUN-identifiable / Strong-Sep: counted by chi(H) AND
                             single-read assignable -- the nice core (O2).
  Tier2 DISTINGUISHABLE-BUT-  unique hap-vector, NO SUN: counted by chi(H), NOT
        UNASSIGNABLE           single-read assignable (needs >=2-PSV co-observation)
                             -- the advisor's "important" unassignable copies.
  Tier3 IDENTICAL / collapsed share a hap-vector, no distinguishing PSV: chi(H)
                             COLLAPSES them to one color -> chi(H) UNDER-counts ->
                             countable only by read-DEPTH / DNA (K=0 frontier, parCN/O4).

Sources (all keyed by the integer family id, consistent across the three):
  bench/psv_graph_genomewide.json   per-family read-level K, group_sizes, cls
  bench/sun_identifiability.json    per-family copyonly_K, tier1/2/3 counts
  bench/sun_identifiability.tsv     per-copy tier + group_size (copy-consensus grouping)

chi(H) is taken on the COPY-CONSENSUS conflict graph (copyonly_K), which is
tier-consistent (singleton groups == Tier1+Tier2) and is the correct O1 COUNTING
object -- it counts copies distinguishable in principle from the assembled copy
sequences, independent of read sampling depth. The read-level psv_graph K is the
ASSIGNMENT-realized count (what THIS read set resolves); it is carried as a
cross-check column and diverges on 11/154 families (documented below).

Run from repo root:  python bench/copy_number_catalog.py
Writes bench/copy_number_catalog.tsv and bench/copy_number_catalog.json (override
with COPY_NUM_OUT / COPY_NUM_JSON).
"""
import json, csv, os
from collections import defaultdict, Counter

PG_JSON  = "bench/psv_graph_genomewide.json"
SUN_JSON = "bench/sun_identifiability.json"
SUN_TSV  = "bench/sun_identifiability.tsv"
OUT_TSV  = os.environ.get("COPY_NUM_OUT",  "bench/copy_number_catalog.tsv")
OUT_JSON = os.environ.get("COPY_NUM_JSON", "bench/copy_number_catalog.json")


def load():
    pg   = {f["family"]: f for f in json.load(open(PG_JSON))["families"]}
    sunj = {f["family"]: f for f in json.load(open(SUN_JSON))["families"]}
    suntsv = defaultdict(list)
    with open(SUN_TSV) as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            suntsv[int(row["family"])].append(row)
    return pg, sunj, suntsv


def family_metrics(fam, pg, sunj, suntsv):
    """Return the copy-number metric row for one family, with the exact identities
    asserted so the algebra is machine-checked (fails loudly if a source drifts)."""
    sj, f, copies = sunj[fam], pg[fam], suntsv[fam]

    n_loci = sj["n_copies"]                 # co-located reference copies fed to the graph
    t1, t2, t3 = sj["tier1"], sj["tier2"], sj["tier3"]
    assert t1 + t2 + t3 == n_loci, (fam, "tiers", t1, t2, t3, n_loci)

    chi_H = sj["copyonly_K"]                # chi(H) on the copy-consensus conflict graph
    chi_H_readlevel = f["K"]               # read-level psv_graph realization (cross-check)

    # group structure from per-copy group_size (copy-consensus / copyonly grouping)
    gsz = Counter(int(c["group_size"]) for c in copies)     # size -> #copies with that size
    n_singleton_groups = gsz.get(1, 0)
    n_nonsingleton_groups = sum(cnt // g for g, cnt in gsz.items() if g > 1)
    collapsed_excess = sum(cnt - cnt // g for g, cnt in gsz.items() if g > 1)  # Sum(size-1)

    # ---- exact identities (Section: metric algebra) ----
    assert n_singleton_groups == t1 + t2, (fam, "singleton", n_singleton_groups, t1, t2)
    assert n_singleton_groups + n_nonsingleton_groups == chi_H, (fam, "K")
    assert t3 == n_loci - n_singleton_groups, (fam, "t3")
    assert collapsed_excess == t3 - n_nonsingleton_groups, (fam, "excess")

    n_resolvable = t1                                       # Tier1 Strong-Sep, assignable
    n_counted_unassignable = chi_H - n_resolvable          # = t2 + n_nonsingleton_groups
    assert n_counted_unassignable == t2 + n_nonsingleton_groups, (fam, "unassignable")
    true_copy_lower_bound = chi_H + collapsed_excess        # = n_loci (reference-resolved)
    assert true_copy_lower_bound == n_loci, (fam, "tclb", true_copy_lower_bound, n_loci)

    return {
        "family": fam,
        "gene": sj["gene"],
        "cls": f["cls"],
        "n_loci": n_loci,
        "chi_H": chi_H,
        "chi_H_readlevel": chi_H_readlevel,
        "readlevel_disagrees": int(chi_H_readlevel != chi_H),
        "n_resolvable": n_resolvable,               # Tier1  (counted AND assignable)
        "n_counted_unassignable": n_counted_unassignable,  # Tier2 + collapsed-group reps
        "tier1": t1, "tier2": t2, "tier3": t3,
        "n_nonsingleton_groups": n_nonsingleton_groups,
        "collapsed_excess": collapsed_excess,       # identical copies chi(H) misses
        "true_copy_lower_bound": true_copy_lower_bound,
    }


def main():
    pg, sunj, suntsv = load()
    rows = [family_metrics(fam, pg, sunj, suntsv) for fam in sorted(pg)]

    cols = ["family", "gene", "cls", "n_loci", "chi_H", "chi_H_readlevel",
            "readlevel_disagrees", "n_resolvable", "n_counted_unassignable",
            "tier1", "tier2", "tier3", "n_nonsingleton_groups",
            "collapsed_excess", "true_copy_lower_bound"]
    with open(OUT_TSV, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=cols, delimiter="\t")
        w.writeheader()
        for r in rows:
            w.writerow(r)

    agg = defaultdict(int)
    for r in rows:
        for k in ("n_loci", "chi_H", "n_resolvable", "n_counted_unassignable",
                  "tier1", "tier2", "tier3", "collapsed_excess", "true_copy_lower_bound"):
            agg[k] += r[k]
    disagree = [r for r in rows if r["readlevel_disagrees"]]

    summary = {
        "n_families": len(rows),
        "totals": dict(agg),
        "bound_chain_sums": {
            "n_resolvable_Tier1": agg["n_resolvable"],
            "chi_H": agg["chi_H"],
            "true_copy_lower_bound_eq_n_loci": agg["true_copy_lower_bound"],
            "statement": "Sum n_resolvable(338) <= Sum chi(H)(361) <= Sum n_loci=true_lb(412)",
        },
        "readlevel_divergence": {
            "n_families_psvgraphK_ne_copyonlyK": len(disagree),
            "families": [{"family": r["family"], "gene": r["gene"],
                          "copyonly_chi_H": r["chi_H"], "readlevel_K": r["chi_H_readlevel"],
                          "direction": ("read_splits_more" if r["chi_H_readlevel"] > r["chi_H"]
                                        else "consensus_splits_more")}
                         for r in disagree],
        },
    }
    json.dump({"summary": summary, "families": rows}, open(OUT_JSON, "w"), indent=1)

    print(f"wrote {OUT_TSV} and {OUT_JSON}  ({len(rows)} families)")
    print(f"bound chain (genome-wide sums): "
          f"n_resolvable(T1)={agg['n_resolvable']} <= chi(H)={agg['chi_H']} "
          f"<= n_loci=true_copy_lower_bound={agg['n_loci']}")
    print(f"  chi(H) = Tier1({agg['tier1']}) + Tier2({agg['tier2']}) + "
          f"collapsed-group-reps({agg['chi_H']-agg['tier1']-agg['tier2']})")
    print(f"  collapsed_excess (identical copies chi(H) misses) = {agg['collapsed_excess']} "
          f"= Tier3({agg['tier3']}) - nonsingleton-groups")
    print(f"  read-level psv_graph K disagrees on {len(disagree)}/154 families "
          f"(both directions; count-vs-assign at the graph level)")


if __name__ == "__main__":
    main()
