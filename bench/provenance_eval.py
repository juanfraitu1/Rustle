#!/usr/bin/env python3
"""Prereg Addendum X: provenance-typed scores for de novo catalogs run with `RUSTLE_ER_EDGE_DUMP`.

usage: provenance_eval.py (--truth truth.tsv | --guided clusters.tsv --contigs c1,c2,...) --expr expr.tsv
                          NAME=DUMP_PREFIX:CATALOG_PREFIX ...
DUMP_PREFIX.nodes.tsv / .edges.tsv are the first E_r call (TX: exon-sum substrate); DUMP_PREFIX.call2.edges.tsv, when the
run had `RUSTLE_ER_UNION_GENOMIC_SPAN`, is the genomic-span call (SPAN); DUMP_PREFIX.genomic_span_rejected_edges.tsv
lists span-only pairs restricted mode rejected. CATALOG_PREFIX.copies.tsv is the emitted catalog.

1. Per-type edge precision. Each rep takes its best-overlapping truth locus (most overlapping bp); an edge whose two reps
   both have one is TRUE iff those loci share a truth family (human: `family`; gorilla: guided cluster), else FALSE;
   other edges are unscored. Types: TX-only, BOTH, SPAN-only admitted, SPAN-only rejected.
2. Presence per truth locus: R = some rep's span lies >= 50% inside the locus; T = no R rep but u >= 3 (primary MAPQ >= 1
   reads with aligned bases in the locus, from --expr); N = u < 3. With the fraction of each category that an emitted
   copy overlaps (locus in an emitted family).
"""
import argparse
import collections
import csv
import os

GATE = 3


def ov(a0, a1, b0, b1):
    return max(0, min(a1, b1) - max(a0, b0))


def read_tsv(path):
    return list(csv.DictReader(open(path), delimiter="\t"))


def load_truth(a):
    """[(name, family, chrom, start0, end)]"""
    if a.truth:
        return [(r["name"], r["family"], r["chrom"], int(r["start0"]), int(r["end"])) for r in read_tsv(a.truth)]
    cs = set(a.contigs.split(","))
    return [(f"{r['cluster_id']}:{r['chrom']}:{r['start']}", r["cluster_id"], r["chrom"], int(r["start"]) - 1, int(r["end"]))
            for r in read_tsv(a.guided) if r["chrom"] in cs]


def best_locus(truth_by_chrom, chrom, s, e):
    hits = [(ov(s, e, t[3], t[4]), -k, k) for k, t in truth_by_chrom.get(chrom, []) if ov(s, e, t[3], t[4]) > 0]
    return max(hits)[2] if hits else None


def pairs(path):
    return {(min(int(r["rep_i"]), int(r["rep_j"])), max(int(r["rep_i"]), int(r["rep_j"]))) for r in read_tsv(path)}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--truth")
    ap.add_argument("--guided")
    ap.add_argument("--contigs")
    ap.add_argument("--expr", required=True)
    ap.add_argument("arms", nargs="+")
    a = ap.parse_args()
    assert bool(a.truth) != bool(a.guided), "give exactly one of --truth / --guided"
    truth = load_truth(a)
    by_chrom = collections.defaultdict(list)
    for k, t in enumerate(truth):
        by_chrom[t[2]].append((k, t))
    expr = {(r["chrom"], int(r["start"]), int(r["end"])): int(r["u"]) for r in read_tsv(a.expr)}
    u = [expr.get((t[2], t[3], t[4])) for t in truth]
    rows0 = {(r["chrom"], int(r["start"]), int(r["end"])): int(r["m0"]) for r in read_tsv(a.expr) if "m0" in r}
    m0 = [rows0.get((t[2], t[3], t[4]), 0) for t in truth] if rows0 else None
    assert all(x is not None for x in u), f"{sum(x is None for x in u)} truth loci missing from --expr"
    print(f"truth: {len(truth)} loci in {len({t[1] for t in truth})} families; u>={GATE}: {sum(x >= GATE for x in u)}")

    for spec in a.arms:
        name, rest = spec.split("=", 1)
        dump, cat = rest.split(":", 1)
        nodes = read_tsv(f"{dump}.nodes.tsv")
        assert [int(r["idx"]) for r in nodes] == list(range(len(nodes)))
        rep_locus = [best_locus(by_chrom, r["chrom"], int(r["start"]), int(r["end"])) for r in nodes]
        tx = pairs(f"{dump}.edges.tsv")
        span_path, rej_path = f"{dump}.call2.edges.tsv", f"{dump}.genomic_span_rejected_edges.tsv"
        span = pairs(span_path) if os.path.exists(span_path) else set()
        rejected = pairs(rej_path) if os.path.exists(rej_path) else set()
        if span:  # the second call must be the genomic-span substrate over the same reps
            n2 = read_tsv(f"{dump}.call2.nodes.tsv")
            assert [r["node_key"] for r in n2] == [r["node_key"] for r in nodes], "call2 reps differ from call1 reps"
        assert rejected <= span - tx, "rejected pairs must be span-only"
        types = {"TX-only": tx - span, "BOTH": tx & span, "SPAN-only admitted": span - tx - rejected,
                 "SPAN-only rejected": rejected}
        print(f"\n## {name}: {len(nodes)} reps, TX {len(tx)} edges, SPAN {len(span)} edges")
        print(f"{'type':20s} {'edges':>6s} {'scored':>6s} {'TRUE':>5s} {'FALSE':>5s} {'same-locus':>10s} {'precision':>9s}")
        for tname, es in types.items():
            sc = [(rep_locus[i], rep_locus[j]) for i, j in es if rep_locus[i] is not None and rep_locus[j] is not None]
            tp = sum(truth[p][1] == truth[q][1] for p, q in sc)
            same = sum(p == q for p, q in sc)
            prec = f"{tp / len(sc):.3f}" if sc else "nan"
            print(f"{tname:20s} {len(es):6d} {len(sc):6d} {tp:5d} {len(sc) - tp:5d} {same:10d} {prec:>9s}")

        copies = read_tsv(f"{cat}.copies.tsv")
        cop_chrom = collections.defaultdict(list)
        for c in copies:
            cop_chrom[c["chrom"]].append((int(c["start"]), int(c["end"])))
        rep_chrom = collections.defaultdict(list)
        for r in nodes:
            rep_chrom[r["chrom"]].append((int(r["start"]), int(r["end"])))
        cat_of = collections.Counter()
        in_fam = collections.Counter()
        per_family = collections.defaultdict(collections.Counter)
        notable = []
        for k, (tn, fam, c, s, e) in enumerate(truth):
            own = any(ov(s, e, x, y) >= 0.5 * (y - x) for x, y in rep_chrom.get(c, []) if y > x)
            cat_k = "R" if own else ("T" if u[k] >= GATE else "N")
            fam_k = any(ov(s, e, x, y) > 0 for x, y in cop_chrom.get(c, []))
            cat_of[cat_k] += 1
            in_fam[cat_k] += fam_k
            if a.truth and (cat_k != "R" or not fam_k):
                notable.append(f"{tn}:{cat_k}{'+fam' if fam_k else '-fam'}")
            if a.truth:
                per_family[fam][cat_k] += 1
                per_family[fam][cat_k + "_fam"] += fam_k
        print("presence: " + "  ".join(f"{x} {cat_of[x]} (in family {in_fam[x]}/{cat_of[x]})" for x in "RTN"))
        if notable:
            print("  not R, or R outside every family: " + " ".join(notable))
        if m0 is not None:  # disclosure, not pre-registered: T/N with MAPQ-0 reads also counted
            t0 = collections.Counter()
            for k, (tn, fam, c, s, e) in enumerate(truth):
                if not any(ov(s, e, x, y) >= 0.5 * (y - x) for x, y in rep_chrom.get(c, []) if y > x):
                    t0["T" if u[k] + m0[k] >= GATE else "N"] += 1
            print(f"  (disclosure, MAPQ-0 reads counted) T {t0['T']}  N {t0['N']}")
        for fam in sorted(per_family):
            pf = per_family[fam]
            print(f"  {fam:8s} " + "  ".join(f"{x} {pf[x]} (in family {pf[x + '_fam']})" for x in "RTN"))


if __name__ == "__main__":
    main()
