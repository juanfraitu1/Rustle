#!/usr/bin/env python3
"""Adjudicate `RUSTLE_ER_COVERAGE_LONGER_FLOOR` against the decision rule of ledger §6bu.

Prices the copies the clause deletes. Compares an OFF and an ON catalog produced by the same
binary on the same substrate, and characterises the LOST copies against the RETAINED ones.

⚠ **D3 is deliberately labelled weak.** A short single-exon node aligning into a long partner has
low `cov_longer` BY CONSTRUCTION, so "the losses are stubs" partly restates the selector rather
than testing it. Reported as consistency, never as validation.
⚠ **D4 is confounded by length** — annotation presence correlates with node size, which correlates
with `cov_longer` — so it is reported WITHIN length strata, per §6bu.
⭐ **D5 is the load-bearing one**: does the clause trim families, or delete them wholesale?

Usage:
  adjudicate_covlonger.py <off_dir> <on_dir> [--truth npip31.regions] [--gff GGO_genomic.gff]
"""
import csv
import sys
import os
import math
import collections


def wilson(k, n, z=1.96):
    if n == 0:
        return (0.0, 1.0)
    p = k / n
    d = 1 + z * z / n
    c = (p + z * z / (2 * n)) / d
    h = z * math.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / d
    return (max(0.0, c - h), min(1.0, c + h))


def load_copies(d):
    return list(csv.DictReader(open(os.path.join(d, "cat.copies.tsv")), delimiter="\t"))


def load_nodes(d):
    """node_key -> row, from the edge-dump node table (has reads/exons/stub/len)."""
    p = os.path.join(d, "dump", "e.nodes.tsv")
    out = {}
    if os.path.exists(p):
        for r in csv.DictReader(open(p), delimiter="\t"):
            out[(r["chrom"], int(r["start"]), int(r["end"]))] = r
    return out


def key(r):
    return (r["chrom"], int(r["start"]), int(r["end"]))


def load_truth(p):
    out = []
    for line in open(p):
        line = line.strip()
        if not line:
            continue
        c, rng = line.split(":")
        s, e = rng.split("-")
        out.append((c, int(s), int(e)))
    return out


def in_truth(k, truth):
    c, s, e = k
    return any(c == tc and s < te and ts < e for tc, ts, te in truth)


def main():
    off_d, on_d = sys.argv[1], sys.argv[2]
    truth_p = sys.argv[sys.argv.index("--truth") + 1] if "--truth" in sys.argv else None
    gff_p = sys.argv[sys.argv.index("--gff") + 1] if "--gff" in sys.argv else None

    off, on = load_copies(off_d), load_copies(on_d)
    nodes = load_nodes(off_d)
    off_k = {key(r): r for r in off}
    on_k = {key(r): r for r in on}
    lost = [k for k in off_k if k not in on_k]
    kept = [k for k in off_k if k in on_k]
    print(f"OFF copies {len(off)}  ON copies {len(on)}  LOST {len(lost)}  KEPT {len(kept)}")

    # ---------- D5: families trimmed, or deleted wholesale?
    foff = collections.defaultdict(list)
    for r in off:
        foff[r["family_id"]].append(key(r))
    fon = collections.defaultdict(list)
    for r in on:
        fon[r["family_id"]].append(key(r))
    on_members = set(on_k)
    deleted_fams = trimmed_fams = intact_fams = 0
    deleted_copies = trimmed_copies = 0
    for fid, mem in foff.items():
        surv = [m for m in mem if m in on_members]
        if not surv:
            deleted_fams += 1
            deleted_copies += len(mem)
        elif len(surv) < len(mem):
            trimmed_fams += 1
            trimmed_copies += len(mem) - len(surv)
        else:
            intact_fams += 1
    print(f"\n=== D5 (REQUIRED): families trimmed vs deleted ===")
    print(f"  OFF families {len(foff)} -> ON families {len(fon)}")
    print(f"  families fully DELETED : {deleted_fams}  (carrying {deleted_copies} copies)")
    print(f"  families TRIMMED       : {trimmed_fams}  (losing  {trimmed_copies} copies)")
    print(f"  families INTACT        : {intact_fams}")
    if lost:
        print(f"  ⟹ share of lost copies that died with their whole family: "
              f"{deleted_copies}/{len(lost)} = {deleted_copies/len(lost):.4f}")

    # ---------- D1: truth
    if truth_p:
        truth = load_truth(truth_p)
        lt = sum(1 for k in lost if in_truth(k, truth))
        kt = sum(1 for k in kept if in_truth(k, truth))
        print(f"\n=== D1: overlap with the truth set ({len(truth)} regions) ===")
        print(f"  truth-overlapping copies KEPT: {kt}   LOST: {lt}")
        if lt:
            for k in lost:
                if in_truth(k, truth):
                    print(f"     ⚠ LOST TRUE-SET COPY {k[0]}:{k[1]}-{k[2]}")

    # ---------- D3 (weak, by construction): node properties
    def props(ks):
        reads, exons, lens, stubs, n = [], [], [], 0, 0
        for k in ks:
            r = nodes.get(k)
            if not r:
                continue
            n += 1
            reads.append(int(r["n_reads"]))
            exons.append(int(r["n_exon"]))
            lens.append(int(r["exon_sum_len"]))
            stubs += 1 if r["stub"] == "true" else 0
        med = lambda v: sorted(v)[len(v) // 2] if v else float("nan")
        return dict(n=n, med_reads=med(reads), med_exons=med(exons), med_len=med(lens),
                    single_exon=sum(1 for e in exons if e == 1), stub=stubs,
                    two_read=sum(1 for x in reads if x <= 2))
    pl, pk = props(lost), props(kept)
    print(f"\n=== D3 (⚠ WEAK — correlated with the selector; consistency only) ===")
    print(f"  {'':22s} {'LOST':>10s} {'KEPT':>10s}")
    for f, lab in [("n", "matched nodes"), ("med_reads", "median reads"),
                   ("med_exons", "median exons"), ("med_len", "median exon-sum bp")]:
        print(f"  {lab:22s} {pl[f]:>10} {pk[f]:>10}")
    for f, lab in [("single_exon", "single-exon"), ("stub", "stub=true"), ("two_read", "<=2 reads")]:
        a = pl[f] / pl["n"] if pl["n"] else float("nan")
        b = pk[f] / pk["n"] if pk["n"] else float("nan")
        print(f"  {lab:22s} {pl[f]:>5} ={a:6.3f} {pk[f]:>5} ={b:6.3f}")

    # ---------- D4: annotation, WITHIN length strata
    if gff_p:
        sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
        src = open(os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                "er_identity_band_external.py")).read()
        ns = {"__name__": "_x"}
        exec(src.split("def main()")[0], ns)
        idx = ns["load_gff"](gff_p)
        ann = ns["annotate"]
        lens_all = sorted(int(nodes[k]["exon_sum_len"]) for k in off_k if k in nodes)
        import bisect
        q = [lens_all[int(len(lens_all) * f)] for f in (0.25, 0.5, 0.75)] if lens_all else [0, 0, 0]

        def stratum(k):
            L = int(nodes[k]["exon_sum_len"])
            return bisect.bisect_left(q, L)
        tally = collections.defaultdict(lambda: collections.Counter())
        for grp, ks in (("LOST", lost), ("KEPT", kept)):
            for k in ks:
                if k not in nodes:
                    continue
                c = tally[(grp, stratum(k))]
                c["n"] += 1
                c["ann"] += 1 if ann(idx, k[0], k[1], k[2]) else 0
        print(f"\n=== D4: reciprocal RefSeq match, WITHIN length quartile (⚠ never pooled) ===")
        print(f"  quartile cuts (exon-sum bp): {q}")
        print(f"  {'stratum':>8s} {'LOST n':>8s} {'ann':>6s} {'rate':>7s}   {'KEPT n':>8s} {'ann':>6s} {'rate':>7s}")
        for s in range(4):
            a, b = tally[("LOST", s)], tally[("KEPT", s)]
            ra = a["ann"] / a["n"] if a["n"] else float("nan")
            rb = b["ann"] / b["n"] if b["n"] else float("nan")
            print(f"  {'Q'+str(s+1):>8s} {a['n']:>8} {a['ann']:>6} {ra:>7.3f}   "
                  f"{b['n']:>8} {b['ann']:>6} {rb:>7.3f}")


if __name__ == "__main__":
    main()
