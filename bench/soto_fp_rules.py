#!/usr/bin/env python3
"""Search for false-merge rules against SOTO-LABELLED pairs (ledger §6by).

This is the first FP-rule search in the project with an **external** positive/negative labelling
of co-membership assertions. Every earlier search used a proxy (CDS concordance) or NPIP.

⛔⛔ **THE CIRCULAR RULE, NAMED SO IT IS NEVER PROPOSED.** Soto's families are **SD98 — >=98%
identity by construction**. Therefore "raise the E_r identity floor to 0.90" scores brilliantly
against this truth set **for the same reason the truth set exists**, not because it is a good rule.
Any rule whose discriminating variable is edge IDENTITY is unadjudicable here. ⟹ **every rule below
is evaluated WITHIN the [0.90,1.00) band**, where Soto's own criterion is satisfied on both sides
and identity can no longer do the separating. That is the same within-stratum discipline §6bn used.

⚠ Precision is UNDERSTATED throughout (Soto is CAT-bounded), so a rule that trims "FP" may be
trimming real copies CAT missed. Rules are ranked by their TP cost, never by precision alone.

Usage: soto_fp_rules.py <bed> <arm_dir> [--floor 0.0001|0.50]
"""
import csv
import sys
import os
import math
import collections
import itertools


def wilson(k, n, z=1.96):
    if n == 0:
        return (0.0, 1.0)
    p = k / n
    d = 1 + z * z / n
    c = (p + z * z / (2 * n)) / d
    h = z * math.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / d
    return (max(0.0, c - h), min(1.0, c + h))


def main():
    bed, arm = sys.argv[1], sys.argv[2]
    floor = float(sys.argv[sys.argv.index("--floor") + 1]) if "--floor" in sys.argv else 0.0001

    genes, fam = {}, {}
    for line in open(bed):
        f = line.rstrip("\n").split("\t")
        genes[f[3]] = (f[0], int(f[1]), int(f[2]))
        fam[f[3]] = f[3].split("|")[1]

    copies = list(csv.DictReader(open(os.path.join(arm, "cat.copies.tsv")), delimiter="\t"))
    a = {}
    for i, r in enumerate(copies):
        c, s, e = r["chrom"], int(r["start"]), int(r["end"])
        best = None
        for n, (gc, gs, ge) in genes.items():
            if gc != c:
                continue
            ov = min(ge, e) - max(gs, s)
            if ov <= 0:
                continue
            fr = ov / max(1, ge - gs)
            if fr >= floor and (best is None or fr > best[0]):
                best = (fr, n)
        if best:
            a[i] = best[1]

    byfam = collections.defaultdict(list)
    for g, f in fam.items():
        byfam[f].append(g)
    truth = {tuple(sorted(p)) for gs in byfam.values() for p in itertools.combinations(sorted(gs), 2)}

    # our asserted pairs, with the family they came from
    pred = collections.defaultdict(set)
    for i, r in enumerate(copies):
        if i in a:
            pred[r["family_id"]].add(a[i])
    pair_fam, ours = {}, set()
    for fid, gs in pred.items():
        for p in itertools.combinations(sorted(gs), 2):
            k = tuple(sorted(p))
            ours.add(k)
            pair_fam[k] = fid
    famsize = {f: len(v) for f, v in pred.items()}

    # direct-edge table: pair -> best edge record
    g_of = {}
    for i, r in enumerate(copies):
        if i in a:
            g_of[(r["chrom"], int(r["start"]), int(r["end"]))] = a[i]
    edge = {}
    ef = os.path.join(arm, "dump", "e.edges.tsv")
    for r in csv.DictReader(open(ef), delimiter="\t"):
        ki = (r["chrom_i"], int(r["start_i"]), int(r["end_i"]))
        kj = (r["chrom_j"], int(r["start_j"]), int(r["end_j"]))
        gi, gj = g_of.get(ki), g_of.get(kj)
        if not gi or not gj or gi == gj:
            continue
        k = tuple(sorted((gi, gj)))
        rec = dict(identity=float(r["identity"]), coverage=float(r["coverage"]),
                   cov_longer=float(r.get("cov_longer") or "nan"),
                   same_chrom=r["chrom_i"] == r["chrom_j"])
        if k not in edge or rec["identity"] > edge[k]["identity"]:
            edge[k] = rec

    def report(label, sub):
        tp = sum(1 for p in sub if p in truth)
        n = len(sub)
        lo, hi = wilson(tp, n) if n else (float("nan"),) * 2
        print(f"  {label:44s} n={n:5d}  TP={tp:5d}  FP={n-tp:5d}  "
              f"P={tp/n if n else float('nan'):.4f} [{lo:.4f},{hi:.4f}]")
        return tp, n

    print(f"arm={os.path.basename(arm.rstrip('/'))}  gene-overlap floor={floor}")
    print(f"asserted pairs={len(ours)}  with a DIRECT edge={sum(1 for p in ours if p in edge)}")

    print("\n=== R1  DIRECT EDGE vs TRANSITIVE CLOSURE (identity plays no part) ===")
    direct = [p for p in ours if p in edge]
    trans = [p for p in ours if p not in edge]
    tpd, nd = report("pair HAS a direct E_r edge", direct)
    tpt, nt = report("pair is TRANSITIVE ONLY (no direct edge)", trans)
    if nd and nt:
        print(f"  ⟹ dropping transitive-only pairs: precision {(tpd+tpt)/(nd+nt):.4f} -> {tpd/nd:.4f}"
              f"   TP cost {tpt} of {tpd+tpt} = {tpt/(tpd+tpt):.1%}")

    print("\n=== WITHIN [0.90,1.00) ONLY — where Soto's own criterion cannot do the separating ===")
    hi_band = [p for p in direct if 0.90 <= edge[p]["identity"] < 1.001]
    report("baseline: direct edge, identity >= 0.90", hi_band)
    print()
    for name, fn in [
        ("R2  coverage(shorter) >= 0.80", lambda p: edge[p]["coverage"] >= 0.80),
        ("R3  coverage(shorter) >= 0.90", lambda p: edge[p]["coverage"] >= 0.90),
        ("R4  cov_longer >= 0.30", lambda p: edge[p]["cov_longer"] >= 0.30),
        ("R5  cov_longer >= 0.50", lambda p: edge[p]["cov_longer"] >= 0.50),
        ("R6  same chromosome", lambda p: edge[p]["same_chrom"]),
        ("R7  cross chromosome", lambda p: not edge[p]["same_chrom"]),
        ("R8  predicted family size <= 5", lambda p: famsize.get(pair_fam[p], 99) <= 5),
        ("R9  predicted family size <= 10", lambda p: famsize.get(pair_fam[p], 99) <= 10),
        ("R10 family size <= 5 AND cov_longer >= 0.30",
         lambda p: famsize.get(pair_fam[p], 99) <= 5 and edge[p]["cov_longer"] >= 0.30),
    ]:
        keep = [p for p in hi_band if fn(p)]
        tpk, nk = report(name, keep)
        base_tp = sum(1 for p in hi_band if p in truth)
        if nk:
            print(f"  {'':44s}   -> keeps {tpk}/{base_tp} = {tpk/base_tp:.1%} of the band's TP")


if __name__ == "__main__":
    main()
