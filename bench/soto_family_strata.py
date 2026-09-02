#!/usr/bin/env python3
"""Are some families easy and others hard — and does the right rule differ by kind? (ledger §6ca)

⚠⚠ **THE TRAP THIS SCRIPT IS BUILT AROUND.** Finding categories in the errors and then fitting a
rule per category is the garden of forking paths, and §6bz already declined to build
`size in [6,15] AND lambda >= 2` for exactly that reason. So **every category here comes from SOTO's
OWN metadata** (`soto_famCN_S1C.tsv`), fixed before this project existed and computed from their data,
never from our error rate:

    No. Family Members          family size
    Family MAD                  THEIR within-family copy-number dispersion — a difficulty label
    No. Protein Coding / No. Unprocessed Pseudogene   composition
    Median famCN                copy number

**ALL strata are reported, never the best one**, and any per-category rule is stated as a hypothesis
with the split that would test it — not as a result.

⚠ Precision is UNDERSTATED (Soto is CAT-bounded); ratios between strata are the robust quantity.

Usage: soto_family_strata.py <bed> <arm_dir> <famcn_tsv>
"""
import csv
import sys
import os
import math
import collections
import itertools


def wilson(k, n, z=1.96):
    if n == 0:
        return (float("nan"), float("nan"))
    p = k / n
    d = 1 + z * z / n
    c = (p + z * z / (2 * n)) / d
    h = z * math.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / d
    return (max(0.0, c - h), min(1.0, c + h))


def num(x):
    try:
        return float(x)
    except (TypeError, ValueError):
        return float("nan")


def main():
    bed, arm, famcn = sys.argv[1], sys.argv[2], sys.argv[3]

    genes, fam = {}, {}
    for line in open(bed):
        f = line.rstrip("\n").split("\t")
        genes[f[3]] = (f[0], int(f[1]), int(f[2]))
        fam[f[3]] = f[3].split("|")[1]
    byfam = collections.defaultdict(list)
    for g, f in fam.items():
        byfam[f].append(g)
    truth = {tuple(sorted(p)) for gs in byfam.values() for p in itertools.combinations(sorted(gs), 2)}
    fam_of = {g: f for g, f in fam.items()}

    # ---- Soto's own per-family metadata (categories fixed outside this project)
    meta = {}
    for r in csv.DictReader(open(famcn), delimiter="\t"):
        fid = r.get("Family ID")
        if fid and fid not in meta:
            meta[fid] = dict(
                size=num(r.get("No. Family Members")),
                mad=num(r.get("Family MAD")),
                coding=num(r.get("No. Protein Coding")),
                pseudo=num(r.get("No. Unprocessed Pseudogene")),
                famcn=num(r.get("Median famCN")),
            )

    # ---- our prediction
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
            if best is None or fr > best[0]:
                best = (fr, n)
        if best:
            a[i] = best[1]
    detected = set(a.values())

    key = lambda r: (r["chrom"], int(r["start"]), int(r["end"]))
    adj_all = collections.defaultdict(set)
    for r in csv.DictReader(open(os.path.join(arm, "dump", "e.edges.tsv")), delimiter="\t"):
        ki = (r["chrom_i"], int(r["start_i"]), int(r["end_i"]))
        kj = (r["chrom_j"], int(r["start_j"]), int(r["end_j"]))
        adj_all[ki].add(kj)
        adj_all[kj].add(ki)

    pred = collections.defaultdict(list)
    for i, r in enumerate(copies):
        pred[r["family_id"]].append((key(r), a.get(i)))

    # asserted gene pairs, with the distance backing each
    asserted = {}
    for fid, mem in pred.items():
        S = {k for k, _ in mem}
        adj = {k: (adj_all[k] & S) for k in S}
        for (ku, gu), (kv, gv) in itertools.combinations(mem, 2):
            if gu is None or gv is None or gu == gv:
                continue
            seen = {ku: 0}
            q = collections.deque([ku])
            while q:
                x = q.popleft()
                for y in adj[x]:
                    if y not in seen:
                        seen[y] = seen[x] + 1
                        q.append(y)
            d = seen.get(kv)
            if d is None:
                continue
            p = tuple(sorted((gu, gv)))
            if p not in asserted or d < asserted[p]:
                asserted[p] = d

    # ---- per-SOTO-family scorecard
    per = {}
    for f, gs in byfam.items():
        if len(gs) < 2:
            continue
        tp_all = {tuple(sorted(p)) for p in itertools.combinations(sorted(gs), 2)}
        det = [g for g in gs if g in detected]
        got = {p for p in tp_all if p in asserted}
        got_d1 = {p for p in got if asserted[p] == 1}
        # pairs we asserted that TOUCH this family but leave it (contamination)
        foreign = 0
        for p, d in asserted.items():
            f0, f1 = fam_of[p[0]], fam_of[p[1]]
            if (f0 == f) != (f1 == f):
                foreign += 1
        per[f] = dict(
            members=len(gs), detected=len(det),
            truth_pairs=len(tp_all), recovered=len(got), recovered_d1=len(got_d1),
            foreign=foreign, **meta.get(f, {}),
        )

    def strat(label, keyfn, buckets):
        print(f"\n=== stratified by {label} (SOTO's own metadata) ===")
        print(f"  {'stratum':>12s} {'fams':>5s} {'members':>8s} {'detect':>7s} "
              f"{'pair recall':>12s} {'|d=1 share':>11s} {'foreign/fam':>12s}")
        for lab, test in buckets:
            fs = [v for v in per.values() if test(keyfn(v))]
            if not fs:
                print(f"  {lab:>12s}     0")
                continue
            m = sum(v["members"] for v in fs)
            dd = sum(v["detected"] for v in fs)
            tp = sum(v["truth_pairs"] for v in fs)
            rc = sum(v["recovered"] for v in fs)
            r1 = sum(v["recovered_d1"] for v in fs)
            fo = sum(v["foreign"] for v in fs)
            print(f"  {lab:>12s} {len(fs):>5d} {m:>8d} {dd/m:>7.3f} "
                  f"{rc/tp if tp else float('nan'):>12.3f} {r1/rc if rc else float('nan'):>11.3f} "
                  f"{fo/len(fs):>12.1f}")

    print(f"scored {len(per)} multi-member Soto families "
          f"({sum(v['members'] for v in per.values())} members, "
          f"{sum(v['truth_pairs'] for v in per.values())} truth pairs)")

    strat("family SIZE (No. Family Members)", lambda v: v.get("size", float("nan")),
          [("2-3", lambda x: x <= 3), ("4-5", lambda x: 4 <= x <= 5),
           ("6-9", lambda x: 6 <= x <= 9), ("10+", lambda x: x >= 10)])
    strat("Family MAD (their CN dispersion)", lambda v: v.get("mad", float("nan")),
          [("MAD=0", lambda x: x == 0), ("0<MAD<=0.3", lambda x: 0 < x <= 0.3),
           ("MAD>0.3", lambda x: x > 0.3), ("N/A", lambda x: x != x)])
    strat("composition (protein-coding share)",
          lambda v: (v.get("coding", 0) or 0) / max(1e-9, (v.get("coding", 0) or 0) + (v.get("pseudo", 0) or 0)),
          [("all pseudo", lambda x: x == 0), ("mixed", lambda x: 0 < x < 1),
           ("all coding", lambda x: x >= 1)])
    strat("Median famCN", lambda v: v.get("famcn", float("nan")),
          [("<=2", lambda x: x <= 2), ("2-4", lambda x: 2 < x <= 4),
           (">4", lambda x: x > 4), ("N/A", lambda x: x != x)])

    print("\n=== the per-family scatter (are there EASY and HARD families at all?) ===")
    rec = sorted((v["recovered"] / v["truth_pairs"], f) for f, v in per.items() if v["truth_pairs"])
    perfect = sum(1 for r, _ in rec if r >= 0.999)
    zero = sum(1 for r, _ in rec if r == 0)
    print(f"  families with pair recall 1.00: {perfect}/{len(rec)}   "
          f"with 0.00: {zero}/{len(rec)}   middle: {len(rec)-perfect-zero}")
    print(f"  hardest: {[f for _, f in rec[:6]]}")
    print(f"  easiest: {[f for _, f in rec[-6:]]}")


if __name__ == "__main__":
    main()
