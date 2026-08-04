#!/usr/bin/env python3
"""Prove-removability: does the READS-ONLY data catalog (gw_family_catalog) group the same
loci into the same families as the ORACLE-fed catalog (family_define/family_rna_refine)?

Two partitions of read-assembled DN loci:
  ORACLE = bench/family_rna_refine.tsv   (family_id | ... | member_dn | chrom | start | end)
  DATA   = GGO_gwcat.copies.tsv          (family_id | copy_idx | tid | chrom | start | end | ...)

Both label loci with the same DN_<chrom>_<pos>_<nexon> scheme (member_dn == tid), so the PRIMARY
join is EXACT locus identity; an overlap join is reported as a robustness cross-check.

On the shared loci we report co-membership agreement (do the two catalogs group them the same):
  pair_precision = co-grouped in both / co-grouped in DATA   (data introduces no spurious merges vs oracle)
  pair_recall    = co-grouped in both / co-grouped in ORACLE (data recovers the oracle's groupings)
  ARI            = adjusted Rand index
  family purity  = fraction of oracle families whose shared members land in ONE data family (and reverse)
Coverage (how many loci each catalog contributes, how many shared) is reported honestly first.
"""
import csv, math
from collections import defaultdict, Counter
from itertools import combinations

ORACLE = "/mnt/c/Users/jfris/Desktop/Rustle/bench/family_rna_refine.tsv"
DATA   = "/mnt/linuxdisk/home/juanfraitu/winloci_data/GGO_gwcat.copies.tsv"

def load(path, fam_prefix, dn_col):
    rows = []  # dict: dn, chrom, start, end, fam
    with open(path) as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            try:
                rows.append(dict(dn=row[dn_col], chrom=row["chrom"],
                                 start=int(row["start"]), end=int(row["end"]),
                                 fam=fam_prefix + row["family_id"]))
            except (KeyError, ValueError):
                continue
    return rows

def pair_metrics(labels_o, labels_d):
    n = len(labels_o)
    both = o_only = d_only = 0
    for i, j in combinations(range(n), 2):
        so = labels_o[i] == labels_o[j]
        sd = labels_d[i] == labels_d[j]
        if so and sd: both += 1
        elif so: o_only += 1
        elif sd: d_only += 1
    prec = both/(both+d_only) if (both+d_only) else float("nan")
    rec  = both/(both+o_only) if (both+o_only) else float("nan")
    f1 = 2*prec*rec/(prec+rec) if (prec and rec and prec==prec and rec==rec) else float("nan")
    # ARI
    cont = Counter(zip(labels_o, labels_d)); a = Counter(labels_o); b = Counter(labels_d)
    C2 = lambda x: x*(x-1)//2
    sij = sum(C2(v) for v in cont.values()); sa = sum(C2(v) for v in a.values()); sb = sum(C2(v) for v in b.values())
    exp = sa*sb/C2(n) if C2(n) else 0; mx = 0.5*(sa+sb)
    ari = (sij-exp)/(mx-exp) if (mx-exp) else float("nan")
    return both, o_only, d_only, prec, rec, f1, ari

def family_purity(pairs_o, pairs_d):
    o2d = defaultdict(list)
    for lo, ld in zip(pairs_o, pairs_d):
        o2d[lo].append(ld)
    pure = tot = 0
    for ds in o2d.values():
        if len(ds) < 2: continue
        tot += 1; pure += (len(set(ds)) == 1)
    return pure, tot

def report(lo, ld, coverage_note):
    print(coverage_note)
    if len(lo) < 2:
        print("  (too few shared loci)"); print(); return
    both, oo, dd, prec, rec, f1, ari = pair_metrics(lo, ld)
    pure, tot = family_purity(lo, ld)
    rpure, rtot = family_purity(ld, lo)
    print(f"  co-membership pairs: both={both}  oracle-only={oo}  data-only={dd}")
    print(f"  pair_precision (data merges match oracle) = {prec:.3f}")
    print(f"  pair_recall    (oracle groupings recovered)= {rec:.3f}")
    print(f"  pair_F1={f1:.3f}   adjusted_Rand={ari:.3f}")
    if tot:  print(f"  oracle families (>=2 shared members) landing in ONE data family : {pure}/{tot} = {pure/tot:.1%}")
    if rtot: print(f"  data   families (>=2 shared members) landing in ONE oracle family: {rpure}/{rtot} = {rpure/rtot:.1%}")
    print()

def recip(a, b):
    if a["chrom"] != b["chrom"]: return 0.0
    s, e = max(a["start"], b["start"]), min(a["end"], b["end"])
    if e <= s: return 0.0
    inter = e - s; la = a["end"]-a["start"]; lb = b["end"]-b["start"]
    return min(inter/la, inter/lb) if la > 0 and lb > 0 else 0.0

def main():
    oracle = load(ORACLE, "O:", "member_dn")
    data   = load(DATA,   "D:", "tid")
    print(f"ORACLE members: {len(oracle)}  families: {len(set(r['fam'] for r in oracle))}")
    print(f"DATA   copies : {len(data)}  families: {len(set(r['fam'] for r in data))}")
    print()

    # ---- PRIMARY: exact DN-id join ----
    o_by_dn = {r["dn"]: r for r in oracle}
    d_by_dn = {r["dn"]: r for r in data}
    shared = sorted(set(o_by_dn) & set(d_by_dn))
    lo = [o_by_dn[dn]["fam"] for dn in shared]
    ld = [d_by_dn[dn]["fam"] for dn in shared]
    note = (f"### EXACT DN-id join (headline)\n"
            f"  shared loci: {len(shared)}   coverage: {len(shared)/len(oracle):.1%} of oracle members, "
            f"{len(shared)/len(data):.1%} of data copies")
    report(lo, ld, note)

    # ---- CROSS-CHECK: reciprocal-overlap join (in case DN ids drifted) ----
    for thr in (0.5, 0.1):
        by_chrom = defaultdict(list)
        for di, d in enumerate(data): by_chrom[d["chrom"]].append(di)
        cands = []
        for oi, o in enumerate(oracle):
            for di in by_chrom.get(o["chrom"], []):
                ro = recip(o, data[di])
                if ro >= thr: cands.append((ro, oi, di))
        cands.sort(reverse=True)
        ou = du = set(); ou, du = set(), set(); pairs = []
        for ro, oi, di in cands:
            if oi in ou or di in du: continue
            ou.add(oi); du.add(di); pairs.append((oi, di))
        loo = [oracle[oi]["fam"] for oi, di in pairs]
        ldd = [data[di]["fam"] for oi, di in pairs]
        note = (f"### reciprocal-overlap join >= {thr} (cross-check)\n"
                f"  matched loci: {len(pairs)}   coverage: {len(pairs)/len(oracle):.1%} oracle, {len(pairs)/len(data):.1%} data")
        report(loo, ldd, note)

if __name__ == "__main__":
    main()
