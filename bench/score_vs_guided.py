#!/usr/bin/env python3
"""Prereg Addenda E/K: score catalogs against guided `gw_units_v3` families.

usage: score_vs_guided.py --contigs include|exclude:<c1,c2,...> [--expr expr.tsv] name=copies.tsv ...
Guided loci are restricted to the contig set (and, with --expr, to loci with u >= 3), then to clusters with >= 2 loci.
A catalog name ending in `+expr` keeps only copies whose own interval (atom key `atom:<idx>` matched by coordinates)
has u >= 3 and then only families with >= 2 such copies (the K1 hybrid). Reports, per catalog:
  any-family:   R_G = guided co-clustered pairs whose loci share ANY family / all guided pairs; P_G = those / all locus
                pairs sharing any family (Addendum E scorer);
  best-overlap: each locus takes the single family of its best-overlapping copy (robustness).
"""
import argparse
import bisect
import collections
import csv

GATE = 3
R = "/mnt/linuxdisk/home/juanfraitu"
ap = argparse.ArgumentParser()
ap.add_argument("--contigs", required=True)
ap.add_argument("--expr")
ap.add_argument("catalogs", nargs="+")
a = ap.parse_args()
mode, cs = a.contigs.split(":", 1)
cset = set(cs.split(","))
keep_chrom = (lambda c: c in cset) if mode == "include" else (lambda c: c not in cset)

expr = {}
if a.expr:
    for r in csv.DictReader(open(a.expr), delimiter="\t"):
        expr[(r["chrom"], int(r["start"]), int(r["end"]))] = int(r["u"])

g = [(r["cluster_id"], r["chrom"], int(r["start"]) - 1, int(r["end"]))
     for r in csv.DictReader(open(f"{R}/mcl_ann/gw_units_v3.clusters.tsv"), delimiter="\t") if keep_chrom(r["chrom"])]
if a.expr:
    g = [x for x in g if expr.get((x[1], x[2], x[3]), 0) >= GATE]
cnt = collections.Counter(x[0] for x in g)
guided = [x for x in g if cnt[x[0]] >= 2]
by_cluster = collections.defaultdict(list)
for li, x in enumerate(guided):
    by_cluster[x[0]].append(li)
nG = sum(len(v) * (len(v) - 1) // 2 for v in by_cluster.values())
print(f"guided: {len(guided)} loci in {len(by_cluster)} clusters (>=2 loci), {nG} co-clustered pairs"
      + (" [expressed loci only]" if a.expr else ""))

for spec in a.catalogs:
    name, path = spec.split("=", 1)
    copies = [r for r in csv.DictReader(open(path), delimiter="\t") if keep_chrom(r["chrom"])]
    if name.endswith("+expr"):
        copies = [r for r in copies if expr.get((r["chrom"], int(r["start"]), int(r["end"])), 0) >= GATE]
        fc = collections.Counter(r["family_id"] for r in copies)
        copies = [r for r in copies if fc[r["family_id"]] >= 2]
    bych = collections.defaultdict(list)
    for r in copies:
        bych[r["chrom"]].append((int(r["start"]), int(r["end"]), r["family_id"]))
    for c in bych:
        bych[c].sort()
    st = {c: [x[0] for x in v] for c, v in bych.items()}
    mx = {c: max(e - s for s, e, _ in v) for c, v in bych.items()}

    def hits(c, s, e):
        if c not in bych:
            return []
        lo = bisect.bisect_left(st[c], s - mx[c])
        hi = bisect.bisect_left(st[c], e)
        return [(min(e, x[1]) - max(s, x[0]), x[2]) for x in bych[c][lo:hi] if x[1] > s]

    H = [hits(c, s, e) for _, c, s, e in guided]
    fs = [{f for _, f in h} for h in H]
    best = [max(h)[1] if h else None for h in H]
    covered = sum(1 for h in H if h)
    rec_any = sum(1 for lis in by_cluster.values() for p in range(len(lis)) for q in range(p + 1, len(lis))
                  if fs[lis[p]] & fs[lis[q]])
    fam_loci = collections.defaultdict(set)
    for li, f in enumerate(fs):
        for x in f:
            fam_loci[x].add(li)
    pred = set()
    for ls in fam_loci.values():
        ls = sorted(ls)
        for p in range(len(ls)):
            for q in range(p + 1, len(ls)):
                pred.add((ls[p], ls[q]))
    tp_any = sum(1 for p, q in pred if guided[p][0] == guided[q][0])
    bf = collections.defaultdict(list)
    for li, b in enumerate(best):
        if b is not None:
            bf[b].append(li)
    tp_b = fp_b = 0
    for ls in bf.values():
        for p in range(len(ls)):
            for q in range(p + 1, len(ls)):
                if guided[ls[p]][0] == guided[ls[q]][0]:
                    tp_b += 1
                else:
                    fp_b += 1
    nf = len({r["family_id"] for r in copies})
    print(f"{name:26s} families={nf:5d} copies={len(copies):5d} loci_with_family_copy={covered:5d}  "
          f"any: R_G={rec_any / nG:.4f} P_G={tp_any / len(pred) if pred else float('nan'):.4f}  "
          f"best: R_G={tp_b / nG:.4f} P_G={tp_b / (tp_b + fp_b) if tp_b + fp_b else float('nan'):.4f}")
