#!/usr/bin/env python3
"""Prereg Addendum AJ-2: does E1 repair the RefSeq truth (keep HGNC same-group pairs, drop cross-group pairs)?

usage: aj_hgnc.py --hgnc hgnc_complete_set.txt (--gff refseq.gff | --nodes nodes.tsv) --e0 E0.clusters.tsv --e1 E1.clusters.tsv
  A locus's symbols: --gff = Names of gene/pseudogene records whose 1-based coordinates equal the clusters.tsv row;
  --nodes = the node's record name (`<nodes>.names.tsv`). H = pairs whose symbols share >= 1 gene_group_id; X = pairs
  where both have a non-empty gene_group_id and share none. REPAIR iff H(E1) >= 0.90 H(E0) and X(E1) <= 0.50 X(E0).
"""
import argparse
import collections
import csv
import itertools


def load_hgnc(path):
    groups = {}
    for r in csv.DictReader(open(path), delimiter="\t"):
        groups[r["symbol"]] = {g for g in r.get("gene_group_id", "").split("|") if g}
    return groups


def locus_symbols(a):
    sym = collections.defaultdict(set)
    if a.gff:
        for line in open(a.gff):
            f = line.rstrip("\n").split("\t")
            if len(f) > 8 and f[2] in ("gene", "pseudogene"):
                for kv in f[8].split(";"):
                    if kv.startswith("Name="):
                        sym[(f[0], int(f[3]), int(f[4]))].add(kv[5:])
    else:
        names = {r["idx"]: r["name"] for r in csv.DictReader(open(a.nodes + ".names.tsv"), delimiter="\t")}
        for r in csv.DictReader(open(a.nodes), delimiter="\t"):
            sym[(r["chrom"], int(r["start"]) + 1, int(r["end"]))].add(names[r["idx"]])
    return sym


def count(clusters, sym, groups):
    by = collections.defaultdict(list)
    for r in csv.DictReader(open(clusters), delimiter="\t"):
        by[r["cluster_id"]].append((r["chrom"], int(r["start"]), int(r["end"])))
    h = x = n = joined = total = 0
    for loci in by.values():
        gs = []
        for k in loci:
            total += 1
            s = sym.get(k, set())
            joined += bool(s)
            g = set().union(*(groups.get(y, set()) for y in s)) if s else set()
            gs.append(g)
        for g1, g2 in itertools.combinations(gs, 2):
            n += 1
            if g1 & g2:
                h += 1
            elif g1 and g2:
                x += 1
    return h, x, n, joined, total


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--hgnc", required=True)
    ap.add_argument("--gff")
    ap.add_argument("--nodes")
    ap.add_argument("--e0", required=True)
    ap.add_argument("--e1", required=True)
    a = ap.parse_args()
    groups = load_hgnc(a.hgnc)
    sym = locus_symbols(a)
    res = {}
    for tag, path in (("E0", a.e0), ("E1", a.e1)):
        h, x, n, joined, total = count(path, sym, groups)
        res[tag] = (h, x)
        print(f"{tag}: pairs {n}  H (same HGNC group) {h}  X (both grouped, disjoint) {x}  "
              f"loci with a symbol {joined}/{total}")
    (h0, x0), (h1, x1) = res["E0"], res["E1"]
    keep, drop = h1 / max(h0, 1), x1 / max(x0, 1)
    verdict = "REPAIR" if keep >= 0.90 and drop <= 0.50 else "NOT REPAIR"
    print(f"H kept {keep:.3f} (>= 0.90), X kept {drop:.3f} (<= 0.50) -> {verdict}")


if __name__ == "__main__":
    main()
