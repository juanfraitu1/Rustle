#!/usr/bin/env python3
"""Prereg Addendum AG: the RNA-level ground truth (definition clause 6) and the goal metrics.

build: rna_truth.py build --clusters guided.clusters.tsv --expr expr.tsv --graph graph.tsv --contigs c1,c2 --out rna.clusters.tsv
  Within each guided cluster, the expressed loci (u >= 3) are joined iff connected by guided-graph edges whose nodes all
  overlap expressed loci of that cluster; RNA truth clusters are those connected components (>= 2 loci).
score: rna_truth.py score --clusters truth.clusters.tsv --expr expr.tsv --contigs c1,c2 [--expressed-only] name=copies.tsv ...
  Best-overlap assignment of every truth locus (in clusters >= 2) to a catalog family; pairwise sensitivity/precision
  and bipartite micro recall/precision/F (`guided_pipeline.pairwise` / `bipartite`).
"""
import argparse
import bisect
import collections
import csv
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import guided_pipeline as gp  # noqa: E402


def load_expr(path):
    return {(r["chrom"], int(r["start"]), int(r["end"])): int(r["u"]) for r in csv.DictReader(open(path), delimiter="\t")}


def cmd_build(a):
    C = set(a.contigs.split(","))
    expr = load_expr(a.expr)
    rows = [r for r in csv.DictReader(open(a.clusters), delimiter="\t") if r["chrom"] in C]
    by_cluster = collections.defaultdict(list)
    for r in rows:
        k = (r["chrom"], int(r["start"]) - 1, int(r["end"]))
        if expr.get(k, 0) >= 3 and k not in by_cluster[r["cluster_id"]]:
            by_cluster[r["cluster_id"]].append(k)
    # graph nodes indexed per contig
    nodes = {}
    idx = collections.defaultdict(list)

    def node_id(name):
        if name not in nodes:
            c, rng = name.rsplit(":", 1)
            s, e = rng.split("-")
            nodes[name] = (c, int(s) - 1, int(e))
        return nodes[name]
    edges = []
    for line in open(a.graph):
        f = line.rstrip("\n").split("\t")
        if len(f) < 3:
            continue
        u, v = node_id(f[0]), node_id(f[1])
        if u[0] in C and v[0] in C:
            edges.append((u, v))
    locus_index = collections.defaultdict(list)  # contig -> [(start, end, cluster, locus_i)]
    for cid, loci in by_cluster.items():
        for i, (c, s, e) in enumerate(loci):
            locus_index[c].append((s, e, cid, i))
    for c in locus_index:
        locus_index[c].sort()
    starts = {c: [x[0] for x in v] for c, v in locus_index.items()}
    maxlen = {c: max(x[1] - x[0] for x in v) for c, v in locus_index.items()}
    memo = {}

    def loci_of(n):
        if n not in memo:
            c, s, e = n
            out = []
            if c in locus_index:
                lo = bisect.bisect_left(starts[c], s - maxlen[c])
                hi = bisect.bisect_left(starts[c], e)
                out = [(cid, i) for x0, x1, cid, i in locus_index[c][lo:hi] if x1 > s and x0 < e]
            memo[n] = out
        return memo[n]
    parent = {(cid, i): (cid, i) for cid, loci in by_cluster.items() for i in range(len(loci))}

    def find(x):
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x

    def union(x, y):
        rx, ry = find(x), find(y)
        if rx != ry:
            parent[rx] = ry
    for n in set(nodes.values()):  # a node overlapping several expressed loci of one cluster joins them
        by = collections.defaultdict(list)
        for cid, i in loci_of(n):
            by[cid].append(i)
        for cid, ii in by.items():
            for i in ii[1:]:
                union((cid, ii[0]), (cid, i))
    n_used = 0
    for u, v in edges:
        lu, lv = loci_of(u), loci_of(v)
        if not lu or not lv:
            continue
        for cu, iu in lu:
            for cv, iv in lv:
                if cu == cv:
                    union((cu, iu), (cv, iv))
                    n_used += 1
    comps = collections.defaultdict(list)
    for cid, loci in by_cluster.items():
        for i in range(len(loci)):
            comps[find((cid, i))].append((cid, i))
    n_out = n_loci = 0
    with open(a.out, "w") as fh:
        fh.write("cluster_id\tchrom\tstart\tend\tguided_cluster\n")
        for k, (root, members) in enumerate(sorted(comps.items(), key=lambda kv: (kv[0][0], min(m[1] for m in kv[1])))):
            if len(members) < 2:
                continue
            n_out += 1
            for cid, i in sorted(members, key=lambda m: m[1]):
                c, s, e = by_cluster[cid][i]
                fh.write(f"RNA{k}\t{c}\t{s + 1}\t{e}\t{cid}\n")
                n_loci += 1
    n_guided = sum(1 for v in by_cluster.values() if len(v) >= 2)
    pairs_dna = sum(len(v) * (len(v) - 1) // 2 for v in by_cluster.values())
    pairs_rna = sum(len(m) * (len(m) - 1) // 2 for m in comps.values())
    print(f"guided clusters with >= 2 expressed loci: {n_guided} ({pairs_dna} expressed pairs); RNA-level truth: {n_out} clusters, "
          f"{n_loci} loci, {pairs_rna} pairs; graph edges used {n_used}")


def cmd_score(a):
    C = set(a.contigs.split(","))
    expr = load_expr(a.expr) if a.expr else None
    g = [(r["cluster_id"], r["chrom"], int(r["start"]) - 1, int(r["end"]))
         for r in csv.DictReader(open(a.clusters), delimiter="\t") if r["chrom"] in C]
    if expr is not None and a.expressed_only:
        g = [x for x in g if expr.get((x[1], x[2], x[3]), 0) >= 3]
    cnt = collections.Counter(x[0] for x in g)
    loci = [x for x in g if cnt[x[0]] >= 2]
    true = [x[0] for x in loci]
    print(f"truth: {len(loci)} loci in {len(set(true))} clusters, "
          f"{sum(v * (v - 1) // 2 for v in collections.Counter(true).values())} pairs")
    print(f"{'catalog':14s} {'pair_sens':>9s} {'pair_prec':>9s} {'bip_R':>6s} {'bip_P':>6s} {'bip_F':>6s}")
    for spec in a.catalogs:
        name, path = spec.split("=", 1)
        by = collections.defaultdict(list)
        for r in csv.DictReader(open(path), delimiter="\t"):
            if r["chrom"] in C:
                by[r["chrom"]].append((int(r["start"]), int(r["end"]), r["family_id"]))
        pred = []
        for i, x in enumerate(loci):
            h = [(min(x[3], e) - max(x[2], s), f) for s, e, f in by[x[1]] if s < x[3] and x[2] < e]
            pred.append(max(h)[1] if h else f"none:{i}")
        ps, pp = gp.pairwise(pred, true)
        br, bp = gp.bipartite(pred, true)
        f = 2 * br * bp / (br + bp) if br + bp else float("nan")
        print(f"{name:14s} {ps:9.3f} {pp:9.3f} {br:6.3f} {bp:6.3f} {f:6.3f}")


def main():
    ap = argparse.ArgumentParser()
    sub = ap.add_subparsers(dest="cmd", required=True)
    p = sub.add_parser("build")
    for k in ("--clusters", "--expr", "--graph", "--contigs", "--out"):
        p.add_argument(k, required=True)
    p = sub.add_parser("score")
    p.add_argument("--clusters", required=True)
    p.add_argument("--expr")
    p.add_argument("--contigs", required=True)
    p.add_argument("--expressed-only", action="store_true")
    p.add_argument("catalogs", nargs="+")
    a = ap.parse_args()
    {"build": cmd_build, "score": cmd_score}[a.cmd](a)


if __name__ == "__main__":
    main()


def scoped_truth(clusters, contigs, seeds_tsv, out):
    """Guided-mode scope: keep guided clusters (on the contigs, >= 2 loci) that contain >= 1 locus overlapping a seed
    gene; all loci of those clusters are kept (hidden members must still be found)."""
    C = set(contigs.split(","))
    seeds = collections.defaultdict(list)
    for line in open(seeds_tsv):
        if line.startswith("name\t"):
            continue
        n, c, s, e, st, sd, ex = line.rstrip("\n").split("\t")
        if sd == "1":
            seeds[c].append((int(s), int(e)))
    rows = [r for r in csv.DictReader(open(clusters), delimiter="\t") if r["chrom"] in C]
    cnt = collections.Counter(r["cluster_id"] for r in rows)
    has_seed = set()
    for r in rows:
        s, e = int(r["start"]) - 1, int(r["end"])
        if any(a < e and s < b for a, b in seeds[r["chrom"]]):
            has_seed.add(r["cluster_id"])
    with open(out, "w") as fh:
        fh.write("cluster_id\tchrom\tstart\tend\n")
        for r in rows:
            if cnt[r["cluster_id"]] >= 2 and r["cluster_id"] in has_seed:
                fh.write(f"{r['cluster_id']}\t{r['chrom']}\t{r['start']}\t{r['end']}\n")
    return len({r["cluster_id"] for r in rows if cnt[r["cluster_id"]] >= 2}), len(has_seed & {r["cluster_id"] for r in rows if cnt[r["cluster_id"]] >= 2})
