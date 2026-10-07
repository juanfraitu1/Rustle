#!/usr/bin/env python3
"""Level C of Amendment 2 of docs/PREREG_locus_units_2026-10-06.md: the connected components of the pre-MCL homology graph (mcl_families --dump-graph) as clusters, in the registered cluster-file format,
so the unchanged scorer (bench/ideal_expression/score.py, fam_score.py) can be run with the components in place of the MCL clusters.

    hier_clusters.py --loci PREFIX.fam.loci.gff3 --graph GRAPH.tsv --out OUTPREFIX
Writes OUTPREFIX.fam.clusters.tsv (components of >= 2 locus keys) and links PREFIX.fam.loci.gff3, .fam.loci.tsv and .families.gtf to OUTPREFIX when --prefix PREFIX is given.
Nodes of the graph that are not loci of the arm (no gene row with that span in the loci file) are ignored; a locus span shared by two loci makes the arm INVALID (the scorer's own rule).
"""
import argparse
import os
import sys

HEADER = ["cluster_id", "size", "density", "frac_in", "corroborated", "chrom", "start", "end"]


def parse_node(s):
    chrom, _, rest = s.rpartition(":")
    a, _, b = rest.partition("-")
    return chrom, int(a), int(b)


def loci_keys(gff3):
    keys = {}
    for ln in open(gff3):
        if ln.startswith("#"):
            continue
        f = ln.rstrip("\n").split("\t")
        if len(f) >= 9 and f[2] == "gene":
            keys.setdefault((f[0], int(f[3]), int(f[4])), []).append(f[8])
    dup = [k for k, v in keys.items() if len(v) > 1]
    if dup:
        sys.exit(f"INVALID: loci share a span: {dup[:3]}")
    return set(keys)


def components(loci_gff3, graph_tsv):
    keys = loci_keys(loci_gff3)
    parent = {}

    def find(x):
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x

    for ln in open(graph_tsv):
        f = ln.rstrip("\n").split("\t")
        if len(f) != 3:
            continue
        a, b = parse_node(f[0]), parse_node(f[1])
        if a not in keys or b not in keys:
            continue
        for x in (a, b):
            parent.setdefault(x, x)
        ra, rb = find(a), find(b)
        if ra != rb:
            parent[ra] = rb
    groups = {}
    for x in parent:
        groups.setdefault(find(x), []).append(x)
    rows = []
    for i, g in enumerate(sorted(groups.values(), key=lambda g: (-len(g), min(g)))):
        if len(g) < 2:
            continue
        for chrom, s, e in sorted(g):
            rows.append((f"CC{i}", chrom, s, e))
    return rows


def write_clusters(rows, path):
    size = {}
    for cid, *_ in rows:
        size[cid] = size.get(cid, 0) + 1
    with open(path, "w") as fh:
        fh.write("\t".join(HEADER) + "\n")
        for cid, chrom, s, e in rows:
            fh.write("\t".join([cid, str(size[cid]), "NA", "NA", "NA", chrom, str(s), str(e)]) + "\n")


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--loci", required=True)
    ap.add_argument("--graph", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--prefix")
    a = ap.parse_args()
    rows = components(a.loci, a.graph)
    write_clusters(rows, a.out + ".fam.clusters.tsv")
    if a.prefix:
        for suf in (".fam.loci.gff3", ".fam.loci.tsv", ".families.gtf"):
            dst = a.out + suf
            if os.path.lexists(dst):
                os.remove(dst)
            os.symlink(os.path.abspath(a.prefix + suf), dst)
    n = len({r[0] for r in rows})
    print(f"[hier_clusters] {len(rows)} loci in {n} components of >= 2 -> {a.out}.fam.clusters.tsv")


if __name__ == "__main__":
    main()
