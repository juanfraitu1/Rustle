#!/usr/bin/env python3
"""Soto 2025 families as a truth table with exons: S1C (Gene ID -> Family ID) joined to CAT v4 (CHM13 v1.0 BED12) by
Gene ID (column 19), never by name; optionally lifted to CHM13 v2.0 through the T2T v1.0 -> v1.1 chain.

    soto_truth.py S1C.tsv CAT_V4.bed OUT.tsv [--chain CHAIN] [--report R.json]

OUT columns: label gene_id chrom start end strand exons — GFF 1-based closed; `exons` = the union of the exon blocks of
every transcript of the gene, `s-e` comma-joined; start/end = their hull. `Unassigned_*` labels are not families and
are dropped. With --chain an exon with an end in a chain gap is dropped and counted, and a gene left without exons is
dropped and counted.
"""
import argparse
import collections
import csv
import json
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from v1_to_v2 import Lifter, read_chain  # noqa: E402


def merge(blocks):
    """0-based half-open blocks -> their union, sorted."""
    out = []
    for s, e in sorted(blocks):
        if out and s <= out[-1][1]:
            out[-1][1] = max(out[-1][1], e)
        else:
            out.append([s, e])
    return [tuple(b) for b in out]


def cat_genes(lines):
    """CAT v4 BED12 -> {gene_id: (chrom, strand, [0-based exon blocks, unioned])}."""
    acc = {}
    for line in lines:
        f = line.rstrip("\n").split("\t")
        if len(f) < 19:
            continue
        chrom, s, strand, gid = f[0], int(f[1]), f[5], f[18]
        sizes = [int(x) for x in f[10].rstrip(",").split(",")]
        offs = [int(x) for x in f[11].rstrip(",").split(",")]
        g = acc.setdefault(gid, (chrom, strand, []))
        assert g[0] == chrom and g[1] == strand, f"{gid}: transcripts on different chromosomes or strands"
        g[2].extend((s + o, s + o + z) for o, z in zip(offs, sizes))
    return {gid: (c, st, merge(b)) for gid, (c, st, b) in acc.items()}


def s1c_rows(path):
    """-> [(Family ID, Gene ID)] for the families proper (no `Unassigned_*`)."""
    with open(path) as fh:
        rows = list(csv.reader(fh, delimiter="\t"))
    h = rows[0]
    gi, fi = h.index("Gene ID"), h.index("Family ID")
    return [(r[fi], r[gi]) for r in rows[1:] if len(r) > fi and r[fi] and not r[fi].startswith("Unassigned_")]


def build(rows, genes, lifter=None):
    out, c = [], collections.Counter()
    for label, gid in sorted(set(rows)):
        if gid not in genes:
            c["genes_not_in_cat"] += 1
            continue
        chrom, strand, ex = genes[gid]
        if lifter is not None:
            lifted = [lifter.lift(chrom, s, e) for s, e in ex]
            c["exons_unliftable"] += sum(r is None for r in lifted)
            ex = merge((r[0], r[1]) for r in lifted if r is not None)
            if not ex:
                c["genes_unliftable"] += 1
                continue
        out.append((label, gid, chrom, ex[0][0] + 1, ex[-1][1], strand, ",".join(f"{s + 1}-{e}" for s, e in ex)))
        c["rows"] += 1
    c["labels"] = len({r[0] for r in out})
    return out, dict(sorted(c.items()))


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("s1c"); ap.add_argument("cat"); ap.add_argument("out")
    ap.add_argument("--chain"); ap.add_argument("--report")
    a = ap.parse_args(argv)
    with open(a.cat) as fh:
        genes = cat_genes(fh)
    lifter = Lifter(read_chain(open(a.chain))) if a.chain else None
    out, counts = build(s1c_rows(a.s1c), genes, lifter)
    with open(a.out, "w") as fh:
        fh.write("label\tgene_id\tchrom\tstart\tend\tstrand\texons\n")
        fh.writelines("\t".join(map(str, r)) + "\n" for r in out)
    rep = {"inputs": vars(a), "counts": counts}
    if a.report:
        with open(a.report, "w") as fh:
            json.dump(rep, fh, indent=1)
    print(json.dumps(rep, indent=1))


if __name__ == "__main__":
    sys.exit(main())
