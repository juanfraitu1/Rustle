#!/usr/bin/env python3
"""Post hoc (docs/SEED_POOL_REAL_READS_2026-10-07.md section 5): whose reads are the secondary alignments that carry a copy's annotated introns?

    paralog_support.py --copies copies.tsv --truth truth.gtf --bam reads.bam --support support.copies.tsv --family NPIP \
                       [--as-table molecules.tsv] --out paralog_support.json

For the copies expressed under Amendment A (ann_expressed == 1 in --support, as in secondary_support.py): k = min(2, the copy's longest annotated intron
chain); an alignment counts when it is on the copy's strand, overlaps the copy's exon union, is not supplementary and carries >= k annotated introns
inside the union's extent. Reported:
  share        of the SECONDARY alignments that count, how many belong to a molecule whose PRIMARY alignment lies in the territory of ANOTHER copy of the family
               (other_copy), in the territory of the same copy (same_copy), or anywhere else (elsewhere). Territory = terr_lo0 .. terr_hi of the copy table;
  good         per copy, how many of those secondary alignments have AS >= 0.98 x the molecule's genome-wide best AS (column 2 of --as-table);
  chains       per copy, the largest number of PRIMARY alignments with an identical junction tuple: the whole read ('whole') and the junctions inside the
               copy's extent ('inside'), each also over MAPQ > 0 only;
  e_sensitivity  E (copies with >= 2 primary alignments that count), with MAPQ > 0 required, and with a floor of 3 primaries.
It reads the BAM only: nothing is written but --out and the text on stdout.
"""
import argparse
import collections
import csv
import json
import os
import statistics
import sys

import pysam

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))
import copy_support as cs  # noqa: E402

GOOD_RATIO = 0.98


def overlaps(chrom, lo, hi, c):
    return chrom == c["chrom"] and lo < c["hi"] and hi > c["lo"]


def where(loc, cid, spans):
    """Where the primary alignment `loc` = (chrom, start, end) of a molecule lies relative to the territories of the family's copies."""
    if loc is None:
        return "elsewhere"
    if any(c2 != cid and overlaps(loc[0], loc[1], loc[2], s) for c2, s in spans.items()):
        return "other_copy"
    if overlaps(loc[0], loc[1], loc[2], spans[cid]):
        return "same_copy"
    return "elsewhere"


def collect(a):
    copies = [r for r in csv.DictReader(open(a.copies), delimiter="\t") if r["family"] == a.family]
    tx = cs.gtf_transcripts(a.truth, {r["isoform_gene"] for r in copies} | {r["cid"] for r in copies})
    expressed = {r["name"] for r in csv.DictReader(open(a.support), delimiter="\t") if r["ann_expressed"] == "1"}
    spans = {c["cid"]: dict(chrom=c["chrom"], lo=int(c["terr_lo0"]), hi=int(c["terr_hi"])) for c in copies}
    bam = pysam.AlignmentFile(a.bam)
    rows = {}
    for c in copies:
        name = c.get("refseq_name") or c.get("cat_name") or c["name"]
        if name not in expressed:
            continue
        texons = tx.get(c["cid"]) or tx.get(c["isoform_gene"], {})
        ann_all = set().union(*[set(cs.introns_of(ex)) for ex in texons.values()])
        k = min(2, max(len(cs.introns_of(ex)) for ex in texons.values()))
        union = cs.merge([e for ex in texons.values() for e in ex])
        lo, hi = union[0][0], union[-1][1]
        prim, sec = [], []
        for rd in bam.fetch(c["chrom"], lo, hi):
            if rd.is_unmapped or rd.is_supplementary or ("-" if rd.is_reverse else "+") != c["strand"]:
                continue
            blocks, juncs = cs.read_blocks_junctions(rd)
            if cs.inter(blocks, union) == 0:
                continue
            inside = tuple(j for j in juncs if j[0] >= lo and j[1] <= hi)
            if sum(1 for j in inside if j in ann_all) < k:
                continue
            rec = dict(name=rd.query_name, AS=(rd.get_tag("AS") if rd.has_tag("AS") else 0), mapq=rd.mapping_quality, whole=tuple(juncs), inside=inside)
            (sec if rd.is_secondary else prim).append(rec)
        rows[c["cid"]] = dict(name=name, k=k, prim=prim, sec=sec)
    return spans, rows, bam


def best_as(path, names):
    best = {}
    with open(path) as fh:
        for ln in fh:
            if ln[0] == "#":
                continue
            i = ln.find("\t")
            if ln[:i] in names:
                best[ln[:i]] = int(ln.split("\t")[1])
    return best


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--copies", required=True)
    ap.add_argument("--truth", required=True)
    ap.add_argument("--bam", required=True)
    ap.add_argument("--support", required=True, help="copy_support.py's PREFIX.copies.tsv (its ann_expressed column picks the copies)")
    ap.add_argument("--family", required=True)
    ap.add_argument("--as-table", help="the best-AS table (as_table); without it the 'good' counts are not reported")
    ap.add_argument("--out", required=True)
    a = ap.parse_args()
    spans, rows, bam = collect(a)
    need = collections.defaultdict(list)
    for cid, r in rows.items():
        for x in r["sec"]:
            need[x["name"]].append((cid, x))
    prim_loc = {}
    for s in spans.values():
        for rd in bam.fetch(s["chrom"], s["lo"], s["hi"]):
            if not (rd.is_unmapped or rd.is_supplementary or rd.is_secondary) and rd.query_name in need:
                prim_loc[rd.query_name] = (s["chrom"], rd.reference_start, rd.reference_end)
    best = best_as(a.as_table, set(need)) if a.as_table else {}
    share, good, nsec = collections.Counter(), collections.Counter(), collections.Counter()
    for nm, lst in need.items():
        for cid, x in lst:
            share[where(prim_loc.get(nm), cid, spans)] += 1
            nsec[cid] += 1
            if nm in best and x["AS"] >= GOOD_RATIO * best[nm]:
                good[cid] += 1
    total = sum(share.values())

    def maxchain(r, key, mapq0=False):
        c = collections.Counter(x[key] for x in r["prim"] if (x["mapq"] > 0 or not mapq0))
        return max(c.values()) if c else 0

    per_copy = {}
    for cid, r in rows.items():
        per_copy[r["name"]] = dict(k=r["k"], primary=len(r["prim"]), secondary=len(r["sec"]), good=(good.get(cid, 0) if best else None),
                                   chain_whole=maxchain(r, "whole"), chain_inside=maxchain(r, "inside"),
                                   chain_whole_mapq0=maxchain(r, "whole", True), chain_inside_mapq0=maxchain(r, "inside", True),
                                   primaries_mapq_pos=sum(1 for x in r["prim"] if x["mapq"] > 0))
    out = dict(family=a.family, expressed=len(rows), secondary_alignments=total, share=dict(share),
               other_copy_fraction=(share["other_copy"] / total if total else None),
               e_sensitivity=dict(E=len(rows), mapq_positive=sum(1 for v in per_copy.values() if v["primaries_mapq_pos"] >= 2),
                                  floor_3=sum(1 for v in per_copy.values() if v["primary"] >= 3),
                                  exactly_2=sum(1 for v in per_copy.values() if v["primary"] == 2)),
               per_copy=per_copy)
    with open(a.out, "w") as fh:
        json.dump(out, fh, indent=1, sort_keys=True)
    print(f"{a.family}: {len(rows)} expressed copies; {total} secondary alignments carry >= k annotated introns; primary of the molecule: {dict(share)}"
          + (f"; other-copy share {out['other_copy_fraction']:.3f}" if total else ""))
    print(f"E sensitivity: E = {out['e_sensitivity']['E']}, MAPQ > 0 required {out['e_sensitivity']['mapq_positive']}, floor of 3 primaries "
          f"{out['e_sensitivity']['floor_3']}, copies with exactly 2 primaries {out['e_sensitivity']['exactly_2']}")
    print(f"{'copy':14}{'k':>2}{'prim':>6}{'sec':>6}{'good':>6}{'chain whole/inside':>20}")
    for nm, v in per_copy.items():
        print(f"{nm:14}{v['k']:>2}{v['primary']:>6}{v['secondary']:>6}{('-' if v['good'] is None else v['good']):>6}{v['chain_whole']:>10}/{v['chain_inside']}")
    if per_copy:
        print(f"median primary {statistics.median(v['primary'] for v in per_copy.values()):g}, median secondary {statistics.median(v['secondary'] for v in per_copy.values()):g}")


if __name__ == "__main__":
    main()
