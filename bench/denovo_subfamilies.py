#!/usr/bin/env python3
"""De novo O1 evaluation against a truth table (prereg Addendum V): family-level recovery of a `gw_family_catalog`
catalog and subfamily clades from exon and intron trees of its copies (D2).

usage: denovo_subfamilies.py --catalog PREFIX --truth truth.tsv --genome genome.fa --iqtree iqtree3 --outdir DIR [--tag T]
PREFIX.copies.tsv / PREFIX.copies.fa are the catalog outputs (same row order). Each truth record is represented by its
best-overlapping copy; per emitted family holding truth records, exon sequence = the copy's spliced sequence, intron
sequence = the copy's genomic span minus its exon blocks (transcript orientation); trees and clade calls are those of
`bench/guided_pipeline.py` (reference-projected alignment, IQ-TREE, SH-aLRT > 75).
"""
import argparse
import collections
import csv
import os
import sys

import pysam

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import guided_pipeline as gp  # noqa: E402


def span_minus_exons(genome, c):
    """A copy's genomic span with its exon blocks removed, in transcript orientation."""
    pieces, cur = [], c["start"]
    for a, b in gp.merge(c["blocks"]):
        if a > cur:
            pieces.append((cur, a))
        cur = max(cur, b)
    if c["end"] > cur:
        pieces.append((cur, c["end"]))
    seq = "".join(genome.fetch(c["chrom"], a, b) for a, b in pieces).upper()
    return gp.rc(seq) if c["strand"] == "-" else seq


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--catalog", required=True)
    ap.add_argument("--truth", required=True)
    ap.add_argument("--genome", required=True)
    ap.add_argument("--iqtree", required=True)
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--tag", default="dn")
    ap.add_argument("--threads", type=int, default=4)
    a = ap.parse_args()
    os.makedirs(a.outdir, exist_ok=True)
    genome = pysam.FastaFile(a.genome)
    truth = list(csv.DictReader(open(a.truth), delimiter="\t"))
    for t in truth:
        t["start0"], t["end"] = int(t["start0"]), int(t["end"])
    rec = {t["name"]: t for t in truth}
    copies = list(csv.DictReader(open(f"{a.catalog}.copies.tsv"), delimiter="\t"))
    seqs = [l.strip() for l in open(f"{a.catalog}.copies.fa") if not l.startswith(">")]
    assert len(seqs) == len(copies), (len(seqs), len(copies))
    for k, c in enumerate(copies):
        c["k"], c["start"], c["end"], c["seq"] = k, int(c["start"]), int(c["end"]), seqs[k]
        c["blocks"] = [tuple(map(int, x.split("-"))) for x in c["exons"].split(",") if x]

    # ---------------- family level ----------------
    best, nfam, ncop = {}, {}, {}
    for t in truth:
        h = [(gp.ov(t["start0"], t["end"], c["start"], c["end"]), -c["k"], c) for c in copies if c["chrom"] == t["chrom"]]
        h = [x for x in h if x[0] > 0]
        if h:
            best[t["name"]] = max(h, key=lambda x: (x[0], x[1]))[2]
            nfam[t["name"]] = len({x[2]["family_id"] for x in h})
            ncop[t["name"]] = len(h)
    shared = collections.Counter(c["k"] for c in best.values())
    families = sorted({t["family"] for t in truth})
    print(f"## {a.tag}: family level ({len(copies)} copies, {len({c['family_id'] for c in copies})} families)")
    for fam in families + ["ALL"]:
        R = [t for t in truth if fam in ("ALL", t["family"])]
        P = [t for t in R if t["name"] in best]
        coll = sum(1 for t in P if shared[best[t["name"]]["k"]] > 1)
        mf = sum(nfam[t["name"]] for t in P) / len(P) if P else float("nan")
        mc = sum(ncop[t["name"]] for t in P) / len(P) if P else float("nan")
        pred = [best[t["name"]]["family_id"] if t["name"] in best else f"missing:{t['name']}" for t in R]
        truef = [t["family"] for t in R]
        ps, pp = gp.pairwise(pred, truef)
        br, bp = gp.bipartite(pred, truef)
        print(f"{fam:6s} present {len(P)}/{len(R)} | mean families/record {mf:.2f} copies/record {mc:.2f} | collapsed {coll} | "
              f"family-level pairwise sens {ps:.3f} prec {pp:.3f} | bipartite micro R {br:.3f} P {bp:.3f} | "
              f"missing {[t['name'] for t in R if t['name'] not in best]}")

    # ---------------- D2 clades ----------------
    groups = gp.literature_groups(truth)
    positional = {"TBC1D3": {t["name"] for t in truth if t["family"] == "TBC1D3" and t["level1"] == "cluster1"}}
    by_family = collections.defaultdict(list)
    for name, c in best.items():
        by_family[c["family_id"]].append(name)
    print(f"\n## {a.tag}: subfamily clades inside emitted families (SH-aLRT > 75)")
    summary = collections.defaultdict(collections.Counter)
    for fid, names in sorted(by_family.items()):
        tf = collections.Counter(rec[n]["family"] for n in names).most_common(1)[0][0]
        members = [c for c in copies if c["family_id"] == fid]
        label = {}
        for n in names:
            if rec[n]["family"] != tf:
                continue
            k = best[n]["k"]
            prev = label.get(k)
            if prev is None or gp.ov(rec[n]["start0"], rec[n]["end"], best[n]["start"], best[n]["end"]) > \
                    gp.ov(rec[prev]["start0"], rec[prev]["end"], best[n]["start"], best[n]["end"]):
                label[k] = n
        ex, it, lab = {}, {}, {}
        for c in members:
            key = f"c{c['k']}"
            ex[key] = c["seq"]
            it[key] = span_minus_exons(genome, c)
            if c["k"] in label:
                lab[key] = label[c["k"]]
        present_truth = sorted(set(lab.values()))
        calls = {}
        for cls, sq in (("exon", ex), ("intron", it)):
            sq = {k: v for k, v in sq.items() if len(v) >= 100}
            if len(sq) < 4 or len({lab[k] for k in sq if k in lab}) < 2:
                print(f"{fid} [{tf}] {cls}: {len(sq)} members, {len(present_truth)} truth -> not treed")
                continue
            tree, ref, kept, dropped = gp.projected_tree(f"{a.tag}_{fid}_{cls}", sq, a.outdir, a.iqtree, a.threads)
            if tree is None:
                print(f"{fid} [{tf}] {cls}: dropped {dropped} -> not treed")
                continue
            res, pos = gp.clade_calls(*tree, {k: v for k, v in lab.items() if k in sq and k not in dropped}, groups[tf],
                                      positional.get(tf))
            calls[cls] = (res, pos)
            print(f"{fid} [{tf}] {cls}: {len(sq) - len(dropped)} leaves ({len(present_truth)} truth: {' '.join(present_truth)}), "
                  f"{kept} columns | " + "; ".join(f"{g}: {v}" for g, v in res.items()) + ("" if pos is None else f" | positional supported: {pos}"))
        for g in groups[tf]:
            cc = {cls: v[0][g] for cls, v in calls.items()}
            if not cc or all(x == "absent" for x in cc.values()):
                continue
            for cls, x in cc.items():
                summary[(tf, g, cls)][x] += 1
            summary[(tf, g, "either")]["RECOVERED" if "RECOVERED" in cc.values() else "not"] += 1
        for cls, v in calls.items():
            if v[1] is not None:
                summary[(tf, "positional split", cls)]["WRONG" if v[1] else "CORRECT"] += 1
    print(f"\n## {a.tag}: clade summary (counts over emitted families holding >= 2 truth records of a group)")
    for k in sorted(summary):
        print(f"{k[0]:6s} {k[1]:28s} {k[2]:6s} {dict(summary[k])}")


if __name__ == "__main__":
    main()
