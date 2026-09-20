#!/usr/bin/env python3
"""Does getting a locus's SIZE right make its FAMILY call better? Or is size cosmetic?

Usage: size_vs_accuracy.py <copies.tsv> <gff> [chr1,chr15]

The size work is only worth presenting if size accuracy CAUSES family accuracy. This measures that
directly instead of asserting it. Each emitted copy is given a size class against the annotation:

    truncated      predicted span < 0.5x the gene
    in-band        0.5x to 2x
    over-extended  > 2x

then every pair of labelled copies is scored two ways, stratified by the size classes of its two members:

    PRECISION  of pairs the pipeline PUT IN ONE FAMILY, how many are truly co-family. If mis-sized loci
               generate the false edges, precision falls away from the in-band/in-band cell.
    RECALL     of pairs that ARE truly co-family, how many the pipeline grouped. If mis-sizing costs
               members, recall falls the same way.

Note what this can and cannot show. A correlation between size class and pair correctness is consistent
with size accuracy CAUSING family accuracy, but also with a third factor (a well-assembled locus is both
correctly sized and easy to group). The honest claim is association, and the strength of it.
"""
import gzip
import re
import sys
from collections import defaultdict
from itertools import combinations

COP, GFF = sys.argv[1], sys.argv[2]
CHROMS = set(sys.argv[3].split(",")) if len(sys.argv) > 3 else {"chr1", "chr15"}
BED = "bench/soto/80_fams.gene_preferred.bed"

genes = defaultdict(list)
opener = gzip.open if GFF.endswith(".gz") else open
with opener(GFF, "rt") as fh:
    for ln in fh:
        if ln[0] == "#":
            continue
        f = ln.rstrip("\n").split("\t")
        if len(f) < 9 or f[2] != "gene" or f[0] not in CHROMS:
            continue
        m = (re.search(r"Name=([^;]+)", f[8]) or re.search(r"source_gene_common_name=([^;]+)", f[8])
             or re.search(r"gene_id=([^;]+)", f[8]))
        if m:
            genes[f[0]].append((int(f[3]) - 1, int(f[4])))
for c in genes:
    genes[c].sort()

copies = []
for i, ln in enumerate(open(COP)):
    if i == 0:
        continue
    f = ln.rstrip("\n").split("\t")
    if len(f) < 7 or f[3] not in CHROMS:
        continue
    s, e = int(f[4]), int(f[5])
    ov = [(min(e, ge) - max(s, gs), ge - gs) for (gs, ge) in genes.get(f[3], []) if s < ge and e > gs]
    if not ov:
        continue
    ratio = (e - s) / max(ov, key=lambda x: x[0])[1]
    copies.append({"fam": f[0], "chrom": f[3], "s": s, "e": e, "ratio": ratio})

members = []
for ln in open(BED):
    c, s, e, nm, *_ = ln.rstrip("\n").split("\t")
    if c in CHROMS:
        members.append((c, int(s), int(e), nm.split("|")[1]))
for cp in copies:
    h = [(min(cp["e"], me) - max(cp["s"], ms), fam) for (c, ms, me, fam) in members
         if c == cp["chrom"] and cp["s"] < me and cp["e"] > ms]
    cp["truth"] = max(h)[1] if h else None


def cls(r):
    return "trunc" if r < 0.5 else ("in-band" if r <= 2 else "over")


lab = [c for c in copies if c["truth"]]
print(f"labelled copies: {len(lab)}   "
      f"({sum(1 for c in lab if cls(c['ratio'])=='trunc')} truncated, "
      f"{sum(1 for c in lab if cls(c['ratio'])=='in-band')} in-band, "
      f"{sum(1 for c in lab if cls(c['ratio'])=='over')} over-extended)\n")

prec = defaultdict(lambda: [0, 0])   # size-class pair -> [true, total grouped]
rec = defaultdict(lambda: [0, 0])    # size-class pair -> [grouped, total true]
for a, b in combinations(lab, 2):
    key = tuple(sorted((cls(a["ratio"]), cls(b["ratio"]))))
    same_fam = a["fam"] == b["fam"]
    same_truth = a["truth"] == b["truth"]
    if same_fam:
        prec[key][1] += 1
        prec[key][0] += same_truth
    if same_truth:
        rec[key][1] += 1
        rec[key][0] += same_fam

order = [("in-band", "in-band"), ("in-band", "over"), ("in-band", "trunc"),
         ("over", "over"), ("over", "trunc"), ("trunc", "trunc")]
print(f"{'size classes of the pair':<28}{'grouped':>9}{'precision':>11}"
      f"{'true pairs':>12}{'recall':>9}")
for k in order:
    kk = tuple(sorted(k))
    p, t = prec[kk], rec[kk]
    ps = f"{p[0]/p[1]:.3f}" if p[1] else "-"
    rs = f"{t[0]/t[1]:.3f}" if t[1] else "-"
    print(f"{' + '.join(k):<28}{p[1]:>9}{ps:>11}{t[1]:>12}{rs:>9}")

both = [k for k in prec if k == ("in-band", "in-band")]
any_bad = [k for k in prec if k != ("in-band", "in-band")]
gp = sum(prec[k][0] for k in both), sum(prec[k][1] for k in both)
bp = sum(prec[k][0] for k in any_bad), sum(prec[k][1] for k in any_bad)
gr = sum(rec[k][0] for k in both), sum(rec[k][1] for k in both)
br = sum(rec[k][0] for k in any_bad), sum(rec[k][1] for k in any_bad)
print(f"\n{'both members in-band':<28}{gp[1]:>9}{gp[0]/max(gp[1],1):>11.3f}"
      f"{gr[1]:>12}{gr[0]/max(gr[1],1):>9.3f}")
print(f"{'at least one mis-sized':<28}{bp[1]:>9}{bp[0]/max(bp[1],1):>11.3f}"
      f"{br[1]:>12}{br[0]/max(br[1],1):>9.3f}")
