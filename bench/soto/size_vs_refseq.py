#!/usr/bin/env python3
"""Are predicted loci the size they should be? Scored against RefSeq gene spans, no Soto file involved.

Usage: size_vs_refseq.py <gff.gz> <copies.tsv> [copies.tsv ...]

Every gene on the covered chromosomes that overlaps an emitted copy is scored as predicted_span /
RefSeq_gene_span. Reported as the IN-BAND fraction (0.5x to 2x), not the median: the median only says
where the mass of the distribution sits, and an earlier version of this analysis reported a median
improvement 0.55 -> 1.03 that was a pure distribution shift -- the shipped representative was actually
MORE often correct (47% in band vs 37%).

Broken out by whether the representative is spliced, because that split has been the dominant one all
along: stub-represented loci run ~0.29x of truth and spliced ones ~0.75x, and the same gene comes out at
0.86-1.04 with a spliced model and 0.03-0.52 with a stub.

A gene is counted once, against the copy that overlaps it most.
"""
import gzip
import re
import sys
from collections import defaultdict

GFF = sys.argv[1]
CATS = sys.argv[2:]

genes = {}
with gzip.open(GFF, "rt") as fh:
    for ln in fh:
        if ln[0] == "#":
            continue
        f = ln.rstrip("\n").split("\t")
        if len(f) < 9 or f[2] != "gene":
            continue
        m = re.search(r"Name=([^;]+)", f[8])
        if m:
            genes.setdefault(m.group(1), (f[0], int(f[3]) - 1, int(f[4])))
by_chrom = defaultdict(list)
for g, (c, s, e) in genes.items():
    by_chrom[c].append((s, e, g))
for c in by_chrom:
    by_chrom[c].sort()
print(f"RefSeq genes loaded: {len(genes)}\n")


def band(v):
    return 100 * sum(1 for x in v if 0.5 <= x <= 2) / len(v)


print(f"{'catalog':<12}{'genes':>7}{'in 0.5-2x':>11}{'<0.5x':>8}{'>2x':>7}{'median':>9}"
      f"{'  |  spliced rep':>17}{'stub rep':>12}")
for cat in CATS:
    copies = []
    for i, ln in enumerate(open(cat)):
        if i == 0:
            continue
        f = ln.rstrip("\n").split("\t")
        if len(f) >= 7:
            copies.append((f[3], int(f[4]), int(f[5]), int(f[6])))
    seen = {}
    for (c, s, e, nx) in copies:
        for (gs, ge, g) in by_chrom.get(c, []):
            if gs >= e:
                break
            if s < ge and e > gs:
                ov = min(e, ge) - max(s, gs)
                if g not in seen or ov > seen[g][0]:
                    seen[g] = (ov, (e - s) / (ge - gs), nx)
    v = sorted(x[1] for x in seen.values())
    sp = sorted(x[1] for x in seen.values() if x[2] > 1)
    st = sorted(x[1] for x in seen.values() if x[2] == 1)
    n = len(v)
    tag = cat.split("/")[-1].split(".")[0]
    print(f"{tag:<12}{n:>7}{band(v):>10.0f}%{100*sum(1 for x in v if x<0.5)/n:>7.0f}%"
          f"{100*sum(1 for x in v if x>2)/n:>6.0f}%{v[n//2]:>9.2f}"
          f"{f'{band(sp):.0f}% (n={len(sp)})' if sp else '-':>17}"
          f"{f'{band(st):.0f}% (n={len(st)})' if st else '-':>12}")
