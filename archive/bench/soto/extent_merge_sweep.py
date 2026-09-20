#!/usr/bin/env python3
"""Control over-extension by how alignment blocks are combined into one locus boundary.

Usage: extent_merge_sweep.py <copies.tsv> [chr1,chr15]

Any-locus bounding lifted median size from 0.55 to 1.03 but pushed >2x from 5% to 37%. The suspected cause
is the COMBINING step, not the bounding idea: blocks were merged across gaps up to 500 bp and the largest
merged run was taken, so an extent can grow along a whole duplicated CLUSTER. AMY2A (0.52 -> 5.22) is the
amylase cluster -- its neighbours are themselves duplicated, so "is this duplicated?" keeps them.

Two dials, swept together:
  merge_gap  -- 0 means take the largest SINGLE alignment block (one copy-correspondence, no growth);
                larger values tolerate gapped alignments at the cost of letting neighbours join.
  max_part   -- optionally use only the best N partner loci rather than every locus that aligns. A real
                copy correspondence should be carried by its best partner; needing many partners to cover
                the extent is the signature of a cluster walk.

Over-extension is the dangerous direction: coverage is aligned/min(qlen,tlen), so an oversized locus
inflates its own denominator and LOSES edges -- the mechanism that cost exon-union 20 recall points. A
setting that trades a little median accuracy for far less >2x is therefore the better trade here.

Reuses the cached read extents and PAF; no BAM or aligner access.
"""
import sys
from collections import defaultdict

COP = sys.argv[1]
CHROMS = set(sys.argv[2].split(",")) if len(sys.argv) > 2 else None
BED = "bench/soto/80_fams.gene_preferred.bed"
CACHE = "/tmp/claude-1000/-mnt-c-Users-jfris-Desktop/161b66ca-e17e-4e70-b1cf-8b479d6c328d/scratchpad"
MIN_ID = 0.90

copies = []
for i, ln in enumerate(open(COP)):
    if i == 0:
        continue
    f = ln.rstrip("\n").split("\t")
    if len(f) < 7 or (CHROMS and f[3] not in CHROMS):
        continue
    copies.append({"chrom": f[3], "s": int(f[4]), "e": int(f[5])})

ext = {}
for ln in open(f"{CACHE}/read_extents.tsv"):
    k, a, b = ln.split()
    ext[k] = (int(a), int(b))
name = {}
for i, c in enumerate(copies):
    c["rlo"], c["rhi"] = ext[f"{c['chrom']}:{c['s']}-{c['e']}"]
    name[f"{c['chrom']}:{c['rlo']}-{c['rhi']}"] = i

recs = defaultdict(list)
for ln in open(f"{CACHE}/anyloc.paf"):
    f = ln.split("\t")
    def back(n):
        ch, rng = n.rsplit(":", 1)
        a, b = rng.split("-")
        return f"{ch}:{int(a)-1}-{b}"
    q, t = back(f[0]), back(f[5])
    if q not in name or t not in name or q == t:
        continue
    ident = 1.0 - float(next((x[5:] for x in f[12:] if x.startswith("de:f:")), "1"))
    if ident < MIN_ID:
        continue
    recs[name[q]].append((name[t], int(f[2]), int(f[3])))

members = []
for ln in open(BED):
    ch, s, e, nm, *_ = ln.rstrip("\n").split("\t")
    if CHROMS and ch not in CHROMS:
        continue
    members.append((ch, int(s), int(e), nm.split("|")[0]))
best_copy = {}
for (ch, ms, me, g) in members:
    ov = [(min(me, c["e"]) - max(ms, c["s"]), i) for i, c in enumerate(copies)
          if c["chrom"] == ch and ms < c["e"] and me > c["s"]]
    if ov:
        best_copy[g] = (max(ov)[1], me - ms)


def bounds(merge_gap, max_part):
    out = {}
    for i, c in enumerate(copies):
        v = recs.get(i, [])
        if not v:
            out[i] = (c["s"], c["e"])
            continue
        if max_part:
            span = defaultdict(int)
            for (t, a, b) in v:
                span[t] += b - a
            keep = {t for t, _ in sorted(span.items(), key=lambda x: -x[1])[:max_part]}
            v = [x for x in v if x[0] in keep]
        iv = sorted((a, b) for (_, a, b) in v)
        merged = [list(iv[0])]
        for a, b in iv[1:]:
            if a <= merged[-1][1] + merge_gap:
                merged[-1][1] = max(merged[-1][1], b)
            else:
                merged.append([a, b])
        bst = max(merged, key=lambda m: m[1] - m[0])
        out[i] = (c["rlo"] + bst[0], c["rlo"] + bst[1])
    return out


def stats(bd):
    v = sorted((bd[i][1] - bd[i][0]) / T for (i, T) in best_copy.values())
    n = len(v)
    return (v[n // 2], 100 * sum(1 for x in v if x >= 0.8) / n,
            100 * sum(1 for x in v if x < 0.5) / n, 100 * sum(1 for x in v if x > 2) / n,
            100 * sum(1 for x in v if 0.5 <= x <= 2) / n)


base = {i: (c["s"], c["e"]) for i, c in enumerate(copies)}
print(f"evaluated truth members: {len(best_copy)}\n")
print(f"{'setting':<34}{'median':>8}{'>=0.8x':>8}{'<0.5x':>7}{'>2x':>6}{'in 0.5-2x':>11}")
m, a, b, c, d = stats(base)
print(f"{'representative (shipped)':<34}{m:>8.2f}{a:>7.0f}%{b:>6.0f}%{c:>5.0f}%{d:>10.0f}%")
for gap in (0, 100, 500, 5000):
    for mp in (0, 1, 3):
        m, a, b, c, d = stats(bounds(gap, mp))
        tag = f"gap={gap} partners={'all' if mp==0 else mp}"
        print(f"{tag:<34}{m:>8.2f}{a:>7.0f}%{b:>6.0f}%{c:>5.0f}%{d:>10.0f}%")

bd = bounds(0, 1)
print(f"\nbest-single-block, best partner only:\n{'gene':<14}{'rep':>7}{'union500':>10}{'single':>9}")
u = bounds(500, 0)
for g in ("AMY2A", "SRGAP2C", "SRGAP2D", "NOTCH2NLC", "GOLGA8K", "SRGAP2", "NOTCH2NLR"):
    if g in best_copy:
        i, T = best_copy[g]
        print(f"{g:<14}{(copies[i]['e']-copies[i]['s'])/T:>7.2f}"
              f"{(u[i][1]-u[i][0])/T:>10.2f}{(bd[i][1]-bd[i][0])/T:>9.2f}")
