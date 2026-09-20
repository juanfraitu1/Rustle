#!/usr/bin/env python3
"""Bound a locus by asking whether its sequence is DUPLICATED AT ALL, not by asking which family it is in.

Usage: anylocus_extent.py <copies.tsv> <bam> <genome.fa> [chr1,chr15]

WHY NOT THE SIBLING VERSION. Bounding a locus by its predicted family siblings needs correct families to
fix the boundaries that feed the substrate that produces the families. That is a feedback loop, not just a
chicken-and-egg: a wrong family bounds a locus by the wrong siblings, shifting its edges and further
cementing the wrong family. Measured: SRGAP2's four members sit in THREE predicted families, and three of
the four fell back to their representative span, so the family-bounded median of 0.91 was largely measuring
the fallback rather than the method.

THE REFORMULATION. Read extent overshoots because reads spill into NEIGHBOURING transcription, and that
neighbouring sequence is typically unique. A real locus is a copy, so it is duplicated somewhere; the spill
is not. Asking "is this duplicated anywhere in the catalog?" needs no family assignment at all, so the
feedback loop disappears -- and it is the question O1 asks anyway.

GUARD. With no family restriction, common repeats could bound everything. MIN_BLOCK sweeps the minimum
aligned block length so repeat sensitivity is visible rather than assumed.

Read extents are cached to a TSV, so re-running any variant after the first call costs seconds.
"""
import os
import subprocess
import sys
from collections import defaultdict

COP, BAM, GENOME = sys.argv[1], sys.argv[2], sys.argv[3]
CHROMS = set(sys.argv[4].split(",")) if len(sys.argv) > 4 else None
BED = "bench/soto/80_fams.gene_preferred.bed"
CACHE = "/tmp/claude-1000/-mnt-c-Users-jfris-Desktop/161b66ca-e17e-4e70-b1cf-8b479d6c328d/scratchpad"
MIN_ID = 0.90
EXT_TSV = f"{CACHE}/read_extents.tsv"

copies = []
for i, ln in enumerate(open(COP)):
    if i == 0:
        continue
    f = ln.rstrip("\n").split("\t")
    if len(f) < 7 or (CHROMS and f[3] not in CHROMS):
        continue
    copies.append({"fam": f[0], "chrom": f[3], "s": int(f[4]), "e": int(f[5]), "nexon": int(f[6])})

# ---- read extent, cached ---------------------------------------------------------------------------
cached = {}
if os.path.exists(EXT_TSV):
    for ln in open(EXT_TSV):
        a, b, c = ln.split()
        cached[a] = (int(b), int(c))
dirty = False
for c in copies:
    k = f"{c['chrom']}:{c['s']}-{c['e']}"
    if k in cached:
        c["rlo"], c["rhi"] = cached[k]
        continue
    dirty = True
    out = subprocess.run(["samtools", "view", "-F", "2308", BAM, f"{c['chrom']}:{c['s']}-{c['e']}"],
                         capture_output=True, text=True).stdout
    starts, ends, spans = defaultdict(int), defaultdict(int), []
    for ln in out.splitlines():
        f = ln.split("\t")
        p = int(f[3])
        q, n = p, 0
        for ch in f[5]:
            if ch.isdigit():
                n = n * 10 + ord(ch) - 48
            else:
                if ch in "MDN=X":
                    q += n
                n = 0
        spans.append((p, q))
        starts[p // 100] += 1
        ends[q // 100] += 1
    keep = [(a, b) for (a, b) in spans if starts[a // 100] >= 3 or ends[b // 100] >= 3]
    c["rlo"], c["rhi"] = (min(a for a, _ in keep), max(b for _, b in keep)) if keep else (c["s"], c["e"])
    cached[k] = (c["rlo"], c["rhi"])
if dirty:
    with open(EXT_TSV, "w") as fh:
        for k, (a, b) in cached.items():
            fh.write(f"{k}\t{a}\t{b}\n")

fa, paf = f"{CACHE}/anyloc.fa", f"{CACHE}/anyloc.paf"
name = {}
if not os.path.exists(paf):
    with open(f"{CACHE}/anyloc.txt", "w") as fh:
        for c in copies:
            fh.write(f"{c['chrom']}:{c['rlo']+1}-{c['rhi']}\n")
    subprocess.run(f"samtools faidx {GENOME} -r {CACHE}/anyloc.txt > {fa} 2>/dev/null",
                   shell=True, check=True)
    with open(paf, "w") as fh:
        fh.write(subprocess.run(
            ["minimap2", "-x", "asm20", "-c", "--eqx", "-X", "--no-long-join", "-t", "4",
             "-N", "50", "-p", "0.1", fa, fa], capture_output=True, text=True, check=True).stdout)
for i, c in enumerate(copies):
    name[f"{c['chrom']}:{c['rlo']}-{c['rhi']}"] = i

records = []
for ln in open(paf):
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
    records.append((name[q], name[t], int(f[2]), int(f[3])))
print(f"copies {len(copies)}   alignment blocks at id>={MIN_ID}: {len(records)}\n")

members = []
for ln in open(BED):
    ch, s, e, nm, *_ = ln.rstrip("\n").split("\t")
    if CHROMS and ch not in CHROMS:
        continue
    members.append((ch, int(s), int(e), nm.split("|")[0]))


def bound(min_block, same_family_only):
    iv = defaultdict(list)
    for (qi, ti, a, b) in records:
        if b - a < min_block:
            continue
        if same_family_only and copies[qi]["fam"] != copies[ti]["fam"]:
            continue
        iv[qi].append((a, b))
    lo_hi, fallback = {}, 0
    for i, c in enumerate(copies):
        v = sorted(iv.get(i, []))
        if not v:
            lo_hi[i] = (c["s"], c["e"])
            fallback += 1
            continue
        merged = [list(v[0])]
        for a, b in v[1:]:
            if a <= merged[-1][1] + 500:
                merged[-1][1] = max(merged[-1][1], b)
            else:
                merged.append([a, b])
        best = max(merged, key=lambda m: m[1] - m[0])
        lo_hi[i] = (c["rlo"] + best[0], c["rlo"] + best[1])
    return lo_hi, fallback


def evaluate(lo_hi):
    out = []
    for (ch, ms, me, g) in members:
        ov = [(min(me, c["e"]) - max(ms, c["s"]), i) for i, c in enumerate(copies)
              if c["chrom"] == ch and ms < c["e"] and me > c["s"]]
        if not ov:
            continue
        i = max(ov)[1]
        lo, hi = lo_hi[i]
        out.append(((hi - lo) / (me - ms), g, i))
    return out


def med(v):
    v = sorted(v)
    return v[len(v) // 2] if v else 0.0


base = evaluate({i: (c["s"], c["e"]) for i, c in enumerate(copies)})
reads = evaluate({i: (c["rlo"], c["rhi"]) for i, c in enumerate(copies)})
print(f"{'estimator':<34}{'median':>8}{'>=0.8x':>8}{'<0.5x':>7}{'>2x':>6}{'fallback':>10}")


def row(lab, res, fb=None):
    v = [r[0] for r in res]
    print(f"{lab:<34}{med(v):>8.2f}{100*sum(1 for x in v if x>=0.8)/len(v):>7.0f}%"
          f"{100*sum(1 for x in v if x<0.5)/len(v):>6.0f}%{100*sum(1 for x in v if x>2)/len(v):>5.0f}%"
          f"{('' if fb is None else f'{fb}/{len(copies)}'):>10}")


row("representative (shipped)", base)
row("read extent, >=3 support", reads)
for mb in (0, 500, 1000, 2000, 5000):
    lh, fb = bound(mb, False)
    row(f"ANY-LOCUS bounded, block>={mb}", evaluate(lh), fb)
lh, fb = bound(1000, True)
row("family-bounded, block>=1000", evaluate(lh), fb)

lh, _ = bound(1000, False)
res = {g: r for r, g, i in evaluate(lh)}
rb = {g: r for r, g, i in base}
print(f"\n{'gene':<14}{'rep':>7}{'any-locus':>11}")
for g in ("SRGAP2", "SRGAP2B", "SRGAP2C", "SRGAP2D", "NOTCH2NLC", "NOTCH2NLR", "GOLGA8K", "AMY2A"):
    if g in res:
        print(f"{g:<14}{rb[g]:>7.2f}{res[g]:>11.2f}")
