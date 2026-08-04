#!/usr/bin/env python3
"""Soto's criterion, applied honestly: is an exon shared IN FULL, or only in part?

Usage: full_exon_share.py <copies.tsv> <genome.fa> [chr1,chr15]

WHY THIS MATTERS. Soto's step is `bedtools intersect -f 0.99` -- an exon counts as shared only when it is
shared over essentially its whole length, and that full-length requirement is a large part of why the
resulting families are believed. Our `RUSTLE_SHARED_EXON` rule instead accepts on ABSOLUTE aligned length
(`nucleotide_edges_indexed`, default 100 bp), so a 3 kb exon matching over 300 bp of a shared domain or an
exonised repeat passes. That is a materially weaker claim than the one Soto verified.

TWO QUESTIONS, ONE MEASUREMENT.
  1. Does requiring FULL sharing discriminate true co-family pairs from false ones better than partial?
  2. Does a locus representative actually represent its locus? If a representative is a fragment or the
     wrong isoform, its exons will not be shared in full with its true siblings even when the loci really
     are copies -- so the full-share rate over TRUE pairs is a direct check on representative quality,
     independent of any grouping decision.

An exon pair is shared at fraction f when the aligned block covers >= f of the SHORTER exon at >= min_id
identity. f = 0.99 is Soto's rule; lower values show what the relaxation buys and costs.
"""
import os
import subprocess
import sys
from collections import defaultdict

COP, GENOME = sys.argv[1], sys.argv[2]
CHROMS = set(sys.argv[3].split(",")) if len(sys.argv) > 3 else None
BED = "bench/soto/80_fams.gene_preferred.bed"
CACHE = "/tmp/claude-1000/-mnt-c-Users-jfris-Desktop/161b66ca-e17e-4e70-b1cf-8b479d6c328d/scratchpad"
MIN_EXON = 80          # below this an "exon" is mostly alignment noise
TAG = os.path.basename(COP).split(".")[0]

copies = []
for i, ln in enumerate(open(COP)):
    if i == 0:
        continue
    f = ln.rstrip("\n").split("\t")
    if len(f) < 10 or (CHROMS and f[3] not in CHROMS):
        continue
    exons = [tuple(int(x) for x in b.split("-")) for b in f[9].split(",")]
    copies.append({"chrom": f[3], "s": int(f[4]), "e": int(f[5]), "exons": exons})

members = []
for ln in open(BED):
    c, s, e, nm, *_ = ln.rstrip("\n").split("\t")
    if CHROMS and c not in CHROMS:
        continue
    members.append((c, int(s), int(e), nm.split("|")[1]))
for cp in copies:
    hit = [(min(cp["e"], me) - max(cp["s"], ms), fam)
           for (c, ms, me, fam) in members if c == cp["chrom"] and cp["s"] < me and cp["e"] > ms]
    cp["truth"] = max(hit)[1] if hit else None

labelled = [i for i, c in enumerate(copies) if c["truth"]]
print(f"{TAG}: {len(copies)} copies, {len(labelled)} with a truth label")

# ---- one FASTA of every exon of every labelled locus -------------------------------------------------
fa = f"{CACHE}/exons_{TAG}.fa"
owner, exlen = [], []
if not os.path.exists(fa):
    regions = []
    for li in labelled:
        for (a, b) in copies[li]["exons"]:
            if b - a >= MIN_EXON:
                regions.append(f"{copies[li]['chrom']}:{a+1}-{b}")
    open(f"{CACHE}/exreg_{TAG}.txt", "w").write("\n".join(regions) + "\n")
    subprocess.run(f"samtools faidx {GENOME} -r {CACHE}/exreg_{TAG}.txt > {fa} 2>/dev/null",
                   shell=True, check=True)
name2owner = {}
for li in labelled:
    for (a, b) in copies[li]["exons"]:
        if b - a >= MIN_EXON:
            name2owner[f"{copies[li]['chrom']}:{a+1}-{b}"] = li
print(f"exons >= {MIN_EXON} bp: {len(name2owner)}")

paf = f"{CACHE}/exons_{TAG}.paf"
if not os.path.exists(paf):
    with open(paf, "w") as fh:
        fh.write(subprocess.run(
            ["minimap2", "-x", "asm20", "-c", "--eqx", "-X", "--no-long-join", "-t", "3",
             "-N", "50", "-p", "0.1", fa, fa], capture_output=True, text=True, check=True).stdout)

# best shared fraction for each locus pair, at each identity floor
best = defaultdict(lambda: defaultdict(float))
for ln in open(paf):
    f = ln.split("\t")
    qa, ta = name2owner.get(f[0]), name2owner.get(f[5])
    if qa is None or ta is None or qa == ta:
        continue
    ident = 1.0 - float(next((x[5:] for x in f[12:] if x.startswith("de:f:")), "1"))
    frac = (int(f[3]) - int(f[2])) / min(int(f[1]), int(f[6]))
    k = (min(qa, ta), max(qa, ta))
    for mid in (0.90, 0.95, 0.98):
        if ident >= mid and frac > best[k][mid]:
            best[k][mid] = frac

pairs = [(a, b) for ii, a in enumerate(labelled) for b in labelled[ii + 1:]]
print(f"labelled locus pairs: {len(pairs)}\n")
print(f"{'rule':<34}{'TRUE linked':>13}{'FALSE linked':>14}{'precision':>11}")
for mid in (0.90, 0.95, 0.98):
    for fr in (0.20, 0.50, 0.80, 0.99):
        t = f = 0
        for (a, b) in pairs:
            if best.get((a, b), {}).get(mid, 0.0) >= fr:
                if copies[a]["truth"] == copies[b]["truth"]:
                    t += 1
                else:
                    f += 1
        lab = f"id>={mid:.2f}, exon shared >= {int(fr*100):>2}%"
        print(f"{lab:<34}{t:>13}{f:>14}{t/max(t+f,1):>11.3f}")

# Representative quality: over pairs that ARE co-family, how fully is an exon shared?
print("\nrepresentative check -- best shared exon fraction over TRUE co-family pairs (id >= 0.95):")
v = sorted(best.get((a, b), {}).get(0.95, 0.0) for (a, b) in pairs
           if copies[a]["truth"] == copies[b]["truth"])
n = len(v)
if n:
    print(f"  n={n}  median={v[n//2]:.2f}  "
          f">=0.99 (Soto's bar): {100*sum(1 for x in v if x>=0.99)/n:.0f}%  "
          f">=0.80: {100*sum(1 for x in v if x>=0.80)/n:.0f}%  "
          f"no shared exon at all: {100*sum(1 for x in v if x==0)/n:.0f}%")
    print("  A true pair failing the full-share bar means at least one representative is a fragment or the")
    print("  wrong isoform -- the loci are copies, so their exons should be shared in full.")
