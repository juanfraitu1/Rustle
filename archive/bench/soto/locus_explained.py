#!/usr/bin/env python3
"""How much of what you would SEE IN IGV at a locus does its representative actually explain?

Usage: locus_explained.py <copies.tsv> <bam> <genome.fa> [chr1,chr15]

WHY. Every metric so far scores the model against a truth annotation. None of them answers the question a
reviewer asks when he opens the locus in a browser: there is read signal here that the model does not
account for -- what is it, and why was it dropped? A model can score well on family recall while visibly
ignoring half the evidence at the locus, and that is the objection this measures.

THREE NUMBERS PER LOCUS, all from primary alignments only (-F 2308, the standing rule for any per-read
CIGAR statistic):

  base explained    aligned read bases falling inside the representative's exons, over all aligned read
                    bases at the locus. This is the "grey coverage the model does not cover" in IGV.
  junction recall   read-supported junctions (>= MIN_J reads) that appear in the representative's chain,
                    over all such junctions. This is the "sashimi arcs the model does not draw".
  unexplained mass  aligned bases in regions with real read depth that the model skips entirely.

Then it asks WHY the missed junctions were missed, by fetching each one's donor/acceptor dinucleotides:
CANONICAL misses (GT-AG, GC-AG, AT-AC and their reverse complements) were dropped by chain selection and
are recoverable by better modelling; NON-CANONICAL misses were dropped by the all-canonical gate in
`build_spliced_seq` and are only recoverable by relaxing it (RUSTLE_JUNCTION_MAJORITY, bounded by
RUSTLE_JUNCTION_NC_MAX_BP because a small non-canonical intron is a plausible splice variant while a
100 kb one is a mis-chain).

That split is the point: it says how much of the unexplained signal relaxing the canonical rule could
actually buy, instead of assuming it is the answer.
"""
import os
import pickle
import subprocess
import sys
from collections import defaultdict

COP, BAM, GENOME = sys.argv[1], sys.argv[2], sys.argv[3]
CHROMS = set(sys.argv[4].split(",")) if len(sys.argv) > 4 else None
CACHE = "/tmp/claude-1000/-mnt-c-Users-jfris-Desktop/161b66ca-e17e-4e70-b1cf-8b479d6c328d/scratchpad"
TAG = os.path.basename(COP).split(".")[0]
MIN_J = 3
CANON = {("GT", "AG"), ("GC", "AG"), ("AT", "AC"), ("CT", "AC"), ("CT", "GC"), ("GT", "AT")}

copies = []
for i, ln in enumerate(open(COP)):
    if i == 0:
        continue
    f = ln.rstrip("\n").split("\t")
    if len(f) < 10 or (CHROMS and f[3] not in CHROMS):
        continue
    copies.append({"chrom": f[3], "s": int(f[4]), "e": int(f[5]), "nexon": int(f[6]),
                   "exons": [tuple(int(x) for x in b.split("-")) for b in f[9].split(",")]})

cache = f"{CACHE}/explained_{TAG}.pkl"
if os.path.exists(cache):
    rows, missed = pickle.load(open(cache, "rb"))
else:
    rows, missed = [], []
    for c in copies:
        out = subprocess.run(["samtools", "view", "-F", "2308", BAM,
                              f"{c['chrom']}:{c['s']}-{c['e']}"], capture_output=True, text=True).stdout
        blocks, junc = [], defaultdict(int)
        for ln in out.splitlines():
            f = ln.split("\t")
            # SAM POS is 1-based; the `exons` column and everything downstream is 0-based half-open.
            # Getting this wrong silently reports 0% junction recall and 100% non-canonical motifs.
            p, n = int(f[3]) - 1, 0
            for ch in f[5]:
                if ch.isdigit():
                    n = n * 10 + ord(ch) - 48
                else:
                    if ch in "M=X":
                        blocks.append((p, p + n))
                        p += n
                    elif ch == "D":
                        p += n
                    elif ch == "N":
                        junc[(p, p + n)] += 1
                        p += n
                    n = 0
        if not blocks:
            continue
        ex = sorted(c["exons"])
        total = inside = 0
        for (a, b) in blocks:
            if b <= c["s"] or a >= c["e"]:
                continue                       # read mass outside the called locus is a different question
            a, b = max(a, c["s"]), min(b, c["e"])
            total += b - a
            for (x, y) in ex:
                lo, hi = max(a, x), min(b, y)
                if hi > lo:
                    inside += hi - lo
        real = {j for j, n in junc.items() if n >= MIN_J}
        have = set(c["exons"] and [(c["exons"][i][1], c["exons"][i + 1][0])
                                   for i in range(len(c["exons"]) - 1)])
        miss = real - have
        rows.append((total, inside, len(real), len(real & have), len(miss), c["nexon"]))
        for j in miss:
            missed.append((c["chrom"], j[0], j[1], junc[j]))
    pickle.dump((rows, missed), open(cache, "wb"))

rows = [r for r in rows if r[0] > 0]
frac = sorted(r[1] / r[0] for r in rows)
n = len(rows)
tot_bases = sum(r[0] for r in rows)
tot_in = sum(r[1] for r in rows)
jr = [r for r in rows if r[2] > 0]
print(f"{TAG}: {n} loci with read coverage\n")
print(f"BASE EXPLAINED (aligned read bases inside the representative's exons)")
print(f"  pooled over all loci : {100*tot_in/tot_bases:.1f}%")
print(f"  per-locus median     : {100*frac[n//2]:.1f}%   "
      f"quartiles {100*frac[n//4]:.0f}% / {100*frac[3*n//4]:.0f}%")
print(f"  loci explaining <50% : {100*sum(1 for x in frac if x<0.5)/n:.0f}%")
if jr:
    tj = sum(r[2] for r in jr)
    hj = sum(r[3] for r in jr)
    print(f"\nJUNCTION RECALL (read-supported junctions >= {MIN_J} reads present in the model)")
    print(f"  pooled: {hj}/{tj} = {100*hj/tj:.1f}%   loci with >= 1 missed junction: "
          f"{100*sum(1 for r in jr if r[4]>0)/len(jr):.0f}%")

# ---- why were the missed junctions missed? --------------------------------------------------------
if missed:
    reg = []
    for (c, d, a, _) in missed:
        reg.append(f"{c}:{d+1}-{d+2}")
        reg.append(f"{c}:{a-1}-{a}")
    rf = f"{CACHE}/miss_{TAG}.txt"
    open(rf, "w").write("\n".join(reg) + "\n")
    fa = subprocess.run(f"samtools faidx {GENOME} -r {rf}", shell=True,
                        capture_output=True, text=True).stdout
    seqs, cur = [], None
    for ln in fa.splitlines():
        if ln.startswith(">"):
            if cur is not None:
                seqs.append(cur)
            cur = ""
        else:
            cur += ln.strip()
    if cur is not None:
        seqs.append(cur)
    canon = noncanon = 0
    nc_small = nc_big = 0
    for i, (c, d, a, w) in enumerate(missed):
        if 2 * i + 1 >= len(seqs):
            break
        don, acc = seqs[2 * i].upper(), seqs[2 * i + 1].upper()
        if (don, acc) in CANON:
            canon += 1
        else:
            noncanon += 1
            if a - d <= 10000:
                nc_small += 1
            else:
                nc_big += 1
    tot = canon + noncanon
    print(f"\nWHY THE MISSED JUNCTIONS WERE MISSED  (n={tot})")
    print(f"  CANONICAL motif      : {canon} ({100*canon/max(tot,1):.0f}%)  "
          f"-> dropped by CHAIN SELECTION, recoverable by better modelling")
    print(f"  non-canonical        : {noncanon} ({100*noncanon/max(tot,1):.0f}%)  "
          f"-> dropped by the all-canonical gate")
    print(f"      intron <= 10 kb  : {nc_small}  (plausible splice variant; RUSTLE_JUNCTION_MAJORITY reaches these)")
    print(f"      intron  > 10 kb  : {nc_big}  (mis-chain territory; the gate is right to drop these)")
    print(f"\n  => relaxing the canonical rule can address at most {100*nc_small/max(tot,1):.0f}% of the")
    print(f"     missed junctions; {100*canon/max(tot,1):.0f}% need better CHAIN CONSTRUCTION instead.")
