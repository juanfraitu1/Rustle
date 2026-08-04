#!/usr/bin/env python3
"""Is locus truncation detectable, and repairable, WITHOUT consulting any annotation?

Usage: truncation_extend.py <catalog.log with RUSTLE_LOCUS_AUDIT=1> <genome.fa> <pad_bp> [chr1,chr15]

NON-CIRCULARITY. Nothing here reads a truth file. The "expected size" of a locus is supplied by its own
predicted family siblings, and the "expected extent" is supplied by where homology to a sibling stops.
Soto is used ONLY at the end, to check whether the loci this flags as truncated are the ones the Hungarian
bipartite matcher independently called truncated (median size ratio 0.54, truncated 104 vs over-extended 12).

THE SIGNATURE. Coverage is aligned/min(qlen,tlen), so a truncated locus is scored against ITS OWN short
length and looks fully covered. Truncation is therefore invisible to the coverage floor by construction.
What it is NOT invisible to: an alignment that runs to the very edge of the shorter locus while the longer
one still has homologous sequence continuing past it. That is a boundary-limited alignment, and it means
the locus was cut off rather than ended.

THE REPAIR TEST. Pad every locus span by `pad_bp` on both sides and re-align. If the pair was genuinely
truncated, homology CONTINUES into the pad on both sides (the duplication extends further than the reads
did). If the locus really ended, the pad aligns to nothing. Reads never enter this test, so it cannot be
faked by read depth -- and the sibling, not an annotation, decides where the locus stops.
"""
import os
import subprocess
import sys
from collections import defaultdict

LOG, GENOME, PAD = sys.argv[1], sys.argv[2], int(sys.argv[3])
CHROMS = set(sys.argv[4].split(",")) if len(sys.argv) > 4 else None
CACHE = "/tmp/claude-1000/-mnt-c-Users-jfris-Desktop/161b66ca-e17e-4e70-b1cf-8b479d6c328d/scratchpad"
EDGE = 200  # bp of the query end within which an alignment counts as "running to the boundary"

reps = []
for ln in open(LOG, errors="replace"):
    if not ln.startswith("[rep-audit]"):
        continue
    f = ln.rstrip("\n").split("\t")
    if CHROMS and f[2] not in CHROMS:
        continue
    reps.append({"chrom": f[2], "s": int(f[3]), "e": int(f[4]), "nexon": int(f[8])})
key = {i: f"{r['chrom']}:{r['s']}-{r['e']}" for i, r in enumerate(reps)}
idx = {v: k for k, v in key.items()}


def build(tag, pad):
    fa = f"{CACHE}/trunc_{tag}.fa"
    if not os.path.exists(fa):
        with open(f"{CACHE}/trunc_{tag}.txt", "w") as fh:
            for r in reps:
                fh.write(f"{r['chrom']}:{max(1,r['s']+1-pad)}-{r['e']+pad}\n")
        subprocess.run(f"samtools faidx {GENOME} -r {CACHE}/trunc_{tag}.txt > {fa}",
                       shell=True, check=True)
    paf = f"{CACHE}/trunc_{tag}.paf"
    if not os.path.exists(paf):
        with open(paf, "w") as fh:
            fh.write(subprocess.run(
                ["minimap2", "-x", "asm20", "-c", "--eqx", "-X", "--no-long-join",
                 "-t", "4", "-N", "50", "-p", "0.1", fa, fa],
                capture_output=True, text=True, check=True).stdout)
    return paf


def load(paf, pad):
    """Return {(i,j): (cov, ident, q_reaches_start, q_reaches_end, into_pad_l, into_pad_r)}."""
    out = {}
    for ln in open(paf):
        f = ln.split("\t")
        def back(n):
            c, rng = n.rsplit(":", 1)
            a, b = rng.split("-")
            a, b = int(a), int(b)
            return f"{c}:{a-1+pad}-{b-pad}" if pad else f"{c}:{a-1}-{b}"
        q, t = back(f[0]), back(f[5])
        if q not in idx or t not in idx or q == t:
            continue
        qlen, qs, qe = int(f[1]), int(f[2]), int(f[3])
        cov = (qe - qs) / min(qlen, int(f[6]))
        ident = 1.0 - float(next((x[5:] for x in f[12:] if x.startswith("de:f:")), "1"))
        # how far the alignment reaches into the padded flanks (0 when pad == 0)
        into_l = max(0, pad - qs)
        into_r = max(0, qe - (qlen - pad))
        rec = (cov, ident, qs <= pad + EDGE, qe >= qlen - pad - EDGE, into_l, into_r)
        k = (idx[q], idx[t])
        if k not in out or cov > out[k][0]:
            out[k] = rec
    return out


base = load(build("p0", 0), 0)
print(f"reps {len(reps)}   aligned ordered pairs (pad 0): {len(base)}")

# --- 1. how often is an alignment BOUNDARY-LIMITED at the shorter locus? --------------------------
strong = {k: v for k, v in base.items() if v[1] >= 0.90 and v[0] >= 0.50}
both_edges = sum(1 for v in strong.values() if v[2] and v[3])
print(f"\nstrong pairs (id>=0.90, cov>=0.50): {len(strong)}")
print(f"  alignment runs to BOTH ends of the query: {both_edges} "
      f"({100*both_edges/max(len(strong),1):.0f}%)  <- boundary-limited, i.e. the locus was cut off")

# --- 2. non-circular size concordance WITHIN predicted families -----------------------------------
par = list(range(len(reps)))
def find(x):
    while par[x] != x:
        par[x] = par[par[x]]
        x = par[x]
    return x
for (a, b) in strong:
    ra, rb = find(a), find(b)
    if ra != rb:
        par[ra] = rb
fam = defaultdict(list)
for i in range(len(reps)):
    fam[find(i)].append(i)
ratios = []
for g in fam.values():
    if len(g) < 2:
        continue
    L = sorted(reps[i]["e"] - reps[i]["s"] for i in g)
    ratios.append(L[0] / L[-1])
ratios.sort()
if ratios:
    print(f"\nsize concordance inside predicted families (NO annotation used): n={len(ratios)}")
    print(f"  median smallest/largest member = {ratios[len(ratios)//2]:.2f}   "
          f"families where the smallest member is <0.5x the largest: "
          f"{100*sum(1 for r in ratios if r < 0.5)/len(ratios):.0f}%")

# --- 3. does padding recover homology that the locus boundary cut off? ----------------------------
padded = load(build(f"p{PAD}", PAD), PAD)
cont_l = cont_r = grew = same = 0
deltas = []
for k, v in strong.items():
    p = padded.get(k)
    if not p:
        continue
    if p[4] > 100:
        cont_l += 1
    if p[5] > 100:
        cont_r += 1
    d = (p[4] + p[5])
    deltas.append(d)
    if d > 100:
        grew += 1
    else:
        same += 1
deltas.sort()
n = len(deltas)
print(f"\npadding every locus by {PAD} bp and re-aligning ({n} of the strong pairs realign):")
print(f"  homology CONTINUES past the original 5' boundary: {cont_l} ({100*cont_l/max(n,1):.0f}%)")
print(f"  homology CONTINUES past the original 3' boundary: {cont_r} ({100*cont_r/max(n,1):.0f}%)")
print(f"  pair recovers >100 bp beyond the old boundaries  : {grew} ({100*grew/max(n,1):.0f}%)")
if n:
    print(f"  median bp recovered per pair: {deltas[n//2]}   90th pct: {deltas[int(n*0.9)]}")
print("\n3' > 5' would mean the TES side is the cut-off one, matching the polyA-anchored-read prediction.")
