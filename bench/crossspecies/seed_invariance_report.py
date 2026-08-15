#!/usr/bin/env python3
"""P1: does every seed in the family return the same family?

Per seed, rebuild V(s) with the shipped discovery rule (block >=1kb & id>=0.80, merge at 10kb,
keep intervals with >=5kb aligned seed and >=5kb span), then compare the 19 resulting locus sets.

Locus sets are compared two ways, because interval boundaries legitimately shift with the seed:
  - by ANNOTATED GENE recovered (exact, annotation-anchored)
  - by RECIPROCAL-OVERLAP >=50% matching of the intervals themselves (annotation-free)
"""
import sys
from collections import defaultdict

OUT = sys.argv[1]
MIN_BLOCK, MIN_ID, MERGE_D, MIN_ALN, MIN_SPAN = 1000, 0.80, 10000, 5000, 5000

hits = defaultdict(list)
for line in open(f"{OUT}/allseeds_hits.paf"):
    f = line.rstrip("\n").split("\t")
    if len(f) < 12:
        continue
    q, t, ts, te, nm, bl = f[0], f[5], int(f[7]), int(f[8]), int(f[9]), int(f[10])
    if bl >= MIN_BLOCK and nm / bl >= MIN_ID:
        hits[q].append((t, ts, te, nm))


def loci_for(seed):
    by_chr = defaultdict(list)
    for t, ts, te, nm in hits.get(seed, []):
        by_chr[t].append((ts, te, nm))
    out = []
    for c, v in by_chr.items():
        v.sort()
        cs, ce, aln = None, None, 0
        for s, e, nm in v:
            if cs is None:
                cs, ce, aln = s, e, nm
            elif s <= ce + MERGE_D:
                ce = max(ce, e)
                aln += nm
            else:
                if aln >= MIN_ALN and ce - cs >= MIN_SPAN:
                    out.append((c, cs, ce))
                cs, ce, aln = s, e, nm
        if cs is not None and aln >= MIN_ALN and ce - cs >= MIN_SPAN:
            out.append((c, cs, ce))
    return sorted(out)


bed = [l.rstrip("\n").split("\t") for l in open(f"{OUT}/hs_npip.bed")]
seeds = [g[3] for g in bed]
allg = set(seeds)


def genes_of(loci):
    hit = set()
    for c, a, b in loci:
        for g in bed:
            if g[0] == c and int(g[1]) < b and a < int(g[2]):
                hit.add(g[3])
    return hit


def recip_match(A, B):
    m = 0
    for c, a, b in A:
        for c2, a2, b2 in B:
            if c != c2:
                continue
            ov = min(b, b2) - max(a, a2)
            if ov > 0 and ov >= 0.5 * (b - a) and ov >= 0.5 * (b2 - a2):
                m += 1
                break
    return m


res = {s: loci_for(s) for s in seeds}
print("\n" + "=" * 72)
print("P1 SEED-INVARIANCE -- one seed per annotated NPIP gene")
print("=" * 72)
print(f"{'seed':<10}{'loci':>6}{'NPIP genes':>12}   missed")
for s in seeds:
    g = genes_of(res[s])
    miss = sorted(allg - g)
    print(f"{s:<10}{len(res[s]):>6}{len(g):>8}/19   " + (", ".join(miss) if miss else "-"))

ref = res["NPIPB11"]
print(f"\nreference = F(NPIPB11): {len(ref)} loci")
print(f"\n{'seed':<10}{'|F(s)|':>7}{'matched':>9}{'jaccard':>9}   identical gene set?")
jac = []
for s in seeds:
    A = res[s]
    m = recip_match(A, ref)
    j = m / max(len(set(A)) + len(ref) - m, 1)
    jac.append(j)
    same = genes_of(A) == genes_of(ref)
    print(f"{s:<10}{len(A):>7}{m:>9}{j:>9.3f}   {'YES' if same else 'NO'}")

full = [s for s in seeds if genes_of(res[s]) == allg]
print(f"\n  seeds recovering ALL 19 annotated NPIP genes: {len(full)}/19")
print(f"  seeds whose gene set equals F(NPIPB11)'s:      "
      f"{sum(1 for s in seeds if genes_of(res[s])==genes_of(ref))}/19")
jac.sort()
print(f"  interval Jaccard vs F(NPIPB11): median {jac[len(jac)//2]:.3f}  "
      f"min {jac[0]:.3f}  max {jac[-1]:.3f}")
sizes = sorted(len(v) for v in res.values())
print(f"  |F(s)| across seeds: median {sizes[len(sizes)//2]}  min {sizes[0]}  max {sizes[-1]}")
verdict = ("HOLDS" if len(full) == 19 and jac[0] >= 0.9 else
           "HOLDS ON GENE SETS, intervals vary" if len(full) == 19 else "FAILS -- F must be a closure")
print(f"\n  VERDICT: P1 {verdict}")
