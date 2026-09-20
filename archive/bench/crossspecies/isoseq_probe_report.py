#!/usr/bin/env python3
"""Score probe B (read-supported exons) against the 19 annotated NPIP genes."""
import sys

OUT = sys.argv[1]
bed = [l.rstrip("\n").split("\t") for l in open(f"{OUT}/hs_npip.bed")]
loci = [l.rstrip("\n").split("\t") for l in open(f"{OUT}/probeB_loci.bed")]

hit, extra = set(), []
for f in loci:
    c, a, b = f[0], int(f[1]), int(f[2])
    ov = [g[3] for g in bed if g[0] == c and int(g[1]) < b and a < int(g[2])]
    if ov:
        hit.update(ov)
    else:
        extra.append(f"{c}:{a+1}-{b}")

allg = {g[3] for g in bed}
print("=" * 66)
print("PROBE COMPARISON -- all three probes built at the SAME single seed locus")
print("=" * 66)
print(f"  A  RefSeq transcript   6,117 bp  ->   9/19 genes  (2 components)")
print(f"  B  read-supported exons          ->  {len(hit):2d}/19 genes  "
      f"({len(loci)} loci recovered)")
print(f"  C  full genomic span  25,154 bp  ->  19/19 genes  (1 component, density 1.000)")
if allg - hit:
    print(f"\n  B missed: {', '.join(sorted(allg - hit))}")
print(f"  B loci with no annotated NPIP gene: {len(extra)}/{len(loci)}")
for e in extra:
    print(f"     {e}")
