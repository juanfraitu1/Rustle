#!/usr/bin/env python3
"""How many IsoSeq reads at NPIP loci place alignments on MORE THAN ONE NPIP gene?

⚠ This is a CERTIFICATE, not a family definition. Read-conflict (E_c) is the ambiguity oracle and
  belongs to O2 (assignment). O1 defines families by sequence homology (E_r) alone. This number is
  evidence to SHOW that the components were tried against each other -- it must never become the
  edge rule, or O1 and O2 stop being independent.

Counts primary AND secondary alignments (the BAM was built with --secondary=yes -N 50), grouped by
read name, then asks how many DISTINCT NPIP genes each read touches.
"""
import subprocess
import sys
from collections import defaultdict

bam, bed = sys.argv[1], sys.argv[2]
genes = []
for line in open(bed):
    f = line.rstrip("\n").split("\t")
    if len(f) >= 4:
        genes.append((f[0], int(f[1]), int(f[2]), f[3]))

touch = defaultdict(set)
mapq0 = defaultdict(int)
for c, s, e, g in genes:
    # flag 2308 = unmapped(4) + supplementary(2048) + qcfail/dup(256 is SECONDARY, kept on purpose)
    # here we KEEP secondary (256) because multimapping is the whole point; drop 4+2048 only.
    p = subprocess.run(["samtools", "view", "-F", "2052", bam, f"{c}:{s+1}-{e}"],
                       capture_output=True, text=True)
    for line in p.stdout.splitlines():
        f = line.split("\t", 5)
        if len(f) < 5:
            continue
        touch[f[0]].add(g)
        if f[4] == "0":
            mapq0[f[0]] += 1

n = len(touch)
multi = {r: v for r, v in touch.items() if len(v) > 1}
print(f"  reads with >=1 alignment on an NPIP gene : {n}")
print(f"  reads touching >1 NPIP gene              : {len(multi)}  ({len(multi)/max(n,1):.1%})")
if multi:
    sizes = sorted(len(v) for v in multi.values())
    print(f"  genes touched per multi-read: median {sizes[len(sizes)//2]}  max {sizes[-1]}")
print(f"  reads with a MAPQ=0 alignment            : {len(mapq0)}  ({len(mapq0)/max(n,1):.1%})")

pair = defaultdict(int)
for v in multi.values():
    vs = sorted(v)
    for i in range(len(vs)):
        for j in range(i + 1, len(vs)):
            pair[(vs[i], vs[j])] += 1
allg = [g[3] for g in genes]
print(f"\n  gene PAIRS linked by >=1 shared read: {len(pair)} of "
      f"{len(allg)*(len(allg)-1)//2} possible")
seen = set()
for (a, b) in pair:
    seen.add(a)
    seen.add(b)
print(f"  NPIP genes participating in >=1 shared read: {len(seen)}/{len(allg)}")
if set(allg) - seen:
    print(f"  genes with NO shared read: {', '.join(sorted(set(allg)-seen))}")
print("\n  top linked pairs:")
for (a, b), c in sorted(pair.items(), key=lambda x: -x[1])[:12]:
    print(f"    {a:<9} {b:<9} {c} reads")
