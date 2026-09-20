#!/usr/bin/env python3
"""Compare the genomic-span seed against the mRNA seed (bench/crossspecies/seed_cds.sh)."""
import sys
from collections import defaultdict

OUT = sys.argv[1]
ID_FLOOR, COV_FLOOR, GAMMA = 0.80, 0.50, 0.40


def parse_region(n):
    c, r = n.rsplit(":", 1)
    a, b = r.split("-")
    return c, int(a) - 1, int(b)


def build_edges(paf, loci):
    edges = set()
    for line in open(paf):
        f = line.rstrip("\n").split("\t")
        if len(f) < 12:
            continue
        q, ql, qs, qe, t, tl = f[0], int(f[1]), int(f[2]), int(f[3]), f[5], int(f[6])
        nm, bl = int(f[9]), int(f[10])
        if q == t or bl == 0 or nm / bl < ID_FLOOR:
            continue
        if (qe - qs) / max(min(ql, tl), 1) >= COV_FLOOR:
            edges.add(frozenset((q, t)))
    return {e for e in edges if set(e) <= set(loci)}


def comps(nodes, edges):
    p = {n: n for n in nodes}

    def f(x):
        while p[x] != x:
            p[x] = p[p[x]]
            x = p[x]
        return x
    for e in edges:
        a, b = tuple(e)
        ra, rb = f(a), f(b)
        if ra != rb:
            p[ra] = rb
    g = defaultdict(list)
    for n in nodes:
        g[f(n)].append(n)
    return sorted(g.values(), key=len, reverse=True)


def load_bed(path):
    out = []
    for line in open(path):
        x = line.rstrip("\n").split("\t")
        if len(x) >= 3:
            out.append((x[0], int(x[1]), int(x[2]), x[3] if len(x) > 3 else ""))
    return out


for tag in ("GGO", "HS"):
    loci = [l.split("\t")[0] for l in open(f"{OUT}/{tag}.cdsloci.fa.fai")]
    e = build_edges(f"{OUT}/{tag}.cds_ava.paf", loci)
    cs = comps(loci, e)
    big = cs[0]
    n = len(big)
    ein = sum(1 for x in e if set(x) <= set(big))
    d = 2 * ein / (n * (n - 1)) if n > 1 else float("nan")
    span = [parse_region(x) for x in loci]
    print(f"\n{'='*66}\n{tag}: mRNA-seeded  ->  {len(loci)} loci\n{'='*66}")
    print(f"  locus span median {sorted(b-a for _,a,b in span)[len(span)//2]:,} bp   "
          f"total {sum(b-a for _,a,b in span)/1e6:.2f} Mb")
    print(f"  edges {len(e)}  components {len(cs)}  largest {n}/{len(loci)}  "
          f"density {d:.3f}  (gamma={GAMMA}, margin {d/GAMMA:.2f}x)")
    prev = sum(1 for _ in open(f"{OUT}/{tag}.loci.bed"))
    print(f"  genomic-span seed gave {prev} loci  ->  mRNA seed gives {len(loci)}")

# human: does the mRNA seed still recover all 19 annotated NPIP genes?
hs = [l.split("\t")[0] for l in open(f"{OUT}/HS.cdsloci.fa.fai")]
bed = load_bed(f"{OUT}/hs_npip.bed")
hit = set()
extra = []
for L in hs:
    c, a, b = parse_region(L)
    ov = [g for g in bed if g[0] == c and g[1] < b and a < g[2]]
    if ov:
        hit.update(g[3] for g in ov)
    else:
        extra.append(L)
allg = {g[3] for g in bed}
print(f"\n  HUMAN CHECK -- mRNA seed vs the 19 annotated NPIP genes:")
print(f"     recovered {len(hit)}/{len(allg)}"
      + (f"   missed: {', '.join(sorted(allg-hit))}" if allg - hit else ""))
print(f"     loci with no annotated NPIP gene: {len(extra)}/{len(hs)}")
for L in extra:
    print(f"        {L}")
