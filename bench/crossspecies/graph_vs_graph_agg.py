#!/usr/bin/env python3
"""DNA vs spliced-RNA graph -- SINGLE-RECORD rule vs AGGREGATED rule, both substrates.

WHY: concatenating exons manufactures artificial junctions. A single alignment record cannot span a
concatenation boundary where two copies differ in exon structure, so the shipped "identity AND coverage
ON ONE PAF RECORD" rule penalises the spliced substrate specifically. Evidence from the first run: the
RNA all-vs-all produced 2,662 records against DNA's 1,966 over the SAME 27 nodes -- more records, each
shorter. The 0.547 RNA density may therefore measure concatenation, not biology.

FIX: aggregate per pair -- coverage = |union of aligned query intervals| / min(qlen,tlen),
identity = sum(nmatch)/sum(blocklen) over those records. Applied IDENTICALLY to both substrates, so it
cannot flatter one of them.

Both rules are reported side by side. If DNA barely moves and RNA jumps, the single-record rule was the
artifact; if both move together, the substrate difference is real.
"""
import sys
from collections import defaultdict

OUT = sys.argv[1]
ID_FLOOR, COV_FLOOR = 0.80, 0.50

nodes = [l.split("\t")[0] for l in open(f"{OUT}/HS.loci.fa.fai")]
nodeset = set(nodes)
bed = [l.rstrip("\n").split("\t") for l in open(f"{OUT}/hs_npip.bed")]


def gene_of(n):
    c, rng = n.rsplit(":", 1)
    a, b = rng.split("-")
    a, b = int(a) - 1, int(b)
    g = [x[3] for x in bed if x[0] == c and int(x[1]) < b and a < int(x[2])]
    return ",".join(g) if g else "-"


def load(paf):
    recs = defaultdict(list)
    for line in open(paf):
        f = line.rstrip("\n").split("\t")
        if len(f) < 12:
            continue
        q, ql, qs, qe, t, tl = f[0], int(f[1]), int(f[2]), int(f[3]), f[5], int(f[6])
        nm, bl = int(f[9]), int(f[10])
        if q == t or bl == 0 or q not in nodeset or t not in nodeset:
            continue
        recs[(q, t)].append((qs, qe, nm, bl, ql, tl))
    return recs


def union_len(iv):
    iv = sorted(iv)
    tot = cs = ce = 0
    started = False
    for s, e in iv:
        if not started:
            cs, ce, started = s, e, True
        elif s > ce:
            tot += ce - cs
            cs, ce = s, e
        else:
            ce = max(ce, e)
    if started:
        tot += ce - cs
    return tot


def edges_single(recs):
    E = set()
    for (q, t), v in recs.items():
        for qs, qe, nm, bl, ql, tl in v:
            if nm / bl >= ID_FLOOR and (qe - qs) / max(min(ql, tl), 1) >= COV_FLOOR:
                E.add(frozenset((q, t)))
                break
    return E


def edges_agg(recs):
    E = set()
    for (q, t), v in recs.items():
        keep = [r for r in v if r[2] / r[3] >= ID_FLOOR]
        if not keep:
            continue
        ql, tl = keep[0][4], keep[0][5]
        cov = union_len([(r[0], r[1]) for r in keep]) / max(min(ql, tl), 1)
        idy = sum(r[2] for r in keep) / sum(r[3] for r in keep)
        if idy >= ID_FLOOR and cov >= COV_FLOOR:
            E.add(frozenset((q, t)))
    return E


def comps(ns, E):
    p = {n: n for n in ns}

    def f(x):
        while p[x] != x:
            p[x] = p[p[x]]
            x = p[x]
        return x
    for e in E:
        a, b = tuple(e)
        ra, rb = f(a), f(b)
        if ra != rb:
            p[ra] = rb
    g = defaultdict(list)
    for n in ns:
        g[f(n)].append(n)
    return sorted(g.values(), key=len, reverse=True)


N = len(nodes)
POSS = N * (N - 1) // 2
rd, rr = load(f"{OUT}/gvg_dna.paf"), load(f"{OUT}/gvg_rna.paf")

print("=" * 74)
print(f"SINGLE-RECORD vs AGGREGATED rule   |V|={N}, {POSS} possible edges")
print("=" * 74)
fr_d = sum(len(v) for v in rd.values()) / max(len(rd), 1)
fr_r = sum(len(v) for v in rr.values()) / max(len(rr), 1)
print(f"  records per aligned pair:  DNA {fr_d:.1f}   RNA {fr_r:.1f}   "
      f"(fragmentation ratio {fr_r/fr_d:.2f}x)")

G = {}
for lbl, recs in (("DNA", rd), ("RNA", rr)):
    for rule, fn in (("single", edges_single), ("agg", edges_agg)):
        E = fn(recs)
        G[(lbl, rule)] = E
        cs = comps(nodes, E)
        print(f"  {lbl}/{rule:<6} edges {len(E):>4}/{POSS}   density {2*len(E)/(N*(N-1)):.3f}   "
              f"components {len(cs)}   largest {cs[0].__len__()}")

for rule in ("single", "agg"):
    Ed, Er = G[("DNA", rule)], G[("RNA", rule)]
    inter, od, orr = Ed & Er, Ed - Er, Er - Ed
    print(f"\n  --- {rule} rule ---")
    print(f"     shared {len(inter)}   DNA-only {len(od)}   RNA-only {len(orr)}   "
          f"Jaccard {len(inter)/max(len(Ed|Er),1):.3f}")
    print(f"     RNA subset of DNA? {'YES' if not orr else f'NO ({len(orr)})'}")

Es, Ea = G[("RNA", "single")], G[("RNA", "agg")]
print(f"\n  RNA edges recovered by aggregating: {len(Ea - Es)}  "
      f"({len(Es)} -> {len(Ea)})")
Ds, Da = G[("DNA", "single")], G[("DNA", "agg")]
print(f"  DNA edges recovered by aggregating: {len(Da - Ds)}  ({len(Ds)} -> {len(Da)})")
print("\n  READ: if RNA moves a lot and DNA barely moves, the single-record rule was the artifact.")

Ea_, Da_ = G[("RNA", "agg")], G[("DNA", "agg")]
deg_d, deg_r = defaultdict(int), defaultdict(int)
for e in Da_:
    for n in e:
        deg_d[n] += 1
for e in Ea_:
    for n in e:
        deg_r[n] += 1
miss = [n for n in nodes if deg_r[n] < deg_d[n]]
if miss:
    print(f"\n  nodes still short of full degree under the aggregated rule (max {N-1}):")
    for n in sorted(miss, key=lambda x: deg_r[x] - deg_d[x])[:10]:
        print(f"    {n:<26}{gene_of(n)[:11]:<12} DNA {deg_d[n]:>3}  RNA {deg_r[n]:>3}")
