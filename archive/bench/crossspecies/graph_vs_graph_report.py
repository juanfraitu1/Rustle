#!/usr/bin/env python3
"""Compare the DNA graph and the spliced-RNA graph over the same node set."""
import json
import os
import sys
from collections import defaultdict

sys.path.insert(0, os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "soto"))
import rustlib as er_tier  # noqa: E402  -- bench/soto/rustlib.py is THE definition of the E_r rule

OUT = sys.argv[1]
# The PAF is now the SHIPPED SENSITIVE tier, whose identity floor is `sensitive_identity` = 0.60,
# NOT `min_identity` = 0.80 (0.80 pairs with the asm20 tier, which the shipped default skips).
ID_FLOOR, COV_FLOOR = er_tier.SENSITIVE_IDENTITY, er_tier.MIN_COVERAGE

nodes = [l.split("\t")[0] for l in open(f"{OUT}/HS.loci.fa.fai")]
nodeset = set(nodes)
nlen = {}
for l in open(f"{OUT}/gvg_nodelen.tsv"):
    r, d, x = l.rstrip("\n").split("\t")
    nlen[r] = (int(d), int(x))

bed = [l.rstrip("\n").split("\t") for l in open(f"{OUT}/hs_npip.bed")]


def gene_of(n):
    c, rng = n.rsplit(":", 1)
    a, b = rng.split("-")
    a, b = int(a) - 1, int(b)
    g = [x[3] for x in bed if x[0] == c and int(x[1]) < b and a < int(x[2])]
    return ",".join(g) if g else "-"


# ⚠ TIER/RULE CORRECTION 2026-08-10 (defects B1 + M1). This module re-implemented the edge rule and
# got two things wrong: identity was nmatch/blocklen (the binary uses 1-de, with nmatch/blocklen only
# as a FALLBACK), and coverage was `(qe-qs)/min(ql,tl)` -- a QUERY-axis numerator over a possibly
# TARGET denominator, which is not a fraction and exceeded 1.0 on 8.2% of records. Both now come from
# bench/soto/rustlib.py, the single Python mirror of `nucleotide_edges_scored`.
def best(paf):
    """per unordered pair, the highest SHIPPED-RULE coverage among records clearing identity"""
    return er_tier.paf_pairs(paf, min_identity=ID_FLOOR, min_coverage=COV_FLOOR, denom="min")


def diagnose(paf):
    """UNFILTERED per-pair view for the failure-mode breakdown ONLY.

    ⚠ DEFECT FOUND 2026-08-11: `why()` used to be handed `best()`, whose floors are applied PER
    RECORD *inside* paf_pairs, so a pair that fails either floor is simply ABSENT from the dict.
    Every failing pair therefore fell into the `v is None` branch and the breakdown printed
    "no alignment record at all: <all of them>" whatever the truth was -- and the two comparison
    branches were dead code (they also indexed a dict as a tuple, so they would have raised).
    Measured on the 27-node shipped-tier run: it printed 91/91 "no record", while all 351 pairs
    have records and the true mode is 91/91 COVERAGE misses at identity 0.971-1.000.

    Returns {(a,b): (max identity over all records, best M1-coverage among records with id>=floor
    or -1 if none, n_records)}.
    """
    out = {}
    for ln in open(paf):
        f = ln.rstrip("\n").split("\t")
        if len(f) < 12:
            continue
        q, t = f[0], f[5]
        if q == t or q not in nodeset or t not in nodeset:
            continue
        ql, qs, qe = int(f[1]), int(f[2]), int(f[3])
        tl, ts, te = int(f[6]), int(f[7]), int(f[8])
        de = None
        for x in f[12:]:
            if x.startswith("de:f:"):
                de = float(x[5:])
                break
        idt = (1.0 - de) if de is not None else int(f[9]) / max(int(f[10]), 1)
        # M1: the numerator's axis follows the denominator (same rule as rustlib.paf_pairs).
        cov = (qe - qs) / max(ql, 1) if ql <= tl else (te - ts) / max(tl, 1)
        k = (q, t) if q < t else (t, q)
        d = out.setdefault(k, [0.0, -1.0, 0])
        d[2] += 1
        d[0] = max(d[0], idt)
        if idt >= ID_FLOOR:
            d[1] = max(d[1], cov)
    return out


def edges(b):
    return set(b)


def comps(ns, E):
    p = {n: n for n in ns}

    def f(x):
        while p[x] != x:
            p[x] = p[p[x]]
            x = p[x]
        return x
    for e in E:
        a, b_ = tuple(e)
        ra, rb = f(a), f(b_)
        if ra != rb:
            p[ra] = rb
    g = defaultdict(list)
    for n in ns:
        g[f(n)].append(n)
    return sorted(g.values(), key=len, reverse=True)


bd, br = best(f"{OUT}/gvg_dna.paf"), best(f"{OUT}/gvg_rna.paf")
Ed, Er = edges(bd), edges(br)
N = len(nodes)
POSS = N * (N - 1) // 2

print("=" * 74)
print(f"DNA GRAPH vs SPLICED-RNA GRAPH   |V| = {N} (identical node set), {POSS} possible edges")
print("=" * 74)
dl = sorted(v[0] for v in nlen.values())
rl = sorted(v[1] for v in nlen.values())
print(f"  node sequence length  DNA median {dl[len(dl)//2]:,} bp   "
      f"RNA median {rl[len(rl)//2]:,} bp   ratio {rl[len(rl)//2]/dl[len(dl)//2]:.2f}")
empty = [n for n in nodes if nlen.get(n, (0, 0))[1] == 0]
if empty:
    print(f"  ⚠ nodes with NO read support (absent from the RNA graph): {len(empty)}")
    for n in empty:
        print(f"      {n}  [{gene_of(n)}]")

for lbl, E in (("DNA", Ed), ("RNA", Er)):
    cs = comps(nodes, E)
    d = 2 * len(E) / (N * (N - 1))
    print(f"\n  {lbl:<4} edges {len(E):>4}/{POSS}   density {d:.3f}   "
          f"components {len(cs)}   largest {len(cs[0])}   singletons "
          f"{sum(1 for c in cs if len(c)==1)}")

inter, only_d, only_r = Ed & Er, Ed - Er, Er - Ed
uni = Ed | Er
print(f"\n  shared      {len(inter):>4}")
print(f"  DNA only    {len(only_d):>4}")
print(f"  RNA only    {len(only_r):>4}")
print(f"  Jaccard     {len(inter)/max(len(uni),1):>7.3f}")
print(f"  RNA subset of DNA? {'YES' if not only_r else f'NO ({len(only_r)} RNA-only edges)'}")
print(f"  DNA subset of RNA? {'YES' if not only_d else f'NO ({len(only_d)} DNA-only edges)'}")

# WHY do the disagreements happen -- identity or coverage?
def why(pairs, other_diag, lbl):
    if not pairs:
        return
    ni = nc = nn = 0
    covs = []
    for k in pairs:
        v = other_diag.get(k)
        if v is None:
            nn += 1
        elif v[1] < 0:          # records exist, none of them clears the identity floor
            ni += 1
        else:                   # identity clears, coverage does not
            nc += 1
            covs.append(v[1])
    print(f"\n  {lbl}: {len(pairs)} edges. In the other graph they fail because —")
    print(f"     no alignment record at all      : {nn}")
    print(f"     records exist, ALL id < {ID_FLOOR}    : {ni}")
    print(f"     identity ok, coverage < {COV_FLOOR}    : {nc}")
    if covs:
        covs.sort()
        print(f"     their best id-passing coverage: {covs[0]:.4f}–{covs[-1]:.4f} "
              f"(median {covs[len(covs)//2]:.4f})")


dd, dr = diagnose(f"{OUT}/gvg_dna.paf"), diagnose(f"{OUT}/gvg_rna.paf")
why(only_d, dr, "DNA-only")
why(only_r, dd, "RNA-only")

deg_d, deg_r = defaultdict(int), defaultdict(int)
for e in Ed:
    for n in e:
        deg_d[n] += 1
for e in Er:
    for n in e:
        deg_r[n] += 1
print(f"\n  per-node degree (max {N-1}):")
print(f"  {'locus':<26}{'gene':<12}{'DNA':>5}{'RNA':>5}{'Δ':>6}{'rna_bp':>9}")
for n in sorted(nodes, key=lambda x: deg_r[x] - deg_d[x]):
    print(f"  {n:<26}{gene_of(n)[:11]:<12}{deg_d[n]:>5}{deg_r[n]:>5}"
          f"{deg_r[n]-deg_d[n]:>6}{nlen.get(n,(0,0))[1]:>9,}")

json.dump({"nodes": nodes, "gene": {n: gene_of(n) for n in nodes},
           "dna": [sorted(e) for e in Ed], "rna": [sorted(e) for e in Er],
           "len": {n: nlen.get(n, (0, 0)) for n in nodes}},
          open(f"{OUT}/gvg_graphs.json", "w"))
print(f"\n  wrote {OUT}/gvg_graphs.json")
