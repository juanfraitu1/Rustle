#!/usr/bin/env python3
"""Does the definition hold on families defined INDEPENDENTLY of anything we did today?

TRUTH: `bench/soto/80_fams.gene_preferred.bed` -- Soto's curated GENE|ID_NNN membership, 362 members
in 80 families on CHM13. Not our output, not name prefixes. Using the subfamily structure our own run
discovered would be circular; this is not.

DEFINITION v2 -- the two repairs today's failures demanded, applied BEFORE looking at the result:
  (1) RELATIVE discovery floor. v1 used an absolute 5 kb of aligned seed, which excludes 14/18 PCDHB
      members and returns ZERO loci for TUBA1A/TUBA1B. Floor is now max(1000, 0.30 * seed length).
      Third absolute threshold to fail today; scale to the seed, never to a constant.
  (2) CLOSURE, not a single pass. F(s) = fixed point of seed-and-extend: discover from s, then
      re-discover from everything found, until nothing new appears. v1's single pass failed P1 on
      all five families tested, while the union over seeds recovered 4/5 -- the closure is that union
      reached from ONE seed instead of from all of them.

Edge rule unchanged and as shipped: asm20 id>=0.80 UNION sensitive(-k11 -w5) id>=0.60, coverage>=0.50
AGGREGATED over records.

REPORTED PER FAMILY:
  P1  do all member seeds reach the same closure?      (the anti-arbitrariness property)
  P2  density / components of the resulting graph      (is it a quasi-clique, with what margin)
  P3  how much of the CURATED membership is recovered  (non-circular recall)

⚠ Both genome passes are batched -- one minimap2 index build per iteration, not one per seed.
"""
import os
import subprocess
import sys
from collections import defaultdict
sys.path.insert(0, os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "soto"))
import rustlib as er_tier  # noqa: E402  -- bench/soto/rustlib.py is THE definition of the E_r rule


W = "/home/juanfra/winloci_scratch/seedfam/cur"
HS = "/mnt/c/Users/jfris/Desktop/Reference/chm13v2.0.fa"
TRUTH = "/mnt/c/Users/jfris/Desktop/Rustle/bench/soto/80_fams.gene_preferred.bed"
T = "4"
FRAC_FLOOR, ABS_MIN, MIN_BLOCK, MERGE_D = 0.30, 1000, 500, 10000
FAMS = sys.argv[1].split(",") if len(sys.argv) > 1 else [
    "ID_131", "ID_431", "ID_391", "ID_395", "ID_400", "ID_116", "ID_8", "ID_468", "ID_154"]
os.makedirs(W, exist_ok=True)


def sh(cmd, out=None):
    if out:
        with open(out, "w") as fh:
            return subprocess.run(cmd, stdout=fh, stderr=subprocess.DEVNULL).returncode
    return subprocess.run(cmd, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL).returncode


members = defaultdict(list)
for line in open(TRUTH):
    f = line.rstrip("\n").split("\t")
    if len(f) < 4 or "|" not in f[3]:
        continue
    g, fid = f[3].split("|", 1)
    if fid in FAMS:
        members[fid].append((f[0], int(f[1]), int(f[2]), g))
for k in FAMS:
    members[k].sort()
    print(f"{k:<8} {len(members[k]):>3} curated members: "
          + ", ".join(x[3] for x in members[k][:6]) + (" ..." if len(members[k]) > 6 else ""))


def write_fa(path, items):
    with open(path, "w") as fh:
        for name, c, s, e in items:
            out = subprocess.run(["samtools", "faidx", HS, f"{c}:{s+1}-{e}"],
                                 capture_output=True, text=True).stdout.splitlines()
            if len(out) > 1:
                fh.write(f">{name}\n" + "\n".join(out[1:]) + "\n")


def discover(paf, seedlen):
    """PAF -> {seed: [loci]} using the RELATIVE floor."""
    hits = defaultdict(list)
    for line in open(paf):
        f = line.rstrip("\n").split("\t")
        if len(f) < 12:
            continue
        q, t, ts, te, nm, bl = f[0], f[5], int(f[7]), int(f[8]), int(f[9]), int(f[10])
        if bl >= MIN_BLOCK and nm / bl >= 0.80:
            hits[q].append((t, ts, te, nm))
    out = {}
    for q, v in hits.items():
        floor = max(ABS_MIN, FRAC_FLOOR * seedlen.get(q, 0))
        by = defaultdict(list)
        for t, ts, te, nm in v:
            by[t].append((ts, te, nm))
        L = []
        for c, w in by.items():
            w.sort()
            cs = ce = None
            aln = 0
            for s, e, nm in w:
                if cs is None:
                    cs, ce, aln = s, e, nm
                elif s <= ce + MERGE_D:
                    ce = max(ce, e)
                    aln += nm
                else:
                    if aln >= floor and ce - cs >= floor:
                        L.append((c, cs, ce))
                    cs, ce, aln = s, e, nm
            if cs is not None and aln >= floor and ce - cs >= floor:
                L.append((c, cs, ce))
        out[q] = sorted(L)
    return out


# ---------- iteration 1: every curated member as its own seed -----------------
items, seedlen = [], {}
for fid in FAMS:
    for c, s, e, g in members[fid]:
        n = f"{fid}|{g}"
        items.append((n, c, s, e))
        seedlen[n] = e - s
write_fa(f"{W}/it1.fa", items)
print(f"\niteration 1: {len(items)} seeds -> one genome pass")
sh(["minimap2", "-x", "asm20", "-c", "--eqx", "-N", "200", "-p", "0.02", "-I", "2G", "-t", T,
    HS, f"{W}/it1.fa"], out=f"{W}/it1.paf")
L1 = discover(f"{W}/it1.paf", seedlen)

# ---------- iteration 2: re-seed from everything found ------------------------
found = sorted({x for v in L1.values() for x in v})
items2, seedlen2 = [], {}
for c, s, e in found:
    n = f"{c}:{s+1}-{e}"
    items2.append((n, c, s, e))
    seedlen2[n] = e - s
write_fa(f"{W}/it2.fa", items2)
print(f"iteration 2: {len(items2)} discovered loci re-seeded -> one genome pass")
sh(["minimap2", "-x", "asm20", "-c", "--eqx", "-N", "200", "-p", "0.02", "-I", "2G", "-t", T,
    HS, f"{W}/it2.fa"], out=f"{W}/it2.paf")
L2 = discover(f"{W}/it2.paf", seedlen2)


def closure(seed):
    """F(seed) = fixed point: loci from seed, plus loci from each of those."""
    cur = set(L1.get(seed, []))
    for _ in range(3):
        nxt = set(cur)
        for c, s, e in cur:
            nxt |= set(L2.get(f"{c}:{s+1}-{e}", []))
        if nxt == cur:
            break
        cur = nxt
    return cur


def ulen(iv):
    iv = sorted(iv)
    tot = 0
    cs = ce = None
    for s, e in iv:
        if cs is None:
            cs, ce = s, e
        elif s > ce:
            tot += ce - cs
            cs, ce = s, e
        else:
            ce = max(ce, e)
    return tot + (ce - cs if cs is not None else 0)


# ⚠ TIER/RULE CORRECTION 2026-08-10 (defects B1 + M1 + B2). This function re-implemented the edge
# rule and diverged from the binary on FOUR axes: the asm20 UNION sensitive tiering (shipped default
# is sensitive-ONLY), identity = nmatch/blocklen (binary uses 1-de first), the UNION-of-collinear-
# records coverage (= RUSTLE_ER_SUM_COVERAGE, DEFAULT OFF, and it costs cross-family precision
# 0.9038), and the M1 query-axis-over-target-denominator coverage. All of it now comes from
# bench/soto/rustlib.py, the single Python mirror of `nucleotide_edges_scored`.
def er(prefix, ns):
    return er_tier.er_edges([(f"{prefix}.paf", er_tier.SENSITIVE_IDENTITY)], nodes=set(ns),
                           min_coverage=er_tier.MIN_COVERAGE, denom="min")


def comps(ns, E):
    p = {n: n for n in ns}

    def fi(x):
        while p[x] != x:
            p[x] = p[p[x]]
            x = p[x]
        return x
    for e in E:
        a, b = tuple(e)
        ra, rb = fi(a), fi(b)
        if ra != rb:
            p[ra] = rb
    g = defaultdict(list)
    for n in ns:
        g[fi(n)].append(n)
    return sorted(g.values(), key=len, reverse=True)


print("\n" + "=" * 88)
print("P1 CLOSURE-INVARIANCE + P3 CURATED RECALL")
print("=" * 88)
print(f"{'family':<9}{'n':>4}{'|F| min-max':>14}{'all seeds same F?':>19}"
      f"{'curated recovered':>20}{'seeds w/ full set':>19}")
keep = {}
for fid in FAMS:
    cl = {g: closure(f"{fid}|{g}") for _, _, _, g in members[fid]}
    sizes = sorted(len(v) for v in cl.values())
    same = len({frozenset(v) for v in cl.values()}) == 1

    def rec(L):
        return {g for c, a, b in L for cc, s, e, g in members[fid]
                if cc == c and s < b and a < e}
    got = {g: rec(L) for g, L in cl.items()}
    full = sum(1 for v in got.values() if len(v) == len(members[fid]))
    union = set().union(*got.values()) if got else set()
    keep[fid] = cl
    print(f"{fid:<9}{len(members[fid]):>4}{f'{sizes[0]}-{sizes[-1]}':>14}"
          f"{('YES' if same else 'NO'):>19}{f'{len(union)}/{len(members[fid])}':>20}"
          f"{f'{full}/{len(cl)}':>19}")

print("\n" + "=" * 88)
print("P2 GRAPH STRUCTURE of the closure (shipped rule, aggregated coverage)")
print("=" * 88)
print(f"{'family':<9}{'|V|':>5}{'edges':>7}{'density':>9}{'comps':>7}{'largest':>9}{'margin':>9}")
for fid in FAMS:
    U = sorted(set().union(*keep[fid].values())) if keep[fid] else []
    if len(U) < 2:
        print(f"{fid:<9}{len(U):>5}   (too few loci)")
        continue
    with open(f"{W}/{fid}.regions", "w") as fh:
        for c, a, b in U:
            fh.write(f"{c}:{a+1}-{b}\n")
    sh(["samtools", "faidx", HS, "-r", f"{W}/{fid}.regions"], out=f"{W}/{fid}.fa")
    sh(["samtools", "faidx", f"{W}/{fid}.fa"])
    ns = [l.split("\t")[0] for l in open(f"{W}/{fid}.fa.fai")]
    er_tier.ava(f"{W}/{fid}.fa", f"{W}/{fid}.paf", threads=int(T))
    E = er(f"{W}/{fid}", set(ns))
    cs = comps(ns, E)
    big = cs[0]
    d = (2 * sum(1 for x in E if set(x) <= set(big)) / (len(big) * (len(big) - 1))
         if len(big) > 1 else float("nan"))
    print(f"{fid:<9}{len(ns):>5}{len(E):>7}{d:>9.3f}{len(cs):>7}{len(big):>9}{d/0.40:>8.2f}x")
print("\nDONE")
