#!/usr/bin/env python3
"""Is the read-supported exon BLOCK SET a better locus representation than the genomic span?

2x2: SUBSTRATE (genomic span | read-supported exon blocks) x DENOMINATOR (min | max).

WHY THIS SHAPE. Union-of-exons was already measured as a representative and LOST 20 recall points,
while today the same object worked as a discovery PROBE (19/19 vs the annotated transcript's 9/19).
The proposed explanation is the denominator: coverage on min(|u|,|v|) punishes a longer rep by
inflating its own denominator, so a richer representation is penalised for being richer. If that is
right, blocks should lose under min() and win (or stop losing) under max(). If blocks lose under BOTH,
the explanation is wrong and the block set is simply a worse representative.

⚠ SCORED AT THE FAMILY LEVEL, against Soto's curated membership. Changing the representative changes
  what a NODE IS, and node-level metrics have failed to detect that three times end-to-end.
⚠ NODE SET IS IDENTICAL ACROSS ALL FOUR ARMS -- discovered once from genomic seeds. Only the sequence
  each node contributes changes. Otherwise substrate and node set are confounded.
⚠ PRECISION IS A LOWER BOUND: PPIAL4G and PPIAL4H are real family members absent from Soto's list, so
  some "extra" loci are discoveries, not errors. Reported as such, never as a clean precision.
⚠ Discovery floor is max(300, 0.30*seed) -- the 1000 bp clamp silently zeroed PPIAL4 (760 bp members)
  and would zero 10 of 83 curated families.
"""
import os
import subprocess
import sys
from collections import defaultdict
sys.path.insert(0, os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "soto"))
import rustlib as er_tier  # noqa: E402  -- bench/soto/rustlib.py is THE definition of the E_r rule


W = "/home/juanfra/winloci_scratch/seedfam/rep2x2"
CUR = "/home/juanfra/winloci_scratch/seedfam/cur"
HS = "/mnt/c/Users/jfris/Desktop/Reference/chm13v2.0.fa"
BAM = "/mnt/linuxdisk/home/juanfraitu/winloci_data/A119b.t2t.bam"
TRUTH = "/mnt/c/Users/jfris/Desktop/Rustle/bench/soto/80_fams.gene_preferred.bed"
HERE = os.path.dirname(os.path.abspath(__file__))
T = "4"
MIN_BLOCK, MERGE_D, FRAC, ABS = 500, 10000, 0.30, 300
MIN_READS, MAX_LOCI = 3, 120
FAMS = sys.argv[1].split(",") if len(sys.argv) > 1 else [
    "ID_431", "ID_131", "ID_391", "ID_395", "ID_400", "ID_116", "ID_8", "ID_468", "ID_154"]
os.makedirs(W, exist_ok=True)


def sh(cmd, out=None):
    if out:
        with open(out, "w") as fh:
            return subprocess.run(cmd, stdout=fh, stderr=subprocess.DEVNULL).returncode
    return subprocess.run(cmd, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL).returncode


members = defaultdict(list)
for line in open(TRUTH):
    f = line.rstrip("\n").split("\t")
    if len(f) >= 4 and "|" in f[3]:
        g, fid = f[3].split("|", 1)
        if fid in FAMS:
            members[fid].append((f[0], int(f[1]), int(f[2]), g))

seedlen, name, n = {}, None, 0
for l in open(f"{CUR}/it1.fa"):
    if l[0] == ">":
        if name:
            seedlen[name] = n
        name, n = l[1:].strip(), 0
    else:
        n += len(l.strip())
if name:
    seedlen[name] = n

hits = defaultdict(list)
for line in open(f"{CUR}/it1.paf"):
    f = line.rstrip("\n").split("\t")
    if len(f) < 12:
        continue
    q, t, ts, te, nm, bl = f[0], f[5], int(f[7]), int(f[8]), int(f[9]), int(f[10])
    if bl >= MIN_BLOCK and nm / bl >= 0.80:
        hits[q].append((t, ts, te, nm))


def discover(q):
    floor = max(ABS, FRAC * seedlen.get(q, 0))
    by = defaultdict(list)
    for t, ts, te, nm in hits.get(q, []):
        by[t].append((ts, te, nm))
    L = []
    for c, v in by.items():
        v.sort()
        cs = ce = None
        aln = 0
        for s, e, nm in v:
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
    return L


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
def er(prefix, ns, denom):
    """`denom` selects the coverage DENOMINATOR ablation this script exists to run:
    "min" = coverage-of-shorter (shipped), "max" = coverage-of-longer
    (RUSTLE_ER_COVERAGE_LONGER). Under M1 the numerator's axis follows the denominator in
    both arms, so the ablation now compares two well-formed fractions instead of two
    mixed-axis quantities."""
    return er_tier.er_edges([(f"{prefix}.paf", er_tier.SENSITIVE_IDENTITY)], nodes=set(ns),
                           min_coverage=er_tier.MIN_COVERAGE, denom=denom)


def comps(ns, E):
    p = {x: x for x in ns}

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
    for x in ns:
        g[fi(x)].append(x)
    return sorted(g.values(), key=len, reverse=True)


results = defaultdict(dict)
for fid in FAMS:
    U = sorted({x for _, _, _, g in members[fid] for x in discover(f"{fid}|{g}")})
    capped = len(U) > MAX_LOCI
    if capped:
        U = U[:MAX_LOCI]
    regs = [f"{c}:{a+1}-{b}" for c, a, b in U]
    with open(f"{W}/{fid}.regions", "w") as fh:
        fh.write("\n".join(regs) + "\n")
    sh(["samtools", "faidx", HS, "-r", f"{W}/{fid}.regions"], out=f"{W}/{fid}_dna.fa")

    # RNA node = read-supported exon blocks inside the SAME interval
    with open(f"{W}/{fid}_rna.fa", "w") as fh:
        for (c, a, b), reg in zip(U, regs):
            p = subprocess.run(["samtools", "view", "-F", "2308", BAM, f"{c}:{a+1}-{b}"],
                               capture_output=True, text=True)
            ex = subprocess.run(["python3", f"{HERE}/read_exons.py"], input=p.stdout,
                                capture_output=True, text=True).stdout
            cnt = defaultdict(int)
            for line in ex.splitlines():
                x = line.split("\t")
                if len(x) >= 3:
                    cnt[(int(x[1]), int(x[2]))] += 1
            blocks, cur = [], None
            for (s, e) in sorted(cnt):
                s2, e2 = max(s, a), min(e, b)
                if e2 - s2 < 20:
                    continue
                if cur and s2 <= cur[1]:
                    cur = (cur[0], max(cur[1], e2), cur[2] + cnt[(s, e)])
                else:
                    if cur and cur[2] >= MIN_READS:
                        blocks.append(cur[:2])
                    cur = (s2, e2, cnt[(s, e)])
            if cur and cur[2] >= MIN_READS:
                blocks.append(cur[:2])
            if not blocks:
                continue
            seq = []
            for s, e in blocks:
                o = subprocess.run(["samtools", "faidx", HS, f"{c}:{s+1}-{e}"],
                                   capture_output=True, text=True).stdout.splitlines()
                seq.append("".join(o[1:]))
            fh.write(f">{reg}\n" + "".join(seq) + "\n")

    for sub in ("dna", "rna"):
        src = f"{W}/{fid}_{sub}.fa"
        if os.path.getsize(src) == 0:
            continue
        er_tier.ava(src, f"{W}/{fid}_{sub}.paf", threads=int(T))

    have_rna = [l[1:].strip() for l in open(f"{W}/{fid}_rna.fa") if l[0] == ">"]
    for sub in ("dna", "rna"):
        ns = regs if sub == "dna" else have_rna
        if len(ns) < 2:
            continue
        for denom in ("min", "max"):
            E = er(f"{W}/{fid}_{sub}", set(ns), denom)
            cs = comps(ns, E)

            def curated_in(comp):
                s = set()
                for reg in comp:
                    c, rng = reg.rsplit(":", 1)
                    a, b = rng.split("-")
                    a, b = int(a) - 1, int(b)
                    for cc, gs, ge, g in members[fid]:
                        if cc == c and gs < b and a < ge:
                            s.add(g)
                return s
            best = max(cs, key=lambda comp: len(curated_in(comp)))
            got = curated_in(best)
            rc = len(got) / len(members[fid])
            pr = len(got) / len(best)
            results[fid][(sub, denom)] = (rc, pr, len(best), len(ns), len(cs))
    print(f"{fid} done ({len(regs)} loci{' CAPPED' if capped else ''}, "
          f"{len(have_rna)} with reads)", flush=True)

print("\n" + "=" * 92)
print("2x2  substrate x denominator, scored at the FAMILY level vs Soto curated membership")
print("=" * 92)
hdr = ["dna/min", "dna/max", "rna/min", "rna/max"]
print(f"{'family':<9}{'n':>4}" + "".join(f"{h:>20}" for h in hdr))
print(f"{'':<13}" + "".join(f"{'recall / prec*':>20}" for _ in hdr))
agg = defaultdict(list)
for fid in FAMS:
    row = f"{fid:<9}{len(members[fid]):>4}"
    for sub in ("dna", "rna"):
        for denom in ("min", "max"):
            v = results[fid].get((sub, denom))
            if not v:
                row += f"{'--':>20}"
                continue
            rc, pr, sz, ntot, nc = v
            agg[(sub, denom)].append((rc, pr))
            row += f"{f'{rc:.2f} / {pr:.2f} (n={sz})':>20}"
    print(row)
print()
for k in (("dna", "min"), ("dna", "max"), ("rna", "min"), ("rna", "max")):
    v = agg[k]
    if not v:
        continue
    mr = sum(x[0] for x in v) / len(v)
    mp = sum(x[1] for x in v) / len(v)
    f1 = 2 * mr * mp / (mr + mp) if mr + mp else 0
    print(f"  {k[0]}/{k[1]:<4}  mean recall {mr:.3f}   mean precision* {mp:.3f}   F1 {f1:.3f}"
          f"   ({len(v)} families)")
print("\n  *precision is a LOWER BOUND -- PPIAL4G/H are real members missing from Soto's list.")
print("  PREDICTION ON RECORD: blocks lose under min(); if the denominator is the cause, "
      "rna/max should recover.")
print("DONE")
