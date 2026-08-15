#!/usr/bin/env python3
"""Do extra PAF columns let us LOWER the identity/coverage floors without polluting the graph?

MOTIVATION, measured on PPIAL4: at the shipped floors (id>=0.80, cov>=0.50) **100% of the false pairs
pass** -- 451 of 451 pairs involving a PPIA processed retrocopy. Median false identity is 0.893 and
median false coverage 0.995, so neither shipped signal does any work against that failure mode.

CANDIDATE COLUMNS, all already emitted by minimap2, none requiring new computation:
    cm/kb        chained MINIMIZERS per kb of alignment  -- shared seeds, not shared bases
    nblk         number of '=' runs in the cg CIGAR      -- is the alignment solid or shattered
    covlong      aligned query span / max(qlen,tlen)     -- coverage of the LONGER sequence
    rlfrac       rl / blocklen                           -- repetitive-seed fraction
    s2/s1        second-best chain score / best          -- how distinctive the chain is

⚠ SCALE-FREE ONLY. On PPIAL4 the best-looking cut was blocklen >= 700, which is 100% true / 4% false --
  but PPIAL4's members are 760 bp, so that is coverage-of-the-seed wearing an absolute threshold. An
  absolute threshold has failed FIVE separate times in this project. Every candidate here is a ratio
  or a density.
⚠ The signals were chosen after seeing them separate on ONE family. This run is the out-of-sample test.
⚠ TRUTH IS CONSERVATIVE AND BIASED AGAINST US: positive = the locus overlaps a curated member of this
  family; everything else is negative. We already know some "negatives" are real members the catalog
  omits (PPIAL4G/H, and 8 RefSeq NPIP genes), so measured false-positive rates are UPPER bounds.
"""
import os
import re
import subprocess
import sys
import statistics as st
from collections import defaultdict
sys.path.insert(0, os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "soto"))
import rustlib as er_tier  # noqa: E402  -- bench/soto/rustlib.py is THE definition of the E_r rule


W = "/home/juanfra/winloci_scratch/seedfam/rep2x2"
OUT = "/home/juanfra/winloci_scratch/seedfam/pafcols"
TRUTH = "/mnt/c/Users/jfris/Desktop/Rustle/bench/soto/80_fams.gene_preferred.bed"
FAMS = sys.argv[1].split(",") if len(sys.argv) > 1 else [
    "ID_431", "ID_131", "ID_391", "ID_395", "ID_400", "ID_116", "ID_8", "ID_468"]
os.makedirs(OUT, exist_ok=True)

members = defaultdict(list)
for line in open(TRUTH):
    f = line.rstrip("\n").split("\t")
    if len(f) >= 4 and "|" in f[3]:
        g, fid = f[3].split("|", 1)
        members[fid].append((f[0], int(f[1]), int(f[2]), g))


def is_member(fid, reg):
    c, rng = reg.rsplit(":", 1)
    a, b = rng.split("-")
    a, b = int(a) - 1, int(b)
    return any(cc == c and gs < b and a < ge for cc, gs, ge, _ in members[fid])


rows = []
for fid in FAMS:
    fa = f"{W}/{fid}_dna.fa"
    if not os.path.exists(fa):
        print(f"  {fid}: no loci fasta, skipped", flush=True)
        continue
    paf = f"{OUT}/{fid}.paf"
    if not os.path.exists(paf) or os.path.getsize(paf) == 0:
        # ⚠ TIER CORRECTION 2026-08-10 (defect B1). This all-vs-all pass ran `-x asm20 -c --eqx
        # -N 200 -p 0.02` with NO `-X`. The comment justified it as "LOW floors on purpose: -p 0.02
        # keeps secondary chains" -- but `-N`/`-p` are INERT at the shipped tier, and this script
        # exists to characterise the columns of the PAF the BINARY produces. Now the shipped command.
        er_tier.ava(fa, paf, threads=2)
    n = 0
    for line in open(paf):
        p = line.rstrip("\n").split("\t")
        if len(p) < 12:
            continue
        q, t = p[0], p[5]
        if q == t:
            continue
        ql, qs, qe, tl, nm, bl = int(p[1]), int(p[2]), int(p[3]), int(p[6]), int(p[9]), int(p[10])
        if bl == 0:
            continue
        tag = {}
        for f in p[12:]:
            k, _, v = f.split(":", 2)
            tag[k] = v
        blocks = [int(x) for x in re.findall(r"(\d+)=", tag.get("cg", ""))] or [0]
        rows.append(dict(
            fid=fid, true=(is_member(fid, q) and is_member(fid, t)),
            idy=nm / bl,
            cov=(qe - qs) / max(min(ql, tl), 1),
            covlong=(qe - qs) / max(max(ql, tl), 1),
            cmkb=int(tag.get("cm", 0)) / max(bl, 1) * 1000,
            nblk=len(blocks),
            blkfrac=max(blocks) / max(bl, 1),
            rlfrac=int(tag.get("rl", 0)) / max(bl, 1),
            s2s1=int(tag.get("s2", 0)) / max(int(tag.get("s1", 1)), 1)))
        n += 1
    print(f"  {fid}: {n} pair records", flush=True)

T = [r for r in rows if r["true"]]
F = [r for r in rows if not r["true"]]
print(f"\nrecords: TRUE (both curated members) {len(T)}   FALSE (at least one not) {len(F)}")


def auc(key, rev=False):
    A = sorted(r[key] for r in T)
    B = sorted(r[key] for r in F)
    if not A or not B:
        return float("nan")
    import bisect
    w = 0
    for x in A:
        lo = bisect.bisect_left(B, x)
        hi = bisect.bisect_right(B, x)
        w += (lo + (hi - lo) * 0.5) if not rev else (len(B) - hi + (hi - lo) * 0.5)
    return w / (len(A) * len(B))


print(f"\n{'signal':<10}{'TRUE med':>10}{'FALSE med':>11}{'AUC':>8}")
for key, rev in (("idy", False), ("cov", False), ("covlong", False), ("cmkb", False),
                 ("nblk", True), ("blkfrac", False), ("rlfrac", True), ("s2s1", True)):
    print(f"{key:<10}{st.median([r[key] for r in T]):>10.4f}"
          f"{st.median([r[key] for r in F]):>11.4f}{auc(key, rev):>8.3f}")

print("\nEdge rule variants — TRUE kept / FALSE admitted (FALSE is an UPPER bound):")


def score(pred, lab):
    tk = sum(1 for r in T if pred(r))
    fk = sum(1 for r in F if pred(r))
    print(f"  {lab:<52} {tk}/{len(T)} = {tk/max(len(T),1):>5.0%}   {fk}/{len(F)} = {fk/max(len(F),1):>5.0%}")


score(lambda r: r["idy"] >= 0.80 and r["cov"] >= 0.50, "SHIPPED  id>=0.80  cov>=0.50")
score(lambda r: r["idy"] >= 0.60 and r["cov"] >= 0.30, "relaxed  id>=0.60  cov>=0.30")
for thr in (0.70, 0.80, 0.90):
    score(lambda r, t=thr: r["idy"] >= 0.60 and r["cov"] >= 0.30 and r["covlong"] >= t,
          f"relaxed + covlong>={thr}")
for thr in (40, 60, 80):
    score(lambda r, t=thr: r["idy"] >= 0.60 and r["cov"] >= 0.30 and r["cmkb"] >= t,
          f"relaxed + cm/kb>={thr}")
score(lambda r: r["idy"] >= 0.60 and r["cov"] >= 0.30 and r["covlong"] >= 0.80 and r["cmkb"] >= 40,
      "relaxed + covlong>=0.80 + cm/kb>=40")
print("\nPer-family, best variant vs shipped:")
print(f"{'family':<9}{'pairs T/F':>12}{'shipped T/F':>16}{'relaxed+covlong0.8 T/F':>26}")
for fid in FAMS:
    t = [r for r in T if r["fid"] == fid]
    f = [r for r in F if r["fid"] == fid]
    if not t and not f:
        continue
    s = (sum(1 for r in t if r["idy"] >= 0.80 and r["cov"] >= 0.50),
         sum(1 for r in f if r["idy"] >= 0.80 and r["cov"] >= 0.50))
    v = (sum(1 for r in t if r["idy"] >= 0.60 and r["cov"] >= 0.30 and r["covlong"] >= 0.80),
         sum(1 for r in f if r["idy"] >= 0.60 and r["cov"] >= 0.30 and r["covlong"] >= 0.80))
    print(f"{fid:<9}{f'{len(t)}/{len(f)}':>12}{f'{s[0]}/{s[1]}':>16}{f'{v[0]}/{v[1]}':>26}")
