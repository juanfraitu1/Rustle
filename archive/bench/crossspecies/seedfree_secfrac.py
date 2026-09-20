#!/usr/bin/env python3
"""Candidates from sec_frac -> E_r -> components. Sweeps the top-N cut; no seed, no annotation.

⚠⚠ TIER FIX, 2026-08-11. Was running its own `minimap2 -x asm20 -k11 -w5 -c --eqx -N 200 -p 0.02`
  with `nmatch/blocklen >= 0.60` and coverage `(qe-qs)/min(ql,tl)` -- a FOURTH tier on all three
  counts. Alignment, identity (1 - de) and the M1 coverage form now come from bench/soto/rustlib.py,
  the single mirror of denovo_pipeline.rs. PAFs are written as `cand_N.er.paf`: the panel-tier
  `cand_N.paf` files still on disk must never be served to the shipped rule under the same name.

⚠⚠ THE `best` COLUMNS ARE ORACLE-SELECTED. `max(C, key=genes_in)` chooses the component by the truth
  it is scored against -- the shape that killed the 0.237 blind-DNA purity -- and the column is
  MISLABELLED "largest" (it is the best-scoring one, which need not be the largest). Continuity only.
  The whole-partition TSV written per N is the quotable substrate.
"""
import os
import subprocess
import sys

sys.path.insert(0, os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "soto"))
import rustlib as er_tier  # noqa: E402  -- bench/soto/rustlib.py is THE definition of the E_r rule

OUT, HS = sys.argv[1], sys.argv[2]
BIN = 5000
MIN_PRIMARY = 5          # a bin needs some real coverage before its ratio means anything
TOPN = (200, 500, 1000, 2000)
NPIP = "/home/juanfra/winloci_scratch/seedfam/hs_npip.bed"

rows = []
for i, line in enumerate(open(f"{OUT}/bins.tsv")):
    if i == 0:
        continue
    b, pr, se, sf = line.split("\t")
    if int(pr) >= MIN_PRIMARY:
        rows.append((float(sf), int(b), int(pr), int(se)))
rows.sort(reverse=True)
print(f"  bins with >={MIN_PRIMARY} primary reads: {len(rows):,}")

truth = []
for line in open(NPIP):
    f = line.rstrip("\n").split("\t")
    if len(f) >= 4 and f[0] == "chr16":
        truth.append((int(f[1]), int(f[2]), f[3]))
print(f"  chr16 NPIP genes (scored against, NOT used to build): {len(truth)}\n")


def merge(bins):
    s = sorted(bins)
    out, cs, ce = [], None, None
    for b in s:
        if cs is None:
            cs = ce = b
        elif b <= ce + 1:
            ce = b
        else:
            out.append((cs * BIN, (ce + 1) * BIN))
            cs = ce = b
    if cs is not None:
        out.append((cs * BIN, (ce + 1) * BIN))
    return out


print(f"{'topN':>6}{'intervals':>11}{'>=10kb':>8}{'edges':>8}{'comps':>7}"
      f"{'ORACLE|C|':>10}{'NPIP genes in it':>18}{'density':>9}")
for N in TOPN:
    iv = [x for x in merge([b for _, b, _, _ in rows[:N]]) if x[1] - x[0] >= 10000]
    if len(iv) < 2:
        print(f"{N:>6}{len(iv):>11}   (too few)")
        continue
    reg = f"{OUT}/cand_{N}.regions"
    with open(reg, "w") as fh:
        for a, b in iv:
            fh.write(f"chr16:{a+1}-{b}\n")
    fa = f"{OUT}/cand_{N}.fa"
    if not os.path.exists(fa) or os.path.getsize(fa) == 0:
        with open(fa, "w") as fh:
            subprocess.run(["samtools", "faidx", HS, "-r", reg], stdout=fh, stderr=subprocess.DEVNULL)
    paf = f"{OUT}/cand_{N}.er.paf"      # `.er.` = shipped tier. NEVER reuse the old `.paf` name.
    if not os.path.exists(paf) or os.path.getsize(paf) == 0:
        er_tier.ava(fa, paf, threads=2)          # the SHIPPED all-vs-all tier, from one place
    names = [l[1:].strip() for l in open(fa) if l[0] == ">"]
    ns = set(names)
    E = er_tier.er_edges([(paf, er_tier.SENSITIVE_IDENTITY)], nodes=ns,
                         min_coverage=er_tier.MIN_COVERAGE)
    C = er_tier.components(names, E)
    with open(f"{OUT}/cand_{N}.partition.tsv", "w") as fh:
        fh.write("# tier\t" + er_tier.er_provenance(threads=2) + "\n")
        fh.write(f"# identity>={er_tier.SENSITIVE_IDENTITY}\tcoverage>={er_tier.MIN_COVERAGE}"
                 f"\tcovform=M1-axis\tdenom=min\ttopN={N}\n")
        fh.write("comp_idx\tcomp_size\tnode\n")
        for i, c in enumerate(C):
            for nd in c:
                fh.write(f"{i}\t{len(c)}\t{nd}\n")

    def genes_in(comp):
        g = set()
        for r in comp:
            a, b = r.split(":")[1].split("-")
            a, b = int(a) - 1, int(b)
            for gs, ge, nm2 in truth:
                if gs < b and a < ge:
                    g.add(nm2)
        return g
    # ⚠⚠ ORACLE-SELECTED — see the header. Continuity only; not quotable.
    best = max(C, key=lambda c: len(genes_in(c)))
    n = len(best)
    d = er_tier.density(best, E)
    print(f"{N:>6}{len(iv):>11}{sum(1 for a,b in iv if b-a>=10000):>8}{len(E):>8}{len(C):>7}"
          f"{n:>10}{f'{len(genes_in(best))}/{len(truth)}':>18}{d:>9.3f}")
