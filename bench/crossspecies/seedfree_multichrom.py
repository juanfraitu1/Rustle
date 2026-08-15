#!/usr/bin/env python3
"""Seed-free sec_frac -> E_r on several chromosomes: is the NPIP result a pattern or a one-off?

NPIP/chr16 gave: 157 candidates -> 78 components; largest NPIP-bearing component 36 loci, 18/19 genes,
density 0.530, PURITY 0.556 -- against blind DNA self-alignment's 0.136. This asks whether that
reproduces on chromosomes carrying different curated families.

NOTHING FROM THE CATALOG ENTERS THE PIPELINE. Candidates come from multimapper presence alone; the
curated families are read only to SCORE the components afterwards.

⚠ SERIAL by chromosome (WSL2: one heavy job at a time). Each chromosome is one streaming BAM pass plus
  one all-vs-all.
⚠ The top-N cut is the one chosen on chr16 (1000). Applying it unchanged elsewhere is the point -- a cut
  retuned per chromosome would prove nothing. Reported as fixed, not optimised.
⚠ PURITY IS A LOWER BOUND. On chr16, 16 "impure" loci included PDXDC1 (a real co-duplicated passenger),
  LOC105371131 (= CAT's NPIPB5), and 3 loci with no annotation at all -- catalog incompleteness, not
  method error.

⚠⚠ TIER FIX, 2026-08-11. This script used to run its own `minimap2 -x asm20 -k11 -w5 -c --eqx -N 200
  -p 0.02` and threshold on `nmatch/blocklen >= 0.60` with coverage `(qe-qs)/min(ql,tl)`. That is a
  FOURTH tier -- neither the aligner nor the identity metric nor the coverage form was the shipped one,
  so every number this script ever printed was at it. Alignment, identity (1 - de) and coverage (M1,
  axis-following) now ALL come from bench/soto/rustlib.py, the single mirror of denovo_pipeline.rs.
  The PAF is written under a NEW `.er.paf` name on purpose: the old `.paf` files still exist and a
  same-name cache hit would have silently served panel-tier alignments to the shipped rule.

⚠⚠ THE `best` COLUMNS BELOW ARE ORACLE-SELECTED AND MUST NOT BE QUOTED. `max(comps, key=hits)` picks
  the component by the answer being scored. That is the exact shape that killed the 0.237 blind-DNA
  purity. They are left in place ONLY as a continuity check against the 08-07 run; the quotable
  numbers come from the whole-partition TSV this script now writes next to each PAF, scored by the
  pre-registered scorer over ALL components and ALL truth families.
"""
import os
import subprocess
import sys
from collections import defaultdict

sys.path.insert(0, os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "soto"))
import rustlib as er_tier  # noqa: E402  -- bench/soto/rustlib.py is THE definition of the E_r rule

BAM = "/mnt/linuxdisk/home/juanfraitu/winloci_data/A119b.t2t.bam"
HS = "/mnt/c/Users/jfris/Desktop/Reference/chm13v2.0.fa"
TRUTH = "/mnt/c/Users/jfris/Desktop/Rustle/bench/soto/80_fams.gene_preferred.bed"
OUT = "/home/juanfra/winloci_scratch/seedfam/seedfree_multi"
BIN, MIN_PRIMARY, TOPN, MIN_LEN = 5000, 5, 1000, 10000
CHROMS = sys.argv[1].split(",") if len(sys.argv) > 1 else ["chr17", "chr15", "chr7", "chr1"]
os.makedirs(OUT, exist_ok=True)

fam = defaultdict(list)
for line in open(TRUTH):
    f = line.rstrip("\n").split("\t")
    if len(f) >= 4 and "|" in f[3]:
        g, fid = f[3].split("|", 1)
        fam[(fid, f[0])].append((int(f[1]), int(f[2]), g))

print(f"{'chrom':<7}{'fam':<9}{'members':>8}{'cands':>7}{'comps':>7}{'|C|':>6}"
      f"{'genes':>9}{'purity':>8}{'density':>9}")
for chrom in CHROMS:
    binf = f"{OUT}/{chrom}.bins.tsv"
    if not os.path.exists(binf):
        prim, sec = defaultdict(int), defaultdict(int)
        p = subprocess.Popen(["samtools", "view", "-F", "2052", BAM, chrom],
                             stdout=subprocess.PIPE, text=True, bufsize=1 << 20)
        for line in p.stdout:
            i = line.find('\t'); j = line.find('\t', i + 1)
            k = line.find('\t', j + 1); m = line.find('\t', k + 1)
            flag = int(line[i + 1:j]); pos = int(line[k + 1:m])
            (sec if flag & 256 else prim)[pos // BIN] += 1
        p.wait()
        with open(binf, "w") as fh:
            for b in sorted(set(prim) | set(sec)):
                fh.write(f"{b}\t{prim[b]}\t{sec[b]}\n")
    rows = []
    for line in open(binf):
        b, pr, se = line.split("\t")
        pr, se = int(pr), int(se)
        if pr >= MIN_PRIMARY and pr + se:
            rows.append((se / (pr + se), int(b)))
    rows.sort(reverse=True)
    top = sorted(b for _, b in rows[:TOPN])
    iv, cs, ce = [], None, None
    for b in top:
        if cs is None:
            cs = ce = b
        elif b <= ce + 1:
            ce = b
        else:
            iv.append((cs * BIN, (ce + 1) * BIN)); cs = ce = b
    if cs is not None:
        iv.append((cs * BIN, (ce + 1) * BIN))
    iv = [x for x in iv if x[1] - x[0] >= MIN_LEN]
    if len(iv) < 2:
        print(f"{chrom:<7}(too few candidates)")
        continue
    reg, fa = f"{OUT}/{chrom}.regions", f"{OUT}/{chrom}.fa"
    paf = f"{OUT}/{chrom}.er.paf"     # `.er.` = shipped tier. NEVER reuse the old `.paf` name.
    with open(reg, "w") as fh:
        for a, b in iv:
            fh.write(f"{chrom}:{a+1}-{b}\n")
    if not os.path.exists(fa) or os.path.getsize(fa) == 0:
        with open(fa, "w") as fh:
            subprocess.run(["samtools", "faidx", HS, "-r", reg], stdout=fh, stderr=subprocess.DEVNULL)
    if not os.path.exists(paf) or os.path.getsize(paf) == 0:
        er_tier.ava(fa, paf, threads=2)          # the SHIPPED all-vs-all tier, from one place
    names = [f"{chrom}:{a+1}-{b}" for a, b in iv]
    ns = set(names)
    E = er_tier.er_edges([(paf, er_tier.SENSITIVE_IDENTITY)], nodes=ns,
                         min_coverage=er_tier.MIN_COVERAGE)
    comps = er_tier.components(names, E)
    # Whole-partition dump: EVERY component including singletons, for the pre-registered scorer.
    # Nothing here is selected, ranked or filtered.
    with open(f"{OUT}/{chrom}.partition.tsv", "w") as fh:
        fh.write("# tier\t" + er_tier.er_provenance(threads=2) + "\n")
        fh.write(f"# identity>={er_tier.SENSITIVE_IDENTITY}\tcoverage>={er_tier.MIN_COVERAGE}"
                 f"\tcovform=M1-axis\tdenom=min\n")
        fh.write("comp_idx\tcomp_size\tnode\n")
        for i, c in enumerate(comps):
            for nd in c:
                fh.write(f"{i}\t{len(c)}\t{nd}\n")
    with open(f"{OUT}/{chrom}.edges.tsv", "w") as fh:
        for a, b in sorted(E):
            fh.write(f"{a}\t{b}\n")

    def span(r):
        a, b = r.split(":")[1].split("-")
        return int(a) - 1, int(b)
    for (fid, fc), mem in sorted(fam.items(), key=lambda kv: -len(kv[1])):
        if fc != chrom or len(mem) < 3:
            continue
        def hits(c):
            g = set()
            for r in c:
                a, b = span(r)
                for gs, ge, nm2 in mem:
                    if gs < b and a < ge:
                        g.add(nm2)
            return g
        # ⚠⚠ ORACLE-SELECTED — see the header. Continuity only; not quotable.
        best = max(comps, key=lambda c: len(hits(c)))
        n = len(best)
        d = er_tier.density(best, E)
        pure = sum(1 for r in best if any(gs < span(r)[1] and span(r)[0] < ge for gs, ge, _ in mem))
        print(f"{chrom:<7}{fid:<9}{len(mem):>8}{len(iv):>7}{len(comps):>7}{n:>6}"
              f"{f'{len(hits(best))}/{len(mem)}':>9}{pure/n:>8.3f}{d:>9.3f}")
