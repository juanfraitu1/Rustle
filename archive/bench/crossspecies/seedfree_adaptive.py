#!/usr/bin/env python3
"""Does N0's 75% miss rate survive an ADAPTIVE candidate budget?

THE CONFOUND THIS REMOVES. The all-chromosome run used a FIXED top-1000 bins per chromosome, and the
miss rate tracked candidates-per-curated-member almost monotonically: chr1 got 1.7 and missed 12 of 13
families; chr17 got 9.4 and worked on 2 of 3. A constant budget is generous where the benchmark is
sparse and starving where it is dense -- and it is dense exactly where the interesting families live.
That is the seventh absolute threshold to fail in this project.

THE FIX. Take every bin whose sec_frac clears a value, instead of a fixed COUNT of bins. The budget then
adapts to how much duplicated sequence each chromosome actually carries, and it does so using only the
signal -- no catalog, no annotation, nothing circular.

⚠ The cut is SWEPT, not chosen. sec_frac has a wide natural gap (duplicated median 0.94 vs single-copy
  0.05 on chr16), so any value in 0.5-0.95 sits in the gap rather than on a knife edge -- but reporting
  the curve is the only honest way to present it.
⚠ Bin counts are cached from the previous run; no BAM is re-streamed.
⚠ PURITY IS A LOWER BOUND -- the catalog is provably incomplete (PPIAL4G/H, 8 RefSeq NPIP genes,
  LOC105371131 = CAT's NPIPB5 all sit in the "impure" bucket).
⚠ Guard: a setting that proposes more than MAX_CAND intervals is reported and SKIPPED rather than run,
  so a runaway all-vs-all cannot silently eat the machine.

⚠⚠ TIER FIX, 2026-08-11. This script used to run its own `minimap2 -x asm20 -k11 -w5 -c --eqx -N 200
  -p 0.02` and threshold on `nmatch/blocklen >= 0.60` with coverage `(qe-qs)/min(ql,tl)` -- a FOURTH
  tier, matching neither the shipped aligner nor the shipped identity metric nor the M1 coverage form.
  Alignment, identity (1 - de) and coverage now all come from bench/soto/rustlib.py. Output PAFs are
  renamed `.er.paf` so the panel-tier `.paf` files already on disk can never be served to the new rule.

⚠⚠ THE PER-FAMILY `best` ROWS ARE ORACLE-SELECTED (`max(comps, key=hits)`) and must not be quoted:
  the component is chosen by the truth it is about to be scored against. Kept for continuity with the
  08-07 run only. The whole-partition TSV written beside each PAF is the quotable substrate.
"""
import os
import subprocess
import sys
from collections import defaultdict
import statistics as st

sys.path.insert(0, os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "soto"))
import rustlib as er_tier  # noqa: E402  -- bench/soto/rustlib.py is THE definition of the E_r rule

HS = "/mnt/c/Users/jfris/Desktop/Reference/chm13v2.0.fa"
TRUTH = "/mnt/c/Users/jfris/Desktop/Rustle/bench/soto/80_fams.gene_preferred.bed"
BINS = "/home/juanfra/winloci_scratch/seedfam/seedfree_multi"
OUT = "/home/juanfra/winloci_scratch/seedfam/seedfree_adaptive"
BIN, MIN_PRIMARY, MIN_LEN, MAX_CAND = 5000, 5, 10000, 700
CUTS = [0.95, 0.90, 0.50]
CHROMS = sys.argv[1].split(",") if len(sys.argv) > 1 else [
    "chr15", "chr1", "chr7", "chr2", "chr9", "chr16", "chr17", "chr5", "chr10", "chr22", "chr8", "chr11"]
os.makedirs(OUT, exist_ok=True)

fam = defaultdict(list)
for line in open(TRUTH):
    f = line.rstrip("\n").split("\t")
    if len(f) >= 4 and "|" in f[3]:
        g, fid = f[3].split("|", 1)
        fam[(fid, f[0])].append((int(f[1]), int(f[2]), g))


def span(r):
    a, b = r.split(":")[1].split("-")
    return int(a) - 1, int(b)


results = defaultdict(list)
for cut in CUTS:
    for chrom in CHROMS:
        binf = f"{BINS}/{chrom}.bins.tsv"
        if not os.path.exists(binf):
            continue
        keep = []
        for line in open(binf):
            b, pr, se = line.split("\t")
            pr, se = int(pr), int(se)
            if pr >= MIN_PRIMARY and pr + se and se / (pr + se) >= cut:
                keep.append(int(b))
        keep.sort()
        iv, cs, ce = [], None, None
        for b in keep:
            if cs is None:
                cs = ce = b
            elif b <= ce + 1:
                ce = b
            else:
                iv.append((cs * BIN, (ce + 1) * BIN)); cs = ce = b
        if cs is not None:
            iv.append((cs * BIN, (ce + 1) * BIN))
        iv = [x for x in iv if x[1] - x[0] >= MIN_LEN]
        tag = f"{chrom}_{int(cut*100)}"
        if len(iv) < 2:
            print(f"  [{tag}] {len(iv)} candidates — too few, skipped", flush=True)
            continue
        if len(iv) > MAX_CAND:
            print(f"  [{tag}] {len(iv)} candidates — ABOVE the {MAX_CAND} guard, skipped "
                  f"(would be a runaway all-vs-all)", flush=True)
            continue
        reg, fa = f"{OUT}/{tag}.regions", f"{OUT}/{tag}.fa"
        paf = f"{OUT}/{tag}.er.paf"     # `.er.` = shipped tier. NEVER reuse the old `.paf` name.
        with open(reg, "w") as fh:
            for a, b in iv:
                fh.write(f"{chrom}:{a+1}-{b}\n")
        if not os.path.exists(fa) or os.path.getsize(fa) == 0:
            with open(fa, "w") as fh:
                subprocess.run(["samtools", "faidx", HS, "-r", reg], stdout=fh,
                               stderr=subprocess.DEVNULL)
        if not os.path.exists(paf) or os.path.getsize(paf) == 0:
            er_tier.ava(fa, paf, threads=2)      # the SHIPPED all-vs-all tier, from one place
        names = [f"{chrom}:{a+1}-{b}" for a, b in iv]
        ns = set(names)
        E = er_tier.er_edges([(paf, er_tier.SENSITIVE_IDENTITY)], nodes=ns,
                             min_coverage=er_tier.MIN_COVERAGE)
        comps = er_tier.components(names, E)
        with open(f"{OUT}/{tag}.partition.tsv", "w") as fh:
            fh.write("# tier\t" + er_tier.er_provenance(threads=2) + "\n")
            fh.write(f"# identity>={er_tier.SENSITIVE_IDENTITY}\tcoverage>={er_tier.MIN_COVERAGE}"
                     f"\tcovform=M1-axis\tdenom=min\tsec_frac_cut={cut}\n")
            fh.write("comp_idx\tcomp_size\tnode\n")
            for i, c in enumerate(comps):
                for nd in c:
                    fh.write(f"{i}\t{len(c)}\t{nd}\n")
        for (fid, fc), mem in fam.items():
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
            pure = sum(1 for r in best
                       if any(gs < span(r)[1] and span(r)[0] < ge for gs, ge, _ in mem))
            results[cut].append(dict(chrom=chrom, fid=fid, n=len(mem), cands=len(iv),
                                     rec=len(hits(best)) / len(mem), pur=pure / n, csz=n))
        print(f"  [{tag}] {len(iv)} candidates -> {len(comps)} components", flush=True)

print("\n" + "=" * 74)
print("ADAPTIVE BUDGET vs the fixed top-1000 baseline (WORKS = recall>=0.7 AND purity>=0.5)")
print("=" * 74)
print(f"{'sec_frac cut':<14}{'pairs':>7}{'WORKS':>8}{'MERGED':>9}{'MISSED':>9}"
      f"{'med recall':>12}{'med cands':>11}")
print(f"{'fixed top-1000':<14}{55:>7}{'6 (11%)':>8}{'8 (15%)':>9}{'41 (75%)':>9}"
      f"{0.00:>12.2f}{132:>11}")
for cut in CUTS:
    v = results[cut]
    if not v:
        continue
    w = sum(1 for r in v if r['rec'] >= 0.7 and r['pur'] >= 0.5)
    m = sum(1 for r in v if r['rec'] >= 0.7 and r['pur'] < 0.5)
    x = len(v) - w - m
    print(f"{cut:<14.2f}{len(v):>7}{f'{w} ({w/len(v):.0%})':>8}{f'{m} ({m/len(v):.0%})':>9}"
          f"{f'{x} ({x/len(v):.0%})':>9}{st.median([r['rec'] for r in v]):>12.2f}"
          f"{st.median([r['cands'] for r in v]):>11.0f}")
