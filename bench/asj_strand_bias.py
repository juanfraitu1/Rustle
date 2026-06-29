#!/usr/bin/env python3
"""Strand-bias (StrandOddsRatio, SOR) filter for allele-specific-junction (ASJ) calls.

Rationale (longcallR, Nat Methods 2026): a junction that is "used" only by reads from a single
strand is a hallmark of an alignment/mapping artifact rather than a true biological splice choice.
longcallR rejects allele-specific junctions with strand-odds-ratio SOR >= 2. Our ASJ pipeline
(bench/allele_specific_junctions.py -> bench/asj_aggregate.py) never tallied read strand at all, so
strand-biased calls could slip through. This script closes that gap on the existing 475 calls.

For each ASJ call we re-fetch every read that SPANS the junction (ref_start <= donor AND
ref_end >= acceptor) and split it two ways:
  - USED  : the (donor, acceptor) intron is present in the read's CIGAR-derived intron chain
  - NOTUSED: it spans the locus but does not use that exact intron
and by strand (pysam read.is_reverse). That gives the 2x2:
      [ (used,    forward)=R+   (used,    reverse)=R- ]
      [ (notused, forward)=A+   (notused, reverse)=A- ]
"used" is treated as the REF allele and "notused" as the ALT allele of the GATK SOR contingency
table (the labelling is symmetric for the symmetricRatio term; the ref/alt ratios just measure
within-class strand balance, so the choice does not change the headline conclusion).

StrandOddsRatio (GATK), with a +1 pseudocount on each of R+, R-, A+, A-:
    R            = (R+ * A-) / (A+ * R-)
    symmetricRatio = R + 1/R
    refRatio     = min(R+, R-) / max(R+, R-)
    altRatio     = min(A+, A-) / max(A+, A-)
    SOR          = ln(symmetricRatio) + ln(refRatio) - ln(altRatio)
Reference: GATK StrandOddsRatioAnnotation (Broad Institute) and Brad Chapman's derivation; the
+1 pseudocount and this exact combination are the standard GATK implementation. Large SOR => the
"used"/"notused" split is strongly confounded with strand => likely artifact.

Cutoffs: GATK's common hard-filter for SNVs is SOR > 3.0 (we KEEP if SOR < 3.0, the default here).
longcallR uses SOR >= 2 as a REJECT for junctions (KEEP if SOR < 2.0); we also report that stricter
pass-count as a cross-check.

NOTE (Rust follow-up): the Rust production path src/bin/asj.rs cannot apply this gate yet because
read strand (is_reverse) is not threaded into AlignedRead in
src/rustle/vg_family/allele_specific_junctions.rs; that is a separate follow-up.

Run with /home/juanfra/miniforge3/bin/python (pysam + scipy). Deterministic.
"""
import argparse
import math
import os
from collections import defaultdict

import pysam
from scipy.stats import fisher_exact

HERE = os.path.dirname(os.path.abspath(__file__))
IN_TSV = os.path.join(HERE, "asj_calls.tsv")
OUT_TSV = os.path.join(HERE, "asj_calls_strandbias.tsv")
DEFAULT_BAM = "/home/juanfra/winloci_scratch/GGO_mm.bam"
SOR_KEEP = 3.0       # GATK SNV hard-filter cutoff: KEEP if SOR < 3.0
SOR_STRICT = 2.0     # longcallR junction REJECT threshold: KEEP if SOR < 2.0
TRANSITION = ({"A", "G"}, {"C", "T"})


def introns(aln):
    """CIGAR-derived intron chain (ref coords of each N block) for an aligned read."""
    out = []
    rpos = aln.reference_start
    for op, length in aln.cigartuples:
        if op == 3:                       # N: skipped region (intron)
            out.append((rpos, rpos + length)); rpos += length
        elif op in (0, 2, 7, 8):          # M/D/=/X consume ref
            rpos += length
    return out


def strand_table(bam, chrom, donor, acceptor):
    """2x2 strand-by-junction-usage counts over reads spanning the (donor, acceptor) junction.

    Returns (r_fwd, r_rev, a_fwd, a_rev) = (used+, used-, notused+, notused-) RAW (no pseudocount).
    """
    r_fwd = r_rev = a_fwd = a_rev = 0
    j = (donor, acceptor)
    # fetch a little past the acceptor so reads that just reach it are included
    for aln in bam.fetch(chrom, donor, acceptor + 1):
        if aln.is_unmapped or aln.is_secondary or aln.is_supplementary or aln.cigartuples is None:
            continue
        if aln.reference_start > donor or aln.reference_end < acceptor:
            continue                      # must SPAN the junction: ref_start <= donor AND ref_end >= acceptor
        used = j in introns(aln)
        rev = aln.is_reverse
        if used:
            if rev:
                r_rev += 1
            else:
                r_fwd += 1
        else:
            if rev:
                a_rev += 1
            else:
                a_fwd += 1
    return r_fwd, r_rev, a_fwd, a_rev


def sor(r_fwd, r_rev, a_fwd, a_rev):
    """GATK StrandOddsRatio with a +1 pseudocount on every cell (see module docstring)."""
    rp = r_fwd + 1.0
    rm = r_rev + 1.0
    ap = a_fwd + 1.0
    am = a_rev + 1.0
    ratio = (rp * am) / (ap * rm)
    sym = ratio + 1.0 / ratio
    ref_ratio = min(rp, rm) / max(rp, rm)
    alt_ratio = min(ap, am) / max(ap, am)
    return math.log(sym) + math.log(ref_ratio) - math.log(alt_ratio)


def fisher_strand_p(r_fwd, r_rev, a_fwd, a_rev):
    """Cross-check: Fisher exact strand-bias p over the raw 2x2 (used/notused x fwd/rev)."""
    try:
        _, p = fisher_exact([[r_fwd, r_rev], [a_fwd, a_rev]])
        return p
    except ValueError:
        return float("nan")


def main():
    ap = argparse.ArgumentParser(description="SOR strand-bias filter for ASJ calls")
    ap.add_argument("--bam", default=DEFAULT_BAM)
    ap.add_argument("--in", dest="in_tsv", default=IN_TSV)
    ap.add_argument("--out", default=OUT_TSV)
    args = ap.parse_args()

    bam = pysam.AlignmentFile(args.bam, "rb")
    refs = set(bam.references)

    with open(args.in_tsv) as fh:
        hdr = fh.readline().rstrip("\n").split("\t")
        rows = [dict(zip(hdr, line.rstrip("\n").split("\t"))) for line in fh if line.strip()]

    out_cols = hdr + ["R_fwd", "R_rev", "A_fwd", "A_rev", "sor", "fisher_strand_p", "sor_pass"]
    n_keep_3 = n_keep_2 = 0
    n_zero_used = 0   # calls with NO junction-using read in this BAM (support absent)
    fails = []   # (gene, sor, table, anchor_type) for SOR>=3.0
    with open(args.out, "w") as out:
        out.write("\t".join(out_cols) + "\n")
        for r in rows:
            chrom = r["chrom"]; donor = int(r["donor"]); acceptor = int(r["acceptor"])
            if donor > acceptor:
                donor, acceptor = acceptor, donor   # normalise (some calls list the anchor side first)
            if chrom not in refs:
                rf = rr = af = ar = 0
                s = float("nan"); fp = float("nan")
            else:
                rf, rr, af, ar = strand_table(bam, chrom, donor, acceptor)
                s = sor(rf, rr, af, ar)
                fp = fisher_strand_p(rf, rr, af, ar)
            if (rf + rr) == 0 and chrom in refs:
                n_zero_used += 1
            passed = (not math.isnan(s)) and s < SOR_KEEP
            if passed:
                n_keep_3 += 1
            if (not math.isnan(s)) and s < SOR_STRICT:
                n_keep_2 += 1
            if (not math.isnan(s)) and s >= SOR_KEEP:
                fails.append((r["gene"], s, (rf, rr, af, ar), r.get("anchor_type", "")))
            out.write("\t".join(str(r[c]) for c in hdr) +
                      f"\t{rf}\t{rr}\t{af}\t{ar}\t{s:.3f}\t{fp:.3e}\t{int(passed)}\n")

    n = len(rows)
    transversion = [r for r in rows if {r["alleleX"], r["alleleY"]} not in TRANSITION]
    # recompute SOR-pass on the transversion core and on PSMD2/DAXX (cheap, re-fetch)
    def call_sor(r):
        chrom = r["chrom"]; d = int(r["donor"]); a = int(r["acceptor"])
        if d > a:
            d, a = a, d
        if chrom not in refs:
            return float("nan")
        return sor(*strand_table(bam, chrom, d, a))

    tv_sor = [(r, call_sor(r)) for r in transversion]
    tv_keep3 = sum(1 for _, s in tv_sor if not math.isnan(s) and s < SOR_KEEP)
    tv_keep2 = sum(1 for _, s in tv_sor if not math.isnan(s) and s < SOR_STRICT)
    flagship = [r for r in rows if r["gene"] in ("PSMD2", "DAXX")]
    flag_sor = [(r, call_sor(r)) for r in flagship]

    print(f"ASJ strand-bias (SOR) filter over {n} calls  [BAM={os.path.basename(args.bam)}]")
    print(f"  KEEP (SOR < {SOR_KEEP}): {n_keep_3}/{n}  ({100*n_keep_3/n:.1f}%)   "
          f"FAIL: {n - n_keep_3}")
    print(f"  KEEP (SOR < {SOR_STRICT}, longcallR-strict): {n_keep_2}/{n}  ({100*n_keep_2/n:.1f}%)   "
          f"FAIL: {n - n_keep_2}")
    print(f"  (note: {n_zero_used} calls have ZERO junction-using reads in this BAM — the original "
          f"calls were made on single-mapping GGO.bam; in multi-mapping GGO_mm.bam that junction's "
          f"support is reassigned/absent, so their SOR reflects spanning-read strand skew, not "
          f"junction-usage strand skew)")
    print(f"  transversion core ({len(transversion)} calls): "
          f"{tv_keep3} pass SOR<{SOR_KEEP}, {tv_keep2} pass SOR<{SOR_STRICT}")
    print("  flagship PSMD2/DAXX:")
    for r, s in flag_sor:
        verdict = "PASS" if (not math.isnan(s) and s < SOR_KEEP) else "FAIL"
        strict = "PASS" if (not math.isnan(s) and s < SOR_STRICT) else "fail"
        print(f"    {r['gene']:6s} {r['donor']}-{r['acceptor']}  SOR={s:.3f}  "
              f"<{SOR_KEEP}:{verdict}  <{SOR_STRICT}:{strict}")
    if fails:
        print(f"  notable strand-bias FAILURES (SOR >= {SOR_KEEP}), top 10 by SOR:")
        for g, s, tbl, at in sorted(fails, key=lambda x: -x[1])[:10]:
            print(f"    {g:14s} SOR={s:.3f}  (R+,R-,A+,A-)={tbl}  [{at}]")
    else:
        print(f"  no call fails at SOR >= {SOR_KEEP}")
    print(f"\n[wrote {args.out}]")


if __name__ == "__main__":
    main()
