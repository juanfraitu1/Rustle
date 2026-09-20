#!/usr/bin/env python
"""Quantify the INTRON-RETENTION rescue potential for the currently-UNASSIGNABLE (Tied) reads.

The K=0 floor reads are tied because the copies are identical over the EXONS the read observes. Copies that are
exon-identical are often intron-DIVERGENT (introns evolve faster), so a read that RETAINS an intron carries an
intronic sequence that could distinguish the copies — an RNA-intrinsic lever the exonic-PSV gate does not use.

This measures, among the Tied reads, what fraction retain >=1 intron (an upper bound on what the lever could
rescue before the divergence factor). Consensus introns are derived from the reads themselves (the splice sites
seen in >= MIN_INTRON_READS reads); a Tied read RETAINS intron I if its reference span fully covers I but it did
NOT splice there (no matching N-gap) = it aligned through the intron.

Run: /home/juanfra/miniforge3/bin/python bench/intron_retention_rescue.py
"""
import os
from collections import Counter, defaultdict

import pysam

SCRATCH = "/home/juanfra/winloci_scratch"
BAM = f"{SCRATCH}/GGO_mm.bam"
REGIONS = f"{SCRATCH}/o2_regions.txt"
ASSIGN = f"{SCRATCH}/o2_chk.assignments.tsv"
MIN_INTRON_READS = 3        # a splice junction is "real" if >=3 reads use it
TOL = 6                     # bp tolerance matching a read's intron to a consensus intron


def read_introns(aln):
    """list of (donor, acceptor) genomic intron coords from N gaps in the CIGAR."""
    out = []
    if not aln.cigartuples:
        return out
    ref = aln.reference_start
    for op, ln in aln.cigartuples:
        if op in (0, 7, 8, 2):      # M/=/X/D consume reference
            ref += ln
        elif op == 3:               # N = intron
            out.append((ref, ref + ln))
            ref += ln
    return out


def load_tied():
    tied = {}
    with open(ASSIGN) as fh:
        fh.readline()
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) >= 5 and f[3] == "tied":
                tied[f[0]] = int(f[4])      # read_name -> n_decisive
    return tied


def main():
    tied = load_tied()
    print(f"{len(tied)} Tied (unassignable) reads in the assignment output")
    regions = []
    for line in open(REGIONS):
        r = line.split()[0]
        c, se = r.split(":")
        s, e = se.split("-")
        regions.append((c, int(s), int(e)))

    bam = pysam.AlignmentFile(BAM, "rb")
    found = 0
    retain = 0                 # tied reads retaining >=1 consensus intron
    retain_ndec0 = 0           # ... among the strict-K=0 (n_decisive==0) tied reads
    nd0_total = 0
    n_introns_seen = Counter()
    for (c, s, e) in regions:
        # 1) consensus introns in this region
        ic = Counter()
        recs = []
        try:
            it = list(bam.fetch(c, s, e))
        except ValueError:
            continue
        for a in it:
            if a.is_unmapped or a.is_secondary or a.is_supplementary:
                continue
            ins = read_introns(a)
            for d, ac in ins:
                ic[(d // TOL, ac // TOL)] += 1     # coarse-binned junction key
            recs.append((a.query_name, a.reference_start, a.reference_end, ins))
        consensus = {k for k, v in ic.items() if v >= MIN_INTRON_READS}
        if not consensus:
            continue
        # 2) classify the Tied reads in this region
        for name, rstart, rend, ins in recs:
            if name not in tied:
                continue
            found += 1
            nd0 = (tied[name] == 0)
            nd0_total += nd0
            own = {(d // TOL, ac // TOL) for d, ac in ins}
            # a consensus intron fully within the read's ref span that the read did NOT splice = retained
            ret = False
            for (dk, ak) in consensus:
                d, ac = dk * TOL, ak * TOL
                if rstart <= d and rend >= ac and (dk, ak) not in own:
                    # confirm the read didn't splice *near* it either
                    if not any(abs(od - dk) <= 1 and abs(oa - ak) <= 1 for (od, oa) in own):
                        ret = True
                        n_introns_seen[name] += 1
            retain += ret
            retain_ndec0 += (ret and nd0)

    print(f"Tied reads located in the regions: {found}/{len(tied)}")
    if found:
        print(f"\nTied reads that RETAIN >=1 intron: {retain}/{found} = {100*retain/found:.1f}%")
        print(f"  (an UPPER BOUND on intron-retention rescue, before the intron-divergence factor)")
    if nd0_total:
        print(f"strict-K=0 Tied reads (n_decisive==0, no exonic distinguishing column): {nd0_total}/{found}")
        print(f"  of those, retain >=1 intron: {retain_ndec0}/{nd0_total} = {100*retain_ndec0/nd0_total:.1f}%")
    print("\nNote: retention is only the FIRST factor. A retained intron rescues the read only if the copies'")
    print("intronic sequence DIFFERS there (factor 2). This number bounds how many Tied reads even carry intronic")
    print("sequence to exploit; multiply by the intron-divergence fraction for the realistic rescue.")


if __name__ == "__main__":
    main()
