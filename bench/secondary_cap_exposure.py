#!/usr/bin/env python3
"""Measure exposure of the read-conflict family definition to minimap2's secondary cap (N=5 -> <=6 placements
per read). For a co-located array, ALL of a read's within-array placements fall in the one region, so counting
alignment records per qname gives the read's multimapping degree directly. A spike at 6 = cap saturation
(reads that may have MORE unreported placements -> potentially missed conflict edges).

Reports per region: degree distribution (records/read), and the fraction of MULTIMAPPING reads (>=2 placements)
that are saturated at the 6-placement cap. Run with /home/juanfra/miniforge3/bin/python."""
import collections
import pysam

BAM = "/mnt/c/Users/jfris/Desktop/GGO.bam"
CAP = 6  # primary + N=5 secondaries

REGIONS = [
    # high-copy tandem arrays (where the cap could bite)
    ("TSPY array", "NC_073248.2", 34688121, 35060917),
    ("RBMY1 array", "NC_073248.2", 19558214, 19739408),
    ("chrY amplicon 25.7Mb", "NC_073248.2", 25690743, 25747730),
    # validated families (<=7 copies)
    ("DSFAM38 (7 copies)", "NC_073242.2", 6296327, 6401327),
    ("APOBEC3 family", "NC_086018.1", 36980802, 37112124),
    ("RABL2A locus", "NC_073235.2", 15131653, 15147533),
    ("RFPL2/3 span", "NC_086018.1", 30208643, 30376055),
]


def degree_hist(chrom, s, e):
    deg = collections.Counter()
    for r in pysam.AlignmentFile(BAM, "rb").fetch(chrom, s, e):
        if r.is_supplementary or r.is_unmapped:
            continue
        deg[r.query_name] += 1
    hist = collections.Counter(deg.values())
    return deg, hist


def main():
    print(f"{'region':24} {'reads':>6} {'multimap':>9} {'saturated(>=6)':>15} {'maxdeg':>7}")
    for name, c, s, e in REGIONS:
        deg, hist = degree_hist(c, s, e)
        nreads = len(deg)
        multimap = sum(1 for v in deg.values() if v >= 2)
        sat = sum(1 for v in deg.values() if v >= CAP)
        maxd = max(deg.values()) if deg else 0
        frac = f"{sat}/{multimap} ({100*sat/multimap:.0f}%)" if multimap else "0/0"
        print(f"{name:24} {nreads:6} {multimap:9} {frac:>15} {maxd:7}")
    print(f"\n[cap] minimap2 -ax splice:hq default N=5 -> max {CAP} placements/read. Saturation at {CAP} means a "
          f"read MAY have more unreported placements (potentially missed conflict edges in >6-copy arrays).")


if __name__ == "__main__":
    main()
