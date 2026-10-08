#!/usr/bin/env python3
"""Side-by-side read view of one run (docs/PREREG_o3_maternal_reference_2026-10-08.md, Amendment 1): the reads of each absent locus on the
reference haplotype (where the copy is missing) and on the other haplotype (where it is present).

    O3_REF=mat side.py     # WR/fate/side.json (needs WR/fate/fate.json, WR/truth/loci.tsv and both W/map/reads.<hap>.all.bam)

Per locus with >= 3 reads: `pairs` = [read, fate on the reference, de ref, MAPQ ref, de other, MAPQ other]; for LARGE loci (>= 20 reads) also
`ref_track` (the nearest reference paralog, where the reads land) and `other_track` (the locus itself): per-bin coverage and mismatch counts of
the primary alignments, the interval cut in NBINS equal bins (positions are fractions of the interval, not shared coordinates)."""
import csv
import json
import os
import statistics
import sys

import pysam

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import common as C  # noqa: E402

NBINS = 60


def add_alignment(cov, mis, ref_start, cigar, lo, hi):
    """Add one alignment to the per-bin coverage / mismatch counts of the interval [lo, hi) cut in len(cov) equal bins.
    cigar: pysam (op, length) tuples with --eqx operators: 7 '=', 8 'X', 0 'M' (covered, mismatch unknown), 1 'I', 2 'D', 3 'N' (intron: nothing).
    A deleted base counts as covered and mismatched; an insertion adds one mismatch at the current reference base and no coverage."""
    n = len(cov)

    def b(p):
        return min(n - 1, (p - lo) * n // (hi - lo))
    pos = ref_start
    for op, ln in cigar:
        if op in (0, 7, 8, 2):
            for p in range(max(pos, lo), min(pos + ln, hi)):
                k = b(p)
                cov[k] += 1
                if op in (8, 2):
                    mis[k] += 1
            pos += ln
        elif op == 3:
            pos += ln
        elif op == 1 and lo <= pos < hi:
            mis[b(pos)] += 1


def pair_rows(reads, recs_ref, recs_other, fate_of):
    """[[read, fate, de ref, MAPQ ref, de other, MAPQ other]] for the reads with a primary record on both haplotypes"""
    out = []
    for n in reads:
        r, o = recs_ref.get(n), recs_other.get(n)
        if r and o and r[0].primary and o[0].primary:
            out.append([n, fate_of[n], r[0].de, r[0].mapq, o[0].de, o[0].mapq])
    return out


def track(bam_path, names, target, al, nbins=NBINS):
    """per-bin coverage and mismatches of the primary alignments of `names` overlapping target = (accession, lo, hi)"""
    cov, mis = [0] * nbins, [0] * nbins
    with pysam.AlignmentFile(bam_path) as bam:
        idx = {C.accession(s["SN"], al): s["SN"] for s in bam.header.to_dict()["SQ"]}
        for rd in bam.fetch(idx[target[0]], target[1], target[2]):
            if rd.is_unmapped or rd.is_secondary or rd.is_supplementary or rd.query_name not in names:
                continue
            add_alignment(cov, mis, rd.reference_start, rd.cigartuples, target[1], target[2])
    return dict(target=list(target), cov=cov, mis=mis)


def median(xs):
    return statistics.median(xs) if xs else None


def main():
    al = C.alias()
    fate = json.load(open(f"{C.WR}/fate/fate.json"))["loci"]
    loci_ = {r["locus"]: r for r in csv.DictReader(open(f"{C.WR}/truth/loci.tsv"), delimiter="\t")}
    bam_ref, bam_oth = f"{C.W}/map/reads.{C.REF}.all.bam", f"{C.W}/map/reads.{C.OTHER}.all.bam"
    recs_ref, recs_oth = C.read_records(bam_ref, al), C.read_records(bam_oth, al)
    out = {}
    for k, v in fate.items():
        if v["n"] < 3:
            continue
        fate_of = {r[0]: r[1] for r in v["reads"]}
        pairs = pair_rows(sorted(fate_of), recs_ref, recs_oth, fate_of)
        row = dict(n=v["n"], kind=v["kind"], pairs=pairs, de_ref_median=median([p[2] for p in pairs if p[2] is not None]),
                   de_other_median=median([p[4] for p in pairs if p[4] is not None]),
                   mapq_ref_median=median([p[3] for p in pairs]), mapq_other_median=median([p[5] for p in pairs]))
        if v["n"] >= 20 and v["paralog"]:
            L = loci_[k]
            row["ref_track"] = track(bam_ref, set(fate_of), tuple(v["paralog"]), al)
            row["other_track"] = track(bam_oth, set(fate_of), (L["chrom"], int(L["start"]), int(L["end"])), al)
        out[k] = row
    json.dump(out, open(f"{C.WR}/fate/side.json", "w"))
    print(f"reference {C.REF}: {len(out)} loci with >= 3 reads; tracks for {sum(1 for r in out.values() if 'ref_track' in r)}")
    for k, r in out.items():
        print(f"{k}\tn={r['n']}\tde {C.REF} {r['de_ref_median']}\tde {C.OTHER} {r['de_other_median']}\tMAPQ {r['mapq_ref_median']} / {r['mapq_other_median']}")


if __name__ == "__main__":
    main()
