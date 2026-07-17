#!/usr/bin/env python3
"""Decompose the 74-member coverage-INDEPENDENT floor (still undetected WITH ideal 40x coverage) into its
distinct causes, to test whether "K=0" is one information floor or a bucket of several gene-family aspects.

Analyzed on the TOP-UP BAM (soto_reads_topup.bam): every non-silent floor member has >=40 reads there, so
whatever stops it seeding is NOT depth. For each floor member we ask WHERE its reads sit and how uniquely:

  own primary reads   = primary alignments sitting mostly (>=50%) on the member span (its OWN reads)
  MAPQ0 fraction      = of those, how many multi-map (aligner cannot place them uniquely) -> ambiguity floor
  MAPQ>=60 fraction   = how many are UNIQUELY placed -> distinguishable reads
  secondary overlap   = alignments whose PRIMARY is elsewhere (a sibling) but that cover this span -> sink

Taxonomy:
  silent                 - 0 reads (no expressed template)
  secondary-sink         - reads exist but their PRIMARY is parked on a sibling (own primaries < seed floor,
                           many secondaries). Distinguishable in principle; a DETECTION-seeding gap, not a floor.
  ambiguous / K0-like    - own primaries present but MAPQ0-dominated: the aligner cannot place them uniquely
                           -> copies too similar in every read-visible region. The real information floor -> DNA.
  unique-but-unseeded    - own primaries present AND uniquely placed (MAPQ>=60) yet NO locus formed. Reads ARE
                           distinguishable; something OTHER than identity/depth stops detection (dosage, isoform-
                           only difference, or a genuine method gap). The "another aspect" bucket.
  dosage-starved         - own uniquely-placed primaries below the seed floor despite a well-covered family.

Run: /home/juanfra/miniforge3/bin/python bench/soto/soto_floor_decompose.py
"""
import csv
from collections import Counter, defaultdict

import pysam

D = "/mnt/linuxdisk/home/juanfraitu/winloci_data"
BAM = f"{D}/topup/soto_reads_topup.bam"
REC = "bench/soto/soto_topup_recovery.tsv"
SEED_FLOOR = 3          # GATE_MIN_READS
OWN_FRAC = 0.50         # a read is "own" if >=50% of its ref span sits on the member
MAPQ_UNIQUE = 60


def read_stats(bam, chrom, s, e):
    """Reads substantially ON the member span. A read is ambiguous (multi-maps to >=2 copies) if its primary
    here is MAPQ0 OR it appears here as a secondary (primary on a sibling). unique_own = primary here, MAPQ>=60."""
    unique_own = ambig_own_mapq0 = n_sec = 0
    for a in bam.fetch(chrom, max(0, s), e):
        if a.is_unmapped or a.is_supplementary:
            continue
        rs, re = a.reference_start, a.reference_end
        if re is None or re <= rs:
            continue
        ov = min(re, e) - max(rs, s)
        if ov <= 0 or ov / (re - rs) < OWN_FRAC:
            continue
        if a.is_secondary:
            n_sec += 1                       # multi-mapper: primary parked on a sibling
        elif a.mapping_quality == 0:
            ambig_own_mapq0 += 1             # primary here but tied -> multi-mapper
        elif a.mapping_quality >= MAPQ_UNIQUE:
            unique_own += 1                  # primary here, uniquely placed
    return unique_own, ambig_own_mapq0, n_sec


def classify(unique_own, ambig0, n_sec, gid, near_copy):
    ambiguous = ambig0 + n_sec               # all multi-mapper reads on the member
    if unique_own + ambiguous == 0:
        return "silent (no reads on member)"
    if unique_own >= SEED_FLOOR:
        # distinguishable OWN reads exist yet no own locus -> NOT an information floor
        if near_copy:
            return "distinguishable-but-MERGED (unique reads absorbed into a co-located locus)"
        return "distinguishable-but-unseeded (unique reads, no nearby locus -> detection gap)"
    # reads are dominated by multi-mappers -> exonic sequence not uniquely placeable
    if gid is not None and gid >= 0.999:
        return "expressed-K=0: young identical (exon AND genome ~identical)"
    return "expressed-K=0: exon-homogenized (reads tie, genome <99.9% -> gene-conversion candidate)"


def load_copies():
    idx = defaultdict(list)
    r = csv.reader(open(f"{D}/soto_topup.copies.tsv"), delimiter="\t"); next(r, None)
    for row in r:
        if len(row) >= 6:
            try:
                idx[row[3]].append((int(row[4]), int(row[5])))
            except ValueError:
                pass
    return idx


def near_detected(copies, chrom, s, e, win=100_000):
    m = (s + e) // 2
    return any(not (a - win > m or b + win < m) for (a, b) in copies.get(chrom, ()))


def main():
    bam = pysam.AlignmentFile(BAM, "rb")
    copies = load_copies()
    rows = []
    tax = Counter()
    for r in csv.DictReader(open(REC), delimiter="\t"):
        if r["topup_detected"] == "Y":
            continue  # floor = still missing with ideal coverage
        chrom, s, e = r["chrom"], int(r["start"]), int(r["end"])
        gid = None if r["best_sibling_id"] == "NA" else float(r["best_sibling_id"])
        unique_own, ambig0, n_sec = read_stats(bam, chrom, s, e)
        near = near_detected(copies, chrom, s, e, win=5_000_000)  # detector's co-located --win
        cls = classify(unique_own, ambig0, n_sec, gid, near)
        tax[cls] += 1
        rows.append([r["family_id"], r["gene"], chrom, s, e, r["primary_reads"], r["best_sibling_id"],
                     unique_own, ambig0, n_sec, "Y" if near else "N", cls])

    rows.sort(key=lambda x: x[11])
    with open("bench/soto/soto_floor_decomposition.tsv", "w", newline="") as f:
        w = csv.writer(f, delimiter="\t")
        w.writerow(["family_id", "gene", "chrom", "start", "end", "orig_primary_reads", "best_sibling_id",
                    "unique_own_reads", "ambig_own_mapq0", "secondaries", "detected_copy_nearby", "cause"])
        w.writerows(rows)

    print(f"\n=== 74-member floor decomposition (WITH ideal 40x coverage; corrected: secondaries = ambiguity) ===")
    total = sum(tax.values())
    for cls, n in sorted(tax.items(), key=lambda kv: -kv[1]):
        print(f"  {n:3d}  {cls}")
    print(f"  ---  {total} total")
    info = sum(n for c, n in tax.items() if c.startswith("expressed-K=0") or c.startswith("silent"))
    other = sum(n for c, n in tax.items() if c.startswith("distinguishable"))
    print(f"\n  expressed-K=0 information floor (+silent), needs DNA/PSV: {info}")
    print(f"  distinguishable -> NOT an info floor (mergeable / gap):   {other}")
    print(f"  wrote bench/soto/soto_floor_decomposition.tsv")


if __name__ == "__main__":
    main()
