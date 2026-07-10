#!/usr/bin/env python3
"""Validate the candidate rule for rejecting readthrough transcripts admitted as COPIES.

A single-exon de-novo transcript spanning 30-250 kb is unspliced pre-mRNA / intronic pileup, not an mRNA.
When one is admitted as a family COPY it corrupts the copy set, and aligning reads against a 128 kb
"transcript" is quadratic (the RFPL hang). Real intronless genes (retrocopies, TSPYL1, GSPT2, JUND, histones)
must survive whatever rule we adopt.

Four candidate rules were tested. THREE FAILED; the fourth is the one this script measures.

  R1  "T overlaps any spliced transcript"                     -> too loose, never measured seriously
  R2  "some ASSEMBLED spliced transcript has an intron inside T's span"
      FAILS (sensitivity): misses the 128 kb RFPL giant, because that window is too sparsely expressed
      for a spliced transcript to be assembled. The rule depends on the assembler having already
      succeeded, which is exactly what fails there.
  R3  "any READ carries a junction fully inside T's span"
      FAILS (specificity): drops 21.2% of expressed annotated intronless genes (TSPYL1, GSPT2, DERPC,
      ATXN7L3B, JUND ...). A few stray spliced reads cross essentially any locus.
  R4  "T's span contains >= MIN_DISTINCT distinct junctions, each supported by >= 2 reads"   <-- ADOPTED
      A readthrough engulfs whole gene structures, so it contains MANY distinct junctions; an intronless
      gene catches at most a stray one or two. Junction evidence is read-level, so it does not require the
      assembler to have built a spliced transcript.

Also tested and rejected: span / longest-contained-read ratio ("no molecule attests a 128 kb exon"). No
separation -- real intronless genes reach ratio 7.1 (low expression => short reads) while one artifact sits
at 1.04 (a single 238 kb pre-mRNA read).

Usage:
    python bench/readthrough_rule_validation.py --bam GGO_mm.bam --gff GGO_genomic.gff \
        --artifact-gtf 'txdump/*.gtf' --controls single_exon_genes.tsv
"""
import argparse
import collections
import glob
import re
import statistics
import subprocess

MIN_SUPPORT = 2   # reads required before a junction counts as real (guards the alignment-noise floor)
MIN_DISTINCT = 5  # >= this many distinct junctions inside the span => readthrough. See the margin below.


def distinct_junctions(bam, chrom, start, end, min_support=MIN_SUPPORT):
    """Distinct splice junctions, each with >= min_support reads, lying ENTIRELY inside [start, end).

    Containment matters: a gene nested inside another gene's intron sees that host's junctions FLANKING it,
    never inside it, so a nested intronless gene is not penalised for its host's splicing.
    """
    out = subprocess.run(
        ["samtools", "view", "-F", "2308", bam, f"{chrom}:{start}-{end}"],
        capture_output=True, text=True,
    ).stdout.splitlines()
    j = collections.Counter()
    for line in out:
        f = line.split("\t")
        ref = int(f[3])
        for n, op in re.findall(r"(\d+)([MIDNSHP=X])", f[5]):
            n = int(n)
            if op in "M=XD":
                ref += n
            elif op == "N":
                if start <= ref and ref + n <= end:
                    j[(ref, ref + n)] += 1
                ref += n
    return sum(1 for v in j.values() if v >= min_support)


def load_gtf(path):
    tx = collections.defaultdict(lambda: {"exons": [], "reads": 0, "chrom": None})
    for line in open(path):
        f = line.rstrip("\n").split("\t")
        if len(f) < 9:
            continue
        m = re.search(r'transcript_id "([^"]*)"', f[8])
        if not m:
            continue
        t = tx[m.group(1)]
        t["chrom"] = f[0]
        r = re.search(r'reads "(\d+)"', f[8])
        if r:
            t["reads"] = int(r.group(1))
        if f[2] == "exon":
            t["exons"].append((int(f[3]), int(f[4])))
    for t in tx.values():
        t["exons"].sort()
    return tx


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--bam", required=True)
    ap.add_argument("--artifact-gtf", required=True, help="glob of copy_assign --gtf dumps")
    ap.add_argument("--controls", required=True, help="TSV chrom/start/end/strand/gene of intronless genes")
    ap.add_argument("--max-controls", type=int, default=120)
    args = ap.parse_args()

    # --- sensitivity: every single-exon de-novo transcript found in the dumps ---
    singles = {}
    for g in glob.glob(args.artifact_gtf):
        for _, v in load_gtf(g).items():
            if len(v["exons"]) == 1:
                s, e = v["exons"][0]
                singles[(v["chrom"], s, e)] = v
    art = sorted(
        (distinct_junctions(args.bam, c, s, e), (e - s) / 1000, f"{c}:{s}-{e}")
        for (c, s, e) in singles
    )
    print(f"SENSITIVITY -- single-exon de-novo transcripts: n={len(art)}")
    for j, kb, loc in art:
        print(f"  distinct_junctions={j:4d}  span={kb:7.1f}kb  {loc}")
    caught = sum(1 for j, _, _ in art if j >= MIN_DISTINCT)
    print(f"  flagged at >= {MIN_DISTINCT}: {caught}/{len(art)}   min observed = {art[0][0] if art else 'n/a'}\n")

    # --- specificity: expressed annotated intronless genes ---
    rows = [l.split("\t") for l in open(args.controls).read().splitlines()]
    vals = []
    for c, s, e, _strand, gene in rows:
        if len(vals) >= args.max_controls:
            break
        s, e = int(s), int(e)
        n = subprocess.run(["samtools", "view", "-c", "-F", "2308", args.bam, f"{c}:{s}-{e}"],
                           capture_output=True, text=True).stdout.strip()
        if not n or int(n) < 5:
            continue
        vals.append((distinct_junctions(args.bam, c, s, e), gene))
    vals.sort(reverse=True)
    fp = sum(1 for v, _ in vals if v >= MIN_DISTINCT)
    print(f"SPECIFICITY -- expressed annotated intronless genes: n={len(vals)}")
    print(f"  median distinct_junctions = {statistics.median([v for v, _ in vals])}")
    print(f"  worst: " + ", ".join(f"{g}={v}" for v, g in vals[:5]))
    print(f"  false positives at >= {MIN_DISTINCT}: {fp}/{len(vals)}\n")

    if art and vals:
        print(f"MARGIN: controls max {vals[0][0]}  <  artifacts min {art[0][0]}")


if __name__ == "__main__":
    main()
