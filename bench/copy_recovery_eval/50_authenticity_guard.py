#!/usr/bin/env python3
"""U5: tool-independent authenticity guard. For each recovered copy, count reads
whose PRIMARY alignment overlaps the copy's exons (SAM flags only). A recovery is
authentic iff that count >= K."""
import argparse, csv, sys
import pysam


def primary_support(bam_path, chrom, exons):
    """Count distinct reads with a PRIMARY alignment overlapping any exon."""
    n = 0
    seen = set()
    with pysam.AlignmentFile(bam_path, "rb") as bam:
        for (start, end) in exons:
            for read in bam.fetch(chrom, start, end):
                if read.is_secondary or read.is_supplementary or read.is_unmapped:
                    continue
                if read.query_name in seen:
                    continue
                seen.add(read.query_name)
                n += 1
    return n


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--bam", required=True)
    ap.add_argument("--recovery", required=True, help="rustle-VG recovery.tsv (ref_transcript,family_id,fsm,ism)")
    ap.add_argument("--exons", required=True, help="TSV: transcript_id,chrom,exon_starts(csv),exon_ends(csv)")
    ap.add_argument("--k", type=int, default=1)
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

    exon_map = {}
    with open(args.exons) as f:
        for row in csv.DictReader(f, delimiter="\t"):
            starts = [int(x) for x in row["exon_starts"].split(",") if x]
            ends = [int(x) for x in row["exon_ends"].split(",") if x]
            exon_map[row["transcript_id"]] = (row["chrom"], list(zip(starts, ends)))

    with open(args.recovery) as f, open(args.out, "w", newline="") as o:
        rows = list(csv.DictReader(f, delimiter="\t"))
        w = csv.writer(o, delimiter="\t")
        w.writerow(["ref_transcript", "family_id", "primary_support", "authentic"])
        out_rows = []
        for r in rows:
            tx = r["ref_transcript"]
            if tx not in exon_map:
                sys.stderr.write(f"[guard] no exons for {tx}, skipping\n"); continue
            chrom, exons = exon_map[tx]
            n = primary_support(args.bam, chrom, exons)
            out_rows.append((tx, r["family_id"], n, "true" if n >= args.k else "false"))
        out_rows.sort(key=lambda x: (x[1], x[0]))
        for row in out_rows:
            w.writerow(row)


if __name__ == "__main__":
    main()
