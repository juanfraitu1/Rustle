#!/usr/bin/env python3
"""Find the Rustle bundle containing a given genomic locus.

Reads RUSTLE_PARITY_PARTITION_TSV (or any TSV with chrom/start/end/strand
columns) and reports bundles whose interval overlaps the requested locus.

Usage:
    identify_locus_bundle.py --partition-tsv path.tsv --locus chrom:start-end
"""
from __future__ import annotations
import argparse, csv, sys
from pathlib import Path

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--partition-tsv", required=True, type=Path)
    ap.add_argument("--locus", required=True)
    args = ap.parse_args()

    chrom, rng = args.locus.split(":")
    lstart, lend = map(int, rng.split("-"))

    with args.partition_tsv.open() as f:
        for row in csv.DictReader(f, delimiter="\t"):
            if row.get("chrom") != chrom:
                continue
            bs, be = int(row["start"]), int(row["end"])
            if be >= lstart and bs <= lend:
                print(f"{row['chrom']}\t{bs}\t{be}\t{row.get('strand','.')}\t{row.get('signature','-')[:80]}")

if __name__ == "__main__":
    sys.exit(main())
