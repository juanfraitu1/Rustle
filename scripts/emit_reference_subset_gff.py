#!/usr/bin/env python3
"""
Emit a subset of reference GFF/GTF containing only listed transcript IDs.

Typical use after investigate_absent_ref_bam.py --emit-candidates:
  merge junction-backed reference isoforms that are absent from gffcompare tmap into
  evaluation or downstream tooling.

Usage:
  python3 scripts/emit_reference_subset_gff.py \\
    --reference GGO_clusterd.gff \\
    --ids absent_junction_candidates.tsv \\
    --out rare_ref_subset.gff

TSV: first column must be transcript_id (header row optional if it starts with transcript_id).
Alternatively --ids-one-per-line plain text.
"""

from __future__ import annotations

import argparse
import csv
import re
from pathlib import Path
from typing import Set


def load_ids_from_tsv(path: Path) -> Set[str]:
    ids: Set[str] = set()
    with path.open() as f:
        r = csv.reader(f, delimiter="\t")
        rows = list(r)
    if not rows:
        return ids
    start = 0
    if rows[0] and rows[0][0].lower() in ("transcript_id", "ref_id", "id"):
        start = 1
    for row in rows[start:]:
        if row and row[0].strip():
            ids.add(row[0].strip())
    return ids


def load_ids_from_txt(path: Path) -> Set[str]:
    return {ln.strip() for ln in path.open() if ln.strip() and not ln.startswith("#")}


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--reference", required=True, type=Path)
    g = ap.add_mutually_exclusive_group(required=True)
    g.add_argument("--ids", type=Path, help="TSV with transcript_id in column 1")
    g.add_argument("--ids-one-per-line", type=Path, help="Plain text IDs")
    ap.add_argument("--out", required=True, type=Path)
    args = ap.parse_args()

    if args.ids:
        want = load_ids_from_tsv(args.ids)
    else:
        want = load_ids_from_txt(args.ids_one_per_line)

    tid_re = re.compile(r'transcript_id\s+"([^"]+)"')
    args.out.parent.mkdir(parents=True, exist_ok=True)
    n = 0
    with args.reference.open() as inf, args.out.open("w") as out:
        for line in inf:
            if line.startswith("#"):
                out.write(line)
                continue
            if not line.strip():
                continue
            m = tid_re.search(line)
            if m and m.group(1) in want:
                out.write(line)
                n += 1
    print(f"Wrote {n} feature lines for {len(want)} target IDs → {args.out}")


if __name__ == "__main__":
    main()
