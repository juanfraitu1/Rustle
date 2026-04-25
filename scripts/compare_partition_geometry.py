#!/usr/bin/env python3
"""Compare partition geometry TSVs from StringTie (PARITY_PARTITION_TSV) and Rustle (RUSTLE_PARITY_PARTITION_TSV).

Both use the same schema from `parity_partition_dump` / `parity_partition_emit_stringtie`:
  kind, source, chrom, start, end, strand, signature

Join key: (chrom, int(start), int(end)). Report match rate and first mismatches.

Usage:
  python3 compare_partition_geometry.py stringtie_partition.tsv rustle_partition.tsv
  python3 compare_partition_geometry.py a.tsv b.tsv --strict   # exit 1 if any mismatch
  python3 compare_partition_geometry.py a.tsv b.tsv --mismatch-tsv mismatches.tsv --mismatch-max 25
  python3 compare_partition_geometry.py a.tsv b.tsv --segment-slop 3   # multiset segment coords ±slop
"""

from __future__ import annotations

import argparse
import csv
import sys
from pathlib import Path


def _parse_flat_segments(signature: str) -> list[tuple[int, int]]:
    """All (start, end) segments across chains, in canonical chain order (TSV `|` then `+`)."""
    if not signature or not signature.strip():
        return []
    segs: list[tuple[int, int]] = []
    for chain in signature.split("|"):
        chain = chain.strip()
        if not chain:
            continue
        for part in chain.split("+"):
            part = part.strip()
            if not part or "-" not in part:
                continue
            a, b = part.split("-", 1)
            segs.append((int(a), int(b)))
    return segs


def segment_multiset_fuzzy_match(sig_a: str, sig_b: str, slop: int) -> bool:
    """Same number of segments and pairwise match after sorting by (start, end), within ±slop."""
    if slop <= 0:
        return sig_a == sig_b
    sa = sorted(_parse_flat_segments(sig_a))
    sb = sorted(_parse_flat_segments(sig_b))
    if len(sa) != len(sb):
        return False
    for (x0, x1), (y0, y1) in zip(sa, sb):
        if abs(x0 - y0) > slop or abs(x1 - y1) > slop:
            return False
    return True


def load_rows(path: Path) -> dict[tuple[str, int, int], tuple[str, str]]:
    """key -> (source, signature)"""
    out: dict[tuple[str, int, int], tuple[str, str]] = {}
    with path.open(newline="") as f:
        r = csv.DictReader(f, delimiter="\t")
        for row in r:
            if row.get("kind") != "partition_geometry_v1":
                continue
            k = (row["chrom"], int(row["start"]), int(row["end"]))
            out[k] = (row["source"], row["signature"])
    return out


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("stringtie_tsv", type=Path)
    ap.add_argument("rustle_tsv", type=Path)
    ap.add_argument("--strict", action="store_true")
    ap.add_argument("--show", type=int, default=12, help="max mismatch examples")
    ap.add_argument(
        "--mismatch-tsv",
        type=Path,
        default=None,
        help="write joined key mismatches (chrom, start, end, st_signature, ru_signature)",
    )
    ap.add_argument(
        "--mismatch-max",
        type=int,
        default=0,
        help="cap rows written to --mismatch-tsv (0 = all mismatches)",
    )
    ap.add_argument(
        "--segment-slop",
        type=int,
        default=0,
        help="if >0, also report multiset segment match allowing ±slop bp on each segment end "
        "(start and displayed right end as in the TSV); chain topology is not relaxed",
    )
    args = ap.parse_args()

    st = load_rows(args.stringtie_tsv)
    ru = load_rows(args.rustle_tsv)
    keys = sorted(set(st) & set(ru))
    matched = sum(1 for k in keys if st[k][1] == ru[k][1])
    slop = max(0, int(args.segment_slop))
    fuzzy_matched = (
        sum(
            1
            for k in keys
            if segment_multiset_fuzzy_match(st[k][1], ru[k][1], slop)
        )
        if slop > 0
        else matched
    )
    missing_st = sorted(set(ru) - set(st))
    missing_ru = sorted(set(st) - set(ru))
    mism = [k for k in keys if st[k][1] != ru[k][1]]

    print(f"Joined bundles: {len(keys)}")
    print(f"Signature match: {matched} / {len(keys)} ({100.0 * matched / len(keys) if keys else 0:.2f}%)")
    if slop > 0:
        print(
            f"Multiset segment match (±{slop} bp): {fuzzy_matched} / {len(keys)} "
            f"({100.0 * fuzzy_matched / len(keys) if keys else 0:.2f}%)"
        )
    print(f"Only in StringTie file: {len(missing_ru)}")
    print(f"Only in Rustle file: {len(missing_st)}")
    if mism and args.show:
        print(f"First mismatches (up to {args.show}):")
        for k in mism[: args.show]:
            print(f"  {k}")
            print(f"    ST: {st[k][1][:120]}{'...' if len(st[k][1]) > 120 else ''}")
            print(f"    RU: {ru[k][1][:120]}{'...' if len(ru[k][1]) > 120 else ''}")

    if args.mismatch_tsv is not None:
        cap = args.mismatch_max if args.mismatch_max > 0 else len(mism)
        to_write = mism[:cap]
        with args.mismatch_tsv.open("w", newline="") as f:
            w = csv.writer(f, delimiter="\t")
            w.writerow(["chrom", "start", "end", "st_signature", "ru_signature"])
            for chrom, start, end in to_write:
                k = (chrom, start, end)
                w.writerow([chrom, start, end, st[k][1], ru[k][1]])
        print(f"Wrote {len(to_write)} mismatch row(s) to {args.mismatch_tsv}")

    if args.strict and (mism or missing_st or missing_ru):
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
