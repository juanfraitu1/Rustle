#!/usr/bin/env python3
"""Merge per-family (or per-copy) chrom:start-end windows into non-overlapping per-contig windows
suitable for `copy_assign --regions` (see src/bin/copy_assign.rs:1241 validate_no_overlapping_regions
and :1447 the --families containment check, which requires every supplied family's full span to be
CONTAINED in exactly one --regions entry -- not split across two, not uncovered).

Merging every input interval that overlaps or touches another on the same contig guarantees this:
each output window is the union of a chain of overlapping/touching input intervals, so any single
input interval (e.g. one family's own span) that fed the merge is fully contained in the output
window that swallowed it.

Input: one or more files, each line whitespace/tab-separated with a `chrom:start-end` token
somewhere on it (e.g. `gw_units_v3.units.regions`'s `family_id<TAB>chrom:start-end` shape, or a
plain one-`chrom:start-end`-per-line file). The FIRST token matching `chrom:start-end` on each line
is used; blank lines and lines with no such token are skipped with a warning.

Output: one merged `chrom:start-end` per line, contigs in first-seen order, windows sorted by start.

usage: merge_regions.py <in1> [<in2> ...] --out <merged.txt>
"""
import argparse
import re
import sys
from collections import defaultdict

REGION_RE = re.compile(r"^([^\s:]+):(\d+)-(\d+)$")


def extract_region(line):
    for tok in line.split():
        m = REGION_RE.match(tok)
        if m:
            return m.group(1), int(m.group(2)), int(m.group(3))
    return None


def merge_intervals(intervals):
    """Merge overlapping OR touching [start, end) intervals (sorted by start)."""
    merged = []
    for lo, hi in sorted(intervals):
        if merged and lo <= merged[-1][1]:
            merged[-1] = (merged[-1][0], max(merged[-1][1], hi))
        else:
            merged.append((lo, hi))
    return merged


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("inputs", nargs="+")
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

    by_contig = defaultdict(list)
    contig_order = []
    n_lines = 0
    n_skipped = 0
    for path in args.inputs:
        with open(path) as fh:
            for line in fh:
                line = line.rstrip("\n")
                if not line.strip():
                    continue
                n_lines += 1
                r = extract_region(line)
                if r is None:
                    n_skipped += 1
                    continue
                chrom, lo, hi = r
                if chrom not in by_contig:
                    contig_order.append(chrom)
                by_contig[chrom].append((lo, hi))

    total_in = sum(len(v) for v in by_contig.values())
    with open(args.out, "w") as out:
        total_merged = 0
        total_span = 0
        for chrom in contig_order:
            merged = merge_intervals(by_contig[chrom])
            total_merged += len(merged)
            for lo, hi in merged:
                out.write(f"{chrom}:{lo}-{hi}\n")
                total_span += hi - lo

    print(f"read {n_lines} lines ({n_skipped} skipped, no chrom:start-end token found)")
    print(f"{total_in} input windows on {len(contig_order)} contigs -> {total_merged} merged windows")
    print(f"total merged span: {total_span} bp ({total_span / 1e6:.1f} Mb)")
    print(f"wrote {args.out}")


if __name__ == "__main__":
    main()
