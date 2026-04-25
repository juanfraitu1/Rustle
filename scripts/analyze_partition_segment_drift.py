#!/usr/bin/env python3
"""Analyze per-segment boundary drift between StringTie and Rustle partition signatures.

Focuses on bundles where junction sets already match (or are accepted as matched), then measures
segment-level coordinate deltas to guide Batch C boundary-normalization fixes.
"""

from __future__ import annotations

import argparse
import csv
from collections import Counter, defaultdict
from pathlib import Path


def parse_segments(signature: str) -> list[tuple[int, int]]:
    out: list[tuple[int, int]] = []
    s = (signature or "").strip()
    if not s:
        return out
    for chain in s.split("|"):
        chain = chain.strip()
        if not chain:
            continue
        for seg in chain.split("+"):
            seg = seg.strip()
            if "-" not in seg:
                continue
            a, b = seg.split("-", 1)
            out.append((int(a), int(b)))
    return out


def load_partition(path: Path) -> dict[tuple[str, int, int], str]:
    d: dict[tuple[str, int, int], str] = {}
    with path.open(newline="") as f:
        r = csv.DictReader(f, delimiter="\t")
        for row in r:
            if row.get("kind") != "partition_geometry_v1":
                continue
            k = (row["chrom"], int(row["start"]), int(row["end"]))
            d[k] = row.get("signature") or ""
    return d


def load_disambig(path: Path) -> dict[tuple[str, int, int], dict[str, str]]:
    d: dict[tuple[str, int, int], dict[str, str]] = {}
    with path.open(newline="") as f:
        r = csv.DictReader(f, delimiter="\t")
        for row in r:
            k = (row["chrom"], int(row["start"]), int(row["end"]))
            d[k] = row
    return d


def nearest_pairs(st: list[tuple[int, int]], ru: list[tuple[int, int]]) -> list[tuple[tuple[int, int], tuple[int, int]]]:
    # Greedy nearest matching for diagnostics (not used for hard truth claims).
    left = st[:]
    right = ru[:]
    pairs: list[tuple[tuple[int, int], tuple[int, int]]] = []
    used = [False] * len(right)
    for s in left:
        best_j = -1
        best_score = None
        for j, r in enumerate(right):
            if used[j]:
                continue
            score = abs(s[0] - r[0]) + abs(s[1] - r[1])
            if best_score is None or score < best_score:
                best_score = score
                best_j = j
        if best_j >= 0:
            used[best_j] = True
            pairs.append((s, right[best_j]))
    return pairs


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--partition-st", type=Path, required=True)
    ap.add_argument("--partition-ru", type=Path, required=True)
    ap.add_argument("--disambig", type=Path, required=True)
    ap.add_argument("--top-list", type=Path, default=None, help="optional TSV with columns chrom,start,end (e.g. batchC_top23.tsv)")
    ap.add_argument("-o", "--output", type=Path, required=True)
    args = ap.parse_args()

    st = load_partition(args.partition_st)
    ru = load_partition(args.partition_ru)
    ds = load_disambig(args.disambig)

    subset: set[tuple[str, int, int]] | None = None
    if args.top_list is not None:
        subset = set()
        with args.top_list.open(newline="") as f:
            r = csv.DictReader(f, delimiter="\t")
            for row in r:
                subset.add((row["chrom"], int(row["start"]), int(row["end"])))

    keys = sorted(set(st) & set(ru))
    rows: list[dict[str, str]] = []
    dstart = Counter()
    dend = Counter()
    dlen = Counter()

    for k in keys:
        if subset is not None and k not in subset:
            continue
        meta = ds.get(k)
        if not meta:
            continue
        # Batch C focus: junction already aligned.
        if meta.get("junction_exact_string") != "1":
            continue
        if (st[k] or "") == (ru[k] or ""):
            continue

        st_seg = sorted(parse_segments(st[k]))
        ru_seg = sorted(parse_segments(ru[k]))
        pairs = nearest_pairs(st_seg, ru_seg)
        for s, r in pairs:
            ds0 = r[0] - s[0]
            de0 = r[1] - s[1]
            dl0 = (r[1] - r[0]) - (s[1] - s[0])
            dstart[ds0] += 1
            dend[de0] += 1
            dlen[dl0] += 1

        rows.append(
            {
                "chrom": k[0],
                "start": str(k[1]),
                "end": str(k[2]),
                "st_segments": str(len(st_seg)),
                "ru_segments": str(len(ru_seg)),
                "pair_count": str(len(pairs)),
                "hint": meta.get("attribution_hint", ""),
                "class": meta.get("partition_mismatch_class", ""),
            }
        )

    with args.output.open("w", newline="") as f:
        w = csv.DictWriter(
            f,
            fieldnames=[
                "chrom",
                "start",
                "end",
                "st_segments",
                "ru_segments",
                "pair_count",
                "hint",
                "class",
            ],
            delimiter="\t",
        )
        w.writeheader()
        for r in rows:
            w.writerow(r)

    md = args.output.with_suffix(".md")
    with md.open("w") as f:
        f.write("# Partition Segment Drift Analysis\n\n")
        f.write(f"- Candidate bundles analyzed: **{len(rows)}**\n")
        f.write("- Filter: junction exact == 1 and partition signature mismatch.\n\n")
        f.write("## Dominant start deltas (RU - ST)\n")
        for k, v in dstart.most_common(10):
            f.write(f"- {k}: **{v}**\n")
        f.write("\n## Dominant end deltas (RU - ST)\n")
        for k, v in dend.most_common(10):
            f.write(f"- {k}: **{v}**\n")
        f.write("\n## Dominant segment-length deltas (RU - ST)\n")
        for k, v in dlen.most_common(10):
            f.write(f"- {k}: **{v}**\n")
        f.write("\n")

    print(args.output)
    print(md)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

