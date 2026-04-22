#!/usr/bin/env python3
"""
Cluster reference GTF/GFF transcripts by intron-chain signature (chrom + strand + ordered
introns as 1-based donor/acceptor pairs). Single-exon transcripts use (chrom, strand, span)
so they do not all collapse into one bucket.

Use case: Iso-Seq collapse references often contain many transcript_ids with the *same* splicing
pattern; gffcompare assigns each query row to at most one ref, inflating "absent" counts for
sibling duplicate models. Collapsing (or reporting) duplicates helps evaluation and optional
reference slimming.

Outputs (via --out-prefix PATH):
  PATH.groups.tsv            All equivalence classes (one row per key)
  PATH.duplicates.tsv         Only groups with 2+ members
  PATH.representatives.tsv     One chosen transcript_id per group + all_members column
  PATH.summary.txt             Counts and examples

Optional:
  --junction-tolerance BP   Quantize intron coords to (2*tol+1) bp buckets (0 = exact)
  --emit-keep-list PATH     One transcript_id per line (representatives only)

Usage:
  python3 scripts/collapse_reference_by_chain.py \\
    --reference GGO_clusterd.gff \\
    --out-prefix GGO_clusterd_chain_clusters

  # Fuzzy junction grouping (e.g. ±1 bp noise)
  python3 scripts/collapse_reference_by_chain.py \\
    --reference GGO_clusterd.gff \\
    --junction-tolerance 1 \\
    --out-prefix GGO_clusterd_chain_fuzzy1
"""

from __future__ import annotations

import argparse
import sys
from collections import defaultdict
from pathlib import Path
from typing import Dict, List, Set, Tuple

_scripts = Path(__file__).resolve().parent
if str(_scripts) not in sys.path:
    sys.path.insert(0, str(_scripts))

from ref_junction_support import (  # noqa: E402
    parse_gff_exons_by_transcript,
    ref_introns_from_exons,
)


def norm_strand(s: str) -> str:
    if s in ("+", "-"):
        return s
    return "+"


def quantize_introns(
    introns: List[Tuple[int, int]], tol: int
) -> Tuple[Tuple[int, int], ...]:
    if tol <= 0:
        return tuple(introns)
    w = 2 * tol + 1
    return tuple((d // w, a // w) for d, a in introns)


def pick_representative(members: List[Tuple[str, int, int, int]]) -> str:
    """
    members: (tid, span_len, n_exons, min_start) — prefer longer span, then more exons, then id.
    """
    if not members:
        return ""
    members.sort(key=lambda x: (-x[1], -x[2], x[3], x[0]))
    return members[0][0]


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--reference", required=True, type=Path)
    ap.add_argument("--out-prefix", required=True, type=Path)
    ap.add_argument(
        "--junction-tolerance",
        type=int,
        default=0,
        help="bp quantization for intron coords (0 = exact match)",
    )
    ap.add_argument(
        "--emit-keep-list",
        type=Path,
        help="Write one representative transcript_id per line",
    )
    args = ap.parse_args()

    exon_map, _ = parse_gff_exons_by_transcript(args.reference)

    # key -> list of (tid, span_len, n_exons, min_start)
    groups: Dict[Tuple, List[Tuple[str, int, int, int]]] = defaultdict(list)

    for tid, (chrom, strand, exons) in exon_map.items():
        if not exons:
            continue
        strand_n = norm_strand(strand)
        starts = [a for a, _ in exons]
        ends = [b for _, b in exons]
        lo, hi = min(starts), max(ends)
        span = hi - lo + 1
        n_ex = len(exons)
        intr = ref_introns_from_exons(exons, strand)
        qintr = quantize_introns(intr, args.junction_tolerance)
        if not qintr:
            key = ("SE", chrom, strand_n, lo, hi)
        else:
            key = ("ME", chrom, strand_n, qintr)
        groups[key].append((tid, span, n_ex, lo))

    prefix = args.out_prefix
    prefix.parent.mkdir(parents=True, exist_ok=True)

    dup_groups = {k: v for k, v in groups.items() if len(v) > 1}
    n_tx = sum(len(v) for v in groups.values())
    n_dup_tx = sum(len(v) for v in dup_groups.values())
    n_singleton = sum(1 for v in groups.values() if len(v) == 1)

    lines_summary: List[str] = [
        "=== Reference clustering by intron chain ===\n",
        f"reference: {args.reference}",
        f"junction_tolerance (quantize width): {args.junction_tolerance}",
        f"unique equivalence classes (keys): {len(groups)}",
        f"transcript models parsed: {n_tx}",
        f"classes with duplicates (≥2 members): {len(dup_groups)}",
        f"transcripts in duplicated classes: {n_dup_tx}",
        f"singleton classes: {n_singleton}",
        "",
        "Representative choice: longest genomic span, then exon count, then transcript_id.",
        "",
        "Largest duplicate groups (member count):",
    ]

    by_size = sorted(dup_groups.items(), key=lambda kv: -len(kv[1]))[:25]
    for key, members in by_size:
        rep = pick_representative(members)
        tids = [m[0] for m in members]
        preview = ", ".join(tids[:5])
        if len(tids) > 5:
            preview += f", ... (+{len(tids) - 5})"
        lines_summary.append(f"  n={len(members):5d}  rep={rep}  e.g. {preview}")

    (prefix.with_suffix(".summary.txt")).write_text("\n".join(lines_summary) + "\n")

    groups_path = prefix.with_suffix(".groups.tsv")
    dups_path = prefix.with_suffix(".duplicates.tsv")
    rep_path = prefix.with_suffix(".representatives.tsv")

    seen_keep: List[str] = []

    with groups_path.open("w") as fg, dups_path.open("w") as fd, rep_path.open("w") as fr:
        fg.write(
            "kind\tchrom\tstrand\tchain\tmember_count\trepresentative_tid\tmember_ids_semicolon\n"
        )
        fd.write(
            "kind\tchrom\tstrand\tchain\tmember_count\trepresentative_tid\tmember_ids_semicolon\n"
        )
        fr.write("representative_tid\tmember_count\tall_members_semicolon\n")

        for key in sorted(groups.keys(), key=lambda k: (-len(groups[k]), str(k))):
            members = groups[key]
            rep = pick_representative(members)
            kind, chrom, strand_n = key[0], key[1], key[2]
            if kind == "ME":
                chain = ";".join(f"{d}-{a}" for d, a in key[3])
            else:
                chain = f"{key[3]}-{key[4]}"
            tids_sem = ";".join(m[0] for m in sorted(members, key=lambda x: x[0]))
            row = (
                f"{kind}\t{chrom}\t{strand_n}\t{chain}\t{len(members)}\t{rep}\t{tids_sem}\n"
            )
            fg.write(row)
            if len(members) > 1:
                fd.write(row)
            fr.write(f"{rep}\t{len(members)}\t{tids_sem}\n")

    if args.emit_keep_list:
        args.emit_keep_list.parent.mkdir(parents=True, exist_ok=True)
        done: Set[str] = set()
        for key in sorted(groups.keys()):
            rep = pick_representative(groups[key])
            if rep not in done:
                seen_keep.append(rep)
                done.add(rep)
        args.emit_keep_list.write_text("\n".join(seen_keep) + "\n")

    print("\n".join(lines_summary))
    print(f"\nWrote {groups_path}")
    print(f"Wrote {dups_path}")
    print(f"Wrote {rep_path}")
    print(f"Wrote {prefix.with_suffix('.summary.txt')}")
    if args.emit_keep_list:
        print(f"Wrote {args.emit_keep_list} ({len(seen_keep)} representatives)")


if __name__ == "__main__":
    main()
