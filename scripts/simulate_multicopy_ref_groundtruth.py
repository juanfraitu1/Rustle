#!/usr/bin/env python3
"""
Simulate multi-copy / paralog-style reference groups from a GFF alone (no BAM).

Definition (splicing-level paralog hypothesis):
  - Multi-exon transcripts are keyed by (strand, tuple of intron donor/acceptor pairs).
  - Single-exon: keyed by (strand, span_len_bucket) — excluded from cross-locus chain groups by default.

Modes (--mode):
  distant   — Multi-locus paralog proxy: same chain and (different chrom OR midpoint separation
              >= --min-separation-bp on same chrom+strand). **Often zero** if the reference mostly
              holds near-identical duplicate *isoform IDs* at one locus rather than far paralogs.
  redundant — Any chain shared by >= 2 multi-exon transcripts (reference redundancy / gffcompare
              competing-FSM proxy); reports span spread.

A group is "multicopy_candidate" in **distant** mode if it passes the separation rule above.

Outputs (--out-prefix):
  PATH.summary.txt       Counts + interpretation
  PATH.groups.tsv        group_id, n_members, n_chroms, chroms..., transcript_ids

This is a *reference structural* proxy for gene families, not copy-number truth.

Usage:
  python3 scripts/simulate_multicopy_ref_groundtruth.py \\
    --reference GGO_clusterd.gff \\
    --out-prefix GGO_multicopy_sim \\
    --min-separation-bp 100000
"""

from __future__ import annotations

import argparse
import csv
from collections import defaultdict
from pathlib import Path
from typing import Dict, List, Set, Tuple


def parse_transcripts(path: Path) -> Dict[str, Tuple[str, str, int, int, List[Tuple[int, int]]]]:
    """tid -> (chrom, strand, start, end, sorted exons)"""
    import re

    meta: Dict[str, Tuple[str, str, int, int]] = {}
    exons: Dict[str, List[Tuple[int, int]]] = defaultdict(list)
    with path.open() as f:
        for line in f:
            if line.startswith("#") or not line.strip():
                continue
            p = line.rstrip("\n").split("\t")
            if len(p) < 9:
                continue
            feat = p[2]
            chrom, _, _, s, e, _, strand, _, attr = p[:9]
            m = re.search(r'transcript_id\s+"([^"]+)"', attr)
            if not m:
                continue
            tid = m.group(1)
            a, b = int(s), int(e)
            if feat == "transcript":
                meta[tid] = (chrom, strand, a, b)
            elif feat == "exon":
                exons[tid].append((a, b))
    out: Dict[str, Tuple[str, str, int, int, List[Tuple[int, int]]]] = {}
    for tid, (chrom, strand, ts, te) in meta.items():
        el = exons.get(tid, [])
        if el:
            el = sorted(el)
            lo = min(x for x, _ in el)
            hi = max(y for _, y in el)
            out[tid] = (chrom, strand, lo, hi, el)
        else:
            out[tid] = (chrom, strand, ts, te, [])
    return out


def intron_chain(ex: List[Tuple[int, int]], strand: str) -> Tuple[Tuple[int, int], ...]:
    if len(ex) < 2:
        return tuple()
    se = sorted(ex)
    introns = []
    for i in range(len(se) - 1):
        donor = se[i][1]
        acceptor = se[i + 1][0]
        introns.append((donor, acceptor))
    return tuple(introns)


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--reference", required=True, type=Path)
    ap.add_argument("--out-prefix", required=True, type=Path)
    ap.add_argument(
        "--min-separation-bp",
        type=int,
        default=100_000,
        help="Min distance between span midpoints on same chrom+strand (distant mode)",
    )
    ap.add_argument(
        "--mode",
        choices=("distant", "redundant"),
        default="distant",
        help="distant=paralog-style separation; redundant=shared-chain reference pairs (any distance)",
    )
    args = ap.parse_args()

    tx = parse_transcripts(args.reference)
    chain_to_ids: Dict[Tuple[str, Tuple[Tuple[int, int], ...]], List[str]] = defaultdict(list)
    for tid, (chrom, strand, lo, hi, ex) in tx.items():
        ch = intron_chain(ex, strand)
        if not ch:
            continue
        chain_to_ids[(strand, ch)].append(tid)

    multicopy_groups: List[Tuple[str, List[str], Set[str], str]] = []
    gid = 0
    for (strand, ch), ids in chain_to_ids.items():
        if len(ids) < 2:
            continue
        entries = []
        chroms: Set[str] = set()
        mids: List[int] = []
        for tid in ids:
            c, st, lo, hi, _ = tx[tid]
            mid = (lo + hi) // 2
            entries.append((tid, c, st, lo, hi, mid))
            chroms.add(c)
            mids.append(mid)
        mids_sorted = sorted(mids)
        spread = mids_sorted[-1] - mids_sorted[0] if mids_sorted else 0

        if args.mode == "redundant":
            gid += 1
            lbl = f"red_{gid:05d}_spread{spread}"
            multicopy_groups.append((lbl, list(ids), chroms, str(spread)))
            continue

        multilocus = len(chroms) >= 2
        if not multilocus:
            by_ch: Dict[Tuple[str, str], List[int]] = defaultdict(list)
            for _tid, c, st, _lo, _hi, mid in entries:
                by_ch[(c, st)].append(mid)
            for ms in by_ch.values():
                ms_sort = sorted(ms)
                for i in range(len(ms_sort)):
                    for j in range(i + 1, len(ms_sort)):
                        if ms_sort[j] - ms_sort[i] >= args.min_separation_bp:
                            multilocus = True
                            break
                    if multilocus:
                        break
                if multilocus:
                    break

        if not multilocus:
            continue
        gid += 1
        multicopy_groups.append((f"mc_{gid:05d}", list(ids), chroms, str(spread)))

    args.out_prefix.parent.mkdir(parents=True, exist_ok=True)
    summ = args.out_prefix.with_suffix(".summary.txt")
    tsv = args.out_prefix.with_suffix(".groups.tsv")

    with summ.open("w") as f:
        f.write("=== Simulated multi-copy candidates (reference intron-chain paralogs) ===\n\n")
        f.write(f"Reference: {args.reference}\n")
        f.write(f"Multi-exon transcripts: {sum(1 for t in tx.values() if len(t[4])>=2)}\n")
        f.write(f"Distinct intron chains (multi-exon): {len(chain_to_ids)}\n")
        f.write(f"Mode: {args.mode}\n")
        if args.mode == "distant":
            f.write(
                f"Multicopy_candidate groups (>=2 loci by chrom or >= {args.min_separation_bp} bp): {len(multicopy_groups)}\n"
            )
        else:
            f.write(f"Redundant_chain groups (>=2 ref tx sharing intron chain): {len(multicopy_groups)}\n")
        refs_in_mc = set()
        for row in multicopy_groups:
            refs_in_mc.update(row[1])
        f.write(f"Reference transcripts covered: {len(refs_in_mc)}\n")

    with tsv.open("w", newline="") as f:
        w = csv.writer(f, delimiter="\t")
        w.writerow(["group_id", "n_members", "n_chroms", "midpoint_spread", "chroms", "transcript_ids"])
        for gname, ids, chroms, spread in multicopy_groups:
            w.writerow([gname, len(ids), len(chroms), spread, ",".join(sorted(chroms)), ",".join(sorted(ids))])

    print(summ.read_text())
    print(f"Wrote {tsv}")


if __name__ == "__main__":
    main()
