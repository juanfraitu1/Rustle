#!/usr/bin/env python3
"""
Investigate reference transcripts that never appear in gffcompare .tmap (no ref_id row).

With `gffcompare -R`, references are evaluated in ref-reduced space; "absent from tmap"
usually means gffcompare did not assign that ref to any query row — often because no
query transcript's locus overlaps enough (not necessarily zero reads).

This script checks:
  (1) Primary-alignment read counts overlapping the reference transcript span (pysam).
  (2) Whether any *assembled* query transcript (GTF) overlaps that span.
  (3) Optional: splice junction support — counts of BAM junctions matching each ref intron
      (same convention as the assembler's bam_parser: CIGAR N / long D, 1-based donor/acceptor).

Requires: pysam, indexed BAM, Python 3.8+

Example:
  python3 scripts/investigate_absent_ref_bam.py \\
    --reference GGO_clusterd.gff \\
    --tmap GGO_clusterd_fixed_cmp.GGO_clusterd_assembly_fixed.gtf.tmap \\
    --bam GGO_clusterd_aln.bam \\
    --query-gtf GGO_clusterd_assembly_fixed.gtf \\
    --junction-slop 0 \\
    --emit-candidates rare_junction_backed.tsv \\
    --out absent_investigation.txt
"""

from __future__ import annotations

import argparse
import csv
import random
import re
import sys
from collections import defaultdict
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List, Set, Tuple

_scripts_dir = Path(__file__).resolve().parent
if str(_scripts_dir) not in sys.path:
    sys.path.insert(0, str(_scripts_dir))

from ref_junction_support import (
    RefJunctionStats,
    build_junction_support_counter_from_path,
    count_junction_support_for_transcript,
    ref_introns_from_exons,
)

try:
    import pysam
except ImportError as e:
    raise SystemExit("Install pysam: pip install pysam / conda install pysam") from e


@dataclass
class Interval:
    chrom: str
    start: int  # 1-based inclusive
    end: int
    strand: str
    num_exons: int


def parse_gff_transcripts(path: Path) -> Dict[str, Interval]:
    exons: Dict[str, List[Tuple[int, int]]] = defaultdict(list)
    meta: Dict[str, Tuple[str, str]] = {}

    with path.open() as f:
        for line in f:
            if line.startswith("#") or not line.strip():
                continue
            p = line.rstrip("\n").split("\t")
            if len(p) < 9:
                continue
            feat = p[2]
            if feat not in ("transcript", "exon"):
                continue
            chrom, start_s, end_s, strand = p[0], p[3], p[4], p[6]
            m = re.search(r'transcript_id\s+"([^"]+)"', p[8])
            if not m:
                continue
            tid = m.group(1)
            if feat == "transcript":
                meta[tid] = (chrom, strand)
            else:
                exons[tid].append((int(start_s), int(end_s)))

    out: Dict[str, Interval] = {}
    for tid, (chrom, strand) in meta.items():
        ev = exons.get(tid, [])
        if not ev:
            continue
        lo = min(a for a, _ in ev)
        hi = max(b for _, b in ev)
        out[tid] = Interval(chrom, lo, hi, strand, len(ev))
    return out


def parse_gff_exons_map(path: Path) -> Dict[str, Tuple[str, str, List[Tuple[int, int]]]]:
    exons: Dict[str, List[Tuple[int, int]]] = defaultdict(list)
    meta: Dict[str, Tuple[str, str]] = {}
    with path.open() as f:
        for line in f:
            if line.startswith("#") or not line.strip():
                continue
            p = line.rstrip("\n").split("\t")
            if len(p) < 9:
                continue
            if p[2] not in ("transcript", "exon"):
                continue
            chrom, strand = p[0], p[6]
            m = re.search(r'transcript_id\s+"([^"]+)"', p[8])
            if not m:
                continue
            tid = m.group(1)
            if p[2] == "transcript":
                meta[tid] = (chrom, strand)
            else:
                exons[tid].append((int(p[3]), int(p[4])))
    out: Dict[str, Tuple[str, str, List[Tuple[int, int]]]] = {}
    for tid, (chrom, strand) in meta.items():
        ev = exons.get(tid, [])
        if ev:
            out[tid] = (chrom, strand, ev)
    return out


def load_tmap_ref_ids(tmap: Path) -> Set[str]:
    s: Set[str] = set()
    with tmap.open() as f:
        r = csv.DictReader(f, delimiter="\t")
        for row in r:
            if row.get("ref_id"):
                s.add(row["ref_id"])
    return s


def build_chrom_intervals(txs: Dict[str, Interval]) -> Dict[str, List[Tuple[int, int]]]:
    """chrom -> sorted non-overlapping merge not needed; list of (start,end) 1-based."""
    by: Dict[str, List[Tuple[int, int]]] = defaultdict(list)
    for iv in txs.values():
        by[iv.chrom].append((iv.start, iv.end))
    for chrom in by:
        by[chrom].sort(key=lambda x: x[0])
    return by


def interval_overlaps_sorted(intervals: List[Tuple[int, int]], lo: int, hi: int) -> bool:
    """Intervals sorted by start; any [s,e] overlapping [lo,hi] (1-based inclusive)."""
    for s, e in intervals:
        if s > hi:
            break
        if e >= lo:
            return True
    return False


def query_overlaps_absent(
    query_by_chrom: Dict[str, List[Tuple[int, int]]], iv: Interval
) -> bool:
    return interval_overlaps_sorted(query_by_chrom.get(iv.chrom, []), iv.start, iv.end)


def count_reads_span(bam: pysam.AlignmentFile, iv: Interval) -> int:
    """Primary alignments overlapping [start,end] (1-based -> pysam  half-open 0-based)."""
    start0 = iv.start - 1
    end1 = iv.end
    try:
        return bam.count(iv.chrom, start0, end1)
    except ValueError:
        return -1


def median(xs: List[float]) -> float:
    if not xs:
        return 0.0
    ys = sorted(x for x in xs if x >= 0)
    if not ys:
        return 0.0
    m = len(ys) // 2
    if len(ys) % 2:
        return float(ys[m])
    return (ys[m - 1] + ys[m]) / 2.0


def pct_ge(xs: List[int], thr: int) -> float:
    xs = [x for x in xs if x >= 0]
    if not xs:
        return 0.0
    return 100.0 * sum(1 for x in xs if x >= thr) / len(xs)


def norm_strand(s: str) -> str:
    if s in ("+", "-"):
        return s
    return "+"


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--reference", required=True, type=Path)
    ap.add_argument("--tmap", required=True, type=Path)
    ap.add_argument("--bam", required=True, type=Path)
    ap.add_argument("--query-gtf", type=Path, help="Assembled GTF (transcript spans)")
    ap.add_argument("--out", required=True, type=Path)
    ap.add_argument("--seed", type=int, default=42)
    ap.add_argument(
        "--junction-slop",
        type=int,
        default=0,
        help="Quantize donor/acceptor coords to (2*slop+1) bp buckets when indexing (0=exact)",
    )
    ap.add_argument(
        "--skip-junction-scan",
        action="store_true",
        help="Do not scan BAM for junction support (faster if only span counts needed)",
    )
    ap.add_argument(
        "--max-absent",
        type=int,
        default=0,
        help="Subsample at most this many absent transcripts (0=all), stratified by exon bin",
    )
    ap.add_argument(
        "--emit-candidates",
        type=Path,
        help="Write TSV: absent transcripts good targets to merge from reference "
        "(full junction chain in BAM, no query GTF span overlap when --query-gtf set)",
    )
    args = ap.parse_args()
    random.seed(args.seed)

    ref_tx = parse_gff_transcripts(args.reference)
    exon_map = parse_gff_exons_map(args.reference)
    in_tmap = load_tmap_ref_ids(args.tmap)
    absent_ids = [tid for tid in ref_tx if tid not in in_tmap]
    present_ids = [tid for tid in ref_tx if tid in in_tmap]

    def exon_bin(n: int) -> str:
        if n <= 1:
            return "1"
        if n <= 4:
            return "2-4"
        if n <= 12:
            return "5-12"
        return "13+"

    if args.max_absent and args.max_absent < len(absent_ids):
        absent_ids = random.sample(absent_ids, args.max_absent)

    query_by_chrom: Dict[str, List[Tuple[int, int]]] = {}
    if args.query_gtf and args.query_gtf.exists():
        qtx = parse_gff_transcripts(args.query_gtf)
        query_by_chrom = build_chrom_intervals(qtx)

    abs_by_bin: Dict[str, List[str]] = defaultdict(list)
    pres_by_bin: Dict[str, List[str]] = defaultdict(list)
    for tid in absent_ids:
        abs_by_bin[exon_bin(ref_tx[tid].num_exons)].append(tid)
    for tid in present_ids:
        pres_by_bin[exon_bin(ref_tx[tid].num_exons)].append(tid)

    matched_present: List[str] = []
    for b, a_ids in abs_by_bin.items():
        pool = pres_by_bin.get(b, []) or present_ids
        k = len(a_ids)
        if k == 0 or not pool:
            continue
        if len(pool) >= k:
            matched_present.extend(random.sample(pool, k))
        else:
            matched_present.extend(random.choices(pool, k=k))

    jctr = None
    if not args.skip_junction_scan:
        jctr = build_junction_support_counter_from_path(args.bam, args.junction_slop)

    bam = pysam.AlignmentFile(str(args.bam), "rb")
    abs_counts: List[int] = []
    pr_counts: List[int] = []
    absent_query_hit = 0
    absent_zero_reads = 0
    absent_bad_chrom = 0
    absent_multi = 0
    absent_full_j = 0
    control_full_j: List[int] = []

    detail_lines: List[str] = []
    candidate_rows: List[str] = []

    for tid in absent_ids:
        iv = ref_tx[tid]
        c = count_reads_span(bam, iv)
        abs_counts.append(c)
        if c < 0:
            absent_bad_chrom += 1
        elif c == 0:
            absent_zero_reads += 1
        oq = False
        if query_by_chrom:
            oq = query_overlaps_absent(query_by_chrom, iv)
            if oq:
                absent_query_hit += 1

        jstats = RefJunctionStats(0, 0, 0, True)
        if jctr is not None and tid in exon_map:
            chrom, strand, ex = exon_map[tid]
            introns = ref_introns_from_exons(ex, strand)
            jstats = count_junction_support_for_transcript(
                jctr, chrom, norm_strand(strand), introns, args.junction_slop
            )
        if jstats.num_introns > 0:
            absent_multi += 1
            if jstats.full_chain_supported:
                absent_full_j += 1

        if (
            args.emit_candidates
            and jctr is not None
            and jstats.num_introns > 0
            and jstats.full_chain_supported
            and not oq
        ):
            candidate_rows.append(
                f"{tid}\t{iv.chrom}\t{iv.start}\t{iv.end}\t{iv.strand}\t{iv.num_exons}\t"
                f"{c}\t{jstats.num_introns}\t{jstats.introns_with_support}\t{jstats.min_junction_reads}\n"
            )

        detail_lines.append(
            f"{tid}\t{iv.chrom}\t{iv.start}\t{iv.end}\t{iv.strand}\t{iv.num_exons}\t{c}\t{int(oq)}\t"
            f"{jstats.num_introns}\t{jstats.introns_with_support}\t{jstats.min_junction_reads}\t{int(jstats.full_chain_supported)}\n"
        )

    for tid in matched_present:
        pr_counts.append(count_reads_span(bam, ref_tx[tid]))
        if jctr is not None and tid in exon_map:
            chrom, strand, ex = exon_map[tid]
            introns = ref_introns_from_exons(ex, strand)
            jst = count_junction_support_for_transcript(
                jctr, chrom, norm_strand(strand), introns, args.junction_slop
            )
            if jst.num_introns > 0 and jst.full_chain_supported:
                control_full_j.append(1)
            elif jst.num_introns > 0:
                control_full_j.append(0)

    bam.close()

    lines = [
        "# Absent-from-tmap vs BAM / assembly overlap",
        f"# reference transcripts in GTF: {len(ref_tx)}",
        f"# ref_id seen in tmap: {len(in_tmap)}",
        f"# absent analyzed (no tmap row for ref_id): {len(absent_ids)}",
        "",
        "## Read counts over full ref transcript span (pysam.AlignmentFile.count; primary reads in region)",
        f"absent: n={len(abs_counts)} median={median([float(x) for x in abs_counts]):.1f}",
        f"  pct >=1 reads: {pct_ge(abs_counts,1):.1f}%  >=5: {pct_ge(abs_counts,5):.1f}%  >=20: {pct_ge(abs_counts,20):.1f}%",
        f"  zero_count_or_missing_chrom: {absent_zero_reads} bad_chrom: {absent_bad_chrom}",
        f"present_control (exon-bin matched): n={len(pr_counts)} median={median([float(x) for x in pr_counts]):.1f}",
        f"  pct >=1 reads: {pct_ge(pr_counts,1):.1f}%  >=5: {pct_ge(pr_counts,5):.1f}%",
        "",
    ]
    if query_by_chrom:
        lines.append(
            "## Query GTF overlap (any assembled transcript span intersects absent ref span)"
        )
        lines.append(
            f"absent with overlap: {absent_query_hit}/{len(absent_ids)} ({100*absent_query_hit/max(1,len(absent_ids)):.1f}%)"
        )
        lines.append(
            "If overlap% is high: reads+assembly often cover locus but gffcompare still omits ref — check -R super-locus logic, strand, or ID sets."
        )
        lines.append(
            "If overlap% is low: assembler did not place isoforms there (sensitivity) despite possible reads."
        )
        lines.append("")

    if jctr is not None:
        lines.extend(
            [
                "## Junction support (BAM primary alignments vs reference introns)",
                f"  junction_slop: {args.junction_slop}",
                f"absent multi-exon: {absent_multi}",
                f"  full intron chain has >=1 matching read junction per intron: {absent_full_j} / {absent_multi} "
                f"({100*absent_full_j/max(1,absent_multi):.1f}%) of multi-exon absent",
            ]
        )
        mctrl = len(control_full_j)
        if mctrl:
            fc = sum(control_full_j)
            lines.append(
                f"present_control multi-exon with full junction chain support: {fc}/{mctrl} ({100*fc/mctrl:.1f}%)"
            )
        lines.append(
            "Interpretation: span-level read counts can be >0 from unspliced or mis-spliced reads; "
            "low full-chain support implies the reference isoform is not clearly present as spliced alignments."
        )
        lines.append("")

    if absent_ids:
        nz = sum(1 for x in abs_counts if x > 0)
        lines.append(
            f"# Among absent: {nz} ({100*nz/len(absent_ids):.1f}%) have >=1 overlapping primary read in span."
        )

    args.out.parent.mkdir(parents=True, exist_ok=True)
    summary_path = args.out
    with summary_path.open("w") as f:
        f.write("\n".join(lines))

    detail_path = args.out.with_suffix(".absent_read_counts.tsv")
    with detail_path.open("w") as f:
        f.write(
            "transcript_id\tchrom\tstart\tend\tstrand\tnum_exons\tread_count\toverlaps_query_gtf\t"
            "num_introns\tintrons_with_support\tmin_junction_reads\tfull_junction_chain\n"
        )
        f.writelines(detail_lines)

    if args.emit_candidates:
        args.emit_candidates.parent.mkdir(parents=True, exist_ok=True)
        with args.emit_candidates.open("w") as f:
            f.write(
                "transcript_id\tchrom\tstart\tend\tstrand\tnum_exons\tread_count\t"
                "num_introns\tintrons_with_support\tmin_junction_reads\n"
            )
            f.writelines(candidate_rows)
        print(f"Wrote {len(candidate_rows)} junction-backed absent rows (no query span overlap) → {args.emit_candidates}")

    print("\n".join(lines))
    print(f"Wrote {summary_path}")
    print(f"Wrote {detail_path}")


if __name__ == "__main__":
    main()
