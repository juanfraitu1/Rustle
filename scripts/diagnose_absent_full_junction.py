#!/usr/bin/env python3
"""
Drill into refs in bucket `absent_full_junction_overlap_query`:
  full ref intron-chain support in BAM, query GTF span overlap, but no tmap row.

Determines whether the gap is primarily:
  A) Assembly: no overlapping query transcript has the reference intron chain (exact/fuzzy)
  B) gffcompare: an overlapping query *does* have that chain (or IS the best structural match)
     but maps to a different ref_id (competing isoform at locus)

Inputs
  --reference      Ref GFF (same as gffcompare -r)
  --query-gtf      Assembled query GTF
  --tmap           gffcompare .tmap (must share ref_id / qry_id namespace with GFFs)
  Optional:
  --ref-buckets-tsv  Output from analyze_missing_ref_recovery.py (filters bucket column)
  --junction-tol     bp tolerance for intron coords [default 3]
  --max-refs         cap for speed [default all matching rows]

Output
  --out-summary      Text report with counts + recommended levers
  --out-tsv          Per-ref diagnosis rows

Usage:
  python3 scripts/diagnose_absent_full_junction.py \\
    --reference GGO_clusterd.gff \\
    --query-gtf GGO_se_tune_default.gtf \\
    --tmap gffcmp_GGO_se_tune_default.GGO_se_tune_default.gtf.tmap \\
    --ref-buckets-tsv missing_recovery_GGO_tune_default.ref_buckets.tsv \\
    --out-summary afj_diagnosis.summary.txt \\
    --out-tsv afj_diagnosis.per_ref.tsv
"""

from __future__ import annotations

import argparse
import csv
import sys
from collections import Counter, defaultdict
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
    return s if s in ("+", "-") else "+"


def chain_match(
    ref_i: List[Tuple[int, int]],
    qry_i: List[Tuple[int, int]],
    tol: int,
) -> str:
    """exact | fuzzy | len_mismatch | incompatible"""
    if len(ref_i) != len(qry_i):
        return "len_mismatch"
    if not ref_i:
        return "exact" if len(qry_i) == 0 else "len_mismatch"
    for (a, b), (c, d) in zip(ref_i, qry_i):
        if abs(a - c) > tol or abs(b - d) > tol:
            return "incompatible"
    exact = all(
        abs(a - c) == 0 and abs(b - d) == 0 for (a, b), (c, d) in zip(ref_i, qry_i)
    )
    return "exact" if exact else "fuzzy"


def spans_overlap(
    lo1: int, hi1: int, lo2: int, hi2: int
) -> bool:
    return lo2 <= hi1 and hi2 >= lo1


def load_tmap_qry_index(tmap: Path) -> Dict[str, Tuple[str, str]]:
    """qry_id -> (ref_id, class_code)"""
    out: Dict[str, Tuple[str, str]] = {}
    with tmap.open() as f:
        r = csv.DictReader(f, delimiter="\t")
        for row in r:
            qid = row.get("qry_id", "").strip()
            rid = row.get("ref_id", "").strip()
            cc = (row.get("class_code") or "").strip()
            if qid and rid:
                out[qid] = (rid, cc)
    return out


def load_bucket_refs(tsv: Path) -> Set[str]:
    want: Set[str] = set()
    with tsv.open() as f:
        r = csv.DictReader(f, delimiter="\t")
        for row in r:
            if row.get("bucket") == "absent_full_junction_overlap_query":
                tid = row.get("transcript_id", "").strip()
                if tid:
                    want.add(tid)
    return want


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--reference", required=True, type=Path)
    ap.add_argument("--query-gtf", required=True, type=Path)
    ap.add_argument("--tmap", required=True, type=Path)
    ap.add_argument("--ref-buckets-tsv", type=Path)
    ap.add_argument("--junction-tol", type=int, default=3)
    ap.add_argument("--max-refs", type=int, default=0)
    ap.add_argument("--out-summary", required=True, type=Path)
    ap.add_argument("--out-tsv", required=True, type=Path)
    args = ap.parse_args()

    ref_exon_map, _ = parse_gff_exons_by_transcript(args.reference)
    qry_exon_map, _ = parse_gff_exons_by_transcript(args.query_gtf)
    qry_index = load_tmap_qry_index(args.tmap)

    want_refs = load_bucket_refs(args.ref_buckets_tsv) if args.ref_buckets_tsv else set()

    # Without buckets TSV: derive same set from recovery script logic is heavy; require TSV or empty
    if not want_refs:
        print(
            "No refs loaded; provide --ref-buckets-tsv from analyze_missing_ref_recovery.py",
            file=sys.stderr,
        )
        args.out_summary.parent.mkdir(parents=True, exist_ok=True)
        args.out_summary.write_text("No input refs.\n")
        return

    ref_ids = sorted(want_refs)
    if args.max_refs > 0:
        ref_ids = ref_ids[: args.max_refs]

    # Index queries by chrom/strand for overlap queries — simple scan per ref (OK ~20k * 130k = slow)
    # Build chrom/strand -> list of (qid, lo, hi, introns)
    by_loc: Dict[Tuple[str, str], List[Tuple[str, int, int, List[Tuple[int, int]]]]] = defaultdict(list)
    for qid, (chrom, strand, exons) in qry_exon_map.items():
        if not exons:
            continue
        lo = min(a for a, _ in exons)
        hi = max(b for _, b in exons)
        intr = ref_introns_from_exons(exons, strand)
        by_loc[(chrom, norm_strand(strand))].append((qid, lo, hi, intr))

    for key in by_loc:
        by_loc[key].sort(key=lambda x: x[1])

    diag_counts: Counter[str] = Counter()
    competitor_counts: Counter[str] = Counter()
    lines: List[str] = []

    for tid in ref_ids:
        if tid not in ref_exon_map:
            diag_counts["ref_parse_missing"] += 1
            continue
        chrom, strand, exons = ref_exon_map[tid]
        strand_n = norm_strand(strand)
        lo = min(a for a, _ in exons)
        hi = max(b for _, b in exons)
        ref_intr = ref_introns_from_exons(exons, strand)

        candidates = by_loc.get((chrom, strand_n), [])
        overlapping: List[Tuple[str, int, int, List[Tuple[int, int]]]] = []
        if candidates:
            for qid, ql, qh, qintr in candidates:
                if ql > hi:
                    break
                if spans_overlap(lo, hi, ql, qh):
                    overlapping.append((qid, ql, qh, qintr))

        best_code = "no_overlap_queries"
        best_detail = ""
        best_comp = ""

        exact_q: List[str] = []
        fuzzy_q: List[str] = []
        stealers: Counter[str] = Counter()

        for qid, _ql, _qh, qintr in overlapping:
            if qid not in qry_index:
                continue
            maps_ref, _cc = qry_index[qid]
            stealers[maps_ref] += 1
            m = chain_match(ref_intr, qintr, args.junction_tol)
            if m == "exact":
                exact_q.append(qid)
            elif m == "fuzzy":
                fuzzy_q.append(qid)

        if exact_q:
            best_code = "assembly_has_exact_chain"
            maps = [qry_index[q][0] for q in exact_q]
            if tid in maps:
                best_detail = "inconsistent_marked_absent_check_inputs"
            else:
                best_detail = "gffcompare_pairs_query_with_other_ref"
                for ref_alt, k in Counter(maps).most_common(5):
                    competitor_counts[ref_alt] += k
        elif fuzzy_q:
            best_code = "assembly_fuzzy_chain_only"
            best_detail = "within_tol_but_not_exact"
        elif overlapping:
            best_code = "assembly_chain_mismatch"
            best_detail = f"n_overlap_queries={len(overlapping)}"
        else:
            best_code = "no_stranded_queries"
            best_detail = "check_strand_meta"

        diag_counts[best_code] += 1
        lines.append(
            f"{tid}\t{chrom}\t{strand_n}\t{len(ref_intr)}\t{best_code}\t{best_detail}\t"
            f"{len(exact_q)}\t{len(fuzzy_q)}\t{len(overlapping)}\t"
            f"{stealers.most_common(1)[0][0] if stealers else ''}\n"
        )

    # Summary
    slines: List[str] = []
    slines.append("=== absent_full_junction_overlap_query — diagnosis ===\n")
    slines.append(f"refs analyzed: {len(ref_ids)}  junction_tol_bp: {args.junction_tol}\n")
    slines.append("--- Primary situation (per ref) ---")
    for k, v in diag_counts.most_common():
        slines.append(f"  {v:6d}  {k}")

    slines.append("\n--- Interpretation ---")
    slines.append(
        "assembly_chain_mismatch: overlapping query span but no query with the ref intron chain "
        "(within tol) → prioritize assembler: collinearity, emit boundary combos, deep_chain_recovery, "
        "min_reads on rare chains, preserve_alt_splice_sites vs over-merge."
    )
    slines.append(
        "assembly_fuzzy_chain_only: nearly-matching junctions → try small junction_tolerance / "
        "scrub_junctions off / snap internal splice sites if genome-guided."
    )
    slines.append(
        "assembly_has_exact_chain + gffcompare_maps_other_ref: correct chain exists in output "
        "but gffcompare pairs those transcripts with another reference isoform (competing FSM) → "
        "reference isoform redundancy at locus, gffcompare -e / tie rules, or reduce duplicate ref models; "
        "not fixed by another boundary pass alone."
    )
    slines.append(
        "assembly_has_exact_chain (no competitor line): verify tmap/ref_id consistency for edge cases."
    )

    if competitor_counts:
        slines.append("\n--- Top ref_ids that \"own\" overlapping queries' tmap rows (competitors) ---")
        for ref_alt, v in competitor_counts.most_common(15):
            slines.append(f"  {v:6d} hits  {ref_alt}")

    txt = "\n".join(slines) + "\n"
    args.out_summary.parent.mkdir(parents=True, exist_ok=True)
    args.out_summary.write_text(txt)

    args.out_tsv.parent.mkdir(parents=True, exist_ok=True)
    with args.out_tsv.open("w") as f:
        f.write(
            "transcript_id\tchrom\tstrand\tnum_introns\tdiagnosis\tdetail\t"
            "n_exact_chain_queries\tn_fuzzy_chain_queries\tn_span_overlap_queries\t"
            "top_stolen_by_ref\n"
        )
        f.writelines(lines)

    print(txt)
    print(f"Wrote {args.out_summary}")
    print(f"Wrote {args.out_tsv}")


if __name__ == "__main__":
    main()
