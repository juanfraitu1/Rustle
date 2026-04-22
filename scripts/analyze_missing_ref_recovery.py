#!/usr/bin/env python3
"""
Analyze reference transcripts that still lack a gffcompare '=' (exact) association — or are
entirely absent from the .tmap — and relate them to empirical read/junction/query overlap
signals to suggest *which* assembler knobs are most likely to help which buckets.

Inputs
  --reference   Reference GFF/GTF (same as gffcompare -r)
  --tmap        gffcompare .tmap for a *consistent* query run
  Optional:
  --query-gtf   Assembled query GTF (transcript spans) for locus overlap
  --bam         Aligned BAM for span counts + junction support (requires pysam)

Outputs (prefix.*)
  .summary.txt              Human-readable tables + ranked lever hints
  .ref_buckets.tsv          One row per reference transcript + bucket + features
  .subopt_by_class.tsv     Reference transcripts in tmap but never '=' (dominant class)

Usage (GGO example — use a tmap that shares transcript_id space with reference):
  python3 scripts/analyze_missing_ref_recovery.py \\
    --reference GGO_clusterd.gff \\
    --tmap gffcmp_GGO_se_tune_default.GGO_se_tune_default.gtf.tmap \\
    --query-gtf GGO_se_tune_default.gtf \\
    --bam GGO_clusterd_aln.bam \\
    --out-prefix missing_recovery_GGO_default

Requires: Python 3.8+; pysam optional (for --bam).
"""

from __future__ import annotations

import argparse
import csv
import re
import sys
from collections import Counter, defaultdict
from pathlib import Path
from typing import Dict, List, Optional, Set, Tuple

_scripts = Path(__file__).resolve().parent
if str(_scripts) not in sys.path:
    sys.path.insert(0, str(_scripts))

from analyze_gffcompare_gap import (  # noqa: E402
    CLASS_HELP,
    RefTx,
    parse_gff_reference,
)

try:
    import pysam
except ImportError:
    pysam = None

from ref_junction_support import (  # noqa: E402
    RefJunctionStats,
    build_junction_support_counter_from_path,
    count_junction_support_for_transcript,
    parse_gff_exons_by_transcript,
    ref_introns_from_exons,
)


def load_tmap_by_ref(path: Path) -> Tuple[Dict[str, List[str]], List[dict]]:
    """ref_id -> list of class_code per tmap row involving that ref."""
    by_ref: Dict[str, List[str]] = defaultdict(list)
    rows: List[dict] = []
    with path.open() as f:
        r = csv.DictReader(f, delimiter="\t")
        for row in r:
            rid = row.get("ref_id", "").strip()
            if not rid:
                continue
            cc = (row.get("class_code") or "?").strip() or "?"
            by_ref[rid].append(cc)
            rows.append(row)
    return by_ref, rows


def ref_category(codes: List[str]) -> Tuple[str, str]:
    """
    Returns (category, detail).
    category: has_exact | subopt_only | (only for union logic)
    detail: dominant non-= code for subopt, or 'has_=' for exact
    """
    if not codes:
        return "absent", ""
    if "=" in codes:
        extras = [c for c in codes if c != "="]
        if not extras:
            return "has_exact", "="
        # has both exact and other match types
        return "has_exact", "<mixed>_" + Counter(extras).most_common(1)[0][0]
    dom = Counter(codes).most_common(1)[0][0]
    return "subopt_only", dom


def build_query_span_index(query_gtf: Path) -> Dict[str, List[Tuple[int, int]]]:
    qmap, _ = parse_gff_exons_by_transcript(query_gtf)
    by_chrom: Dict[str, List[Tuple[int, int]]] = defaultdict(list)
    for _tid, (chrom, _strand, exons) in qmap.items():
        if not exons:
            continue
        lo = min(a for a, _ in exons)
        hi = max(b for _, b in exons)
        by_chrom[chrom].append((lo, hi))
    for chrom in by_chrom:
        by_chrom[chrom].sort(key=lambda x: x[0])
    return by_chrom


def span_overlaps_query(
    by_chrom: Dict[str, List[Tuple[int, int]]], chrom: str, lo: int, hi: int
) -> bool:
    for s, e in by_chrom.get(chrom, []):
        if s > hi:
            break
        if e >= lo:
            return True
    return False


def count_span_reads(bam: "pysam.AlignmentFile", chrom: str, start: int, end: int) -> int:
    try:
        return bam.count(chrom, start - 1, end)
    except (ValueError, TypeError):
        return -1


def norm_strand(s: str) -> str:
    return s if s in ("+", "-") else "+"


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--reference", required=True, type=Path)
    ap.add_argument("--tmap", required=True, type=Path)
    ap.add_argument("--query-gtf", type=Path)
    ap.add_argument("--bam", type=Path)
    ap.add_argument("--junction-slop", type=int, default=0)
    ap.add_argument("--out-prefix", type=Path, default=Path("missing_recovery"))
    args = ap.parse_args()

    ref_by_id: Dict[str, RefTx] = parse_gff_reference(args.reference)
    exon_map, _nex = parse_gff_exons_by_transcript(args.reference)
    tmap_by_ref, tmap_rows = load_tmap_by_ref(args.tmap)

    if args.bam and pysam is None:
        raise SystemExit("Install pysam for --bam, or omit --bam")

    jctr = None
    bam = None
    if args.bam:
        jctr = build_junction_support_counter_from_path(args.bam, args.junction_slop)
        bam = pysam.AlignmentFile(str(args.bam), "rb")

    query_by_chrom: Dict[str, List[Tuple[int, int]]] = {}
    if args.query_gtf and args.query_gtf.exists():
        query_by_chrom = build_query_span_index(args.query_gtf)

    # --- Per-ref features ---
    bucket_counts: Counter[str] = Counter()
    by_recovery: Counter[str] = Counter()
    subopt_class_counts: Counter[str] = Counter()
    absent_bucket_counts: Counter[str] = Counter()
    lever_volume: Counter[str] = Counter()  # which lever text lines get weight
    bucket_examples: Dict[str, List[str]] = defaultdict(list)

    tsv_header = (
        "transcript_id\trecovery_category\tdominant_class\tbucket\t"
        "chrom\tstrand\tnum_exons\tlen_bp\tspan_reads\t"
        "num_introns\tintrons_with_read_support\tmin_junction_reads\tfull_junction_chain\t"
        "overlaps_query_span\n"
    )
    tsv_lines: List[str] = []

    subopt_refs: List[Tuple[str, str, str]] = []  # tid, dominant, codes joined

    for tid, rtx in sorted(ref_by_id.items()):
        codes = tmap_by_ref.get(tid, [])
        cat, detail = ref_category(codes)
        if cat == "has_exact" and detail == "=":
            # Fully satisfied; skip row in "missing" TSV or include with bucket OK?
            continue
        if cat == "has_exact" and detail != "=":
            # exact + subopt rows — still "recovered" exactly somewhere
            continue

        dom = Counter(codes).most_common(1)[0][0] if codes else ""
        if cat == "subopt_only":
            subopt_refs.append((tid, dom, ",".join(sorted(set(codes)))))

        # Length / exons
        ln = max(1, rtx.end - rtx.start + 1)
        oq = 0
        if query_by_chrom:
            oq = 1 if span_overlaps_query(query_by_chrom, rtx.chrom, rtx.start, rtx.end) else 0

        jstats = RefJunctionStats(0, 0, 0, True)
        reads: Optional[int] = None
        if tid in exon_map:
            chrom_e, strand_e, ex = exon_map[tid]
            introns = ref_introns_from_exons(ex, strand_e)
            if jctr is not None:
                jstats = count_junction_support_for_transcript(
                    jctr, chrom_e, norm_strand(strand_e), introns, args.junction_slop
                )
            else:
                ni = len(introns)
                jstats = RefJunctionStats(ni, 0, 0, ni == 0)
        if bam is not None:
            reads = count_span_reads(bam, rtx.chrom, rtx.start, rtx.end)

        # Assign empirical bucket + lever hints
        bucket, levers = assign_bucket_and_levers(
            category=cat,
            dominant_class=dom,
            num_exons=rtx.num_exons,
            span_reads=reads,
            jstats=jstats,
            overlaps_query=bool(oq),
            has_bam=bam is not None,
        )
        bucket_counts[bucket] += 1
        by_recovery[cat] += 1
        if cat == "absent":
            absent_bucket_counts[bucket] += 1
        elif cat == "subopt_only":
            subopt_class_counts[dom] += 1
        for lv in levers:
            lever_volume[lv] += 1
        if len(bucket_examples[bucket]) < 3:
            bucket_examples[bucket].append(tid)

        rs = "NA" if reads is None else str(reads)
        tsv_lines.append(
            f"{tid}\t{cat}\t{dom}\t{bucket}\t{rtx.chrom}\t{rtx.strand}\t{rtx.num_exons}\t{ln}\t"
            f"{rs}\t{jstats.num_introns}\t{jstats.introns_with_support}\t{jstats.min_junction_reads}\t"
            f"{int(jstats.full_chain_supported)}\t{oq}\n"
        )

    if bam is not None:
        bam.close()

    # --- Summary text ---
    lines: List[str] = []
    n_ref = len(ref_by_id)
    n_needy = sum(bucket_counts.values())
    lines.append("=== Missing / subopt reference recovery analysis ===\n")
    lines.append(f"Reference transcripts: {n_ref}")
    lines.append(f"Transcripts needing improvement (absent_from_tmap OR subopt_only): {n_needy}")
    lines.append(
        f"  (Excluded refs that already have at least one '=' row in tmap, including mixed rows.)\n"
    )
    lines.append(f"tmap data rows: {len(tmap_rows)}  refs_with_any_tmap_row: {len(tmap_by_ref)}\n")

    lines.append("--- Split: absent from tmap vs in tmap but never '=' (subopt_only) ---")
    for label, c in by_recovery.most_common():
        pct = 100.0 * c / max(1, n_needy)
        lines.append(f"  {label}: {c} ({pct:.1f}% of needy refs)")
    lines.append("")
    lines.append("Absent-from-tmap buckets (empirical, read-backed when --bam given):")
    for bkt, c in absent_bucket_counts.most_common():
        lines.append(f"  {c:6d}  {bkt}")
    lines.append("")
    lines.append("Subopt_only: dominant gffcompare class_code (one ref can appear in multiple tmap rows;")
    lines.append("we use majority class among rows lacking '='):")
    for cls, c in subopt_class_counts.most_common():
        pct = 100.0 * c / max(1, by_recovery.get("subopt_only", 0))
        ch, _lev = CLASS_HELP.get(cls, ("", []))
        lines.append(f"  {c:6d}  [{cls}] ({pct:.1f}% of subopt_only) — {ch[:90]}")
    lines.append("")

    lines.append("--- Bucket counts (all needy refs; empirical stratification) ---")
    for bkt, n in bucket_counts.most_common():
        ex = ", ".join(bucket_examples.get(bkt, []))
        lines.append(f"  {n:6d}  {bkt}")
        if ex:
            lines.append(f"          e.g. {ex}")

    lines.append("\n--- Suggested levers ranked by *number of refs* that carry that hint ---")
    lines.append("(A ref can contribute to multiple levers; treat as prioritization, not independence.)\n")
    for lv, w in lever_volume.most_common(25):
        lines.append(f"  {w:6d}  {lv}")

    lines.append("\n--- gffcompare class reference (from analyze_gffcompare_gap) ---")
    lines.append("For subopt_only refs, see dominant_class; suggested levers from class dictionary:\n")
    for code in sorted(set(c for _, c, _ in subopt_refs)):
        desc, class_levers = CLASS_HELP.get(code, ("(unknown code)", []))
        lines.append(f"  [{code}] {desc}")
        for x in class_levers[:4]:
            lines.append(f"      - {x}")

    out_prefix = args.out_prefix
    out_prefix.parent.mkdir(parents=True, exist_ok=True)
    summary_path = out_prefix.with_suffix(".summary.txt")
    summary_path.write_text("\n".join(lines) + "\n")

    tsv_path = out_prefix.with_suffix(".ref_buckets.tsv")
    with tsv_path.open("w") as f:
        f.write(tsv_header)
        f.writelines(tsv_lines)

    subopt_path = out_prefix.with_suffix(".subopt_by_class.tsv")
    with subopt_path.open("w") as f:
        f.write("transcript_id\tdominant_class\tall_codes_seen\n")
        for tid, dom, allc in sorted(subopt_refs):
            f.write(f"{tid}\t{dom}\t{allc}\n")

    print("\n".join(lines))
    print(f"\nWrote {summary_path}")
    print(f"Wrote {tsv_path} ({len(tsv_lines)} rows)")
    print(f"Wrote {subopt_path} ({len(subopt_refs)} rows)")


def assign_bucket_and_levers(
    category: str,
    dominant_class: str,
    num_exons: int,
    span_reads: Optional[int],
    jstats: RefJunctionStats,
    overlaps_query: bool,
    has_bam: bool,
) -> Tuple[str, List[str]]:
    """
    Returns (bucket_id, lever_key_strings).
    """
    levers: List[str] = []

    if category == "absent":
        if not has_bam:
            if num_exons <= 1:
                bucket = "absent_no_bam_SE"
            elif overlaps_query:
                bucket = "absent_no_bam_ME_overlap_query"
            else:
                bucket = "absent_no_bam_ME_no_overlap"
            levers.append(
                "INSTRUMENTATION: pass --bam (and ensure chrom names match GFF) for span/junction buckets"
            )
            return bucket, levers
        assert span_reads is not None
        if span_reads == 0:
            bucket = "absent_no_reads"
            levers.append("DATA: no primary reads on ref span — sequencing depth, locus filter, or wrong sample/BAM")
        elif span_reads < 0:
            bucket = "absent_bad_chrom"
            levers.append("DATA: BAM chrom name mismatch vs GFF")
        elif jstats.num_introns == 0:
            bucket = "absent_single_exon_ref"
            levers.append("ASSEMBLY: reference is single-exon — min_single_exon_reads, single_exon_cluster_tolerance, SE vs spliced (class e/c/k) at locus")
            levers.append("COMPARE: tall single-exon query isoforms often drive k/c; see GFFCOMPARE_INVESTIGATION_AND_TUNING.md")
        elif jstats.full_chain_supported:
            if overlaps_query:
                bucket = "absent_full_junction_overlap_query"
                levers.append(
                    "DIAGNOSE: run scripts/diagnose_absent_full_junction.py on this bucket — "
                    "if most rows are assembly_has_exact_chain, the intron chain already exists in "
                    "query GTF and gffcompare mapped those queries to a *different* reference isoform "
                    "(competing FSM). Assembler changes rarely help; consider collapsing redundant ref "
                    "models, alternative chain-level metrics, or accepting one-to-one tmap limits."
                )
                levers.append(
                    "If diagnosis shows assembly_fuzzy_chain_only: junction_tol / scrub_junctions / "
                    "internal splice snap (small cohort)."
                )
                levers.append(
                    "If assembly_chain_mismatch (rare in this bucket): collinearity, "
                    "emit_all_boundary_combos, deep_chain_recovery."
                )
            else:
                bucket = "absent_full_junction_no_query_span"
                levers.append(
                    "SENSITIVITY: junction evidence present but no assembled span overlap — "
                    "lower min_reads / deep_chain_recovery / fragmented_chains; "
                    "review single-exon clouds stealing reads"
                )
                levers.append("COVERAGE: try --no-refine-with-coverage if wide SE (k) dominates nearby")
        else:
            bucket = "absent_partial_junction"
            partial = jstats.introns_with_support
            total = max(1, jstats.num_introns)
            levers.append(
                f"JUNCTION_GRAPH: only {partial}/{total} introns with read junctions — "
                "junction_tolerance (small non-zero), preserve_alt_splice_sites, "
                "deep_chain / stitched recovery for weak introns"
            )
            levers.append("FILTER: min_intron_length / chimeric_filter may be stripping rare splices")
        return bucket, levers

    # subopt_only
    bucket = f"subopt_{dominant_class or 'unknown'}"
    desc, class_levers = CLASS_HELP.get(
        dominant_class, ("(see gffcompare manual)", ["Review class in tmap"])
    )
    levers.append(f"CLASS {dominant_class}: {desc[:100]}")
    for x in class_levers[:5]:
        levers.append(f"LEVER-{dominant_class}: {x}")

    if dominant_class == "j" and jstats.num_introns > 0 and jstats.full_chain_supported:
        levers.append("EMPIRICAL-j: chain fully read-supported — prioritize snap_boundaries (50–80), boundary_tolerance; careful with refine_with_coverage (can raise k)")
    if dominant_class == "c":
        levers.append("EMPIRICAL-c: merge_by_extension (multi-exon only), emit_all_boundary_combos, post_chain min_subset_ratio")
    if dominant_class == "k":
        levers.append("EMPIRICAL-k: single-exon_cluster_tolerance + skip SE merge-by-extension (default); avoid wide SE hulls")
    if dominant_class in ("o", "n", "x", "y"):
        levers.append(
            f"EMPIRICAL-{dominant_class}: alt structure / novel intron — preserve_alt_splice_sites vs over-merge; "
            "fragmented_chains / junction_tolerance tradeoff"
        )

    if (
        not overlaps_query
        and span_reads is not None
        and span_reads > 0
    ):
        levers.append("OVERLAP: no query span overlap despite reads — sensitivity at locus (not only class-specific lever)")

    if not has_bam:
        levers.append(
            "INSTRUMENTATION: repeat with --bam to separate boundary/junction issues (e.g. class j vs o)"
        )

    return bucket, levers


if __name__ == "__main__":
    main()
