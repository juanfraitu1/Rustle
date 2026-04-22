#!/usr/bin/env python3
"""
Analyze gffcompare .tmap + reference GTF to bucket "missing" recovery cases and
suggest empirical assembler levers (de novo — no reference cheating).

Usage:
  python3 scripts/analyze_gffcompare_gap.py \\
    --reference GGO_clusterd.gff \\
    --tmap GGO_no_deep_cmp.GGO_no_deep.gtf.tmap \\
    --out-prefix gap_report

Requires: Python 3.8+ (stdlib only).
"""

from __future__ import annotations

import argparse
import csv
import re
from collections import Counter, defaultdict
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Set, Tuple


# gffcompare class codes (StringTie/gffcompare manual; minor version differences exist).
CLASS_HELP: Dict[str, Tuple[str, List[str]]] = {
    "=": (
        "Exact intron chain + compatible terminals (within -e).",
        ["(none — baseline target)"],
    ),
    "c": (
        "Containment: one transcript fully contains the other (length/terminal mismatch).",
        [
            "Increase boundary_tolerance / TSS-TES clustering width",
            "emit_one_extra_terminal_variant / emit_all_boundary_combos (if too few terminals)",
            "merge_by_extension for contained pairs",
            "post_chain_filter: lower min_subset_ratio slightly if valid variants are dropped",
        ],
    ),
    "j": (
        "Same intron chain but terminal exons differ (junction match, boundary mismatch).",
        [
            "Raise boundary_tolerance for chain-first clustering",
            "snap_boundaries (read-consensus snap, not reference)",
            "emit_all_boundary_combos or extra terminal variant",
            "refine_with_coverage / coverage-boundary refinement if BAM available",
        ],
    ),
    "k": (
        "Cufflinks 'chain' code: reference has fragments matching multiple queries (or vice versa).",
        [
            "Reduce duplicate isoform collapse / post_chain_filter aggressiveness",
            "Review gene_family_mode / multimapping if loci are duplicated",
        ],
    ),
    "e": (
        "Single-exon transfrag vs reference (often SE vs spliced or partial).",
        [
            "single_exon_precision_mode / min_single_exon_reads tuning",
            "TSS/TES tolerances (tss_tolerance, tes_tolerance)",
            "micro-exon / short intron handling if ref is spliced",
        ],
    ),
    "o": (
        "Exon overlap but incompatible intron chain (alternative structure).",
        [
            "preserve_alt_splice_sites + junction_tolerance tradeoff",
            "alt-splice / trie variants if implemented",
            "fragmented_chains / chain_stitch_v2 for partial support",
        ],
    ),
    "u": (
        "No reliable match (often intergenic, very low overlap, or incompatible structure).",
        [
            "sensitive / min_reads floor (singletons)",
            "strict_strand off if antisense artifacts",
            "deep_chain / stitching only if evidence exists (precision risk)",
            "Validate locus overlap: may be absent from query super-loci (-R)",
        ],
    ),
    "m": (
        "Reference skipped (multiple best matches).",
        ["Review gffcompare -M/-N; reduce redundant query isoforms"],
    ),
    "n": (
        "Novel exon/intron structure vs reference.",
        [
            "Junction discovery: junction_tolerance, min_reads on rare junctions",
            "scrub_junctions off; targeted recovery modules",
        ],
    ),
    "i": (
        "Intron retention vs reference.",
        [
            "min_intron_length (allow IR-like short introns carefully)",
            "variation_graph_sensitive / IR detection paths if available",
        ],
    ),
    "x": (
        "Exon overlap, incompatible introns.",
        [
            "Similar to o: alt splice + junction graph exploration",
            "reduce over-merge of junction clusters",
        ],
    ),
    "y": (
        "Fuzzy intron match (near splice sites).",
        [
            "junction_tolerance small non-zero; preserve_alt_splice_sites balance",
            "correct_splice_sites with genome (if not de novo evaluation)",
        ],
    ),
    "s": (
        "Match on opposite strand (antisense).",
        [
            "strict_strand false (default in sensitive); chimeric_filter review",
        ],
    ),
    "p": (
        "Polymerase run-on / rare class depending on version.",
        ["Filter RT artifacts; tune 3' end clustering"],
    ),
    "r": (
        "Repeat / low complexity (version-dependent).",
        ["Repeat masking awareness; mapq / coverage filters"],
    ),
}


@dataclass
class RefTx:
    transcript_id: str
    chrom: str
    strand: str
    start: int
    end: int
    num_exons: int


def parse_gff_reference(path: Path) -> Dict[str, RefTx]:
    """Minimal GTF/GFF parse: transcript lines + exon counts per transcript_id."""
    tid_to_exons: Dict[str, List[Tuple[int, int]]] = defaultdict(list)
    tid_meta: Dict[str, Tuple[str, str, int, int]] = {}

    with path.open() as f:
        for line in f:
            if line.startswith("#") or not line.strip():
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 9:
                continue
            chrom, source, feat, start_s, end_s, score, strand, frame, attr = parts[:9]
            if feat not in ("transcript", "exon"):
                continue
            m = re.search(r'transcript_id\s+"([^"]+)"', attr)
            if not m:
                continue
            tid = m.group(1)
            start, end = int(start_s), int(end_s)
            if feat == "transcript":
                tid_meta[tid] = (chrom, strand, start, end)
            else:
                tid_to_exons[tid].append((start, end))

    out: Dict[str, RefTx] = {}
    for tid, (chrom, strand, _ts, _te) in tid_meta.items():
        exons = tid_to_exons.get(tid, [])
        if not exons:
            # fallback: transcript coords only
            out[tid] = RefTx(tid, chrom, strand, tid_meta[tid][2], tid_meta[tid][3], 1)
        else:
            starts = [a for a, b in exons]
            ends = [b for a, b in exons]
            out[tid] = RefTx(
                tid,
                chrom,
                strand,
                min(starts),
                max(ends),
                len(exons),
            )
    return out


def load_tmap(path: Path) -> Tuple[List[dict], Set[str]]:
    rows: List[dict] = []
    ref_ids: Set[str] = set()
    with path.open() as f:
        reader = csv.DictReader(f, delimiter="\t")
        for r in reader:
            if not r.get("ref_id"):
                continue
            ref_ids.add(r["ref_id"])
            rows.append(r)
    return rows, ref_ids


def bucket_exons(n: str) -> str:
    try:
        v = int(n)
    except (TypeError, ValueError):
        return "unknown"
    return "single_exon" if v == 1 else "multi_exon"


def histogram(vals: Iterable[int], bins: List[int]) -> Counter:
    c: Counter = Counter()
    for v in vals:
        placed = False
        for hi in bins:
            if v <= hi:
                c[f"<={hi}"] += 1
                placed = True
                break
        if not placed:
            c[f">{bins[-1]}"] += 1
    return c


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--reference", required=True, type=Path, help="Reference GTF/GFF")
    ap.add_argument("--tmap", required=True, type=Path, help="gffcompare .tmap file")
    ap.add_argument("--out-prefix", type=Path, default=Path("gap_analysis"), help="Output prefix")
    args = ap.parse_args()

    ref_by_id = parse_gff_reference(args.reference)
    tmap_rows, tmap_ref_ids = load_tmap(args.tmap)

    absent_ids = set(ref_by_id.keys()) - tmap_ref_ids

    # Class counts + exon bucket
    by_class = Counter()
    cross: Dict[str, Counter] = defaultdict(Counter)
    lens_by_class: Dict[str, List[int]] = defaultdict(list)

    for r in tmap_rows:
        cc = r.get("class_code", "?").strip() or "?"
        by_class[cc] += 1
        cross[cc][bucket_exons(r.get("num_exons", ""))] += 1
        le = r.get("len", "")
        try:
            lens_by_class[cc].append(int(float(le)))
        except (TypeError, ValueError):
            pass

    out_lines: List[str] = []
    out_lines.append("=== gffcompare gap analysis (de novo guidance) ===\n")
    out_lines.append(f"Reference transcripts parsed: {len(ref_by_id)}")
    out_lines.append(f"tmap rows: {len(tmap_rows)}  unique ref_id in tmap: {len(tmap_ref_ids)}")
    out_lines.append(f"Reference transcripts with NO tmap row: {len(absent_ids)}")
    out_lines.append(
        "(Usually no genomic overlap with any query transcript under -R, or outside evaluated regions.)\n"
    )

    out_lines.append("--- Class code counts (tmap) ---")
    for c, n in by_class.most_common():
        out_lines.append(f"  {c!r}: {n}")

    out_lines.append("\n--- Cross-tab: class_code vs single_exon / multi_exon (from tmap num_exons) ---")
    for c in sorted(cross.keys()):
        ct = cross[c]
        out_lines.append(
            f"  {c!r}: single_exon={ct['single_exon']} multi_exon={ct['multi_exon']} unknown={ct['unknown']}"
        )

    out_lines.append("\n--- Suggested assembly levers by class (read-side / de novo) ---")
    for c in sorted(by_class.keys()):
        desc, levers = CLASS_HELP.get(c, ("(see gffcompare manual for this code)", ["Review release notes"]))
        out_lines.append(f"\n[{c}] {desc}")
        for lv in levers:
            out_lines.append(f"    - {lv}")

    # Absent-from-tmap characterization
    if absent_ids:
        abs_lens = [max(1, ref_by_id[i].end - ref_by_id[i].start + 1) for i in absent_ids]
        abs_exons = [ref_by_id[i].num_exons for i in absent_ids]
        out_lines.append("\n--- Reference transcripts absent from tmap (no row) ---")
        out_lines.append(f"Count: {len(absent_ids)}")
        out_lines.append(f"Length (bp) histogram (transcript span): {dict(histogram(abs_lens, [500, 2000, 5000, 10000]))}")
        exon_hist = Counter(abs_exons)
        out_lines.append(f"Exon count (top): {exon_hist.most_common(10)}")

    out_lines.append("\n--- Empirical priorities (by volume) ---")
    for c, n in by_class.most_common():
        if c == "=":
            continue
        pct = 100.0 * n / max(1, len(tmap_rows))
        out_lines.append(f"  {c!r}: {n} rows ({pct:.1f}% of tmap) — address with: {CLASS_HELP.get(c, ('',[]))[0][:80]}...")

    out_text = "\n".join(out_lines) + "\n"
    prefix = args.out_prefix
    prefix.parent.mkdir(parents=True, exist_ok=True)
    summary_path = prefix.with_suffix(".summary.txt")
    summary_path.write_text(out_text)
    print(out_text)
    print(f"Wrote {summary_path}")

    # Optional TSV for absent IDs
    absent_path = prefix.with_suffix(".absent_from_tmap.tsv")
    if absent_ids:
        with absent_path.open("w") as f:
            f.write("transcript_id\tchrom\tstrand\tstart\tend\tnum_exons\tlen_bp\n")
            for tid in sorted(absent_ids):
                r = ref_by_id[tid]
                ln = max(1, r.end - r.start + 1)
                f.write(f"{tid}\t{r.chrom}\t{r.strand}\t{r.start}\t{r.end}\t{r.num_exons}\t{ln}\n")
        print(f"Wrote {absent_path} ({len(absent_ids)} rows)")


if __name__ == "__main__":
    main()
