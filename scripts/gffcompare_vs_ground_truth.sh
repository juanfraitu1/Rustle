#!/usr/bin/env bash
# Compare an assembler GTF to the intersected NCBI reference (from build_ground_truth_gtf.sh).
#
# IMPORTANT — naming / scope:
#   Optional [chromosome] limits the reference to ALL transcripts on that chromosome from the
#   **full 80-family** ground-truth GTF. It does **not** restrict to one benchmark BED label
#   (e.g. ID_14 / LRRC37 panel). If your query was built with --gene-family-bed-label ID_14,
#   that does **not** make this comparison "ID_14-only" vs reference.
#   Use ./scripts/gffcompare_vs_ground_truth_over_bed_label.sh when the reference side should
#   be limited to regions for a single panel label.
#   Good prefix examples for THIS script: gffcmp_miniBam_chr941_vs_gt80, gffcmp_fullchr_vs_gt80
#   Avoid misleading prefixes like gffcmp_id14_only unless you used the bed-label script below.
#
# Uses gffcompare -R -Q:
#   -R  reference transcripts are restricted to those overlapping the query (Sn correction)
#   -Q  query transcripts restricted to those overlapping the reference (precision correction; novel loci dropped)
#
# Usage:
#   ./scripts/gffcompare_vs_ground_truth.sh <query.gtf> [output_prefix] [chromosome]
#
# Examples:
#   ./scripts/gffcompare_vs_ground_truth.sh my_assembly.gtf gffcmp_my
#   ./scripts/gffcompare_vs_ground_truth.sh my_assembly.gtf gffcmp_chr941_vs_gt80 NC_060941.1
#
# Requires: run ./scripts/build_ground_truth_gtf.sh first (HSA_genomic.gff ∩ ground_truth_80families.bed).
#   Canonical file: reference_HSA_genomic_x_gt80families.gffread.gtf
#   Legacy symlink:  ground_truth_80families_ref.gffread.gtf
# Override: GFFCOMPARE_REF_GTF=/path/to/ref.gtf

set -euo pipefail

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
if [[ -n "${GFFCOMPARE_REF_GTF:-}" ]]; then
  GTFGT="$GFFCOMPARE_REF_GTF"
else
  GTFGT=""
  for cand in \
    "$REPO_ROOT/reference_HSA_genomic_x_gt80families.gffread.gtf" \
    "$REPO_ROOT/ground_truth_80families_ref.gffread.gtf"; do
    if [[ -f "$cand" ]]; then
      GTFGT="$cand"
      break
    fi
  done
fi

QUERY="${1:?usage: $0 <query.gtf> [out_prefix] [chromosome]}"
PREFIX="${2:-gffcmp_vs_gt}"
CHROM="${3:-}"

if ! command -v gffcompare >/dev/null 2>&1; then
  echo "error: gffcompare not in PATH" >&2
  exit 1
fi
if [[ -z "$GTFGT" || ! -f "$GTFGT" ]]; then
  echo "error: reference GTF not found (expected reference_HSA_genomic_x_gt80families.gffread.gtf or symlink)." >&2
  echo "  Build: ./scripts/build_ground_truth_gtf.sh" >&2
  echo "  Or:   export GFFCOMPARE_REF_GTF=/path/to/ref.gtf" >&2
  exit 1
fi
if [[ ! -f "$QUERY" ]]; then
  echo "error: query GTF not found: $QUERY" >&2
  exit 1
fi

REF_TMP="${PREFIX}.ref_for_compare.gtf"
if [[ -n "$CHROM" ]]; then
  grep "^${CHROM}" "$GTFGT" >"$REF_TMP" || true
  echo "Using reference lines for ${CHROM} only ($(wc -l <"$REF_TMP") lines) -> $REF_TMP"
  echo "Scope: chr-wide slice of 80-family GTF (not limited to one BED label such as ID_14)."
else
  cp "$GTFGT" "$REF_TMP"
  echo "Using full reference GTF -> $REF_TMP"
  echo "Scope: entire 80-family intersect GTF."
fi

gffcompare -r "$REF_TMP" -R -Q -o "$PREFIX" "$QUERY" 2>"${PREFIX}.stderr"
echo ""
echo "Wrote: ${PREFIX}.stats (and .annotated.gtf, .tracking, ...)"
echo "gffcompare stderr: ${PREFIX}.stderr"
cat "${PREFIX}.stats"
