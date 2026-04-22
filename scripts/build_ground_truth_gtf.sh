#!/usr/bin/env bash
# Build a reference annotation for gffcompare: **HSA_genomic.gff ∩ ground_truth_80families.bed**
# (bedtools intersect -wa: any GFF line overlapping a BED interval).
#
# Default output prefix (repo root):
#   reference_HSA_genomic_x_gt80families
# Produces:
#   reference_HSA_genomic_x_gt80families.intersected.gff
#   reference_HSA_genomic_x_gt80families.exons.gff
#   reference_HSA_genomic_x_gt80families.gffread.gtf   (if gffread is installed)
#
# For backward compatibility, when you use the default prefix the script also updates symlinks:
#   ground_truth_80families_ref.intersected.gff  -> same as .intersected.gff above
#   ground_truth_80families_ref.gffread.gtf      -> same as .gffread.gtf above
#
# Usage:
#   ./scripts/build_ground_truth_gtf.sh [HSA_genomic.gff] [ground_truth_80families.bed] [output_prefix]
#
# Requires: bedtools; optional: gffread for GTF (-T)
#
# Example:
#   ./scripts/build_ground_truth_gtf.sh
#   gffcompare -r reference_HSA_genomic_x_gt80families.gffread.gtf -o cmp your_assembly.gtf

set -euo pipefail

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
GFF="${1:-$REPO_ROOT/HSA_genomic.gff}"
BED="${2:-$REPO_ROOT/ground_truth_80families.bed}"
DEFAULT_PREFIX="$REPO_ROOT/reference_HSA_genomic_x_gt80families"
OUT_PREFIX="${3:-$DEFAULT_PREFIX}"
USE_DEFAULT_NAME=false
if [[ "${3:-}" == "" ]]; then
  USE_DEFAULT_NAME=true
fi

if ! command -v bedtools >/dev/null 2>&1; then
  echo "error: bedtools not found in PATH" >&2
  exit 1
fi

if [[ ! -f "$GFF" ]]; then
  echo "error: GFF not found: $GFF" >&2
  exit 1
fi
if [[ ! -f "$BED" ]]; then
  echo "error: BED not found: $BED" >&2
  exit 1
fi

echo "Intersecting annotation with BED (bedtools intersect -wa):"
echo "  GFF: $GFF"
echo "  BED: $BED"
echo "  -> prefix: $OUT_PREFIX"
echo ""

# Any GFF feature overlapping a BED interval (seq names must match BED chrom column).
bedtools intersect -a "$GFF" -b "$BED" -wa >"${OUT_PREFIX}.intersected.gff"

grep $'^[^\t]*\t[^\t]*\texon\t' "${OUT_PREFIX}.intersected.gff" >"${OUT_PREFIX}.exons.gff" || true

if command -v gffread >/dev/null 2>&1; then
  gffread -T -o "${OUT_PREFIX}.gffread.gtf" "${OUT_PREFIX}.intersected.gff"
else
  echo "note: gffread not in PATH; skipped ${OUT_PREFIX}.gffread.gtf" >&2
fi

if [[ "$USE_DEFAULT_NAME" == true && -f "${OUT_PREFIX}.gffread.gtf" ]]; then
  (cd "$REPO_ROOT" && ln -sf "$(basename "${OUT_PREFIX}.intersected.gff")" ground_truth_80families_ref.intersected.gff)
  (cd "$REPO_ROOT" && ln -sf "$(basename "${OUT_PREFIX}.gffread.gtf")" ground_truth_80families_ref.gffread.gtf)
  echo ""
  echo "Symlinks for legacy paths (gffcompare defaults, etc.):"
  echo "  $REPO_ROOT/ground_truth_80families_ref.intersected.gff"
  echo "  $REPO_ROOT/ground_truth_80families_ref.gffread.gtf"
fi

echo ""
echo "Wrote:"
echo "  ${OUT_PREFIX}.intersected.gff"
echo "  ${OUT_PREFIX}.exons.gff"
if [[ -f "${OUT_PREFIX}.gffread.gtf" ]]; then
  echo "  ${OUT_PREFIX}.gffread.gtf"
fi
echo ""
echo "gffcompare (reference = HSA ∩ BED, query = your assembler GTF):"
echo "  gffcompare -r ${OUT_PREFIX}.intersected.gff -o cmp_prefix your_assembly.gtf"
if [[ -f "${OUT_PREFIX}.gffread.gtf" ]]; then
  echo "  gffcompare -r ${OUT_PREFIX}.gffread.gtf -o cmp_prefix your_assembly.gtf"
fi
echo ""
echo "Or set: export GFFCOMPARE_REF_GTF=${OUT_PREFIX}.gffread.gtf"
