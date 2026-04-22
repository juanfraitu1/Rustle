#!/usr/bin/env bash
# gffcompare where the **reference** is restricted to transcripts overlapping BED rows whose
# 4th column ends with |LABEL (same rule as --gene-family-bed-label).
#
# Use output prefixes that reflect this, e.g.:
#   gffcmp_query_vs_gt80_bedID14
# not gffcmp_id14_subset_minimal (which sounds like an unclear scope).
#
# Usage:
#   ./scripts/gffcompare_vs_ground_truth_over_bed_label.sh <query.gtf> <out_prefix> <panel.bed> <LABEL> [chromosome]
#
# [chromosome] if set: first restrict the ground-truth GTF to that seqname, then intersect with BED.
#
# Requires: gffcompare, bedtools, reference_HSA_genomic_x_gt80families.gffread.gtf (from build_ground_truth_gtf.sh) or GFFCOMPARE_REF_GTF.

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

QUERY="${1:?query.gtf}"
PREFIX="${2:?out_prefix}"
BED="${3:?panel.bed}"
LABEL="${4:?LABEL e.g. ID_14}"
CHROM="${5:-}"

WORK="${TMPDIR:-/tmp}/gffc_bedlabel_$$"
mkdir -p "$WORK"
trap 'rm -rf "$WORK"' EXIT

if ! command -v gffcompare >/dev/null 2>&1; then echo "error: gffcompare not in PATH" >&2; exit 1; fi
if ! command -v bedtools >/dev/null 2>&1; then echo "error: bedtools not in PATH" >&2; exit 1; fi
if [[ -z "${GTFGT:-}" || ! -f "$GTFGT" ]]; then
  echo "error: reference GTF not found; run ./scripts/build_ground_truth_gtf.sh or set GFFCOMPARE_REF_GTF" >&2
  exit 1
fi
if [[ ! -f "$QUERY" ]]; then echo "error: query GTF not found: $QUERY" >&2; exit 1; fi
if [[ ! -f "$BED" ]]; then echo "error: BED not found: $BED" >&2; exit 1; fi

awk -F '\t' -v lbl="$LABEL" 'NF >= 4 && $4 ~ "\\|" lbl "$" { print $1 "\t" $2 "\t" $3 }' "$BED" >"$WORK/matched.bed"
N=$(wc -l <"$WORK/matched.bed" | tr -d ' ')
if [[ "$N" -eq 0 ]]; then
  echo "error: no BED rows with name suffix |${LABEL}" >&2
  exit 1
fi
sort -k1,1 -k2,2n "$WORK/matched.bed" | bedtools merge >"$WORK/merged.bed"
echo "BED label |${LABEL}: ${N} row(s) -> $(wc -l <"$WORK/merged.bed" | tr -d ' ') merged interval(s)"

REF_SRC="$WORK/ref_src.gtf"
if [[ -n "$CHROM" ]]; then
  grep "^${CHROM}" "$GTFGT" >"$REF_SRC" || true
  echo "Restricted ground truth to ${CHROM} ($(wc -l <"$REF_SRC" | tr -d ' ') lines)"
else
  cp "$GTFGT" "$REF_SRC"
  echo "Using full ground-truth GTF before BED intersect ($(wc -l <"$REF_SRC" | tr -d ' ') lines)"
fi

# Do not use bedtools -sorted here: GTF mixes feature types; a single-pass chrom/start sort is not
# enough for `intersect -sorted`. Unsorted intersect is fine for panel-sized references.
grep -v '^#' "$REF_SRC" >"$WORK/ref_nohdr.gtf"
bedtools intersect -a "$WORK/ref_nohdr.gtf" -b "$WORK/merged.bed" -wa >"$WORK/ref_bedlabel.gtf" || {
  echo "error: bedtools intersect failed" >&2
  exit 1
}

REF_LINES=$(wc -l <"$WORK/ref_bedlabel.gtf" | tr -d ' ')
if [[ "$REF_LINES" -eq 0 ]]; then
  echo "error: no reference GTF lines overlap merged BED for |${LABEL}" >&2
  exit 1
fi

REF_TMP="${PREFIX}.ref_for_compare.bedlabel_${LABEL}.gtf"
cp "$WORK/ref_bedlabel.gtf" "$REF_TMP"
echo "Reference for gffcompare: ${REF_LINES} lines (transcripts/features overlapping BED |${LABEL}) -> $REF_TMP"
echo "Scope: **BED-label ${LABEL}** slice of ground truth (then -R -Q as usual)."

gffcompare -r "$REF_TMP" -R -Q -o "$PREFIX" "$QUERY" 2>"${PREFIX}.stderr"
echo ""
echo "Wrote: ${PREFIX}.stats (and .annotated.gtf, .tracking, ...)"
echo "gffcompare stderr: ${PREFIX}.stderr"
cat "${PREFIX}.stats"
