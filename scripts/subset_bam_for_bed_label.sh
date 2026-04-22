#!/usr/bin/env bash
# Subset a BAM to reads overlapping BED intervals whose 4th column ends with |LABEL
# (same rule as --gene-family-bed-label, e.g. ID_14 does not match ID_141).
#
# Usage:
#   ./scripts/subset_bam_for_bed_label.sh <in.bam> <panel.bed> <LABEL> <out.bam> [slop_bp]
#
# Example:
#   ./scripts/subset_bam_for_bed_label.sh A119b_gt80_mini.bam ground_truth_80families.bed ID_14 A119b_ID14_only.bam 200000
#
# Naming: this BAM is “reads overlapping BED |LABEL regions”, not “gffcompare ID_14-only”.
# For gffcompare vs reference restricted to those regions, use:
#   ./scripts/gffcompare_vs_ground_truth_over_bed_label.sh
#
# Requires: samtools, awk; optional slop uses bedtools slop (needs chrom sizes from BAM header).

set -euo pipefail

IN_BAM="${1:?in.bam}"
BED="${2:?panel.bed}"
LABEL="${3:?LABEL e.g. ID_14}"
OUT_BAM="${4:?out.bam}"
SLOP="${5:-0}"

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
TMPDIR="${TMPDIR:-/tmp}"
STEM="$(basename "$OUT_BAM" .bam)"
WORK="$TMPDIR/subset_bed_label_${STEM}_$$"
mkdir -p "$WORK"
trap 'rm -rf "$WORK"' EXIT

if [[ ! -f "$IN_BAM" ]]; then echo "error: missing $IN_BAM" >&2; exit 1; fi
if [[ ! -f "$BED" ]]; then echo "error: missing $BED" >&2; exit 1; fi

# Strict suffix: |ID_14$ in awk (label must not contain regex specials; IDs are alphanumeric)
awk -F '\t' -v lbl="$LABEL" '
  NF >= 4 && $4 ~ "\\|" lbl "$" { print $1 "\t" $2 "\t" $3 "\t" $4 }
' "$BED" >"$WORK/matched.bed"

N=$(wc -l <"$WORK/matched.bed" | tr -d ' ')
if [[ "$N" -eq 0 ]]; then
  echo "error: no BED rows with name suffix |${LABEL}" >&2
  exit 1
fi
echo "Matched ${N} BED interval(s) for |${LABEL}"

sort -k1,1 -k2,2n "$WORK/matched.bed" >"$WORK/sorted.bed"
bedtools merge -i "$WORK/sorted.bed" >"$WORK/merged.bed"

REGIONS="$WORK/merged.bed"
if [[ "$SLOP" != "0" && "$SLOP" != "00" ]]; then
  samtools view -H "$IN_BAM" | awk '/^@SQ/ {
    sn=""; ln=""
    for (i=2;i<=NF;i++) {
      if ($i ~ /^SN:/) sn=substr($i,4)
      if ($i ~ /^LN:/) ln=substr($i,4)
    }
    if (sn != "" && ln != "") print sn "\t" ln
  }' >"$WORK/genome.txt"
  bedtools slop -i "$WORK/merged.bed" -g "$WORK/genome.txt" -b "$SLOP" >"$WORK/slopped.bed"
  REGIONS="$WORK/slopped.bed"
  echo "Applied slop ${SLOP} bp (clamped to chrom lengths from BAM header)"
fi

samtools view -b -L "$REGIONS" -o "$OUT_BAM" "$IN_BAM"
samtools index "$OUT_BAM"
echo "Wrote $OUT_BAM and ${OUT_BAM}.bai"
echo "Regions used ($(wc -l <"$REGIONS" | tr -d ' ') line(s)):"
cat "$REGIONS"
