#!/bin/bash
# Extract all BAM alignments overlapping a BED (for testing assemblies on validated regions).
# Requires: samtools, sorted/indexed input BAM (.bai).
#
# Usage: ./extract_mini_bam_from_bed.sh <input.bam> <regions.bed> <output_mini.bam> [threads]
#
# Example:
#   ./extract_mini_bam_from_bed.sh A119b.bam ground_truth_80families.bed A119b_gt80_mini.bam 8

set -euo pipefail

BAM="${1:-}"
BED="${2:-}"
OUT="${3:-}"
THREADS="${4:-4}"

if [[ -z "$BAM" || -z "$BED" || -z "$OUT" ]]; then
    echo "Usage: $0 <input.bam> <regions.bed> <output_mini.bam> [threads]"
    exit 1
fi

if [[ ! -f "$BAM" ]]; then
    echo "Missing BAM: $BAM"
    exit 1
fi
if [[ ! -f "${BAM}.bai" ]]; then
    echo "Need index: ${BAM}.bai"
    exit 1
fi
if [[ ! -f "$BED" ]]; then
    echo "Missing BED: $BED"
    exit 1
fi

SORTED_BED="${BED%.bed}.sorted.bed"
TMP="${OUT%.bam}.unsorted.bam"

sort -k1,1 -k2,2n "$BED" > "$SORTED_BED"

echo "Extracting alignments overlapping $SORTED_BED -> $OUT"
samtools view -b -@"$THREADS" -L "$SORTED_BED" "$BAM" -o "$TMP"
samtools sort -@"$THREADS" -o "$OUT" "$TMP"
rm -f "$TMP"
samtools index "$OUT"

echo "---"
samtools flagstat "$OUT"
