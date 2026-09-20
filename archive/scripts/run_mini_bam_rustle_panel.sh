#!/usr/bin/env bash
# Extract a mini-BAM from a BED panel and print rustle + gffcompare commands for fast locus work.
#
# Usage:
#   ./scripts/run_mini_bam_rustle_panel.sh <input.bam> <panel.bed> <out_mini.bam> [threads]
#
# Recommended iteration loop (one locus at a time):
#   head -1 bench/ggo19_overprediction_panel.bed > bench/my_locus.bed
#   ./scripts/extract_mini_bam_from_bed.sh GGO_19.bam bench/my_locus.bed bench/my_locus_reads.bam 8
#   samtools view -c bench/my_locus_reads.bam   # sanity-check read count
#   ./target/release/rustle bench/my_locus_reads.bam -L --chrom NC_073243.2 -o bench/my_locus_rustle.gtf -p 8
#   gffcompare -r GGO_19.gtf -R -Q -o bench/my_locus_cmp bench/my_locus_rustle.gtf
#
# Multi-row BED panels (many wide loci) can still take as long as a near–full-chromosome rustle run;
# do not pipe rustle to `tail` with `-v` (stderr stays buffered and looks hung). Log instead:
#   ./target/release/rustle ... 2>&1 | tee bench/panel_rustle.log
#
# Single-bundle debugging on the *full* BAM (still scans the chromosome; can be slow):
#   ./target/release/rustle GGO_19.bam -L -o /tmp/locus.gtf -p 8 --chrom NC_073243.2 \
#     --debug-bundle "NC_073243.2:START-END" --only-debug-bundle
#
set -euo pipefail

BAM="${1:-}"
BED="${2:-}"
OUT="${3:-}"
THREADS="${4:-8}"

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO="$(cd "$HERE/.." && pwd)"

if [[ -z "$BAM" || -z "$BED" || -z "$OUT" ]]; then
  echo "Usage: $0 <input.bam> <panel.bed> <out_mini.bam> [threads]"
  exit 1
fi

"$HERE/extract_mini_bam_from_bed.sh" "$BAM" "$BED" "$OUT" "$THREADS"

PREFIX="${OUT%.bam}"
GTF="${PREFIX}_rustle.gtf"
CMP="${PREFIX}_cmp"
# Derive chrom from first BED line if present (column 1).
CHROM="$(awk 'NF>=1 && $1 !~ /^#/ {print $1; exit}' "$BED")"
CHROM_ARG=""
if [[ -n "$CHROM" ]]; then
  CHROM_ARG="--chrom $CHROM "
fi

cat <<EOF

--- Next: run rustle on the mini BAM (add --chrom when the BAM still names a full assembly) ---
cd $REPO
./target/release/rustle $OUT -L ${CHROM_ARG}-o $GTF -p $THREADS

--- Compare to reference at the same genomic loci (recommended for mini-BAM / regional GTFs) ---
# -R: sensitivity uses only reference transcripts overlapping the query (Sn correction).
# -Q: precision uses only query transcripts overlapping some reference (drops novel-only query).
gffcompare -r GGO_19.gtf -R -Q -o $CMP $GTF

--- Rich trace for one locus (same bundle as gffcompare XLOC; trim noise with RUSTLE_TRACE_LOCUS) ---
# RUSTLE_TRACE_LOCUS=36445000-36633000 target/release/rustle ... \\
#   --parity-stage-tsv bench/locus_parity.tsv \\
#   --snapshot-jsonl bench/locus_snap.jsonl \\
#   --trace-reference GGO_19.gtf --trace-output bench/locus_trace.txt

EOF
