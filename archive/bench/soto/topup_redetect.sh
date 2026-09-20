#!/bin/bash
# Re-detect on the TOPPED-UP BAM (real Soto reads + ideal-coverage simulated top-up of each member's OWN
# transcript) with the CURRENT binary, to split the residual not-found members into:
#   COVERAGE-RECOVERABLE  -> a locus now seeds with ideal coverage (real-data miss = depth/mis-chain, fixable)
#   K=0 STRUCTURAL FLOOR  -> still no locus even with ideal coverage (topup is K=0-self-enforcing: an
#                            identical-transcript member gets identical top-up reads that still map MAPQ-0)
# Per-chrom (family-build is single-threaded; co-located siblings share a chrom so K=0 collapse is exercised).
# Outputs to winloci_scratch (NOT /tmp), default flags (= --refine, matches the definitive baseline catalog).
set -u
D=/mnt/linuxdisk/home/juanfraitu/winloci_data
TOPBAM=$D/topup/soto_reads_topup.bam
FA=$D/chm13v2.0.fa
BIN=/mnt/c/Users/jfris/Desktop/Rustle/target/release/gw_family_catalog
SAM=/home/juanfra/miniforge3/bin/samtools
OUT=/home/juanfra/winloci_scratch/soto_cache/topup_redetect
mkdir -p "$OUT"
CHROMS="chr1 chr2 chr3 chr5 chr7 chr8 chr9 chr11 chr14 chr15 chr16 chr17 chr21 chr22"
echo "[$(date +%H:%M:%S)] extracting per-chrom mini-BAMs from topup BAM..."
for c in $CHROMS; do
  if [ ! -f "$OUT/$c.bam.bai" ]; then
    $SAM view -b "$TOPBAM" "$c" -o "$OUT/$c.bam" 2>/dev/null && $SAM index "$OUT/$c.bam" 2>/dev/null
  fi
done
echo "[$(date +%H:%M:%S)] running current-binary detection per chrom (<=4 parallel)..."
for c in $CHROMS; do
  ( timeout 1800 "$BIN" --bam "$OUT/$c.bam" --fasta "$FA" --out "$OUT/$c" > "$OUT/$c.log" 2>&1
    echo "  [$c] copies=$(($(wc -l < "$OUT/$c.copies.tsv" 2>/dev/null || echo 1)-1))" ) &
  while [ "$(jobs -r | wc -l)" -ge 4 ]; do wait -n; done
done
wait
{ echo -e "family_id\tcopy_idx\ttid\tchrom\tstart\tend\tn_exon\tstrand\tn_reads"
  for c in $CHROMS; do [ -f "$OUT/$c.copies.tsv" ] && tail -n +2 "$OUT/$c.copies.tsv"; done; } > "$OUT/topup_catalog.copies.tsv"
echo "[$(date +%H:%M:%S)] DONE topup catalog: $(($(wc -l < "$OUT/topup_catalog.copies.tsv")-1)) copies -> $OUT/topup_catalog.copies.tsv"
