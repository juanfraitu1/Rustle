#!/bin/bash
# CORRECTED top-up test: merge the valid SIMTOPUP reads into the FLANK-BEARING cache BAM (soto_regions.bam,
# the substrate that actually seeds -- the earlier soto_reads.bam is exon-scoped and recovers only 23% of
# real-found members, so its verdicts are invalid). Then re-detect per-chrom with the current binary.
#   member recovers here but not on real data  -> COVERAGE-RECOVERABLE (it "lacked reads"; ideal depth seeds it)
#   member still missing with ideal coverage+flank -> K=0 (identical transcript maps MAPQ-0, never seeds own locus)
# Also emits a power check (control real-found recovery) so we KNOW the test is powered before trusting verdicts.
set -u
D=/mnt/linuxdisk/home/juanfraitu/winloci_data
TOPBAM=$D/topup/soto_reads_topup.bam
CACHE=/home/juanfra/winloci_scratch/soto_cache
FA=$D/chm13v2.0.fa
BIN=/mnt/c/Users/jfris/Desktop/Rustle/target/release/gw_family_catalog
SAM=/home/juanfra/miniforge3/bin/samtools
OUT=$CACHE/topup_flank
mkdir -p "$OUT"
CHROMS="chr1 chr2 chr3 chr5 chr7 chr8 chr9 chr11 chr14 chr15 chr16 chr17 chr21 chr22"

echo "[$(date +%H:%M:%S)] extracting SIMTOPUP reads from topup BAM..."
if [ ! -f "$OUT/topup_only.bam" ]; then
  $SAM view -h "$TOPBAM" 2>/dev/null | awk '/^@/ || /^SIMTOPUP/' | $SAM view -b -o "$OUT/topup_only.bam" - 2>/dev/null
fi
echo "  SIMTOPUP reads: $($SAM view -c "$OUT/topup_only.bam" 2>/dev/null)"
echo "[$(date +%H:%M:%S)] merging cache (flank-bearing) + topup -> cache_topup.bam ..."
if [ ! -f "$OUT/cache_topup.bam.bai" ]; then
  $SAM merge -f -@ 4 "$OUT/cache_topup.bam" "$CACHE/soto_regions.bam" "$OUT/topup_only.bam" 2>/dev/null
  $SAM index "$OUT/cache_topup.bam"
fi
echo "[$(date +%H:%M:%S)] per-chrom detection on cache_topup.bam (<=4 parallel)..."
for c in $CHROMS; do
  if [ ! -f "$OUT/$c.bam.bai" ]; then
    $SAM view -b "$OUT/cache_topup.bam" "$c" -o "$OUT/$c.bam" 2>/dev/null && $SAM index "$OUT/$c.bam" 2>/dev/null
  fi
  ( timeout 1800 "$BIN" --bam "$OUT/$c.bam" --fasta "$FA" --out "$OUT/$c" > "$OUT/$c.log" 2>&1
    echo "  [$c] copies=$(($(wc -l < "$OUT/$c.copies.tsv" 2>/dev/null || echo 1)-1))" ) &
  while [ "$(jobs -r | wc -l)" -ge 4 ]; do wait -n; done
done
wait
{ echo -e "family_id\tcopy_idx\ttid\tchrom\tstart\tend\tn_exon\tstrand\tn_reads"
  for c in $CHROMS; do [ -f "$OUT/$c.copies.tsv" ] && tail -n +2 "$OUT/$c.copies.tsv"; done; } > "$OUT/cache_topup_catalog.copies.tsv"
echo "[$(date +%H:%M:%S)] DONE: $(($(wc -l < "$OUT/cache_topup_catalog.copies.tsv")-1)) copies -> $OUT/cache_topup_catalog.copies.tsv"
