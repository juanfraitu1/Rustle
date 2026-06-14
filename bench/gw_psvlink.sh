#!/usr/bin/env bash
# Genome-wide PSV-linkage attribution harness (SERIAL, resumable). Adds a
# rustle --vg-layer2-psv-linkage output (pl_$C.gtf) per chromosome, on top of the
# existing three-way scorecard outputs (st_$C.gtf / vg_$C.gtf / ref_$C.gff3 from
# gw_threeway.sh). The delta pl vs vg isolates the PSV-linkage channel's marginal
# NCBI-corroborated contribution over the rest of the VG layer.
# Run AFTER gw_threeway.sh. Per-chrom serial (whole-genome OOMs). Run to completion.
set -uo pipefail

BAM=/mnt/c/Users/jfris/Desktop/GGO.bam
FASTA=/mnt/c/Users/jfris/Desktop/GGO.fasta
FAI=/mnt/c/Users/jfris/Desktop/GGO.fasta.fai
RUSTLE=/mnt/c/Users/jfris/Desktop/Rustle/target/release/rustle
GFFCMP=/home/juanfra/miniforge3/bin/gffcompare
SAMTOOLS=/home/juanfra/miniforge3/bin/samtools
OUT=/tmp/gw

FAILED=()
CONTIGS=$(awk '$2 > 1000000 {print $1}' "$FAI")
JOB_START=$(date +%s)

for C in $CONTIGS; do
  if [[ -f "$OUT/pl_$C.done" ]]; then echo "[$C] pl already done (skip)"; continue; fi
  if [[ ! -s "$OUT/st_$C.gtf" || ! -s "$OUT/ref_$C.gff3" ]]; then echo "[$C] no st_/ref_ (not in 3-way), skip"; continue; fi
  echo "===== [$C] PL START $(date '+%H:%M:%S') ====="
  C_START=$(date +%s)

  "$SAMTOOLS" view -b "$BAM" "$C" -o "$OUT/$C.bam" 2>"$OUT/$C.plslice.log"
  if [[ $? -ne 0 || ! -s "$OUT/$C.bam" ]]; then echo "[$C] FAIL: bam slice"; FAILED+=("$C:slice"); rm -f "$OUT/$C.bam"; continue; fi
  "$SAMTOOLS" index "$OUT/$C.bam" 2>>"$OUT/$C.plslice.log"

  # rustle VG + Layer2 + PSV-linkage, guided (same config as vg_$C plus the flag).
  RAYON_NUM_THREADS=4 "$RUSTLE" -L "$OUT/$C.bam" --vg --vg-layer2 --vg-layer2-psv-linkage \
    -G "$OUT/st_$C.gtf" --genome-fasta "$FASTA" -o "$OUT/pl_$C.gtf" 2>"$OUT/pl_$C.log"
  RC=$?
  rm -f "$OUT/$C.bam" "$OUT/$C.bam.bai"
  if [[ $RC -ne 0 || ! -s "$OUT/pl_$C.gtf" ]]; then echo "[$C] FAIL: rustle-psvlink (rc=$RC)"; FAILED+=("$C:rustle_pl"); continue; fi

  "$GFFCMP" -r "$OUT/ref_$C.gff3" -o "$OUT/gc_pl_$C" "$OUT/pl_$C.gtf" 2>"$OUT/gc_pl_$C.log"
  PL_TX=$(grep -c $'\ttranscript\t' "$OUT/pl_$C.gtf" 2>/dev/null)
  VG_TX=$(grep -c $'\ttranscript\t' "$OUT/vg_$C.gtf" 2>/dev/null)
  NPSV=$(grep -c 'PSV-linked' "$OUT/pl_$C.log" 2>/dev/null)
  echo "[$C] PL DONE wall=$(($(date +%s)-C_START))s pl_tx=${PL_TX} vg_tx=${VG_TX} psvlink_log_lines=${NPSV}"
  touch "$OUT/pl_$C.done"
done

echo "===== PL TOTAL WALL = $(($(date +%s)-JOB_START))s ====="
if [[ ${#FAILED[@]} -gt 0 ]]; then echo "PL FAILED CHROMS: ${FAILED[*]}"; else echo "PL FAILED CHROMS: none"; fi
echo "GW_PSVLINK_ALL_DONE"
