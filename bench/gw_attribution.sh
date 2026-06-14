#!/usr/bin/env bash
# 4-way attribution harness (SERIAL, resumable). Adds a rustle GUIDED-BASELINE
# (rustle -G st.gtf, NO --vg) per chromosome to the three-way scorecard, so the
# VG-only-vs-StringTie recoveries can be split into:
#   - "assembler"  = also found by the guided baseline (rustle's assembler beats ST)
#   - "VG-layer"   = found ONLY with --vg --vg-layer2 (the secondary-alignment mechanism)
# Reuses st_$C.gtf, vg_$C.gtf, ref_$C.gff3 produced by gw_threeway.sh.
# Run AFTER gw_threeway.sh. Do NOT run chroms in parallel (OOM). Run to completion.
set -uo pipefail

BAM=/mnt/c/Users/jfris/Desktop/GGO.bam
FAI=/mnt/c/Users/jfris/Desktop/GGO.fasta.fai
RUSTLE=/mnt/c/Users/jfris/Desktop/Rustle/target/release/rustle
GFFCMP=/home/juanfra/miniforge3/bin/gffcompare
SAMTOOLS=/home/juanfra/miniforge3/bin/samtools
OUT=/tmp/gw

FAILED=()
CONTIGS=$(awk '$2 > 1000000 {print $1}' "$FAI")
JOB_START=$(date +%s)

for C in $CONTIGS; do
  if [[ -f "$OUT/gd_$C.done" ]]; then echo "[$C] gd already done (skip)"; continue; fi
  # Only chroms that the three-way already processed (st_ exists = had reads).
  if [[ ! -s "$OUT/st_$C.gtf" ]]; then echo "[$C] no st_ (no reads / not in 3-way), skip"; continue; fi
  echo "===== [$C] GD START $(date '+%H:%M:%S') ====="
  C_START=$(date +%s)

  "$SAMTOOLS" view -b "$BAM" "$C" -o "$OUT/$C.bam" 2>"$OUT/$C.gdslice.log"
  if [[ $? -ne 0 || ! -s "$OUT/$C.bam" ]]; then echo "[$C] FAIL: bam slice"; FAILED+=("$C:slice"); rm -f "$OUT/$C.bam"; continue; fi
  "$SAMTOOLS" index "$OUT/$C.bam" 2>>"$OUT/$C.gdslice.log"

  # rustle GUIDED baseline: -G st.gtf, NO --vg (no genome-fasta needed).
  RAYON_NUM_THREADS=4 "$RUSTLE" -L "$OUT/$C.bam" -G "$OUT/st_$C.gtf" -o "$OUT/gd_$C.gtf" 2>"$OUT/gd_$C.log"
  RC=$?
  rm -f "$OUT/$C.bam" "$OUT/$C.bam.bai"
  if [[ $RC -ne 0 || ! -s "$OUT/gd_$C.gtf" ]]; then echo "[$C] FAIL: rustle-guided (rc=$RC)"; FAILED+=("$C:rustle_gd"); continue; fi

  "$GFFCMP" -r "$OUT/ref_$C.gff3" -o "$OUT/gc_gd_$C" "$OUT/gd_$C.gtf" 2>"$OUT/gc_gd_$C.log"
  GD_TX=$(grep -c $'\ttranscript\t' "$OUT/gd_$C.gtf" 2>/dev/null)
  echo "[$C] GD DONE wall=$(($(date +%s)-C_START))s gd_tx=${GD_TX}"
  touch "$OUT/gd_$C.done"
done

echo "===== GD TOTAL WALL = $(($(date +%s)-JOB_START))s ====="
if [[ ${#FAILED[@]} -gt 0 ]]; then echo "GD FAILED CHROMS: ${FAILED[*]}"; else echo "GD FAILED CHROMS: none"; fi
echo "GW_ATTRIBUTION_ALL_DONE"
