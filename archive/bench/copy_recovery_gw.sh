#!/usr/bin/env bash
# Genome-wide multimapper copy-recovery generation (SERIAL, OOM-safe, resumable).
# Runs rustle with the full recovery stack (fix#1 reference-anchored PSV + fix#2 component
# split + perf-restrict) per chrom, reusing cached StringTie guides /tmp/gw/st_$C.gtf.
# Output: /tmp/gw/rec_$C.gtf. Run AFTER the perf-restrict build. Then merge + SQANTI3.
set -uo pipefail
BAM=/mnt/c/Users/jfris/Desktop/GGO.bam
FAI=/mnt/c/Users/jfris/Desktop/GGO.fasta.fai
FA=/mnt/c/Users/jfris/Desktop/GGO.fasta
RUSTLE=/mnt/c/Users/jfris/Desktop/Rustle/target/release/rustle
SAMTOOLS=/home/juanfra/miniforge3/bin/samtools
OUT=/tmp/gw
# CONTIGS env override (e.g. the 2 chroms carrying all 27 targets) → faster than genome-wide.
CONTIGS="${CONTIGS:-$(awk '$2 > 1000000 {print $1}' "$FAI")}"
JOB_START=$(date +%s)
for C in $CONTIGS; do
  if [[ -f "$OUT/rec_$C.done" ]]; then echo "[$C] rec done (skip)"; continue; fi
  if [[ ! -s "$OUT/st_$C.gtf" ]]; then echo "[$C] no st guide, skip"; continue; fi
  echo "===== [$C] REC START $(date '+%H:%M:%S') ====="
  CS=$(date +%s)
  "$SAMTOOLS" view -b "$BAM" "$C" -o "$OUT/$C.rec.bam" 2>/dev/null
  [[ -s "$OUT/$C.rec.bam" ]] || { echo "[$C] slice fail"; continue; }
  "$SAMTOOLS" index "$OUT/$C.rec.bam"
  # recovery is activated by the RUSTLE_VG_RECOVER_COPIES env var (no CLI flag exists yet).
  RAYON_NUM_THREADS=4 RUSTLE_VG_RECOVER_COPIES=1 "$RUSTLE" -L "$OUT/$C.rec.bam" \
    --vg --vg-layer2 --vg-layer2-psv-linkage \
    -G "$OUT/st_$C.gtf" --genome-fasta "$FA" -o "$OUT/rec_$C.gtf" 2>"$OUT/rec_$C.log"
  RC=$?
  rm -f "$OUT/$C.rec.bam" "$OUT/$C.rec.bam.bai"
  [[ $RC -eq 0 && -s "$OUT/rec_$C.gtf" ]] || { echo "[$C] rustle fail rc=$RC"; continue; }
  echo "[$C] REC DONE wall=$(($(date +%s)-CS))s rec_tx=$(grep -cP '\ttranscript\t' "$OUT/rec_$C.gtf")"
  touch "$OUT/rec_$C.done"
done
echo "REC TOTAL WALL=$(($(date +%s)-JOB_START))s"; echo "COPY_RECOVERY_GW_DONE"
