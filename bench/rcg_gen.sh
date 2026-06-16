#!/usr/bin/env bash
# Task 7 generation: GATED read-coherence (`--read-coherence`) genome-wide, per-chrom (SERIAL,
# OOM-safe), reusing cached StringTie guides /tmp/gw/st_$C.gtf. Mirrors rc_gen.sh but runs the
# gated layer (degrade-collapse + realness gate) instead of the ungated `--read-chain`.
# Per-chrom FASTA is sliced so the gate's genome load stays tiny (vs loading all of GGO.fasta).
# Output: /tmp/gw/rcg_$C.gtf  (gated read-coherence over the -G flow baseline gd_$C.gtf).
set -uo pipefail
BAM=/mnt/c/Users/jfris/Desktop/GGO.bam
FASTA=/mnt/c/Users/jfris/Desktop/GGO.fasta
FAI=/mnt/c/Users/jfris/Desktop/GGO.fasta.fai
RUSTLE=/mnt/c/Users/jfris/Desktop/Rustle/target/release/rustle
SAMTOOLS=/home/juanfra/miniforge3/bin/samtools
OUT=/tmp/gw
CONTIGS=$(awk '$2 > 1000000 {print $1}' "$FAI")
JOB_START=$(date +%s)
for C in $CONTIGS; do
  if [[ -f "$OUT/rcg_$C.done" ]]; then echo "[$C] rcg done (skip)"; continue; fi
  if [[ ! -s "$OUT/st_$C.gtf" ]]; then echo "[$C] no st guide, skip"; continue; fi
  echo "===== [$C] RCG START $(date '+%H:%M:%S') ====="
  CS=$(date +%s)
  "$SAMTOOLS" view -b "$BAM" "$C" -o "$OUT/$C.rcg.bam" 2>/dev/null
  [[ -s "$OUT/$C.rcg.bam" ]] || { echo "[$C] slice fail"; continue; }
  "$SAMTOOLS" index "$OUT/$C.rcg.bam"
  "$SAMTOOLS" faidx "$FASTA" "$C" > "$OUT/$C.rcg.fa" 2>/dev/null
  "$SAMTOOLS" faidx "$OUT/$C.rcg.fa"
  RAYON_NUM_THREADS=4 "$RUSTLE" -L "$OUT/$C.rcg.bam" -G "$OUT/st_$C.gtf" \
      --read-coherence --genome-fasta "$OUT/$C.rcg.fa" -o "$OUT/rcg_$C.gtf" 2>"$OUT/rcg_$C.log"
  RC=$?
  rm -f "$OUT/$C.rcg.bam" "$OUT/$C.rcg.bam.bai" "$OUT/$C.rcg.fa" "$OUT/$C.rcg.fa.fai"
  [[ $RC -eq 0 && -s "$OUT/rcg_$C.gtf" ]] || { echo "[$C] rustle fail rc=$RC"; continue; }
  TX=$(grep -cP '\ttranscript\t' "$OUT/rcg_$C.gtf")
  RCTAG=$(grep -cP 'source "read_coherence"' "$OUT/rcg_$C.gtf" 2>/dev/null || echo 0)
  GATED=$(grep -oP 'global_read_coherence_gate.*' "$OUT/rcg_$C.log" 2>/dev/null | tail -1)
  echo "[$C] RCG DONE wall=$(($(date +%s)-CS))s rcg_tx=$TX rc_tagged=$RCTAG  $GATED"
  touch "$OUT/rcg_$C.done"
done
echo "RCG TOTAL WALL=$(($(date +%s)-JOB_START))s"; echo "RCG_GEN_DONE"
