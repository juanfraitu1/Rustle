#!/usr/bin/env bash
# ⛔ STALE (§6s2): this script invokes a binary named `rustle`, the monolithic assembler
# deleted in 2e046730 ("retire(2/2): delete the StringTie assembler + network-flow island").
# Cargo builds no such bin — see Cargo.toml [[bin]]. It has not been runnable since that commit.
# Kept because documents cite it as the provenance of recorded numbers; DO NOT expect it to run.
# The current pipeline entry point is `copy_assign` (its CLI is NOT a drop-in for rustle's).
# Generate read-coherence (#1) transcriptomes genome-wide: rustle -G st --read-chain,
# per-chrom (SERIAL, OOM-safe), reusing cached StringTie guides /tmp/gw/st_$C.gtf.
# Output: /tmp/gw/rc_$C.gtf (additive read-coherence over the -G flow baseline gd_$C.gtf).
set -uo pipefail
BAM=/mnt/c/Users/jfris/Desktop/GGO.bam
FAI=/mnt/c/Users/jfris/Desktop/GGO.fasta.fai
RUSTLE=/mnt/c/Users/jfris/Desktop/Rustle/target/release/rustle
SAMTOOLS=/home/juanfra/miniforge3/bin/samtools
OUT=/tmp/gw
CONTIGS=$(awk '$2 > 1000000 {print $1}' "$FAI")
JOB_START=$(date +%s)
for C in $CONTIGS; do
  if [[ -f "$OUT/rc_$C.done" ]]; then echo "[$C] rc done (skip)"; continue; fi
  if [[ ! -s "$OUT/st_$C.gtf" ]]; then echo "[$C] no st guide, skip"; continue; fi
  echo "===== [$C] RC START $(date '+%H:%M:%S') ====="
  CS=$(date +%s)
  "$SAMTOOLS" view -b "$BAM" "$C" -o "$OUT/$C.rc.bam" 2>/dev/null
  [[ -s "$OUT/$C.rc.bam" ]] || { echo "[$C] slice fail"; continue; }
  "$SAMTOOLS" index "$OUT/$C.rc.bam"
  RAYON_NUM_THREADS=4 "$RUSTLE" -L "$OUT/$C.rc.bam" -G "$OUT/st_$C.gtf" --read-chain -o "$OUT/rc_$C.gtf" 2>"$OUT/rc_$C.log"
  RC=$?
  rm -f "$OUT/$C.rc.bam" "$OUT/$C.rc.bam.bai"
  [[ $RC -eq 0 && -s "$OUT/rc_$C.gtf" ]] || { echo "[$C] rustle fail rc=$RC"; continue; }
  echo "[$C] RC DONE wall=$(($(date +%s)-CS))s rc_tx=$(grep -cP '\ttranscript\t' "$OUT/rc_$C.gtf")"
  touch "$OUT/rc_$C.done"
done
echo "RC TOTAL WALL=$(($(date +%s)-JOB_START))s"; echo "RC_GEN_DONE"
