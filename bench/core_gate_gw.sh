#!/usr/bin/env bash
# Genome-wide contiguous-core gate measurement: run rustle --vg --vg-layer2 with the gate-trace
# (RUSTLE_VG_CORE_GATE_TRACE=1, gate OFF/default so decisions are unchanged), per-chrom SERIAL
# (OOM-safe; --vg is memory-heavy), capturing every jaccard-passing cross-copy singleton-exon
# merge pair's (jaccard, core_cov, would_gate@0.13). One gate-off pass = the full distribution +
# the would-gate count genome-wide. Per-chrom sliced FASTA keeps the genome load small.
# Output: /tmp/gw/cg_$C.trace (the [CORE_TRACE] lines) + /tmp/gw/cg_$C.gtf (gate-off baseline).
set -uo pipefail
BAM=/home/juanfra/winloci_scratch/GGO.bam
FASTA=/home/juanfra/winloci_scratch/GGO.fasta
FAI=/home/juanfra/winloci_scratch/GGO.fasta.fai
RUSTLE=/mnt/c/Users/jfris/Desktop/Rustle/target/release/rustle
SAM=/home/juanfra/miniforge3/bin/samtools
OUT=/tmp/gw
CONTIGS=$(awk '$2 > 1000000 {print $1}' "$FAI")
JOB_START=$(date +%s)
for C in $CONTIGS; do
  if [[ -f "$OUT/cg_$C.done" ]]; then echo "[$C] cg done (skip)"; continue; fi
  if [[ ! -s "$OUT/st_$C.gtf" ]]; then echo "[$C] no st guide, skip"; continue; fi
  echo "===== [$C] CG START $(date '+%H:%M:%S') ====="
  CS=$(date +%s)
  "$SAM" view -b "$BAM" "$C" -o "$OUT/$C.cg.bam" 2>/dev/null
  [[ -s "$OUT/$C.cg.bam" ]] || { echo "[$C] slice fail"; continue; }
  "$SAM" index "$OUT/$C.cg.bam"
  "$SAM" faidx "$FASTA" "$C" > "$OUT/$C.cg.fa" 2>/dev/null && "$SAM" faidx "$OUT/$C.cg.fa"
  RUSTLE_VG_CORE_GATE_TRACE=1 RAYON_NUM_THREADS=4 "$RUSTLE" --vg --vg-layer2 \
      -L "$OUT/$C.cg.bam" -G "$OUT/st_$C.gtf" --genome-fasta "$OUT/$C.cg.fa" \
      -o "$OUT/cg_$C.gtf" 2>"$OUT/cg_$C.log"
  RC=$?
  grep '\[CORE_TRACE\]' "$OUT/cg_$C.log" > "$OUT/cg_$C.trace" 2>/dev/null
  rm -f "$OUT/$C.cg.bam" "$OUT/$C.cg.bam.bai" "$OUT/$C.cg.fa" "$OUT/$C.cg.fa.fai"
  if [[ $RC -ne 0 ]]; then echo "[$C] rustle fail rc=$RC (see cg_$C.log)"; continue; fi
  echo "[$C] CG DONE wall=$(($(date +%s)-CS))s trace_pairs=$(wc -l < "$OUT/cg_$C.trace" 2>/dev/null || echo 0) gated=$(grep -c 'would_gate=true' "$OUT/cg_$C.trace" 2>/dev/null || echo 0)"
  touch "$OUT/cg_$C.done"
done
echo "CG TOTAL WALL=$(($(date +%s)-JOB_START))s"; echo "CG_GW_DONE"
