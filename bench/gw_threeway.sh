#!/usr/bin/env bash
# Genome-wide three-way gffcompare scorecard harness (SERIAL, resumable).
# Per chrom: StringTie vs rustle-VG vs NCBI annotation + parity (VG vs ST).
# Do NOT run chroms in parallel (OOM). Run to completion.
set -uo pipefail

BAM=/mnt/c/Users/jfris/Desktop/GGO.bam
FASTA=/mnt/c/Users/jfris/Desktop/GGO.fasta
FAI=/mnt/c/Users/jfris/Desktop/GGO.fasta.fai
GFF=/mnt/c/Users/jfris/Desktop/GGO_genomic.gff
ST=/mnt/c/Users/jfris/Desktop/Rustle/tools/stringtie/stringtie
RUSTLE=/mnt/c/Users/jfris/Desktop/Rustle/target/release/rustle
GFFCMP=/home/juanfra/miniforge3/bin/gffcompare
SAMTOOLS=/home/juanfra/miniforge3/bin/samtools
OUT=/tmp/gw

mkdir -p "$OUT"
FAILED=()

# Build list of contigs with size > 1e6 (skip MT and tiny).
CONTIGS=$(awk '$2 > 1000000 {print $1}' "$FAI")

JOB_START=$(date +%s)

for C in $CONTIGS; do
  if [[ -f "$OUT/$C.done" ]]; then
    echo "[$C] already done (skip)"
    continue
  fi
  echo "===== [$C] START $(date '+%H:%M:%S') ====="
  C_START=$(date +%s)

  # 1. Slice BAM for this contig.
  "$SAMTOOLS" view -b "$BAM" "$C" -o "$OUT/$C.bam" 2>"$OUT/$C.slice.log"
  if [[ $? -ne 0 || ! -s "$OUT/$C.bam" ]]; then
    echo "[$C] FAIL: bam slice"
    FAILED+=("$C:slice")
    rm -f "$OUT/$C.bam"
    continue
  fi
  "$SAMTOOLS" index "$OUT/$C.bam" 2>>"$OUT/$C.slice.log"
  NREADS=$("$SAMTOOLS" idxstats "$OUT/$C.bam" 2>/dev/null | awk -v c="$C" '$1==c {print $3}')
  if [[ -z "$NREADS" || "$NREADS" -eq 0 ]]; then
    echo "[$C] SKIP: 0 reads"
    rm -f "$OUT/$C.bam" "$OUT/$C.bam.bai"
    continue
  fi

  # 2. StringTie.
  "$ST" -L -p 4 "$OUT/$C.bam" -o "$OUT/st_$C.gtf" 2>"$OUT/st_$C.log"
  if [[ $? -ne 0 || ! -s "$OUT/st_$C.gtf" ]]; then
    echo "[$C] FAIL: stringtie"
    FAILED+=("$C:stringtie")
    rm -f "$OUT/$C.bam" "$OUT/$C.bam.bai"
    continue
  fi

  # 3. rustle VG guided.
  R_START=$(date +%s)
  RAYON_NUM_THREADS=4 "$RUSTLE" -L "$OUT/$C.bam" --vg --vg-layer2 \
    -G "$OUT/st_$C.gtf" --genome-fasta "$FASTA" \
    -o "$OUT/vg_$C.gtf" 2>"$OUT/vg_$C.log"
  RC=$?
  R_END=$(date +%s)
  R_WALL=$((R_END - R_START))
  if [[ $RC -ne 0 || ! -s "$OUT/vg_$C.gtf" ]]; then
    echo "[$C] FAIL: rustle-vg (rc=$RC, ${R_WALL}s)"
    FAILED+=("$C:rustle")
    rm -f "$OUT/$C.bam" "$OUT/$C.bam.bai"
    continue
  fi
  if [[ $R_WALL -gt 360 ]]; then
    echo "[$C] !!! SPEED WARNING: rustle-VG took ${R_WALL}s (> 6 min)"
  fi

  # 4. Per-chrom NCBI ref.
  grep -P "^$C\t" "$GFF" > "$OUT/ref_$C.gff3"

  # 5. Three gffcompare runs.
  "$GFFCMP" -r "$OUT/ref_$C.gff3" -o "$OUT/gc_st_$C" "$OUT/st_$C.gtf" 2>"$OUT/gc_st_$C.log"
  "$GFFCMP" -r "$OUT/ref_$C.gff3" -o "$OUT/gc_vg_$C" "$OUT/vg_$C.gtf" 2>"$OUT/gc_vg_$C.log"
  "$GFFCMP" -r "$OUT/st_$C.gtf"   -o "$OUT/par_$C"   "$OUT/vg_$C.gtf" 2>"$OUT/par_$C.log"

  # tx counts.
  ST_TX=$(grep -c $'\ttranscript\t' "$OUT/st_$C.gtf" 2>/dev/null)
  VG_TX=$(grep -c $'\ttranscript\t' "$OUT/vg_$C.gtf" 2>/dev/null)

  # 6. Cleanup BAM, record time, sentinel.
  rm -f "$OUT/$C.bam" "$OUT/$C.bam.bai"
  C_END=$(date +%s)
  C_WALL=$((C_END - C_START))
  echo "[$C] DONE wall=${C_WALL}s (rustle=${R_WALL}s) st_tx=${ST_TX} vg_tx=${VG_TX} reads=${NREADS}"
  echo "$C_WALL $R_WALL $ST_TX $VG_TX $NREADS" > "$OUT/$C.timing"
  touch "$OUT/$C.done"
done

JOB_END=$(date +%s)
echo "===== TOTAL WALL = $((JOB_END - JOB_START))s ====="
if [[ ${#FAILED[@]} -gt 0 ]]; then
  echo "FAILED CHROMS: ${FAILED[*]}"
else
  echo "FAILED CHROMS: none"
fi
echo "GW_THREEWAY_ALL_DONE"
