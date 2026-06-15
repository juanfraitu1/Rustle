#!/usr/bin/env bash
# Genome-wide StringTie-EXACT floor harness (SERIAL, resumable).
#
# Establishes the baseline FLOOR the user asked for: rustle in eonly-guided mode
# (`-e -G st.gtf`) must reproduce StringTie's transcript set EXACTLY (same count,
# transcript-level Sn/Pr = 100/100). This is the anchor for the floor-then-add
# framing: VG (`--vg --vg-layer2 -G st.gtf`) starts from this exact floor and only
# ADDS (the assembler "better decisions" + the VG-layer recoveries).
#
# Per chrom (reusing st_$C.gtf from gw_threeway.sh):
#   - run `rustle -L bam -e -G st_$C.gtf -o ee_$C.gtf`  (eonly: emit ONLY ref matches)
#   - gffcompare ee vs StringTie's own GTF  -> parity (want Sn=Pr=100, equal counts)
#   - gffcompare ee vs NCBI                 -> floor's NCBI match (want == StringTie's)
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
  if [[ -f "$OUT/ee_$C.done" ]]; then echo "[$C] ee already done (skip)"; continue; fi
  if [[ ! -s "$OUT/st_$C.gtf" ]]; then echo "[$C] no st_ (no reads / not in 3-way), skip"; continue; fi
  echo "===== [$C] EE START $(date '+%H:%M:%S') ====="
  C_START=$(date +%s)

  "$SAMTOOLS" view -b "$BAM" "$C" -o "$OUT/$C.bam" 2>"$OUT/$C.eeslice.log"
  if [[ $? -ne 0 || ! -s "$OUT/$C.bam" ]]; then echo "[$C] FAIL: bam slice"; FAILED+=("$C:slice"); rm -f "$OUT/$C.bam"; continue; fi
  "$SAMTOOLS" index "$OUT/$C.bam" 2>>"$OUT/$C.eeslice.log"

  # rustle EONLY guided: -e -G st.gtf, NO --vg (emit ONLY transcripts matching the guide).
  RAYON_NUM_THREADS=4 "$RUSTLE" -L "$OUT/$C.bam" -e -G "$OUT/st_$C.gtf" -o "$OUT/ee_$C.gtf" 2>"$OUT/ee_$C.log"
  RC=$?
  rm -f "$OUT/$C.bam" "$OUT/$C.bam.bai"
  if [[ $RC -ne 0 || ! -s "$OUT/ee_$C.gtf" ]]; then echo "[$C] FAIL: rustle-eonly (rc=$RC)"; FAILED+=("$C:rustle_ee"); continue; fi

  # Parity: ee vs StringTie's own GTF (the floor must equal StringTie).
  "$GFFCMP" -r "$OUT/st_$C.gtf"  -o "$OUT/gc_ee_par_$C" "$OUT/ee_$C.gtf" 2>"$OUT/gc_ee_par_$C.log"
  # Floor's NCBI match (should equal StringTie's NCBI-matched count).
  "$GFFCMP" -r "$OUT/ref_$C.gff3" -o "$OUT/gc_ee_$C"     "$OUT/ee_$C.gtf" 2>"$OUT/gc_ee_$C.log"

  EE_TX=$(grep -c $'\ttranscript\t' "$OUT/ee_$C.gtf" 2>/dev/null)
  ST_TX=$(grep -c $'\ttranscript\t' "$OUT/st_$C.gtf" 2>/dev/null)
  echo "[$C] EE DONE wall=$(($(date +%s)-C_START))s ee_tx=${EE_TX} st_tx=${ST_TX}"
  touch "$OUT/ee_$C.done"
done

echo "===== EE TOTAL WALL = $(($(date +%s)-JOB_START))s ====="
if [[ ${#FAILED[@]} -gt 0 ]]; then echo "EE FAILED CHROMS: ${FAILED[*]}"; else echo "EE FAILED CHROMS: none"; fi
echo "GW_EONLY_PARITY_ALL_DONE"
