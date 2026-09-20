#!/usr/bin/env bash
# Three-condition VG borrow benchmark on GOLGA8 region.
# Usage: bash bench/vg_borrow_benchmark.sh
# Requires: target/release/rustle built, gffcompare on PATH
set -euo pipefail

RUSTLE=target/release/rustle
BAM=bench/multi_copy_eval/golga8_region.bam
REF=bench/multi_copy_eval/ref_GOLGA8.gff
FASTA=/mnt/c/Users/jfris/Desktop/GGO.fasta
OUTDIR=bench/vg_borrow_benchmark

mkdir -p "$OUTDIR"/{off,legacy,enhanced}

run_condition() {
    local label=$1
    local extra_var=${2:-}          # optional VAR=value string
    echo "=== Running: $label ==="
    if [[ -n "$extra_var" ]]; then
        env "$extra_var" \
            "$RUSTLE" "$BAM" \
            --vg --vg-solver em --vg-snp \
            --genome-fasta "$FASTA" \
            -o "$OUTDIR/$label/out.gtf" \
            2>"$OUTDIR/$label/rustle.stderr"
    else
        "$RUSTLE" "$BAM" \
            --vg --vg-solver em --vg-snp \
            --genome-fasta "$FASTA" \
            -o "$OUTDIR/$label/out.gtf" \
            2>"$OUTDIR/$label/rustle.stderr"
    fi
    gffcompare -r "$REF" "$OUTDIR/$label/out.gtf" \
        -o "$OUTDIR/$label/cmp" \
        2>/dev/null
    local matching
    matching=$(grep "Matching transcripts:" "$OUTDIR/$label/cmp.stats" \
               | awk '{print $NF}')
    echo "  Matching transcripts: $matching"
    echo "$matching" > "$OUTDIR/$label/matching_transcripts.txt"
}

run_condition "off"      "RUSTLE_VG_NO_BORROW=1"
run_condition "legacy"   "RUSTLE_VG_BORROW_LEGACY=1"
run_condition "enhanced"

OFF=$(cat "$OUTDIR/off/matching_transcripts.txt")
LEG=$(cat "$OUTDIR/legacy/matching_transcripts.txt")
ENH=$(cat "$OUTDIR/enhanced/matching_transcripts.txt")

echo ""
echo "=== GOLGA8 borrow benchmark summary ==="
echo "Baseline (StringTie 3.0): 6 exact matches (historical reference)"
printf "%-12s %s\n" "Condition" "Matching transcripts"
printf "%-12s %s\n" "---------" "--------------------"
printf "%-12s %s\n" "OFF"       "$OFF"
printf "%-12s %s\n" "Legacy"    "$LEG"
printf "%-12s %s\n" "Enhanced"  "$ENH"
