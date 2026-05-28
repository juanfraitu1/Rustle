#!/usr/bin/env bash
# Generate ML feature dump for training.
# Usage: bash bench/run_ml_dump.sh <BAM> <REF_GTF> [output_prefix]
# Outputs: <prefix>_features.jsonl, <prefix>_all_transcripts.gtf, <prefix>_cmp.* in the working directory.
set -euo pipefail

BAM=${1:?Usage: run_ml_dump.sh <BAM> <REF_GTF> [prefix]}
REF_GTF=${2:?}
PREFIX=${3:-ml_training}

RUSTLE=./target/release/rustle
GTF_OUT="${PREFIX}_all_transcripts.gtf"
CMP_OUT="${PREFIX}_cmp"

echo "[run_ml_dump] Running rustle with -f 0 and feature dump..."
RUSTLE_ML_FEATURE_DUMP=1 \
RUSTLE_ML_FEATURE_DUMP_PATH="${PREFIX}_features.jsonl" \
"$RUSTLE" -L "$BAM" -f 0 -o "$GTF_OUT" 2>/dev/null

echo "[run_ml_dump] Running gffcompare..."
gffcompare -RQ -r "$REF_GTF" "$GTF_OUT" -o "$CMP_OUT" 2>/dev/null

echo "[run_ml_dump] Done."
echo "  Features : ${PREFIX}_features.jsonl"
echo "  GTF      : $GTF_OUT"
echo "  tmap     : ${CMP_OUT}.${GTF_OUT##*/}.tmap"
