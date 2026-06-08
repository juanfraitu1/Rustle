#!/usr/bin/env bash
# U2: run SQANTI3 QC for one arm's GTF against the RefSeq reference + genome.
# Usage: 20_run_sqanti3.sh <arm_name> <arm_gtf>
source "$(dirname "$0")/config.sh"
set -euo pipefail
ARM="${1:?arm name}"; GTF="${2:?arm gtf}"
OUT="$OUTDIR/sqanti3/$ARM"; mkdir -p "$OUT"

REF_GTF="$OUTDIR/ref.gtf"
if [ ! -s "$REF_GTF" ]; then
  "$GFFREAD" "$REF_GFF" -T -o "$REF_GTF"
fi

# SQANTI3 v5.5.x uses named args (--isoforms/--refGTF/--refFasta). GTF input is
# auto-detected by the .gtf extension.
mamba run -n "$SQANTI3_ENV" python "$SQANTI3_DIR/sqanti3_qc.py" \
  --isoforms "$GTF" --refGTF "$REF_GTF" --refFasta "$GENOME_FASTA" \
  -d "$OUT" -o "$ARM" --report skip --skipORF 2>&1 | tail -8

ls "$OUT/${ARM}_classification.txt"
echo "[U2] $ARM classification -> $OUT/${ARM}_classification.txt"
