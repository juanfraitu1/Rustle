#!/usr/bin/env bash
# Byte-identity gate for deliverable B. FOREGROUND, serial, winloci_scratch (crash rule).
# Usage: byte_identity_gate.sh freeze   # write baseline
#        byte_identity_gate.sh check    # diff current md5s vs baseline, non-zero on mismatch
set -euo pipefail
MODE="${1:?freeze|check}"
RUSTLE=/mnt/c/Users/jfris/Desktop/Rustle
BIN=$RUSTLE/target/release
SCRATCH=/home/juanfra/winloci_scratch
BAM=$SCRATCH/GGO_mm.bam
FASTA=$SCRATCH/GGO.fasta
BASELINE=$RUSTLE/bench/mechanism/byte_identity_baseline.txt
WORK=$SCRATCH/bgate
mkdir -p "$WORK"; cd "$WORK"

# --- the fixed corpus (region-restricted, foreground) ---
"$BIN/copy_assign" --bam "$BAM" --fasta "$FASTA" \
  --region NC_073224.2:129160000-129230000 --homology-primary --min-copies 2 \
  --dump-psv --out gstm >/dev/null 2>&1
# gw_family_catalog has NO --region/--regions flag (genome-wide only by design); region-restrict
# by pre-subsetting the BAM to the gate regions with samtools (foreground, <1s, tiny output) so
# the catalog run itself stays fast and serial.
GATE_BAM="$WORK/gate_regions.bam"
if [ ! -f "$GATE_BAM" ] || [ "$BAM" -nt "$GATE_BAM" ]; then
  samtools view -b "$BAM" $(cat "$RUSTLE/bench/mechanism/gate_regions.txt") -o "$GATE_BAM"
  samtools index "$GATE_BAM"
fi
"$BIN/gw_family_catalog" --bam "$GATE_BAM" --fasta "$FASTA" --out cat >/dev/null 2>&1 || true

# --- collect md5s of the scientific outputs ---
collect(){
  for f in gstm.assignments.tsv gstm.families.tsv gstm.quant.tsv gstm.psv_cols.tsv \
           cat.families.tsv cat.copies.tsv; do
    [ -f "$f" ] && printf '%s  %s\n' "$(md5sum < "$f" | cut -d' ' -f1)" "$f"
  done
}

if [ "$MODE" = freeze ]; then
  { echo "# baseline @ $(git -C "$RUSTLE" rev-parse --short HEAD)"; collect; } > "$BASELINE"
  echo "froze $(grep -vc '^#' "$BASELINE") md5s -> $BASELINE"
else
  cur=$(collect)
  base=$(grep -v '^#' "$BASELINE")
  if [ "$cur" = "$base" ]; then echo "GATE PASS: all md5 identical to baseline"; exit 0
  else echo "GATE FAIL: md5 drift vs baseline:"; diff <(echo "$base") <(echo "$cur") || true; exit 1; fi
fi
