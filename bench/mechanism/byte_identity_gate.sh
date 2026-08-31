#!/usr/bin/env bash
# Byte-identity gate for deliverable B. FOREGROUND, serial, winloci_scratch (crash rule).
# Usage: byte_identity_gate.sh freeze   # write baseline
#        byte_identity_gate.sh check    # diff current md5s vs baseline, non-zero on mismatch
set -euo pipefail
MODE="${1:?freeze|check}"
RUSTLE=/mnt/c/Users/jfris/Desktop/Rustle
# The mandated build dir is CARGO_TARGET_DIR=/mnt/linuxdisk/home/juanfraitu/rustle_target (the vhdx-backed
# in-repo target/ is not used and $RUSTLE/target/release DOES NOT EXIST); overridable for A/B-ing two
# binaries. A missing binary is a hard error here — a gate that silently produces no outputs "passes".
BIN=${BIN:-/mnt/linuxdisk/home/juanfraitu/rustle_target/release}
SCRATCH=/home/juanfra/winloci_scratch
BAM=$SCRATCH/GGO_mm.bam
FASTA=$SCRATCH/GGO.fasta
BASELINE=$RUSTLE/bench/mechanism/byte_identity_baseline.txt
WORK=$SCRATCH/bgate
# CORPUS 2 outputs are ~170 MB, which belongs on the PHYSICAL disk, not inside the sparse vhdx.
PSV_WORK=/mnt/linuxdisk/home/juanfraitu/bgate_psv
for b in copy_assign gw_family_catalog; do
  [ -x "$BIN/$b" ] || { echo "GATE ABORT: $BIN/$b missing (build with CARGO_TARGET_DIR=/mnt/linuxdisk/home/juanfraitu/rustle_target --release)"; exit 2; }
done
mkdir -p "$WORK" "$PSV_WORK"; cd "$WORK"

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

# --- CORPUS 2: the MULTI-FAMILY corpus (added 2026-08-31) -------------------------------------------
# WHY: corpus 1 (gstm) is ONE family and therefore offers ZERO cross-family candidate pairs, so any
# change scoped to interactions BETWEEN families comes back "unchanged" there BY CONSTRUCTION — the
# empty-denominator instrument failure already on record (a human 150-window panel returned 2/150
# UNCHANGED while offering 0 qualifying pairs). This corpus is the `--families` batch from
# /mnt/linuxdisk/home/juanfraitu/mec/run_psv.sh: 12 supplied families (8 admitted), 101 copies, 2 regions
# on npip3.bam, ~4 min, FOREGROUND and serial per the crash rule. It supplies 1,587 molecules genotyped
# in >= 2 families and 517 assigned in >= 2, i.e. a NON-EMPTY denominator.
MEC=/mnt/linuxdisk/home/juanfraitu/mec
if [ -f "$MEC/batch.copies.tsv" ] && [ -f /mnt/linuxdisk/home/juanfraitu/npip_cat/npip3.bam ]; then
  ( cd "$PSV_WORK" && "$BIN/copy_assign" \
      --bam /mnt/linuxdisk/home/juanfraitu/npip_cat/npip3.bam \
      --fasta /mnt/linuxdisk/home/juanfraitu/_from_wsl/winloci_scratch/GGO.fasta \
      --families "$MEC/batch.copies.tsv" \
      --copies-fa /mnt/linuxdisk/home/juanfraitu/npip_cat/arm_f2/cat.copies.fa \
      --regions "$MEC/regions.txt" --dump-psv --out psv >psv.log 2>psv.err )
else
  echo "NOTE: corpus 2 inputs absent; gate covers corpus 1 only (READ THE CANDIDATE COUNT before quoting a verdict)"
fi

# --- collect md5s of the scientific outputs ---
collect(){
  for f in gstm.assignments.tsv gstm.families.tsv gstm.quant.tsv gstm.psv_cols.tsv \
           cat.families.tsv cat.copies.tsv; do
    [ -f "$f" ] && printf '%s  %s\n' "$(md5sum < "$f" | cut -d' ' -f1)" "$f"
  done
  for f in psv.assignments.tsv psv.families.tsv psv.quant.tsv psv.psv_reads.tsv psv.psv_cols.tsv \
           psv.family_join.tsv; do
    [ -f "$PSV_WORK/$f" ] && printf '%s  %s\n' "$(md5sum < "$PSV_WORK/$f" | cut -d' ' -f1)" "$f"
  done
}

if [ "$MODE" = freeze ]; then
  # A re-freeze MUST say what made the old baseline stale, or the next reader cannot tell a legitimate
  # re-freeze from one that papered over a real regression. Pass it in REASON=...
  { echo "# baseline @ $(git -C "$RUSTLE" rev-parse --short HEAD)"
    echo "# reason: ${REASON:-UNSTATED (a re-freeze without a stated reason is not auditable)}"
    collect; } > "$BASELINE"
  echo "froze $(grep -vc '^#' "$BASELINE") md5s -> $BASELINE"
else
  cur=$(collect)
  base=$(grep -v '^#' "$BASELINE")
  if [ "$cur" = "$base" ]; then echo "GATE PASS: all md5 identical to baseline"; exit 0
  else echo "GATE FAIL: md5 drift vs baseline:"; diff <(echo "$base") <(echo "$cur") || true; exit 1; fi
fi
