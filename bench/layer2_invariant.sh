#!/usr/bin/env bash
# Standing Layer-2 invariant harness. Per-chrom serial (whole-genome -L OOMs).
# Run: bench/layer2_invariant.sh
#
# Always-run current-tree invariants, on BOTH chr19 (GGO_19.bam) AND a chrY family
# chrom (bench/fixtures/chrY_family.bam):
#   (A) VG default (--vg, Layer 2 OFF) == baseline (mutual coord-signature superset)
#   (B) VG Layer 2 (--vg --vg-layer2) superset baseline
# Optional leg (C): RUSTLE_PRECISE=1 byte-identical to commit 4705ab1 — runs only if
#   bench/ref/4705ab1_precise_GGO_19.gtf exists. Generate it on a CLEAN tree with
#   LAYER2_GEN_REF=1 (does a git checkout 4705ab1 build dance; refuses on dirty tree).
set -euo pipefail
cd "$(dirname "$0")/.."

BIN=./target/release/rustle
REF_PRECISE=bench/ref/4705ab1_precise_GGO_19.gtf
CHRY_BAM=${CHRY_BAM:-bench/fixtures/chrY_family.bam}
mkdir -p bench/ref /tmp/layer2
export RAYON_NUM_THREADS=1

echo "== build (current tree) =="
cargo build --release

# --- optional: generate the 4705ab1 reference (CLEAN tree only, opt-in) ---
if [ "${LAYER2_GEN_REF:-0}" = "1" ] && [ ! -f "$REF_PRECISE" ]; then
  if [ -n "$(git status --porcelain)" ]; then
    echo "  REFUSING to generate 4705ab1 reference: working tree is dirty."
    echo "  Commit/stash your changes first, then re-run with LAYER2_GEN_REF=1."
    exit 1
  fi
  echo "  generating 4705ab1 reference"
  cur=$(git rev-parse --abbrev-ref HEAD)
  git checkout 4705ab1
  cargo build --release
  RUSTLE_PRECISE=1 "$BIN" -L GGO_19.bam -o "$REF_PRECISE" 2>/dev/null
  git checkout "$cur"
  cargo build --release
fi

# --- per-chromosome current-tree invariant runner ---
check_chrom () {
  local tag="$1" bam="$2"
  echo "== [$tag] baseline (no --vg) =="
  "$BIN" -L "$bam" -o "/tmp/layer2/${tag}_baseline.gtf" 2>/dev/null
  echo "== [$tag] VG default (Layer 2 OFF) == baseline =="
  "$BIN" -L "$bam" --vg -o "/tmp/layer2/${tag}_vg_default.gtf" 2>/dev/null
  python3 scripts/coord_signature_superset.py "/tmp/layer2/${tag}_vg_default.gtf" "/tmp/layer2/${tag}_baseline.gtf"
  python3 scripts/coord_signature_superset.py "/tmp/layer2/${tag}_baseline.gtf" "/tmp/layer2/${tag}_vg_default.gtf" \
    && echo "  OK [$tag]: VG-default == baseline (mutual superset)"
  echo "== [$tag] VG Layer-2 superset baseline =="
  "$BIN" -L "$bam" --vg --vg-layer2 -o "/tmp/layer2/${tag}_vg_layer2.gtf" 2>/dev/null
  python3 scripts/coord_signature_superset.py "/tmp/layer2/${tag}_vg_layer2.gtf" "/tmp/layer2/${tag}_baseline.gtf" \
    && echo "  OK [$tag]: VG superset baseline with Layer 2 on"
}

# --- (C) optional precise byte-identity leg ---
if [ -f "$REF_PRECISE" ]; then
  echo "== (C) RUSTLE_PRECISE byte-identity vs 4705ab1 (chr19) =="
  RUSTLE_PRECISE=1 "$BIN" -L GGO_19.bam -o /tmp/layer2/precise.gtf 2>/dev/null
  diff -q "$REF_PRECISE" /tmp/layer2/precise.gtf \
    && echo "  OK: byte-identical to 4705ab1" \
    || { echo "  FAIL: RUSTLE_PRECISE drifted from 4705ab1"; exit 1; }
else
  echo "== (C) RUSTLE_PRECISE byte-identity vs 4705ab1: SKIPPED (no reference) =="
  echo "   generate on a clean tree: LAYER2_GEN_REF=1 bench/layer2_invariant.sh"
fi

echo "== (A,B) chr19 current-tree invariants =="
check_chrom chr19 GGO_19.bam

echo "== (A,B) chrY current-tree invariants (repeat-rich; required) =="
if [ ! -f "$CHRY_BAM" ]; then
  echo "  FAIL: chrY family BAM ($CHRY_BAM) missing — the chrY superset invariant is REQUIRED."
  exit 1
fi
check_chrom chrY "$CHRY_BAM"

echo "ALL CURRENT-TREE INVARIANTS PASS"
