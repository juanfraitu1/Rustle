#!/usr/bin/env bash
# Standing Layer-2 invariant harness. Per-chrom serial (whole-genome -L OOMs).
# Run: bench/layer2_invariant.sh
#
# Runs on BOTH chr19 (GGO_19.bam) AND a chrY family chrom
# (bench/fixtures/chrY_family.bam). Two tiers of invariant:
#
#   (A) HARD (must pass every milestone until consumers emit, M3.2+):
#       --vg --vg-layer2 output == --vg output. Layer 2 is ADDITIVE — until the
#       union-by-chain emitter lands (M6) it must not change a single chain.
#
#   (B) TARGET (the project goal): baseline coord-signatures ⊆ VG-layer2 output.
#       This is NOT satisfied yet — plain --vg drops some baseline chains
#       (pre-existing VG-vs-baseline regression; see project_vg_baseline_regression).
#       M6's union-by-chain closes it. Reported here as a measured gap; promote to a
#       HARD failure with LAYER2_STRICT=1 (flip on in CI once M6 lands).
#
# Optional leg (C): RUSTLE_PRECISE=1 byte-identical to commit 4705ab1 — runs only if
#   bench/ref/4705ab1_precise_GGO_19.gtf exists. Generate it on a CLEAN tree with
#   LAYER2_GEN_REF=1 (does a git checkout 4705ab1 build dance; refuses on dirty tree).
set -euo pipefail
cd "$(dirname "$0")/.."

BIN=./target/release/rustle
REF_PRECISE=bench/ref/4705ab1_precise_GGO_19.gtf
CHRY_BAM=${CHRY_BAM:-bench/fixtures/chrY_family.bam}
SUPERSET=scripts/coord_signature_superset.py
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

# --- per-chromosome invariant runner ---
check_chrom () {
  local tag="$1" bam="$2"
  echo "== [$tag] generate baseline / VG-default / VG-layer2 =="
  "$BIN" -L "$bam" -o "/tmp/layer2/${tag}_baseline.gtf" 2>/dev/null
  "$BIN" -L "$bam" --vg -o "/tmp/layer2/${tag}_vg_default.gtf" 2>/dev/null
  "$BIN" -L "$bam" --vg --vg-layer2 -o "/tmp/layer2/${tag}_vg_layer2.gtf" 2>/dev/null

  # (A) HARD — additivity vs the VG copy-recovery output: Layer 2 may ADD chains
  #     (the recovered baseline floor + future novel copies) but must NEVER DROP a
  #     chain that --vg already produced. So vg_layer2 ⊇ vg_default.
  if python3 "$SUPERSET" "/tmp/layer2/${tag}_vg_layer2.gtf" "/tmp/layer2/${tag}_vg_default.gtf" >/dev/null; then
    echo "  OK [$tag] (A): VG Layer-2 ⊇ VG default (additive; no VG chain dropped)"
  else
    echo "  FAIL [$tag] (A): --vg-layer2 DROPPED a chain present in --vg — Layer 2 must be additive"
    exit 1
  fi

  # (B) HARD (M-FLOOR established 2026-06-13) — baseline ⊆ VG-layer2. The pulled-forward
  #     union floor guarantees every baseline chain survives VG. Any later milestone that
  #     breaks this is a regression. (LAYER2_STRICT kept as a no-op alias for callers.)
  if python3 "$SUPERSET" "/tmp/layer2/${tag}_vg_layer2.gtf" "/tmp/layer2/${tag}_baseline.gtf"; then
    echo "  OK [$tag] (B): VG Layer-2 ⊇ baseline (floor holds)"
  else
    echo "  FAIL [$tag] (B): baseline NOT subset of VG-layer2 — M-FLOOR regression"
    exit 1
  fi
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

echo "== chr19 invariants =="
check_chrom chr19 GGO_19.bam

echo "== chrY invariants (repeat-rich; required) =="
if [ ! -f "$CHRY_BAM" ]; then
  echo "  FAIL: chrY family BAM ($CHRY_BAM) missing — the chrY invariant is REQUIRED."
  exit 1
fi
check_chrom chrY "$CHRY_BAM"

echo "ALL INVARIANTS PASS: VG Layer-2 is additive (⊇ VG default) AND ⊇ baseline (M-FLOOR floor holds) on chr19 + chrY."
