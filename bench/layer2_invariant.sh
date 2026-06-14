#!/usr/bin/env bash
# Standing Layer-2 invariant harness. Per-chrom serial (whole-genome -L OOMs).
# Run: bench/layer2_invariant.sh
#
# Runs on BOTH chr19 (GGO_19.bam) AND a chrY family chrom
# (bench/fixtures/chrY_family.bam). Two tiers of invariant:
#
#   (A) HARD — additivity: --vg --vg-layer2 ⊇ --vg output. Layer 2 may ADD chains
#       (the recovered baseline floor + future novel copies) but must NEVER DROP a
#       chain --vg already produced.
#
#   (B) HARD — the floor (established by M-FLOOR, 2026-06-13): baseline coord-signatures
#       ⊆ --vg --vg-layer2 output. The pulled-forward union floor recovers every baseline
#       chain VG's family-EM dropped (chr19 RSTL.589.1 + 18 chrY chains). Any later
#       milestone that breaks this is a regression. LAYER2_STRICT kept as a no-op alias.
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

# (5) simulated 2-copy paralog NO-FALSE-POSITIVE check. Both copies already
# assemble in --vg (copyB is NOT starved here — its 5 primary reads suffice), so
# Layer 2 must add NOTHING: layer2 == default. This guards against the chimeric
# cross-copy emission bug (decompose emitting a node's union span instead of
# per-copy coordinates) that this fixture originally exposed. (A genuine
# starved-copy RECOVERY fixture — needing homology coordinate-transfer — is
# separate future work.)
echo "== (5) simulated paralog: Layer-2 emits NO chimera (layer2 == default) =="
SIM_BAM=${SIM_BAM:-bench/fixtures/sim_paralog.bam}
SIM_FA=${SIM_FA:-bench/fixtures/sim_paralog.fa}
if [ -f "$SIM_BAM" ] && [ -f "$SIM_FA" ]; then
  "$BIN" -L "$SIM_BAM" --vg --genome-fasta "$SIM_FA" -o /tmp/layer2/sim_default.gtf 2>/dev/null
  "$BIN" -L "$SIM_BAM" --vg --vg-layer2 --genome-fasta "$SIM_FA" -o /tmp/layer2/sim_layer2.gtf 2>/dev/null
  d=$(grep -c $'\ttranscript\t' /tmp/layer2/sim_default.gtf || true)
  l=$(grep -c $'\ttranscript\t' /tmp/layer2/sim_layer2.gtf || true)
  echo "  sim transcripts: default=$d layer2=$l"
  if python3 "$SUPERSET" /tmp/layer2/sim_layer2.gtf /tmp/layer2/sim_default.gtf >/dev/null \
     && python3 "$SUPERSET" /tmp/layer2/sim_default.gtf /tmp/layer2/sim_layer2.gtf >/dev/null; then
    echo "  OK: layer2 == default on sim (no chimera, no spurious addition)"
  else
    echo "  FAIL: layer2 != default on sim — Layer 2 emitted a spurious/chimeric chain"
    exit 1
  fi
else
  echo "  SKIP: $SIM_BAM / $SIM_FA not present"
fi

echo "== (6) starvation fixture: Layer-2 RECOVERS the dropped copy =="
ST_BAM=${ST_BAM:-bench/fixtures/sim_starved.bam}; ST_FA=${ST_FA:-bench/fixtures/sim_starved.fa}
if [ -f "$ST_BAM" ] && [ -f "$ST_FA" ]; then
  "$BIN" -L "$ST_BAM" --vg --genome-fasta "$ST_FA" -o /tmp/layer2/starved_default.gtf 2>/dev/null
  "$BIN" -L "$ST_BAM" --vg --vg-layer2 --genome-fasta "$ST_FA" -o /tmp/layer2/starved_layer2.gtf 2>/dev/null
  db=$(awk -F'\t' '$3=="transcript" && $4>20000' /tmp/layer2/starved_default.gtf | wc -l)
  lb=$(awk -F'\t' '$3=="transcript" && $4>20000' /tmp/layer2/starved_layer2.gtf | wc -l)
  echo "  copyB transcripts (region>20000): default=$db layer2=$lb"
  if [ "$db" -eq 0 ] && [ "$lb" -ge 1 ]; then echo "  OK: Layer-2 recovered the starved copy B"; else echo "  FAIL: starved copy B not recovered (default=$db layer2=$lb)"; exit 1; fi
else echo "  SKIP: starved fixture absent"; fi

# (7) GENOME-BACKED real chrY recovery leg. The legs above run chrY WITHOUT a
# genome, so run_layer2 (which requires --genome-fasta) is SKIPPED — they only
# test the M-FLOOR floor, not homology-transfer recovery. This leg exercises the
# real genome-backed homology-transfer path on real chrY families (the path that
# panicked on real data before the loci re-sync fix). Env-gated on a full genome
# FASTA covering chrY; SKIPs if absent (large, not in-repo). HARD: must not panic,
# must stay additive (layer2 ⊇ default). Real-data validation confirmed 2 real
# copy recoveries here (RSTL.82.1=LOC129530264) with zero false positives.
echo "== (7) genome-backed chrY: homology-transfer runs + stays additive =="
GGO_GENOME=${GGO_GENOME:-/mnt/c/Users/jfris/Desktop/GGO.fasta}
if [ -f "$GGO_GENOME" ]; then
  "$BIN" -L "$CHRY_BAM" --vg --genome-fasta "$GGO_GENOME" -o /tmp/layer2/chrYg_default.gtf 2>/dev/null
  "$BIN" -L "$CHRY_BAM" --vg --vg-layer2 --genome-fasta "$GGO_GENOME" -o /tmp/layer2/chrYg_layer2.gtf 2>/dev/null \
    || { echo "  FAIL: genome-backed --vg-layer2 errored/panicked on chrY"; exit 1; }
  if python3 "$SUPERSET" /tmp/layer2/chrYg_layer2.gtf /tmp/layer2/chrYg_default.gtf >/dev/null; then
    nv=$(grep -cP 'copy_status "novel"' /tmp/layer2/chrYg_layer2.gtf || true)
    echo "  OK: genome-backed layer2 ⊇ default (additive); $nv homology-recovered copy(ies)"
  else
    echo "  FAIL: genome-backed layer2 DROPPED a chain present in --vg default"; exit 1
  fi
else
  echo "  SKIP: genome $GGO_GENOME absent (set GGO_GENOME to enable)"
fi

echo "ALL INVARIANTS PASS: VG Layer-2 is additive (⊇ VG default) AND ⊇ baseline (M-FLOOR floor holds) on chr19 + chrY."
