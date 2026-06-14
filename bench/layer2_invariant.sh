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

# (5) simulated 2-copy paralog. Originally this asserted strict layer2==default
# (Layer 2 must add NOTHING) to guard the chimeric cross-copy emission bug
# (decompose emitting a node's UNION span instead of per-copy coordinates). That
# strict form is now RELAXED (2026-06-13): Part A (secondary-enriched per-copy
# isoform discovery, default-on under --vg-layer2) legitimately ADDS within-copy
# alt-splice isoforms — real-data-validated on chrY as read-backed novel isoforms
# of real paralog genes (10 j / 1 c / 3 n, ZERO invented junctions, ZERO
# phantoms). On this fixture Part A adds copy B's exon-skip isoform (chrSIM
# 30001-32000) from copyA-skip secondaries shadowing copyB's near-identical locus.
# So we no longer forbid additions; we keep the two invariants that actually
# matter:
#   (i)  ADDITIVE — layer2 ⊇ default (drops nothing VG already produced).
#   (ii) NO CHIMERA — no transcript spans BOTH copies. The two copies sit at
#        chrSIM:10001-12000 and 30001-32000 (each ~2kb, ~18kb apart); the original
#        union-span chimera bridged them (10001-32000, ~22kb). Every legitimate
#        per-copy isoform is bounded by one copy's ~2kb extent, so we cap the max
#        transcript span well below the inter-copy distance (5kb: above any single
#        copy, far below a ~22kb cross-copy union). A reborn chimera blows past it.
echo "== (5) simulated paralog: Layer-2 additive + NO cross-copy chimera =="
SIM_BAM=${SIM_BAM:-bench/fixtures/sim_paralog.bam}
SIM_FA=${SIM_FA:-bench/fixtures/sim_paralog.fa}
SIM_CHIMERA_SPAN_MAX=${SIM_CHIMERA_SPAN_MAX:-5000}
if [ -f "$SIM_BAM" ] && [ -f "$SIM_FA" ]; then
  "$BIN" -L "$SIM_BAM" --vg --genome-fasta "$SIM_FA" -o /tmp/layer2/sim_default.gtf 2>/dev/null
  "$BIN" -L "$SIM_BAM" --vg --vg-layer2 --genome-fasta "$SIM_FA" -o /tmp/layer2/sim_layer2.gtf 2>/dev/null
  d=$(grep -c $'\ttranscript\t' /tmp/layer2/sim_default.gtf || true)
  l=$(grep -c $'\ttranscript\t' /tmp/layer2/sim_layer2.gtf || true)
  echo "  sim transcripts: default=$d layer2=$l"
  # (i) additive: layer2 ⊇ default (Part A may add; it must never drop).
  if ! python3 "$SUPERSET" /tmp/layer2/sim_layer2.gtf /tmp/layer2/sim_default.gtf >/dev/null; then
    echo "  FAIL: layer2 DROPPED a sim chain present in default (not additive)"
    exit 1
  fi
  # (i+) positive Part-A regression guard: Part A MUST recover copy B's exon-skip
  # isoform from the copyA-skip secondaries shadowing copyB's locus, so layer2 has
  # strictly more transcripts than default. If this stops firing, Part A regressed.
  if [ "$l" -le "$d" ]; then
    echo "  FAIL: Part A added nothing on sim ($l <= $d) — multi-isoform recovery regressed"
    exit 1
  fi
  # (ii) no cross-copy chimera: max transcript span must stay within a single copy.
  lmax=$(awk -F'\t' '$3=="transcript"{s=$5-$4; if(s>m)m=s} END{print m+0}' /tmp/layer2/sim_layer2.gtf)
  echo "  sim max transcript span: layer2=$lmax (chimera bound=$SIM_CHIMERA_SPAN_MAX)"
  if [ "$lmax" -le "$SIM_CHIMERA_SPAN_MAX" ]; then
    echo "  OK: layer2 additive + no cross-copy chimera (max span within single-copy bound)"
  else
    echo "  FAIL: layer2 transcript span $lmax > $SIM_CHIMERA_SPAN_MAX — possible cross-copy union chimera"
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

# (8) GHOST fixture: M5 recovers a no-primary copy. Distinct from leg 6's starved
# copy (which has 1 primary → a bundle forms → the existing homology-transfer path
# handles it). Here copy B has ZERO primary reads: Layer 1 builds NO bundle there,
# so its secondaries get locus=None and are dropped — ONLY M5 (the default-off
# --vg-layer2-new-copies path) clusters those locus=None secondaries and recovers
# the copy. So copy B is ABSENT with --vg-layer2 (M5 OFF) and RECOVERED only when
# --vg-layer2-new-copies is added. Must also stay additive (on ⊇ off).
echo "== (8) ghost fixture: M5 recovers a no-primary copy =="
GH_BAM=${GH_BAM:-bench/fixtures/sim_ghost.bam}; GH_FA=${GH_FA:-bench/fixtures/sim_ghost.fa}
if [ -f "$GH_BAM" ] && [ -f "$GH_FA" ]; then
  "$BIN" -L "$GH_BAM" --vg --vg-layer2 --genome-fasta "$GH_FA" -o /tmp/layer2/ghost_off.gtf 2>/dev/null
  "$BIN" -L "$GH_BAM" --vg --vg-layer2 --vg-layer2-new-copies --genome-fasta "$GH_FA" -o /tmp/layer2/ghost_on.gtf 2>/dev/null
  goff=$(awk -F'\t' '$3=="transcript" && $4>20000' /tmp/layer2/ghost_off.gtf | wc -l)
  gon=$(awk -F'\t' '$3=="transcript" && $4>20000' /tmp/layer2/ghost_on.gtf | wc -l)
  # copy A sits at <20000, copy B at >20000, so every exon line >20000 belongs to a
  # recovered copy-B transcript. More copy-B exon lines than copy-B transcripts ⇒ at
  # least one recovered transcript is MULTI-EXON (pigeonhole) — M5 must reconstruct a
  # spliced chain, not a single-exon blob (its gate requires >=2 exons).
  gex=$(awk -F'\t' '$3=="exon" && $4>20000' /tmp/layer2/ghost_on.gtf | wc -l)
  echo "  copyB transcripts (region>20000): M5-off=$goff M5-on=$gon ($gex copyB exon lines)"
  # (i) ghost gate: M5 OFF must leave copy B absent (no bundle ever formed); M5 ON
  #     must recover it MULTI-EXON. If OFF != 0 the existing path already handles it
  #     (wrong fixture); if ON < 1 M5 regressed; if exons <= tx the recovery degenerated
  #     to single-exon blobs.
  if [ "$goff" -eq 0 ] && [ "$gon" -ge 1 ] && [ "$gex" -gt "$gon" ]; then
    echo "  OK: M5 recovered the no-primary ghost copy B, multi-exon (off=0 → on=$gon tx / $gex exons)"
  else
    echo "  FAIL: ghost copy B recovery wrong (M5-off=$goff must be 0, M5-on=$gon must be >=1, exons=$gex must be > tx)"
    exit 1
  fi
  # (ii) additive: M5-on ⊇ M5-off (recovery must never drop a chain VG already had).
  if python3 "$SUPERSET" /tmp/layer2/ghost_on.gtf /tmp/layer2/ghost_off.gtf >/dev/null; then
    echo "  OK: M5-on ⊇ M5-off (additive; no chain dropped)"
  else
    echo "  FAIL: --vg-layer2-new-copies DROPPED a chain present without it (not additive)"
    exit 1
  fi
else echo "  SKIP: ghost fixture absent"; fi

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
