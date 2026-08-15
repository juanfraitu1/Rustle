#!/usr/bin/env bash
# ─────────────────────────────────────────────────────────────────────────────────────────────────
# THE SHIPPED E_r ALIGNMENT TIER — the shell mirror of `ER_TIER_FLAGS` / `ER_SENSITIVE_SEED` in
# src/rustle/vg_family/denovo_pipeline.rs. Source this; do not re-type the flags.
#
# WHY THIS FILE EXISTS (defect B1). Every O1 property number published before 2026-08-10 was computed
# with the WRONG aligner. The binary runs
#       minimap2 -c -X --no-long-join -t N -k 11 -w 5
# but all eight bench/crossspecies panel scripts ran
#       minimap2 -x asm20|-k 11 -w 5 -c --eqx -N 200 -p 0.02 -t N        (NO -X)
# and additionally kept an `-x asm20` leg that the shipped default SKIPS
# (`RUSTLE_ER_SENSITIVE_ONLY` has defaulted to true since 2026-08-07).
#
# Measured over 14 panels on BYTE-IDENTICAL FASTA: the partition differs on 4/14 and the edge count on
# 10/14. `-N` and `-p` are INERT at this tier (they cap secondary alignments, which `-X` all-vs-all
# self-mapping does not produce here); `-X` is the ONE operative difference, because it implies
# `--dual=no` and therefore emits ONE orientation per pair.
#
# ⚠ `--eqx` is also dropped: the binary does not pass it. It expands the CIGAR M ops into =/X and does
# not touch de:f:, coordinates, or block length, so no downstream statistic in these scripts moves --
# but a panel that claims to run "the shipped command" must run the shipped command.
#
# ⚠ THIS IS THE ALL-VS-ALL (E_r) TIER ONLY. Genome PASSES (seed -> genome, locating copies) are a
# different question and legitimately keep `-N 200 -p 0.02` and drop `-X`: there the query and target
# are different files, secondaries ARE the answer, and `-X` would be wrong.

# Usage: er_tier_allvsall <threads> <fasta> <out.paf>
er_tier_allvsall() {
  local T="$1" FA="$2" OUT="$3"
  minimap2 -c -X --no-long-join -t "$T" -k 11 -w 5 "$FA" "$FA" > "$OUT" 2>/dev/null
}

# The literal argv, for logging into a run's provenance file.
ER_TIER_CMDLINE="minimap2 -c -X --no-long-join -t <T> -k 11 -w 5 <fa> <fa>"
