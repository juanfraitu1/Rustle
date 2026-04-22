#!/usr/bin/env bash
# Ranks 6–10 from GGO_19_after_termfix_missing (see bench/ggo19_needy_top15.panel.tsv).
# Per-bundle summary: auto-written as bench/ggo19_needy_next5_run_STRG_*.needy_bundle_summary.tsv
# (omit --parity-stage-tsv so rustle uses needy-locus auto path).
# Run from repo root. Requires: target/release/rustle, GGO_19.bam, gffcompare.
set -euo pipefail
REFDIR=bench/ggo19_needy_top15_refs
OUT=bench/ggo19_needy_next5_run

# STRG.442
RUSTLE_TRACE_LOCUS=80297718-80330444 \
  target/release/rustle GGO_19.bam -L --chrom NC_073243.2 \
  -o "${OUT}_STRG_442.gtf" \
  --debug-bundle "NC_073243.2:80297718-80330444" --only-debug-bundle \
  --trace-reference "${REFDIR}/STRG_442.gtf" \
  --trace-output "${OUT}_STRG_442.trace.txt" \
  2> "${OUT}_STRG_442.stderr.log"
gffcompare -r "${REFDIR}/STRG_442.gtf" -R -Q -o "${OUT}_STRG_442_cmp" "${OUT}_STRG_442.gtf"

# STRG.566
RUSTLE_TRACE_LOCUS=111146297-111177575 \
  target/release/rustle GGO_19.bam -L --chrom NC_073243.2 \
  -o "${OUT}_STRG_566.gtf" \
  --debug-bundle "NC_073243.2:111146297-111177575" --only-debug-bundle \
  --trace-reference "${REFDIR}/STRG_566.gtf" \
  --trace-output "${OUT}_STRG_566.trace.txt" \
  2> "${OUT}_STRG_566.stderr.log"
gffcompare -r "${REFDIR}/STRG_566.gtf" -R -Q -o "${OUT}_STRG_566_cmp" "${OUT}_STRG_566.gtf"

# STRG.445
RUSTLE_TRACE_LOCUS=80975942-81174569 \
  target/release/rustle GGO_19.bam -L --chrom NC_073243.2 \
  -o "${OUT}_STRG_445.gtf" \
  --debug-bundle "NC_073243.2:80975942-81174569" --only-debug-bundle \
  --trace-reference "${REFDIR}/STRG_445.gtf" \
  --trace-output "${OUT}_STRG_445.trace.txt" \
  2> "${OUT}_STRG_445.stderr.log"
gffcompare -r "${REFDIR}/STRG_445.gtf" -R -Q -o "${OUT}_STRG_445_cmp" "${OUT}_STRG_445.gtf"

# STRG.29
RUSTLE_TRACE_LOCUS=17431800-17453950 \
  target/release/rustle GGO_19.bam -L --chrom NC_073243.2 \
  -o "${OUT}_STRG_29.gtf" \
  --debug-bundle "NC_073243.2:17431800-17453950" --only-debug-bundle \
  --trace-reference "${REFDIR}/STRG_29.gtf" \
  --trace-output "${OUT}_STRG_29.trace.txt" \
  2> "${OUT}_STRG_29.stderr.log"
gffcompare -r "${REFDIR}/STRG_29.gtf" -R -Q -o "${OUT}_STRG_29_cmp" "${OUT}_STRG_29.gtf"

# STRG.440
RUSTLE_TRACE_LOCUS=80273347-80288948 \
  target/release/rustle GGO_19.bam -L --chrom NC_073243.2 \
  -o "${OUT}_STRG_440.gtf" \
  --debug-bundle "NC_073243.2:80273347-80288948" --only-debug-bundle \
  --trace-reference "${REFDIR}/STRG_440.gtf" \
  --trace-output "${OUT}_STRG_440.trace.txt" \
  2> "${OUT}_STRG_440.stderr.log"
gffcompare -r "${REFDIR}/STRG_440.gtf" -R -Q -o "${OUT}_STRG_440_cmp" "${OUT}_STRG_440.gtf"
