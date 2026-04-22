#!/usr/bin/env bash
# Run from repo root. Compare panel output to the locus reference slice with -R -Q.
set -euo pipefail

# STRG.251 refs=STRG.251.1,STRG.251.2,STRG.251.4,STRG.251.5,STRG.251.8,STRG.251.9
RUSTLE_TRACE_LOCUS=36443168-36634779 \
  target/release/rustle GGO_19.bam -L --chrom NC_073243.2 \
  -o bench/ggo19_needy_top5_run_STRG_251.gtf \
  --debug-bundle "NC_073243.2:36443168-36634779" --only-debug-bundle \
  --trace-reference bench/ggo19_needy_top5_refs/STRG_251.gtf \
  --trace-output bench/ggo19_needy_top5_run_STRG_251.trace.txt \
  --parity-stage-tsv bench/ggo19_needy_top5_run_STRG_251.parity.tsv \
  --snapshot-jsonl bench/ggo19_needy_top5_run_STRG_251.snap.jsonl

gffcompare -r bench/ggo19_needy_top5_refs/STRG_251.gtf -R -Q -o bench/ggo19_needy_top5_run_STRG_251_cmp bench/ggo19_needy_top5_run_STRG_251.gtf


# STRG.151 refs=STRG.151.1,STRG.151.2,STRG.151.5,STRG.151.6,STRG.151.7,STRG.151.4
RUSTLE_TRACE_LOCUS=23033526-23064635 \
  target/release/rustle GGO_19.bam -L --chrom NC_073243.2 \
  -o bench/ggo19_needy_top5_run_STRG_151.gtf \
  --debug-bundle "NC_073243.2:23033526-23064635" --only-debug-bundle \
  --trace-reference bench/ggo19_needy_top5_refs/STRG_151.gtf \
  --trace-output bench/ggo19_needy_top5_run_STRG_151.trace.txt \
  --parity-stage-tsv bench/ggo19_needy_top5_run_STRG_151.parity.tsv \
  --snapshot-jsonl bench/ggo19_needy_top5_run_STRG_151.snap.jsonl

gffcompare -r bench/ggo19_needy_top5_refs/STRG_151.gtf -R -Q -o bench/ggo19_needy_top5_run_STRG_151_cmp bench/ggo19_needy_top5_run_STRG_151.gtf


# STRG.503 refs=STRG.503.1,STRG.503.4,STRG.503.5,STRG.503.6,STRG.503.2
RUSTLE_TRACE_LOCUS=97540381-97743862 \
  target/release/rustle GGO_19.bam -L --chrom NC_073243.2 \
  -o bench/ggo19_needy_top5_run_STRG_503.gtf \
  --debug-bundle "NC_073243.2:97540381-97743862" --only-debug-bundle \
  --trace-reference bench/ggo19_needy_top5_refs/STRG_503.gtf \
  --trace-output bench/ggo19_needy_top5_run_STRG_503.trace.txt \
  --parity-stage-tsv bench/ggo19_needy_top5_run_STRG_503.parity.tsv \
  --snapshot-jsonl bench/ggo19_needy_top5_run_STRG_503.snap.jsonl

gffcompare -r bench/ggo19_needy_top5_refs/STRG_503.gtf -R -Q -o bench/ggo19_needy_top5_run_STRG_503_cmp bench/ggo19_needy_top5_run_STRG_503.gtf


# STRG.157 refs=STRG.157.6,STRG.157.7,STRG.157.3,STRG.157.5,STRG.157.4
RUSTLE_TRACE_LOCUS=23288512-23308475 \
  target/release/rustle GGO_19.bam -L --chrom NC_073243.2 \
  -o bench/ggo19_needy_top5_run_STRG_157.gtf \
  --debug-bundle "NC_073243.2:23288512-23308475" --only-debug-bundle \
  --trace-reference bench/ggo19_needy_top5_refs/STRG_157.gtf \
  --trace-output bench/ggo19_needy_top5_run_STRG_157.trace.txt \
  --parity-stage-tsv bench/ggo19_needy_top5_run_STRG_157.parity.tsv \
  --snapshot-jsonl bench/ggo19_needy_top5_run_STRG_157.snap.jsonl

gffcompare -r bench/ggo19_needy_top5_refs/STRG_157.gtf -R -Q -o bench/ggo19_needy_top5_run_STRG_157_cmp bench/ggo19_needy_top5_run_STRG_157.gtf


# STRG.453 refs=STRG.453.2,STRG.453.3,STRG.453.4,STRG.453.5
RUSTLE_TRACE_LOCUS=82243221-82437729 \
  target/release/rustle GGO_19.bam -L --chrom NC_073243.2 \
  -o bench/ggo19_needy_top5_run_STRG_453.gtf \
  --debug-bundle "NC_073243.2:82243221-82437729" --only-debug-bundle \
  --trace-reference bench/ggo19_needy_top5_refs/STRG_453.gtf \
  --trace-output bench/ggo19_needy_top5_run_STRG_453.trace.txt \
  --parity-stage-tsv bench/ggo19_needy_top5_run_STRG_453.parity.tsv \
  --snapshot-jsonl bench/ggo19_needy_top5_run_STRG_453.snap.jsonl

gffcompare -r bench/ggo19_needy_top5_refs/STRG_453.gtf -R -Q -o bench/ggo19_needy_top5_run_STRG_453_cmp bench/ggo19_needy_top5_run_STRG_453.gtf

