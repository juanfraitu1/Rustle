#!/usr/bin/env bash
# Run from repo root. Compare panel output to the locus reference slice with -R -Q.
set -euo pipefail

# STRG.251 refs=STRG.251.1,STRG.251.2,STRG.251.4,STRG.251.5,STRG.251.8,STRG.251.9
RUSTLE_TRACE_LOCUS=36443168-36634779 \
  target/release/rustle GGO_19.bam -L --chrom NC_073243.2 \
  -o bench/ggo19_needy_top15_run_STRG_251.gtf \
  --debug-bundle "NC_073243.2:36443168-36634779" --only-debug-bundle \
  --trace-reference bench/ggo19_needy_top15_refs/STRG_251.gtf \
  --trace-output bench/ggo19_needy_top15_run_STRG_251.trace.txt \
  --parity-stage-tsv bench/ggo19_needy_top15_run_STRG_251.parity.tsv \
  --snapshot-jsonl bench/ggo19_needy_top15_run_STRG_251.snap.jsonl

gffcompare -r bench/ggo19_needy_top15_refs/STRG_251.gtf -R -Q -o bench/ggo19_needy_top15_run_STRG_251_cmp bench/ggo19_needy_top15_run_STRG_251.gtf


# STRG.151 refs=STRG.151.1,STRG.151.2,STRG.151.5,STRG.151.6,STRG.151.7,STRG.151.4
RUSTLE_TRACE_LOCUS=23033526-23064635 \
  target/release/rustle GGO_19.bam -L --chrom NC_073243.2 \
  -o bench/ggo19_needy_top15_run_STRG_151.gtf \
  --debug-bundle "NC_073243.2:23033526-23064635" --only-debug-bundle \
  --trace-reference bench/ggo19_needy_top15_refs/STRG_151.gtf \
  --trace-output bench/ggo19_needy_top15_run_STRG_151.trace.txt \
  --parity-stage-tsv bench/ggo19_needy_top15_run_STRG_151.parity.tsv \
  --snapshot-jsonl bench/ggo19_needy_top15_run_STRG_151.snap.jsonl

gffcompare -r bench/ggo19_needy_top15_refs/STRG_151.gtf -R -Q -o bench/ggo19_needy_top15_run_STRG_151_cmp bench/ggo19_needy_top15_run_STRG_151.gtf


# STRG.503 refs=STRG.503.1,STRG.503.4,STRG.503.5,STRG.503.6,STRG.503.2
RUSTLE_TRACE_LOCUS=97540381-97743862 \
  target/release/rustle GGO_19.bam -L --chrom NC_073243.2 \
  -o bench/ggo19_needy_top15_run_STRG_503.gtf \
  --debug-bundle "NC_073243.2:97540381-97743862" --only-debug-bundle \
  --trace-reference bench/ggo19_needy_top15_refs/STRG_503.gtf \
  --trace-output bench/ggo19_needy_top15_run_STRG_503.trace.txt \
  --parity-stage-tsv bench/ggo19_needy_top15_run_STRG_503.parity.tsv \
  --snapshot-jsonl bench/ggo19_needy_top15_run_STRG_503.snap.jsonl

gffcompare -r bench/ggo19_needy_top15_refs/STRG_503.gtf -R -Q -o bench/ggo19_needy_top15_run_STRG_503_cmp bench/ggo19_needy_top15_run_STRG_503.gtf


# STRG.157 refs=STRG.157.6,STRG.157.7,STRG.157.3,STRG.157.5,STRG.157.4
RUSTLE_TRACE_LOCUS=23288512-23308475 \
  target/release/rustle GGO_19.bam -L --chrom NC_073243.2 \
  -o bench/ggo19_needy_top15_run_STRG_157.gtf \
  --debug-bundle "NC_073243.2:23288512-23308475" --only-debug-bundle \
  --trace-reference bench/ggo19_needy_top15_refs/STRG_157.gtf \
  --trace-output bench/ggo19_needy_top15_run_STRG_157.trace.txt \
  --parity-stage-tsv bench/ggo19_needy_top15_run_STRG_157.parity.tsv \
  --snapshot-jsonl bench/ggo19_needy_top15_run_STRG_157.snap.jsonl

gffcompare -r bench/ggo19_needy_top15_refs/STRG_157.gtf -R -Q -o bench/ggo19_needy_top15_run_STRG_157_cmp bench/ggo19_needy_top15_run_STRG_157.gtf


# STRG.453 refs=STRG.453.2,STRG.453.3,STRG.453.4,STRG.453.5
RUSTLE_TRACE_LOCUS=82243221-82437729 \
  target/release/rustle GGO_19.bam -L --chrom NC_073243.2 \
  -o bench/ggo19_needy_top15_run_STRG_453.gtf \
  --debug-bundle "NC_073243.2:82243221-82437729" --only-debug-bundle \
  --trace-reference bench/ggo19_needy_top15_refs/STRG_453.gtf \
  --trace-output bench/ggo19_needy_top15_run_STRG_453.trace.txt \
  --parity-stage-tsv bench/ggo19_needy_top15_run_STRG_453.parity.tsv \
  --snapshot-jsonl bench/ggo19_needy_top15_run_STRG_453.snap.jsonl

gffcompare -r bench/ggo19_needy_top15_refs/STRG_453.gtf -R -Q -o bench/ggo19_needy_top15_run_STRG_453_cmp bench/ggo19_needy_top15_run_STRG_453.gtf


# STRG.442 refs=STRG.442.3,STRG.442.6,STRG.442.9
RUSTLE_TRACE_LOCUS=80297718-80330444 \
  target/release/rustle GGO_19.bam -L --chrom NC_073243.2 \
  -o bench/ggo19_needy_top15_run_STRG_442.gtf \
  --debug-bundle "NC_073243.2:80297718-80330444" --only-debug-bundle \
  --trace-reference bench/ggo19_needy_top15_refs/STRG_442.gtf \
  --trace-output bench/ggo19_needy_top15_run_STRG_442.trace.txt \
  --parity-stage-tsv bench/ggo19_needy_top15_run_STRG_442.parity.tsv \
  --snapshot-jsonl bench/ggo19_needy_top15_run_STRG_442.snap.jsonl

gffcompare -r bench/ggo19_needy_top15_refs/STRG_442.gtf -R -Q -o bench/ggo19_needy_top15_run_STRG_442_cmp bench/ggo19_needy_top15_run_STRG_442.gtf


# STRG.566 refs=STRG.566.12,STRG.566.9,STRG.566.10,STRG.566.11
RUSTLE_TRACE_LOCUS=111146297-111177575 \
  target/release/rustle GGO_19.bam -L --chrom NC_073243.2 \
  -o bench/ggo19_needy_top15_run_STRG_566.gtf \
  --debug-bundle "NC_073243.2:111146297-111177575" --only-debug-bundle \
  --trace-reference bench/ggo19_needy_top15_refs/STRG_566.gtf \
  --trace-output bench/ggo19_needy_top15_run_STRG_566.trace.txt \
  --parity-stage-tsv bench/ggo19_needy_top15_run_STRG_566.parity.tsv \
  --snapshot-jsonl bench/ggo19_needy_top15_run_STRG_566.snap.jsonl

gffcompare -r bench/ggo19_needy_top15_refs/STRG_566.gtf -R -Q -o bench/ggo19_needy_top15_run_STRG_566_cmp bench/ggo19_needy_top15_run_STRG_566.gtf


# STRG.445 refs=STRG.445.2,STRG.445.5,STRG.445.4
RUSTLE_TRACE_LOCUS=80975942-81174569 \
  target/release/rustle GGO_19.bam -L --chrom NC_073243.2 \
  -o bench/ggo19_needy_top15_run_STRG_445.gtf \
  --debug-bundle "NC_073243.2:80975942-81174569" --only-debug-bundle \
  --trace-reference bench/ggo19_needy_top15_refs/STRG_445.gtf \
  --trace-output bench/ggo19_needy_top15_run_STRG_445.trace.txt \
  --parity-stage-tsv bench/ggo19_needy_top15_run_STRG_445.parity.tsv \
  --snapshot-jsonl bench/ggo19_needy_top15_run_STRG_445.snap.jsonl

gffcompare -r bench/ggo19_needy_top15_refs/STRG_445.gtf -R -Q -o bench/ggo19_needy_top15_run_STRG_445_cmp bench/ggo19_needy_top15_run_STRG_445.gtf


# STRG.29 refs=STRG.29.2,STRG.29.3,STRG.29.4
RUSTLE_TRACE_LOCUS=17431800-17453950 \
  target/release/rustle GGO_19.bam -L --chrom NC_073243.2 \
  -o bench/ggo19_needy_top15_run_STRG_29.gtf \
  --debug-bundle "NC_073243.2:17431800-17453950" --only-debug-bundle \
  --trace-reference bench/ggo19_needy_top15_refs/STRG_29.gtf \
  --trace-output bench/ggo19_needy_top15_run_STRG_29.trace.txt \
  --parity-stage-tsv bench/ggo19_needy_top15_run_STRG_29.parity.tsv \
  --snapshot-jsonl bench/ggo19_needy_top15_run_STRG_29.snap.jsonl

gffcompare -r bench/ggo19_needy_top15_refs/STRG_29.gtf -R -Q -o bench/ggo19_needy_top15_run_STRG_29_cmp bench/ggo19_needy_top15_run_STRG_29.gtf


# STRG.440 refs=STRG.440.3,STRG.440.4,STRG.440.5
RUSTLE_TRACE_LOCUS=80273347-80288948 \
  target/release/rustle GGO_19.bam -L --chrom NC_073243.2 \
  -o bench/ggo19_needy_top15_run_STRG_440.gtf \
  --debug-bundle "NC_073243.2:80273347-80288948" --only-debug-bundle \
  --trace-reference bench/ggo19_needy_top15_refs/STRG_440.gtf \
  --trace-output bench/ggo19_needy_top15_run_STRG_440.trace.txt \
  --parity-stage-tsv bench/ggo19_needy_top15_run_STRG_440.parity.tsv \
  --snapshot-jsonl bench/ggo19_needy_top15_run_STRG_440.snap.jsonl

gffcompare -r bench/ggo19_needy_top15_refs/STRG_440.gtf -R -Q -o bench/ggo19_needy_top15_run_STRG_440_cmp bench/ggo19_needy_top15_run_STRG_440.gtf


# STRG.300 refs=STRG.300.1,STRG.300.6,STRG.300.5
RUSTLE_TRACE_LOCUS=44054004-44096254 \
  target/release/rustle GGO_19.bam -L --chrom NC_073243.2 \
  -o bench/ggo19_needy_top15_run_STRG_300.gtf \
  --debug-bundle "NC_073243.2:44054004-44096254" --only-debug-bundle \
  --trace-reference bench/ggo19_needy_top15_refs/STRG_300.gtf \
  --trace-output bench/ggo19_needy_top15_run_STRG_300.trace.txt \
  --parity-stage-tsv bench/ggo19_needy_top15_run_STRG_300.parity.tsv \
  --snapshot-jsonl bench/ggo19_needy_top15_run_STRG_300.snap.jsonl

gffcompare -r bench/ggo19_needy_top15_refs/STRG_300.gtf -R -Q -o bench/ggo19_needy_top15_run_STRG_300_cmp bench/ggo19_needy_top15_run_STRG_300.gtf


# STRG.125 refs=STRG.125.6,STRG.125.9,STRG.125.5
RUSTLE_TRACE_LOCUS=22454887-22475442 \
  target/release/rustle GGO_19.bam -L --chrom NC_073243.2 \
  -o bench/ggo19_needy_top15_run_STRG_125.gtf \
  --debug-bundle "NC_073243.2:22454887-22475442" --only-debug-bundle \
  --trace-reference bench/ggo19_needy_top15_refs/STRG_125.gtf \
  --trace-output bench/ggo19_needy_top15_run_STRG_125.trace.txt \
  --parity-stage-tsv bench/ggo19_needy_top15_run_STRG_125.parity.tsv \
  --snapshot-jsonl bench/ggo19_needy_top15_run_STRG_125.snap.jsonl

gffcompare -r bench/ggo19_needy_top15_refs/STRG_125.gtf -R -Q -o bench/ggo19_needy_top15_run_STRG_125_cmp bench/ggo19_needy_top15_run_STRG_125.gtf


# STRG.210 refs=STRG.210.4,STRG.210.5,STRG.210.1
RUSTLE_TRACE_LOCUS=30118918-30195392 \
  target/release/rustle GGO_19.bam -L --chrom NC_073243.2 \
  -o bench/ggo19_needy_top15_run_STRG_210.gtf \
  --debug-bundle "NC_073243.2:30118918-30195392" --only-debug-bundle \
  --trace-reference bench/ggo19_needy_top15_refs/STRG_210.gtf \
  --trace-output bench/ggo19_needy_top15_run_STRG_210.trace.txt \
  --parity-stage-tsv bench/ggo19_needy_top15_run_STRG_210.parity.tsv \
  --snapshot-jsonl bench/ggo19_needy_top15_run_STRG_210.snap.jsonl

gffcompare -r bench/ggo19_needy_top15_refs/STRG_210.gtf -R -Q -o bench/ggo19_needy_top15_run_STRG_210_cmp bench/ggo19_needy_top15_run_STRG_210.gtf


# STRG.52 refs=STRG.52.2,STRG.52.4
RUSTLE_TRACE_LOCUS=19033042-19077751 \
  target/release/rustle GGO_19.bam -L --chrom NC_073243.2 \
  -o bench/ggo19_needy_top15_run_STRG_52.gtf \
  --debug-bundle "NC_073243.2:19033042-19077751" --only-debug-bundle \
  --trace-reference bench/ggo19_needy_top15_refs/STRG_52.gtf \
  --trace-output bench/ggo19_needy_top15_run_STRG_52.trace.txt \
  --parity-stage-tsv bench/ggo19_needy_top15_run_STRG_52.parity.tsv \
  --snapshot-jsonl bench/ggo19_needy_top15_run_STRG_52.snap.jsonl

gffcompare -r bench/ggo19_needy_top15_refs/STRG_52.gtf -R -Q -o bench/ggo19_needy_top15_run_STRG_52_cmp bench/ggo19_needy_top15_run_STRG_52.gtf


# STRG.95 refs=STRG.95.2,STRG.95.5
RUSTLE_TRACE_LOCUS=20647910-20710598 \
  target/release/rustle GGO_19.bam -L --chrom NC_073243.2 \
  -o bench/ggo19_needy_top15_run_STRG_95.gtf \
  --debug-bundle "NC_073243.2:20647910-20710598" --only-debug-bundle \
  --trace-reference bench/ggo19_needy_top15_refs/STRG_95.gtf \
  --trace-output bench/ggo19_needy_top15_run_STRG_95.trace.txt \
  --parity-stage-tsv bench/ggo19_needy_top15_run_STRG_95.parity.tsv \
  --snapshot-jsonl bench/ggo19_needy_top15_run_STRG_95.snap.jsonl

gffcompare -r bench/ggo19_needy_top15_refs/STRG_95.gtf -R -Q -o bench/ggo19_needy_top15_run_STRG_95_cmp bench/ggo19_needy_top15_run_STRG_95.gtf

