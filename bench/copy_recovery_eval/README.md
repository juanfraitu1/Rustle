# VG Copy-Recovery Evaluation Protocol

This protocol measures whether `rustle --vg` recovers RefSeq-annotated paralog gene copies
that StringTie2 misses, evaluated genome-wide (no locus windows) using SQANTI3 structural
classification and a primary-support authenticity guard (minimum primary-aligned reads per
copy) to distinguish genuine recovery from fabrication.

## Run order

`run_protocol.sh` chains stages 00 through 60 in sequence:
00 (setup/index) → 10 (run arms) → 20 (SQANTI3 classification) → 30 (paralog family table) →
40 (per-copy primary-support guard) → 50 (summary statistics) → 60 (final report).

## Configuration

All pinned paths, tool locations, arm flags, and thresholds live in `config.sh`.
Source it at the top of every stage script: `source "$(dirname "$0")/config.sh"`.
Override any variable via environment before invoking a stage (e.g. `CHROM_SUBSET=chr19`
for smoke tests).

## Headline result

To be recorded here after the genome-wide run completes (Task 10).
