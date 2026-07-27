#!/usr/bin/env bash
# --from-genome benchmark: reproduce the Soto 83-family set from the CHM13 genome alone, score vs Soto.
# Uses the debug binary (the heavy work is external minimap2, so debug orchestration is fine) to avoid a
# slow release build. One minimap2 genome index (~few min) + 362 short window self-alignments.
set -euo pipefail
BIN=target/debug/gw_family_catalog
FA=/mnt/c/Users/jfris/Desktop/Reference/chm13v2.0.fa
BED=bench/soto/80_fams.chr.bed
OUT=/home/juanfra/winloci_scratch/dna_mode
cargo build --bin gw_family_catalog 2>&1 | tail -3
# min_block 300 to admit the smaller Soto members; SD-discovery identity floor 0.90.
RUSTLE_GENOME_MIN_BLOCK=300 "$BIN" --from-genome "$BED" --fasta "$FA" --min-identity 0.90 --out "$OUT"
echo "=== SCORE ==="
python3 bench/soto/score_from_genome.py "$OUT.copies.tsv"
echo "RUN_FROM_GENOME_DONE"
