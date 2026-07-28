#!/usr/bin/env bash
# Anti-over/under-merge lever sweep for --from-genome. Each cell is a full run scored on PARTITION quality
# (over-merge / split / homogeneity / completeness), not just member recall — member recall is partition-blind.
# Levers: min_coverage (repeat/duplicon-bridge defense), min_identity, gamma (quasi-clique density).
set -uo pipefail
BIN=target/debug/gw_family_catalog
FA=/mnt/c/Users/jfris/Desktop/Reference/chm13v2.0.fa
BED=bench/soto/80_fams.chr.bed
OUTD=/home/juanfra/winloci_scratch
run() { # label  identity  coverage  gamma
  local lab=$1 id=$2 cov=$3 g=$4
  local out="$OUTD/dnasweep_$lab"
  RUSTLE_GENOME_MIN_BLOCK=300 RUSTLE_GENOME_MIN_COVERAGE="$cov" RUSTLE_GENOME_GAMMA="$g" \
    "$BIN" --from-genome "$BED" --fasta "$FA" --min-identity "$id" --out "$out" 2>&1 | grep -E "grouping:|reps$|wrote" | sed 's/^/    /'
  printf "%-16s id=%-5s cov=%-5s g=%-5s  " "$lab" "$id" "$cov" "$g"
  python3 bench/soto/partition_score.py "$out.copies.tsv" 2>&1 | tail -1
}
echo "=== ANTI-MERGE LEVER SWEEP (partition quality) ==="
echo "baseline          id=0.90 cov=0.50 g=0.20  (already run)"
python3 bench/soto/partition_score.py "$OUTD/dna_mode.copies.tsv" | sed 's/^/                                          /'
run cov080 0.90 0.80 0.20
run cov090 0.90 0.90 0.20
run id098  0.98 0.50 0.20
run gam050 0.90 0.50 0.50
echo "SWEEP_DONE"
