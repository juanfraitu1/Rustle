#!/usr/bin/env bash
# The held-out test of docs/archive/2026-10/merged/PREREG_o3_psv_clusters_v2_2026-10-15.md: base seeds 61-63, hard seeds 71-75 (--strands), E0 / E2, arm P,
# copy_assign --psv-clusters (the Rust product, scored with score.py --rust), prototype parity on every run.
#   pixi run bash bench/o3_unexplained/heldout_v2.sh OUT_DIR
set -euo pipefail
OUT=$1; HERE=$(cd "$(dirname "$0")" && pwd); BIN=${CARGO_TARGET_DIR:-$HOME/rustle_target}/release
mkdir -p "$OUT"
for spec in base:61 base:62 base:63 hard:71 hard:72 hard:73 hard:74 hard:75; do
  design=${spec%%:*}; s=${spec#*:}; W=$OUT/$design$s
  [ -s "$W/P.fam.copies.tsv" ] || python "$HERE/world.py" build "$W" --seed "$s" --strands --design "$design"
  for v in E0 E2; do
    D=$W/$v.P
    [ -s "$D/A.assignments.tsv" ] || python "$HERE/world.py" run "$W" --variant "$v" --arm P
    [ -s "$D/R.psv_clusters.tsv" ] || "$BIN/copy_assign" --bam "$W/$v.reads.bam" --fasta "$W/genome.ref.fa" --regions "$W/regions.txt" \
        --families "$W/P.fam.copies.tsv" --copies-fa "$W/P.fam.copies.fa" --threads 4 --flag-missing-copies --no-as-tied-only \
        --psv-clusters --out "$D/R" > "$D/R.log" 2>&1
    cmp -s "$D/A.assignments.tsv" "$D/R.assignments.tsv" || echo "ASSIGNMENTS DIFFER $D"
    [ -s "$D/psvclust.tsv" ] || python "$HERE/psvclust.py" "$W" "$v" P > "$D/psvclust.out"
    python "$HERE/parity.py" "$D" "$D/R" || true
    python "$HERE/score.py" rows "$W" "$v" P "$OUT/rows.tsv" --rust
  done
done
python "$HERE/score.py" summary "$OUT/rows.tsv" | tee "$OUT/summary.txt"
