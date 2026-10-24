#!/usr/bin/env bash
# The held-out test of docs/archive/2026-10/merged/PREREG_o3_unexplained_2026-10-15.md: seeds 41-45 (--strands), variants E0 and E2, arms P and T.
#   pixi run bash bench/o3_unexplained/heldout.sh OUT_DIR
set -euo pipefail
OUT=$1; HERE=$(cd "$(dirname "$0")" && pwd)
mkdir -p "$OUT"
for s in 41 42 43 44 45; do
  W=$OUT/s$s
  [ -s "$W/P.fam.copies.tsv" ] || python "$HERE/world.py" build "$W" --seed "$s" --strands
  for v in E0 E2; do for arm in P T; do
    [ -s "$W/$v.$arm/A.assignments.tsv" ] || python "$HERE/world.py" run "$W" --variant "$v" --arm "$arm"
    [ -s "$W/$v.$arm/psvclust.tsv" ] || python "$HERE/psvclust.py" "$W" "$v" "$arm" > "$W/$v.$arm/psvclust.out"
    python "$HERE/score.py" rows "$W" "$v" "$arm" "$OUT/rows.tsv"
  done; done
done
python "$HERE/score.py" summary "$OUT/rows.tsv" | tee "$OUT/summary.txt"
