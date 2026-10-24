#!/usr/bin/env bash
# The held-out test of docs/archive/2026-10/merged/PREREG_o3_structure_projection_2026-10-14.md: seeds 21-25 (--strands), variants E0 and E2, the Rust product
# scored by score.py, IsoCon as the reconstruction comparator. Run inside the default pixi environment:
#   pixi run bash bench/o3_structure/heldout.sh OUT_DIR
set -euo pipefail
OUT=$1; HERE=$(cd "$(dirname "$0")" && pwd)
mkdir -p "$OUT"
for s in 21 22 23 24 25; do
  W=$OUT/s$s
  [ -s "$W/P.fam.copies.tsv" ] || python "$HERE/world.py" build "$W" --seed "$s" --strands
  for v in E0 E2; do
    [ -s "$W/$v/P.cand.structure.tsv" ] || python "$HERE/world.py" run "$W" --variant "$v" 2>&1 | grep -v bam_sort
    python "$HERE/score.py" score "$W" "$v" "$W/$v/P.cand" --rows "$OUT/rows.tsv"
    python "$HERE/isocon_cmp.py" "$W" "$v" "$OUT/isocon.tsv"
  done
done
python "$HERE/score.py" summary "$OUT/rows.tsv" --isocon "$OUT/isocon.tsv" | tee "$OUT/summary.txt"
