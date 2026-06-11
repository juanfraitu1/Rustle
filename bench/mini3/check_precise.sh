#!/usr/bin/env bash
# Invariant: RUSTLE_PRECISE=1 must byte-match the frozen 4705ab1 behavior.
set -euo pipefail
cd "$(git rev-parse --show-toplevel)"
OUT=$(mktemp --suffix=.gtf)
RUSTLE_PRECISE=1 target/release/rustle -L bench/mini3/mini3.bam -o "$OUT" 2>/dev/null
if diff <(grep -v '^#' bench/mini3/expected_precise.gtf) <(grep -v '^#' "$OUT") >/dev/null; then
  echo "ESCAPE-HATCH OK: RUSTLE_PRECISE=1 byte-matches 4705ab1"
else
  echo "ESCAPE-HATCH DRIFT — RUSTLE_PRECISE=1 diverged from 4705ab1"; diff <(grep -v '^#' bench/mini3/expected_precise.gtf) <(grep -v '^#' "$OUT") | head; exit 1
fi
rm -f "$OUT"
