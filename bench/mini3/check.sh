#!/usr/bin/env bash
# Fast ST-faithful convergence check on the 3-loci mini fixture (~0.1s).
# Usage: bench/mini3/check.sh            -> rustle (default) vs StringTie
#        RUSTLE_PRECISE=1 bench/mini3/check.sh  -> escape-hatch (must match commit 4705ab1)
set -euo pipefail
cd "$(git rev-parse --show-toplevel)"
BAM=bench/mini3/mini3.bam
EXP=bench/mini3/expected_st.gtf
OUT=$(mktemp --suffix=.gtf)
target/release/rustle -L "$BAM" -o "$OUT" 2>/dev/null
echo "rustle-vs-StringTie on mini3 (target: 0 / 0):"
python3 bench/gtf_chain_diff.py "$OUT" "$EXP" 2>/dev/null | grep -E "in both|Rustle-only|ST-only" | head -1
rm -f "$OUT"
