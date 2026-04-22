#!/usr/bin/env bash
# Full-BAM trace vs precomputed StringTie GTF (classify NotExtracted vs Filter for miss refs).
# Usage: from repo root, after building release rustle:
#   ./scripts/run_full_bam_stringtie_trace.sh /path/to/GGO_19.bam /path/to/GGO_19_stringtie.gtf out_dir
set -euo pipefail
BAM="${1:?Usage: $0 BAM STRINGTIE_GTF OUT_DIR}"
ST_GTF="${2:?}"
OUT_DIR="${3:?}"
mkdir -p "$OUT_DIR"
RUSTLE="${RUSTLE:-./target/release/rustle}"
exec "$RUSTLE" -L -p "${RUSTLE_THREADS:-8}" \
  --trace-reference "$ST_GTF" \
  --trace-output "$OUT_DIR/trace_vs_stringtie.txt" \
  "$BAM" \
  >"$OUT_DIR/rustle_trace.stdout.log" 2>"$OUT_DIR/rustle_trace.stderr.log"
