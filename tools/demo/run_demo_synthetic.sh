#!/usr/bin/env bash
# Tier 1 — synthetic_family walkthrough.
# Produces all dumps + V1 + V3 + V4 PNGs for the 95-read synthetic dataset.
set -euo pipefail

cd "$(dirname "$0")/../.."

OUT="tools/demo/out/synthetic"
mkdir -p "$OUT"

BAM="test_data/synthetic_family/reads_sorted.bam"
GTF_TRUTH="test_data/synthetic_family/truth.gtf"

GTF_BASELINE="$OUT/baseline.gtf"
GTF_VG="$OUT/vg.gtf"

NODES_TSV="$OUT/nodes.tsv"
EDGES_TSV="$OUT/edges.tsv"
FLOW_TSV="$OUT/flow.tsv"
READS_TSV="$OUT/reads.tsv"
PARTITION_TSV="$OUT/partition.tsv"
VG_LOG="$OUT/vg_trace.log"

# clean prior dumps so headers are fresh
rm -f "$NODES_TSV" "$EDGES_TSV" "$FLOW_TSV" "$READS_TSV" "$PARTITION_TSV"

echo "== Running Rustle baseline (no VG) =="
RUSTLE_PARITY_GRAPH_TSV="$NODES_TSV" \
RUSTLE_PARITY_GRAPH_EDGES_TSV="$EDGES_TSV" \
RUSTLE_FLOW_ITER_TSV="$FLOW_TSV" \
RUSTLE_PARITY_READ_DUMP_TSV="$READS_TSV" \
RUSTLE_PARITY_PARTITION_TSV="$PARTITION_TSV" \
./target/release/rustle -L -o "$GTF_BASELINE" "$BAM"

echo "== Running Rustle --vg (capturing trace) =="
RUSTLE_VG_TRACE=1 \
./target/release/rustle -L --vg -o "$GTF_VG" "$BAM" 2> "$VG_LOG"

echo "== Per-stage dump line counts =="
wc -l "$READS_TSV" "$PARTITION_TSV" "$NODES_TSV" "$EDGES_TSV" "$FLOW_TSV" "$GTF_BASELINE" "$GTF_VG"

echo "== Rendering V1 splice graph (first bundle) =="
FIRST_LOCUS=$(awk -F'\t' 'NR>1 {print $2":"$3"-"$4; exit}' "$NODES_TSV")
echo "  Locus = $FIRST_LOCUS"
python3 tools/demo/render_splice_graph.py \
  --nodes-tsv "$NODES_TSV" --edges-tsv "$EDGES_TSV" \
  --locus "$FIRST_LOCUS" \
  --output-png "$OUT/v1_splice_graph.png" \
  --output-dot "$OUT/v1_splice_graph.dot" \
  --title "Tier 1 — synthetic_family bundle 1"

echo "== Rendering V3 flow iterations =="
python3 tools/demo/render_flow_iterations.py \
  --flow-tsv "$FLOW_TSV" --nodes-tsv "$NODES_TSV" --edges-tsv "$EDGES_TSV" \
  --locus "$FIRST_LOCUS" \
  --output-png "$OUT/v3_flow_iterations.png"

echo "== Rendering V4 multi-mapper pinning (if real-data trace has per-iter rows) =="
# The current rustle binary emits only [VG] summary lines, not per-EM-iter
# weights, so this typically fails — that's expected. The renderer will exit
# with a clean error message.
READ_ID=$(grep -m1 -oE 'read_id=[^ ]+' "$VG_LOG" | head -1 | cut -d= -f2 || echo "")
if [[ -n "$READ_ID" ]]; then
  python3 tools/demo/render_multimapper_pinning.py \
    --vg-trace-log "$VG_LOG" --read-id "$READ_ID" \
    --output-png "$OUT/v4_multimapper_${READ_ID}.png" 2>&1 || \
    echo "  (V4 skipped — no per-EM-iter rows in current binary trace; fall back to fixture)"
else
  echo "  (V4 skipped — no read_id matches in VG trace; current binary emits summary only)"
fi

echo "== Tier 1 complete. Outputs in $OUT =="
ls -lh "$OUT"
