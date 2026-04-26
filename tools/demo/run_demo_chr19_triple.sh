#!/usr/bin/env bash
# Tier 2 — chr19 GOLGA6L paralog triple from GGO_19.bam.
# Three loci on NC_073243.2 — all annotated "golgin subfamily A member
# 6-like protein 7", all + strand, all 8 introns:
#   LOC115931294 104789647-104796276
#   LOC134757625 104830536-104837094
#   LOC101137218 104871356-104877901
set -euo pipefail

cd "$(dirname "$0")/../.."

OUT="tools/demo/out/chr19_triple"
mkdir -p "$OUT"

BAM="/mnt/c/Users/jfris/Desktop/GGO_19.bam"
LOCUS_REGION="NC_073243.2:104780000-104900000"

LOCUS_A="NC_073243.2:104789647-104796276"   # LOC115931294
LOCUS_B="NC_073243.2:104830536-104837094"   # LOC134757625
LOCUS_C="NC_073243.2:104871356-104877901"   # LOC101137218

# slice BAM to the region for fast runs
REGION_BAM="$OUT/region.bam"
echo "== Slicing BAM to $LOCUS_REGION =="
samtools view -bh "$BAM" "$LOCUS_REGION" > "$REGION_BAM"
samtools index "$REGION_BAM"
samtools view -c "$REGION_BAM"  # report read count

NODES_TSV="$OUT/nodes.tsv"
EDGES_TSV="$OUT/edges.tsv"
FLOW_TSV="$OUT/flow.tsv"
READS_TSV="$OUT/reads.tsv"
PARTITION_TSV="$OUT/partition.tsv"
VG_LOG="$OUT/vg_trace.log"

GTF_BASELINE="$OUT/baseline.gtf"
GTF_VG="$OUT/vg.gtf"

rm -f "$NODES_TSV" "$EDGES_TSV" "$FLOW_TSV" "$READS_TSV" "$PARTITION_TSV"

echo "== Rustle baseline on chr19 triple region =="
RUSTLE_PARITY_GRAPH_TSV="$NODES_TSV" \
RUSTLE_PARITY_GRAPH_EDGES_TSV="$EDGES_TSV" \
RUSTLE_FLOW_ITER_TSV="$FLOW_TSV" \
RUSTLE_PARITY_READ_DUMP_TSV="$READS_TSV" \
RUSTLE_PARITY_PARTITION_TSV="$PARTITION_TSV" \
./target/release/rustle -L -o "$GTF_BASELINE" "$REGION_BAM"

echo "== Rustle --vg on same region =="
RUSTLE_VG_TRACE=1 \
./target/release/rustle -L --vg -o "$GTF_VG" "$REGION_BAM" 2> "$VG_LOG"

echo "== Bundles overlapping each paralog locus =="
for L in "$LOCUS_A" "$LOCUS_B" "$LOCUS_C"; do
  echo "-- $L"
  python3 tools/demo/identify_locus_bundle.py --partition-tsv "$PARTITION_TSV" --locus "$L" || true
done

echo "== Rendering V1 individual splice graphs =="
python3 tools/demo/render_splice_graph.py --nodes-tsv "$NODES_TSV" --edges-tsv "$EDGES_TSV" \
  --locus "$LOCUS_A" --output-png "$OUT/v1_locA.png" --output-dot "$OUT/v1_locA.dot" \
  --title "LOC115931294" || echo "  (V1 locA failed; locus may not have a bundle)"
python3 tools/demo/render_splice_graph.py --nodes-tsv "$NODES_TSV" --edges-tsv "$EDGES_TSV" \
  --locus "$LOCUS_B" --output-png "$OUT/v1_locB.png" --output-dot "$OUT/v1_locB.dot" \
  --title "LOC134757625" || echo "  (V1 locB failed)"
python3 tools/demo/render_splice_graph.py --nodes-tsv "$NODES_TSV" --edges-tsv "$EDGES_TSV" \
  --locus "$LOCUS_C" --output-png "$OUT/v1_locC.png" --output-dot "$OUT/v1_locC.dot" \
  --title "LOC101137218" || echo "  (V1 locC failed)"

echo "== Rendering V2 three-way isomorphism =="
python3 tools/demo/compare_graphs_isomorphism.py \
  --nodes-tsv "$NODES_TSV" --edges-tsv "$EDGES_TSV" \
  --locus "$LOCUS_A" --locus "$LOCUS_B" --locus "$LOCUS_C" \
  --output-json "$OUT/v2_isomorphism.json" \
  --output-png "$OUT/v2_isomorphism.png" || echo "  (V2 failed; check that all 3 loci have bundles)"

if [[ -f "$OUT/v2_isomorphism.json" ]]; then
  echo "  V2 result:"
  python3 -c "
import json
d = json.load(open('$OUT/v2_isomorphism.json'))
print('  isomorphic =', d['isomorphic'])
print('  pairwise =', list(d['pairwise'].keys()))
for s in d.get('graph_summary', []):
    print('   ', s)
"
fi

echo "== Rendering V3 flow iterations for locus A =="
python3 tools/demo/render_flow_iterations.py \
  --flow-tsv "$FLOW_TSV" --nodes-tsv "$NODES_TSV" --edges-tsv "$EDGES_TSV" \
  --locus "$LOCUS_A" --output-png "$OUT/v3_flow_locA.png" || echo "  (V3 failed)"

echo "== Rendering V4 multi-mapper pinning (best-effort) =="
READ_ID=$(grep -m1 -oE 'read_id=\S+' "$VG_LOG" | head -1 | cut -d= -f2 || echo "")
if [[ -n "$READ_ID" ]]; then
  python3 tools/demo/render_multimapper_pinning.py \
    --vg-trace-log "$VG_LOG" --read-id "$READ_ID" \
    --output-png "$OUT/v4_multimapper_${READ_ID}.png" 2>&1 || \
    echo "  (V4 skipped — no per-EM-iter rows in current binary trace)"
else
  echo "  (V4 skipped — current binary emits summary lines only)"
fi

echo "== Side-by-side baseline-vs-VG transcript counts =="
echo "Baseline transcripts: $(grep -cE '^[^#].*\ttranscript\t' "$GTF_BASELINE" 2>/dev/null || echo 0)"
echo "VG transcripts:       $(grep -cE '^[^#].*\ttranscript\t' "$GTF_VG" 2>/dev/null || echo 0)"

echo "== Tier 2 complete. Outputs in $OUT =="
ls -lh "$OUT"
