#!/usr/bin/env bash
# Genome-wide downstream on/off for the contiguous-core family-merge gate (the deferred
# validation from core_gate_pipeline.md). Per-chrom SERIAL (OOM-safe, --vg is heavy), gate
# OFF (default) vs ON @0.13, coordinate-keyed diff (bench/gate_keyed_diff.py) so rayon
# line-order noise is ignored. Confirms genome-wide whether flipping the gate ON loses any
# real transcript (byte-identical-safety prerequisite for a default flip).
# Output: /tmp/gw/onoff_summary.tsv (one row/contig) + /tmp/gw/onoff_$C.diff (details).
set -uo pipefail
BAM=/home/juanfra/winloci_scratch/GGO.bam
FASTA=/home/juanfra/winloci_scratch/GGO.fasta
FAI=/home/juanfra/winloci_scratch/GGO.fasta.fai
RUSTLE=/mnt/c/Users/jfris/Desktop/Rustle/target/release/rustle
SAM=/home/juanfra/miniforge3/bin/samtools
PY=/home/juanfra/miniforge3/bin/python
DIFF=/mnt/c/Users/jfris/Desktop/Rustle/bench/gate_keyed_diff.py
OUT=/tmp/gw
SUM=$OUT/onoff_summary.tsv
echo -e "chrom\ttx_off\ttx_on\tlost\tgained\treal_attr\tcosmetic\tgated" > "$SUM"
CONTIGS=$(awk '$2 > 1000000 {print $1}' "$FAI")
T0=$(date +%s)
for C in $CONTIGS; do
  [[ -f "$OUT/onoff_$C.done" ]] && { echo "[$C] done (skip)"; continue; }
  [[ -s "$OUT/st_$C.gtf" ]] || { echo "[$C] no guide, skip"; continue; }
  CS=$(date +%s)
  "$SAM" view -b "$BAM" "$C" -o "$OUT/$C.bam" 2>/dev/null
  RN=$("$SAM" view -c "$OUT/$C.bam" 2>/dev/null || echo 0)
  if [[ "$RN" -eq 0 ]]; then echo "[$C] 0 reads, skip"; rm -f "$OUT/$C.bam"; continue; fi
  "$SAM" index "$OUT/$C.bam"
  "$SAM" faidx "$FASTA" "$C" > "$OUT/$C.fa" 2>/dev/null && "$SAM" faidx "$OUT/$C.fa"
  RAYON_NUM_THREADS=4 "$RUSTLE" --vg --vg-layer2 --genome-fasta "$OUT/$C.fa" \
      -G "$OUT/st_$C.gtf" -L "$OUT/$C.bam" -o "$OUT/off_$C.gtf" 2>/dev/null
  RAYON_NUM_THREADS=4 RUSTLE_VG_CORE_GATE_TRACE=1 RUSTLE_VG_FAMILY_MIN_CORE_COVERAGE=0.13 \
      "$RUSTLE" --vg --vg-layer2 --genome-fasta "$OUT/$C.fa" \
      -G "$OUT/st_$C.gtf" -L "$OUT/$C.bam" -o "$OUT/on_$C.gtf" 2>"$OUT/on_$C.log"
  GATED=$(grep -c 'would_gate=true' "$OUT/on_$C.log" 2>/dev/null || echo 0)
  ROW=$("$PY" "$DIFF" "$OUT/off_$C.gtf" "$OUT/on_$C.gtf" 2>"$OUT/onoff_$C.diff")
  echo -e "${ROW}\t${GATED}" >> "$SUM"
  echo "[$C] reads=$RN gated=$GATED  ->  $ROW  (wall=$(($(date +%s)-CS))s)"
  rm -f "$OUT/$C.bam" "$OUT/$C.bam.bai" "$OUT/$C.fa" "$OUT/$C.fa.fai"
  touch "$OUT/onoff_$C.done"
done
echo "=== GENOME-WIDE on/off DONE wall=$(($(date +%s)-T0))s ==="
"$PY" - <<'PY'
import csv
rows=list(csv.DictReader(open("/tmp/gw/onoff_summary.tsv"),delimiter="\t"))
S=lambda k:sum(int(r[k]) for r in rows)
print(f"contigs={len(rows)} tx_off={S('tx_off')} tx_on={S('tx_on')} "
      f"LOST={S('lost')} GAINED={S('gained')} real_attr_changed={S('real_attr')} "
      f"cosmetic_renumber={S('cosmetic')} gated_merges={S('gated')}")
PY
echo "ONOFF_GW_DONE"
