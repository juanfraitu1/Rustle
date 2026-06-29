#!/usr/bin/env bash
# Full protocol driver: U1 -> U3 -> (U2 + U4 per arm) -> U5 -> U6.
# Deterministic; no LLM in the loop. Python stages run via $PY (an interpreter
# with pysam/pandas; base env lacks them).
source "$(dirname "$0")/config.sh"
set -euo pipefail
D="$(dirname "$0")"

# U1 — transcriptomes (StringTie2 + rustle arms).
bash "$D/10_generate_transcriptomes.sh"

# U3 — paralog universe (tool-independent; once). Produces universe.tsv + exons.tsv.
$PY "$D/30_build_universe.py" --ref-gff "$REF_GFF" --genome "$GENOME_FASTA" \
  --min-identity "$MIN_IDENTITY" --min-cov-frac "$MIN_COV_FRAC" \
  --minimap2 "$MINIMAP2" --gffread "$GFFREAD" --outdir "$OUTDIR"

# U2 + U4 per arm (StringTie2 baseline + rustle-VG headline).
for arm in stringtie rustle_vg; do
  GTF="$OUTDIR/gtf/${arm}.gtf"
  [ -s "$GTF" ] || { echo "[driver] skip $arm (no GTF)"; continue; }
  bash "$D/20_run_sqanti3.sh" "$arm" "$GTF"
  $PY "$D/40_recovery_sets.py" \
    --classification "$OUTDIR/sqanti3/$arm/${arm}_classification.txt" \
    --universe "$OUTDIR/universe.tsv" --out "$OUTDIR/${arm}_recovery.tsv"
done

# U5 — authenticity guard (rustle-VG recoveries only): decisive own-copy evidence
# vs sister copies -> authentic / unresolvable / phantom.
$PY "$D/50_authenticity_guard.py" --bam "$BAM" \
  --recovery "$OUTDIR/rustle_vg_recovery.tsv" --exons "$OUTDIR/exons.tsv" \
  --universe "$OUTDIR/universe.tsv" \
  --k-decisive "$GUARD_K_DECISIVE" --k-tied "$GUARD_K_TIED" \
  --out "$OUTDIR/authenticity.tsv"

# U6 — head-to-head + report.
$PY "$D/60_headtohead.py" \
  --rustle-recovery "$OUTDIR/rustle_vg_recovery.tsv" \
  --stringtie-recovery "$OUTDIR/stringtie_recovery.tsv" \
  --authenticity "$OUTDIR/authenticity.tsv" --outdir "$OUTDIR"

echo "[DONE] see $OUTDIR/REPORT.md and $OUTDIR/headtohead.tsv"
