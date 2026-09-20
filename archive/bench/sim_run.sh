#!/usr/bin/env bash
# Reproduce the fully-simulated ground-truth benchmark end-to-end:
#   sim_genome.py (plant genome+reads) -> minimap2 -> gw_family_catalog (O1) -> copy_assign (O2) -> sim_eval.py
# Everything synthetic and labelled => NON-circular: the read name carries the TRUE family/copy.
set -euo pipefail

OUT=/home/juanfra/winloci_scratch
PY=/home/juanfra/miniforge3/bin/python
MM2=/home/juanfra/miniforge3/bin/minimap2
SAM=/home/juanfra/miniforge3/bin/samtools
BENCH="$(cd "$(dirname "$0")" && pwd)"
REL="$BENCH/../target/release"

echo "== 1. plant genome + labelled reads =="
"$PY" "$BENCH/sim_genome.py"

echo "== 2. map (keep secondaries so co-located copies multimap => MAPQ-0) =="
"$MM2" -ax splice:hq -uf --eqx -Y -N 50 -p 0.1 --secondary=yes "$OUT/simgw.fasta" "$OUT/simgw_reads.fastq" 2>/dev/null \
  | "$SAM" sort -o "$OUT/simgw.bam"
"$SAM" index "$OUT/simgw.bam"

# The simulated copies sit on two separate contigs (simA/simB), so this exercises the dispersed case. That
# needs no flag: the default homology (E_r) catalog never restricts membership by chromosome.
echo "== 3. O1: family catalog (homology/E_r, >=2 copies) =="
printf 'simA:0-200000\nsimB:0-200000\n' > "$OUT/simgw_regions.txt"
"$REL/gw_family_catalog" --bam "$OUT/simgw.bam" --fasta "$OUT/simgw.fasta" \
  --min-copies 2 --out "$OUT/simgw_fam" 2>/dev/null

echo "== 4. O2: per-read copy assignment (de-novo families + posterior/zone) =="
"$REL/copy_assign" --bam "$OUT/simgw.bam" --fasta "$OUT/simgw.fasta" \
  --regions "$OUT/simgw_regions.txt" --min-copies 2 --posterior \
  --out "$OUT/simgw_as" 2>/dev/null

echo "== 5. score vs PLANTED truth =="
"$PY" "$BENCH/sim_eval.py"
