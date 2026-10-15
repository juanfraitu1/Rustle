#!/usr/bin/env bash
# missing_copy_flag positive control (d = 0.02, 0.01; seeds 1-5) and negative control (template reads only, d = 0.02),
# docs/PREREG_missing_copy_posctl_2026-10-14.md. Run inside the pixi environment (minimap2 2.30, samtools 1.22.1):
#   pixi run bash bench/o3_posctl.sh REF_DIR OUT_DIR [BIN_DIR]
# REF_DIR holds chm13_chr20.fa, chr20_ref.gtf, chm13_chr20.RefSeq.gff, chm13_chr20.splice.mmi (see the prereg).
set -euo pipefail
REF=$1; OUT=$2; BIN=${3:-${CARGO_TARGET_DIR:-$HOME/rustle_target}/release}
HERE=$(cd "$(dirname "$0")" && pwd)
MM2=$(python3 -c "import sys; sys.path.insert(0, '$HERE'); import sim; print(sim.MM2)")
mkdir -p "$OUT"; cd "$OUT"
flag() {  # flag BAM SIM_PREFIX RUN_PREFIX
  "$BIN/missing_copy_flag" --bam "$1" --fasta "$REF/chm13_chr20.fa" --loci "$REF/chr20_ref.gtf" \
    --index "$REF/chm13_chr20.splice.mmi" --gff "$REF/chm13_chr20.RefSeq.gff" --confirm copies="$2.copies.mmi" \
    --threads 4 --out "$3" 2> "$3.log"
}
for d in 0.02 0.01; do
  for s in 1 2 3 4 5; do
    sim=sim_d${d}_s${s}; run=run_d${d}_s${s}
    [ -s "$sim.bam" ] || python3 -I "$HERE/sim.py" missing-copy genomic "$REF/chr20_ref.gtf" "$REF/chm13_chr20.fa" chr20 40 "$d" \
      "$sim" "$s" --threads 4 > "$sim.log" 2>&1
    [ -s "$run.missing_copy.tsv" ] || flag "$sim.bam" "$sim" "$run"
    python3 -I "$HERE/score.py" missing-copy-control --sim "$sim" --run "$run" --label "$run" --per-gene "$run.per_gene.tsv"
    if [ "$d" = 0.02 ]; then  # N: the extra copy's reads removed
      neg=neg_d${d}_s${s}
      if [ ! -s "$neg.bam" ]; then
        awk 'NR % 4 == 1 { keep = ($0 ~ /^@template\|/) } keep' "$sim.fq" > "$neg.fq"
        $MM2 -t 4 "$REF/chm13_chr20.fa" "$neg.fq" 2>/dev/null | samtools sort -@2 -o "$neg.bam" - && samtools index "$neg.bam"
      fi
      [ -s "$neg.missing_copy.tsv" ] || flag "$neg.bam" "$sim" "$neg"
      python3 -I "$HERE/score.py" missing-copy-control --sim "$sim" --run "$neg" --label "$neg" --per-gene "$neg.per_gene.tsv"
    fi
  done
done
