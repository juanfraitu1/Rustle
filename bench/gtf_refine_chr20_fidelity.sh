#!/bin/bash
# chr20 fidelity anchors for --gtf-refine and the dedup fix
# (docs/superpowers/specs/2026-09-16-gtf-refine-and-dedup-fix-design.md, "Fidelity anchors"). Serial, foreground.
#
# Usage: bash bench/gtf_refine_chr20_fidelity.sh [arm_label]
#   With no argument, runs all seven arms in sequence (each is a whole-chr20 run,
#   can take several minutes). With an arm_label, runs only that one arm (for
#   restarting a single long run under the tool's background/foreground split).
set -euo pipefail
W=/mnt/linuxdisk/home/juanfraitu/bakeoff/human_chr20
F=$W/fidelity
BIN=/mnt/linuxdisk/home/juanfraitu/rustle_target/release/copy_assign
mkdir -p "$F"
arm() { # $1 label, $2 legacy(0/1), $3 refine list ('' = none)
  local label=$1 leg=$2 refine=$3
  mkdir -p "$F/$label"
  local extra=()
  [ -n "$refine" ] && extra=(--gtf-refine "$refine")
  ( cd "$F/$label"
    if [ "$leg" = 1 ]; then export RUSTLE_LEGACY_PLACEMENT_DEDUP=1; fi
    "$BIN" --gtf --bam "$W/chr20.bam" --fasta "$W/chr20.fa" --region chr20:1-66210255 "${extra[@]}" --out ours \
      > ours.stdout 2> ours.stderr )
  ( cd "$F/$label" && gffcompare -r "$W/chr20_ref.gtf" -o gffc ours.gtf > /dev/null 2>&1 )
  echo "== $label"; sed -n '/Query mRNAs/p;/Transcript level/p;/Intron chain level/p;/Locus level/p;/Matching transcripts/p' "$F/$label/gffc.stats"
}
run_arm() {
  case "$1" in
    legacy_none) arm legacy_none 1 '' ;;
    fixed_none) arm fixed_none 0 '' ;;
    fixed_fragsupport) arm fixed_fragsupport 0 fragsupport ;;
    fixed_tss) arm fixed_tss 0 tss ;;
    legacy_subset) arm legacy_subset 1 subset ;;
    legacy_strand) arm legacy_strand 1 strand ;;
    legacy_strand_subset_mono) arm legacy_strand_subset_mono 1 strand,subset,mono ;;
    fixed_all) arm fixed_all 0 all ;;
    *) echo "unknown arm: $1" >&2; exit 1 ;;
  esac
}
ALL_ARMS=(legacy_none fixed_none fixed_fragsupport legacy_subset legacy_strand legacy_strand_subset_mono fixed_all)
if [ "$#" -ge 1 ]; then
  run_arm "$1"
else
  for a in "${ALL_ARMS[@]}"; do run_arm "$a"; done
  cmp -s "$F/legacy_none/ours.gtf" "$W/ours/ours.gtf" && echo "legacy_none GTF byte-identical to the 2026-09-15 ours.gtf" || echo "legacy_none GTF DIFFERS from the 2026-09-15 ours.gtf"
fi
