#!/bin/bash
# rustle_reassemble.sh — the closed loop's pass 2, one command per step (docs/archive/2026-09/PREREG_tied_read_loop_2026-09-25.md).
#
#   union   copy assignment of tied reads with ONE test per read over every candidate copy (the union certificate),
#           on the default families' units            copy_assign --families UNITS --union-certificate -> PREFIX.loop.*
#   home    each ASSIGNED read's home copy (tied / ambiguous reads are never moved; a read with no alignment at its
#           home is left out and counted)             bench/loop_home.py union               -> PREFIX.loop.home.tsv
#   pass2   the assemble stage again, each listed read taken only at its home copy (all its other alignments
#           dropped, the primary included; every other read as in pass 1)
#                                                     copy_assign --assemble-only ... with RUSTLE_READ_HOME_TABLE
#                                                     -> PREFIX.pass2.gtf (--arm S, base = pass 1's seeding)
#                                                     -> PREFIX.pass2p.gtf (--arm P, base = primaries only)
#   g0      the pre-registered stop-rule counts (records moved) from the pass-2 log            -> PREFIX.<tag>.g0.tsv
#   all     union, home, pass2, g0
#
# usage: tools/rustle_reassemble.sh STAGE --bam B --fasta G --out PREFIX [--units U.tsv --units-fa U.fa]
#        [--arm S|P] [--home HOME.tsv] [--tag NAME] [--threads N] [--bin DIR] [--python PY]
#
# Inputs from the pipeline driver (tools/rustle_pipeline.sh, same PREFIX): PREFIX.molecules.tsv (the genome-wide
# best-AS table of `assemble`) and the families' units PREFIX.fam.units.tsv/.fa (the `families` stage run with --bam,
# which writes units in the copy_assign --families contract). --units/--units-fa name another roster (the legacy
# catalog PREFIX.cat.copies.tsv/.fa is for development mechanics only; the prereg forbids a verdict on it).
# --home replaces the union table (e.g. the ORACLE table of a simulation, bench/loop_home.py oracle); --tag names the
# pass-2 product (default pass2 / pass2p). The families are frozen from pass 1: nothing here recomputes them.
# Heavy steps: `union` is one copy_assign over every contig (the driver's `assign` cost: est. minutes to hours per
# sample); `pass2` costs about one pass-1 assembly. Run each under the machine lock, foreground.
set -euo pipefail
STAGE=${1:-all}; shift || true
HERE="$(cd "$(dirname "$0")" && pwd)"
BAM=""; FASTA=""; OUT=""; UNITS=""; UNITS_FA=""; ARM=S; HOME_TSV=""; TAG=""; THREADS=4
BIN="$HERE/../target/release"; PY=python3
while [ $# -gt 0 ]; do
  case "$1" in
    --bam) BAM=$2; shift 2;; --fasta) FASTA=$2; shift 2;; --out) OUT=$2; shift 2;;
    --units) UNITS=$2; shift 2;; --units-fa) UNITS_FA=$2; shift 2;; --arm) ARM=$2; shift 2;;
    --home) HOME_TSV=$2; shift 2;; --tag) TAG=$2; shift 2;; --threads) THREADS=$2; shift 2;;
    --bin) BIN=$2; shift 2;; --python) PY=$2; shift 2;;
    *) echo "unknown argument $1" >&2; exit 2;;
  esac
done
[ -n "$BAM" ] && [ -n "$FASTA" ] && [ -n "$OUT" ] || { echo "need --bam, --fasta, --out" >&2; exit 2; }
case "$ARM" in S|P) ;; *) echo "--arm is S (pass-1 seeding as base) or P (primaries only as base)" >&2; exit 2;; esac
UNITS=${UNITS:-$OUT.fam.units.tsv}; UNITS_FA=${UNITS_FA:-$OUT.fam.units.fa}
[ -n "$TAG" ] || { [ "$ARM" = S ] && TAG=pass2 || TAG=pass2p; }
export TMPDIR=${TMPDIR:-/tmp}
# the pass-1 assembly's polish flags, read from the driver so the two passes can never drift apart
POLISH=$(sed -n 's/^POLISH="\(.*\)"$/\1/p' "$HERE/rustle_pipeline.sh")
[ -n "$POLISH" ] || { echo "could not read POLISH from $HERE/rustle_pipeline.sh" >&2; exit 2; }
say() { echo "[rustle_reassemble] $(date +%H:%M:%S) $*" >&2; }

stage_union() {
  for f in "$UNITS" "$UNITS_FA"; do
    [ -s "$f" ] || { echo "missing $f: run the families stage with --bam (units in the copy_assign contract), or pass --units/--units-fa" >&2; exit 2; }
  done
  say "union: copy assignment with the union certificate on $UNITS"
  samtools view -H "$BAM" | awk '$1=="@SQ"{sub("SN:","",$2); sub("LN:","",$3); print $2":1-"$3}' > "$OUT.loop.regions.txt"
  "$BIN/copy_assign" --bam "$BAM" --fasta "$FASTA" --regions "$OUT.loop.regions.txt" \
    --families "$UNITS" --copies-fa "$UNITS_FA" --union-certificate --out "$OUT.loop" > "$OUT.loop.log" 2>&1
  say "union: $(grep -o '\[union\] --union-certificate: .*' "$OUT.loop.log" | head -1 | cut -c1-220)"
}
stage_home() {
  [ -s "$OUT.loop.union_certificate.tsv" ] || { echo "missing $OUT.loop.union_certificate.tsv: run the union step" >&2; exit 2; }
  say "home: assigned reads -> home copies"
  "$PY" "$HERE/../bench/loop_home.py" union --union "$OUT.loop.union_certificate.tsv" --families "$UNITS" \
    --bam "$BAM" --out "$OUT.loop.home.tsv" --summary "$OUT.loop.home.summary.tsv" 2> "$OUT.loop.home.log"
  say "home: $(awk -F'\t' '$1=="A_assigned"||$1=="H_no_record_at_home"||$1=="home_table_molecules"{printf "%s=%s ", $1, $2}' "$OUT.loop.home.summary.tsv")"
}
stage_pass2() {
  local table=${HOME_TSV:-$OUT.loop.home.tsv}
  [ -f "$table" ] || { echo "missing $table: run the home step (or pass --home)" >&2; exit 2; }
  local seed_env=()
  if [ "$ARM" = S ]; then
    [ -s "$OUT.molecules.tsv" ] || { echo "missing $OUT.molecules.tsv: run the driver's assemble stage first (same PREFIX)" >&2; exit 2; }
    seed_env=(RUSTLE_GTF_SECONDARY=1 RUSTLE_GTF_SECONDARY_AS_RATIO=0.98 "RUSTLE_GTF_SECONDARY_AS_TABLE=$OUT.molecules.tsv")
  fi
  say "pass2 ($ARM): re-assembly with $table -> $OUT.$TAG.gtf"
  env "${seed_env[@]}" "RUSTLE_READ_HOME_TABLE=$table" "$BIN/copy_assign" --assemble-only --genome-wide \
    --assembly-junctions strict $POLISH --gtf-tpm --bam "$BAM" --fasta "$FASTA" --out "$OUT.$TAG" > "$OUT.$TAG.assemble.log" 2>&1
  say "pass2: $(awk -F'\t' '$3=="transcript"' "$OUT.$TAG.gtf" | wc -l) transcripts"
}
stage_g0() {
  "$PY" "$HERE/../bench/loop_home.py" g0 --log "$OUT.$TAG.assemble.log" --home "${HOME_TSV:-$OUT.loop.home.tsv}" > "$OUT.$TAG.g0.tsv"
  say "g0: $(tr '\n' ' ' < "$OUT.$TAG.g0.tsv")"
}
case "$STAGE" in
  union) stage_union;; home) stage_home;; pass2) stage_pass2;; g0) stage_g0;;
  all) stage_union; stage_home; stage_pass2; stage_g0;;
  *) echo "unknown stage $STAGE" >&2; exit 2;;
esac
say "done"
