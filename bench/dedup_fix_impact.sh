#!/bin/bash
# Dedup-fix impact measurement (docs/superpowers/specs/2026-09-16-gtf-refine-and-dedup-fix-design.md, Part 1).
# Three arms per substrate: base (pre-fix binary), fixed (new binary, default), legacy (new binary,
# RUSTLE_LEGACY_PLACEMENT_DEDUP=1). Serial, foreground, outputs under /mnt/linuxdisk.
#
# Usage:
#   bash bench/dedup_fix_impact.sh                # run all substrates x all arms
#   bash bench/dedup_fix_impact.sh SUBSTRATE       # run one substrate, all arms
#   bash bench/dedup_fix_impact.sh SUBSTRATE ARM   # run exactly one (substrate, arm) job
# SUBSTRATE in {o2_families, denovo, chr20_gtf}; ARM in {base, fixed, legacy}.
# Each (substrate, arm) run is independent -- callers under a wall-clock cap (e.g. a 10-minute
# tool call) should invoke this once per (substrate, arm) pair rather than the bare no-arg form.
set -euo pipefail
O=/mnt/linuxdisk/home/juanfraitu/dedupfix
BASE=$O/copy_assign.base
NEW=/mnt/linuxdisk/home/juanfraitu/rustle_target/release/copy_assign
BAM=/mnt/linuxdisk/home/juanfraitu/npip_cat/npip3.bam
FA=/mnt/linuxdisk/home/juanfraitu/_from_wsl/winloci_scratch/GGO.fasta
CAT=/mnt/linuxdisk/home/juanfraitu/mec/batch.copies.tsv
CFA=/mnt/linuxdisk/home/juanfraitu/npip_cat/arm_f2/cat.copies.fa
REG=/mnt/linuxdisk/home/juanfraitu/mec/regions.txt
C20=/mnt/linuxdisk/home/juanfraitu/bakeoff/human_chr20

run() { # $1=substrate $2=arm $3=binary $4=legacy(0/1) rest=args
  local sub=$1 arm=$2 bin=$3 leg=$4; shift 4
  mkdir -p "$O/$sub/$arm"
  if [ "$leg" = 1 ]; then
    ( cd "$O/$sub/$arm" && RUSTLE_LEGACY_PLACEMENT_DEDUP=1 "$bin" "$@" --out run > run.stdout 2> run.stderr )
  else
    ( cd "$O/$sub/$arm" && "$bin" "$@" --out run > run.stdout 2> run.stderr )
  fi
  echo "$sub/$arm exit=$?"
}

run_one() { # $1=substrate $2=arm
  local sub=$1 arm=$2
  local bin=$NEW leg=0
  [ "$arm" = base ] && bin=$BASE
  [ "$arm" = legacy ] && leg=1
  case "$sub" in
    o2_families) run o2_families "$arm" "$bin" "$leg" --bam $BAM --fasta $FA --families $CAT --copies-fa $CFA --regions $REG --dump-psv ;;
    denovo)      run denovo "$arm" "$bin" "$leg" --bam $BAM --fasta $FA --regions $REG ;;
    chr20_gtf)   run chr20_gtf "$arm" "$bin" "$leg" --gtf --bam $C20/chr20.bam --fasta $C20/chr20.fa --region chr20:1-66210255 ;;
    *) echo "unknown substrate: $sub" >&2; exit 2 ;;
  esac
}

SUBSTRATES="o2_families denovo chr20_gtf"
ARMS="base fixed legacy"

if [ $# -eq 2 ]; then
  run_one "$1" "$2"
elif [ $# -eq 1 ]; then
  for arm in $ARMS; do run_one "$1" "$arm"; done
elif [ $# -eq 0 ]; then
  for sub in $SUBSTRATES; do
    for arm in $ARMS; do run_one "$sub" "$arm"; done
  done
else
  echo "usage: $0 [SUBSTRATE [ARM]]" >&2
  exit 2
fi
