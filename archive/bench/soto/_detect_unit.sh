#!/bin/bash
# Detect one unit (a per-chrom BED, or the cross-chrom BED) with the given flags. Extract its reads from
# the cached subset BAM (once, cached), then run gw_family_catalog. Called in parallel by recompute_perchrom.sh.
BED=$1; FLAGS=$2
CACHE=${SOTO_CACHE:-/home/juanfra/winloci_scratch/soto_cache}; PC=$CACHE/perchrom
# The repo-local target/ DOES NOT EXIST, so this default used to resolve to nothing: every unit died
# with rc=127 while the caller carried on and scored stale outputs (o1_ledger.md §6am). The mandated
# build dir is CARGO_TARGET_DIR=/mnt/linuxdisk/home/juanfraitu/rustle_target; overridable via GWCAT.
BIN=${GWCAT:-/mnt/linuxdisk/home/juanfraitu/rustle_target/release/gw_family_catalog}
# Aborts before the rm -f of this unit's previous outputs, so a missing binary cannot destroy them. (§6am)
[ -x "$BIN" ] || { echo "[$(basename "$BED" .bed)] FATAL: gw_family_catalog not executable at $BIN (set GWCAT)" >&2; exit 2; }
# The big-data disk is a manual WSL mount and is often absent; fall back to the Desktop Reference copy
# rather than failing every unit with "failed to open FASTA".
FA=${SOTO_FASTA:-}
if [ -z "$FA" ]; then
  for cand in /mnt/linuxdisk/home/juanfraitu/winloci_data/chm13v2.0.fa \
              /mnt/c/Users/jfris/Desktop/Reference/chm13v2.0.fa; do
    [ -f "$cand" ] && FA="$cand" && break
  done
fi
[ -f "$FA" ] || { echo "[$(basename "$1" .bed)] FATAL: no CHM13 FASTA found (set SOTO_FASTA)" >&2; exit 2; }
SAM=${SAMTOOLS:-/home/juanfra/miniforge3/bin/samtools}
NAME=$(basename "$BED" .bed)
if [ ! -f "$PC/${NAME}.bam.bai" ]; then
  "$SAM" view -b -M -L "$BED" "$CACHE/soto_regions.bam" -o "$PC/${NAME}.bam" 2>/dev/null
  "$SAM" index "$PC/${NAME}.bam" 2>/dev/null
fi
# Stale outputs from a previous run must not survive a failure: the combine step downstream globs
# *.copies.tsv and cannot tell a fresh result from a leftover one. A whole recompute once "succeeded" in
# 17s -- every unit had died instantly on a missing FASTA and the old files were silently re-combined.
rm -f "$PC/${NAME}.copies.tsv" "$PC/${NAME}.families.tsv" "$PC/${NAME}.copies.fa"
timeout "${SOTO_TIMEOUT:-5400}" "$BIN" --bam "$PC/${NAME}.bam" --fasta "$FA" $FLAGS --out "$PC/${NAME}" > "$PC/${NAME}.log" 2>&1
RC=$?   # capture BEFORE any other command overwrites $?
if [ $RC -ne 0 ]; then
  echo "[$NAME] FAILED rc=$RC -- $(tail -1 "$PC/${NAME}.log" 2>/dev/null | cut -c1-100)" >&2
  exit $RC
fi
echo "[$NAME] copies=$(($(wc -l < "$PC/${NAME}.copies.tsv" 2>/dev/null || echo 1)-1)) rc=$RC"
