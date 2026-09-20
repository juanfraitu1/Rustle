#!/bin/bash
# FAST parallel per-family Soto-validation recompute. Each of the 83 families is a tiny mini-BAM;
# gw_family_catalog runs on them IN PARALLEL (wall-clock ~ slowest family, not the serial sum) — the fix
# for gw_family_catalog's single-threaded family-build phase. Cache: soto_regions.bam + per-family BEDs
# (build_soto_cache.sh + the per-family BED builder). Per-family mini-BAMs are extracted once and cached.
#
# Usage: bench/soto/recompute_perfam.sh ["<gw_family_catalog flags>"] [nproc]
#   default: "--cross-chrom --refine"  6
set -u
CACHE=${SOTO_CACHE:-/home/juanfra/winloci_scratch/soto_cache}; PF=$CACHE/perfam
HERE="$(cd "$(dirname "$0")" && pwd)"
export SOTO_CACHE SOTO_FLAGS="${1:---cross-chrom --refine}"
NP="${2:-6}"
[ -d "$PF/beds" ] || { echo "per-family BEDs missing — build them first"; exit 1; }
FAMS=$(ls "$PF/beds"/*.bed | sed 's#.*/##;s/\.bed//')
NF=$(echo "$FAMS" | wc -l)
echo "parallel per-family recompute: $NF families, flags=$SOTO_FLAGS, P=$NP  (binary $(date -r "${GWCAT:-/mnt/c/Users/jfris/Desktop/Rustle/target/release/gw_family_catalog}" '+%m-%d %H:%M'))"
t0=$(date +%s 2>/dev/null || echo 0)
echo "$FAMS" | xargs -P "$NP" -I{} bash "$HERE/_perfam_one.sh" {}
t1=$(date +%s 2>/dev/null || echo 0)
echo "=== all families done in $((t1-t0))s ==="
# combine per-family catalogs
{ echo -e "family_id\tcopy_idx\ttid\tchrom\tstart\tend\tn_exon\tstrand\tn_reads"
  for f in "$PF"/*.copies.tsv; do [ -f "$f" ] && tail -n +2 "$f"; done; } > "$CACHE/perfam_catalog.copies.tsv"
echo "combined catalog: $(($(wc -l < "$CACHE/perfam_catalog.copies.tsv")-1)) copies"
python3 "$HERE/soto_cache_score.py" "$CACHE/perfam_catalog.copies.tsv"
