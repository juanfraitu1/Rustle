#!/bin/bash
# FAST + CORRECT parallel Soto-validation recompute. gw_family_catalog's default (non-cross-chrom) mode
# processes each CHROMOSOME independently, so per-chrom parallel jobs reproduce the genome-wide within-chrom
# detection EXACTLY (unlike per-family, which isolates co-located families and over-detects). The 18
# cross-chrom families are handled by ONE overlapping --cross-chrom pass on just their regions. Cache:
# soto_regions.bam + perchrom/beds/*.bed + perchrom/xchrom.bed. Per-unit mini-BAMs are extracted once + cached.
#
# Usage: bench/soto/recompute_perchrom.sh [nproc]   (default 5 for the per-chrom pool; +1 for xchrom)
set -u
CACHE=${SOTO_CACHE:-/home/juanfra/winloci_scratch/soto_cache}; PC=$CACHE/perchrom
HERE="$(cd "$(dirname "$0")" && pwd)"
export SOTO_CACHE
NP="${1:-5}"
[ -d "$PC/beds" ] || { echo "per-chrom BEDs missing"; exit 1; }
NCH=$(ls "$PC/beds"/*.bed | wc -l)
echo "per-chrom recompute: $NCH chroms (P=$NP) + 1 cross-chrom pass  (binary $(date -r "${GWCAT:-/mnt/c/Users/jfris/Desktop/Rustle/target/release/gw_family_catalog}" '+%m-%d %H:%M' 2>/dev/null))"
t0=$(date +%s 2>/dev/null || echo 0)
# cross-chrom pass in background (overlaps the per-chrom pool). Its BAM is the 18 cross-chrom families'
# regions (~1.9 GB) and --cross-chrom is heavy, so give it a long timeout (the default 1800s truncated it).
SOTO_TIMEOUT="${XCHROM_TIMEOUT:-7200}" bash "$HERE/_detect_unit.sh" "$PC/xchrom.bed" "--cross-chrom --refine" &
XPID=$!
# per-chrom pool (NO --cross-chrom = chrom-independent = genome-wide within-chrom behavior)
ls "$PC/beds"/*.bed | xargs -P "$NP" -I{} bash "$HERE/_detect_unit.sh" {} "--refine"
wait $XPID
t1=$(date +%s 2>/dev/null || echo 0)
echo "=== all units done in $((t1-t0))s ==="
# combine: per-chrom copies + cross-chrom copies (union)
{ echo -e "family_id\tcopy_idx\ttid\tchrom\tstart\tend\tn_exon\tstrand\tn_reads"
  for f in "$PC"/chr*.copies.tsv "$PC"/xchrom.copies.tsv; do [ -f "$f" ] && tail -n +2 "$f"; done; } > "$CACHE/perchrom_catalog.copies.tsv"
echo "combined catalog: $(($(wc -l < "$CACHE/perchrom_catalog.copies.tsv")-1)) copies"
python3 "$HERE/soto_cache_score.py" "$CACHE/perchrom_catalog.copies.tsv"
