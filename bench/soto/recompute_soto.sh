#!/bin/bash
# FAST cached Soto-validation recompute. Run after ANY binary/pipeline change to re-measure member-detection
# recall/precision on the cached subset BAM (built once by build_soto_cache.sh). Minutes, not hours.
#
# Usage:  bench/soto/recompute_soto.sh ["<gw_family_catalog flags>"]
#   default flags: --cross-chrom --refine   (the genome-wide Soto catalog recipe)
# Writes <cache>/catalog.{copies,families}.tsv and prints recall vs the committed 76.2% baseline.
set -euo pipefail
CACHE=${SOTO_CACHE:-/home/juanfra/winloci_scratch/soto_cache}
BIN=${GWCAT:-/mnt/c/Users/jfris/Desktop/Rustle/target/release/gw_family_catalog}
FA=${SOTO_FASTA:-/mnt/linuxdisk/home/juanfraitu/winloci_data/chm13v2.0.fa}
FLAGS=${1:---cross-chrom --refine}
HERE="$(cd "$(dirname "$0")" && pwd)"
[ -f "$CACHE/soto_regions.bam" ] || { echo "cache missing — run bench/soto/build_soto_cache.sh first"; exit 1; }

echo "recompute: gw_family_catalog $FLAGS  (binary $(date -r "$BIN" '+%Y-%m-%d %H:%M'))"
stdbuf -oL -eL /usr/bin/time -v "$BIN" --bam "$CACHE/soto_regions.bam" --fasta "$FA" $FLAGS \
   --threads 8 --out "$CACHE/catalog" > "$CACHE/catalog.log" 2>&1
echo "  catalog: $(($(wc -l < "$CACHE/catalog.copies.tsv")-1)) copies, $(($(wc -l < "$CACHE/catalog.families.tsv")-1)) families  ($(grep -oE 'wall clock.*' "$CACHE/catalog.log" | tail -1))"
python3 "$HERE/soto_cache_score.py" "$CACHE/catalog.copies.tsv"
