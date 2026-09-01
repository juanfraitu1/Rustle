#!/bin/bash
# FAST + CORRECT parallel Soto-validation recompute.
#
# Each per-chrom unit runs `--cross-chrom`, because THAT IS THE ALGORITHM THE HEADLINE CATALOG USED. The
# `definitive` catalog (863 copies / 280 families / the 76.2% member figure) is ONE `--cross-chrom` pass over
# the whole soto_regions.bam -- Louvain community split on a global conflict graph
# (`detect_conflict_catalog_genome_wide_xchrom`). gw_family_catalog's DEFAULT mode is a different algorithm:
# per-chrom connected components emitted as "clean families (same-strand, disjoint-loci, >= 2 copies)"
# (`gw_reps_and_catalog`). It is strictly more conservative, and running it here is what made this recipe
# non-comparable to the headline.
#
# Earlier versions of this header asserted that the default mode "reproduces the genome-wide within-chrom
# detection EXACTLY". That was false and it silently halved the benchmark: the recipe scored 95/362 = 26.2%
# against a 76.2% headline. Measured per chromosome, default vs --cross-chrom vs definitive's own breakdown:
#   chr1   45 -> 180 (definitive 163)      chr15  18 -> 132 (123)      chr9    6 -> 98 (97)
# Cross-chrom families still need the separate overlapping pass, since a per-chromosome split cannot form an
# edge between loci on different chromosomes.
#
# Cache: soto_regions.bam + perchrom/beds/*.bed + perchrom/xchrom.bed. Per-unit mini-BAMs extracted once.
#
# Usage: bench/soto/recompute_perchrom.sh [nproc]   (default 5 for the per-chrom pool; +1 for xchrom)
set -u
CACHE=${SOTO_CACHE:-/home/juanfra/winloci_scratch/soto_cache}; PC=$CACHE/perchrom
HERE="$(cd "$(dirname "$0")" && pwd)"
export SOTO_CACHE
NP="${1:-3}"   # --cross-chrom units are heavier than the old default-mode ones (WSL2 memory)
[ -d "$PC/beds" ] || { echo "per-chrom BEDs missing"; exit 1; }
NCH=$(ls "$PC/beds"/*.bed | wc -l)
echo "per-chrom recompute: $NCH chroms (P=$NP) + 1 cross-chrom pass  (binary $(date -r "${GWCAT:-/mnt/linuxdisk/home/juanfraitu/rustle_target/release/gw_family_catalog}" '+%m-%d %H:%M' 2>/dev/null))"
t0=$(date +%s 2>/dev/null || echo 0)
# cross-chrom pass in background (overlaps the per-chrom pool). Its BAM is the 18 cross-chrom families'
# regions (~1.9 GB) and --cross-chrom is heavy, so give it a long timeout (the default 1800s truncated it).
SOTO_TIMEOUT="${XCHROM_TIMEOUT:-7200}" bash "$HERE/_detect_unit.sh" "$PC/xchrom.bed" "--cross-chrom --refine" &
XPID=$!
# per-chrom pool. --cross-chrom is the headline algorithm (see header); within a single-chromosome BAM it
# can only form same-chrom families, so this stays chrom-independent and parallel-safe.
ls "$PC/beds"/*.bed | xargs -P "$NP" -I{} bash "$HERE/_detect_unit.sh" {} "--cross-chrom --refine"
wait $XPID
t1=$(date +%s 2>/dev/null || echo 0)
echo "=== all units done in $((t1-t0))s ==="
# combine: per-chrom copies + cross-chrom copies (union)
{ echo -e "family_id\tcopy_idx\ttid\tchrom\tstart\tend\tn_exon\tstrand\tn_reads"
  for f in "$PC"/chr*.copies.tsv "$PC"/xchrom.copies.tsv; do [ -f "$f" ] && tail -n +2 "$f"; done; } > "$CACHE/perchrom_catalog.copies.tsv"
echo "combined catalog: $(($(wc -l < "$CACHE/perchrom_catalog.copies.tsv")-1)) copies"
python3 "$HERE/soto_cache_score.py" "$CACHE/perchrom_catalog.copies.tsv"
