#!/bin/bash
# FASTEST Soto recompute: per-CO-LOCATED-CLUSTER parallel. Splits the genome-wide catalog into 34 co-located
# 5 Mb clusters + 18 small per-cross-chrom-family units, all run in parallel. Removes the single-threaded
# big-chrom / 1.9 GB cross-chrom TAIL that bounds recompute_perchrom (biggest unit now ~3.3 Mb, not a whole
# chromosome). Correctness = per-chrom (co-located conflict resolution preserved within each 5 Mb cluster);
# the cross-chrom per-family units isolate those 18 families (mild over-detection on that subset only, same as
# the xchrom pass). Cache: soto_regions.bam + percluster/beds + manifest.json (build with the cluster builder).
#
# Usage: bench/soto/recompute_percluster.sh [nproc]   (default = all cores)
#   fast dev tier: SOTO_UNIT_GLOB='cl0*' bench/soto/recompute_percluster.sh   (run a subset of clusters)
set -u
HERE="$(cd "$(dirname "$0")" && pwd)"
CACHE=${SOTO_CACHE:-/home/juanfra/winloci_scratch/soto_cache}; PCL=$CACHE/percluster
export SOTO_CACHE
NP="${1:-$(nproc)}"
[ -f "$PCL/manifest.json" ] || { echo "cluster manifest missing — build the cluster cache first"; exit 1; }
NU=$(python3 -c "import json;print(len(json.load(open('$PCL/manifest.json'))))")
echo "per-cluster recompute: $NU units, P=$NP  (binary $(date -r "${GWCAT:-/mnt/c/Users/jfris/Desktop/Rustle/target/release/gw_family_catalog}" '+%m-%d %H:%M' 2>/dev/null))"
t0=$(date +%s 2>/dev/null || echo 0)
# task list = "name|flags" per unit; run in parallel (newline-delimited so flags with spaces survive)
python3 -c "import json;[print(u['name']+'|'+u['flags']) for u in json.load(open('$PCL/manifest.json'))]" \
  | xargs -P "$NP" -d '\n' -I{} bash "$HERE/_cluster_one.sh" "{}"
t1=$(date +%s 2>/dev/null || echo 0)
echo "=== all $NU clusters done in $((t1-t0))s ==="
{ echo -e "family_id\tcopy_idx\ttid\tchrom\tstart\tend\tn_exon\tstrand\tn_reads"
  for f in "$PCL"/cl*.copies.tsv "$PCL"/xf_*.copies.tsv; do [ -f "$f" ] && tail -n +2 "$f"; done; } > "$CACHE/percluster_catalog.copies.tsv"
echo "combined catalog: $(($(wc -l < "$CACHE/percluster_catalog.copies.tsv")-1)) copies"
python3 "$HERE/soto_cache_score.py" "$CACHE/percluster_catalog.copies.tsv"
