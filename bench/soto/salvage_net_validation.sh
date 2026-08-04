#!/bin/bash
# NET verdict for mis-chain salvage: full per-chrom Soto recompute OFF vs ON (same binary), scored for member
# recall AND copy count (FP proxy). The isolated CDH12 test showed a regression (conflict-graph over-connection);
# this measures whether the net genome-wide effect is positive, neutral, or negative. Controller-run (crash rule).
set -u
HERE=/mnt/c/Users/jfris/Desktop/Rustle/bench/soto
CACHE=/home/juanfra/winloci_scratch/soto_cache
ncopies(){ echo $(($(wc -l < "$1" 2>/dev/null || echo 1)-1)); }
recall(){ python3 "$HERE/soto_cache_score.py" "$1" 2>/dev/null | grep -oE "NEW: [0-9]+/362 = [0-9.]+%" | sed 's/NEW: //'; }

echo "[$(date '+%H:%M:%S')] === OFF (default, salvage off) full per-chrom recompute ==="
rm -f "$CACHE/perchrom/"chr*.copies.tsv "$CACHE/perchrom/xchrom.copies.tsv"
bash "$HERE/recompute_perchrom.sh" 5 > "$CACHE/salv_net_off.log" 2>&1
cp "$CACHE/perchrom_catalog.copies.tsv" "$CACHE/salv_off_catalog.copies.tsv"
OFFC=$(ncopies "$CACHE/salv_off_catalog.copies.tsv"); OFFR=$(recall "$CACHE/salv_off_catalog.copies.tsv")
echo "[$(date '+%H:%M:%S')] OFF: copies=$OFFC recall=$OFFR"

echo "[$(date '+%H:%M:%S')] === ON (RUSTLE_MISCHAIN_SALVAGE=1) full per-chrom recompute ==="
rm -f "$CACHE/perchrom/"chr*.copies.tsv "$CACHE/perchrom/xchrom.copies.tsv"
RUSTLE_MISCHAIN_SALVAGE=1 bash "$HERE/recompute_perchrom.sh" 5 > "$CACHE/salv_net_on.log" 2>&1
cp "$CACHE/perchrom_catalog.copies.tsv" "$CACHE/salv_on_catalog.copies.tsv"
ONC=$(ncopies "$CACHE/salv_on_catalog.copies.tsv"); ONR=$(recall "$CACHE/salv_on_catalog.copies.tsv")
echo "[$(date '+%H:%M:%S')] ON:  copies=$ONC recall=$ONR"

echo "=== NET VERDICT ==="
echo "OFF: copies=$OFFC recall=$OFFR"
echo "ON : copies=$ONC recall=$ONR"
echo "delta copies=$((ONC-OFFC))  (recall: OFF $OFFR -> ON $ONR)"
echo "SHIP iff ON recall > OFF recall AND copy inflation is real recoveries, not FP over-connection."
