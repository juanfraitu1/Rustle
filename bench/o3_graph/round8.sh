#!/bin/bash
# round8.sh OUTROOT FAMILIES SET SEED MASKDIR SCENARIO  (docs/archive/2026-10/merged/PREREG_o3_graph_round8_2026-10-16.md)
# One world, end to end: build (alleles, --het 0.001) + map + the current method (copy_assign O2 + its --psv-clusters O3) +
# the integrated O1 -> O2 -> O3 (o1seed.py, o2graph.py run W integ 2, detect.py arm I). REMOVED maps to MASKDIR's masked
# genome and index (realworld.py mask, then minimap2 -x splice:hq -d); CROSS projects by liftover.
set -u
cd "$(dirname "$0")/../.."
OUT=$1; FJ=$2; SET=$3; SEED=$4; MASK=$5; sc=$6; W=$OUT/$sc; L=$OUT/$sc.log
IDX=/mnt/linuxdisk/home/juanfraitu/npip_ladder/idx/target.splice.mmi; G=""; X=""
if [ "$sc" = REMOVED ]; then G="--genome $MASK/masked.fa"; IDX=$MASK/masked.splice.mmi; fi
if [ "$sc" = CROSS ]; then X="--projection liftover"; fi
mkdir -p "$OUT"
step() { local t=$1; shift; RLOCK_TIMEOUT=$t bash tools/rlock.sh heavy "$@" >> "$L" 2>&1; }
: > "$L"
step 200 python3 bench/o3_graph/realworld.py build $W $FJ --set $SET --scenario $sc --seed $SEED --het 0.001 $G $X || { echo "FAIL build $W"; exit 1; }
step 400 python3 bench/o3_graph/run.py map $W --index $IDX || { echo "FAIL map $W"; exit 1; }
step 300 python3 bench/o3_graph/run.py cur $W /mnt/linuxdisk/tmp/o3_fb_build/copy_assign.fb || { echo "FAIL cur $W"; exit 1; }
step 200 python3 bench/o3_graph/o1seed.py $W || { echo "FAIL o1seed $W"; exit 1; }
step 300 python3 bench/o3_graph/o2graph.py run $W integ 2 || { echo "FAIL integ $W"; exit 1; }
step 500 python3 bench/o3_graph/detect.py $W || { echo "FAIL detect $W"; exit 1; }
echo "$W: $(grep -o 'o1seed.*' $L | cut -c1-120) | $(grep 'o3_graph\]' $L | grep -o 'I [0-9]* calls')"
