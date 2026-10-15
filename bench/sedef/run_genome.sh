#!/bin/bash
# One whole-genome `sd` stage into the persistent cache, under a 23 GiB address-space cap, timed
# (plan docs/superpowers/plans/2026-10-13-o1-layers-phase-1-2.md, Task 8). An allocation over the cap aborts sd-scan
# instead of waking the kernel's OOM killer.
#   run_genome.sh NAME FASTA   ->  /mnt/linuxdisk/tmp/o1_layers/NAME.sd.{prefix,log,time}
set -euo pipefail
NAME=$1; FA=$2
O=/mnt/linuxdisk/tmp/o1_layers; mkdir -p "$O"
export RUSTLE_SD_CACHE_DIR=/mnt/linuxdisk/home/juanfraitu/sd_cache
ulimit -v 24117248
cd "$(dirname "$0")/../.."
/usr/bin/time -f $'wall_s\t%e\npeak_rss_kb\t%M' -o "$O/$NAME.sd.time" \
  bash tools/rustle_pipeline.sh sd --fasta "$FA" --out "$O/$NAME" --threads 4 --bin /mnt/linuxdisk/home/juanfraitu/rustle_target/release
cat "$O/$NAME.sd.time"
