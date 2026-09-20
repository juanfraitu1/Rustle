#!/bin/bash
# Detect one cluster unit "name|flags": extract its reads from the cached subset BAM (once, cached) then run
# gw_family_catalog with the unit's flags. Called in parallel by recompute_percluster.sh.
IFS='|' read -r NAME FLAGS <<< "$1"
CACHE=${SOTO_CACHE:-/home/juanfra/winloci_scratch/soto_cache}; PCL=$CACHE/percluster
BIN=${GWCAT:-/mnt/c/Users/jfris/Desktop/Rustle/target/release/gw_family_catalog}
FA=${SOTO_FASTA:-/mnt/linuxdisk/home/juanfraitu/winloci_data/chm13v2.0.fa}
SAM=${SAMTOOLS:-/home/juanfra/miniforge3/bin/samtools}
if [ ! -f "$PCL/${NAME}.bam.bai" ]; then
  "$SAM" view -b -M -L "$PCL/beds/${NAME}.bed" "$CACHE/soto_regions.bam" -o "$PCL/${NAME}.bam" 2>/dev/null
  "$SAM" index "$PCL/${NAME}.bam" 2>/dev/null
fi
timeout "${SOTO_TIMEOUT:-1200}" "$BIN" --bam "$PCL/${NAME}.bam" --fasta "$FA" $FLAGS --out "$PCL/${NAME}" > "$PCL/${NAME}.log" 2>&1
echo "[$NAME] copies=$(($(wc -l < "$PCL/${NAME}.copies.tsv" 2>/dev/null || echo 1)-1)) rc=$?"
