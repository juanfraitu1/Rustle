#!/bin/bash
# Detect one unit (a per-chrom BED, or the cross-chrom BED) with the given flags. Extract its reads from
# the cached subset BAM (once, cached), then run gw_family_catalog. Called in parallel by recompute_perchrom.sh.
BED=$1; FLAGS=$2
CACHE=${SOTO_CACHE:-/home/juanfra/winloci_scratch/soto_cache}; PC=$CACHE/perchrom
BIN=${GWCAT:-/mnt/c/Users/jfris/Desktop/Rustle/target/release/gw_family_catalog}
FA=${SOTO_FASTA:-/mnt/linuxdisk/home/juanfraitu/winloci_data/chm13v2.0.fa}
SAM=${SAMTOOLS:-/home/juanfra/miniforge3/bin/samtools}
NAME=$(basename "$BED" .bed)
if [ ! -f "$PC/${NAME}.bam.bai" ]; then
  "$SAM" view -b -M -L "$BED" "$CACHE/soto_regions.bam" -o "$PC/${NAME}.bam" 2>/dev/null
  "$SAM" index "$PC/${NAME}.bam" 2>/dev/null
fi
timeout "${SOTO_TIMEOUT:-1800}" "$BIN" --bam "$PC/${NAME}.bam" --fasta "$FA" $FLAGS --out "$PC/${NAME}" > "$PC/${NAME}.log" 2>&1
echo "[$NAME] copies=$(($(wc -l < "$PC/${NAME}.copies.tsv" 2>/dev/null || echo 1)-1)) rc=$?"
