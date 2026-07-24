#!/bin/bash
# Detect ONE Soto family: extract its reads from the cached subset BAM (once, cached) then run
# gw_family_catalog on that tiny per-family BAM. Called in parallel by recompute_perfam.sh.
FID=$1
CACHE=${SOTO_CACHE:-/home/juanfra/winloci_scratch/soto_cache}; PF=$CACHE/perfam
BIN=${GWCAT:-/mnt/c/Users/jfris/Desktop/Rustle/target/release/gw_family_catalog}
FA=${SOTO_FASTA:-/mnt/linuxdisk/home/juanfraitu/winloci_data/chm13v2.0.fa}
SAM=${SAMTOOLS:-/home/juanfra/miniforge3/bin/samtools}
FLAGS=${SOTO_FLAGS:---cross-chrom --refine}
if [ ! -f "$PF/${FID}.bam.bai" ]; then
  "$SAM" view -b -M -L "$PF/beds/${FID}.bed" "$CACHE/soto_regions.bam" -o "$PF/${FID}.bam" 2>/dev/null
  "$SAM" index "$PF/${FID}.bam" 2>/dev/null
fi
timeout 300 "$BIN" --bam "$PF/${FID}.bam" --fasta "$FA" $FLAGS --out "$PF/${FID}" > "$PF/${FID}.log" 2>&1
echo "[$FID] copies=$(($(wc -l < "$PF/${FID}.copies.tsv" 2>/dev/null || echo 1)-1)) rc=$?"
