#!/bin/bash
# Build the CACHED Soto-region subset BAM. ONE-TIME — rebuild only when the alignment (A119b.t2t.bam)
# changes. All Soto-validation recomputes then run on this ~5 GB subset via recompute_soto.sh (minutes),
# instead of the ~96 GB full BAM genome-wide (hours). The advisor accepts only Soto, so we recompute often.
#
# Design: extract every read overlapping a Soto family member locus (+/- PAD) via the BAM index (-L), so
# the fetch touches only the ~33 Mb of member neighborhoods, not the whole genome. Co-located siblings and
# cross-chrom partners are all included (every member is in the BED), so family detection is preserved.
set -euo pipefail
BAM=${SOTO_BAM:-/mnt/linuxdisk/home/juanfraitu/winloci_data/A119b.t2t.bam}
BED_SRC="$(cd "$(dirname "$0")" && pwd)/80_fams.chr.bed"
CACHE=${SOTO_CACHE:-/home/juanfra/winloci_scratch/soto_cache}
SAM=${SAMTOOLS:-/home/juanfra/miniforge3/bin/samtools}
PAD=${PAD:-50000}
mkdir -p "$CACHE"

# merged region BED = member intervals +/- PAD, sorted + merged
awk -F'\t' -v p="$PAD" '{s=$2-p; if(s<0)s=0; print $1"\t"s"\t"$3+p}' "$BED_SRC" \
  | sort -k1,1 -k2,2n \
  | awk 'BEGIN{OFS="\t"} {if($1==c && $2<=e){if($3>e)e=$3} else {if(c!="")print c,s,e; c=$1;s=$2;e=$3}} END{if(c!="")print c,s,e}' \
  > "$CACHE/soto_regions.bed"
echo "merged regions: $(wc -l < "$CACHE/soto_regions.bed")  span: $(awk '{s+=$3-$2}END{printf "%.1f Mb\n",s/1e6}' "$CACHE/soto_regions.bed")"

# indexed fetch of reads overlapping the regions (multi-region iterator)
"$SAM" view -b -M -L "$CACHE/soto_regions.bed" -@ 8 "$BAM" -o "$CACHE/soto_regions.bam"
"$SAM" index "$CACHE/soto_regions.bam"
echo "cache built: $CACHE/soto_regions.bam  ($("$SAM" view -c "$CACHE/soto_regions.bam") reads, $(du -h "$CACHE/soto_regions.bam" | cut -f1))"
echo "now recompute with: bench/soto/recompute_soto.sh"
