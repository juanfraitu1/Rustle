#!/bin/bash
# map_reads.sh <mat|pat> <reads.fa> <out.bam>: the fibroblast BAM's own minimap2 command (@PG of GCA_029281585.2_flnc_mm.bam)
# against a haplotype splice index. Run under tools/rlock.sh heavy (loads a 13 GB index, ~15 GB RSS).
set -euo pipefail
hap=${1:?mat|pat}; fa=${2:?reads.fa}; out=${3:?out.bam}
IDX=/mnt/linuxdisk/home/juanfraitu/winloci_data/mGorGor1.$hap.splice.mmi
[ -s "$IDX" ] || { echo "missing $IDX" >&2; exit 2; }
minimap2 -ax splice:hq -uf --eqx -Y -N 50 -p 0.1 --secondary=yes -K 100M -t 4 "$IDX" "$fa" 2> "$out.log" \
  | samtools sort -@ 1 -m 500M -o "$out" -
samtools index "$out"
echo "$out: $(samtools view -c "$out") records, $(samtools view -c -F 2308 "$out") primaries, $(samtools view -c -f 4 "$out") unmapped"
