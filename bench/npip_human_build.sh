#!/bin/bash
# Human NPIP substrate (PREREG_npip_human_2026-09-07, md5 dfbe88b6): chr16+chr18 gene spans, CHM13.
set -u
D=/mnt/linuxdisk/home/juanfraitu; W=$D/soto_mcl/npip_hsa; mkdir -p $W
GFF=$D/soto_mcl/hsa.gff; FA=$D/winloci_data/chm13v2.0.fa; BAM=$D/soto_adj/soto.bam
MF=$D/rustle_target/release/mcl_families
t0=$(date +%s)
awk -F'\t' '($1=="chr16"||$1=="chr18") && ($3=="gene"||$3=="pseudogene"){print $1":"$4"-"$5}' $GFF | sort -u > $W/spans.regions
echo "spans: $(wc -l < $W/spans.regions)"
samtools faidx $FA -r $W/spans.regions > $W/spans.fa || exit 1
echo "fasta: $(grep -c '>' $W/spans.fa) seqs $(du -h $W/spans.fa | cut -f1); minimap2 ..."
minimap2 -x asm20 -c -X -N 50 -p 0.1 -t 4 $W/spans.fa $W/spans.fa > $W/hsa.paf 2> $W/mm2.log || exit 1
echo "paf: $(wc -l < $W/hsa.paf) records, $(( $(date +%s) - t0 )) s"
$MF --paf $W/hsa.paf --gff $GFF --min-exonic-bp 1 --merge-overlapping-loci --core-refine --core-from-paf \
    --emit-units --bam $BAM --fasta $FA --out $W/cat > $W/cat.log 2>&1
echo "mcl_families exit=$? units=$(($(wc -l < $W/cat.units.tsv 2>/dev/null || echo 1)-1)) total $(( $(date +%s) - t0 )) s"
