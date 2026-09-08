#!/bin/bash
# Human chrY ampliconic substrate (PREREG_yags_daz_2026-09-07, md5 be38f3f2).
set -u
D=/mnt/linuxdisk/home/juanfraitu; W=$D/soto_mcl/yag; mkdir -p $W
GFF=$D/soto_mcl/hsa.gff; FA=$D/winloci_data/chm13v2.0.fa; BAM=$D/winloci_data/A119b.t2t.bam
SED=$D/winloci_data/HSA_sedef_pairs.bed; MF=$D/rustle_target/release/mcl_families
t0=$(date +%s)
# AMENDMENT 1: euchromatic MSY only (start < 28 Mb). The Yq12 satellite block carries several 0.8-1.1 Mb
# models whose all-vs-all is quadratic in anchors; no alignment parameter is changed, only the node set.
awk -F'\t' '$1=="chrY" && ($3=="gene"||$3=="pseudogene") && $4 < 28000000 {print $1":"$4"-"$5}' $GFF | sort -u > $W/spans.regions
echo "spans: $(wc -l < $W/spans.regions)"
samtools faidx $FA -r $W/spans.regions > $W/spans.fa || exit 1
minimap2 -x asm20 -c -X -N 50 -p 0.1 -t 4 $W/spans.fa $W/spans.fa > $W/y.paf 2> $W/mm2.log || exit 1
echo "paf: $(wc -l < $W/y.paf) records, $(( $(date +%s) - t0 )) s"
$MF --paf $W/y.paf --gff $GFF --min-exonic-bp 1 --merge-overlapping-loci --core-refine --sedef $SED \
    --emit-units --emit-readthrough-units --bam $BAM --fasta $FA --out $W/cat > $W/cat.log 2>&1
echo "mcl_families exit=$? units=$(($(wc -l < $W/cat.units.tsv 2>/dev/null || echo 1)-1)) total $(( $(date +%s) - t0 )) s"
grep -iE "readthrough|clusters" $W/cat.log | tail -3
