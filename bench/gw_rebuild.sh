#!/bin/bash
# Genome-wide gorilla catalog under the 2026-09-05 defaults (fold within clusters, exon-less span, exon-to-exon,
# units follow reads, min_size 2), with SEDEF cores, the curated repeat library, and the genome-wide downsampled
# IsoSeq BAM for units. Output: mcl_ann/gw_units_v1.*  (clusters, loci, cores, blocks, units).
set -u
M=/mnt/linuxdisk/home/juanfraitu/mcl_ann; W=/mnt/linuxdisk/home/juanfraitu/winloci_data
MF=/mnt/linuxdisk/home/juanfraitu/rustle_target/release/mcl_families
cd $M && ( ulimit -v 20000000; /usr/bin/time -f "TIME %e s rss %M kB" $MF --paf allgenes_gw.asm20.paf --gff $W/GGO_genomic.gff --min-exonic-bp 1 \
  --sedef $W/GGO_sedef_final.bed --rmsk $W/rmsk/GCF_029281585.2.repeatMasker.out \
  --bam $W/GGO_ds.bam --fasta /home/juanfra/winloci_scratch/GGO.fasta --out $M/gw_units_v1 ) > $M/adj/gw_units_v1.log 2>&1
echo "gw exit=$?"; tail -4 $M/adj/gw_units_v1.log
