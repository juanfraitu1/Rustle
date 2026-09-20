#!/bin/bash
# One arm of the annotation-degradation ablation (PREREG_annotation_ablation_2026-09-07, md5 65efc5b2).
# usage: o1_ablation_run.sh <ARM>   (A0 A0r A1 A2 A3 A4 A5)  — run ONE arm per call, foreground.
set -u
ARM=$1
D=/mnt/linuxdisk/home/juanfraitu; W=$D/mcl_ann/adj/ablation; mkdir -p $W
GFF=$D/winloci_data/GGO_genomic.gff; FA=$D/npip_cat/npip3_contigs.fa; BAM=$D/npip_cat/npip3.bam
SEDEF=$D/winloci_data/GGO_sedef_final.bed; RMSK=$D/winloci_data/rmsk/substrate3.rm.out
MF=$D/rustle_target/release/mcl_families; C=NC_073241.2,NC_073242.2,NC_073244.2
BENCH=/mnt/c/Users/jfris/Desktop/Rustle/bench
t0=$(date +%s)
if [ "$ARM" = "A0" ]; then AG=$GFF; else
  AG=$W/$ARM.gff; python3 $BENCH/o1_degrade_gff.py $GFF ${ARM%r} $AG $C || exit 1
fi
SP=$W/$ARM.regions
awk -F'\t' -v c=$C 'BEGIN{split(c,a,",");for(i in a)k[a[i]]=1} ($1 in k)&&($3=="gene"||$3=="pseudogene"){print $1":"$4"-"$5}' $AG | sort -u > $SP
echo "[$ARM] spans: $(wc -l < $SP)"
PAF=$W/$ARM.paf
if [ "$ARM" = "A0r" ] && [ -s $W/A0.paf ]; then PAF=$W/A0.paf; echo "[$ARM] reusing A0.paf"
else
  samtools faidx $FA -r $SP > $W/$ARM.fa || exit 1
  echo "[$ARM] fasta $(grep -c '>' $W/$ARM.fa) seqs, $(du -h $W/$ARM.fa | cut -f1); minimap2 ..."
  minimap2 -x asm20 -c -X -N 50 -p 0.1 -t 4 $W/$ARM.fa $W/$ARM.fa > $PAF 2> $W/$ARM.mm2.log || exit 1
  echo "[$ARM] paf $(wc -l < $PAF) records, $(( $(date +%s) - t0 )) s"
fi
$MF --paf $PAF --gff $AG --min-exonic-bp 1 --merge-overlapping-loci --core-refine --emit-units \
    --sedef $SEDEF --bam $BAM --fasta $FA --rmsk $RMSK --out $W/$ARM > $W/$ARM.log 2>&1
echo "[$ARM] mcl_families exit=$? units=$(($(wc -l < $W/$ARM.units.tsv 2>/dev/null || echo 1)-1)) total $(( $(date +%s) - t0 )) s"
