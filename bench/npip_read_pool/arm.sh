#!/bin/bash
# arm.sh ARM (P|GOOD|ALL): frozen copy_assign --assemble-only on chr16 (A119b) with the shipped polish; then loci, all-vs-all, MCL.
set -uo pipefail
A=$1; D=/mnt/linuxdisk/tmp/readpool_npip; export TMPDIR=$D/tmp; mkdir -p $TMPDIR
BIN=/mnt/linuxdisk/tmp/rustle_figures/cc_bin_frozen; M2=/mnt/linuxdisk/home/juanfraitu/rustle_target_m2/release
FASTA=/mnt/linuxdisk/home/juanfraitu/winloci_data/chm13v2.0.fa; BAM=/mnt/linuxdisk/home/juanfraitu/winloci_data/A119b.t2t.bam
TAB=/mnt/linuxdisk/tmp/rustle_figures/runs/human_A119b/human_A119b.molecules.tsv
POLISH="--assembly-polish full --polish-isoform-fraction 0.02 --polish-mono-shadow --polish-mono-quantile 0.82 --polish-ism-ratio 0.7 --polish-retained-ratio 10"
case $A in P) ENV=();; GOOD) ENV=(RUSTLE_GTF_SECONDARY=1 RUSTLE_GTF_SECONDARY_AS_RATIO=0.98 RUSTLE_GTF_SECONDARY_AS_TABLE=$TAB);;
  ALL) ENV=(RUSTLE_GTF_SECONDARY=1);; *) echo bad; exit 2;; esac
echo "chr16:0-96330374" > $D/chr16.txt
step=${2:-assemble}
if [ "$step" = assemble ]; then
  t0=$(date +%s)
  env "${ENV[@]}" /usr/bin/time -v $BIN/copy_assign --assemble-only --regions $D/chr16.txt --assembly-junctions strict $POLISH \
    --bam $BAM --fasta $FASTA --out $D/$A > $D/$A.log 2>&1; rc=$?
  echo "$A assemble rc=$rc $(( $(date +%s)-t0 )) s $(grep 'Maximum resident' $D/$A.log)"
else
  t0=$(date +%s)
  python3 /mnt/linuxdisk/tmp/gw22/sec/loci_from_gtf.py $D/$A.gtf $FASTA $D/$A > $D/$A.loci.log 2>&1
  /usr/bin/time -v minimap2 -x asm20 -c -X -N 50 -p 0.1 --secondary=yes -t 4 $D/$A.loci.fa $D/$A.loci.fa > $D/$A.paf 2> $D/$A.mm2.log
  $M2/mcl_families --paf $D/$A.paf --gff $D/$A.gff3 --min-exonic-bp 1 --min-shared-exon-frac 0.60 --out $D/$A.fam > $D/$A.mcl.log 2>&1
  echo "$A family rc=$? $(( $(date +%s)-t0 )) s loci=$(grep -c '^>' $D/$A.loci.fa) paf=$(wc -l < $D/$A.paf) $(grep 'Maximum resident' $D/$A.mm2.log)"
fi
