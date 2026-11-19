#!/bin/bash
#SBATCH --account=open
#SBATCH --partition=basic
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=2:00:00
#SBATCH --job-name=r1829_l1
# PREREG r1829 layer 1: near-tie tiers per sample from each BAM's as_table
source /storage/group/kdm16/default/jxi21/fritz/env.sh
cd $FRITZ/Rustle
W=$FRITZ/ws/r1829_near_tie_census; mkdir -p $W/layer1
R=$HOME/reads; C=$HOME/data/chm13
S="A119b:HSA/A119b.t2t.bam HG002:HSA/HG002_kinnex.t2t.bam OR6737:GGO/GGO_mm.bam KB3781:GGO/GGO_KB3781_fib_mm.bam
PTR_testis:PTR/PTR_mm.bam PTR_AG18354_lcl:PTR/PTR_AG18354_lcl_mm.bam PAB_testis:PAB/PAB_mm.bam PAB_AG06213_fib:PAB/PAB_AG06213_fib_mm.bam
PPA_testis:PPA/PPA_mm.bam PPA_PR00251_fib:PPA/PPA_PR00251_fib_mm.bam PPY_testis:PPY/PPY_mm.bam PPY_AG05252_fib:PPY/PPY_AG05252_fib_mm.bam"
for x in $S; do
  n=${x%%:*}; b=$R/${x#*:}; t=$b.as_table.tsv
  # the table must describe this BAM (header bam_bytes = the BAM's size)
  want=$(head -1 $t | tr '\t' '\n' | grep '^bam_bytes=' | cut -d= -f2); have=$(stat -Lc %s $b)
  echo "$n table_bytes=$want bam_bytes=$have $([ "$want" = "$have" ] && echo MATCH || echo MISMATCH)" >> $W/layer1/tables.check
  extra=""
  case $n in
    A119b) extra="--sd $C/chm13v2.0.fa.sedef.bed --exclude $C/chm13v2.0.fa.exclude.bed --genes x";;
    HG002) extra="--sd $C/chm13v2.0.fa.sedef.bed --exclude $C/chm13v2.0.fa.exclude.bed --genes x --family-names $HOME/ws/r1776_hg002_room/hg.sel.fa --hash-pct 25";;
  esac
  pixi run -q python bench/o2_truth/near_tie_census.py $n $t $extra > $W/layer1/$n.census.tsv 2> $W/layer1/$n.err &
done
wait
cat $W/layer1/tables.check
cat $W/layer1/*.census.tsv | awk 'NR==1 || $1!="sample"' > $W/census.tsv
wc -l $W/census.tsv; cat $W/layer1/*.err | head
