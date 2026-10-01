#!/bin/bash
# SQANTI3 QC + rules filter against the CAT reference GTF; same arms, inputs and flags as the registered RefSeq runs (sq3/ -> sq3_cat/)
set -u
B=/mnt/linuxdisk/home/juanfraitu/bakeoff; S=/mnt/linuxdisk/home/juanfraitu/_from_wsl/tools/SQANTI3
PY=/home/juanfra/miniforge3/envs/sqanti3/bin/python; export PATH=/home/juanfra/miniforge3/envs/sqanti3/bin:$PATH
c=$1
declare -A IN=([SHIP]=ship/SHIP.gtf [stringtie]=stringtie/st.gtf [flair]=flair/flair.isoforms.gtf)
D=$B/human_chr$c; O=$D/sq3_cat; mkdir -p $O
for a in SHIP stringtie flair; do
  [ -s $O/filt_$a/${a}_RulesFilter_result_classification.txt ] && { echo "chr$c $a done"; continue; }
  $PY $S/sqanti3_qc.py --isoforms $D/${IN[$a]} --refGTF $D/chr${c}_ref_cat.gtf --refFasta $D/chr$c.fa --report skip -t 4 -d $O/$a -o $a > $O/$a.log 2>&1 || { echo "chr$c $a QC FAILED"; exit 1; }
  $PY $S/sqanti3_filter.py rules --sqanti_class $O/$a/${a}_classification.txt --filter_gtf $O/$a/${a}_corrected.gtf --skip_report -d $O/filt_$a -o $a > $O/filt_$a.log 2>&1 || { echo "chr$c $a FILTER FAILED"; exit 1; }
  echo "chr$c $a ok"
done
