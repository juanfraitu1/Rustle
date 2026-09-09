#!/bin/bash
# flair arm of the tool bakeoff (docs/PREREG_tool_bakeoff_2026-09-08.md).
# FOREGROUND, one heavy run at a time -- flair align is a genome-wide minimap2 pass.
# flair at its DOCUMENTED DEFAULTS; the point of the arm is what the defaults do.
set -euo pipefail
W=/mnt/linuxdisk/home/juanfraitu/bakeoff
BAM=/mnt/linuxdisk/home/juanfraitu/npip_cat/npip3.bam
FA=/mnt/linuxdisk/home/juanfraitu/npip_cat/npip3_contigs.fa
COP=/mnt/linuxdisk/home/juanfraitu/mcl_ann/sweep_v19/fam_MCL1_073242/copies.tsv
source /home/juanfra/miniforge3/etc/profile.d/conda.sh; conda activate flair
mkdir -p "$W/flair"; cd "$W/flair"

# 1. reads. flair aligns from FASTQ, so the molecules must leave the BAM.
#    -F 2308 keeps primary/mapped/non-supplementary only, so both tools start from the SAME molecules.
if [ ! -s reads.fq ]; then
  samtools fastq -F 2308 -@ 4 "$BAM" > reads.fq 2> fastq.log
fi
echo "reads: $(( $(wc -l < reads.fq) / 4 ))"

# 2. flair align  -- NOTE minimap2 runs with --secondary=no inside flair (flair_align.py:150).
#    Nothing here forces that; it is flair's own default and is the finding under test.
[ -s flair.align.bed ] || flair align -g "$FA" -r reads.fq -o flair --threads 4 2>&1 | tail -5

# 3. flair correct -- no annotation given, so this is the unguided form (no --gtf/--shortread).
[ -s flair_all_corrected.bed ] || flair correct -q flair.bed -g "$FA" -o flair --threads 4 2>&1 | tail -5

# 4. flair collapse -- emits the isoform GTF that the common scorer reads.
[ -s flair.isoforms.gtf ] || flair collapse -g "$FA" -q flair_all_corrected.bed -r reads.fq \
     -o flair --threads 4 --generate_map 2>&1 | tail -5

ls -la "$W/flair"
