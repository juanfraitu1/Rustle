#!/bin/bash
#SBATCH -p singlenode -c 72 -t 04:00:00 -J r1694cat
#SBATCH -o /anvme/workspace/v138da13-rustle/r1694_exgain/log/catalog.%j.log
export PATH=$HOME/.pixi/detached-envs/rustle-7434162673922408322/envs/default/bin:$PATH
D=$HOME/work/jfim/rustle_data; W=/anvme/workspace/v138da13-rustle/r1694_exgain
# 3 h stop per the task: a catalog that runs longer is reported, not replaced
timeout 3h bash $HOME/Rustle/tools/rustle_pipeline.sh catalog --bam $D/reads/GCA_029281585.2_flnc_mm.bam --fasta $D/genomes/GGO.fa \
  --out $W/cat/kb --threads 72 --bin /anvme/workspace/v138da13-rustle/audit/g2/r76/bin
echo "catalog exit $? after ${SECONDS}s"
