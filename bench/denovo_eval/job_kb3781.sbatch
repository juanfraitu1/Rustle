#!/bin/bash
#SBATCH --job-name=kb3781_fam
#SBATCH --partition=singlenode
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=72
#SBATCH --time=10:00:00
#SBATCH --output=/home/hpc/v138da/v138da13/runs/gorilla_KB3781/slurm-%j.out
# Gorilla KB3781 fibroblast Iso-Seq (GCA_029281585.2_flnc_mm.bam on GCF_029281585.2): assemble -> de novo families
# (baseline, exon-gated all-vs-all) -> the same families with --shared-exon-columns, replaying the cached PAF.
set -euo pipefail
export PATH=$HOME/.pixi/detached-envs/rustle-7434162673922408322/envs/default/bin:$PATH   # minimap2-gate (gate6), samtools
DATA=$HOME/work/jfim/rustle_data
BAM=$DATA/reads/GCA_029281585.2_flnc_mm.bam
FASTA=$DATA/genomes/GGO.fa
BIN=$HOME/rustle_target/release
RUN=$HOME/runs/gorilla_KB3781
T=${THREADS:-36}
cd $HOME/Rustle
echo "[job] $(date) host=$(hostname) commit=$(git rev-parse --short HEAD) minimap2-gate=$(minimap2-gate --version)"
STAGES=${STAGES:-assemble families columns}
for stage in $STAGES; do
  case $stage in
    assemble)
      echo "[job] $(date) assemble (t+${SECONDS}s)"
      bash tools/rustle_pipeline.sh assemble --bam $BAM --fasta $FASTA --out $RUN/kb --threads $T --bin $BIN
      ;;
    families)
      echo "[job] $(date) families (baseline) (t+${SECONDS}s)"
      bash tools/rustle_pipeline.sh families --bam $BAM --fasta $FASTA --out $RUN/kb --threads $T --bin $BIN
      ;;
    columns)
      echo "[job] $(date) families --shared-exon-columns (replaying $RUN/kb.cache) (t+${SECONDS}s)"
      RUSTLE_CACHE_DIR=$RUN/kb.cache $BIN/mcl_families --from-gtf $RUN/kb.families.gtf --fasta $FASTA --threads $T \
        --min-exonic-bp 1 --min-shared-exon-frac 0.60 --emit-units --shared-exon-columns --out $RUN/kbcols.fam \
        > $RUN/kbcols.families.log 2>&1
      ;;
    columns-union)
      echo "[job] $(date) families --shared-exon-columns --shared-exon-union (ungated all-vs-all) (t+${SECONDS}s)"
      RUSTLE_CACHE_DIR=$RUN/kb.cache $BIN/mcl_families --from-gtf $RUN/kb.families.gtf --fasta $FASTA --threads $T \
        --min-exonic-bp 1 --min-shared-exon-frac 0.60 --emit-units --shared-exon-columns --shared-exon-union --out $RUN/kbcolsu.fam \
        > $RUN/kbcolsu.families.log 2>&1
      ;;
  esac
done
echo "[job] $(date) done"
