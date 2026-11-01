#!/bin/bash
#SBATCH -J o2fast
#SBATCH -p spr1tb
#SBATCH -t 01:30:00
#SBATCH -o /home/hpc/v138da/v138da13/o2levers_fast/fast.%j.out
#SBATCH -e /home/hpc/v138da/v138da13/o2levers_fast/fast.%j.err
# O2 LEVERS fast panel (register 1481 amendment): base and --union-certificate on the same 57-family panel, side by side
set -u; cd $HOME/o2levers_fast
export PATH=$HOME/.pixi/detached-envs/rustle-7434162673922408322/envs/default/bin:$PATH
BIN=$HOME/rustle_target/release; R=$HOME/runs/human_A119b; S=$HOME/o2levers_A119b
BAM=$HOME/work/jfim/rustle_data/reads/A119b.t2t.bam; FASTA=$HOME/data/chm13/chm13v2.0.fa
run() { tag=$1; shift; st=$(date +%s)
  $BIN/copy_assign --bam $BAM --fasta $FASTA --regions panel.regions.txt --families $S/a1.fam.copies.tsv --copies-fa $S/a1.fam.copies.fa \
    --only-families panel.families.txt --region-threads 8 --skip-poa-diagnostic --tied-prefilter "$@" --out $tag > $tag.log 2>&1
  echo "[job] $(date) $tag exit $? $(( $(date +%s) - st )) s"; }
run base &
run union --union-certificate &
wait
python3 $HOME/Rustle/bench/o2levers/levers_score.py base.assignments.tsv union.assignments.tsv --names $HOME/o2diag/a1.lt5.names \
  --base-log base.log --log union.log > union.score.tsv 2> union.score.err; cat union.score.tsv
echo "[job] $(date) done"
