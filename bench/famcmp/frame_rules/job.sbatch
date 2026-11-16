#!/bin/bash
#SBATCH -p singlenode -N1 -n1 -c 72 -t 04:00:00 -J r1766 -o /anvme/workspace/v138da13-rustle/r1766_allmiss/job.%j.out
# bambu de novo (the arm r1763 lacks) + ours at a 1-read floor with today's defaults, PTR and PPA testis
export PATH=$HOME/.pixi/detached-envs/rustle-7434162673922408322/envs/default/bin:$HOME/.pixi/bin:$PATH
W=/anvme/workspace/v138da13-rustle; D=$W/r1766_allmiss; E=$W/bambu_env/pixi.toml; S=$HOME/Rustle/bench/famcmp
for t in PTR PPA; do n=${t,,}
  (s=$SECONDS; pixi run --manifest-path $E Rscript $S/bambu_run.R $W/reads/$t/${t}_mm.bam $W/genomes/$t/$t.fa $D/$n.bambu 20 $TMPDIR/rc_$n > $D/$n.bambu.log 2>&1; echo "$t bambu exit $? $((SECONDS-s)) s") &
done
for t in PTR PPA; do n=${t,,}; mkdir -p $D/floor1
  RUSTLE_MIN_CHAIN_READS=1 bash $D/snap/tools/rustle_pipeline.sh assemble --bam $W/reads/$t/${t}_mm.bam --fasta $W/genomes/$t/$t.fa --out $D/floor1/$n --threads 30 --bin $HOME/rustle_target/release > $D/floor1/$n.out 2>&1; echo "$t floor1 exit $? $SECONDS s $(awk '$3=="transcript"' $D/floor1/$n.gtf | wc -l) transcripts"
done
wait; echo done $SECONDS s
