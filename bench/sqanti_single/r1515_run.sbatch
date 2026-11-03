#!/bin/bash
#SBATCH -J r1515iso
#SBATCH -p singlenode
#SBATCH -t 05:00:00
#SBATCH -o /home/hpc/v138da/v138da13/isocon_only_A119b/run.%j.out
#SBATCH -e /home/hpc/v138da/v138da13/isocon_only_A119b/run.%j.err
# PREREG ISOCON-ONLY TRANSCRIPTS r1515 (docs/PREREG_isocon_only_2026-11-03.md): our chain on the 370 single-copy pools, all consensuses
set -u; D=$HOME/isocon_only_A119b; W=/anvme/workspace/v138da13-rustle/isocon_only; SQD=$HOME/sqanti_sc_A119b; cd $D
RUST=$HOME/.pixi/detached-envs/rustle-7434162673922408322/envs/default/bin; SQ=$(ls -d $HOME/.pixi/detached-envs/sqanti3-*/envs/default)/bin
export PATH=$HOME/.pixi/bin:$RUST:$PATH TMPDIR=$W/tmp; mkdir -p $TMPDIR $W/cache
BIN=$HOME/rustle_target_test/release; T0=$(date +%s); stamp() { echo "[job] $(date +%T) $1 (t+$(( $(date +%s) - T0 )) s)"; }
BAM=$HOME/work/jfim/rustle_data/reads/A119b.t2t.bam; FA=$HOME/data/chm13/chm13v2.0.fa; PY=$HOME/Rustle/bench/sqanti_single
IDX=$W/chm13.splicehq.mmi; [ -s $IDX ] || { minimap2 -x splice:hq -I 8G -t 32 -d $IDX $FA > $W/index.log 2>&1; stamp "index exit $?"; }
st=$(date +%s)
RUSTLE_CACHE_DIR=$W/cache timeout 7200 $BIN/candidate_copies --bam $BAM --fasta $FA --copies $W/sc.copies.tsv --copies-fa $W/sc.copies.fa \
  --index $IDX --max-reads 5000 --threads 64 --out $W/sc.cand --all-clusters $W/all_clusters.fa > $W/cand.log 2>&1
stamp "candidate_copies exit $? wall $(( $(date +%s) - st )) s consensuses $(grep -c '>' $W/all_clusters.fa 2>/dev/null) $(grep -o 'pass A .*' $W/cand.log | head -1)"
python3 $PY/isocon_only.py gtf $W/all_clusters.fa $IDX $D/arm.gtf
python3 $PY/sc_panel.py cut $SQD/targets.tsv $D/arm.gtf $D/arm.cut.gtf
stamp "gtf: $(grep -c '	transcript	' $D/arm.gtf) placed, cut $(grep -c '	transcript	' $D/arm.cut.gtf)"
export PATH=$SQ:$PATH
sqanti3_qc.py --isoforms $D/arm.cut.gtf --refGTF $SQD/chr20_ref_cat.gtf --refFasta $SQD/chr20.fa --report skip -t 16 -d $D/sq3/arm -o arm > $D/sq3_arm.log 2>&1 || echo "QC FAILED"
sqanti3_filter.py rules --sqanti_class $D/sq3/arm/arm_classification.txt --filter_gtf $D/sq3/arm/arm_corrected.gtf --skip_report -d $D/sq3/filt_arm -o arm > $D/filt_arm.log 2>&1 || echo "FILTER FAILED"
stamp "sqanti done"
export PATH=$RUST:$PATH
python3 $PY/sc_panel.py table $SQD/targets.tsv $(for a in ours stringtie flair isocon; do f=$(ls $SQD/sq3/filt_$a/*classification* | head -1); echo "$a=$SQD/sq3/$a/${a}_classification.txt:$f"; done) arm=$D/sq3/arm/arm_classification.txt:$(ls $D/sq3/filt_arm/*classification* | head -1) > $D/table.tsv 2> $D/table.err
cat $D/table.tsv
python3 $PY/isocon_only.py recover $D/only.tsv $D/sq3/arm/arm_classification.txt $D/recover.tsv
stamp done
