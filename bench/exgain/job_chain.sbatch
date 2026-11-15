#!/bin/bash
#SBATCH -p singlenode -c 72 -t 04:00:00 -J r1694chain
#SBATCH -o /anvme/workspace/v138da13-rustle/r1694_exgain/log/chain.%j.log
# r1694 rebuild chain; each stage is skipped when its output exists, each gate is a hard stop (Amendment 1).
set -euo pipefail
export PATH=$HOME/.pixi/detached-envs/rustle-7434162673922408322/envs/default/bin:$PATH
S=$HOME/Rustle/bench/exgain; W=/anvme/workspace/v138da13-rustle/r1694_exgain; D=$HOME/work/jfim/rustle_data
BAM=$D/reads/GCA_029281585.2_flnc_mm.bam; SEX=$S/SEX.tsv; PYE=$HOME/Rustle/.pixi/envs/bench/bin/python3  # bench env: edlib
declare -A FA=([ref]=$D/genomes/GGO.fa [mat]=$D/genomes/KB3781.mat.fa.gz [pat]=$D/genomes/KB3781.pat.fa.gz)
cd $W
for g in ref mat pat; do   # step 2: copy transcripts on each assembly
  [ -s hapcn_$g.paf ] || { minimap2 -c -x splice:hq -uf -N 50 -p 0.5 -I 16G -t 72 ${FA[$g]} cat/kb.cat.copies.fa > hapcn_$g.paf.tmp 2> log/hapcn_$g.mm2.log && mv hapcn_$g.paf.tmp hapcn_$g.paf; }
done
[ -s exons.fa ] || python3 $S/ex_units.py cat/kb.cat.copies.tsv hapcn_ref.paf $SEX ${FA[ref]} $W | tee exons.log
# --secondary=yes: minimap2 2.30's sr preset prints no secondary records, so -N 500 had no effect and no copy could be
# counted (run 4470947, kept in run1_nosecondary/); added after that run's G1/G2 were seen, disclosed in the report.
for g in ref mat pat; do   # step 3: exons on each assembly
  [ -s ex_sr_$g.paf ] || { minimap2 -c -x sr --secondary=yes -N 500 -p 0.5 -I 16G -t 72 ${FA[$g]} exons.fa > ex_sr_$g.paf.tmp 2> log/ex_sr_$g.mm2.log && mv ex_sr_$g.paf.tmp ex_sr_$g.paf; }
done
grep -h "loaded/built the index" log/*.mm2.log | sort | uniq -c   # one index part per run expected
python3 $S/ex_count.py $W $SEX | tee gates12.txt   # G1, G2: exit 3 stops the chain
[ -s units.tsv ] || python3 $S/exhap_units.py $W $SEX | tee units.log
if [ ! -s net.tsv ]; then
  for k in $(seq 0 35); do python3 $S/exhap_net.py $W $BAM $k 36 & done; wait
  { printf 'family\tread\tused\trole\tvia\tclass\tin_reads_tsv\n'; cat net_{0..35}.tsv.tmp | sort -s -k1,1; } > net.tsv.tmp && mv net.tsv.tmp net.tsv && rm net_*.tsv.tmp
fi
[ -s exceil.tsv ] && grep -q PASS gate3.txt || $PYE $S/exceil.py $W $SEX $BAM ${FA[ref]} ${FA[mat]} ${FA[pat]} | tee gate3.txt   # G3: exit 3 stops
[ -s exhap.families.tsv ] || { ./copy_phasing run --bam $BAM --copies copies.tsv --genome ${FA[ref]} --net net.tsv --out exhap.tmp --threads 72 > log/copy_phasing.log 2>&1
  for f in exhap.tmp.*; do mv $f exhap.${f#exhap.tmp.}; done; }
python3 $S/join.py $W exhap | tee join.txt
echo "chain done ${SECONDS}s"
