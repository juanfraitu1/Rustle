#!/bin/bash
#SBATCH --account=open
#SBATCH --partition=basic
#SBATCH --cpus-per-task=64
#SBATCH --mem=240G
#SBATCH --time=4:00:00
#SBATCH --job-name=r1829_score
# PREREG r1829 layer 2, part B. Usage: sbatch l2_score.sbatch [SEL.bam OUTDIR]   (default: the full selection)
source /storage/group/kdm16/default/jxi21/fritz/env.sh
W=$FRITZ/ws/r1829_near_tie_census/layer2; cd $W
SEL=${1:-$W/sel.bam}; O=${2:-$W/truth}; mkdir -p $O
P="pixi run -q --manifest-path $FRITZ/Rustle/pixi.toml python"; T=$FRITZ/Rustle/bench/o2_truth
REF=$HOME/data/chm13/chm13v2.0.fa; MAT=$HOME/ws/r1776_hg002_room/mat.mmi; PAT=$HOME/ws/r1776_hg002_room/pat.mmi
B=$HOME/reads/HSA/HG002_kinnex.t2t.bam; FAMFA=$HOME/ws/r1776_hg002_room/hg.sel.fa
export NT_REPAIR=1 NT_THREADS=32 NT_CAP=30000 TS_STRATA=S1,S2 AG_GUARD=1000000 OP_THREADS=32 TMPDIR=${TMPDIR:-/tmp}
set -o pipefail
$P $T/neutral_truth.py score $REF $MAT $PAT $O $SEL > $O/score.log 2>&1 || { echo "score failed"; tail $O/score.log; exit 1; }
echo "score done $SECONDS s"; cat $O/score.log | tail -8
# read -> family name, for the tie split (family reads = neutral_truth's selection of residues 0-24)
$P - $FAMFA $O/reads.tsv $O/family.tsv <<'EOF'
import sys
fam = {l[1:].split()[0] for l in open(sys.argv[1]) if l[0] == '>'}
with open(sys.argv[3], 'w') as w:
    for l in open(sys.argv[2]):
        q = l.split('\t')[0]
        if q != 'read': w.write(f"{q}\t{'family' if q in fam else 'other'}\n")
EOF
$P $T/tie_split.py $O/reads.tsv $O/reads.fa $O/placements.fa $MAT $PAT $O/family.tsv $O/tie_split.tsv > $O/tie_split.log 2>&1 || { echo "tie split failed"; tail $O/tie_split.log; exit 1; }
echo "tie split done $SECONDS s"; cat $O/tie_split.log
( cd $T && $P allele_graph.py regions $O/reads.tsv $O/regions.bed ) || { echo "regions failed"; exit 1; }
OP_CALLS=$O/own_calls.tsv $P $T/own_psv.py $O/reads.tsv $SEL $B $O/regions.bed > $O/own_psv.log 2>&1 || { echo "own psv failed"; tail $O/own_psv.log; exit 1; }
echo "own psv done $SECONDS s"
ANN=$HOME/ws/r1797_famsim_mimic/annot/chm13v2.0_CAT_Liftoff.slim.gff3.gz
$P $T/near_tie_classes.py $O/reads.tsv $SEL $O/tie_split.tsv $O/own_calls.tsv $ANN $REF $FAMFA $O/classes || { echo "classes failed"; exit 1; }
echo "all done $SECONDS s"
