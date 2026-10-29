#!/bin/bash
#SBATCH --job-name=hg002_u1
#SBATCH --partition=singlenode
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=72
#SBATCH --time=24:00:00
#SBATCH --export=NONE
unset SLURM_EXPORT_ENV
set -euo pipefail
# PREREG U1-FRESH-HG002-INTRON (docs/LEDGER_O1.md, 2026-10-29). STAGES env selects stages.
export PATH=$HOME/.pixi/detached-envs/rustle-7434162673922408322/envs/default/bin:$PATH
RUN=$HOME/runs/hg002_kinnex; BIN=$HOME/rustle_target/release; REPO=$HOME/Rustle
FASTA=$HOME/data/hg002/hg002.pat.fa; SRC_BAM=$HOME/work/jfim/rustle_data/reads/HG002_kinnex.t2t.bam
BAM=$RUN/hg002.pat.bam; T=72
STAGES=${STAGES:-align assemble families templates lift}
cd $RUN
for st in $STAGES; do
  echo "[job] $(date) $st (t+${SECONDS}s)"
  case $st in
    align)
      samtools fastq -@ 8 -F 0x900 $SRC_BAM 2> $RUN/fastq.log | minimap2 -ax splice:hq -uf --eqx -Y -N 50 -p 0.1 --secondary=yes -t $T $FASTA - 2> $RUN/mm2.log \
        | samtools sort -@ 16 -m 2G -o $BAM - && samtools index $BAM
      samtools flagstat $BAM > $RUN/hg002.pat.flagstat ;;
    assemble)
      bash $REPO/tools/rustle_pipeline.sh assemble --bam $BAM --fasta $FASTA --out $RUN/hg --threads $T --bin $BIN > $RUN/hg.assemble.log 2>&1 ;;
    families)
      bash $REPO/tools/rustle_pipeline.sh families --bam $BAM --fasta $FASTA --out $RUN/hg --threads $T --bin $BIN > $RUN/hg.families.log 2>&1 ;;
    templates)
      mkdir -p $RUN/edgeh/hg
      python3 $REPO/bench/denovo_eval/liftmode.py templates $RUN/hg.families.gtf $FASTA $RUN/edgeh/hg > $RUN/edgeh/templates.log 2>&1
      python3 - $RUN/edgeh/hg <<'PY'
import json, sys, collections
d = sys.argv[1]; M = json.load(open(f'{d}/meta.json'))['templates']
by = collections.defaultdict(list)
for n, t in M.items(): by[t['locus']].append(n)
keep = {max(ns, key=lambda n: (M[n]['reads'], M[n]['len'])) for ns in by.values()}
ok = False
with open(f'{d}/reps.fa', 'w') as w:
    for line in open(f'{d}/templates.fa'):
        if line.startswith('>'): ok = line[1:].strip() in keep
        if ok: w.write(line)
print('reps', len(keep))
PY
      cd $RUN/edgeh/hg
      minimap2 -c -x splice -N 100 -p 0.1 --secondary=yes -t 36 $FASTA reps.fa > S_rep.paf 2> S_rep.log &
      minimap2 -c -x asm20 -N 100 -p 0.1 --secondary=yes -t 36 $FASTA reps.fa > A_rep.paf 2> A_rep.log &
      wait; cd $RUN ;;
    lift)
      cd $RUN/edgeh/hg
      LIFT_MASK=rep LIFT_RETRO=1 LIFT_RECIPROCAL=1 python3 $REPO/bench/denovo_eval/liftmode.py families . S_rep.paf A_rep.paf > lift.log 2>&1
      awk '$3=="gene"{match($9,/ID=gene-([^;]+)/,m); print m[1]"\t"$1":"$4"-"$5}' $RUN/hg.fam.loci.gff3 > keymap.tsv
      python3 - <<'PY'
key = dict(l.rstrip('\n').split('\t') for l in open('keymap.tsv'))
with open('edges.u1.tsv', 'w') as w:
    for l in open('edges.tsv'):
        a, b, wt, ident = l.rstrip('\n').split('\t'); w.write(f'{key[a]}\t{key[b]}\t{wt}\t{ident}\n')
PY
      cd $RUN ;;
  esac
done
echo "[job] $(date) done (t+${SECONDS}s)"
