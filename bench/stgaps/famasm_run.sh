#!/bin/bash
#SBATCH -J famasm
#SBATCH -p singlenode
#SBATCH -t 12:00:00
#SBATCH -o famasm.%j.out
#SBATCH -e famasm.%j.err
# FAMASM driver for one sample: famasm_run.sh SAMPLE PREFIX BAM GENOME RUN_DIR SPECIES
#   RUN_DIR/ref/{reftx.tsv,reftx.paf,ref.gtf} must exist (Task 1 of the plan). Env: THREADS (72), CONTIG (pool smoke only,
#   no scoring), OBSERVED=0 skips the observed.py BAM pass, SKIP_ASSIGN=1 refuses to run the assign stage, SKIP_POOL=1 reuses an existing RUN_DIR/pool/.
# Products: RUN_DIR/{pool,lift,admit}/, famasm.gtf, strata.txt, famasm.fam.*, mainstays.txt, pairs.txt, verdict.txt.
set -u
export PATH=$HOME/.pixi/detached-envs/gff-17846079296155688661/envs/default/bin:$HOME/.pixi/detached-envs/rustle-7434162673922408322/envs/default/bin:$PATH
export DENOVO_CAT_GFF=$HOME/data/denovo_eval/chm13v2.0_gencode.genes.gff3 DENOVO_CAT_TSV=$HOME/work/jfim/rustle_data/annotation/chm13v2.0_CAT_Liftoff.genes.tsv DENOVO_PARALOGS=$HOME/data/denovo_eval/paralogs.tsv DENOVO_HGNC=$HOME/data/denovo_eval/hgnc_complete_set.txt DENOVO_GGO_GFF=$HOME/data/denovo_eval/GGO_genomic.gff DENOVO_HG002_GFF=$HOME/data/hg002/genes.gff3 DENOVO_HG002_XREF=$HOME/data/hg002/xref.tsv
SAMPLE=$1; PREFIX=$2; BAM=$3; GENOME=$4; RUN=$5; SPECIES=$6
THREADS=${THREADS:-72}; CONTIG=${CONTIG:-}; OBSERVED=${OBSERVED:-1}; export FAMASM_THREADS=${FAMASM_THREADS:-$THREADS}
B=$HOME/Rustle/bench/stgaps; D=$HOME/Rustle/bench/denovo_eval; BIN=$HOME/rustle_target/release
COPIES=$PREFIX.fam.copies.tsv; BASE=$PREFIX.families.gtf; ASG=$PREFIX.assign.assignments.tsv
mkdir -p "$RUN"; cd "$RUN" || exit 2
printf 'SAMPLE=%s\nSPECIES=%s\nPREFIX=%s\nBASE_GTF=%s\nCOPIES=%s\nBAM=%s\nGENOME=%s\n' "$SAMPLE" "$SPECIES" "$PREFIX" "$BASE" "$COPIES" "$BAM" "$GENOME" > famasm.env
step() { echo "[famasm $(date +%H:%M:%S)] $*"; }

if [ -n "$CONTIG" ]; then   # smoke: pool one contig, nothing else
  step "pool --contig $CONTIG"; python3 $B/famasm.py pool "$COPIES" "$BAM" "$PREFIX.molecules.tsv" "$RUN" --contig "$CONTIG"; echo "pool exit $?"
  step "lift"; python3 $B/famasm.py lift "$COPIES" "$GENOME" "$RUN"; echo "lift exit $?"
  if [ "$(wc -l < "$ASG" 2>/dev/null || echo 0)" -gt 1 ]; then
    step "admit"; python3 $B/famasm.py admit "$RUN" "$ASG" "$COPIES"; echo "admit exit $?"
  fi
  exit 0
fi

# 1. certified reads: the assign stage on the families' copy table (skipped when its product is already there)
if [ "$(wc -l < "$ASG" 2>/dev/null || echo 0)" -le 1 ]; then
  [ "${SKIP_ASSIGN:-0}" = 1 ] && { echo "$ASG missing and SKIP_ASSIGN=1" >&2; exit 2; }
  step "assign stage -> $ASG"
  bash $HOME/Rustle/tools/rustle_pipeline.sh assign --bam "$BAM" --fasta "$GENOME" --out "$PREFIX" --threads "$THREADS" --bin "$BIN" > assign.log 2>&1; rc=$?; echo "assign exit $rc"
  [ "$rc" = 0 ] && [ "$(wc -l < "$ASG" 2>/dev/null || echo 0)" -gt 1 ] || { echo "assign stage failed or wrote no $ASG; see assign.log" >&2; exit 2; }
fi
# observed-chain counts (descriptive), in parallel with the pools
if [ "$OBSERVED" = 1 ] && [ ! -s ref/counts.tsv ]; then
  step "observed count (background)"; python3 $B/observed.py count "$BAM" ref/counts.tsv 2> ref/counts.err &
fi
# 2. FAMASM
fail() { echo "[famasm] $1 failed (exit $2): no verdict" >&2; echo "verdict INVALID: $1 failed" > verdict.txt; exit 2; }
if [ "${SKIP_POOL:-0}" = 1 ] && [ -s pool.stats.tsv ]; then echo "pool reused (SKIP_POOL=1)"; else
step "pool"; python3 $B/famasm.py pool "$COPIES" "$BAM" "$PREFIX.molecules.tsv" "$RUN" || fail pool $?; echo "pool exit 0"; fi
step "lift"; python3 $B/famasm.py lift "$COPIES" "$GENOME" "$RUN" || fail lift $?; echo "lift exit 0"
step "admit"; python3 $B/famasm.py admit "$RUN" "$ASG" "$COPIES" || fail admit $?; echo "admit exit 0"
step "gtf"; python3 $B/famasm.py gtf "$RUN" "$BASE" "$COPIES" famasm.gtf || fail gtf $?; echo "gtf exit 0"; cat gtf.stats.tsv
# 3. strata (M1-M3)
step "stratify"; python3 $B/stratify.py ref/reftx.tsv ref/reftx.paf ref/ref.gtf strata base="$BASE" famasm=famasm.gtf > strata.txt 2> stratify.err || fail stratify $?; echo "stratify exit 0"
# 4. families on the FAMASM GTF (M4)
step "mcl_families --from-gtf famasm.gtf"
RUSTLE_CACHE_DIR=$PREFIX.cache $BIN/mcl_families --from-gtf famasm.gtf --fasta "$GENOME" --threads "$THREADS" --min-exonic-bp 1 --min-shared-exon-frac 0.60 --emit-units --out famasm.fam > famasm.families.log 2>&1 || fail mcl_families $?; echo "mcl exit 0"
python3 $D/mainstay_score.py --species "$SPECIES" "$BASE" "$COPIES" famasm.fam.copies.tsv > mainstays.txt 2> mainstays.err || fail mainstay_score $?; echo "mainstay exit 0"
python3 $D/pairscore.py --species "$SPECIES" --gtf "$BASE" "$COPIES" famasm.fam.copies.tsv > pairs.txt 2> pairs.err || fail pairscore $?; echo "pairscore exit 0"
wait
if [ "$OBSERVED" = 1 ] && [ -s ref/counts.tsv ]; then
  python3 $B/observed.py score ref/ref.gtf strata/stratum.tsv ref/counts.tsv base=strata/gc_base.base.gtf.tmap famasm=strata/gc_famasm.famasm.gtf.tmap > observed.txt 2> observed.err; echo "observed exit $?"
fi
# 5. verdict
step "score"; python3 $B/famasm.py score "$RUN" > verdict.txt || fail score $?; cat verdict.txt
echo done
