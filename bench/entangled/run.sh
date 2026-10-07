#!/bin/bash
# bench/entangled/run.sh — runner of docs/PREREG_entangled_baseline_2026-10-06.md. FAMILY = NPIP | TBC1D3, REP = 1 | 2 (the ideal-expression read sets, bench/ideal_expression).
#   run.sh tools FAMILY REP    StringTie (-L, no -G) and FLAIR (bam2bed + collapse --trust_ends, no correct, no annotation) on the ideal BAM, the lab's recipes (benchmark_collapse/run_*.sbatch)
#   run.sh units FAMILY REP    the driver's assemble and families stages with RUSTLE_BRIDGE_REGROUP=f1units (set for these two commands only)
#   run.sh levers FAMILY REP   Amendment 1: the driver's assemble and families stages with one existing opt-in change each: P = --no-seed-secondaries, C = RUSTLE_POLISH_SUBCHAIN=drop, PC = both (set for these commands only)
#   run.sh units2 FAMILY REP   docs/PREREG_locus_units_2026-10-06.md: arms Q (Rule 1), S_D (Rule 2 on all), S_Q (Rules 1+2) = locus_units.py on the families input, then mcl_families with the driver's flags, family-level scoring (score.py), Compara F (fam_score.py)
#   run.sh f1 FAMILY REP       Amendment 2: the driver with RUSTLE_BRIDGE_REGROUP=f1 (assemble + families), then family-level scoring (arm F1)
#   run.sh f1q FAMILY REP      Amendment 2: F1's families input with Rule 1 (primary chains of the primaries-only assembly), mcl_families, scoring (arm F1Q)
#   run.sh hier FAMILY REP     Amendment 2: for every arm, the pre-MCL graph (mcl_families --dump-graph, clusters must equal the registered ones), its components as clusters (hier_clusters.py), score.py and fam_score.py at level C
#   run.sh score FAMILY REP    arms_score.py over every arm whose GTF exists
#   run.sh pool                pool.py over the four runs
# Environment: RS_BIN (release dir), RS_IDEAL (the ideal-expression products), RS_WORK (this study's products). Nothing is read from a RUSTLE_* variable of the calling shell.
set -euo pipefail
HERE=$(cd "$(dirname "$0")" && pwd); REPO=$(cd "$HERE/../.." && pwd)
BIN=${RS_BIN:-/mnt/linuxdisk/home/juanfraitu/rustle_target_m2/release}
IDEAL=${RS_IDEAL:-/mnt/linuxdisk/tmp/ideal_expression_2026-10-06}
W=${RS_WORK:-/mnt/linuxdisk/tmp/entangled_2026-10-06}
FA=/mnt/linuxdisk/home/juanfraitu/winloci_data/chm13v2.0.fa
ST=/home/juanfra/miniforge3/bin/stringtie
FLENV=/home/juanfra/miniforge3/envs/flair/bin
export TMPDIR=$W/tmp; mkdir -p "$TMPDIR"
export RLOCK_TIMEOUT=${RLOCK_TIMEOUT:-550}
leak=$(env | grep '^RUSTLE_' || true)
[ -z "$leak" ] || { echo "[entangled] refusing: RUSTLE_* set in the calling environment: $leak" >&2; exit 2; }
light() { bash "$REPO/tools/rlock.sh" light "$@"; }
heavy() { bash "$REPO/tools/rlock.sh" heavy "$@"; }
cmd=${1:?tools|units|levers|units2|f1|f1q|hier|score|pool}
if [ "$cmd" = pool ]; then python3 "$HERE/pool.py" "$W"; exit 0; fi
F=${2:?family}; R=${3:?rep}
D=$IDEAL/$F/rep$R; O=$W/$F/rep$R; mkdir -p "$O/tmp"
case "$cmd" in
tools)
  { echo "date	$(date -Is)"; echo "stringtie	$($ST --version)"; echo "flair	$($FLENV/flair --version 2>&1 | head -1)"; echo "minimap2	$(minimap2 --version)"; echo "bam	$D/reads.bam"; echo "reads	$D/reads.fq"; } > "$O/tools.run.log"
  [ -s "$O/st.gtf" ] || light $ST -L -p 4 -o "$O/st.gtf.tmp" -A "$O/st.abund" "$D/reads.bam" > "$O/st.log" 2>&1
  [ -s "$O/st.gtf" ] || mv "$O/st.gtf.tmp" "$O/st.gtf"
  export PATH="$HERE/flair_shims:$FLENV:$PATH"
  [ -s "$O/fl.bed" ] || { light "$HERE/flair_shims/flair_bam2bed.py" -b "$D/reads.bam" -o "$O/fl.tmp" > "$O/fl.bam2bed.log" 2>&1; mv "$O/fl.tmp.bed" "$O/fl.bed"; }
  [ -s "$O/fl.isoforms.gtf" ] || heavy flair collapse -q "$O/fl.bed" -g "$FA" -r "$D/reads.fq" -o "$O/fl" -t 4 --trust_ends --generate_map --temp_dir "$O/tmp" > "$O/fl.collapse.log" 2>&1
  echo "[entangled] $F rep$R: stringtie $(awk '$3=="transcript"' "$O/st.gtf" | wc -l) transcripts; flair $(awk '$3=="transcript"' "$O/fl.isoforms.gtf" | wc -l) isoforms" ;;
units)
  A=$O/units; { echo "date	$(date -Is)"; echo "bam	$D/reads.bam"; for b in copy_assign mcl_families as_table; do echo "$b	$(sha1sum "$BIN/$b" | cut -d' ' -f1)"; done; echo "mode	RUSTLE_BRIDGE_REGROUP=f1units"; } > "$A.run.log"
  rm -f "$A".fam.* "$A.families.log"
  /usr/bin/time -v env RUSTLE_BRIDGE_REGROUP=f1units bash "$REPO/tools/rlock.sh" heavy bash "$REPO/tools/rustle_pipeline.sh" assemble --bam "$D/reads.bam" --fasta "$FA" --out "$A" --bin "$BIN" --threads 4 > "$A.assemble.driver.log" 2> "$A.assemble.driver.stderr"
  /usr/bin/time -v env RUSTLE_BRIDGE_REGROUP=f1units bash "$REPO/tools/rlock.sh" heavy bash "$REPO/tools/rustle_pipeline.sh" families --bam "$D/reads.bam" --fasta "$FA" --out "$A" --bin "$BIN" --threads 4 > "$A.families.driver.log" 2> "$A.families.driver.stderr"
  echo "[entangled] units $F rep$R: $(grep -h 'assemble:' "$A.assemble.driver.stderr" | tail -1 | sed 's/.*assemble: //')" ;;
levers)
  for arm in P C PC; do
    A=$O/lever_$arm; extra=(); envs=()
    case $arm in P) extra=(--no-seed-secondaries);; C) envs=(RUSTLE_POLISH_SUBCHAIN=drop);; PC) extra=(--no-seed-secondaries); envs=(RUSTLE_POLISH_SUBCHAIN=drop);; esac
    { echo "date	$(date -Is)"; echo "bam	$D/reads.bam"; for b in copy_assign mcl_families as_table; do echo "$b	$(sha1sum "$BIN/$b" | cut -d' ' -f1)"; done; echo "arm	$arm"; echo "extra	${extra[*]:-}"; echo "env	${envs[*]:-}"; } > "$A.run.log"
    rm -f "$A".fam.* "$A.families.log"
    /usr/bin/time -v env "${envs[@]}" bash "$REPO/tools/rlock.sh" heavy bash "$REPO/tools/rustle_pipeline.sh" assemble --bam "$D/reads.bam" --fasta "$FA" --out "$A" --bin "$BIN" --threads 4 "${extra[@]}" > "$A.assemble.driver.log" 2> "$A.assemble.driver.stderr"
    /usr/bin/time -v env "${envs[@]}" bash "$REPO/tools/rlock.sh" heavy bash "$REPO/tools/rustle_pipeline.sh" families --bam "$D/reads.bam" --fasta "$FA" --out "$A" --bin "$BIN" --threads 4 "${extra[@]}" > "$A.families.driver.log" 2> "$A.families.driver.stderr"
    echo "[entangled] lever $arm $F rep$R: $(grep -h 'assemble:' "$A.assemble.driver.stderr" | tail -1 | sed 's/.*assemble: //' | cut -c1-110)"
  done ;;
units2)
  PFAM=$O/lever_P.families.gtf; DFAM=$D/asm.families.gtf
  [ -s "$PFAM" ] || { echo "run levers first (the primaries-only families input)" >&2; exit 2; }
  case "$F" in NPIP) CH=chr16;; *) CH=chr17;; esac
  { echo "date	$(date -Is)"; echo "locus_units.py	$(sha1sum "$HERE/locus_units.py" | cut -d' ' -f1)"; echo "mcl_families	$(sha1sum "$BIN/mcl_families" | cut -d' ' -f1)"; echo "primary_gtf	$PFAM"; echo "base_gtf	$DFAM"; } > "$O/units2.run.log"
  python3 "$HERE/locus_units.py" --base "$DFAM" --out "$O/q.families.gtf" --primary-gtf "$PFAM" --rule1
  python3 "$HERE/locus_units.py" --base "$DFAM" --out "$O/sd.families.gtf" --rule2 all --side "$O/sd.separators.tsv"
  python3 "$HERE/locus_units.py" --base "$DFAM" --out "$O/sq.families.gtf" --primary-gtf "$PFAM" --rule1 --rule2 primary --side "$O/sq.separators.tsv"
  for arm in q sd sq; do
    A=$O/$arm; rm -f "$A".fam.* "$A.families.log"
    heavy "$BIN/mcl_families" --from-gtf "$A.families.gtf" --fasta "$FA" --threads 4 --min-exonic-bp 1 --min-shared-exon-frac 0.60 --emit-units --out "$A.fam" > "$A.families.log" 2>&1
    light python3 "$REPO/bench/ideal_expression/score.py" --truth "$D/reads" --strata "$D/strata" --asm "$A" --single "$D/strata.single_copy.tsv" --out "$O/score_$arm" > "$O/score_$arm.log" 2>&1
    light python3 "$REPO/bench/ideal_expression/fam_score.py" --fs "$BIN/family_score" --clusters "$A.fam.clusters.tsv" --contig "$CH" --windows "$D/reads.windows.tsv" --label "${F}_rep${R}_$arm" --out "$O/famscore_$arm.json" --work "$O/famscore_work_$arm" > "$O/famscore_$arm.log" 2>&1
    echo "[entangled] units2 $arm $F rep$R: $(grep -h 'families:' "$A.families.log" | head -1 | cut -c1-100) | $(tail -1 "$O/score_$arm.log" | cut -c1-120)"
  done ;;
f1)
  A=$O/f1; case "$F" in NPIP) CH=chr16;; *) CH=chr17;; esac
  { echo "date	$(date -Is)"; echo "bam	$D/reads.bam"; for b in copy_assign mcl_families as_table; do echo "$b	$(sha1sum "$BIN/$b" | cut -d' ' -f1)"; done; echo "mode	RUSTLE_BRIDGE_REGROUP=f1"; } > "$A.run.log"
  rm -f "$A".fam.* "$A.families.log"
  /usr/bin/time -v env RUSTLE_BRIDGE_REGROUP=f1 bash "$REPO/tools/rlock.sh" heavy bash "$REPO/tools/rustle_pipeline.sh" assemble --bam "$D/reads.bam" --fasta "$FA" --out "$A" --bin "$BIN" --threads 4 > "$A.assemble.driver.log" 2> "$A.assemble.driver.stderr"
  /usr/bin/time -v env RUSTLE_BRIDGE_REGROUP=f1 bash "$REPO/tools/rlock.sh" heavy bash "$REPO/tools/rustle_pipeline.sh" families --bam "$D/reads.bam" --fasta "$FA" --out "$A" --bin "$BIN" --threads 4 > "$A.families.driver.log" 2> "$A.families.driver.stderr"
  light python3 "$REPO/bench/ideal_expression/score.py" --truth "$D/reads" --strata "$D/strata" --asm "$A" --single "$D/strata.single_copy.tsv" --out "$O/score_f1" > "$O/score_f1.log" 2>&1
  light python3 "$REPO/bench/ideal_expression/fam_score.py" --fs "$BIN/family_score" --clusters "$A.fam.clusters.tsv" --contig "$CH" --windows "$D/reads.windows.tsv" --label "${F}_rep${R}_f1" --out "$O/famscore_f1.json" --work "$O/famscore_work_f1" > "$O/famscore_f1.log" 2>&1
  echo "[entangled] f1 $F rep$R: $(grep -h 'assemble:' "$A.assemble.driver.stderr" | tail -1 | sed 's/.*assemble: //' | cut -c1-100)" ;;
f1q)
  A=$O/f1q; case "$F" in NPIP) CH=chr16;; *) CH=chr17;; esac
  python3 "$HERE/locus_units.py" --base "$O/f1.families.gtf" --out "$A.families.gtf" --primary-gtf "$O/lever_P.families.gtf" --rule1
  rm -f "$A".fam.* "$A.families.log"
  heavy "$BIN/mcl_families" --from-gtf "$A.families.gtf" --fasta "$FA" --threads 4 --min-exonic-bp 1 --min-shared-exon-frac 0.60 --emit-units --out "$A.fam" > "$A.families.log" 2>&1
  light python3 "$REPO/bench/ideal_expression/score.py" --truth "$D/reads" --strata "$D/strata" --asm "$A" --single "$D/strata.single_copy.tsv" --out "$O/score_f1q" > "$O/score_f1q.log" 2>&1
  light python3 "$REPO/bench/ideal_expression/fam_score.py" --fs "$BIN/family_score" --clusters "$A.fam.clusters.tsv" --contig "$CH" --windows "$D/reads.windows.tsv" --label "${F}_rep${R}_f1q" --out "$O/famscore_f1q.json" --work "$O/famscore_work_f1q" > "$O/famscore_f1q.log" 2>&1
  echo "[entangled] f1q $F rep$R done" ;;
hier)
  case "$F" in NPIP) CH=chr16;; *) CH=chr17;; esac
  for arm in ${HIER_ARMS:-D P PC q sd sq f1 f1q}; do
    case $arm in D) PRE=$D/asm;; P|PC|C) PRE=$O/lever_$arm;; *) PRE=$O/$arm;; esac
    [ -s "$PRE.fam.clusters.tsv" ] || { echo "[entangled] hier $arm $F rep$R: no arm"; continue; }
    G=$O/hier_$arm.graph.tsv
    heavy "$BIN/mcl_families" --from-gtf "$PRE.families.gtf" --fasta "$FA" --threads 4 --min-exonic-bp 1 --min-shared-exon-frac 0.60 --emit-units --dump-graph "$G" --out "$O/hier_$arm.hg" > "$O/hier_$arm.log" 2>&1
    if ! cmp -s "$O/hier_$arm.hg.fam.clusters.tsv" "$PRE.fam.clusters.tsv"; then echo "[entangled] hier $arm $F rep$R: INVALID (the re-run clusters differ from the registered ones)"; continue; fi
    python3 "$HERE/hier_clusters.py" --loci "$PRE.fam.loci.gff3" --graph "$G" --out "$O/hierC_$arm" --prefix "$PRE" > /dev/null
    light python3 "$REPO/bench/ideal_expression/score.py" --truth "$D/reads" --strata "$D/strata" --asm "$O/hierC_$arm" --single "$D/strata.single_copy.tsv" --out "$O/scoreC_$arm" > "$O/scoreC_$arm.log" 2>&1
    light python3 "$REPO/bench/ideal_expression/fam_score.py" --fs "$BIN/family_score" --clusters "$O/hierC_$arm.fam.clusters.tsv" --contig "$CH" --windows "$D/reads.windows.tsv" --label "${F}_rep${R}_C_$arm" --out "$O/famscoreC_$arm.json" --work "$O/famscoreC_work_$arm" > "$O/famscoreC_$arm.log" 2>&1
    echo "[entangled] hier $arm $F rep$R: $(grep -h 'graph:' "$O/hier_$arm.log" | head -1 | cut -c1-80) | $(tail -1 "$O/scoreC_$arm.log" | cut -c1-60)"
  done ;;
score)
  args=(--arm "D_asm=$D/asm.gtf" --arm "D_fam=$D/asm.families.gtf")
  [ -s "$O/st.gtf" ] && args+=(--arm "S=$O/st.gtf")
  [ -s "$O/fl.isoforms.gtf" ] && args+=(--arm "F=$O/fl.isoforms.gtf")
  [ -s "$O/units.gtf" ] && args+=(--arm "U_asm=$O/units.gtf")
  [ -s "$O/units.families.gtf" ] && args+=(--arm "U_fam=$O/units.families.gtf")
  for arm in q sd sq; do
    [ -s "$O/$arm.families.gtf" ] && args+=(--arm "${arm}_fam=$O/$arm.families.gtf")
  done
  for arm in P C PC; do
    [ -s "$O/lever_$arm.gtf" ] && args+=(--arm "${arm}_asm=$O/lever_$arm.gtf")
    [ -s "$O/lever_$arm.families.gtf" ] && args+=(--arm "${arm}_fam=$O/lever_$arm.families.gtf")
  done
  light python3 "$HERE/arms_score.py" --truth "$D/reads" "${args[@]}" --out "$O/arms" | tee "$O/arms.log" ;;
*) echo "usage: run.sh tools|units|levers|units2|f1|f1q|hier|score FAMILY REP | pool" >&2; exit 2 ;;
esac
