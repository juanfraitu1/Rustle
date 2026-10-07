#!/bin/bash
# bench/entangled/run.sh — runner of docs/PREREG_entangled_baseline_2026-10-06.md. FAMILY = NPIP | TBC1D3, REP = 1 | 2 (the ideal-expression read sets, bench/ideal_expression).
#   run.sh tools FAMILY REP    StringTie (-L, no -G) and FLAIR (bam2bed + collapse --trust_ends, no correct, no annotation) on the ideal BAM, the lab's recipes (benchmark_collapse/run_*.sbatch)
#   run.sh units FAMILY REP    the driver's assemble and families stages with RUSTLE_BRIDGE_REGROUP=f1units (set for these two commands only)
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
cmd=${1:?tools|units|score|pool}
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
score)
  args=(--arm "D_asm=$D/asm.gtf" --arm "D_fam=$D/asm.families.gtf")
  [ -s "$O/st.gtf" ] && args+=(--arm "S=$O/st.gtf")
  [ -s "$O/fl.isoforms.gtf" ] && args+=(--arm "F=$O/fl.isoforms.gtf")
  [ -s "$O/units.gtf" ] && args+=(--arm "U_asm=$O/units.gtf")
  [ -s "$O/units.families.gtf" ] && args+=(--arm "U_fam=$O/units.families.gtf")
  light python3 "$HERE/arms_score.py" --truth "$D/reads" "${args[@]}" --out "$O/arms" | tee "$O/arms.log" ;;
*) echo "usage: run.sh tools|units|score FAMILY REP | pool" >&2; exit 2 ;;
esac
