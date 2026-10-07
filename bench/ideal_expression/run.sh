#!/bin/bash
# bench/ideal_expression/run.sh — the runner of docs/PREREG_ideal_expression_2026-10-06.md. FAMILY = NPIP | TBC1D3, REP = 1 | 2 (read seeds 20261006 | 20261007).
#   run.sh truth FAMILY REP     truth tables + the ideal reads (sim_windows.py), then G1 (verify_g1.py)          -> W/FAMILY/repREP/reads.*
#   run.sh map FAMILY REP       map the reads (whole-genome splice index, read-disjoint parts, one part per call; re-run until "mapping complete")
#   run.sh strata FAMILY REP    E0, the reachable stratum R and its sha1, G2 — from the BAM alone, BEFORE the pipeline
#   run.sh asm FAMILY REP       the driver's `assemble` and `families` stages with the HEAD binaries
#   run.sh anno FAMILY REP      G6 positive control: the canonicalized annotation as loci through the driver's `families` stage
#   run.sh score FAMILY REP     score.py (+ G7) for the arm; `score-anno` for the control
#   run.sh both FAMILY          verdict.py over the two replicates
# Environment: RS_BIN (HEAD release dir), RS_WORK (products), RS_PAD (smoke tests only), RS_PARTS. Nothing is read from a RUSTLE_* variable of the calling shell.
set -euo pipefail
HERE=$(cd "$(dirname "$0")" && pwd); REPO=$(cd "$HERE/../.." && pwd)
BIN=${RS_BIN:-/mnt/linuxdisk/home/juanfraitu/rustle_target_m2/release}
W=${RS_WORK:-/mnt/linuxdisk/tmp/ideal_expression_2026-10-06}
HB=/mnt/linuxdisk/home/juanfraitu
FA=$HB/winloci_data/chm13v2.0.fa
INDEX=$HB/npip_ladder/idx/target.splice.mmi
export TMPDIR=$W/tmp; mkdir -p "$TMPDIR"
export RLOCK_TIMEOUT=${RLOCK_TIMEOUT:-550}
leak=$(env | grep '^RUSTLE_' || true)
[ -z "$leak" ] || { echo "[ideal_expression] refusing: RUSTLE_* set in the calling environment: $leak" >&2; exit 2; }
sha() { sha1sum "$1" | cut -d' ' -f1; }
stamp() { for b in copy_assign mcl_families as_table family_score; do echo "$b	$(sha "$BIN/$b")"; done; for f in sim_windows.py strata.py score.py anno_loci.py verify_g1.py map_reads.py verdict.py; do echo "$f	$(sha "$HERE/$f")"; done; echo "head	$(git -C "$REPO" rev-parse --short HEAD)"; }
light() { bash "$REPO/tools/rlock.sh" light "$@"; }
heavy() { bash "$REPO/tools/rlock.sh" heavy "$@"; }
cmd=${1:?truth|map|strata|asm|anno|score|score-anno|both}; F=${2:?family}
R=${3:-1}
case "$R" in 1) SEED=20261006;; 2) SEED=20261007;; *) echo "REP must be 1 or 2" >&2; exit 2;; esac
D=$W/$F/rep$R; mkdir -p "$D"
P=$D/reads          # truth tables, reads, BAM
S=$D/strata         # E0, R, G2
A=$D/asm            # the driver's prefix for the arm
C=$D/anno           # the driver's prefix for the G6 control
case "$cmd" in
truth)
  { echo "date	$(date -Is)"; echo "seed	$SEED"; stamp; } > "$D/truth.run.log"
  light python3 "$HERE/sim_windows.py" --family "$F" --out "$P" --seed "$SEED" ${RS_PAD:+--pad "$RS_PAD"} | tee "$D/sim.log"
  light python3 "$HERE/verify_g1.py" "$P" | tee "$D/g1.log" ;;
map)
  case "$F" in NPIP) PARTS=${RS_PARTS:-6};; *) PARTS=${RS_PARTS:-5};; esac
  { echo "date	$(date -Is)"; echo "parts	$PARTS"; stamp; } >> "$P.map.run.log"
  heavy python3 "$HERE/map_reads.py" "$P" "$INDEX" --parts "$PARTS" --threads 4 --max-parts-per-call 1 2>&1 | tail -2 ;;
strata)
  [ -s "$P.bam" ] || { echo "map first" >&2; exit 2; }
  light python3 "$HERE/strata.py" --truth "$P" --bam "$P.bam" --out "$S" | tee "$D/strata.log"
  echo "[ideal_expression] R sha1: $(cat "$S.R.sha1")" ;;
asm)
  [ -s "$S.R.sha1" ] || { echo "run strata first (R is fixed before the pipeline)" >&2; exit 2; }
  { echo "date	$(date -Is)"; echo "R_sha1	$(cat "$S.R.sha1")"; echo "bam	$P.bam"; stamp; } > "$A.run.log"
  rm -f "$A".fam.* "$A.families.log"
  /usr/bin/time -v bash "$REPO/tools/rlock.sh" heavy bash "$REPO/tools/rustle_pipeline.sh" assemble --bam "$P.bam" --fasta "$FA" --out "$A" --bin "$BIN" --threads 4 > "$A.assemble.driver.log" 2> "$A.assemble.driver.stderr"
  /usr/bin/time -v bash "$REPO/tools/rlock.sh" heavy bash "$REPO/tools/rustle_pipeline.sh" families --bam "$P.bam" --fasta "$FA" --out "$A" --bin "$BIN" --threads 4 > "$A.families.driver.log" 2> "$A.families.driver.stderr"
  echo "[ideal_expression] asm $F rep$R: $(grep -h 'assemble:' "$A.assemble.driver.stderr" | tail -1 | sed 's/.*assemble: //'); $(grep -h 'families:' "$A.families.driver.stderr" | head -2 | sed 's/.*families: //' | tr '\n' ';')" ;;
anno)
  light python3 "$HERE/anno_loci.py" "$P" "$C"
  rm -f "$C".fam.* "$C.families.log"
  /usr/bin/time -v bash "$REPO/tools/rlock.sh" heavy bash "$REPO/tools/rustle_pipeline.sh" families --bam "$P.bam" --fasta "$FA" --out "$C" --bin "$BIN" --threads 4 > "$C.families.driver.log" 2> "$C.families.driver.stderr"
  echo "[ideal_expression] G6 control $F: $(grep -h 'families:' "$C.families.driver.stderr" | head -2 | sed 's/.*families: //' | tr '\n' ';')" ;;
score)
  light python3 "$HERE/score.py" --truth "$P" --strata "$S" --asm "$A" --single "$S.single_copy.tsv" --out "$D/score" | tee "$D/score.log" ;;
score-anno)
  light python3 "$HERE/score.py" --truth "$P" --strata "$S" --asm "$C" --single "$S.single_copy.tsv" --out "$D/score_anno" | tee "$D/score_anno.log" ;;
both)
  python3 "$HERE/verdict.py" --rep1 "$W/$F/rep1/score" --rep2 "$W/$F/rep2/score" --out "$W/$F/verdict.json" ;;
*) echo "usage: run.sh truth|map|strata|asm|anno|score|score-anno|both FAMILY [REP]" >&2; exit 2 ;;
esac
