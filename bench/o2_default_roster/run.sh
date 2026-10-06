#!/bin/bash
# bench/o2_default_roster/run.sh — the runner of docs/PREREG_o2_default_roster_2026-10-06.md: O2 with read-level truth on the DEFAULT families roster.
#
#   run.sh roster CONTIG   HEAD assemble (the driver's `assemble` flags, --region CONTIG) + the driver's `families` stage -> W/CONTIG/roster.fam.copies.{tsv,fa,regions};
#                          cmp against the stored e163d955 products where they exist
#   run.sh sim CONTIG      bench/sim.py copies on the roster against the whole-genome splice index, seed 20260925, 2 mapping parts, one part per call:
#                          re-run until it prints "sim complete"
#   run.sh assign CONTIG   the catalog (roster minus copies without reads), the driver's `assign` stage on the simulated BAM (arm O2) and
#                          copy_assign --union-certificate (arm U2, beside)
#   run.sh score CONTIG    bench/o2_default_roster/report.py (CONTIG: chr16 = NPIP, chr17 = TBC1D3, chrY = Y)
# Environment: RS_BIN (HEAD release dir), RS_WORK (products). Nothing is read from a RUSTLE_* variable of the calling shell.
set -euo pipefail
HERE=$(cd "$(dirname "$0")" && pwd); REPO=$(cd "$HERE/../.." && pwd)
BIN=${RS_BIN:-/mnt/linuxdisk/home/juanfraitu/rustle_target_m2/release}
W=${RS_WORK:-/mnt/linuxdisk/tmp/o2_default_roster_2026-10-06}
HB=/mnt/linuxdisk/home/juanfraitu
BAM=$HB/winloci_data/A119b.t2t.bam; FA=$HB/winloci_data/chm13v2.0.fa
INDEX=$HB/npip_ladder/idx/target.splice.mmi
MOL=/mnt/linuxdisk/tmp/rustle_figures/runs/human_A119b/human_A119b.molecules.tsv
CH=/mnt/linuxdisk/tmp/rustle_figures_dev/container_headroom/data/human_A119b
SEED=20260925
POLISH="--assembly-polish full --polish-isoform-fraction 0.02 --polish-mono-shadow --polish-mono-quantile 0.82 --polish-ism-ratio 0.7 --polish-retained-ratio 10"
export TMPDIR=$W/tmp; mkdir -p "$TMPDIR"
export RLOCK_TIMEOUT=${RLOCK_TIMEOUT:-550}
leak=$(env | grep '^RUSTLE_' || true)
[ -z "$leak" ] || { echo "[o2_default_roster] refusing: RUSTLE_* set in the calling environment: $leak" >&2; exit 2; }
grep -qF -- "POLISH=\"$POLISH\"" "$REPO/tools/rustle_pipeline.sh" || { echo "[o2_default_roster] the driver's POLISH differs from this runner's" >&2; exit 2; }
sha() { sha1sum "$1" | cut -d' ' -f1; }
stamp() { for b in copy_assign mcl_families; do echo "$b	$(sha "$BIN/$b")"; done; echo "sim.py	$(sha "$REPO/bench/sim.py")"; echo "score.py	$(sha "$REPO/bench/score.py")"; echo "head	$(git -C "$REPO" rev-parse --short HEAD)"; }
cmd=${1:?roster|sim|assign|score}; C=${2:?contig (chr16|chr17|chrY)}
D=$W/$C; mkdir -p "$D"
R=$D/roster          # the roster: HEAD default assemble + families on the contig
S=$D/sim             # the simulation
A=$D/asg             # the assignment (the driver's assign stage prefix)
case "$cmd" in
roster)
  mkdir -p "$W/mol"; [ -e "$W/mol/mol.tsv" ] || ln -s "$MOL" "$W/mol/mol.tsv"
  len=$(awk -v c="$C" '$1==c{print $2}' "$FA.fai"); [ -n "$len" ] || { echo "contig $C not in $FA.fai" >&2; exit 2; }
  { echo "date	$(date -Is)"; echo "contig	$C"; stamp; } > "$R.run.log"
  /usr/bin/time -v env RUSTLE_GTF_SECONDARY=1 RUSTLE_GTF_SECONDARY_AS_RATIO=0.98 "RUSTLE_GTF_SECONDARY_AS_TABLE=$W/mol/mol.tsv" \
    bash "$REPO/tools/rlock.sh" heavy "$BIN/copy_assign" --assemble-only --region "$C:0-$len" --assembly-junctions strict $POLISH --bridge-regroup f1v2 --gtf-tpm \
    --bam "$BAM" --fasta "$FA" --out "$R" > "$R.assemble.log" 2> "$R.assemble.stderr"
  rm -f "$R".fam.* "$R.families.log"
  /usr/bin/time -v bash "$REPO/tools/rlock.sh" heavy bash "$REPO/tools/rustle_pipeline.sh" families --bam "$BAM" --fasta "$FA" --out "$R" --bin "$BIN" --threads 4 \
    > "$R.driver.log" 2> "$R.driver.stderr"
  echo "[o2_default_roster] roster $C: $(awk 'NR>1' "$R.fam.copies.tsv" | wc -l) copies in $(awk 'NR>1{print $1}' "$R.fam.copies.tsv" | sort -u | wc -l) families"
  case "$C" in
    chr16|chr17) SA=$CH/asm_f1v2/$C/human_A119b.$C; SF=$CH/fam/$C/D;;
    *) SA=""; SF=$CH/repfam/$C/D;;
  esac
  {
    echo "provenance vs the e163d955 container-headroom products (cmp):"
    if [ -n "$SA" ]; then for p in "gtf|$R.gtf|$SA.gtf" "families.gtf|$R.families.gtf|$SA.families.gtf"; do IFS='|' read -r l x y <<< "$p"; if cmp -s "$x" "$y"; then echo "  $l: identical"; else echo "  $l: DIFFERS"; fi; done; fi
    for p in "clusters|$R.fam.clusters.tsv|$SF.fam.clusters.tsv" "loci.gff3|$R.fam.loci.gff3|$SF.fam.loci.gff3" "copies|$R.fam.copies.tsv|$SF.fam.copies.tsv"; do
      IFS='|' read -r l x y <<< "$p"; if [ ! -e "$y" ]; then echo "  $l: no stored product"; elif cmp -s "$x" "$y"; then echo "  $l: identical"; else echo "  $l: DIFFERS"; fi
    done
  } | tee "$R.provenance.txt" ;;
sim)
  [ -s "$R.fam.copies.tsv" ] || { echo "run roster first" >&2; exit 2; }
  { echo "date	$(date -Is)"; echo "seed	$SEED"; stamp; } >> "$S.run.log"
  # the stage's genome-wide mapping, one part per call (the index load, ~1 min, is paid per part); a finished simulation is marked by S.done
  bash "$REPO/tools/rlock.sh" heavy python3 "$REPO/bench/sim.py" copies "$R.fam.copies.tsv" "$R.fam.copies.fa" "$INDEX" "$S" "$SEED" \
    --threads 4 --parts 2 --max-parts-per-call 1 --reuse-fastq --sibling none >> "$S.log" 2>&1
  if [ -e "$S.done" ]; then echo "[o2_default_roster] sim complete: $(cat "$S.done")"; else echo "[o2_default_roster] sim: parts remain, re-run: $(tail -1 "$S.log")"; fi ;;
assign)
  [ -e "$S.done" ] || { echo "the simulation is not complete" >&2; exit 2; }
  python3 - "$REPO" "$R" "$S" "$D" <<'PY'
import os, shutil, sys
repo, R, S, D = sys.argv[1:5]
sys.path.insert(0, os.path.join(repo, "figures"))
import _o2
cat = _o2.derive_catalog(f"{R}.fam.copies.tsv", f"{R}.fam.copies.fa", S, f"{D}/cat", force=True)
shutil.copy(cat["tsv"], f"{D}/asg.fam.copies.tsv"); shutil.copy(cat["fa"], f"{D}/asg.fam.copies.fa"); shutil.copy(f"{R}.fam.copies.regions", f"{D}/asg.fam.copies.regions")
print(f"[o2_default_roster] catalog: {sum(1 for _ in open(cat['tsv'])) - 1} copies kept, {sum(1 for _ in open(cat['dropped'])) - 1} dropped (no read over their span)")
PY
  rm -f "$A".assign.*
  /usr/bin/time -v bash "$REPO/tools/rlock.sh" heavy bash "$REPO/tools/rustle_pipeline.sh" assign --bam "$S.bam" --fasta "$FA" --out "$A" --bin "$BIN" --threads 4 \
    > "$A.driver.log" 2> "$A.driver.stderr"
  grep -h "assign:" "$A.driver.stderr" | tail -2
  /usr/bin/time -v bash "$REPO/tools/rlock.sh" heavy "$BIN/copy_assign" --bam "$S.bam" --fasta "$FA" --regions "$A.regions.txt" --families "$A.fam.copies.tsv" \
    --copies-fa "$A.fam.copies.fa" --union-certificate --out "$A.union" > "$A.union.log" 2> "$A.union.stderr"
  echo "[o2_default_roster] assign $C: O2 $(grep Elapsed "$A.driver.stderr" | awk '{print $NF}'), U2 $(grep Elapsed "$A.union.stderr" | awk '{print $NF}')" ;;
score)
  case "$C" in chr16) T=NPIP;; chr17) T=TBC1D3;; chrY) T=Y;; *) T=none;; esac
  light() { bash "$REPO/tools/rlock.sh" light "$@"; }
  light python3 "$HERE/report.py" --work "$D" --contig "$C" --target "$T" --sim "$S" --catalog "$A.fam.copies.tsv" --o2 "$A.assign" --u2 "$A.union" --logs "$D/score" | tee "$D/report.txt" ;;
*) echo "usage: run.sh roster|sim|assign|score CONTIG" >&2; exit 2 ;;
esac
