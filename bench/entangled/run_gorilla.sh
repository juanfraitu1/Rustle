#!/bin/bash
# bench/entangled/run_gorilla.sh — runner of docs/PREREG_gorilla_overlap_2026-10-07.md (gorilla OR6737 testis, genome-wide).
#   run_gorilla.sh arms          our assembler D (default), P (--no-seed-secondaries), PC (P + RUSTLE_POLISH_SUBCHAIN=drop): the driver's assemble stage, one after the other
#   run_gorilla.sh annotation    valid annotated chains (GFF + FASTA)                         -> W/truth.annotation.pkl
#   run_gorilla.sh support CONTIGS   read support of the annotated chains on those contigs (comma list; empty = all)   -> W/truth.support.*.pkl
#   run_gorilla.sh tables        truth.chains.tsv, truth.genes.tsv, truth.stats.json
#   run_gorilla.sh score         real_score.py: D P PC and the lab's S F I, on the P1 and on the P2 denominator
# Environment: RS_BIN (release dir), RS_WORK. Nothing is read from a RUSTLE_* variable of the calling shell.
set -euo pipefail
HERE=$(cd "$(dirname "$0")" && pwd); REPO=$(cd "$HERE/../.." && pwd)
BIN=${RS_BIN:-/mnt/linuxdisk/home/juanfraitu/rustle_target_m2/release}
W=${RS_WORK:-/mnt/linuxdisk/tmp/gorilla_overlap_2026-10-07}
HB=/mnt/linuxdisk/home/juanfraitu
BAM=$HB/winloci_data/GGO_mm.bam; FA=$HB/_from_wsl/winloci_scratch/GGO.fasta; GFF=$HB/winloci_data/GGO_genomic.gff
LAB=/mnt/c/Users/jfris/Desktop
S_GTF=$LAB/benchmark_collapse/stringtie_GGO/GGO.stringtie.gtf; F_GTF=$LAB/benchmark_collapse/flair_GGO/GGO.flair.isoforms.gtf; I_GFF=$LAB/isoseq_upload/isoseq_GGO_OR6737/GGO_OR6737.collapsed.gff.gz
mkdir -p "$W" "$W/tmp"; export TMPDIR=$W/tmp
export RLOCK_TIMEOUT=${RLOCK_TIMEOUT:-550}
leak=$(env | grep '^RUSTLE_' || true)
[ -z "$leak" ] || { echo "[gorilla_overlap] refusing: RUSTLE_* set in the calling environment: $leak" >&2; exit 2; }
heavy() { bash "$REPO/tools/rlock.sh" heavy "$@"; }
cmd=${1:?arms|annotation|support|tables|score}
case "$cmd" in
arms)
  for arm in D P PC; do
    A=$W/$arm; extra=(); envs=()
    case $arm in P) extra=(--no-seed-secondaries);; PC) extra=(--no-seed-secondaries); envs=(RUSTLE_POLISH_SUBCHAIN=drop);; esac
    if [ "$arm" != D ] && [ -s "$W/D.molecules.tsv" ]; then cp -n "$W/D.molecules.tsv" "$A.molecules.tsv"; [ -s "$W/D.molecules.tsv.asbin" ] && cp -n "$W/D.molecules.tsv.asbin" "$A.molecules.tsv.asbin"; fi
    [ -s "$A.gtf" ] && { echo "[gorilla_overlap] $arm exists"; continue; }
    { echo "date	$(date -Is)"; for b in copy_assign as_table mcl_families; do echo "$b	$(sha1sum "$BIN/$b" | cut -d' ' -f1)"; done; echo "arm	$arm"; echo "extra	${extra[*]:-}"; echo "env	${envs[*]:-}"; } > "$A.run.log"
    /usr/bin/time -v env "${envs[@]}" bash "$REPO/tools/rlock.sh" heavy bash "$REPO/tools/rustle_pipeline.sh" assemble --bam "$BAM" --fasta "$FA" --out "$A" --bin "$BIN" --threads 4 "${extra[@]}" > "$A.assemble.driver.log" 2> "$A.assemble.driver.stderr"
    echo "[gorilla_overlap] $arm: $(grep -h 'assemble:' "$A.assemble.driver.stderr" | tail -1 | sed 's/.*assemble: //' | cut -c1-120); $(grep -h 'Elapsed' "$A.assemble.driver.stderr" | head -1)"
  done ;;
annotation)
  heavy python3 "$HERE/real_truth.py" --stage annotation --gff "$GFF" --fasta "$FA" --bam "$BAM" --as-table "$W/D.molecules.tsv" --out "$W/truth" ;;
support)
  heavy python3 "$HERE/real_truth.py" --stage support --gff "$GFF" --fasta "$FA" --bam "$BAM" --as-table "$W/D.molecules.tsv" --out "$W/truth" --contigs "${2:-}" ;;
tables)
  python3 "$HERE/real_truth.py" --stage tables --gff "$GFF" --fasta "$FA" --bam "$BAM" --as-table "$W/D.molecules.tsv" --out "$W/truth" ;;
score)
  [ -s "$W/I.collapsed.gff" ] || zcat "$I_GFF" > "$W/I.collapsed.gff"
  for den in P1 P2; do
    python3 "$HERE/real_score.py" --truth "$W/truth" --denominator $den --arm "D=$W/D.gtf" --arm "P=$W/P.gtf" --arm "PC=$W/PC.gtf" --arm "S=$S_GTF" --arm "F=$F_GTF" --arm "I=$W/I.collapsed.gff" --out "$W/score_$den" | tee "$W/score_$den.log"
  done ;;
*) echo "usage: run_gorilla.sh arms|annotation|support CONTIGS|tables|score" >&2; exit 2 ;;
esac
