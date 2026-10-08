#!/bin/bash
# bench/default_rescore/run.sh — the runner of docs/PREREG_default_rescore_npip_2026-10-06.md: the current default pipeline (HEAD binaries)
# on human A119b chr16, scored on the 25 CAT/Liftoff NPIP copies (strict FOUND, Amendment E) and on the U2 / Compara / Soto family truths.
#
#   run.sh gates         G0 (HEAD scorer on the frozen read-pool arms), G1 (own-node rule vs the stored pagedata flags), G3 (HEAD family_score
#                        on the stored e163d955 default clusters) -> W/g0.json, g1.json, g3.json
#   run.sh asm ARM       chr16 assembly with the HEAD copy_assign, the driver's `assemble` command with --region in place of --genome-wide
#                        (equal to the genome-wide run on the contig: docs/archive/2026-09/CONTAINER_HEADROOM_2026-09-30.md gates); ARM = DEF (the driver's
#                        default, --bridge-regroup f1v2) | PRE (--bridge-regroup off, the pre-flip control)
#   run.sh fam ARM       the driver's `families` stage on that assembly (PRE: RUSTLE_BRIDGE_REGROUP=off RUSTLE_MIN_COV_SHORTER=0)
#   run.sh score         own nodes, copy_support (DEF and PRE; twice, under two PYTHONHASHSEEDs), family_score (DEF, PRE), provenance diffs
#   run.sh verdict       bench/default_rescore/verdict.py
# Environment: RS_BIN (HEAD release dir), RS_WORK (products). Nothing is read from a RUSTLE_* variable of the calling shell.
set -euo pipefail
HERE=$(cd "$(dirname "$0")" && pwd); REPO=$(cd "$HERE/../.." && pwd)
BIN=${RS_BIN:-/mnt/linuxdisk/home/juanfraitu/rustle_target_m2/release}
W=${RS_WORK:-/mnt/linuxdisk/tmp/rescore_2026-10-06}
HB=/mnt/linuxdisk/home/juanfraitu
BAM=$HB/winloci_data/A119b.t2t.bam; FA=$HB/winloci_data/chm13v2.0.fa
MOL=/mnt/linuxdisk/tmp/rustle_figures/runs/human_A119b/human_A119b.molecules.tsv   # the genome-wide best-AS table of the 09-25 run
FROZEN=/mnt/linuxdisk/tmp/readpool_npip                                             # P / GOOD / ALL arms of 2026-10-01 and their stored scores
ANN=/mnt/linuxdisk/tmp/rustle_figures_dev/copy_recovery_tools_cat/ann               # the 25 CAT/Liftoff NPIP copies and their truth
CH=/mnt/linuxdisk/tmp/rustle_figures_dev/container_headroom                         # the e163d955 default products (D = f1v2 + cov .70, B0 = pre-flip)
POLISH="--assembly-polish full --polish-isoform-fraction 0.02 --polish-mono-shadow --polish-mono-quantile 0.82 --polish-ism-ratio 0.7 --polish-retained-ratio 10"
export TMPDIR=$W/tmp; mkdir -p "$TMPDIR"
export RLOCK_TIMEOUT=${RLOCK_TIMEOUT:-550}
leak=$(env | grep '^RUSTLE_' || true)
[ -z "$leak" ] || { echo "[default_rescore] refusing: RUSTLE_* set in the calling environment: $leak" >&2; exit 2; }
grep -qF -- "POLISH=\"$POLISH\"" "$REPO/tools/rustle_pipeline.sh" || { echo "[default_rescore] the driver's POLISH differs from this runner's" >&2; exit 2; }
sha() { sha1sum "$1" | cut -d' ' -f1; }
heavy() { bash "$REPO/tools/rlock.sh" heavy "$@"; }
light() { bash "$REPO/tools/rlock.sh" light "$@"; }
stamp() { for b in copy_assign mcl_families family_score; do echo "$b	$(sha "$BIN/$b")"; done; echo "copy_support.py	$(sha "$REPO/bench/copy_support.py")"; }
LOCI_FROZEN=(--loci "P=$FROZEN/P.gff3,$FROZEN/P.gtf" --loci "GOOD=$FROZEN/GOOD.gff3,$FROZEN/GOOD.gtf" --loci "ALL=$FROZEN/ALL.gff3,$FROZEN/ALL.gtf")
cmd=${1:?gates|asm|fam|score|verdict}; shift
mkdir -p "$W"
case "$cmd" in
gates)
  { echo "date	$(date -Is)"; echo "head	$(git -C "$REPO" rev-parse --short HEAD)"; stamp; } > "$W/gates.log"
  light python3 "$REPO/bench/copy_support.py" --copies "$ANN/copies.hsa.tsv" --truth "$ANN/truth.hsa.gtf" --bam "$BAM" --family NPIP --out "$W/g0" \
    "${LOCI_FROZEN[@]}" --nodes "$FROZEN/pagedata.json" > "$W/g0.stdout"
  # G1 on P and GOOD: the frozen ALL arm has two loci on one span (DN_chr16_33611806_2 and _3), which pagedata.py's span join resolved silently and nodes.py refuses
  python3 "$HERE/nodes.py" --copies "$ANN/copies.hsa.tsv" --truth "$ANN/truth.hsa.gtf" --family NPIP --exons-json "$FROZEN/npip_read_pool.json" \
    --arm "P=$FROZEN/P.gff3,$FROZEN/P.fam.clusters.tsv" --arm "GOOD=$FROZEN/GOOD.gff3,$FROZEN/GOOD.fam.clusters.tsv" \
    --out "$W/g1.nodes.json" --check "$FROZEN/pagedata.json" --check-out "$W/g1.json"
  python3 "$HERE/score_families.py" --fs "$BIN/family_score" --clusters "$CH/data/human_A119b/fam/chr16/D.fam.clusters.tsv" --contig chr16 --label g3_D --out "$W/g3.json" --work "$W/fs"
  echo "[default_rescore] gates written to $W/{g0,g1,g3}.json" ;;
asm)
  ARM=${1:?DEF|PRE}
  case "$ARM" in DEF) MODE=f1v2;; PRE) MODE=off;; *) echo "ARM must be DEF or PRE" >&2; exit 2;; esac
  D=$W/$ARM; mkdir -p "$D" "$W/mol"; [ -e "$W/mol/mol.tsv" ] || ln -s "$MOL" "$W/mol/mol.tsv"   # the .asbin sidecar lands here, not in the stored run
  len=$(awk '$1=="chr16"{print $2}' "$FA.fai"); P=$D/human_A119b.chr16
  { echo "date	$(date -Is)"; echo "head	$(git -C "$REPO" rev-parse --short HEAD)"; echo "arm	$ARM	bridge-regroup $MODE"; stamp; } > "$P.run.log"
  /usr/bin/time -v env RUSTLE_GTF_SECONDARY=1 RUSTLE_GTF_SECONDARY_AS_RATIO=0.98 "RUSTLE_GTF_SECONDARY_AS_TABLE=$W/mol/mol.tsv" \
    bash "$REPO/tools/rlock.sh" heavy "$BIN/copy_assign" --assemble-only --region "chr16:0-$len" --assembly-junctions strict $POLISH --bridge-regroup "$MODE" --gtf-tpm \
    --bam "$BAM" --fasta "$FA" --out "$P" > "$P.assemble.log" 2> "$P.assemble.stderr"
  echo "[default_rescore] asm $ARM: $(grep Elapsed "$P.assemble.stderr" | awk '{print $NF}') $(grep 'Maximum resident' "$P.assemble.stderr" | awk '{print $NF}') KB, $(awk -F'\t' '$3=="transcript"' "$P.gtf" | wc -l) transcripts" ;;
fam)
  ARM=${1:?DEF|PRE}
  case "$ARM" in DEF) ENVV=();; PRE) ENVV=(RUSTLE_BRIDGE_REGROUP=off RUSTLE_MIN_COV_SHORTER=0);; *) echo "ARM must be DEF or PRE" >&2; exit 2;; esac
  D=$W/$ARM; P=$D/human_A119b.chr16
  rm -f "$P".fam.* "$P.families.log"
  /usr/bin/time -v env "${ENVV[@]}" bash "$REPO/tools/rlock.sh" heavy bash "$REPO/tools/rustle_pipeline.sh" families --bam "$BAM" --fasta "$FA" --out "$P" --bin "$BIN" --threads 4 \
    > "$P.driver.log" 2> "$P.driver.stderr"
  echo "[default_rescore] fam $ARM: $(grep Elapsed "$P.driver.stderr" | awk '{print $NF}'); $(grep -h 'families:' "$P.driver.stderr" | sed 's/.*families: //' | tr '\n' ';')" ;;
score)
  PD=$W/DEF/human_A119b.chr16; PP=$W/PRE/human_A119b.chr16
  for f in "$PD.fam.loci.gff3" "$PD.families.gtf" "$PP.fam.loci.gff3" "$PP.gtf"; do [ -s "$f" ] || { echo "missing $f: run asm and fam for both arms" >&2; exit 2; }; done
  python3 "$HERE/nodes.py" --copies "$ANN/copies.hsa.tsv" --truth "$ANN/truth.hsa.gtf" --family NPIP --exons-json "$FROZEN/npip_read_pool.json" \
    --arm "DEF=$PD.fam.loci.gff3,$PD.fam.clusters.tsv" --arm "PRE=$PP.fam.loci.gff3,$PP.fam.clusters.tsv" --out "$W/nodes.json"
  for seed in 0 1; do
    PYTHONHASHSEED=$seed light python3 "$REPO/bench/copy_support.py" --copies "$ANN/copies.hsa.tsv" --truth "$ANN/truth.hsa.gtf" --bam "$BAM" --family NPIP \
      --out "$W/support_s$seed" --loci "DEF=$PD.fam.loci.gff3,$PD.families.gtf" --loci "PRE=$PP.fam.loci.gff3,$PP.gtf" --nodes "$W/nodes.json" > "$W/support_s$seed.stdout"
  done
  python3 - "$W" <<'PY'
import json, sys
w = sys.argv[1]
a, b = (json.load(open(f"{w}/support_s{s}.json")) for s in (0, 1))
print("scorer deterministic under PYTHONHASHSEED 0 vs 1:", a == b)
if a != b:
    sys.exit("copy_support.py summary differs between the two seeds")
PY
  cp "$W/support_s0.json" "$W/support.json"; cp "$W/support_s0.copies.tsv" "$W/support.copies.tsv"
  python3 "$HERE/score_families.py" --fs "$BIN/family_score" --clusters "$PD.fam.clusters.tsv" --contig chr16 --label DEF --out "$W/fs_DEF.json" --work "$W/fs"
  python3 "$HERE/score_families.py" --fs "$BIN/family_score" --clusters "$PP.fam.clusters.tsv" --contig chr16 --label PRE --out "$W/fs_PRE.json" --work "$W/fs"
  { echo "provenance: HEAD products vs the e163d955 products of the container-headroom run (cmp; 'identical' or the first difference)"
    for pair in "DEF gtf|$PD.gtf|$CH/data/human_A119b/asm_f1v2/chr16/human_A119b.chr16.gtf" "DEF families.gtf|$PD.families.gtf|$CH/data/human_A119b/asm_f1v2/chr16/human_A119b.chr16.families.gtf" \
                "DEF clusters|$PD.fam.clusters.tsv|$CH/data/human_A119b/fam/chr16/D.fam.clusters.tsv" "DEF loci.gff3|$PD.fam.loci.gff3|$CH/data/human_A119b/fam/chr16/D.fam.loci.gff3" \
                "DEF copies|$PD.fam.copies.tsv|$CH/data/human_A119b/fam/chr16/D.fam.copies.tsv" \
                "PRE gtf|$PP.gtf|$CH/data/human_A119b/asm_off/chr16/human_A119b.chr16.gtf" "PRE clusters|$PP.fam.clusters.tsv|$CH/data/human_A119b/fam/chr16/B0.fam.clusters.tsv" \
                "PRE loci.gff3|$PP.fam.loci.gff3|$CH/data/human_A119b/fam/chr16/B0.fam.loci.gff3"; do
      IFS='|' read -r label x y <<< "$pair"
      if cmp -s "$x" "$y"; then echo "  $label: identical"; else echo "  $label: DIFFERS ($(diff <(cut -f1-9 "$x") <(cut -f1-9 "$y") | grep -c '^[<>]') differing lines; $(cmp "$x" "$y" 2>&1 | head -1))"; fi
    done; } | tee "$W/provenance.txt" ;;
verdict)
  python3 "$HERE/verdict.py" --work "$W" ;;
*) echo "usage: run.sh gates | asm ARM | fam ARM | score | verdict" >&2; exit 2 ;;
esac
