#!/bin/bash
# bench/rep_rule/run.sh — the runner of docs/archive/2026-10/PREREG_locus_representative_rule_2026-10-04.md: the de novo locus representative
# R_M (most-reads, the shipped rule) vs R_J (most-junctions, `mcl_families --representative most-junctions`).
#
#   bench/rep_rule/run.sh families SPECIES CONTIG ARM   ARM = R_M | R_J: the driver's `families` stage on the contig's Figure 7
#                                                       de novo GTF, ONE heavy call (tools/rlock.sh heavy, /usr/bin/time -v)
#   bench/rep_rule/run.sh h3 SPECIES CONTIG             H3: copies = the contig's annotated protein-coding genes (score.py
#                                                       h3-inputs: human CAT/Liftoff v2.0, gorilla RefSeq), bench/copy_support.py
#                                                       on both arms' loci in one call (rlock heavy: a whole contig's BAM)
#   bench/rep_rule/run.sh npip                          H3 on human chr16: the 25 NPIP copies of
#                                                       docs/archive/2026-10/SPLICED_COPY_SUPPORT_2026-10-04.md (same copies and truth), both arms
# SPECIES = human | gorilla. Scoring (H1/H2/H3 tables) and the registered decision: bench/rep_rule/score.py.
#
# Inputs, read-only: ${FIG7}/<species>_<contig>.denovo.gtf (the genome-wide de novo assembly restricted to the contig, the
# Figure 7 development tables of 2026-09-25), copied once to ${W}/<species>_<contig>/; each arm's PREFIX is
# ${W}/<species>_<contig>/<ARM> and the stage reads PREFIX.gtf (a copy of it). BAM / FASTA from figures/inputs.local.tsv.
# Environment (defaults = the run of 2026-10-04): REP_WORK (products), REP_FIG7 (the read-only Figure 7 cache), REP_NPIP_ANN
# (the NPIP copies + truth of `npip`), REP_CAT_GFF (read by score.py h3-inputs); all exported to score.py with REP_BIN.
# REP_BIN (default the m2 release dir) holds the ONE binary both arms use; every log records the sha1 of its mcl_families and
# family_score before and after the call (score.py refuses a contig whose arms differ).
# Arms: R_M runs with RUSTLE_REPRESENTATIVE unset (`env -u`), R_J with RUSTLE_REPRESENTATIVE=most-junctions on that one command.
# ⚠ RUSTLE_BRIDGE_REGROUP=off on BOTH arms: the Figure 7 GTFs were assembled on 2026-09-25, before f1v2 became the default
# (2026-09-29); they hold no bridge relation (0 `fusion_of`) and no PREFIX.families.gtf, so the driver's guard refuses its
# default f1v2 on them ("set RUSTLE_BRIDGE_REGROUP=off ... on a GTF assembled with off"). Every other families setting is the
# stage's shipped default (--min-exonic-bp 1 --min-shared-exon-frac 0.60 --emit-units; --min-cov-shorter = the binary's 0.70).
set -euo pipefail
HERE=$(cd "$(dirname "$0")" && pwd); REPO=$(cd "$HERE/../.." && pwd)
BIN=${REP_BIN:-/mnt/linuxdisk/home/juanfraitu/rustle_target_m2/release}
W=${REP_WORK:-/mnt/linuxdisk/tmp/rep_rule}
FIG7=${REP_FIG7:-/mnt/linuxdisk/tmp/rustle_figures/fig7/current}
NPIP_ANN=${REP_NPIP_ANN:-/mnt/linuxdisk/tmp/rustle_figures_dev/copy_recovery_tools_cat/ann}
export REP_BIN=$BIN REP_WORK=$W REP_FIG7=$FIG7   # score.py reads the same values
THREADS=4
export TMPDIR=$W/tmp
mkdir -p "$TMPDIR"
cfgval() { (cd "$REPO/figures" && python3 -c 'import figlib, sys; print(figlib.load_inputs()[sys.argv[1]])' "$1"); }
sha() { sha1sum "$1" | cut -d' ' -f1; }
stamp() {   # binary provenance, one line each (score.py reads them)
  echo "mcl_families_sha1	$(sha "$BIN/mcl_families")"
  echo "family_score_sha1	$(sha "$BIN/family_score")"
}
# no RUSTLE_* knob may leak into an arm from the calling shell (R_J's is set on its one command below)
leak=$(env | grep '^RUSTLE_' || true)
[ -z "$leak" ] || { echo "[rep_rule] refusing: RUSTLE_* set in the calling environment: $leak" >&2; exit 2; }

cmd=${1:?families|h3|npip}; shift
case "$cmd" in
families)
  SP=${1:?species}; C=${2:?contig}; ARM=${3:?R_M|R_J}
  case "$ARM" in R_M) REP=(env -u RUSTLE_REPRESENTATIVE);; R_J) REP=(env RUSTLE_REPRESENTATIVE=most-junctions);;
    *) echo "ARM must be R_M or R_J" >&2; exit 2;; esac
  D=$W/${SP}_${C}; mkdir -p "$D"
  src=$FIG7/${SP}_${C}.denovo.gtf; gtf=$D/${SP}_${C}.denovo.gtf
  [ -s "$gtf" ] || cp "$src" "$gtf"
  cmp -s "$src" "$gtf" || { echo "[rep_rule] $gtf differs from $src" >&2; exit 2; }
  cp "$gtf" "$D/$ARM.gtf"
  rm -f "$D/$ARM".fam.* "$D/$ARM".families.log
  log=$D/$ARM.run.log
  {
    echo "date	$(date -Is)"; echo "arm	$ARM"; echo "species	$SP"; echo "contig	$C"
    echo "rule	$([ "$ARM" = R_J ] && echo 'RUSTLE_REPRESENTATIVE=most-junctions' || echo 'RUSTLE_REPRESENTATIVE unset (most-reads)')"
    echo "bridge_regroup	off (the input GTF predates f1v2: 0 fusion_of, no PREFIX.families.gtf)"
    echo "gtf_source	$src"; echo "gtf_sha1	$(sha "$gtf")"; echo "bin	$BIN"; stamp
    echo "--- command"
  } > "$log"
  BAM=$(cfgval "${SP}_bam"); FASTA=$(cfgval "${SP}_fasta")
  set +e
  "${REP[@]}" RUSTLE_BRIDGE_REGROUP=off bash "$REPO/tools/rlock.sh" heavy /usr/bin/time -v \
    bash "$REPO/tools/rustle_pipeline.sh" families --bam "$BAM" --fasta "$FASTA" --out "$D/$ARM" --bin "$BIN" \
    --threads "$THREADS" >> "$log" 2>&1
  rc=$?
  set -e
  { echo "--- after"; echo "exit	$rc"; stamp | sed 's/^/after_/'; } >> "$log"
  # the representative row: present (last line) in R_J's params.tsv, absent in R_M's
  if [ "$rc" = 0 ]; then
    last=$(tail -1 "$D/$ARM.fam.params.tsv")
    if [ "$ARM" = R_J ]; then
      [ "$last" = "$(printf 'representative\tmost-junctions')" ] || { echo "[rep_rule] R_J params.tsv lacks the representative row" | tee -a "$log" >&2; exit 3; }
    else
      ! grep -q '^representative' "$D/$ARM.fam.params.tsv" || { echo "[rep_rule] R_M params.tsv has a representative row" | tee -a "$log" >&2; exit 3; }
    fi
    echo "params_representative_row	$([ "$ARM" = R_J ] && echo present || echo absent)" >> "$log"
  fi
  grep -E 'Elapsed \(wall clock\)|Maximum resident|Exit status' "$log" | sed 's/^\s*//'
  echo "[rep_rule] $SP $C $ARM exit $rc ($log)"
  exit "$rc" ;;
h3)
  SP=${1:?species}; C=${2:?contig}
  D=$W/${SP}_${C}; gtf=$D/${SP}_${C}.denovo.gtf
  for a in R_M R_J; do [ -s "$D/$a.fam.loci.gff3" ] || { echo "[rep_rule] $D/$a.fam.loci.gff3 missing: run the families stage" >&2; exit 2; }; done
  python3 "$HERE/score.py" h3-inputs --species "$SP" --contig "$C" --out "$D/h3"
  BAM=$(cfgval "${SP}_bam")
  log=$D/h3.run.log
  { echo "date	$(date -Is)"; echo "species	$SP"; echo "contig	$C"; echo "bam	$BAM"; stamp; echo "copy_support_sha1	$(sha "$REPO/bench/copy_support.py")"; echo "--- command"; } > "$log"
  set +e
  bash "$REPO/tools/rlock.sh" heavy /usr/bin/time -v python3 "$REPO/bench/copy_support.py" --copies "$D/h3.copies.tsv" \
    --truth "$D/h3.truth.gtf" --bam "$BAM" --family ALL --out "$D/h3.support" \
    --loci "R_M=$D/R_M.fam.loci.gff3,$gtf" --loci "R_J=$D/R_J.fam.loci.gff3,$gtf" > "$D/h3.support.stdout" 2>> "$log"
  rc=$?
  set -e
  echo "exit	$rc" >> "$log"
  grep -E 'Elapsed \(wall clock\)|Maximum resident|Exit status' "$log" | sed 's/^\s*//'
  echo "[rep_rule] h3 $SP $C exit $rc ($log)"
  exit "$rc" ;;
npip)
  D=$W/human_chr16; gtf=$D/human_chr16.denovo.gtf
  ANN=$NPIP_ANN
  BAM=$(cfgval human_bam)
  log=$D/npip.run.log
  { echo "date	$(date -Is)"; echo "copies	$ANN/copies.hsa.tsv	$(sha "$ANN/copies.hsa.tsv")"; echo "truth	$ANN/truth.hsa.gtf	$(sha "$ANN/truth.hsa.gtf")"; stamp; echo "--- command"; } > "$log"
  set +e
  bash "$REPO/tools/rlock.sh" light /usr/bin/time -v python3 "$REPO/bench/copy_support.py" --copies "$ANN/copies.hsa.tsv" \
    --truth "$ANN/truth.hsa.gtf" --bam "$BAM" --family NPIP --out "$D/npip.support" \
    --loci "R_M=$D/R_M.fam.loci.gff3,$gtf" --loci "R_J=$D/R_J.fam.loci.gff3,$gtf" > "$D/npip.support.stdout" 2>> "$log"
  rc=$?
  set -e
  echo "exit	$rc" >> "$log"
  grep -E 'Elapsed \(wall clock\)|Maximum resident|Exit status' "$log" | sed 's/^\s*//'
  echo "[rep_rule] npip exit $rc ($log)"
  exit "$rc" ;;
*) echo "usage: run.sh families SPECIES CONTIG ARM | h3 SPECIES CONTIG | npip" >&2; exit 2 ;;
esac
