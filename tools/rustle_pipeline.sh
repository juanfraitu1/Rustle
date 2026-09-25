#!/bin/bash
# rustle_pipeline.sh — the whole pipeline, one command per stage, shipped defaults (2026-09-23).
#
#   assemble  reads -> loci -> isoform GTF                              copy_assign --assemble-only --genome-wide
#   families  gene families from the de novo loci (all-vs-all -> MCL)   mcl_families --from-gtf
#   catalog   the copy catalog the assignment stage consumes            gw_family_catalog
#   assign    per-read copy assignment on the catalog (assign/abstain)  copy_assign --families
#   flag      copies the reference does not contain, from RNA alone     missing_copy_flag --scan-only / --from-scan
#             (+ optional DNA confirmation against --confirm genomes)
#   all       every stage in order
# (thesis record: families/catalog = O1, assign = O2, flag = O3)
#
# usage: tools/rustle_pipeline.sh STAGE --bam B --fasta G --out PREFIX [--index G.splice.mmi] [--gff ANNOT.gff]
#        [--confirm NAME=X.mmi ...] [--foreign NAME=X.mmi ...] [--threads N] [--bin DIR] [--no-seed-secondaries]
# Loci are seeded from primaries PLUS secondaries within 2% of the molecule's GENOME-WIDE best alignment score
#   (one `as_table` pass over the BAM -> PREFIX.molecules.tsv, reused if present): on the held-out gorilla contig
#   this finds 30 more loci (97% annotated) and doubles the >=90%-identity referee pairs joined, at pair precision
#   1.000 (register rows 1060/1100/1101); the streaming assembler applies the filter at no memory cost (validated
#   identical to the buffered path). --no-seed-secondaries restores primaries-only seeding (2026-09-24 default flip).
#   The table is only as genome-wide as the BAM: on a region SLICE it admits secondaries whose real best lies
#   outside the slice (tes44: 4,210 transcripts vs 3,911 with the full-BAM table) — give the driver the full BAM.
# The splice index is needed by `flag` (home search); the annotation by `flag` (IG/TR screen) and, when given, by
# `families` as the guided locus set instead of the de novo one. Every product carries the PREFIX.
set -euo pipefail
STAGE=${1:-all}; shift || true
BAM=""; FASTA=""; OUT=""; INDEX=""; GFF=""; THREADS=4; BIN="$(dirname "$0")/../target/release"; CONFIRM=(); FOREIGN=(); SEED_SEC=1
while [ $# -gt 0 ]; do
  case "$1" in
    --bam) BAM=$2; shift 2;; --fasta) FASTA=$2; shift 2;; --out) OUT=$2; shift 2;; --index) INDEX=$2; shift 2;;
    --gff) GFF=$2; shift 2;; --threads) THREADS=$2; shift 2;; --bin) BIN=$2; shift 2;;
    --confirm) CONFIRM+=(--confirm "$2"); shift 2;; --foreign) FOREIGN+=(--foreign "$2"); shift 2;;
    --seed-secondaries) SEED_SEC=1; shift;; --no-seed-secondaries) SEED_SEC=0; shift;;
    *) echo "unknown argument $1" >&2; exit 2;;
  esac
done
[ -n "$BAM" ] && [ -n "$FASTA" ] && [ -n "$OUT" ] || { echo "need --bam, --fasta, --out" >&2; exit 2; }
export TMPDIR=${TMPDIR:-/tmp}
POLISH="--assembly-polish full --polish-isoform-fraction 0.02 --polish-mono-shadow --polish-mono-quantile 0.82 --polish-ism-ratio 0.7 --polish-retained-ratio 10"
say() { echo "[rustle_pipeline] $(date +%H:%M:%S) $*" >&2; }

stage_assemble() {
  say "assemble: $BAM -> $OUT.gtf"
  local seed_env=()
  if [ "$SEED_SEC" = 1 ]; then
    if ! head -1 "$OUT.molecules.tsv" 2>/dev/null | grep -qF "bam=$(readlink -f "$BAM")$(printf '\t')"; then   # absent, truncated or from another BAM
      say "assemble: genome-wide best-AS table -> $OUT.molecules.tsv (one pass over the BAM)"
      "$BIN/as_table" --bam "$BAM" --out "$OUT.molecules.tsv" --threads "$THREADS" > "$OUT.as_table.log" 2>&1
    fi
    seed_env=(RUSTLE_GTF_SECONDARY=1 RUSTLE_GTF_SECONDARY_AS_RATIO=0.98 "RUSTLE_GTF_SECONDARY_AS_TABLE=$OUT.molecules.tsv")
  fi
  env "${seed_env[@]}" "$BIN/copy_assign" --assemble-only --genome-wide --assembly-junctions strict $POLISH --gtf-tpm \
    --bam "$BAM" --fasta "$FASTA" --out "$OUT" > "$OUT.assemble.log" 2>&1
  say "assemble: $(awk -F'\t' '$3=="transcript"' "$OUT.gtf" | wc -l) transcripts"
}
stage_families() {
  say "families: gene families on the de novo loci of $OUT.gtf"
  "$BIN/mcl_families" --from-gtf "$OUT.gtf" --fasta "$FASTA" --threads "$THREADS" \
    --min-exonic-bp 1 --min-shared-exon-frac 0.60 --out "$OUT.fam" > "$OUT.families.log" 2>&1
  say "families: $(awk 'NR>1' "$OUT.fam.clusters.tsv" | cut -f1 | sort -u | wc -l) clusters ($OUT.fam.clusters.tsv)"
}
stage_catalog() {
  say "catalog: gw_family_catalog on $BAM"
  "$BIN/gw_family_catalog" --bam "$BAM" --fasta "$FASTA" --threads "$THREADS" --out "$OUT.cat" > "$OUT.catalog.log" 2>&1
  say "catalog: $(awk 'NR>1' "$OUT.cat.copies.tsv" | wc -l) copies in $(awk 'NR>1 && $2>=2' "$OUT.cat.families.tsv" | wc -l) multi-copy families"
}
stage_assign() {
  say "assign: per-read copy assignment on $OUT.cat"
  samtools view -H "$BAM" | awk '$1=="@SQ"{sub("SN:","",$2); sub("LN:","",$3); print $2":1-"$3}' > "$OUT.regions.txt"
  "$BIN/copy_assign" --bam "$BAM" --fasta "$FASTA" --regions "$OUT.regions.txt" \
    --families "$OUT.cat.copies.tsv" --copies-fa "$OUT.cat.copies.fa" --out "$OUT.assign" > "$OUT.assign.log" 2>&1
  say "assign: $(awk -F'\t' 'NR>1 && $4=="assigned"' "$OUT.assign.assignments.tsv" | wc -l) assigned rows of $(awk 'NR>1' "$OUT.assign.assignments.tsv" | wc -l) (one row per read x family)"
}
stage_flag() {
  [ -n "$INDEX" ] || { echo "flag needs --index (splice .mmi of the primary genome)" >&2; exit 2; }
  local LOCI=${GFF:-$OUT.gtf}
  say "flag: scan $BAM on $LOCI"
  "$BIN/missing_copy_flag" --bam "$BAM" --fasta "$FASTA" --loci "$LOCI" ${GFF:+--gff "$GFF"} --index x --threads "$THREADS" \
    --out "$OUT.flag_scan" --scan-only > "$OUT.flag_scan.log" 2>&1
  say "flag: align + verdict"
  "$BIN/missing_copy_flag" --bam "$BAM" --fasta "$FASTA" --loci "$LOCI" --index "$INDEX" --threads "$THREADS" \
    "${CONFIRM[@]}" "${FOREIGN[@]}" --out "$OUT.flag" --from-scan "$OUT.flag_scan" > "$OUT.flag.log" 2>&1
  say "flag: $(grep -o 'verdicts: .*' "$OUT.flag.log")"
}
case "$STAGE" in
  assemble) stage_assemble;; families) stage_families;; catalog) stage_catalog;; assign) stage_assign;; flag) stage_flag;;
  all) stage_assemble; stage_families; stage_catalog; stage_assign; stage_flag;;
  *) echo "unknown stage $STAGE" >&2; exit 2;;
esac
say "done"
