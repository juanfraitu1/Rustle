#!/bin/bash
# rustle_pipeline.sh — the whole pipeline, one command per stage, shipped defaults (2026-09-23).
#
#   assemble  reads -> loci -> isoform GTF                              copy_assign --assemble-only --genome-wide
#   families  gene families from the de novo loci (all-vs-all -> MCL)   mcl_families --from-gtf --emit-units
#             = THE default de novo family definition (user decision 2026-09-25): one copy per member locus =
#             its representative transcript (PREFIX.fam.copies.tsv/.fa, the copy table copy assignment consumes)
#   catalog   LEGACY copy catalog (gw_family_catalog; kept, not the default definition)
#   assign    per-read copy assignment on the catalog (assign/abstain)  copy_assign --families
#   flag      copies the reference does not contain, from RNA alone     missing_copy_flag --scan-only / --from-scan
#             (+ optional DNA confirmation against --confirm genomes)
#   all       every stage in order
# (thesis record: families/catalog = O1, assign = O2, flag = O3)
#
# usage: tools/rustle_pipeline.sh STAGE --bam B --fasta G --out PREFIX [--index G.splice.mmi] [--gff ANNOT.gff]
#        [--confirm NAME=X.mmi ...] [--foreign NAME=X.mmi ...] [--threads N] [--bin DIR] [--no-seed-secondaries]
#        [--no-cache] [--inspect] [--piecewise [--max-pieces N] [--budget-s S] [--piece-records R] [--piece LABEL]]
#   tools/rustle_pipeline.sh cache-ls --out PREFIX      list what PREFIX.cache holds (cache-clear: delete it)
# --piecewise (catalog only; needs the cache): the catalog's representatives (both BAM passes + the span-overlap
#   collapse and its POA, none of which compares two contigs) are built one CONTIG per piece, each cached in
#   PREFIX.cache/reps/; when all are cached they are merged into the one-run order and the catalog continues as one
#   run (the k11 all-vs-all of ALL representatives, edges, families; families still cross contigs). Same products
#   as without it (cmp-checked 2026-09-25). --max-pieces N / --budget-s S bound one call: it computes that much,
#   appends to PREFIX.catalog.log and EXITS 75 while work remains (pieces pending, or the merge still to do) — call
#   it again until it exits 0. --piece-records R also cuts a contig of more than R BAM records into pieces of ~R
#   at positions no record crosses (still exact; human chr1 alone exceeds 10 min as one piece); --piece LABEL
#   computes that one piece only (as the log labels it, e.g. chr13:65273-15760749), for a known-heavy piece. The all-vs-all
#   itself can be sharded through RUSTLE_MINIMAP2 (tools/mm2_shard.sh).
# CACHE (default on, --no-cache turns it off): families and catalog keep their expensive intermediates in
#   PREFIX.cache/ (RUSTLE_CACHE_DIR): the catalog's collapsed representatives (reps/<key>/reps.tsv + reps.fa, both
#   BAM passes and the locus collapse) and every all-vs-all PAF (paf/<key>/out.paf). A re-run that changes only
#   downstream settings (the E_r edge rule, gamma, coverage split) replays them: human chr16 358 s -> 0.9 s,
#   byte-identical outputs. Keys cover the binary, BAM/FASTA (+ indexes) and every upstream RUSTLE_* setting.
#   The families PAF is keyed on every byte of PREFIX.fam.loci.fa (hashed as it is written) and replayed as a HARD
#   LINK to PREFIX.fam.loci.paf, not a copy (chimp: the hit step 2.6-13 s -> 0.01 s); a write through that link
#   invalidates the entry (mtime/inode/sampled-content pins), and RUSTLE_CACHE_VERIFY=1 re-hashes it on every hit.
# --inspect: also write the analyst dumps — catalog edge tables (PREFIX.cache/inspect/catalog.*: reps.fa, PAF,
#   edges.tsv, nodes.tsv, rule.tsv, params.tsv), per-round collapse statistics in the catalog log, and the
#   assignment evidence (PREFIX.assign.psv_*.tsv, PREFIX.assign.posterior.tsv).
# Loci are seeded from primaries PLUS secondaries within 2% of the molecule's GENOME-WIDE best alignment score
#   (one `as_table` pass over the BAM -> PREFIX.molecules.tsv + its binary load sidecar PREFIX.molecules.tsv.asbin,
#   reused if present): on gorilla NC_073244.2 (the
#   seeding pre-registration's verdict contig) this finds 30 more loci (97% annotated) and joins more
#   >=90%-identity referee pairs (21 -> 26 of 83 on genes with >= 2 exonic primary reads; all of the gain is one
#   tandem array) at pair precision 1.000 (rows 1060/1100/1101/1116); the streaming assembler applies the filter
#   at no memory cost (validated identical to the buffered path). --no-seed-secondaries restores primaries-only seeding (2026-09-24 default flip).
#   The table is only as genome-wide as the BAM: on a region SLICE it admits secondaries whose real best lies
#   outside the slice (tes44: 4,210 transcripts vs 3,911 with the full-BAM table) — give the driver the full BAM.
# The splice index is needed by `flag` (home search); the annotation (--gff) by `flag` only (IG/TR screen).
# Environment: RUSTLE_POLISH_SUBCHAIN=tag|drop adds `--polish-subchain` to `assemble` (default unset = off);
#   RUSTLE_POLISH_TSS=tag|rescue|split adds `--polish-tss` to `assemble` (default unset = off);
#   RUSTLE_POLISH_TES=tag|pas-end adds `--polish-tes` to `assemble` (default unset = off);
#   RUSTLE_POLISH_JUNCTION_SNAP=equiv|reads adds `--polish-junction-snap` to `assemble` (default unset = off);
#   RUSTLE_GTF_REGROUP=1 adds `--gtf-regroup` to `assemble` (RG3 regroup after polish; default unset = off);
#   RUSTLE_FAMILY_CONTAINER=1 adds `--emit-container` to `families` (PREFIX.fam.container*.tsv; default unset = off).
# `families` is the DE NOVO mode (loci from the assembled GTF). The GUIDED mode (loci = the annotation's gene and
# pseudogene bodies, PREREG_heldout_families_2026-09-20 §2) is not a driver stage: figures/_o1_recovery.py
# (guided_families) runs its recipe step by step. Every product carries the PREFIX.
set -euo pipefail
STAGE=${1:-all}; shift || true
BAM=""; FASTA=""; OUT=""; INDEX=""; GFF=""; THREADS=4; BIN="$(dirname "$0")/../target/release"; CONFIRM=(); FOREIGN=(); SEED_SEC=1; CACHE=1; INSPECT=0
PIECEWISE=0; MAX_PIECES=0; BUDGET_S=0; PIECE_RECORDS=0; PIECE=""
while [ $# -gt 0 ]; do
  case "$1" in
    --bam) BAM=$2; shift 2;; --fasta) FASTA=$2; shift 2;; --out) OUT=$2; shift 2;; --index) INDEX=$2; shift 2;;
    --gff) GFF=$2; shift 2;; --threads) THREADS=$2; shift 2;; --bin) BIN=$2; shift 2;;
    --confirm) CONFIRM+=(--confirm "$2"); shift 2;; --foreign) FOREIGN+=(--foreign "$2"); shift 2;;
    --seed-secondaries) SEED_SEC=1; shift;; --no-seed-secondaries) SEED_SEC=0; shift;;
    --cache) CACHE=1; shift;; --no-cache) CACHE=0; shift;; --inspect) INSPECT=1; shift;;
    --piecewise) PIECEWISE=1; shift;; --max-pieces) MAX_PIECES=$2; shift 2;; --budget-s) BUDGET_S=$2; shift 2;;
    --piece-records) PIECE_RECORDS=$2; shift 2;; --piece) PIECE=$2; shift 2;;
    *) echo "unknown argument $1" >&2; exit 2;;
  esac
done
if [ "$STAGE" = cache-clear ]; then
  [ -n "$OUT" ] || { echo "cache-clear needs --out PREFIX" >&2; exit 2; }
  [ -d "$OUT.cache" ] && du -sh "$OUT.cache" && rm -rf "$OUT.cache"
  exit 0
fi
if [ "$STAGE" = cache-ls ]; then
  [ -n "$OUT" ] || { echo "cache-ls needs --out PREFIX" >&2; exit 2; }
  for d in "$OUT".cache/*/*/; do
    [ -f "$d/DONE" ] || continue
    printf '%s\t%s\t%s\n' "$(du -sh "$d" | cut -f1)" "${d%/}" "$(head -1 "$d/key.tsv")"
  done
  exit 0
fi
[ -n "$BAM" ] && [ -n "$FASTA" ] && [ -n "$OUT" ] || { echo "need --bam, --fasta, --out" >&2; exit 2; }
export TMPDIR=${TMPDIR:-/tmp}
if [ "$CACHE" = 1 ]; then export RUSTLE_CACHE_DIR="$OUT.cache"; else unset RUSTLE_CACHE_DIR; fi
INSPECT_CAT=(); INSPECT_ASSIGN=()
if [ "$INSPECT" = 1 ]; then
  mkdir -p "$OUT.cache/inspect"
  INSPECT_CAT=("RUSTLE_ER_EDGE_DUMP=$OUT.cache/inspect/catalog" RUSTLE_COLLAPSE_STATS=1)
  INSPECT_ASSIGN=(--dump-psv --posterior)
fi
POLISH="--assembly-polish full --polish-isoform-fraction 0.02 --polish-mono-shadow --polish-mono-quantile 0.82 --polish-ism-ratio 0.7 --polish-retained-ratio 10"
# RUSTLE_POLISH_SUBCHAIN=tag|drop (opt-in; unset = off, the same command): copy_assign --polish-subchain, which tags
# (or drops) each transcript that is an end-compatible contiguous sub-chain of a longer emitted transcript of the same
# locus with >= 1/2 its reads (`subchain_of` / `subchain_missing`; drop is a documented dev-only trade, see its --help)
case "${RUSTLE_POLISH_SUBCHAIN:-}" in
  "") ;;
  off|tag|drop) POLISH="$POLISH --polish-subchain $RUSTLE_POLISH_SUBCHAIN" ;;
  *) echo "[rustle_pipeline] RUSTLE_POLISH_SUBCHAIN must be off, tag or drop (got '$RUSTLE_POLISH_SUBCHAIN')" >&2; exit 2 ;;
esac
# RUSTLE_POLISH_TSS=tag|rescue|split (opt-in; unset = off, the same command): copy_assign --polish-tss, the read-proven
# TSS (tag: `tss_clusters` only; rescue: also keeps proven short-TSS forms the polish dropped; split: also emits a chain
# once per proven TSS cluster). Dev-only, in-sample evidence: see its --help
case "${RUSTLE_POLISH_TSS:-}" in
  "") ;;
  off|tag|rescue|split) POLISH="$POLISH --polish-tss $RUSTLE_POLISH_TSS" ;;
  *) echo "[rustle_pipeline] RUSTLE_POLISH_TSS must be off, tag, rescue or split (got '$RUSTLE_POLISH_TSS')" >&2; exit 2 ;;
esac
# RUSTLE_POLISH_TES=tag|pas-end (opt-in; unset = off, the same command): copy_assign --polish-tes, the read + genome
# proven TES (tag: `tes_clusters` / `tes_pas` / `tes_primed` only; pas-end: also moves an internally primed 3' end to
# the most-3' PAS-proven cluster of the transcript's own reads). Dev-only, in-sample evidence: see its --help
case "${RUSTLE_POLISH_TES:-}" in
  "") ;;
  off|tag|pas-end) POLISH="$POLISH --polish-tes $RUSTLE_POLISH_TES" ;;
  *) echo "[rustle_pipeline] RUSTLE_POLISH_TES must be off, tag or pas-end (got '$RUSTLE_POLISH_TES')" >&2; exit 2 ;;
esac
# RUSTLE_POLISH_JUNCTION_SNAP=equiv|reads (opt-in; unset = off, the same command): copy_assign --polish-junction-snap,
# the fuzzy-junction merge that keeps read-proven splice sites (reads: a junction within 10 bp of a better-supported one
# in its locus moves onto it only when its own reads do not prove it and the other's reads do). Dev-only, in-sample
# evidence: see its --help
case "${RUSTLE_POLISH_JUNCTION_SNAP:-}" in
  "") ;;
  off|equiv|reads) POLISH="$POLISH --polish-junction-snap $RUSTLE_POLISH_JUNCTION_SNAP" ;;
  *) echo "[rustle_pipeline] RUSTLE_POLISH_JUNCTION_SNAP must be off, equiv or reads (got '$RUSTLE_POLISH_JUNCTION_SNAP')" >&2; exit 2 ;;
esac
# RUSTLE_GTF_REGROUP=1 (opt-in; unset or 0 = off, the same command): copy_assign --gtf-regroup, RG3 — after every polish
# step, split a gene_id whose surviving transcripts share no same-strand exonic base (a dropped readthrough bridge, or
# two pre-polish components that collided on one base tid); split-only, intron chains unchanged, the deeper piece keeps
# the name and the others become <gene_id>.rg<k>. See its --help
case "${RUSTLE_GTF_REGROUP:-}" in
  ""|0) ;;
  1) POLISH="$POLISH --gtf-regroup" ;;
  *) echo "[rustle_pipeline] RUSTLE_GTF_REGROUP must be 0 or 1 (got '$RUSTLE_GTF_REGROUP')" >&2; exit 2 ;;
esac
# RUSTLE_FAMILY_CONTAINER=1 (opt-in; unset or 0 = off, the same command): mcl_families --emit-container, the container
# of each family member's extra pieces (docs/PREREG_fusion_container_sim_2026-09-28.md §1): every clustered locus's
# all-transcript exon blocks, core (aligned to an exon base of another member of its family) or accessory, and the
# other families each accessory block aligns to -> PREFIX.fam.container.tsv / .container_relations.tsv /
# .container_summary.tsv. It never changes a family. See its --help
FAM_EXTRA=()
case "${RUSTLE_FAMILY_CONTAINER:-}" in
  ""|0) ;;
  1) FAM_EXTRA=(--emit-container) ;;
  *) echo "[rustle_pipeline] RUSTLE_FAMILY_CONTAINER must be 0 or 1 (got '$RUSTLE_FAMILY_CONTAINER')" >&2; exit 2 ;;
esac
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
# families: the de novo families AND their copy table (--emit-units with --from-gtf: PREFIX.fam.copies.tsv/.fa/.regions,
# the gw_family_catalog copies contract, one copy per member locus = its representative transcript and its spliced exon
# sum). clusters.tsv / loci.* are byte-identical to a run without it (cmp-checked 2026-09-25,
# docs/PREREG_families_copy_table_2026-09-25.md); the copy table is a new product. A binary older than the copy table
# (its --help does not name <out>.copies.tsv) still writes the families, with a warning and no copy table.
stage_families() {
  say "families: gene families on the de novo loci of $OUT.gtf"
  local copies=() help
  help=$("$BIN/mcl_families" --help 2>&1 || true)
  case "$help" in
    *'<out>.copies.tsv'*) copies=(--emit-units);;
    *) say "families: WARNING $BIN/mcl_families predates the families copy table; rebuild it to write $OUT.fam.copies.tsv";;
  esac
  if [ ${#FAM_EXTRA[@]} -gt 0 ] && [[ "$help" != *'--emit-container'* ]]; then
    echo "[rustle_pipeline] RUSTLE_FAMILY_CONTAINER=1 but $BIN/mcl_families predates --emit-container; rebuild it" >&2; exit 2
  fi
  "$BIN/mcl_families" --from-gtf "$OUT.gtf" --fasta "$FASTA" --threads "$THREADS" \
    --min-exonic-bp 1 --min-shared-exon-frac 0.60 "${copies[@]}" "${FAM_EXTRA[@]}" --out "$OUT.fam" > "$OUT.families.log" 2>&1
  say "families: $(awk 'NR>1' "$OUT.fam.clusters.tsv" | cut -f1 | sort -u | wc -l) clusters ($OUT.fam.clusters.tsv)"
  if [ ${#copies[@]} -gt 0 ]; then
    say "families: $(awk 'NR>1' "$OUT.fam.copies.tsv" | wc -l) copies (locus representatives) in $(awk 'NR>1' "$OUT.fam.copies.tsv" | cut -f1 | sort -u | wc -l) families ($OUT.fam.copies.tsv)"
  fi
  if [ ${#FAM_EXTRA[@]} -gt 0 ]; then
    say "families: container $(awk -F'\t' '$1=="blocks"{b=$2} $1=="accessory_blocks"{a=$2} $1=="family_relations_directed"{r=$2} END{print b" blocks, "a" accessory, "r" directed family relations"}' "$OUT.fam.container_summary.tsv") ($OUT.fam.container.tsv)"
  fi
}
# catalog: LEGACY (2026-09-25). gw_family_catalog's copy catalog (primaries only, span-aware POA collapse, exon-sum k11
# edges, gamma-quasi-clique) was the O2 roster before the families stage wrote its own copy table; kept runnable for
# comparison and for the runs made with it. The default de novo family definition is the families stage.
stage_catalog() {
  say "catalog: gw_family_catalog on $BAM"
  if [ "$PIECEWISE" = 1 ]; then
    [ "$CACHE" = 1 ] || { echo "catalog --piecewise keeps its pieces in PREFIX.cache: drop --no-cache" >&2; exit 2; }
    local rc=0
    echo "=== [rustle_pipeline] $(date '+%F %T') catalog --piecewise --max-pieces $MAX_PIECES --budget-s $BUDGET_S --piece-records $PIECE_RECORDS" >> "$OUT.catalog.log"
    env "${INSPECT_CAT[@]}" "$BIN/gw_family_catalog" --bam "$BAM" --fasta "$FASTA" --threads "$THREADS" --out "$OUT.cat" \
      --piecewise --max-pieces "$MAX_PIECES" --budget-s "$BUDGET_S" --piece-records "$PIECE_RECORDS" ${PIECE:+--piece "$PIECE"} \
      >> "$OUT.catalog.log" 2>&1 || rc=$?
    if [ "$rc" = 75 ]; then
      say "catalog: $(grep -oE '[0-9]+ of [0-9]+ pieces cached[^;]*|all-vs-all aligner stopped with its progress kept' "$OUT.catalog.log" | tail -1); call again to continue (exit 75)"
      exit 75
    fi
    [ "$rc" = 0 ] || { say "catalog failed (exit $rc), see $OUT.catalog.log"; exit "$rc"; }
  else
    env "${INSPECT_CAT[@]}" "$BIN/gw_family_catalog" --bam "$BAM" --fasta "$FASTA" --threads "$THREADS" --out "$OUT.cat" > "$OUT.catalog.log" 2>&1
  fi
  say "catalog: $(awk 'NR>1' "$OUT.cat.copies.tsv" | wc -l) copies in $(awk 'NR>1 && $2>=2' "$OUT.cat.families.tsv" | wc -l) multi-copy families"
}
stage_assign() {
  say "assign: per-read copy assignment on $OUT.cat"
  samtools view -H "$BAM" | awk '$1=="@SQ"{sub("SN:","",$2); sub("LN:","",$3); print $2":1-"$3}' > "$OUT.regions.txt"
  "$BIN/copy_assign" --bam "$BAM" --fasta "$FASTA" --regions "$OUT.regions.txt" \
    --families "$OUT.cat.copies.tsv" --copies-fa "$OUT.cat.copies.fa" "${INSPECT_ASSIGN[@]}" --out "$OUT.assign" > "$OUT.assign.log" 2>&1
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
