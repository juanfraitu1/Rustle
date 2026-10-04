#!/bin/bash
# rustle_pipeline.sh — the whole pipeline, one command per stage, shipped defaults (2026-09-23).
#
#   assemble  reads -> loci -> isoform GTF                              copy_assign --assemble-only --genome-wide
#   families  gene families from the de novo loci (all-vs-all -> MCL)   mcl_families --from-gtf --emit-units
#             = THE default de novo family definition (user decision 2026-09-25): one copy per member locus =
#             its representative transcript (PREFIX.fam.copies.tsv/.fa, the copy table copy assignment consumes)
#   candidates  reference-absent copies from each family's own reads    o3_candidates + tools/o3_augment.py + minimap2
#             OPT-IN (ruling R14, 2026-10-02: its pre-registered acceptance, Amendment 12, FAILED —
#             docs/O3_CANDIDATES_ACCEPTANCE_2026-10-02.md; 2026-10-03: the re-run, Amendment 13, passed, but the no-deletion
#             control, Amendment 14, failed, so the default flip was reverted — docs/O3_CANDIDATES_CONTROL_A14_2026-10-03.md):
#             naming the stage runs it, `all` runs it only with --candidates.
#             (read net -> read clusters -> consensus -> flag/link/merge -> one exon-union contig per flagged
#             candidate, PREFIX.cand.*), then the augmentation (PREFIX.aug.fa = genome + the contigs, PREFIX.aug.copies.*
#             = the copy table + one row per candidate) and the patch realignment of the candidate families' reads to
#             PREFIX.aug.fa (PREFIX.aug.bam): a read whose best alignment is on a candidate contig is PLACED there by the
#             realignment (ruling R13). --delta D (0.00958), --cand-max-reads N (1000). Cost: PREFIX.aug.fa is a copy of
#             the whole genome, and the realignment indexes it on every run (about 400 s and 19 GB on a human genome).
#             The stage's own cost on a full BAM (spec §9b, R23): on the 23-GB gorilla fibroblast BAM one batch of 50 families
#             did not finish in a 10-minute call (nets phase 397-469 s: the whole-BAM sweep ~150 s, the attribution alignment
#             against the batch's 2.6 GB of net reads 213-250 s; peak RSS 10.7 GB at the stop).
#   catalog   LEGACY copy catalog (gw_family_catalog; kept, not the default definition)
#   assign    per-read copy assignment (assign/abstain) on the families'  copy_assign --families
#             copy table PREFIX.fam.copies.*; with --candidates and flagged candidates, two runs (the candidate families
#             on PREFIX.aug.*, each assigned as a family that includes its candidate copies; the others on the originals)
#             concatenated into PREFIX.assign.*. O2 assigns AS-tied molecules only (the gate is unchanged): a read the
#             realignment places uniquely on a candidate has no row (R13). Without --candidates the candidates stage's
#             products are not used. --legacy-catalog: on the legacy catalog PREFIX.cat.* instead (the stage as it ran
#             before 2026-10-02; refused together with --candidates)
#   flag      copies the reference does not contain, from RNA alone     missing_copy_flag --scan-only / --from-scan
#             (+ optional DNA confirmation against --confirm genomes; with --candidates, + the `o3_candidate` column
#             naming the flagged candidate whose nearest locus is the row's)
#   all       assemble, families, assign, flag; with --candidates the candidates stage before assign, with
#             --legacy-catalog the catalog there instead
# (thesis record: families/catalog = O1, candidates = O3, assign = O2, flag = O3; in order O1 families -> O3 candidates ->
#  O2 assign -> O3 flag)
#
# usage: tools/rustle_pipeline.sh STAGE --bam B --fasta G --out PREFIX [--index G.splice.mmi] [--gff ANNOT.gff]
#        [--confirm NAME=X.mmi ...] [--foreign NAME=X.mmi ...] [--threads N] [--bin DIR] [--no-seed-secondaries]
#        [--no-cache] [--inspect] [--piecewise [--max-pieces N] [--budget-s S] [--piece-records R] [--piece LABEL]]
#        [--candidates | --no-candidates] [--delta D] [--cand-max-reads N] [--legacy-catalog]
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
# The splice index is needed by `candidates` (the consensus sequences' genome hits) and `flag` (home search); the
# annotation (--gff) by `flag` only (IG/TR screen).
# Environment: RUSTLE_POLISH_SUBCHAIN=tag|drop adds `--polish-subchain` to `assemble` (default unset = off);
#   RUSTLE_POLISH_TSS=tag|rescue|split adds `--polish-tss` to `assemble` (default unset = off);
#   RUSTLE_POLISH_TES=tag|pas-end adds `--polish-tes` to `assemble` (default unset = off);
#   RUSTLE_POLISH_JUNCTION_SNAP=equiv|reads adds `--polish-junction-snap` to `assemble` (default unset = off);
#   RUSTLE_GTF_REGROUP=1 adds `--gtf-regroup` to `assemble` (RG3 regroup after polish; default unset = off; needs
#     RUSTLE_BRIDGE_REGROUP=off, whose default contains RG3's split);
#   RUSTLE_BRIDGE_REGROUP=off|f1|f1v2 sets `--bridge-regroup` on `assemble` (unset = f1v2, THE DEFAULT since 2026-09-29;
#     off = the 2026-09-25 products, byte for byte). With f1 or f1v2, `families`, and `flag` without --gff, read
#     PREFIX.families.gtf, the assembly without its bridge transcripts (bridges are relations, not loci, for every
#     downstream stage); every stage derives the mode the same way (BRIDGE_MODE below);
#   RUSTLE_MIN_COV_SHORTER=C sets `mcl_families --min-cov-shorter C` on `families` (unset = the binary's own default,
#     0.70 since 2026-09-29, nothing passed; 0 = the 2026-09-25 edge weights);
#   RUSTLE_FAMILY_CONTAINER=1 adds `--emit-container` to `families` (PREFIX.fam.container*.tsv; default unset = off).
# `families` is the DE NOVO mode (loci from the assembled GTF). The GUIDED mode (loci = the annotation's gene and
# pseudogene bodies, PREREG_heldout_families_2026-09-20 §2) is not a driver stage: figures/_o1_recovery.py
# (guided_families) runs its recipe step by step. Every product carries the PREFIX.
set -euo pipefail
STAGE=${1:-all}; shift || true
BAM=""; FASTA=""; OUT=""; INDEX=""; GFF=""; THREADS=4; BIN="$(dirname "$0")/../target/release"; CONFIRM=(); FOREIGN=(); SEED_SEC=1; CACHE=1; INSPECT=0
PIECEWISE=0; MAX_PIECES=0; BUDGET_S=0; PIECE_RECORDS=0; PIECE=""
LEGACY_CATALOG=0; CANDIDATES=0; DELTA=0.00958; CAND_MAX=1000
while [ $# -gt 0 ]; do
  case "$1" in
    --bam) BAM=$2; shift 2;; --fasta) FASTA=$2; shift 2;; --out) OUT=$2; shift 2;; --index) INDEX=$2; shift 2;;
    --gff) GFF=$2; shift 2;; --threads) THREADS=$2; shift 2;; --bin) BIN=$2; shift 2;;
    --confirm) CONFIRM+=(--confirm "$2"); shift 2;; --foreign) FOREIGN+=(--foreign "$2"); shift 2;;
    --seed-secondaries) SEED_SEC=1; shift;; --no-seed-secondaries) SEED_SEC=0; shift;;
    --cache) CACHE=1; shift;; --no-cache) CACHE=0; shift;; --inspect) INSPECT=1; shift;;
    --piecewise) PIECEWISE=1; shift;; --max-pieces) MAX_PIECES=$2; shift 2;; --budget-s) BUDGET_S=$2; shift 2;;
    --piece-records) PIECE_RECORDS=$2; shift 2;; --piece) PIECE=$2; shift 2;;
    --legacy-catalog) LEGACY_CATALOG=1; shift;; --candidates) CANDIDATES=1; shift;; --no-candidates) CANDIDATES=0; shift;;
    --delta) DELTA=$2; shift 2;; --cand-max-reads) CAND_MAX=$2; shift 2;;
    *) echo "unknown argument $1" >&2; exit 2;;
  esac
done
# --candidates and --legacy-catalog name two different copy tables for assign: refused together, whatever the stage. Written
# as a `case "$STAGE"` block, as the final dispatch is, so figures/samples.py's driver_stage_code leaves it out of every
# stage's code hash: it refuses a call and changes no product.
case "$STAGE" in
  *) if [ "$CANDIDATES" = 1 ] && [ "$LEGACY_CATALOG" = 1 ]; then
       echo "[rustle_pipeline] --candidates and --legacy-catalog exclude each other (assign reads the families' copy table with its candidates, or the legacy catalog): drop one" >&2; exit 2
     fi;;
esac
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
# the name and the others become <gene_id>.rg<k>. Exclusive with the bridge regroup below, whose default (f1v2) contains
# this split: RUSTLE_GTF_REGROUP=1 needs RUSTLE_BRIDGE_REGROUP=off. See its --help
case "${RUSTLE_GTF_REGROUP:-}" in
  ""|0) ;;
  1) POLISH="$POLISH --gtf-regroup" ;;
  *) echo "[rustle_pipeline] RUSTLE_GTF_REGROUP must be 0 or 1 (got '$RUSTLE_GTF_REGROUP')" >&2; exit 2 ;;
esac
# RUSTLE_BRIDGE_REGROUP=off|f1|f1v2 (unset = f1v2, THE DEFAULT since 2026-09-29 by the user's decision on the F1v2 held-out
# Outcome and its family-level side result, docs/PREREG_o1_cover_growth_2026-09-29.md; off = the 2026-09-25 products,
# byte for byte): copy_assign --bridge-regroup, passed explicitly by `assemble`. A transcript that is the only link
# between two pieces of its gene_id, when reads end at a PAS inside its intron and other reads start there at their own
# promoter (f1v2: and it carries fewer reads than each piece), becomes a relation `<gene_id>.fus<k>` with `fusion_of`,
# and the pieces split as RUSTLE_GTF_REGROUP splits them (so the two are exclusive). With f1 or f1v2, `assemble` also
# writes PREFIX.families.gtf, the GTF without the bridges, and `families` reads it, as the held-out runs did
# (docs/PREREG_f1_bridge_locus_2026-09-28.md, docs/PREREG_f1v2_readshare_2026-09-29.md), as does `flag`'s scan of the
# de novo loci (without --gff): a bridge is a relation, never a locus to scan. See its --help.
# BRIDGE_MODE is the one derivation every stage uses (assemble's flag, families' and flag's input, the guard below).
BRIDGE_MODE=${RUSTLE_BRIDGE_REGROUP:-f1v2}
FAM_GTF=$OUT.gtf
case "$BRIDGE_MODE" in
  off) ;;
  f1|f1v2)
    [ "${RUSTLE_GTF_REGROUP:-0}" = 0 ] || { echo "[rustle_pipeline] RUSTLE_BRIDGE_REGROUP=$BRIDGE_MODE (unset = f1v2) already splits every gene_id as RUSTLE_GTF_REGROUP does: set RUSTLE_BRIDGE_REGROUP=off with RUSTLE_GTF_REGROUP=1" >&2; exit 2; }
    FAM_GTF=$OUT.families.gtf ;;
  *) echo "[rustle_pipeline] RUSTLE_BRIDGE_REGROUP must be off, f1 or f1v2 (got '$RUSTLE_BRIDGE_REGROUP')" >&2; exit 2 ;;
esac
POLISH="$POLISH --bridge-regroup $BRIDGE_MODE"
# RUSTLE_MIN_COV_SHORTER=C (unset = nothing passed: the binary's default, 0.70 since 2026-09-29 by the user's decision,
# register 1006/1014; 0 = the 2026-09-25 edge weights, byte-identical to every catalog built before the flip):
# mcl_families --min-cov-shorter C, the §6x4 containment escape — a pair whose coverage of the LONGER locus fails also
# passes at coverage >= C of the SHORTER locus's exonic length, with that coverage as its edge weight. Known regressions
# (register 1007/1009): NPIP in GUIDED mode, Soto F .833 -> .800; semi-guided SD-region nodes, precision .973 -> .833.
FAM_COV=()
case "${RUSTLE_MIN_COV_SHORTER:-}" in
  "") ;;
  *) [[ "$RUSTLE_MIN_COV_SHORTER" =~ ^(0|1|1\.0+|0?\.[0-9]+)$ ]] || { echo "[rustle_pipeline] RUSTLE_MIN_COV_SHORTER must be a number in [0, 1] (got '$RUSTLE_MIN_COV_SHORTER')" >&2; exit 2; }
     FAM_COV=(--min-cov-shorter "$RUSTLE_MIN_COV_SHORTER") ;;
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
# The de novo loci a stage reads ($FAM_GTF: PREFIX.families.gtf when BRIDGE_MODE is f1 or f1v2, else PREFIX.gtf) must
# come from the assembly on disk: with a bridge mode (the default), PREFIX.families.gtf must exist and be no older than
# PREFIX.gtf; with off, a newer PREFIX.families.gtf means the GTF was assembled with a bridge mode and its bridges would
# be read as loci. $1 = the stage.
fam_gtf_guard() {
  if [ "$FAM_GTF" != "$OUT.gtf" ]; then
    [ -s "$FAM_GTF" ] && [ ! "$FAM_GTF" -ot "$OUT.gtf" ] || { echo "[rustle_pipeline] RUSTLE_BRIDGE_REGROUP=$BRIDGE_MODE (unset = f1v2) but $FAM_GTF is missing or older than $OUT.gtf: run assemble with it, or set RUSTLE_BRIDGE_REGROUP=off for $1 on a GTF assembled with off" >&2; exit 2; }
  elif [ -e "$OUT.families.gtf" ] && [ ! "$OUT.families.gtf" -ot "$OUT.gtf" ]; then
    echo "[rustle_pipeline] $OUT.gtf was assembled with a bridge mode (its bridges are relations, not loci) but RUSTLE_BRIDGE_REGROUP=off: unset it, or set assemble's mode, for $1 too" >&2; exit 2
  fi
}

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
  if [ "$FAM_GTF" != "$OUT.gtf" ]; then
    say "assemble: $(grep -c 'fusion_of "' "$OUT.gtf") bridge transcripts kept as fusion_of relations; families input $FAM_GTF"
  fi
}
# families: the de novo families AND their copy table (--emit-units with --from-gtf: PREFIX.fam.copies.tsv/.fa/.regions,
# the gw_family_catalog copies contract, one copy per member locus = its representative transcript and its spliced exon
# sum). clusters.tsv / loci.* are byte-identical to a run without it (cmp-checked 2026-09-25,
# docs/PREREG_families_copy_table_2026-09-25.md); the copy table is a new product. A binary older than the copy table
# (its --help does not name <out>.copies.tsv) still writes the families, with a warning and no copy table.
stage_families() {
  fam_gtf_guard families
  say "families: gene families on the de novo loci of $FAM_GTF"
  local copies=() help
  help=$("$BIN/mcl_families" --help 2>&1 || true)
  case "$help" in
    *'<out>.copies.tsv'*) copies=(--emit-units);;
    *) say "families: WARNING $BIN/mcl_families predates the families copy table; rebuild it to write $OUT.fam.copies.tsv";;
  esac
  if [ ${#FAM_EXTRA[@]} -gt 0 ] && [[ "$help" != *'--emit-container'* ]]; then
    echo "[rustle_pipeline] RUSTLE_FAMILY_CONTAINER=1 but $BIN/mcl_families predates --emit-container; rebuild it" >&2; exit 2
  fi
  if [ ${#FAM_COV[@]} -gt 0 ] && [[ "$help" != *'--min-cov-shorter'* ]]; then
    echo "[rustle_pipeline] RUSTLE_MIN_COV_SHORTER is set but $BIN/mcl_families predates --min-cov-shorter; rebuild it" >&2; exit 2
  fi
  "$BIN/mcl_families" --from-gtf "$FAM_GTF" --fasta "$FASTA" --threads "$THREADS" \
    --min-exonic-bp 1 --min-shared-exon-frac 0.60 "${copies[@]}" "${FAM_EXTRA[@]}" "${FAM_COV[@]}" --out "$OUT.fam" > "$OUT.families.log" 2>&1
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
# candidates (O3; spec docs/superpowers/specs/2026-10-02-o3-candidates-design.md §4, §7; OPT-IN since 2026-10-02, ruling R14:
# its pre-registered acceptance, Amendment 12, FAILED, docs/O3_CANDIDATES_ACCEPTANCE_2026-10-02.md; the re-run, Amendment 13,
# passed, but its no-deletion control, Amendment 14, failed, so the 2026-10-03 default flip was reverted,
# docs/O3_CANDIDATES_CONTROL_A14_2026-10-03.md — naming the stage runs it, `all` runs it only with --candidates, and assign
# and flag use its products only with --candidates): o3_candidates turns each family's read net (reads with a record on its
# copies, plus the unmapped and the poorly placed reads >= 300 bp that align, map-ont, over >= 50% of their length at
# de <= 0.20 to the run's net reads or copies: Amendment 13b) into read clusters at --delta, one consensus per cluster on its
# structurally central member (Amendments 13d/13e), and candidate copies (clusters beyond delta of every reference locus,
# merged by the significance test, flagged with >= 6 reads), each represented by the exon union of its clusters ->
# PREFIX.cand.{candidates.tsv,contigs.fa,nets.fa,...}. With a flagged candidate: tools/o3_augment.py writes PREFIX.aug.{fa,copies.tsv,copies.fa,regions.txt,
# families.txt} (each candidate a contig `cand_<family>_<k>` of PREFIX.aug.fa and a `member_status candidate` row of its
# family), and the reads of the candidate families (PREFIX.cand.nets.fa, every read of each such family's net) are
# realigned to PREFIX.aug.fa with the pipeline's own minimap2 flags -> PREFIX.aug.bam, the reads `assign --candidates` gives
# those families. That realignment is the placement of the candidates' reads (ruling R13): a read whose best alignment is on
# a candidate contig lands there (the fixture: all 60 reads of the deleted copy, MAPQ 60). An earlier run's PREFIX.aug.* are
# removed first, so a run without a flagged candidate leaves none behind. PREFIX.aug.fa is a whole-genome copy that minimap2
# indexes again on every run (about 400 s and 19 GB on a human genome).
stage_candidates() {
  [ -n "$INDEX" ] || { echo "candidates needs --index (splice .mmi of the primary genome)" >&2; exit 2; }
  [ -s "$OUT.fam.copies.tsv" ] || { echo "candidates needs $OUT.fam.copies.tsv (run the families stage)" >&2; exit 2; }
  rm -f "$OUT".aug.{fa,fa.fai,copies.tsv,copies.fa,regions.txt,families.txt,bam,bam.bai,mm2.log}
  if [ "$(awk 'NR > 1 && NF' "$OUT.fam.copies.tsv" | wc -l)" = 0 ]; then
    say "candidates: $OUT.fam.copies.tsv lists no copy (no family): nothing to do"; return 0
  fi
  say "candidates: o3_candidates on $OUT.fam.copies.tsv (delta $DELTA, at most $CAND_MAX reads per family)"
  "$BIN/o3_candidates" --bam "$BAM" --fasta "$FASTA" --copies "$OUT.fam.copies.tsv" --copies-fa "$OUT.fam.copies.fa" \
    --index "$INDEX" --delta "$DELTA" --max-reads "$CAND_MAX" --threads "$THREADS" --out "$OUT.cand" > "$OUT.candidates.log" 2>&1 \
    || { local rc=$?; say "candidates: o3_candidates failed (exit $rc), see $OUT.candidates.log"; exit "$rc"; }
  local n
  n=$(awk -F'\t' 'NR>1 && $5==1' "$OUT.cand.candidates.tsv" | wc -l)
  say "candidates: $n flagged candidate copies in $(awk -F'\t' 'NR>1 && $5==1' "$OUT.cand.candidates.tsv" | cut -f1 | sort -u | wc -l) families ($OUT.cand.candidates.tsv)"
  [ "$n" -gt 0 ] || return 0
  python3 "$(dirname "$0")/o3_augment.py" --fasta "$FASTA" --copies "$OUT.fam.copies.tsv" --copies-fa "$OUT.fam.copies.fa" \
    --regions "$OUT.fam.copies.regions" --cand "$OUT.cand" --out "$OUT.aug" || exit $?
  samtools faidx "$OUT.aug.fa"
  say "candidates: realigning the $(grep -c '^>' "$OUT.cand.nets.fa") reads of the candidate families to $OUT.aug.fa"
  if ! "${RUSTLE_MINIMAP2:-minimap2}" -ax splice:hq -uf --eqx -Y -N 50 -p 0.1 --secondary=yes -t "$THREADS" "$OUT.aug.fa" "$OUT.cand.nets.fa" 2> "$OUT.aug.mm2.log" \
       | samtools sort -@ 2 -o "$OUT.aug.bam" -; then
    say "candidates: the patch realignment failed, see $OUT.aug.mm2.log"; rm -f "$OUT.aug.bam"; exit 1
  fi
  samtools index "$OUT.aug.bam"
  say "candidates: $OUT.aug.bam: $(samtools view -F 2308 "$OUT.aug.bam" | awk -F'\t' 'FILENAME == ARGV[1] { if (FNR > 1 && $5 == 1) c[$2] = 1; next } { n++; k += ($3 in c) } END { print n + 0 " primary records, " k + 0 " of them on a candidate contig" }' "$OUT.cand.candidates.tsv" -)"
}
# assign (O2): per-read copy assignment on the families' copy table (PREFIX.fam.copies.tsv/.fa, the default O1 output),
# swept over PREFIX.regions.txt (fam_regions). With --candidates and flagged candidates (PREFIX.aug.families.txt from the
# candidates stage of this families run): assign_with_candidates; candidate products present without --candidates are not
# used (one line says so). --legacy-catalog: the legacy catalog PREFIX.cat.* over whole contigs, as the stage ran before
# 2026-10-02. The three helpers are defined inside the stage (they are its code alone: figures/samples.py's driver_stage_code
# then keeps them out of the other stages' code hashes).
stage_assign() {
  # fam_regions all|only|skip [LIST]: copy_assign's --regions for the families of PREFIX.fam.copies.regions (`{fid}\t{chrom}:
  # {lo}-{hi}`, the copies' hull +- 5 kb per family and contig, as mcl_families writes it): every row, the rows of the
  # families LIST names (one id per line), or the other rows, MERGED per contig into disjoint intervals (an interval that
  # overlaps or touches the previous one joins it). copy_assign refuses overlapping --regions, and families' hulls
  # overlap and nest wherever families interleave (human chr16: 92 overlapping neighbours among 111 regions), so the
  # second column cannot be passed as it is. Every family's own interval lies inside exactly one merged region, which is
  # what binds it; its reads are still gathered around each copy (copy_assign's copy windows), not over the region.
  # tools/o3_augment.py merges the candidate families' rows by the same rule (PREFIX.aug.regions.txt).
  fam_regions() {
    local mode=$1 list=${2:-/dev/null}
    awk -F'\t' -v OFS='\t' -v mode="$mode" -v list="$list" '
      FILENAME == list { if ($1 != "") pick[$1] = 1; next }
      mode == "all" || (mode == "only") == ($1 in pick) {
        p = match($2, /:[0-9]+-[0-9]+$/)
        if (!p) { print "[rustle_pipeline] bad region \"" $2 "\" in " FILENAME > "/dev/stderr"; bad = 1; exit 2 }
        split(substr($2, p + 1), r, "-"); print substr($2, 1, p - 1), r[1], r[2] }
      END { if (bad) exit 2 }' "$list" "$OUT.fam.copies.regions" |
      LC_ALL=C sort -t "$(printf '\t')" -k1,1 -k2,2n -k3,3n |
      awk -F'\t' '$1 != c || $2 > e { if (c != "") print c ":" s "-" e; c = $1; s = $2; e = $3; next }
                  $3 > e { e = $3 }
                  END { if (c != "") print c ":" s "-" e }'
  }
  # cand_ready: true when this run uses the candidates stage's products (--candidates, PREFIX.cand.candidates.tsv present);
  # a table older than PREFIX.fam.copies.tsv was made for other families: an error, not a silent mismatch.
  cand_ready() {
    [ "$CANDIDATES" = 1 ] && [ -s "$OUT.cand.candidates.tsv" ] || return 1
    if [ "$OUT.cand.candidates.tsv" -ot "$OUT.fam.copies.tsv" ]; then
      echo "[rustle_pipeline] $OUT.cand.candidates.tsv is older than $OUT.fam.copies.tsv (families ran again): run the candidates stage again, or drop --candidates" >&2; exit 2
    fi
  }
  # assign_with_candidates: the two-run O2 split (spec §4/§7; v1, §10 names the single run over a merged BAM as the end
  # state). The families PREFIX.aug.families.txt lists are assigned on the augmented inputs (PREFIX.aug.bam = their nets
  # realigned to PREFIX.aug.fa, PREFIX.aug.copies.* = their copies + candidates, PREFIX.aug.regions.txt), every other family
  # on the original BAM and copy table; copy_assign's --only-families / --skip-families select them, so each family is
  # assigned in exactly one run (a run with no family left is not made). Every per-family table (PREFIX.assign_cand.<t>.tsv
  # and PREFIX.assign_rest.<t>.tsv: assignments, families, quant, family_join, famcn_readonly, and with --inspect psv_reads,
  # psv_cols, psv_copies, posterior) is concatenated, one header, into PREFIX.assign.<t>.tsv; the run certificates are not
  # family tables and stay per run (PREFIX.assign_cand.params.tsv, PREFIX.assign_rest.params.tsv); PREFIX.assign.log holds
  # both logs. Ruling R13: each candidate family is assigned as a family that includes its candidate copies, under the
  # unchanged AS-tied gate, so a read the realignment placed uniquely on a candidate is not AS-tied and has no row (its
  # placement is PREFIX.aug.bam's); rows arise only for molecules tied between copies. How the split differs from one run
  # over every family: copy_assign keeps a process-wide registry of the molecules tied outside the catalog (keyed by read
  # name) and demotes a family's verdicts with it, so a row's `status`, not only its `tie_outside_catalog`, can differ from
  # a single run; the cross-family steps see only their own run's families (RUSTLE_XFAM_RECONCILE's reconciliation, and
  # `--union-certificate` in a hand-made split); and in the candidate run, candidate families whose real-chromosome hulls
  # overlap register each other's tied molecules and can demote each other's verdicts.
  assign_with_candidates() {
    [ -s "$OUT.aug.bam.bai" ] && [ ! "$OUT.aug.bam.bai" -ot "$OUT.aug.families.txt" ] \
      || { echo "[rustle_pipeline] $OUT.aug.bam is missing or older than $OUT.aug.families.txt: run the candidates stage again, or drop --candidates" >&2; exit 2; }
    local n_cand n_rest t f
    n_cand=$(grep -c . "$OUT.aug.families.txt")
    n_rest=$(awk -F'\t' 'FILENAME == ARGV[1] { c[$1] = 1; next } FNR > 1 && !($1 in c) { print $1 }' "$OUT.aug.families.txt" "$OUT.fam.copies.tsv" | sort -u | wc -l)
    say "assign: $n_cand candidate families on $OUT.aug.*, the other $n_rest on $OUT.fam.copies.tsv"
    rm -f "$OUT".assign.*.tsv "$OUT".assign_cand.*.tsv "$OUT".assign_rest.*.tsv
    "$BIN/copy_assign" --bam "$OUT.aug.bam" --fasta "$OUT.aug.fa" --regions "$OUT.aug.regions.txt" --families "$OUT.aug.copies.tsv" \
      --copies-fa "$OUT.aug.copies.fa" --only-families "$OUT.aug.families.txt" "${INSPECT_ASSIGN[@]}" --out "$OUT.assign_cand" > "$OUT.assign_cand.log" 2>&1 \
      || { local rc=$?; say "assign: copy_assign failed on the candidate families (exit $rc), see $OUT.assign_cand.log; the likely cause: a real copy whose reads all realigned to its candidate, so it has no read in $OUT.aug.bam (copy_assign requires one on every real copy)"; exit "$rc"; }
    local runs=("$OUT.assign_cand")
    if [ "$n_rest" -gt 0 ]; then
      fam_regions skip "$OUT.aug.families.txt" > "$OUT.regions.txt"
      "$BIN/copy_assign" --bam "$BAM" --fasta "$FASTA" --regions "$OUT.regions.txt" --families "$OUT.fam.copies.tsv" \
        --copies-fa "$OUT.fam.copies.fa" --skip-families "$OUT.aug.families.txt" "${INSPECT_ASSIGN[@]}" --out "$OUT.assign_rest" > "$OUT.assign_rest.log" 2>&1 \
        || { local rc=$?; say "assign: copy_assign failed on the other families (exit $rc), see $OUT.assign_rest.log"; exit "$rc"; }
      runs+=("$OUT.assign_rest")
    fi
    for t in $(for r in "${runs[@]}"; do for f in "$r".*.tsv; do [ -e "$f" ] || continue; f=${f#"$r."}; echo "${f%.tsv}"; done; done | sort -u); do
      [ "$t" = params ] && continue
      local first=""
      : > "$OUT.assign.$t.tsv"
      for r in "${runs[@]}"; do
        [ -e "$r.$t.tsv" ] || continue
        if [ -z "$first" ]; then
          first=$r; cat "$r.$t.tsv" >> "$OUT.assign.$t.tsv"
        elif [ "$(head -1 "$r.$t.tsv")" = "$(head -1 "$first.$t.tsv")" ]; then
          awk 'NR > 1' "$r.$t.tsv" >> "$OUT.assign.$t.tsv"
        else
          echo "[rustle_pipeline] $r.$t.tsv and $first.$t.tsv have different headers: cannot concatenate them into $OUT.assign.$t.tsv" >&2
          rm -f "$OUT".assign.*.tsv; exit 2
        fi
      done
    done
    for r in "${runs[@]}"; do echo "=== [rustle_pipeline] $r.log"; cat "$r.log"; done > "$OUT.assign.log"
    # assigned_copy is the family's copy INDEX: family_join names the copy (copy_tid)
    local on_cand
    on_cand=$(awk -F'\t' 'FNR == 1 { next } FILENAME == ARGV[1] { tid[$1 "\t" $2] = $3; next } $4 == "assigned" && tid[$2 "\t" $3] ~ /^cand_/ { n++ } END { print n + 0 }' \
      "$OUT.assign.family_join.tsv" "$OUT.assign.assignments.tsv")
    say "assign: $(awk -F'\t' 'NR>1 && $4=="assigned"' "$OUT.assign.assignments.tsv" | wc -l) assigned rows of $(awk 'NR>1' "$OUT.assign.assignments.tsv" | wc -l) (one row per read x family; $on_cand on a candidate copy; a read the realignment placed uniquely on a candidate is not AS-tied and has no row)"
  }
  if [ "$LEGACY_CATALOG" = 1 ]; then
    say "assign: per-read copy assignment on $OUT.cat"
    samtools view -H "$BAM" | awk '$1=="@SQ"{sub("SN:","",$2); sub("LN:","",$3); print $2":1-"$3}' > "$OUT.regions.txt"
    "$BIN/copy_assign" --bam "$BAM" --fasta "$FASTA" --regions "$OUT.regions.txt" \
      --families "$OUT.cat.copies.tsv" --copies-fa "$OUT.cat.copies.fa" "${INSPECT_ASSIGN[@]}" --out "$OUT.assign" > "$OUT.assign.log" 2>&1
    say "assign: $(awk -F'\t' 'NR>1 && $4=="assigned"' "$OUT.assign.assignments.tsv" | wc -l) assigned rows of $(awk 'NR>1' "$OUT.assign.assignments.tsv" | wc -l) (one row per read x family)"
    return 0
  fi
  [ -s "$OUT.fam.copies.tsv" ] || { echo "assign needs $OUT.fam.copies.tsv (run the families stage with an mcl_families that writes the copy table), or pass --legacy-catalog" >&2; exit 2; }
  if [ "$(awk 'NR > 1 && NF' "$OUT.fam.copies.tsv" | wc -l)" = 0 ]; then
    say "assign: $OUT.fam.copies.tsv lists no copy (no family): nothing to assign"; return 0
  fi
  [ -s "$OUT.fam.copies.regions" ] || { echo "assign needs $OUT.fam.copies.regions (the families stage writes it beside the copy table)" >&2; exit 2; }
  if cand_ready && [ -s "$OUT.aug.families.txt" ]; then
    assign_with_candidates; return 0
  fi
  if [ "$CANDIDATES" != 1 ] && [ -s "$OUT.aug.families.txt" ] && [ ! "$OUT.aug.families.txt" -ot "$OUT.fam.copies.tsv" ]; then
    say "assign: the candidates stage's products ($OUT.aug.*) are present but unused (the stage is opt-in, ruling R14): pass --candidates to assign the candidate families on them"
  fi
  say "assign: per-read copy assignment on $OUT.fam.copies.tsv"
  rm -f "$OUT".assign.*.tsv "$OUT".assign_cand.* "$OUT".assign_rest.*   # an earlier split run's tables
  fam_regions all > "$OUT.regions.txt"
  "$BIN/copy_assign" --bam "$BAM" --fasta "$FASTA" --regions "$OUT.regions.txt" --families "$OUT.fam.copies.tsv" \
    --copies-fa "$OUT.fam.copies.fa" "${INSPECT_ASSIGN[@]}" --out "$OUT.assign" > "$OUT.assign.log" 2>&1 \
    || { local rc=$?; say "assign: copy_assign failed (exit $rc), see $OUT.assign.log"; exit "$rc"; }
  say "assign: $(awk -F'\t' 'NR>1 && $4=="assigned"' "$OUT.assign.assignments.tsv" | wc -l) assigned rows of $(awk 'NR>1' "$OUT.assign.assignments.tsv" | wc -l) (one row per read x family)"
}
stage_flag() {
  [ -n "$INDEX" ] || { echo "flag needs --index (splice .mmi of the primary genome)" >&2; exit 2; }
  # the loci to scan: the annotation with --gff, else the de novo loci `families` reads (bridges are relations, not loci)
  [ -n "$GFF" ] || fam_gtf_guard flag
  local LOCI=${GFF:-$FAM_GTF} cand=()
  # with --candidates, the two O3 sources corroborate: the verdict table's `o3_candidate` column names a flagged candidate on
  # the row's locus (a candidates table older than PREFIX.fam.copies.tsv was made for other families: an error)
  if [ "$CANDIDATES" = 1 ] && [ -s "$OUT.cand.candidates.tsv" ]; then
    [ ! "$OUT.cand.candidates.tsv" -ot "$OUT.fam.copies.tsv" ] \
      || { echo "[rustle_pipeline] $OUT.cand.candidates.tsv is older than $OUT.fam.copies.tsv (families ran again): run the candidates stage again, or drop --candidates" >&2; exit 2; }
    cand=(--candidates "$OUT.cand.candidates.tsv")
  fi
  say "flag: scan $BAM on $LOCI"
  "$BIN/missing_copy_flag" --bam "$BAM" --fasta "$FASTA" --loci "$LOCI" ${GFF:+--gff "$GFF"} --index x --threads "$THREADS" \
    --out "$OUT.flag_scan" --scan-only > "$OUT.flag_scan.log" 2>&1
  say "flag: align + verdict"
  "$BIN/missing_copy_flag" --bam "$BAM" --fasta "$FASTA" --loci "$LOCI" --index "$INDEX" --threads "$THREADS" \
    "${CONFIRM[@]}" "${FOREIGN[@]}" "${cand[@]}" --out "$OUT.flag" --from-scan "$OUT.flag_scan" > "$OUT.flag.log" 2>&1
  say "flag: $(grep -o 'verdicts: .*' "$OUT.flag.log")"
  if [ ${#cand[@]} -gt 0 ]; then say "flag: $(grep -o '[0-9]* flagged candidates with a nearest locus; [0-9]* rows name one' "$OUT.flag.log") (o3_candidate column)"; fi
}
case "$STAGE" in
  assemble) stage_assemble;; families) stage_families;; candidates) stage_candidates;; catalog) stage_catalog;;
  assign) stage_assign;; flag) stage_flag;;
  # --legacy-catalog: `assign` reads the legacy catalog, so `all` builds it; --candidates: `all` runs the (opt-in) candidates stage
  all) stage_assemble; stage_families
       if [ "$LEGACY_CATALOG" = 1 ]; then stage_catalog
       elif [ "$CANDIDATES" = 1 ]; then stage_candidates
       else say "candidates: skipped (opt-in, --candidates; R14: Amendment 12 failed; R22: Amendment 14 failed)"; fi
       stage_assign; stage_flag;;
  *) echo "unknown stage $STAGE" >&2; exit 2;;
esac
say "done"
