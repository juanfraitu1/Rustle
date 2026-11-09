#!/bin/bash
# rustle_pipeline.sh — the whole pipeline, one command per stage, shipped defaults (2026-09-23).
#
#   assemble  reads -> loci -> isoform GTF                              copy_assign --assemble-only --genome-wide
#   families  gene families from the de novo loci (all-vs-all -> MCL)   mcl_families --from-gtf --emit-units
#             = THE default de novo family definition (user decision 2026-09-25): one copy per member locus =
#             its representative transcript (PREFIX.fam.copies.tsv/.fa, the copy table copy assignment consumes)
#   candidates  reference-absent copies from each family's own reads    candidate_copies + utilities candidate-augment + minimap2
#             OPT-IN (ruling R14, 2026-10-02: its pre-registered acceptance, Amendment 12, FAILED —
#             docs/archive/2026-10/O3_CANDIDATES_ACCEPTANCE_2026-10-02.md; 2026-10-03: the re-run, Amendment 13, passed, but the no-deletion
#             control, Amendment 14, failed, so the default flip was reverted — docs/archive/2026-10/O3_CANDIDATES_CONTROL_A14_2026-10-03.md):
#             naming the stage runs it, `all` runs it only with --candidates.
#             (read net -> read clusters -> consensus -> flag/link/merge -> one exon-union contig per flagged
#             candidate, PREFIX.cand.*), then the augmentation (PREFIX.aug.fa = genome + the contigs, PREFIX.aug.copies.*
#             = the copy table + one row per candidate) and the patch realignment of the candidate families' reads to
#             PREFIX.aug.fa (PREFIX.aug.bam): a read whose best alignment is on a candidate contig is PLACED there by the
#             realignment (ruling R13). --delta D (0.00958), --cand-max-reads N (1000). Cost: PREFIX.aug.fa is a copy of
#             the whole genome, and the realignment indexes it on every run (about 400 s and 19 GB on a human genome).
#             The stage's own cost on a full BAM (spec §9b, R23): on the 23-GB gorilla fibroblast BAM one batch of 50 families
#             did not finish in a 10-minute call (nets phase 397-469 s: the whole-BAM sweep ~150 s, the attribution alignment
#             against the batch's 2.6 GB of net reads 213-250 s; peak RSS 10.7 GB at the stop). The driver runs EVERY family in
#             one call (no batching; the binary's `--families` batches): the attribution targets are then all families' nets, so
#             its peak memory exceeds that batch's 10.7 GB — an OOM risk on this machine.
#   catalog   LEGACY copy catalog (gw_family_catalog; kept, not the default definition)
#   assign    per-read copy assignment on the catalog (assign/abstain)  copy_assign --families
#   all --to STAGE / --from STAGE   stop after / start at a stage (escape hatch; `--to families` answers "what do the loci
#             look like before reads are assigned": PREFIX.fam.copies.tsv + PREFIX.fam.loci.gff3). 2026-10-30.
#   merged    THE genome-wide entry point: assemble + families + assign on the FAMILIES copy table, one command,
#             byte-identical products, with the genome-wide levers ON by default (--skip-poa-diagnostic and
#             --region-threads, both cmp-proven byte-identical) and per-phase resume (a finished phase is never
#             redone; the families all-vs-all PAF replays from PREFIX.cache when the loci FASTA is unchanged).
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
#   sd        the genome's segmental duplications and duplicons            utilities sd-scan (vendored sedef-rs)
#             once per genome, cached by the FASTA's md5 and sd-scan's version key in RUSTLE_SD_CACHE_DIR (default
#             PREFIX.cache)/sd/<key>/genome.sd.*; PREFIX.sd.prefix names that prefix for `mcl_families --emit-hierarchy
#             --sd-cache`. Needs --fasta and --out only; `all` does not run it. An entry is a hit only if its DONE lists
#             every file with its size (an empty DONE, as this stage wrote before 2026-10-14, is a miss); --no-cache scans
#             again, into the same entry.
#             Gates in docs/archive/2026-11/LEDGER_O1.md:4038 (PREREG_sedef_rs_layers_2026-10-13.md) (register 1278-1280): sd-scan's SD calls failed V1/V1b
#             (recall about 0.40 against SEDEF) and its duplicons failed V2 (held-out family score 0.5153 against
#             DupMasker), so `--sd-cache` is an opt-in, unvalidated alternative to `--hierarchy-sd` /
#             `--hierarchy-duplicons`; its SD layer uses every `ok` pair (no identity or length floor). Cost with
#             --threads 4: about 18 min and 10 GiB peak per human genome, about 55 min and 19 GiB for gorilla. The stage
#             sets no memory cap (bench/sedef/run_genome.sh runs it under `ulimit -v`, 23 GiB).
# (thesis record: families/catalog = O1, candidates = O3, assign = O2, flag = O3; in order O1 families -> O3 candidates ->
#  O2 assign -> O3 flag)
#
# usage: tools/rustle_pipeline.sh STAGE --bam B --fasta G --out PREFIX [--index G.splice.mmi] [--gff ANNOT.gff]
#        [--confirm NAME=X.mmi ...] [--foreign NAME=X.mmi ...] [--threads N] [--bin DIR]
#        [--seed-pool primary|good|all [--seed-as-ratio R]] [--as-table FILE] [--contig NAME[,NAME...]] [--no-seed-secondaries]
#        [--no-cache] [--inspect] [--piecewise [--max-pieces N] [--budget-s S] [--piece-records R] [--piece LABEL]]
#        [--candidates | --no-candidates] [--delta D] [--cand-max-reads N] [--legacy-catalog]
#   tools/rustle_pipeline.sh cache-ls --out PREFIX      list what PREFIX.cache holds (cache-clear: delete it)
#   tools/rustle_pipeline.sh sd --fasta G --out PREFIX [--threads N] [--bin DIR] [--no-cache]   the genome's SDs, cached (stage sd)
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
#   (one `as_table` pass over the BAM -> BAM.as_table.tsv next to the BAM, shared by every run on it, + its binary load
#   sidecar .asbin; PREFIX.molecules.tsv links to it; PREFIX.molecules.tsv itself when the BAM's directory is read-only): on gorilla NC_073244.2 (the
#   seeding pre-registration's verdict contig) this finds 30 more loci (97% annotated) and joins more
#   >=90%-identity referee pairs (21 -> 26 of 83 on genes with >= 2 exonic primary reads; all of the gain is one
#   tandem array) at pair precision 1.000 (rows 1060/1100/1101/1116); the streaming assembler applies the filter
#   at no memory cost (validated identical to the buffered path). --no-seed-secondaries restores primaries-only seeding (2026-09-24 default flip).
#   The table is only as genome-wide as the BAM: on a region SLICE it admits secondaries whose real best lies
#   outside the slice (tes44: 4,210 transcripts vs 3,911 with the full-BAM table) — give the driver the full BAM.
# SEEDING POOL, the switch (2026-10-07): which alignments seed the assembly. One flag, three pools, everything else equal:
#   --seed-pool primary   primary alignments only                                   (= --no-seed-secondaries)
#   --seed-pool good      primaries + secondaries with AS >= R x the molecule's genome-wide best AS, R = 0.98
#                         (THE DEFAULT; = --seed-secondaries; needs the as_table pass above)
#   RUSTLE_SEED_UNSPLICED=primary|support (opt-in, passed through to the binary; r1518/r1520): what the pool does with an UNSPLICED
#                         secondary: primary = drop it; support = it adds support to the primary-seeded cluster it overlaps but never
#                         seeds or bridges one (array-wide single-exon blobs gone; 13% of CAT single-exon pseudogene hits lost).
#   --seed-pool all       primaries + EVERY secondary (RUSTLE_GTF_SECONDARY=1, no filter, no table): the 2026-09-21 pool; in the
#                         2026-10-07 study 1.4-1.5x the transcripts of good on every contig, and the aligner
#                         records of good x17.5 (chr16; x22.9 of primary), x14.0 (chr17; x15.0), x3.1 (gorilla NPIP; x3.9)
#                         (human chr16: 14,183 against 9,473 transcripts, 1,648,588 against 94,433 PAF records;
#                         docs/archive/2026-11/LEDGER_O1.md:2696 (SEED_POOL_REAL_READS_2026-10-07.md))
#   --seed-as-ratio R     the width R of the good pool, a number in (0, 1]; 1 = exact ties only. Refused with another pool.
#   RUSTLE_SEED_POOL / RUSTLE_SEED_AS_RATIO set the same two as defaults (a flag wins; a malformed exported ratio is refused,
#   a valid one is ignored by the other pools). The pool the flag names is the pool the binary gets: an exported
#   RUSTLE_GTF_SECONDARY* of the shell is overridden (primary: unset; all: ratio and table unset), so those three variables
#   cannot contaminate a comparison between pools (other exported RUSTLE_* read-pool knobs still reach the binary;
#   bench/seed_pool/run.sh refuses any exported RUSTLE_*). One behaviour change: `--no-seed-secondaries` with an exported
#   RUSTLE_GTF_SECONDARY* used to let it through, and now unsets it (the default and the legacy flags are otherwise
#   byte-identical to the earlier driver, stub-checked). If the best-AS table is empty or from another BAM the binary says
#   so in PREFIX.assemble.log and the good pool fails with status 3 instead of running a wider pool. The
#   pool is announced in the log ("seeding pool: ..."), recorded in PREFIX.merged.env when it is not the default, and
#   changes nothing but that environment: the arguments of copy_assign are the same for every pool. Use a different
#   --out PREFIX per pool (bench/seed_pool/run.sh does, and prints the comparison). --as-table FILE reads (or builds, when absent
#   or made from another BAM) the best-AS table there instead of PREFIX.molecules.tsv, so many pools and widths share one pass over
#   the BAM. --contig NAME[,NAME...] (assemble) assembles those contigs only, copy_assign --region NAME:0-LEN (--regions FILE for
#   several) with LEN from FASTA.fai, in place of --genome-wide: the products equal the genome-wide run's on that contig (docs/archive/2026-09/CONTAINER_HEADROOM_2026-09-30.md) and it
#   makes a pool comparison on one contig a matter of minutes; the families, assign and flag stages read whatever assembly is there.
# The splice index is needed by `candidates` (the consensus sequences' genome hits) and `flag` (home search); the
# annotation (--gff) by `flag` only (IG/TR screen).
# Environment: RUSTLE_POLISH_SUBCHAIN=tag|drop adds `--polish-subchain` to `assemble` (default unset = off);
#   RUSTLE_POLISH_TSS=off|tag|rescue|split sets `--polish-tss` for `assemble` (default unset = rescue since 2026-11-01, register 1485);
#   RUSTLE_POLISH_TSS_END=median adds `--polish-tss-end` to `assemble` (default unset = off; register 1487);
#   RUSTLE_POLISH_TES_EXTEND=on adds `--polish-tes-extend` to `assemble` (default unset = off; register 1488);
#   RUSTLE_MONO_MODEL=cluster adds `--mono-model` to `assemble` (default unset = off; register 1489);
#   RUSTLE_FINAL_SUBCHAIN=drop|off: the sub-chain drop as a fourth final.gtf layer (default drop, register 1498);
#   RUSTLE_ASSEMBLY_LAYERS=on|off: `assemble` also writes PREFIX.final.gtf, the assembly with the three end layers, from a second
#   parallel run (default on; PREFIX.gtf / PREFIX.families.gtf unchanged by construction; registers 1487-1491);
#   RUSTLE_POLISH_TES=tag|pas-end adds `--polish-tes` to `assemble` (default unset = off);
#   RUSTLE_POLISH_MONO_PRIMED=0 adds `--no-polish-mono-primed` to `assemble` (default unset = 1 = on since 2026-10-26);
#   RUSTLE_POLISH_JUNCTION_SNAP=equiv|reads adds `--polish-junction-snap` to `assemble` (default unset = off);
#   RUSTLE_GTF_REGROUP=1 adds `--gtf-regroup` to `assemble` (RG3 regroup after polish; default unset = off; needs
#     RUSTLE_BRIDGE_REGROUP=off, whose default contains RG3's split);
#   RUSTLE_BRIDGE_REGROUP=off|f1|f1v2|f1units sets `--bridge-regroup` on `assemble` (unset = f1v2, THE DEFAULT since
#     2026-09-29; off = the 2026-09-25 products, byte for byte). With f1, f1v2 or f1units, `families`, and `flag` without
#     --gff, read PREFIX.families.gtf (f1 / f1v2: the assembly without its bridge transcripts, bridges are relations, not
#     loci; f1units, OPT-IN: each bridge transcript replaced by its UNITS, the exon runs between its bridge introns, see
#     copy_assign --help); every stage derives the mode the same way (BRIDGE_MODE below);
#   RUSTLE_BRIDGE_UNITS_LIST=FILE (f1units only; no whitespace in the path) adds `--bridge-units-list FILE` to `assemble`:
#     the transcripts and introns to cut come from FILE (any detector, e.g. bench/units_from_annotation.py) and F1's own
#     read evidence is not consulted;
#   RUSTLE_FAMILY_RELATIONS=1 adds `--emit-relations` to `families` (PREFIX.fam.relations.tsv, PREFIX.fam.members_by_locus.tsv:
#     one relation record per split transcript and the members of each family by locus; default unset = off);
#   RUSTLE_MIN_COV_SHORTER=C sets `mcl_families --min-cov-shorter C` on `families` (unset = the binary's own default,
#     0.70 since 2026-09-29, nothing passed; 0 = the 2026-09-25 edge weights);
#   RUSTLE_REPRESENTATIVE=most-reads|most-junctions sets `mcl_families --representative` on `families`, the transcript that
#     represents each de novo locus (loci.gff3, the copy table). Unset = nothing passed: the binary's default, most-reads
#     (the transcript with the most reads), byte for byte the earlier products; most-junctions = the transcript with the most
#     junctions (gaps >= 50 bp between exons), ties to the most reads — OPT-IN: the arm of
#     docs/archive/2026-10/PREREG_locus_representative_rule_2026-10-04.md, whose held-out decision (2026-10-04) kept most-reads: chr6 clause (a)
#     FAILED, Compara F 0.444 -> 0.372, on pre-f1v2 GTFs (docs/archive/2026-10/LOCUS_REPRESENTATIVE_RULE_2026-10-04.md); an mcl_families
#     without the flag is refused;
#   RUSTLE_FAMILY_CONTAINER=1 adds `--emit-container` to `families` (PREFIX.fam.container*.tsv; default unset = off).
#   RUSTLE_FAMILY_EDGE_RULE=share|coverage (default share since 2026-11-09, row 1593): `families` builds edges by exon sharing
#     (`--min-shared-exon-frac 0.0 --min-exon-share ${RUSTLE_MIN_EXON_SHARE:-0.5} --inflation ${RUSTLE_INFLATION:-2.0}`);
#     coverage = the earlier rule (`--min-shared-exon-frac ${RUSTLE_MIN_SHARED_EXON_FRAC:-0.60}`, inflation 2.8), also chosen when
#     RUSTLE_MIN_SHARED_EXON_FRAC is set. The share rule aligns more pairs (exon gate at fraction 0): about 2x the families wall;
#   RUSTLE_SEDEF_BED=FILE turns on SD damping in `families` (unset: FASTA.sedef.bed when that sidecar exists; empty: off) (`--sd-damp-bed FILE --sd-damp ${RUSTLE_SD_DAMP:-0.25}`): edges
#     between loci no SEDEF pair links are down-weighted before MCL (row 1580: HG002 held-out +0.160 precision at equal
#     recall; A119b .802/.740 -> .925/.718). Set it whenever the genome has a SEDEF build (human CHM13, gorilla; paths in
#     docs/DATA.md); RUSTLE_SD_DAMP=1 keeps the BED but turns damping off. Unset = off, the earlier products.
# `families` is the DE NOVO mode (loci from the assembled GTF). The GUIDED mode (loci = the annotation's gene and
# pseudogene bodies, PREREG_heldout_families_2026-09-20 §2) is not a driver stage: figures/_o1_recovery.py
# (guided_families) runs its recipe step by step. Every product carries the PREFIX.
set -euo pipefail
# Run from a snapshot of this file (round-9 audit): bash reads a script as it executes it, so an edit or `git pull` of the
# repo during a multi-hour run spliced old and new code (seen 2026-11-09: "unexpected EOF" after a 10-minute assemble).
# The snapshot keeps size and mtime (cp -p), so the resume records' driver line is unchanged; RP_DRIVER_DIR (no RUSTLE_
# prefix: they must not enter the resume records) keeps the repo's tools/ directory for the sibling scripts and --bin.
# A stage call of `all` runs the same snapshot (its $0 is RP_SNAPSHOT); a call of the repo file (a user, rustle_samples.sh) snapshots.
if [ "${RP_SNAPSHOT:-}" != "$0" ]; then
  export RP_DRIVER_DIR RP_SNAPSHOT RP_OWNER_PID=$$; RP_DRIVER_DIR=$(cd "$(dirname "$0")" && pwd)
  RP_SNAPSHOT=$(mktemp "${TMPDIR:-/tmp}/rustle_pipeline.XXXXXX") && cp -p "$0" "$RP_SNAPSHOT" && exec bash "$RP_SNAPSHOT" "$@"
fi
# the process that made the snapshot (exec keeps its pid) removes it at exit; stage calls (other pids) leave it to their parent
if [ "${RP_OWNER_PID:-}" = "$$" ]; then trap 'rm -f "$RP_SNAPSHOT"' EXIT; fi
STAGE=${1:-all}; shift || true
ORIG_ARGS=("$@")
# Multi-sample runs (a comma in --bam): tools/rustle_samples.sh takes the whole call (pool, then split). A `case "$STAGE"` block
# at column 0, so figures/samples.py's driver_stage_code leaves it out of every stage's hash; one BAM never enters it.
case "$STAGE" in
  *) for ((i = 0; i + 1 < ${#ORIG_ARGS[@]}; i++)); do
       if [ "${ORIG_ARGS[$i]}" = --bam ] && [[ "${ORIG_ARGS[$((i + 1))]}" == *,* ]]; then
         exec bash "$RP_DRIVER_DIR/rustle_samples.sh" "$STAGE" "${ORIG_ARGS[@]}"
       fi
     done;;
esac
BAM=""; FASTA=""; OUT=""; INDEX=""; GFF=""; THREADS=4; THREADS_SET=0; BIN="$RP_DRIVER_DIR/../target/release"; CONFIRM=(); FOREIGN=(); SEED_POOL=${RUSTLE_SEED_POOL:-good}; SEED_RATIO=${RUSTLE_SEED_AS_RATIO:-0.98}; SEED_RATIO_FLAG=0; AS_TABLE=""; CONTIG=""; CONTIG_FLAG=0; CACHE=1; INSPECT=0
PIECEWISE=0; MAX_PIECES=0; BUDGET_S=0; PIECE_RECORDS=0; PIECE=""
LEGACY_CATALOG=0; CANDIDATES=0; DELTA=0.00958; CAND_MAX=1000
# `all --to STAGE` stops after STAGE (the escape hatch: `--to families` = the loci and their copy table before any read is
# assigned, 2026-10-30 NEXT TARGETS item 2); `--from STAGE` starts there (the earlier products must exist).
TO_STAGE=; FROM_STAGE=
INDEX_GB=8; LOW_MEM=0; NO_FLAG=0; RESUME=1; MAX_WINDOW_RECORDS=${RUSTLE_MAX_WINDOW_RECORDS:-0}
while [ $# -gt 0 ]; do
  case "$1" in
    --bam) BAM=$2; shift 2;; --fasta) FASTA=$2; shift 2;; --out) OUT=$2; shift 2;; --index) INDEX=$2; shift 2;;
    --gff) GFF=$2; shift 2;; --threads) THREADS=$2; THREADS_SET=1; shift 2;; --bin) BIN=$2; shift 2;;
    --confirm) CONFIRM+=(--confirm "$2"); shift 2;; --foreign) FOREIGN+=(--foreign "$2"); shift 2;;
    --seed-pool) SEED_POOL=$2; shift 2;; --seed-as-ratio) SEED_RATIO=$2; SEED_RATIO_FLAG=1; shift 2;; --as-table) AS_TABLE=$2; shift 2;;
    --contig) CONTIG=$2; CONTIG_FLAG=1; shift 2;;
    --seed-secondaries) SEED_POOL=good; shift;; --no-seed-secondaries) SEED_POOL=primary; shift;;
    --cache) CACHE=1; shift;; --no-cache) CACHE=0; shift;; --inspect) INSPECT=1; shift;;
    --piecewise) PIECEWISE=1; shift;; --max-pieces) MAX_PIECES=$2; shift 2;; --budget-s) BUDGET_S=$2; shift 2;;
    --piece-records) PIECE_RECORDS=$2; shift 2;; --piece) PIECE=$2; shift 2;;
    --legacy-catalog) LEGACY_CATALOG=1; shift;; --candidates) CANDIDATES=1; shift;; --no-candidates) CANDIDATES=0; shift;;
    --to) TO_STAGE=$2; shift 2;; --from) FROM_STAGE=$2; shift 2;;
    --delta) DELTA=$2; shift 2;; --cand-max-reads) CAND_MAX=$2; shift 2;;
    --index-gb) INDEX_GB=$2; shift 2;; --low-mem) LOW_MEM=1; shift;; --no-flag) NO_FLAG=1; shift;; --no-resume) RESUME=0; shift;;
    --max-window-records) MAX_WINDOW_RECORDS=$2; shift 2;;
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
# The seeding pool (see the SEEDING POOL note above): a misspelling must not quietly run another pool. Same `case "$STAGE"`
# form as the check above.
case "$STAGE" in
  *) case "$SEED_POOL" in
       primary|good|all) ;;
       *) echo "[rustle_pipeline] --seed-pool (RUSTLE_SEED_POOL) must be primary, good or all (got '$SEED_POOL')" >&2; exit 2;;
     esac
     [[ "$SEED_RATIO" =~ ^(1(\.0+)?|0?\.[0-9]*[1-9][0-9]*)$ ]] || { echo "[rustle_pipeline] --seed-as-ratio (RUSTLE_SEED_AS_RATIO) must be a number in (0, 1] (got '$SEED_RATIO')" >&2; exit 2; }
     if [ "$SEED_RATIO_FLAG" = 1 ] && [ "$SEED_POOL" != good ]; then
       echo "[rustle_pipeline] --seed-as-ratio sets the width of the good pool and only the good pool has one (--seed-pool $SEED_POOL): drop one" >&2; exit 2
     fi
     if [ "$CONTIG_FLAG" = 1 ] && [ -z "$CONTIG" ]; then
       echo "[rustle_pipeline] --contig is empty (an empty shell variable?): name a contig, or drop the flag for a genome-wide run" >&2; exit 2
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
# sd (docs/superpowers/specs/2026-10-13-o1-layered-hierarchy-design.md): the genome's segmental duplications with the
# vendored sedef-rs (`utilities sd-scan`), once per genome: cached in RUSTLE_SD_CACHE_DIR (default PREFIX.cache) under
# sd/<md5 of the FASTA>-<md5 of sd-scan's version key, 8 hex>/genome.sd.*. PREFIX.sd.prefix names the cached prefix
# that `mcl_families --emit-hierarchy --sd-cache` takes; it is removed first and written only after a hit or a finished
# scan, so a failed run leaves none. Needs --fasta and --out only (no BAM), so it is a stage function plus a
# `case "$STAGE"` block that exits before the shared --bam check: figures/samples.py's driver_stage_code leaves every
# other stage's code hash as it was.
# The entry is marked done by DONE, written last: one `name<TAB>bytes` line per file of the entry (every genome.sd.* file
# and key.tsv), as src/family.rs's cache entries do. A hit needs the cache on (--no-cache always scans, into the same
# entry), a non-empty DONE and every listed file present with its listed size; anything else is a miss (an empty DONE,
# as this stage wrote before 2026-10-14, included). A miss removes DONE, scans, writes key.tsv, syncs the files, writes
# DONE.tmp and renames it to DONE.
stage_sd() {
  [ -n "$FASTA" ] && [ -n "$OUT" ] || { echo "[rustle_pipeline] sd needs --fasta and --out" >&2; exit 2; }
  [ -s "$FASTA" ] || { echo "[rustle_pipeline] sd: $FASTA is missing or empty" >&2; exit 2; }
  rm -f "$OUT.sd.prefix"
  local key md5 dir hit=0 f="" n="" b=""
  key=$("$BIN/utilities" sd-scan --version-key) || { echo "[rustle_pipeline] sd: $BIN/utilities has no sd-scan; rebuild it" >&2; exit 2; }
  md5=$(md5sum < "$FASTA" | cut -c1-32)
  dir=${RUSTLE_SD_CACHE_DIR:-$OUT.cache}/sd/$md5-$(printf '%s' "$key" | md5sum | cut -c1-8)
  # one stage per cache entry at a time (held until this stage exits): two misses on a shared cache used to rescan into
  # the same files, truncating the tables the run that finished first was reading
  # (opened for append, read-only when another user's lock file is not ours to write: a shared-cache hit then still runs;
  # NFS needs the writable fd for an exclusive lock, a read-only one locks shared)
  mkdir -p "${dir%/*}"
  # (the braces keep the 2>/dev/null to this open: on a bare exec it would silence the rest of the stage)
  local ro=0
  if { exec 8>> "$dir.lock"; } 2>/dev/null; then flock 8; else exec 8< "$dir.lock" && flock -s 8; ro=1; fi
  if [ "$CACHE" = 1 ] && [ -s "$dir/DONE" ]; then
    hit=1
    while IFS=$'\t' read -r n b || [ -n "$n" ]; do
      [ -n "$n" ] && [ -f "$dir/$n" ] && [ "$(wc -c < "$dir/$n")" = "$b" ] || { hit=0; break; }
    done < "$dir/DONE"
  fi
  if [ "$hit" = 1 ]; then
    echo "[rustle_pipeline] $(date +%H:%M:%S) sd: cached ($dir)" >&2
  else
    # a shared lock (the lock file is not ours to write) does not keep a second scanner out: scan only under our own lock
    [ "$ro" = 0 ] || { echo "[rustle_pipeline] sd: $dir is not a cache hit and $dir.lock is not writable by this user: scan with your own RUSTLE_SD_CACHE_DIR" >&2; exit 2; }
    mkdir -p "$dir"
    rm -f "$dir/DONE" "$dir/DONE.tmp"
    echo "[rustle_pipeline] $(date +%H:%M:%S) sd: utilities sd-scan on $FASTA -> $dir" >&2
    "$BIN/utilities" sd-scan --genome "$FASTA" --out "$dir/genome" --threads "$THREADS" > "$OUT.sd.log" 2>&1 \
      || { echo "[rustle_pipeline] sd: sd-scan failed (exit $?); see $OUT.sd.log" >&2; exit 1; }
    printf 'sd\t%s\t%s\t%s\n' "$md5" "$key" "$FASTA" > "$dir/key.tsv"
    sync "$dir"/genome.sd.* "$dir/key.tsv" || { echo "[rustle_pipeline] sd: could not sync the files of $dir; not marked done" >&2; exit 1; }
    for f in "$dir"/genome.sd.* "$dir/key.tsv"; do
      if [ -f "$f" ]; then printf '%s\t%s\n' "${f##*/}" "$(wc -c < "$f")"; fi
    done > "$dir/DONE.tmp"
    sync "$dir/DONE.tmp" && mv "$dir/DONE.tmp" "$dir/DONE" || { echo "[rustle_pipeline] sd: could not write $dir/DONE" >&2; exit 1; }
  fi
  printf '%s\n' "$dir/genome" > "$OUT.sd.prefix"
  echo "[rustle_pipeline] $(date +%H:%M:%S) sd: $(awk -F'\t' '!/^#/ && $11=="ok"' "$dir/genome.sd.pairs.tsv" | wc -l) SD pairs (ok), $(awk '!/^#/' "$dir/genome.sd.families.tsv" | wc -l) duplication families; --sd-cache $dir/genome ($OUT.sd.prefix)" >&2
}
case "$STAGE" in
  sd) stage_sd; exit 0;;
esac
if [ "$STAGE" = index ]; then
  [ -n "$FASTA" ] && [ -n "$OUT" ] || { echo "index needs --fasta and --out" >&2; exit 2; }
else
  [ -n "$BAM" ] && [ -n "$FASTA" ] && [ -n "$OUT" ] || { echo "need --bam, --fasta, --out" >&2; exit 2; }
fi
# --low-mem (docs/archive/2026-11/LEDGER_O1.md:3502 (GENOMEWIDE_FAST_2026-10-12.md)): the splice index in 1-Gbp parts (flag holds one part at a time) and one
# assign region at a time; every other setting unchanged
[[ "$MAX_WINDOW_RECORDS" =~ ^[0-9]+$ ]] || { echo "[rustle_pipeline] --max-window-records must be a whole number (got '$MAX_WINDOW_RECORDS')" >&2; exit 2; }
if [ "$LOW_MEM" = 1 ]; then INDEX_GB=1; export RUSTLE_REGION_THREADS=1; [ "$MAX_WINDOW_RECORDS" -gt 0 ] || MAX_WINDOW_RECORDS=500000; fi
[[ "$INDEX_GB" =~ ^[1-9][0-9]*$ ]] || { echo "[rustle_pipeline] --index-gb must be a positive integer (got '$INDEX_GB')" >&2; exit 2; }
# skip lists that `merged` and the --candidates assign do not apply: refused here, before assemble and families run (round-23 audit);
# the stage functions refuse them again. Same `case "$STAGE"` form as the checks above.
case "$STAGE" in
  *) if [ "$MAX_WINDOW_RECORDS" -gt 0 ] || [ -n "${RUSTLE_ANNOTATED_ONLY:-}" ] || [ -n "${RUSTLE_EXCLUDE_BED:-}" ]; then
       if [ "$STAGE" = merged ]; then
         echo "[rustle_pipeline] merged applies no --max-window-records / --low-mem / RUSTLE_ANNOTATED_ONLY / RUSTLE_EXCLUDE_BED: use \`all\` for them" >&2; exit 2
       elif [ "$CANDIDATES" = 1 ] && { [ "$STAGE" = all ] || [ "$STAGE" = assign ]; }; then
         echo "[rustle_pipeline] --candidates does not apply --max-window-records / --low-mem / RUSTLE_ANNOTATED_ONLY / RUSTLE_EXCLUDE_BED: drop one of them" >&2; exit 2
       fi
     fi;;
esac
export TMPDIR=${TMPDIR:-/tmp}
if [ "$CACHE" = 1 ]; then export RUSTLE_CACHE_DIR="$OUT.cache"; else unset RUSTLE_CACHE_DIR; fi
INSPECT_CAT=(); INSPECT_ASSIGN=()
# ASSIGN_FAST (docs/archive/2026-11/LEDGER_O1.md:3502 (GENOMEWIDE_FAST_2026-10-12.md)): every copy_assign of the assign stage skips the POA homology DIAGNOSTIC
# (byte-identical families / assignments / abundance; ~85% of the wall time on dense families; only the diagnostic
# .fallback.tsv and its log counts are lost) and runs RUSTLE_REGION_THREADS regions at once (byte-identical). Defaults since
# r1543 (docs/PREREG_pipeline_speed_2026-11-05.md gw7, A119b genome-wide: 6,589 s vs 19,358 s, peak 199 GB on a 250 GB node):
# 8 regions, 4 family threads x 16 star threads, 4M-record read budget. LOW_MEM=1 keeps regions sequential.
# RUSTLE_POA_DIAGNOSTIC=1 restores the diagnostic pass.
case "${RUSTLE_REGION_THREADS:-8}" in ''|*[!0-9]*|0) echo "[rustle_pipeline] RUSTLE_REGION_THREADS must be a positive integer (got '${RUSTLE_REGION_THREADS:-}')" >&2; exit 2;; esac
ASSIGN_FAST=(--region-threads "${RUSTLE_REGION_THREADS:-8}")
export RUSTLE_FAMILY_THREADS="${RUSTLE_FAMILY_THREADS:-4}" RUSTLE_STAR_THREADS="${RUSTLE_STAR_THREADS:-16}" RUSTLE_REGION_READ_BUDGET="${RUSTLE_REGION_READ_BUDGET:-4000000}"
[ "${RUSTLE_POA_DIAGNOSTIC:-0}" = 1 ] || ASSIGN_FAST+=(--skip-poa-diagnostic)
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
# RUSTLE_POLISH_TSS=off|tag|rescue|split: copy_assign --polish-tss, the read-proven TSS (tag: `tss_clusters` only; rescue:
# also keeps proven short-TSS forms the polish dropped; split: also emits a chain once per proven TSS cluster). Unset = rescue,
# the default since 2026-11-01 (register 1485, PREREG ASSEMBLY ENDS E1: held-out A119b chr11 / chr7 chains +7 / +9 at
# transcript precision -0.07 / -0.01 points); RUSTLE_POLISH_TSS=off restores the earlier output (`--polish-tss off` is byte-identical to no flag).
POLISH_TSS=${RUSTLE_POLISH_TSS:-rescue}
case "$POLISH_TSS" in
  off|tag|rescue|split) POLISH="$POLISH --polish-tss $POLISH_TSS" ;;
  *) echo "[rustle_pipeline] RUSTLE_POLISH_TSS must be off, tag, rescue or split (got '$RUSTLE_POLISH_TSS')" >&2; exit 2 ;;
esac
# RUSTLE_POLISH_TSS_END=median (opt-in; unset = off, the same command): copy_assign --polish-tss-end, each multi-exon transcript's
# 5' end at the lower median of its exact-chain reads' starts (register 1487, PASSED held-out chr11 / chr7: starts within 50 bp
# of CAT 40.6 -> 64.6% / 39.3 -> 63.4%, chains identical); needs --polish-tss (the default rescue). Default pending the user.
case "${RUSTLE_POLISH_TSS_END:-}" in
  "") ;;
  off|median) POLISH="$POLISH --polish-tss-end $RUSTLE_POLISH_TSS_END" ;;
  *) echo "[rustle_pipeline] RUSTLE_POLISH_TSS_END must be off or median (got '$RUSTLE_POLISH_TSS_END')" >&2; exit 2 ;;
esac
# RUSTLE_POLISH_TES_EXTEND=on (opt-in; unset = off, the same command): copy_assign --polish-tes-extend, a multi-exon model's 3' end
# extended to a PAS-proven cluster of unspliced reads over >= 2-read coverage (register 1488, PASSED held-out chr11 / chr7: base
# sensitivity +2.1 / +2.4, precision unchanged, chains identical); needs --polish-tss (the default rescue). Default pending the user.
case "${RUSTLE_POLISH_TES_EXTEND:-}" in
  "") ;;
  off|on) POLISH="$POLISH --polish-tes-extend $RUSTLE_POLISH_TES_EXTEND" ;;
  *) echo "[rustle_pipeline] RUSTLE_POLISH_TES_EXTEND must be off or on (got '$RUSTLE_POLISH_TES_EXTEND')" >&2; exit 2 ;;
esac
# RUSTLE_MONO_MODEL=cluster (opt-in; unset = off, the same command): copy_assign --mono-model, single-exon models rebuilt from
# end-clustered unspliced reads (register 1489, PASSED held-out chr11 / chr7: base precision +17.2 / +16.5, transcript precision
# +0.2 / +0.4, chains identical); needs --polish-tss (the default rescue). Default pending the user.
case "${RUSTLE_MONO_MODEL:-}" in
  "") ;;
  off|cluster) POLISH="$POLISH --mono-model $RUSTLE_MONO_MODEL" ;;
  *) echo "[rustle_pipeline] RUSTLE_MONO_MODEL must be off or cluster (got '$RUSTLE_MONO_MODEL')" >&2; exit 2 ;;
esac
# RUSTLE_POLISH_TES=tag|pas-end (opt-in; unset = off, the same command): copy_assign --polish-tes, the read + genome
# proven TES (tag: `tes_clusters` / `tes_pas` / `tes_primed` only; pas-end: also moves an internally primed 3' end to
# the most-3' PAS-proven cluster of the transcript's own reads). Dev-only, in-sample evidence: see its --help
case "${RUSTLE_POLISH_TES:-}" in
  "") ;;
  off|tag|pas-end) POLISH="$POLISH --polish-tes $RUSTLE_POLISH_TES" ;;
  *) echo "[rustle_pipeline] RUSTLE_POLISH_TES must be off, tag or pas-end (got '$RUSTLE_POLISH_TES')" >&2; exit 2 ;;
esac
# RUSTLE_POLISH_MONO_PRIMED (default ON since 2026-10-26, user decision; unset or 1 = on, the same command; 0 = off):
# copy_assign drops a single-exon transcript when >= 80% of its unspliced reads are internally primed (StringTie's
# long-read poly-A filter, genomic form; held-out PASS bench/stgaps/PREREG_mono_primed_2026-10-26.md)
case "${RUSTLE_POLISH_MONO_PRIMED:-}" in
  ""|1) ;;
  0) POLISH="$POLISH --no-polish-mono-primed" ;;
  *) echo "[rustle_pipeline] RUSTLE_POLISH_MONO_PRIMED must be 0 or 1 (got '$RUSTLE_POLISH_MONO_PRIMED')" >&2; exit 2 ;;
esac
# RUSTLE_COMPAT_SUPPORT=1 (opt-in; unset or 0 = off, the same command): copy_assign --compat-support, compatible-read support in
# the assembler's chain admission (a chain short of the read floor is admitted when every junction is carried by enough reads
# whose chains are contiguous, end-compatible sub-chains of it; docs/archive/2026-11/LEDGER_O1.md:3055 (PREREG_assembler_compatible_support_2026-10-11.md))
case "${RUSTLE_COMPAT_SUPPORT:-}" in
  ""|0) ;;
  1) POLISH="$POLISH --compat-support" ;;
  *) echo "[rustle_pipeline] RUSTLE_COMPAT_SUPPORT must be 0 or 1 (got '$RUSTLE_COMPAT_SUPPORT')" >&2; exit 2 ;;
esac
# RUSTLE_DIRECT_TRANSFRAGS=1|union (opt-in; unset or 0 = off, the same command): copy_assign --direct-transfrags (or, with union,
# --direct-transfrags-union: plus every chain the default admits), transcripts from
# maximal transfrags with the reads of the transfrags they contain transferred to them, no flow
# (docs/archive/2026-11/LEDGER_O1.md:3131 (PREREG_direct_transfrags_2026-10-11.md)); exclusive with RUSTLE_COMPAT_SUPPORT
case "${RUSTLE_DIRECT_TRANSFRAGS:-}" in
  ""|0) ;;
  1) POLISH="$POLISH --direct-transfrags" ;;
  union) POLISH="$POLISH --direct-transfrags-union" ;;
  *) echo "[rustle_pipeline] RUSTLE_DIRECT_TRANSFRAGS must be 0, 1 or union (got '$RUSTLE_DIRECT_TRANSFRAGS')" >&2; exit 2 ;;
esac
# Opt-in levers toward FLAIR's exclusive transcripts (docs/archive/2026-11/LEDGER_GENERAL.md:511 (PENDING_2026-10-11.md); unset = off, the same command):
#   RUSTLE_TRANSFER_CREDIT=full|best copy_assign --transfer-credit full|best (with RUSTLE_DIRECT_TRANSFRAGS)
#   RUSTLE_READ_JUNCTION_SNAP=W    copy_assign --read-junction-snap W (bp; read junctions snap to better-supported ones)
#   RUSTLE_MIN_CHAIN_READS=N       copy_assign --min-chain-reads N (1 = single-read spliced chains)
#   RUSTLE_READ_STRAND=0           assemble: single-exon strand back to the '+' placeholder (default 1 = read-orientation vote, r1525)
#   RUSTLE_DUP_ADMIT=1             pass 1: coordinate-duplicate spliced records count toward the floor for a chain that would otherwise
#                                  fall under it (r1524: rescues 2-read family chains; genome-wide too broad, 68,749 chains at 1% CAT match)
#   RUSTLE_CHECK_SPLICE=1          copy_assign --check-splice (FLAIR --check_splice on each spliced read)
#   RUSTLE_TERMINAL_EXON_ANCHOR=N  copy_assign --terminal-exon-anchor N (bp; FLAIR --stringent's end rule is 25)
#   RUSTLE_EXTEND_TRANSFRAGS=1     copy_assign --extend-transfrags (StringTie-style greedy extension, no flow)
#   RUSTLE_TRANSFER_ENDS=1         copy_assign --transfer-ends (ends from reads sharing the first / last junctions)
#   RUSTLE_READ_ISOFORM_K=K        copy_assign --read-isoform-k K (also admit a chain whose every junction >= K reads carry;
#                                  2 = StringTie's single-read choice, register 1267; buffered reader, slower)
case "${RUSTLE_TRANSFER_CREDIT:-}" in ""|split) ;; full|best) POLISH="$POLISH --transfer-credit $RUSTLE_TRANSFER_CREDIT" ;;
  *) echo "[rustle_pipeline] RUSTLE_TRANSFER_CREDIT must be split, full or best (got '$RUSTLE_TRANSFER_CREDIT')" >&2; exit 2 ;; esac
case "${RUSTLE_READ_JUNCTION_SNAP:-}" in ""|0) ;; *[!0-9]*) echo "[rustle_pipeline] RUSTLE_READ_JUNCTION_SNAP must be a whole number of bp (got '$RUSTLE_READ_JUNCTION_SNAP')" >&2; exit 2 ;;
  *) POLISH="$POLISH --read-junction-snap $RUSTLE_READ_JUNCTION_SNAP" ;; esac
case "${RUSTLE_MIN_CHAIN_READS:-}" in ""|0) ;; *[!0-9]*) echo "[rustle_pipeline] RUSTLE_MIN_CHAIN_READS must be a whole number (got '$RUSTLE_MIN_CHAIN_READS')" >&2; exit 2 ;;
  *) POLISH="$POLISH --min-chain-reads $RUSTLE_MIN_CHAIN_READS" ;; esac
case "${RUSTLE_CHECK_SPLICE:-}" in ""|0) ;; 1) POLISH="$POLISH --check-splice" ;;
  *) echo "[rustle_pipeline] RUSTLE_CHECK_SPLICE must be 0 or 1 (got '$RUSTLE_CHECK_SPLICE')" >&2; exit 2 ;; esac
case "${RUSTLE_TERMINAL_EXON_ANCHOR:-}" in ""|0) ;; *[!0-9]*) echo "[rustle_pipeline] RUSTLE_TERMINAL_EXON_ANCHOR must be a whole number of bp (got '$RUSTLE_TERMINAL_EXON_ANCHOR')" >&2; exit 2 ;;
  *) POLISH="$POLISH --terminal-exon-anchor $RUSTLE_TERMINAL_EXON_ANCHOR" ;; esac
case "${RUSTLE_EXTEND_TRANSFRAGS:-}" in ""|0) ;; 1) POLISH="$POLISH --extend-transfrags" ;;
  *) echo "[rustle_pipeline] RUSTLE_EXTEND_TRANSFRAGS must be 0 or 1 (got '$RUSTLE_EXTEND_TRANSFRAGS')" >&2; exit 2 ;; esac
case "${RUSTLE_TRANSFER_ENDS:-}" in ""|0) ;; 1) POLISH="$POLISH --transfer-ends" ;;
  *) echo "[rustle_pipeline] RUSTLE_TRANSFER_ENDS must be 0 or 1 (got '$RUSTLE_TRANSFER_ENDS')" >&2; exit 2 ;; esac
case "${RUSTLE_READ_ISOFORM_K:-}" in ""|0) ;; *[!0-9]*) echo "[rustle_pipeline] RUSTLE_READ_ISOFORM_K must be a whole number (got '$RUSTLE_READ_ISOFORM_K')" >&2; exit 2 ;;
  *) POLISH="$POLISH --read-isoform-k $RUSTLE_READ_ISOFORM_K" ;; esac
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
# RUSTLE_BRIDGE_REGROUP=off|f1|f1v2|f1units (unset = f1v2, THE DEFAULT since 2026-09-29 by the user's decision on the F1v2 held-out
# Outcome and its family-level side result, docs/archive/2026-09/PREREG_o1_cover_growth_2026-09-29.md; off = the 2026-09-25 products,
# byte for byte): copy_assign --bridge-regroup, passed explicitly by `assemble`. A transcript that is the only link
# between two pieces of its gene_id, when reads end at a PAS inside its intron and other reads start there at their own
# promoter (f1v2: and it carries fewer reads than each piece), becomes a relation `<gene_id>.fus<k>` with `fusion_of`,
# and the pieces split as RUSTLE_GTF_REGROUP splits them (so the two are exclusive). With f1 or f1v2, `assemble` also
# writes PREFIX.families.gtf, the GTF without the bridges, and `families` reads it, as the held-out runs did
# (docs/archive/2026-09/PREREG_f1_bridge_locus_2026-09-28.md, docs/archive/2026-09/PREREG_f1v2_readshare_2026-09-29.md), as does `flag`'s scan of the
# de novo loci (without --gff): a bridge is a relation, never a locus to scan. See its --help.
# BRIDGE_MODE is the one derivation every stage uses (assemble's flag, families' and flag's input, the guard below).
# f1units (OPT-IN, docs/archive/2026-09/PREREG_container_units_v2_dev_2026-09-30.md Part C): F1's bridges without its read-share rule are
# kept in PREFIX.families.gtf as units instead of being left out (PREFIX.gtf is exactly f1's: only the transcripts that were cut are
# tagged); RUSTLE_BRIDGE_UNITS_LIST names the cuts instead of F1.
BRIDGE_MODE=${RUSTLE_BRIDGE_REGROUP:-f1v2}
FAM_GTF=$OUT.gtf
case "$BRIDGE_MODE" in
  off) ;;
  f1|f1v2|f1units)
    [ "${RUSTLE_GTF_REGROUP:-0}" = 0 ] || { echo "[rustle_pipeline] RUSTLE_BRIDGE_REGROUP=$BRIDGE_MODE (unset = f1v2) already splits every gene_id as RUSTLE_GTF_REGROUP does: set RUSTLE_BRIDGE_REGROUP=off with RUSTLE_GTF_REGROUP=1" >&2; exit 2; }
    FAM_GTF=$OUT.families.gtf ;;
  *) echo "[rustle_pipeline] RUSTLE_BRIDGE_REGROUP must be off, f1, f1v2 or f1units (got '$RUSTLE_BRIDGE_REGROUP')" >&2; exit 2 ;;
esac
POLISH="$POLISH --bridge-regroup $BRIDGE_MODE"
case "${RUSTLE_BRIDGE_UNITS_LIST:-}" in
  "") ;;
  *[[:space:]]*) echo "[rustle_pipeline] RUSTLE_BRIDGE_UNITS_LIST must not contain whitespace (got '$RUSTLE_BRIDGE_UNITS_LIST')" >&2; exit 2 ;;
  *) [ "$BRIDGE_MODE" = f1units ] || { echo "[rustle_pipeline] RUSTLE_BRIDGE_UNITS_LIST names the units of RUSTLE_BRIDGE_REGROUP=f1units (got '$BRIDGE_MODE')" >&2; exit 2; }
     [ -s "$RUSTLE_BRIDGE_UNITS_LIST" ] || { echo "[rustle_pipeline] RUSTLE_BRIDGE_UNITS_LIST=$RUSTLE_BRIDGE_UNITS_LIST is missing or empty" >&2; exit 2; }
     POLISH="$POLISH --bridge-units-list $RUSTLE_BRIDGE_UNITS_LIST" ;;
esac
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
# of each family member's extra pieces (docs/archive/2026-09/PREREG_fusion_container_sim_2026-09-28.md §1): every clustered locus's
# all-transcript exon blocks, core (aligned to an exon base of another member of its family) or accessory, and the
# other families each accessory block aligns to -> PREFIX.fam.container.tsv / .container_relations.tsv /
# .container_summary.tsv. It never changes a family. See its --help
FAM_EXTRA=()
case "${RUSTLE_FAMILY_CONTAINER:-}" in
  ""|0) ;;
  1) FAM_EXTRA=(--emit-container) ;;
  *) echo "[rustle_pipeline] RUSTLE_FAMILY_CONTAINER must be 0 or 1 (got '$RUSTLE_FAMILY_CONTAINER')" >&2; exit 2 ;;
esac
# RUSTLE_FAMILY_RELATIONS=1 (opt-in; unset or 0 = off, the same command): mcl_families --emit-relations, the relation
# record of every unit-split transcript (RUSTLE_BRIDGE_REGROUP=f1units) and the members of each family by locus ->
# PREFIX.fam.relations.tsv / PREFIX.fam.members_by_locus.tsv (docs/archive/2026-09/PREREG_container_units_v2_dev_2026-09-30.md Part C).
# It never changes a family. Without units the relations table is empty and every locus is a `whole` member.
FAM_REL=()
case "${RUSTLE_FAMILY_RELATIONS:-}" in
  ""|0) ;;
  1) FAM_REL=(--emit-relations) ;;
  *) echo "[rustle_pipeline] RUSTLE_FAMILY_RELATIONS must be 0 or 1 (got '$RUSTLE_FAMILY_RELATIONS')" >&2; exit 2 ;;
esac
say() { echo "[rustle_pipeline] $(date +%H:%M:%S) $*" >&2; }
# The de novo loci a stage reads ($FAM_GTF: PREFIX.families.gtf when BRIDGE_MODE is f1, f1v2 or f1units, else PREFIX.gtf) must
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
  # --contig NAME[,NAME...]: these contigs only (copy_assign --region NAME:0-LEN, or --regions FILE for several, LEN from the FASTA
  # index), in place of --genome-wide; the products equal the genome-wide run's on those contigs (docs/archive/2026-09/CONTAINER_HEADROOM_2026-09-30.md gates)
  local scope=(--genome-wide)
  if [ -n "$CONTIG" ]; then
    case "$CONTIG" in *$'\n'*|*$'\r'*) echo "[rustle_pipeline] --contig: a newline or carriage return in the value (the list is comma separated)" >&2; exit 2;; esac
    case ",$CONTIG," in *,,*) echo "[rustle_pipeline] --contig '$CONTIG': an empty name in the list" >&2; exit 2;; esac
    local items=() c clen regions="" nctg=0 seen=","
    IFS=, read -ra items <<< "$CONTIG"
    for c in "${items[@]}"; do
      case "$c" in
        *:*|*[[:space:]]*) echo "[rustle_pipeline] --contig '$c': a name with ':' or whitespace cannot be a region (copy_assign splits at the first ':')" >&2; exit 2;;
      esac
      case "$seen" in *",$c,"*) echo "[rustle_pipeline] --contig lists $c twice" >&2; exit 2;; esac
      seen+="$c,"
      clen=""
      if [ -s "$FASTA.fai" ]; then clen=$(awk -F'\t' -v c="$c" '($1 "")==(c ""){print $2; exit}' "$FASTA.fai"); fi   # as strings: awk would equate 01, 1.0 and 1e0
      [ -n "$clen" ] || { echo "[rustle_pipeline] --contig $c: no such contig in $FASTA.fai (samtools faidx writes it)" >&2; exit 2; }
      regions+="$c:0-$clen"$'\n'; nctg=$((nctg + 1))
    done
    if [ "$nctg" = 1 ]; then
      scope=(--region "${regions%$'\n'}")
      say "assemble: contig $CONTIG only (--region ${regions%$'\n'})"
    else
      printf '%s' "$regions" > "$OUT.assemble_regions.txt"
      scope=(--regions "$OUT.assemble_regions.txt")
      say "assemble: contigs $CONTIG only ($nctg regions in $OUT.assemble_regions.txt)"
    fi
  fi
  # the seeding pool (the SEEDING POOL note above): the three RUSTLE_GTF_SECONDARY* variables the binary reads are set, or
  # unset, for every pool, so a variable exported in the shell cannot change the pool the flag names
  local seed_env=(-u RUSTLE_GTF_SECONDARY -u RUSTLE_GTF_SECONDARY_AS_RATIO -u RUSTLE_GTF_SECONDARY_AS_TABLE)
  if [ "$SEED_POOL" = good ]; then
    say "assemble: seeding pool: good (AS >= $SEED_RATIO x the molecule's genome-wide best)"
    # The table depends on the BAM alone, so by default it is SHARED by every run on that BAM: BAM.as_table.tsv next to
    # the BAM (2026-11-07; the pass is ~200 s of I/O on a 96 GB BAM, register 1567), or PREFIX.molecules.tsv when the BAM's
    # directory is not writable. --as-table FILE names it instead. Reused when its header names this BAM and it is not
    # older than the BAM; otherwise rebuilt under a temporary name and renamed into place (a concurrent run never reads
    # a half-written table; rename keeps the inode and mtime the .asbin sidecar is keyed on).
    local table=${AS_TABLE:-}
    as_table_ok() { head -1 "$1" 2>/dev/null | grep -qF "bam=$(readlink -f "$BAM")$(printf '\t')bam_bytes=$(stat -Lc %s "$BAM")$(printf '\t')" && [ ! "$1" -ot "$BAM" ]; }
    if [ -z "$table" ]; then
      # next to the BAM when a valid table is there or the directory takes a new file (-w alone misses a full file quota);
      # otherwise PREFIX.molecules.tsv, which the check below then reuses like any other table
      local probe="$BAM.as_table.probe.$$"
      table="$OUT.molecules.tsv"
      if as_table_ok "$BAM.as_table.tsv" || { : > "$probe"; } 2>/dev/null; then table="$BAM.as_table.tsv"; fi
      rm -f "$probe"
    fi
    if as_table_ok "$table"; then
      say "assemble: genome-wide best-AS table: reusing $table"
    else   # absent, truncated, from another BAM, or older than the BAM
      say "assemble: genome-wide best-AS table -> $table (one pass over the BAM)"
      local tmp="$table.tmp.$$"
      if "$BIN/as_table" --bam "$BAM" --out "$tmp" --threads "$THREADS" > "$OUT.as_table.log" 2>&1 && [ -s "$tmp" ]; then
        [ -e "$tmp.asbin" ] && mv -f "$tmp.asbin" "$table.asbin"
        mv -f "$tmp" "$table"
      else
        rm -f "$tmp" "$tmp.asbin"
        echo "[rustle_pipeline] as_table failed (see $OUT.as_table.log)" >&2; exit 3
      fi
    fi
    # PREFIX.molecules.tsv stays the per-run name later tools read (tools/rustle_reassemble.sh): a link to the shared table
    if [ "$(readlink -f "$table")" != "$(readlink -f "$OUT.molecules.tsv")" ]; then
      ln -sfn "$(readlink -f "$table")" "$OUT.molecules.tsv"
      ln -sfn "$(readlink -f "$table").asbin" "$OUT.molecules.tsv.asbin"
    fi
    seed_env=(RUSTLE_GTF_SECONDARY=1 RUSTLE_GTF_SECONDARY_AS_RATIO="$SEED_RATIO" "RUSTLE_GTF_SECONDARY_AS_TABLE=$table")
  elif [ "$SEED_POOL" = all ]; then
    say "assemble: seeding pool: all (every secondary alignment, no filter)"
    seed_env=(-u RUSTLE_GTF_SECONDARY_AS_RATIO -u RUSTLE_GTF_SECONDARY_AS_TABLE RUSTLE_GTF_SECONDARY=1)
  else
    say "assemble: seeding pool: primary (primary alignments only)"
  fi
  # r1525 PASS (2026-11-04): single-exon models take their strand from the read-orientation vote (default on here;
  # RUSTLE_READ_STRAND=0 restores the '+' placeholder). A119b: spliced chains identical, CAT single-exon gene hits 1,465 -> 1,490,
  # complete truth families 204 -> 209, 3,031 placeholder-strand models now shadowed by their same-strand spliced model.
  # After the pool branches on purpose: they reassign seed_env wholesale (the first ship, b8d29f0d, was inert for the good pool).
  seed_env+=("RUSTLE_READ_STRAND=${RUSTLE_READ_STRAND:-1}")
  # the exact flags of this run, for tools/rustle_reassemble.sh's pass 2 (it must not re-derive them from the defaults)
  printf '%s\n' "$POLISH" > "$OUT.assemble.polish"
  # ... and its seeding environment and contig scope, one word per line (pass 2 must seed and scope as pass 1 did)
  printf '%s\n' "${seed_env[@]}" > "$OUT.assemble.seed_env"
  printf '%s\n' "${scope[@]}" > "$OUT.assemble.scope"
  # PREFIX.final.gtf (2026-11-02, registers 1487-1491): the assembly with the three end layers (median 5' ends, 3' UTR extension,
  # single-exon models from end clusters), each PASSED on held-out chromosomes, but each changed the families stage when it read
  # them (r1490 / r1491). So the layered assembly is a SECOND run, in parallel, written only as PREFIX.final.gtf: PREFIX.gtf and
  # PREFIX.families.gtf (every later stage's input) are byte-identical to a run without it. RUSTLE_ASSEMBLY_LAYERS=off skips it.
  local layers=()
  case "${RUSTLE_ASSEMBLY_LAYERS:-on}" in
    on) # a layer the environment already set (RUSTLE_POLISH_TSS_END / _TES_EXTEND / RUSTLE_MONO_MODEL) is in $POLISH for both
        # runs, and a second copy of the flag stops clap: the user's value wins (audit 2026-11-07, fourth pass)
        case "$POLISH" in *--polish-tss-end*) ;; *) layers+=(--polish-tss-end median) ;; esac
        case "$POLISH" in *--polish-tes-extend*) ;; *) layers+=(--polish-tes-extend on) ;; esac
        case "$POLISH" in *--mono-model*) ;; *) layers+=(--mono-model cluster) ;; esac
        # RUSTLE_FINAL_SUBCHAIN=drop|off (register 1498 PASS; unset = drop): --polish-subchain drop as a fourth layer of PREFIX.final.gtf only
        case "${RUSTLE_FINAL_SUBCHAIN:-drop}" in
          drop) case "$POLISH" in *--polish-subchain*) ;; *) layers+=(--polish-subchain drop) ;; esac ;;
          off) ;;
          *) echo "[rustle_pipeline] RUSTLE_FINAL_SUBCHAIN must be drop or off (got '$RUSTLE_FINAL_SUBCHAIN')" >&2; exit 2 ;;
        esac ;;
    off) ;;
    *) echo "[rustle_pipeline] RUSTLE_ASSEMBLY_LAYERS must be on or off (got '$RUSTLE_ASSEMBLY_LAYERS')" >&2; exit 2 ;;
  esac
  # the three end layers read the --polish-tss evidence pass: with RUSTLE_POLISH_TSS=off only the sub-chain layer is kept
  # (it needs no --polish-tss; round-23 audit), and with no layer left there is no layered assembly to build
  if [ "${#layers[@]}" -gt 0 ] && [ "$POLISH_TSS" = off ]; then
    local kept=() i
    for (( i=0; i<${#layers[@]}; i+=2 )); do [ "${layers[i]}" = --polish-subchain ] && kept+=("${layers[i]}" "${layers[i+1]}"); done
    if [ "${#kept[@]}" -gt 0 ]; then say "assemble: RUSTLE_POLISH_TSS=off, so $OUT.final.gtf has only ${kept[*]} (the end layers need --polish-tss)"
    else say "assemble: RUSTLE_POLISH_TSS=off, so no layered $OUT.final.gtf (its end layers need --polish-tss); RUSTLE_ASSEMBLY_LAYERS=off silences this"; fi
    layers=("${kept[@]}")
  fi
  # r1564 PASS (2026-11-07): contigs swept in parallel (--region-threads; low-mem mode sets 1 = serial), A119b assemble
  # 2,588 -> 1,685 s, six products byte-identical, RSS 16 GB. The read budget is 20M records here, not assign's 4M:
  # at 4M chr1/chr2/chr13 (5-6M records each) run alone and the sweep stays near serial.
  local fin_pid=""
  rm -f "$OUT.final.gtf" "$OUT".final_run.*
  if [ "${#layers[@]}" -gt 0 ]; then
    env "${seed_env[@]}" RUSTLE_REGION_READ_BUDGET=20000000 "$BIN/copy_assign" --assemble-only --region-threads "${RUSTLE_REGION_THREADS:-8}" "${scope[@]}" --assembly-junctions strict $POLISH "${layers[@]}" --gtf-tpm \
      --bam "$BAM" --fasta "$FASTA" --out "$OUT.final_run" > "$OUT.final.log" 2>&1 &
    fin_pid=$!
  fi
  # set -e ends the stage when the base run fails: the background layered run must not outlive it holding memory and threads
  local rc=0
  env "${seed_env[@]}" RUSTLE_REGION_READ_BUDGET=20000000 "$BIN/copy_assign" --assemble-only --region-threads "${RUSTLE_REGION_THREADS:-8}" "${scope[@]}" --assembly-junctions strict $POLISH --gtf-tpm \
    --bam "$BAM" --fasta "$FASTA" --out "$OUT" > "$OUT.assemble.log" 2>&1 || rc=$?
  if [ "$rc" != 0 ]; then
    if [ -n "$fin_pid" ]; then kill "$fin_pid" 2>/dev/null || true; wait "$fin_pid" 2>/dev/null || true; rm -f "$OUT".final_run.*; fi
    echo "[rustle_pipeline] assemble failed (exit $rc, see $OUT.assemble.log)" >&2; exit "$rc"
  fi
  if [ -n "$fin_pid" ]; then
    if wait "$fin_pid" && [ -s "$OUT.final_run.gtf" ]; then
      mv -f "$OUT.final_run.gtf" "$OUT.final.gtf"
      say "assemble: $(awk -F'\t' '$3=="transcript"' "$OUT.final.gtf" | wc -l) transcripts in the layered assembly $OUT.final.gtf (log $OUT.final.log)"
    else
      say "assemble: WARNING the layered assembly failed (see $OUT.final.log); $OUT.gtf and the families input are unaffected"
    fi
    rm -f "$OUT".final_run.*
  fi
  # the good pool is only the good pool if the binary used the table: an empty table is ignored and the region-local rule applies, a table
  # from another BAM keeps every secondary it does not know; either way the run is wider than --seed-pool good says
  if [ "$SEED_POOL" = good ] && grep -qE 'table IGNORED|the table may come from a different BAM' "$OUT.assemble.log"; then
    echo "[rustle_pipeline] --seed-pool good: copy_assign did not apply the best-AS table ($(grep -m1 -E 'table IGNORED|the table may come from a different BAM' "$OUT.assemble.log" | cut -c1-300)); the assembly is not the good pool: delete $table and $table.asbin, then run again (an absent table is rebuilt; see $OUT.assemble.log)" >&2; exit 3
  fi
  # copy_assign writes PREFIX.families.gtf and PREFIX.gtf a few ms apart; on a filesystem with sub-second mtimes fam_gtf_guard's
  # `-ot` can then call it older (seen on NFS, 2026-10-29). Stamp it with PREFIX.gtf's mtime: same assembly, same age.
  [ "$FAM_GTF" = "$OUT.gtf" ] || [ ! -s "$FAM_GTF" ] || touch -r "$OUT.gtf" "$FAM_GTF"
  say "assemble: $(awk -F'\t' '$3=="transcript"' "$OUT.gtf" | wc -l) transcripts"
  if [ "$BRIDGE_MODE" = f1units ]; then
    say "assemble: $(grep -c 'fusion_of "' "$OUT.gtf") bridge transcripts cut into $(awk -F'\t' '$3=="transcript"' "$FAM_GTF" | grep -c 'fusion_unit "') units in the families input $FAM_GTF ($OUT.bridge_units.tsv)"
  elif [ "$FAM_GTF" != "$OUT.gtf" ]; then
    say "assemble: $(grep -c 'fusion_of "' "$OUT.gtf") bridge transcripts kept as fusion_of relations; families input $FAM_GTF"
  fi
}
# families: the de novo families AND their copy table (--emit-units with --from-gtf: PREFIX.fam.copies.tsv/.fa/.regions,
# the gw_family_catalog copies contract, one copy per member locus = its representative transcript and its spliced exon
# sum). clusters.tsv / loci.* are byte-identical to a run without it (cmp-checked 2026-09-25,
# docs/archive/2026-09/PREREG_families_copy_table_2026-09-25.md); the copy table is a new product. A binary older than the copy table
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
  if [ ${#FAM_REL[@]} -gt 0 ] && [[ "$help" != *'--emit-relations'* ]]; then
    echo "[rustle_pipeline] RUSTLE_FAMILY_RELATIONS=1 but $BIN/mcl_families predates --emit-relations; rebuild it" >&2; exit 2
  fi
  # RUSTLE_REPRESENTATIVE=most-reads|most-junctions (unset = nothing passed: the binary's default, most-reads, the earlier
  # products byte for byte; `most-reads` names that default): mcl_families --representative, which transcript represents each
  # de novo locus (docs/archive/2026-10/PREREG_locus_representative_rule_2026-10-04.md). Read here, in this stage's code alone, so that
  # figures/samples.py's driver_stage_code leaves the other stages' code hashes as they were.
  local rep=()
  case "${RUSTLE_REPRESENTATIVE:-}" in
    "") ;;
    most-reads|most-junctions)
      [[ "$help" == *'--representative'* ]] || { echo "[rustle_pipeline] RUSTLE_REPRESENTATIVE is set but $BIN/mcl_families predates --representative; rebuild it" >&2; exit 2; }
      rep=(--representative "$RUSTLE_REPRESENTATIVE") ;;
    *) echo "[rustle_pipeline] RUSTLE_REPRESENTATIVE must be most-reads or most-junctions (got '$RUSTLE_REPRESENTATIVE')" >&2; exit 2 ;;
  esac
  # Edge rule. DEFAULT since 2026-11-09 (user decision; row 1593 PASS on HG002 held-out): `share` = an edge needs one record on
  # which >= half of the fewer-exon locus's exons match an exon of the other (`--min-exon-share`, RUSTLE_MIN_EXON_SHARE,
  # default 0.5), with the coverage gate off (--min-shared-exon-frac 0.0) and MCL inflation RUSTLE_INFLATION (default 2.0).
  # HG002: plain .545/.396 vs the coverage rule's .466/.345; with RUSTLE_SEDEF_BED .703/.355 vs .719/.321.
  # RUSTLE_FAMILY_EDGE_RULE=coverage = the rule before 2026-11-09: RUSTLE_MIN_SHARED_EXON_FRAC=F (r1528; default 0.60) and the
  # binary's inflation (2.8); setting RUSTLE_MIN_SHARED_EXON_FRAC alone also selects it.
  local edge=()
  case "${RUSTLE_FAMILY_EDGE_RULE:-$([ -n "${RUSTLE_MIN_SHARED_EXON_FRAC:-}" ] && echo coverage || echo share)}" in
    share)
      [[ "$help" == *'--min-exon-share'* ]] || { echo "[rustle_pipeline] the default edge rule needs an mcl_families with --min-exon-share; rebuild it (or RUSTLE_FAMILY_EDGE_RULE=coverage)" >&2; exit 2; }
      edge=(--min-shared-exon-frac 0.0 --min-exon-share "${RUSTLE_MIN_EXON_SHARE:-0.5}" --inflation "${RUSTLE_INFLATION:-2.0}") ;;
    coverage)
      case "${RUSTLE_MIN_SHARED_EXON_FRAC:-0.60}" in
        0.[0-9]*|1|1.0) ;;
        *) echo "[rustle_pipeline] RUSTLE_MIN_SHARED_EXON_FRAC must be a fraction 0-1 (got '$RUSTLE_MIN_SHARED_EXON_FRAC')" >&2; exit 2 ;;
      esac
      edge=(--min-shared-exon-frac "${RUSTLE_MIN_SHARED_EXON_FRAC:-0.60}")
      [ -n "${RUSTLE_INFLATION:-}" ] && edge+=(--inflation "$RUSTLE_INFLATION") ;;
    *) echo "[rustle_pipeline] RUSTLE_FAMILY_EDGE_RULE must be share or coverage (got '$RUSTLE_FAMILY_EDGE_RULE')" >&2; exit 2 ;;
  esac
  say "families: edge rule ${edge[*]}"
  # RUSTLE_SEDEF_BED=FILE (SEDEF pairs BED of THIS genome; human CHM13 and gorilla have one, docs/DATA.md): SD damping ON at
  # RUSTLE_SD_DAMP (default 0.25, row 1580 PASS on HG002 held-out, default per user 2026-11-09; 1 = off). Unset = no damping,
  # the same command as before. The factor is always passed, so an older binary cannot fall back to its own default.
  local sd=()
  # the genome's own SEDEF BED as a sidecar, FASTA.sedef.bed, is used when RUSTLE_SEDEF_BED is unset (round-4 audit: SD damping
  # is the default for genomes that have SEDEF, and an unset variable used to drop it silently); RUSTLE_SEDEF_BED= (empty) turns it off
  if [ -z "${RUSTLE_SEDEF_BED+x}" ] && [ -s "$FASTA.sedef.bed" ]; then RUSTLE_SEDEF_BED=$FASTA.sedef.bed; fi
  if [ -n "${RUSTLE_SEDEF_BED:-}" ]; then
    [ -s "$RUSTLE_SEDEF_BED" ] || { echo "[rustle_pipeline] RUSTLE_SEDEF_BED=$RUSTLE_SEDEF_BED: no such non-empty file" >&2; exit 2; }
    [[ "$help" == *'--sd-damp-bed'* ]] || { echo "[rustle_pipeline] RUSTLE_SEDEF_BED is set but $BIN/mcl_families predates --sd-damp-bed; rebuild it" >&2; exit 2; }
    sd=(--sd-damp-bed "$RUSTLE_SEDEF_BED" --sd-damp "${RUSTLE_SD_DAMP:-0.25}")
    say "families: SD damping ${RUSTLE_SD_DAMP:-0.25} with $RUSTLE_SEDEF_BED"
  fi
  # an earlier run's copy table must not survive a run that writes none (a binary older than the copy table): assign would read it
  rm -f "$OUT".fam.copies.{tsv,fa,regions,merged.tsv}
  # threads: without --threads, every CPU of the job (the all-vs-all is ~98% of this stage; A119b 72 threads 410 s, the same
  # families as 12 threads, 2026-11-09). The products do not depend on it, so the settings record keeps $THREADS.
  local fth=$THREADS
  [ "$THREADS_SET" = 1 ] || fth=${SLURM_CPUS_PER_TASK:-$(nproc)}
  say "families: $fth threads"
  "$BIN/mcl_families" --from-gtf "$FAM_GTF" --fasta "$FASTA" --threads "$fth" \
    --min-exonic-bp 1 "${edge[@]}" "${copies[@]}" "${FAM_EXTRA[@]}" "${FAM_REL[@]}" "${FAM_COV[@]}" "${rep[@]}" "${sd[@]}" --out "$OUT.fam" > "$OUT.families.log" 2>&1
  # without minimap2-gate on PATH (or RUSTLE_MINIMAP2) every locus pair is aligned: same families, ~5x the wall (audit round 1)
  if grep -q 'has no --gate' "$OUT.families.log"; then
    say "families: WARNING no minimap2-gate on PATH: every pair was aligned (~5x slower); build it with tools/minimap2_gate/build.sh"
  fi
  say "families: $(awk 'NR>1' "$OUT.fam.clusters.tsv" | cut -f1 | sort -u | wc -l) clusters ($OUT.fam.clusters.tsv)"
  if [ ${#copies[@]} -gt 0 ]; then
    say "families: $(awk 'NR>1' "$OUT.fam.copies.tsv" | wc -l) copies (locus representatives) in $(awk 'NR>1' "$OUT.fam.copies.tsv" | cut -f1 | sort -u | wc -l) families ($OUT.fam.copies.tsv)"
  fi
  if [ ${#FAM_REL[@]} -gt 0 ]; then
    say "families: relations $(awk 'NR>1' "$OUT.fam.relations.tsv" | wc -l) split transcripts, $(awk -F'\t' 'NR>1 && $14=="cover"' "$OUT.fam.relations.tsv" | wc -l) cover, $(awk 'NR>1' "$OUT.fam.members_by_locus.tsv" | wc -l) members by locus ($OUT.fam.relations.tsv, $OUT.fam.members_by_locus.tsv)"
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
# its pre-registered acceptance, Amendment 12, FAILED, docs/archive/2026-10/O3_CANDIDATES_ACCEPTANCE_2026-10-02.md; the re-run, Amendment 13,
# passed, but its no-deletion control, Amendment 14, failed, so the 2026-10-03 default flip was reverted,
# docs/archive/2026-10/O3_CANDIDATES_CONTROL_A14_2026-10-03.md — naming the stage runs it, `all` runs it only with --candidates, and assign
# and flag use its products only with --candidates): candidate_copies turns each family's read net (reads with a record on its
# copies, plus the unmapped and the poorly placed reads >= 300 bp that align, map-ont, over >= 50% of their length at
# de <= 0.20 to the run's net reads or copies: Amendment 13b) into read clusters at --delta, one consensus per cluster on its
# structurally central member (Amendments 13d/13e), and candidate copies (clusters beyond delta of every reference locus,
# merged by the significance test, flagged with >= 6 reads), each represented by the exon union of its clusters ->
# PREFIX.cand.{candidates.tsv,contigs.fa,nets.fa,...}. With a flagged candidate: `utilities candidate-augment` writes PREFIX.aug.{fa,copies.tsv,copies.fa,regions.txt,
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
  rm -f "$OUT".aug.{fa,fa.fai,copies.tsv,copies.fa,regions.txt,families.txt,bam,bam.bai,mm2.log} "$OUT".cand.candidates.tsv
  if [ "$(awk 'NR > 1 && NF' "$OUT.fam.copies.tsv" | wc -l)" = 0 ]; then
    say "candidates: $OUT.fam.copies.tsv lists no copy (no family): nothing to do"; return 0
  fi
  say "candidates: candidate_copies on $OUT.fam.copies.tsv (delta $DELTA, at most $CAND_MAX reads per family)"
  "$BIN/candidate_copies" --bam "$BAM" --fasta "$FASTA" --copies "$OUT.fam.copies.tsv" --copies-fa "$OUT.fam.copies.fa" \
    --index "$INDEX" --delta "$DELTA" --max-reads "$CAND_MAX" --threads "$THREADS" --out "$OUT.cand" > "$OUT.candidates.log" 2>&1 \
    || { local rc=$?; say "candidates: candidate_copies failed (exit $rc), see $OUT.candidates.log"; exit "$rc"; }
  local n
  n=$(awk -F'\t' 'NR>1 && $5==1' "$OUT.cand.candidates.tsv" | wc -l)
  say "candidates: $n flagged candidate copies in $(awk -F'\t' 'NR>1 && $5==1' "$OUT.cand.candidates.tsv" | cut -f1 | sort -u | wc -l) families ($OUT.cand.candidates.tsv)"
  [ "$n" -gt 0 ] || return 0
  "$BIN/utilities" candidate-augment --fasta "$FASTA" --copies "$OUT.fam.copies.tsv" --copies-fa "$OUT.fam.copies.fa" \
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
  unset RUSTLE_READTHROUGH_JUNCTIONS   # assemble-only: copy_assign refuses it on every assign call (round-23 audit)
  # fam_regions all|only|skip [LIST]: copy_assign's --regions for the families of PREFIX.fam.copies.regions (`{fid}\t{chrom}:
  # {lo}-{hi}`, the copies' hull +- 5 kb per family and contig, as mcl_families writes it): every row, the rows of the
  # families LIST names (one id per line), or the other rows, MERGED per contig into disjoint intervals (an interval that
  # overlaps or touches the previous one joins it). copy_assign refuses overlapping --regions, and families' hulls
  # overlap and nest wherever families interleave (human chr16: 92 overlapping neighbours among 111 regions), so the
  # second column cannot be passed as it is. Every family's own interval lies inside exactly one merged region, which is
  # what binds it; its reads are still gathered around each copy (copy_assign's copy windows), not over the region.
  # `utilities candidate-augment` merges the candidate families' rows by the same rule (PREFIX.aug.regions.txt).
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
      --copies-fa "$OUT.aug.copies.fa" --only-families "$OUT.aug.families.txt" "${ASSIGN_FAST[@]}" "${INSPECT_ASSIGN[@]}" --out "$OUT.assign_cand" > "$OUT.assign_cand.log" 2>&1 \
      || { local rc=$?; say "assign: copy_assign failed on the candidate families (exit $rc), see $OUT.assign_cand.log; the likely cause: a real copy whose reads all realigned to its candidate, so it has no read in $OUT.aug.bam (copy_assign requires one on every real copy)"; exit "$rc"; }
    local runs=("$OUT.assign_cand")
    if [ "$n_rest" -gt 0 ]; then
      fam_regions skip "$OUT.aug.families.txt" > "$OUT.regions.txt"
      "$BIN/copy_assign" --bam "$BAM" --fasta "$FASTA" --regions "$OUT.regions.txt" --families "$OUT.fam.copies.tsv" \
        --copies-fa "$OUT.fam.copies.fa" --skip-families "$OUT.aug.families.txt" "${ASSIGN_FAST[@]}" "${INSPECT_ASSIGN[@]}" --out "$OUT.assign_rest" > "$OUT.assign_rest.log" 2>&1 \
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
    samtools view -H "$BAM" | awk '$1=="@SQ"{sub("SN:","",$2); sub("LN:","",$3); print $2":0-"$3}' > "$OUT.regions.txt"
    "$BIN/copy_assign" --bam "$BAM" --fasta "$FASTA" --regions "$OUT.regions.txt" \
      --families "$OUT.cat.copies.tsv" --copies-fa "$OUT.cat.copies.fa" "${ASSIGN_FAST[@]}" "${INSPECT_ASSIGN[@]}" --out "$OUT.assign" > "$OUT.assign.log" 2>&1
    say "assign: $(awk -F'\t' 'NR>1 && $4=="assigned"' "$OUT.assign.assignments.tsv" | wc -l) assigned rows of $(awk 'NR>1' "$OUT.assign.assignments.tsv" | wc -l) (one row per read x family)"
    return 0
  fi
  [ -s "$OUT.fam.copies.tsv" ] || { echo "assign needs $OUT.fam.copies.tsv (run the families stage with an mcl_families that writes the copy table), or pass --legacy-catalog" >&2; exit 2; }
  if [ "$(awk 'NR > 1 && NF' "$OUT.fam.copies.tsv" | wc -l)" = 0 ]; then
    rm -f "$OUT".assign.*.tsv "$OUT".assign.*.txt "$OUT".assign_cand.* "$OUT".assign_rest.*   # an earlier run's products
    say "assign: $OUT.fam.copies.tsv lists no copy (no family): nothing to assign"; return 0
  fi
  [ -s "$OUT.fam.copies.regions" ] || { echo "assign needs $OUT.fam.copies.regions (the families stage writes it beside the copy table)" >&2; exit 2; }
  # flagged candidates but no augmented families: the candidates stage did not finish; assigning without them would be silent
  if cand_ready && [ ! -s "$OUT.aug.families.txt" ] && awk -F'\t' 'NR > 1 && $5 == 1 { f = 1 } END { exit !f }' "$OUT.cand.candidates.tsv"; then
    echo "[rustle_pipeline] $OUT.cand.candidates.tsv has flagged candidates but $OUT.aug.families.txt is missing (the candidates stage did not finish): run it again, or drop --candidates" >&2; exit 2
  fi
  if cand_ready && [ -s "$OUT.aug.families.txt" ]; then
    # the candidates path applies none of the skip lists below: refused rather than silently ignored (round-10 audit)
    if [ "$MAX_WINDOW_RECORDS" -gt 0 ] || [ -n "${RUSTLE_ANNOTATED_ONLY:-}" ] || [ -n "${RUSTLE_EXCLUDE_BED:-}" ]; then
      echo "[rustle_pipeline] --candidates does not apply --max-window-records / --low-mem / RUSTLE_ANNOTATED_ONLY / RUSTLE_EXCLUDE_BED: drop one of them" >&2; exit 2
    fi
    assign_with_candidates; return 0
  fi
  if [ "$CANDIDATES" != 1 ] && [ -s "$OUT.aug.families.txt" ] && [ ! "$OUT.aug.families.txt" -ot "$OUT.fam.copies.tsv" ]; then
    say "assign: the candidates stage's products ($OUT.aug.*) are present but unused (the stage is opt-in, ruling R14): pass --candidates to assign the candidate families on them"
  fi
  say "assign: per-read copy assignment on $OUT.fam.copies.tsv"
  rm -f "$OUT".assign.*.tsv "$OUT".assign_cand.* "$OUT".assign_rest.*   # an earlier split run's tables
  local mem=()
  # --tied-prefilter (byte-identical): each region is read light first and only the AS-tied molecules' records are loaded in full
  mem+=(--tied-prefilter)
  # --max-window-records N (--low-mem: 500,000): a family whose merged copy windows (copy +- 50 kb clipped to its region, what copy_assign reads) hold more
  # than N records copy_assign would keep (rDNA arrays, acrocentric short arms) is not assigned, and its regions are not swept: listed
  # with the count reached in PREFIX.assign.depth_skipped.tsv (`utilities window-depth`). Only reads such a family shares with others
  # can change: their placement on its copies is no longer seen.
  rm -f "$OUT.assign.depth_skipped.tsv" "$OUT.assign.depth_skipped.txt"
  if [ "$MAX_WINDOW_RECORDS" -gt 0 ]; then
    # clipped to the regions copy_assign sweeps (every family's: the skip lists only shrink them, so still an upper bound)
    fam_regions all > "$OUT.assign.depth_regions.txt"
    "$BIN/utilities" window-depth --bam "$BAM" --copies "$OUT.fam.copies.tsv" --max "$MAX_WINDOW_RECORDS" \
      --regions "$OUT.assign.depth_regions.txt" --out "$OUT.assign.depth_skipped.tsv" 2>> "$OUT.assign.depth.log"
    if [ -s "$OUT.assign.depth_skipped.tsv" ]; then
      cut -f1 "$OUT.assign.depth_skipped.tsv" > "$OUT.assign.depth_skipped.txt"   # --skip-families is passed once below
      say "assign: $(wc -l < "$OUT.assign.depth_skipped.tsv") famil(ies) with a copy window over $MAX_WINDOW_RECORDS records are not assigned ($OUT.assign.depth_skipped.tsv)"
    fi
  fi
  # RUSTLE_ANNOTATED_ONLY=GFF (opt-in, row 1551): families with no copy on a `gene` feature of GFF join the skip list
  # (A119b: 212 of 1,565, among them every chr13 hotspot family that sets the assign wall). Output-changing: not a default.
  rm -f "$OUT.assign.unannotated.txt"
  if [ -n "${RUSTLE_ANNOTATED_ONLY:-}" ]; then
    python3 "$RP_DRIVER_DIR/annotated_families.py" "$OUT.fam.copies.tsv" "$RUSTLE_ANNOTATED_ONLY" > "$OUT.assign.unannotated.txt" 2>> "$OUT.assign.depth.log" || exit 2
    say "assign: annotated-only: $(wc -l < "$OUT.assign.unannotated.txt") famil(ies) with no copy on a gene of $RUSTLE_ANNOTATED_ONLY are not assigned ($OUT.assign.unannotated.txt)"
  fi
  # RUSTLE_EXCLUDE_BED=BED (opt-in, row 1559): families whose copies mostly lie in a confounding element (T2T censat
  # satellites: the chr13 HSat1A array is 618 of A119b's 1,565 families and the assign wall) join the skip list.
  rm -f "$OUT.assign.excluded.txt"
  if [ -n "${RUSTLE_EXCLUDE_BED:-}" ]; then
    python3 "$RP_DRIVER_DIR/annotated_families.py" "$OUT.fam.copies.tsv" "$RUSTLE_EXCLUDE_BED" --majority-in > "$OUT.assign.excluded.txt" 2>> "$OUT.assign.depth.log" || exit 2
    say "assign: exclude-bed: $(wc -l < "$OUT.assign.excluded.txt") famil(ies) mostly inside $RUSTLE_EXCLUDE_BED are not assigned ($OUT.assign.excluded.txt)"
  fi
  # the regions are rebuilt without the skipped families: a region left with no family would be read whole
  { for f in "$OUT.assign.depth_skipped.txt" "$OUT.assign.unannotated.txt" "$OUT.assign.excluded.txt"; do [ -s "$f" ] && cat "$f"; done; true; } | sort -u > "$OUT.assign.skipped.txt"
  if [ -s "$OUT.assign.skipped.txt" ]; then
    mem+=(--skip-families "$OUT.assign.skipped.txt")
    fam_regions skip "$OUT.assign.skipped.txt" > "$OUT.regions.txt"
  else
    fam_regions all > "$OUT.regions.txt"
  fi
  if [ ! -s "$OUT.regions.txt" ]; then
    say "assign: every family is on a skip list ($OUT.assign.skipped.txt): nothing to assign"; return 0
  fi
  "$BIN/copy_assign" --bam "$BAM" --fasta "$FASTA" --regions "$OUT.regions.txt" --families "$OUT.fam.copies.tsv" \
    --copies-fa "$OUT.fam.copies.fa" "${ASSIGN_FAST[@]}" "${mem[@]}" "${INSPECT_ASSIGN[@]}" --out "$OUT.assign" > "$OUT.assign.log" 2>&1 \
    || { local rc=$?; say "assign: copy_assign failed (exit $rc), see $OUT.assign.log"; exit "$rc"; }
  say "assign: $(awk -F'\t' 'NR>1 && $4=="assigned"' "$OUT.assign.assignments.tsv" | wc -l) assigned rows of $(awk 'NR>1' "$OUT.assign.assignments.tsv" | wc -l) (one row per read x family)"
}
# assign on the FAMILIES copy table (the merged stage's phase 3): the same byte-identical levers the genome-wide
# sweep pre-registers recommended — `--skip-poa-diagnostic` (the POA homology pass is purely diagnostic; measured
# 6.8x on the heaviest family) and `--region-threads` (cross-contig rayon pool; byte-identical by construction) —
# are ON here by default instead of opt-in. Products: PREFIX.assign_fam.* (the catalog assign keeps PREFIX.assign.*).
stage_assign_fam() {
  # this phase reads every family: the skip lists of the `assign` stage are not applied, so they are refused (round-10 audit)
  if [ "$MAX_WINDOW_RECORDS" -gt 0 ] || [ -n "${RUSTLE_ANNOTATED_ONLY:-}" ] || [ -n "${RUSTLE_EXCLUDE_BED:-}" ]; then
    echo "[rustle_pipeline] merged applies no --max-window-records / --low-mem / RUSTLE_ANNOTATED_ONLY / RUSTLE_EXCLUDE_BED: use \`all\` for them" >&2; exit 2
  fi
  say "assign: per-read copy assignment on $OUT.fam.copies (families table; --skip-poa-diagnostic --region-threads $THREADS)"
  samtools view -H "$BAM" | awk '$1=="@SQ"{sub("SN:","",$2); sub("LN:","",$3); print $2":0-"$3}' > "$OUT.regions.txt"
  # RUSTLE_READTHROUGH_JUNCTIONS is assemble-only (copy_assign refuses it here): dropped for this call alone, so the
  # merged.env record written after this phase still names it and a rerun resumes (round-24 audit)
  env -u RUSTLE_READTHROUGH_JUNCTIONS "$BIN/copy_assign" --bam "$BAM" --fasta "$FASTA" --regions "$OUT.regions.txt" \
    --families "$OUT.fam.copies.tsv" --copies-fa "$OUT.fam.copies.fa" \
    --skip-poa-diagnostic --region-threads "$THREADS" "${INSPECT_ASSIGN[@]}" --out "$OUT.assign_fam" > "$OUT.assign_fam.log" 2>&1
  say "assign: $(awk -F'\t' 'NR>1 && $4=="assigned"' "$OUT.assign_fam.assignments.tsv" | wc -l) assigned rows of $(awk 'NR>1' "$OUT.assign_fam.assignments.tsv" | wc -l) (one row per read x family)"
}
# The genome-wide entry point: assemble -> families -> assign(fam), one command. Resume: PREFIX.merged.env records
# the exact inputs/settings of a COMPLETED merged run (written only when all three phases finished); a phase is
# skipped only when that file matches AND the phase's product exists, so a stale product from another BAM/genome is
# never silently reused. A crashed run redoes the phases (the families all-vs-all replays from PREFIX.cache when
# the loci FASTA is unchanged — the expensive step is never repeated).
# seed_sec= keeps its two values for the default and primary pools (a finished run of either is still recognised); the all pool and
# a non-default width of the good pool add a line, so a finished run of one pool never stands in for another
merged_env() {
  printf 'bam=%s\nfasta=%s\nthreads=%s\nseed_sec=%s\nbridge=%s\n' "$(readlink -f "$BAM")" "$(readlink -f "$FASTA")" "$THREADS" \
    "$([ "$SEED_POOL" = primary ] && echo 0 || echo 1)" "$BRIDGE_MODE"
  if [ "$SEED_POOL" = all ]; then printf 'seed_pool=all\n'; fi
  if [ "$SEED_POOL" = good ] && [ "$SEED_RATIO" != 0.98 ]; then printf 'seed_ratio=%s\n' "$SEED_RATIO"; fi
  if [ -n "$CONTIG" ]; then printf 'contig=%s\n' "$CONTIG"; fi
  # the SEDEF sidecar the families stage picks up when RUSTLE_SEDEF_BED is unset (round-4 audit)
  if [ -z "${RUSTLE_SEDEF_BED+x}" ] && [ -s "$FASTA.sedef.bed" ]; then stat -L -c "sedef_sidecar %s %Y" "$FASTA.sedef.bed"; fi
  # every RUSTLE_* knob changes some stage's products: a resume must not skip a stage run with other values (audit 2026-11-07)
  env | grep '^RUSTLE_' | LC_ALL=C sort || true
  # and the binaries: a rebuild must rerun the stages it changes, not replay the old products (audit 2026-11-07, second pass)
  local b; for b in copy_assign mcl_families as_table utilities missing_copy_flag candidate_copies gw_family_catalog; do
    [ -n "${BIN:-}" ] && [ -e "$BIN/$b" ] && stat -L -c "bin $b %s %Y" "$BIN/$b"
  done; true
}
stage_merged() {
  local env_ok=0
  if [ -f "$OUT.merged.env" ] && [ "$(cat "$OUT.merged.env")" = "$(merged_env)" ]; then env_ok=1; fi
  # removed before any phase runs (a run that dies midway must not leave the record vouching for its products), and once one
  # phase runs every later phase runs too (its input changed), as in stage_all
  rm -f "$OUT.merged.env"
  if [ "$env_ok" = 1 ] && [ -s "$OUT.gtf" ]; then
    say "merged: phase 1 assemble SKIPPED ($OUT.gtf present, env matches)"
  else
    stage_assemble; env_ok=0
  fi
  if [ "$env_ok" = 1 ] && [ -s "$OUT.fam.copies.tsv" ]; then
    say "merged: phase 2 families SKIPPED ($OUT.fam.copies.tsv present, env matches)"
  else
    stage_families; env_ok=0
  fi
  if [ "$env_ok" = 1 ] && [ -s "$OUT.assign_fam.assignments.tsv" ]; then
    say "merged: phase 3 assign SKIPPED ($OUT.assign_fam.assignments.tsv present, env matches)"
  else
    stage_assign_fam
  fi
  merged_env > "$OUT.merged.env"
  say "merged: done (env recorded in $OUT.merged.env; phases skip on the next identical run)"
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
  # a FASTA as --index (the --low-mem profile) is indexed on the fly in INDEX_GB-Gbp parts; a prebuilt .mmi is used as it is
  local igb=()
  case "$INDEX" in *.mmi) ;; *) igb=(--index-gb "$INDEX_GB");; esac
  "$BIN/missing_copy_flag" --bam "$BAM" --fasta "$FASTA" --loci "$LOCI" --index "$INDEX" --threads "$THREADS" "${igb[@]}" \
    "${CONFIRM[@]}" "${FOREIGN[@]}" "${cand[@]}" --out "$OUT.flag" --from-scan "$OUT.flag_scan" > "$OUT.flag.log" 2>&1
  say "flag: $(grep -o 'verdicts: .*' "$OUT.flag.log")"
  if [ ${#cand[@]} -gt 0 ]; then say "flag: $(grep -o '[0-9]* flagged candidates with a nearest locus; [0-9]* rows name one' "$OUT.flag.log") (o3_candidate column)"; fi
}
# index: the splice index flag aligns its consensus sequences against (shared next to the FASTA, below). `-I INDEX_GB` Gbp per part: the default
# 8 gives one part for a mammalian genome (the figures' index, ~13 GB to hold); --index-gb 1 / --low-mem writes ~1-Gbp parts
# that flag loads one at a time (minimap2 --split-prefix merges their hits: near-identical, not byte-identical).
# The index depends on the FASTA, the part size and minimap2 alone, so it is SHARED by every run on that genome when it can be
# (round-2 audit: each run directory built its own ~13 GB copy): in RUSTLE_INDEX_DIR when set, else next to the FASTA unless
# the FASTA lies under $HOME (large files stay off the home quota), as <FASTA name>.splice.I<GB>G.mmi; otherwise, or when that
# directory takes no new file, PREFIX.splice.mmi. Its .key (FASTA path, size, mtime, part size, minimap2 version) must match
# to be reused; it is built under a temporary name and renamed, so a concurrent run never reads half of it.
index_key() { printf 'fasta=%s\nfasta_stat=%s\nindex_gb=%s\nminimap2=%s\n' "$(readlink -f "$FASTA")" "$(stat -Lc '%s %Y' "$FASTA")" "$INDEX_GB" "$("${RUSTLE_MINIMAP2:-minimap2}" --version 2>/dev/null)"; }
index_path() {
  local f d p probe; f=$(readlink -f "$FASTA"); d=${RUSTLE_INDEX_DIR:-$(dirname "$f")}
  if [ -z "${RUSTLE_INDEX_DIR:-}" ] && [[ "$f" == "$(readlink -f "$HOME")"/* ]]; then echo "$OUT.splice.mmi"; return; fi
  p="$d/$(basename "$f").splice.I${INDEX_GB}G.mmi"; probe="$p.probe.$$"
  if { [ -s "$p" ] && [ "$(cat "$p.key" 2>/dev/null)" = "$(index_key)" ]; } || { : > "$probe"; } 2>/dev/null; then rm -f "$probe"; echo "$p"
  else echo "$OUT.splice.mmi"; fi
}
stage_index() {
  local idx=${INDEX:-$(index_path)} key; key=$(index_key)
  if [ -s "$idx" ] && [ "$(cat "$idx.key" 2>/dev/null)" = "$key" ]; then say "index: $idx present"; return 0; fi
  rm -f "$idx.key"
  say "index: minimap2 -x splice -I ${INDEX_GB}G -d $idx"
  "${RUSTLE_MINIMAP2:-minimap2}" -x splice -I "${INDEX_GB}G" -t "$THREADS" -d "$idx.tmp.$$" "$FASTA" > "$OUT.index.log" 2>&1
  mv -f "$idx.tmp.$$" "$idx"
  printf '%s\n' "$key" > "$idx.key"
}
# stage_all_env STAGE: what a finished `all` stage was run with; a stage is skipped only when it matches (and its product is there)
stage_all_env() {
  merged_env
  printf 'candidates=%s\nlegacy=%s\ngff=%s\nindex=%s\nindex_gb=%s\nregion_threads=%s\npoa=%s\n' "$CANDIDATES" "$LEGACY_CATALOG" \
    "$GFF" "${INDEX:-$OUT.splice.mmi}" "$INDEX_GB" "${RUSTLE_REGION_THREADS:-8}" "${RUSTLE_POA_DIAGNOSTIC:-0}"
  printf 'max_window_records=%s\n' "$MAX_WINDOW_RECORDS"
  # settings only one stage reads go only into that stage's record (a flag-only option re-runs flag alone; round-23 audit),
  # and the driver itself everywhere (a default flipped here reaches no binary stat; round-10 audit)
  case "$1" in
    flag) printf 'confirm=%s\nforeign=%s\n' "${CONFIRM[*]-}" "${FOREIGN[*]-}";;
    candidates) printf 'delta=%s\ncand_max=%s\n' "$DELTA" "$CAND_MAX";;
    catalog|assign) printf 'inspect=%s\n' "$INSPECT";;
  esac
  [ -f "$0" ] && stat -L -c "driver %s %Y" "$0"; true
}
# stage_all: every stage as its own driver call under /usr/bin/time (PREFIX.<stage>.time: wall seconds, peak RSS of the stage's
# process tree, exit). Resume: a stage whose product exists and is not newer than its time record (a direct stage call rewrote
# it), whose time record says exit 0 and whose OWN settings record PREFIX.<stage>.env matches is skipped — until one stage runs,
# after which every later stage runs too. A stage's record is removed before it runs and written only when it succeeds, so a
# failed or killed stage is redone and the stages before it are not; `--from` writes no record for the stages it skips.
stage_all() {
  local plan=(assemble:"$OUT.gtf" families:"$OUT.fam.copies.tsv")
  if [ "$LEGACY_CATALOG" = 1 ]; then plan+=(catalog:"$OUT.cat.copies.tsv")
  elif [ "$CANDIDATES" = 1 ]; then plan+=(candidates:"$OUT.cand.candidates.tsv"); fi
  plan+=(assign:"$OUT.assign.assignments.tsv")
  local env_now ran=0 st prod args=("${ORIG_ARGS[@]}")
  # families, assign, index and flag all run minimap2: a missing one fails at start, not after the assembly (round-13 audit)
  command -v "${RUSTLE_MINIMAP2:-minimap2}" > /dev/null || { echo "[rustle_pipeline] no ${RUSTLE_MINIMAP2:-minimap2} on PATH (families, assign and flag run it); set RUSTLE_MINIMAP2 or PATH" >&2; exit 2; }
  if [ "$NO_FLAG" = 0 ]; then
    # no --index: the default profile builds the shared FASTA.splice.I<GB>G.mmi (index stage); --low-mem aligns against the FASTA itself,
    # indexed on the fly in INDEX_GB-Gbp parts (no index file on disk)
    if [ -z "$INDEX" ] && [ "$LOW_MEM" = 1 ]; then INDEX=$FASTA; args+=(--index "$FASTA")
    elif [ -z "$INDEX" ]; then INDEX=$(index_path); plan+=(index:"$INDEX"); args+=(--index "$INDEX"); fi
    plan+=(flag:"$OUT.flag.missing_copy.tsv")
  fi
  local started=1; [ -n "$FROM_STAGE" ] && started=0
  for item in "${plan[@]}"; do
    st=${item%%:*}; prod=${item#*:}; env_now=$(stage_all_env "$st")
    if [ "$started" = 0 ]; then
      if [ "$st" = "$FROM_STAGE" ]; then started=1; else say "all: $st skipped (--from $FROM_STAGE)"; continue; fi
    fi
    if [ "$RESUME" = 1 ] && [ "$ran" = 0 ] && [ -s "$prod" ] && [ -f "$OUT.$st.env" ] && [ "$(cat "$OUT.$st.env")" = "$env_now" ] \
       && grep -qx $'exit\t0' "$OUT.$st.time" 2>/dev/null && [ ! "$prod" -nt "$OUT.$st.time" ]; then
      say "all: $st SKIPPED ($prod present, settings unchanged)"; continue
    fi
    ran=1
    rm -f "$OUT.$st.env" "$OUT.all.env"   # all.env: the single record of builds before 2026-11-07
    say "all: $st"
    local rc=0
    if [ -x /usr/bin/time ]; then
      /usr/bin/time -f $'wall_s\t%e\npeak_rss_kb\t%M' -o "$OUT.$st.time.tmp" bash "$0" "$st" "${args[@]}" || rc=$?
    else
      # compute nodes without GNU time (r1543 e2e, 2026-11-06): wall only, no peak RSS
      local st0=$SECONDS; bash "$0" "$st" "${args[@]}" || rc=$?; printf 'wall_s\t%s\n' "$(( SECONDS - st0 ))" > "$OUT.$st.time.tmp"
    fi
    { cat "$OUT.$st.time.tmp" 2>/dev/null | grep -E '^(wall_s|peak_rss_kb)'; printf 'exit\t%s\n' "$rc"; } > "$OUT.$st.time"; rm -f "$OUT.$st.time.tmp"
    [ "$rc" = 0 ] || { say "all: $st FAILED (exit $rc); see $OUT.$st.time and the stage's log"; exit "$rc"; }
    printf '%s\n' "$env_now" > "$OUT.$st.env"
    if [ -n "$TO_STAGE" ] && [ "$st" = "$TO_STAGE" ]; then say "all: stopped after $st (--to $TO_STAGE); product $prod"; return 0; fi
  done
  say "all: done. stage, wall, peak memory:"
  for item in "${plan[@]}"; do
    st=${item%%:*}
    awk -F'\t' -v st="$st" '$1=="wall_s"{w=$2} $1=="peak_rss_kb"{m=$2} END{printf "[rustle_pipeline]   %-11s %8.0f s  %s\n", st, w, m=="" ? "  peak n/a (no GNU time)" : sprintf("%6.1f GB", m/1048576)}' "$OUT.$st.time" >&2
  done
}
case "$STAGE" in
  merged) stage_merged;;
  assemble) stage_assemble;; families) stage_families;; candidates) stage_candidates;; catalog) stage_catalog;;
  assign) stage_assign;; flag) stage_flag;; index) stage_index;;
  # --legacy-catalog: `assign` reads the legacy catalog, so `all` builds it; --candidates: `all` runs the (opt-in) candidates stage
  all) stage_all;;
  *) echo "unknown stage $STAGE" >&2; exit 2;;
esac
say "done"
