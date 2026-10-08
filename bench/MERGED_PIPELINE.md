# Merged genome-wide pipeline (`rustle_pipeline.sh merged`)

2026-10-04. The genome-wide entry point: `assemble -> families -> assign on the FAMILIES copy table`
in one driver command, byte-identical products, with the genome-wide levers ON by default and
per-phase resume. Motivation and stage dataflow: `bench/PERFORMANCE_AND_IO.md`; the extraction of the
families stage into `src/rustle/vg_family/fam_from_gtf.rs` (imported back by `mcl_families`, behavior
unchanged) is registered in `docs/MODULE_STATUS.md`.

## Why a driver-level merge (and not one process)

Explored and rejected after the Step-0 timings: `copy_assign`'s `main` is a ~3.9k-line function whose
two modes (assemble / assign) communicate through process-global state (`SKIP_READ_SEQUENCE`, env
vars set mid-run, `OnceLock` caches). Running both modes in one process risks silent cross-phase
leaks against a byte-identical requirement, while the measurable saving of the in-process handoffs
is seconds: the GTF/copies re-parse the merge would skip is ~1–5 s against a families stage of
~3.6k s (gorilla genome-wide, below). The merged stage therefore orchestrates the SAME binaries the
separate stages run, with three real changes:

1. **assign runs on the families copy table with the genome-wide levers ON**
   (`--skip-poa-diagnostic`, `--region-threads N`) — both cmp-proven byte-identical
   (`bench/PERFORMANCE_AND_IO.md`), both opt-in until now. The legacy catalog assign
   (`stage_assign`, on `gw_family_catalog` copies) is untouched.
2. **Per-phase resume**: `PREFIX.merged.env` records the exact inputs/settings of a COMPLETED merged
   run and is written only when all three phases finished; a phase is skipped only when the env
   matches AND its product exists. A crashed run redoes phases, but the families all-vs-all PAF
   replays from `PREFIX.cache` when the loci FASTA is unchanged (content-keyed, pinned) — the
   expensive step is never repeated.
3. **The PAF cache is populated as a side effect of the first run** — the dominant cost of any rerun
   (including parameter sweeps over downstream settings) becomes ~2 min of MCL + copies instead of
   ~1 h of minimap2.

What merging does NOT fix: the all-vs-all itself (minimap2, `-x asm20 -c -X -N 50 -p 0.1`, already
multi-threaded) is CPU-bound; on the 5-core box it is the genome-wide floor. For the human 96 GB
BAM this stage belongs on the second machine (docs/archive/2026-09/HANDOFF_SECOND_MACHINE_2026-09-30.md).

## Step 0 baseline — current pipeline, gorilla genome-wide (GGO_mm.bam 11.7 GB, mGorGor1)

One `--seed-secondaries` run (the completeness config; `as_table` pass included), then the stages
as the driver runs them today. Binaries: pre-extraction build (e163d955 working tree), `--threads 4`.

| stage | wall | peak RSS | notes |
|---|---|---|---|
| as_table (full BAM scan) | 76 s | 0.78 GB | seeds secondaries within 2% of genome-wide best AS |
| assemble (`--assemble-only --genome-wide`) | 112 s | see log | 26 contigs, streaming |
| families (`mcl_families --from-gtf --emit-units`) | **3,568 s** | — | minimap2 all-vs-all ~59 min of it (1 GB loci FASTA); MCL+copies ~2 min; 498 clusters / 1,607 copies |
| assign on fam copies, `--skip-poa-diagnostic --region-threads 4` | 232 s | 8.4 GB | 557 assigned / 9,353 rows (read × family) |
| **total** | **~62.8 min** | | **minimap2 = 90% of wall-clock** |

Assign is NOT the genome-wide bottleneck on gorilla once the diagnostic is skipped and regions are
parallel — the merged stage's defaults are cheap insurance, not the headline. The headline is
resume + cache: a second merged run costs assemble (~2 min) + families replay (~2 min) + assign
(~4 min).

## Verification (byte-identity)

### Step 1 — `mcl_families` extraction (`fam_from_gtf.rs`)

The NEW binary (post-extraction) rerun on the baseline's `gw.families.gtf` with the same `--out gw.fam`
prefix and the warm PAF cache reproduced every families product **byte-identical** (clusters.tsv,
copies.tsv/.fa/.regions/.merged.tsv, loci.gff3/.fa/.paf, params.tsv, loci.tsv; cmp, 2026-10-04). Stage wall
time on a warm cache: **25 s** (the 59-min minimap2 replays from `PREFIX.cache`; `bench/../merged_bench/cmp2/`).
`cargo test --release`: 1041 tests, 0 failed.

### Step 5 — merged stage vs staged pipeline, gorilla genome-wide

`rustle_pipeline.sh merged` on a fresh prefix with a warm PAF cache: **8m 24s total** (as_table 1:16,
assemble 1:52, families 0:40 with the PAF replay, assign 4:36), peak RSS 7.9 GB. Every product cmp'd
against the staged baseline: **all identical** — `gtf`, `families.gtf`, `molecules.tsv`,
`bridge_junctions.tsv`, `bridges.tsv`, `fam.clusters.tsv`, `fam.copies.tsv/.fa/.regions/.merged.tsv`,
`fam.loci.gff3/.fa/.paf`, `assignments.tsv` (557 assigned / 9,353 rows), `assign.quant.tsv`,
`assign.families.tsv`, `assign.famcn_readonly.tsv`, and `fam.params.tsv` after prefix normalization.
Two EXPECTED, content-neutral exceptions, both by design:
- `molecules.tsv.asbin` embeds the TSV's size/mtime/inode/device staleness identity
  (`AsTsvIdentity`, denovo_assemble.rs) — never equal across runs, like the driver's own staleness rule;
- `fam.params.tsv` rows `paf`/`gff` carry the output prefix (same for any two prefixes in the staged
  pipeline; normalized before cmp).

Resume: an identical re-invocation skips all three phases in 0.03 s (`PREFIX.merged.env` match).

### Step 5.3 — human chr20+chr21 slice (A119b), staged vs merged

Deep-data scale check (`A119b.t2t.bam` slice, 3.8 GB, 13,238 transcripts -> 6,705 loci, 126 MB loci FASTA;
`chm13v2.0.fa`). Staged path 5,007 s (as_table+assemble 31 s, **minimap2 all-vs-all ~79 min** — slower than
the ENTIRE gorilla genome's 59 min on 8x less sequence: deep data is disproportionately expensive at
`-N 50 -p 0.1`), assign ~3 min; merged run 4,570 s (its own cold all-vs-all; with a shared cache it
replays like Step 5). **All 14 products byte-identical** (gtf, families.gtf, molecules.tsv,
bridge_junctions.tsv, fam.clusters.tsv, fam.copies.tsv/.fa, fam.loci.gff3/.fa/.paf, fam.params.tsv
prefix-normalized, assignments.tsv, assign.quant.tsv). Human scaling: a whole-genome cold all-vs-all on
this box is a ~15-25 h class job — run it on the second machine (or overnight) once; every rerun then
replays from `PREFIX.cache`.

## How to run it

```sh
tools/rustle_pipeline.sh merged --bam READS.bam --fasta GENOME.fa --out run --threads 4
# rerun after any change: finished phases skip; the all-vs-all replays from run.cache when the loci FASTA is unchanged
```
