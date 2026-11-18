# Rustle

Rust tools for **multi-copy gene families in long-read RNA (IsoSeq/HiFi)**: define the families from the
reads and the genome, assign ambiguous reads to copies, and flag copies the reference does not contain.
Substrate: gorilla (mGorGor1) and human (CHM13) IsoSeq. This is a thesis codebase; the scientific record is
kept next to the code and is pre-registered.

## What it does

| stage | what | where |
|---|---|---|
| **family definition** | multi-copy gene families from the reads and the genome — de novo (reads only) and guided (annotation), with the goal of closing the gap between them | `docs/seeded_family_definition.md`, `mcl_families`, `gw_family_catalog` |
| **copy assignment** | each AS-tied multimapping read (a rival alignment within 10 AS of its best) to one of its very similar copies, or abstention: PSVs, junctions and divergence, never a 1/k split. A read minimap2 placed clearly is not O2's | `docs/copy_assignment_definition.md`, `copy_assign --families` |
| **missing copies** | expressed copies the reference does not contain (diverged or with rearranged exons): detect, characterise, screen, and hand copy number to DNA | `docs/archive/2026-11/O3_STATUS.md`, `docs/archive/2026-09/PREREG_o3_rna_only_2026-09-23.md`, `missing_copy_flag` |

In the thesis record these are objectives O1, O2 and O3; the ledger, register and pre-registrations use those
labels. `docs/THESIS_OBJECTIVES.md` states them in full.

## Build, test, run

```bash
cargo build --release          # binaries in target/release/, where the driver looks (real numbers always come from --release);
                               # after `pixi run build` they are in ~/rustle_target/release: pass --bin ~/rustle_target/release to the driver
cargo test --release           # ~1,600 tests; fixtures are in tests/fixtures and src/testdata
cargo test --profile dev-opt   # the iteration loop: full suite in minutes (opt 2, no LTO; see AGENTS.md §4)
```

Run the whole pipeline on one sample (one run directory, every stage resumable; `--to`/`--from STAGE` stop or start part way):

```bash
tools/minimap2_gate/build.sh               # once: minimap2-gate (minimap2 2.30 + exon gate) into ~/.local/bin; without it families is ~5x slower
tools/rustle_pipeline.sh all --bam reads.bam --fasta genome.fa --out run/s1 --threads 72
```

Needs `minimap2` (2.30) on PATH (or `RUSTLE_MINIMAP2`); the driver checks before the first stage. The BAM must be indexed,
minimap2 `-ax splice:hq -uf --eqx -Y -N 50 -p 0.1 --secondary=yes` (`docs/DATA.md`), and the FASTA must carry a `.fai`. Two
sidecars are shared by every run when their directory takes files: `BAM.as_table.tsv` (else `--as-table PATH`) and
`FASTA.splice.I8G.mmi` (not for a FASTA under $HOME: there each run builds `PREFIX.splice.mmi` unless `RUSTLE_INDEX_DIR` names a shared directory). A119b (96 GB BAM, CHM13) runs in 22 min from scratch on 72 cores; `docs/PERFORMANCE.md`.

`REPRODUCE.md` is the recipe book: every shipped number, the command that produces it, and the expected
output. `docs/DATA.md` lists the BAMs, genomes and annotations those recipes read (large inputs live outside
the repository). `docs/MODULE_STATUS.md` says which modules are reachable at defaults, which are opt-in, and
which only a side binary uses (a test checks that the nine modules of `src/family.rs` are listed; the status column is kept by hand).

Main binaries (`src/`; five more are not pipeline stages: `copy_phasing`, `family_gfa`, `semisec`, `semifam`, `slimbam`): the pipeline stages — `copy_assign` (`--assemble-only --genome-wide`
streams a whole BAM into loci and isoforms; `--families` is the copy assignment; `--flag-missing-copies` the
pairwise missing-copy test), `mcl_families` (family definition; `--from-gtf` runs the de novo stage from an
assembled GTF in one command, with the copy table the assignment consumes), `gw_family_catalog` (the LEGACY copy
catalog, assigned only with the driver's `--legacy-catalog`), `candidate_copies` (the opt-in candidates stage: candidate
copies absent from the reference, from each family's reads), `missing_copy_flag` (missing copies from RNA), `as_table` (one pass over a BAM → each molecule's genome-wide best alignment score; the pipeline's assemble stage seeds loci with the tied secondaries it names, `--no-seed-secondaries` to opt out) — the `utilities` binary (11 subcommands: the converters `mcl-port`, `locus-bed`, `gff-to-gtf`, `parcn`, and `guide`, `sd-scan`, `candidate-augment`, `sample-counts`, `end-artifacts`, `window-depth`, `psv-discover`), plus the standalone comparator `family_score` and filter `readthrough_filter`.
`tools/rustle_pipeline.sh assemble|families|candidates|catalog|assign|merged|flag|all` runs any stage, or all of them (`all` =
assemble, families, assign, flag); `merged` is the genome-wide entry point (assemble + families + assign on the families copy table, byte-identical, resumable — bench/MERGED_PIPELINE.md), with the shipped defaults (since 2026-09-29 these include bridge-aware regrouping,
`copy_assign --bridge-regroup f1v2`, and the containment escape `mcl_families --min-cov-shorter 0.70`; since 2026-10-02
`assign` reads the families' copy table, the legacy catalog only with `--legacy-catalog`, and `candidates` — O3 candidate
copies from each family's reads, `candidate_copies` — is OPT-IN, `--candidates`, because its pre-registered acceptance failed:
`docs/archive/2026-10/O3_CANDIDATES_ACCEPTANCE_2026-10-02.md`; the re-run acceptance passed on 2026-10-03, but its no-deletion control failed
and the default flip was reverted the same day, `docs/archive/2026-10/O3_CANDIDATES_CONTROL_A14_2026-10-03.md`; on a full BAM one batch of 50
families did not finish in a 10-minute call, spec §9b, R23; the driver runs every family in one call, only the binary's `--families` batches, so its peak memory exceeds that batch's 10.7 GB). `RUSTLE_BRIDGE_REGROUP=off RUSTLE_MIN_COV_SHORTER=0 ... --legacy-catalog` is
the earlier pipeline, byte for byte except `assign` on a catalog with cross-chromosome families (2026-10-02, `copy_assign
--help`, `--families`) — `REPRODUCE.md` §5 and its driver section, `bench/ASSEMBLY_POLISH.md` addendum 3). `bench/` holds the
Python files behind the record's analysis commands, among them the library `lib.py`, the subcommand scripts `score.py`,
`sim.py` and `truth.py`, `guided_pipeline.py`, `mcl_port.py`, and the `layer_order/` and `soto/` chains
(`bench/README.md` lists each one, and maps every retired script name to its current command); `tools/` the pipeline
driver and the cleanup-audit script.

### Stages added in November 2026

`all` now also writes two side files after its last stage, each a copy of `PREFIX.gtf` with extra models and nothing
removed: `PREFIX.earlyends.gtf` (`RUSTLE_EARLYENDS=0` turns it off) and `PREFIX.capstarts.gtf` (`RUSTLE_CAPSTARTS=0`;
refused, and reported, on a library without a usable cap G). `families` and `assign` use extra locus representatives
(`famreps`, `RUSTLE_FAMREPS=0` off; needs `python3`). Two sidecars next to the FASTA change the defaults when present:
`FASTA.sedef.bed` turns on SD damping in `families` and the `PREFIX.fam.sdcore.tsv` side table (`RUSTLE_SEDEF_BED=` off,
`RUSTLE_SDCORE=0` off), and `FASTA.exclude.bed` makes `assign` skip the families inside it (`RUSTLE_EXCLUDE_BED=` off).
Everything else below is opt-in, runs only when named, and leaves
the default products untouched:

| stage | product | register |
|---|---|---|
| `assemble --guide G.gtf` / `guided` | `PREFIX.guided.gtf` | 1622 |
| `siblings` | `PREFIX.sibling.gtf` | 1628 |
| `flairlike` | `PREFIX.flairlike.gtf` | 1640 |
| `altends` | `PREFIX.altends.gtf` | 1660 |
| `trunc5` | `PREFIX.trunc5.gtf` | 1677, 1692 |
| `shiftdrop` | `PREFIX.shiftdrop.gtf` | 1728-1739 |
| `chainscore` | `PREFIX.chainscore.gtf` (`PREFIX.gtf` without multi-exon chains below the fixed cut of the frozen chain score; default since 2026-11-16: the same chains are also removed from `PREFIX.final.gtf`, `RUSTLE_CHAINSCORE=0` = off) | 1677, 1760, 1761, 1769. Default weight file since 2026-11-18: `chainscore_r1769.json` (fitted on 11 samples of six species, cut 0.055); `RUSTLE_CHAINSCORE_MODEL=r1677` = the earlier one. |
| `shortexons` | `PREFIX.shortexons.gtf` (default reassembly after reads that skip an exon of <= 15 bp are rewritten through it; needs room for a second copy of the BAM while it runs) | 1746, 1753, 1754, 1757 |
| `oneread` | `PREFIX.oneread.gtf` (one-read chains with every junction in another read and a clean reading frame; `RUSTLE_ONEREAD=1` appends them to `PREFIX.final.gtf`; default test failed its yield bar) | 1766 |
| `assemble --names ANNOT` / `names` | `PREFIX.named.gtf`, `PREFIX.final.named.gtf` (the same models with the annotation's gene and transcript names as attributes; independent of `--guide`) | none |
| `eqclasses` | `PREFIX.eq.classes.tsv`, `PREFIX.eq.loci.tsv` (sets of loci the alignment score does not separate; per locus the reads placed only there, the maximum with every tied read, an EM count between the two) | none |
| `famwise` | `PREFIX.famwise.families.tsv`, `.tx.tsv`, `.classes.tsv`, `.reads.tsv`, `.gtf` (transcripts joined across loci by tied molecules, each molecule counted once; per transcript a count range and a status; reporting only) | none |
| `psi` | `PREFIX.psi.tsv` (PSI of each intron among the introns sharing its donor or acceptor, with the range tied reads allow); after `exonmap` also `PREFIX.fam.psi.tsv` (the same junction across the copies of a family, with added rows for copies whose junction the stage's read floors dropped or that never use it, and the columns `silent`, `n_eligible`, `range_lo`, `range_hi`, `gap`) and `PREFIX.psi.introns.tsv` | none |
| `repeatloci` | `PREFIX.repeat_loci.tsv` | 1698, 1702 |
| `exonmap` | `PREFIX.fam.exonmap.columns.tsv` + one GFA per family | 1595 |
| `relatives` | `PREFIX.fam.relatives.tsv` | 1676, 1737 |

Since 2026-11-16 the transcript lines of `PREFIX.final.gtf` also carry `orf_len` and `ptc` (`ptc "1"` = the longest ORF
stops more than 50 nt before the last junction; a tag, nothing is dropped; `RUSTLE_PTC=0` off; register 1766).

The header of `tools/rustle_pipeline.sh` has each stage's inputs and switches; `docs/LEDGER_O1.md` (update 2026-11-15)
has the verdicts. Several of these stages run Python scripts under `bench/` and `tools/` and need `python3` with pysam.

## The record

- `docs/LEDGER_O1.md`, `docs/LEDGER_O2.md`, `docs/LEDGER_O3.md` — the running ledgers (the O1 notebook to 2026-10, §6xx, is `docs/archive/2026-11/o1_ledger.md`).
- `docs/NEGATIVE_RESULTS_REGISTER.md` — every refuted idea, one row each; check it before proposing anything.
- `docs/PREREG_*.md` while open, then `docs/archive/YYYY-MM/PREREG_*.md` — decision rules written before looking, with their outcomes appended.
- `docs/REFERENCE.md` — glossary and the working set (what is live and what was pruned); open items are in each ledger's `## Current state`.

## Pruned material

Earlier documents, run outputs and retired scripts are not deleted: tracked files are at the git tags
`notebook-2026-09-19`, `notebook-2026-09-20`, `retired-modules-2026-09-20`, `notebook-2026-09-23`,
`notebook-2026-09-23b`, `notebook-2026-09-23c` and `notebook-2026-09-24` (`git checkout <tag> -- <path>`), and
everything moved out of the working tree sits in `~/Desktop/Rustle_attic/` with a manifest per wave.
