# Rustle

Rust tools for **multi-copy gene families in long-read RNA (IsoSeq/HiFi)**: define the families from the
reads and the genome, assign ambiguous reads to copies, and flag copies the reference does not contain.
Substrate: gorilla (mGorGor1) and human (CHM13) IsoSeq. This is a thesis codebase; the scientific record is
kept next to the code and is pre-registered.

## What it does

| stage | what | where |
|---|---|---|
| **family definition** | multi-copy gene families from the reads and the genome — de novo (reads only) and guided (annotation), with the goal of closing the gap between them | `docs/seeded_family_definition.md`, `mcl_families`, `gw_family_catalog` |
| **copy assignment** | each ambiguous read to one copy, or abstention: PSVs, junctions and divergence, never a 1/k split | `docs/copy_assignment_definition.md`, `copy_assign --families` |
| **missing copies** | expressed copies the reference does not contain (diverged or with rearranged exons): detect, characterise, screen, and hand copy number to DNA | `docs/O3_STATUS.md`, `docs/PREREG_o3_rna_only_2026-09-23.md`, `missing_copy_flag` |

In the thesis record these are objectives O1, O2 and O3; the ledger, register and pre-registrations use those
labels. `docs/THESIS_OBJECTIVES.md` states them in full.

## Build, test, run

```bash
cargo build --release          # binaries in target/release/
cargo test --release           # ~860 tests; fixtures are in tests/fixtures and src/rustle/vg_family/testdata
```

`REPRODUCE.md` is the recipe book: every shipped number, the command that produces it, and the expected
output. `docs/DATA.md` lists the BAMs, genomes and annotations those recipes read (large inputs live outside
the repository). `docs/MODULE_STATUS.md` says which modules are reachable at defaults, which are opt-in, and
which only a side binary uses — enforced by a test.

Main binaries (`src/bin/`, ten): the pipeline stages — `copy_assign` (`--assemble-only --genome-wide`
streams a whole BAM into loci and isoforms; `--families` is the copy assignment; `--flag-missing-copies` the
pairwise missing-copy test), `mcl_families` (family definition; `--from-gtf` runs the de novo stage from an
assembled GTF in one command), `gw_family_catalog` (the copy catalog the assignment consumes),
`missing_copy_flag` (missing copies from RNA) — and five comparators and converters: `family_score`, `mcl_port`, `readthrough_filter`, `locus_bed`, `gff_to_gtf`, `parcn`.
`tools/rustle_pipeline.sh assemble|families|catalog|assign|flag|all` runs any stage, or all of them, with the shipped defaults. `bench/` holds the
39 analysis scripts the record cites (`bench/README.md` lists each one); `tools/` the sweep, audit and attic
scripts.

## The record

- `docs/o1_ledger.md` — the running notebook (sections §6xx are what the other documents cite).
- `docs/NEGATIVE_RESULTS_REGISTER.md` — every refuted idea, one row each; check it before proposing anything.
- `docs/PREREG_*.md` — decision rules written before looking, with their outcomes appended.
- `docs/PENDING_*.md` — open items; `docs/ACTIVE_WORKING_SET.md` — what is live and what was pruned.

## Pruned material

Earlier documents, run outputs and retired scripts are not deleted: tracked files are at the git tags
`notebook-2026-09-19`, `notebook-2026-09-20`, `retired-modules-2026-09-20` and `notebook-2026-09-23`
(`git checkout <tag> -- <path>`), and everything moved out of the working tree sits in
`~/Desktop/Rustle_attic/` with a manifest per wave.
