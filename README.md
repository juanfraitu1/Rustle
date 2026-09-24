# Rustle

Rust tools for **multi-copy gene families in long-read RNA (IsoSeq/HiFi)**: define the families from the
reads and the genome, assign ambiguous reads to copies, and flag copies the reference does not contain.
Substrate: gorilla (mGorGor1) and human (CHM13) IsoSeq. This is a thesis codebase; the scientific record is
kept next to the code and is pre-registered.

## The three objectives

| | objective | where it lives |
|---|---|---|
| **O1** | topological RNA multi-copy **family definition** — de novo (reads only) and guided (annotation), with the goal of closing the gap between them | `docs/seeded_family_definition.md`, `mcl_families`, `family_define`, `gw_family_catalog` |
| **O2** | **copy assignment** under MAPQ-0 ambiguity: PSVs + junctions + divergence, assign-or-abstain, never 1/k | `docs/copy_assignment_definition.md`, `copy_assign` |
| **O3** | **reference-absent copies**: detect and flag, hand copy number to DNA | `docs/O3_STATUS.md`, `docs/PREREG_o3_rna_only_2026-09-23.md`, `o3_rna_flag` |

`docs/THESIS_OBJECTIVES.md` states them in full.

## Build, test, run

```bash
cargo build --release          # binaries in target/release/
cargo test --release           # ~860 tests; fixtures are in tests/fixtures and src/rustle/vg_family/testdata
```

`REPRODUCE.md` is the recipe book: every shipped number, the command that produces it, and the expected
output. `docs/DATA.md` lists the BAMs, genomes and annotations those recipes read (large inputs live outside
the repository). `docs/MODULE_STATUS.md` says which modules are reachable at defaults, which are opt-in, and
which only a side binary uses — enforced by a test.

Main binaries (`src/bin/`): `copy_assign` (loci → isoform assembly → GTF; `--assemble-only` is the assembler,
`--genome-wide` streams a whole BAM; O2 assignment and the O3 flag pass live here too), `mcl_families`,
`family_define`, `mcl_refine`, `gw_family_catalog` (O1), `o3_rna_flag` (O3), and the helpers `gff_to_gtf`,
`locus_bed`, `readthrough_filter`, `family_score`, `mcl_port`, `parcn`, `gamma_refine`, `filter_bam_by_as`,
`index_bam`, `bam_header`. `bench/` holds the analysis scripts and per-topic reports that the record cites;
`tools/` the sweep and audit scripts.

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
