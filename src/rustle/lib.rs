//! Rustle — multi-copy gene-family analysis over long-read RNA (great-ape pan-transcriptomics).
//!
//! The thesis lives entirely in `vg_family`: O1 family definition, O2 copy assignment under MAPQ-0
//! ambiguity, O3 reference-absent / missing copies (detect and flag) — over `genome` plus two small
//! foundational modules (`types` for the deterministic hash aliases, `bam` for BAM opening and CIGAR
//! exon blocks). Every binary in `src/bin/` imports only these.
//!
//! The legacy StringTie assembler + network-flow stack (~50k lines, ~55 modules) was RETIRED on
//! 2026-07-14; see `docs/RETIREMENT_AND_MIGRATION.md`. Its residue was removed in three steps:
//! 2026-09-20 (StringTie-exact presets, the VG-HMM rescue cluster, the dropped ASJ objective's modules;
//! tag `retired-modules-2026-09-20`) and 2026-09-24 (the bundle/read types, `RunConfig`, the dead half
//! of `bam`, `util::constants`, and the `minimizers`/`bridge_detector`/`repeat_catalog` modules; tag
//! `notebook-2026-09-24`). Build new work in `vg_family`.

pub mod types; // FixedBuild / DetHashMap / DetHashSet (deterministic FxHash containers)
pub mod bam; // open_bam + CIGAR -> exon blocks
pub mod genome; // genome / chromosome metadata (GenomeIndex, IndexedFasta)
pub mod vg_family; // THESIS LAYER: O1 family definition, O2 copy assignment, O3 missing copies
