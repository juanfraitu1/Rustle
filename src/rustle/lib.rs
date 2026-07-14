//! Rustle — multi-copy gene-family analysis over long-read RNA (great-ape pan-transcriptomics).
//!
//! The thesis lives entirely in `vg_family`: O1 family definition, O2 copy assignment under MAPQ-0
//! ambiguity, O3 allele-specific junctions, O4 reference-absent copies — over `genome` plus the
//! foundational IO/type modules (`util`, `types`, `bam`). The four thesis binaries
//! (`copy_assign`, `gw_family_catalog`, `asj`, `asj_verify`) import only these.
//!
//! The legacy StringTie assembler + network-flow stack (~50k lines, ~55 modules) was RETIRED on
//! 2026-07-14; see `docs/RETIREMENT_AND_MIGRATION.md`. Build new work in `vg_family`.

pub mod util; // bitset, bitvec, constants, coord, hard_counters (crate-wide low-level utilities)
pub mod types; // RunConfig, Bundle, Junction, all shared data types
pub mod bam; // BAM parsing / read ingestion
pub mod genome; // genome / chromosome metadata (GenomeIndex)
pub mod vg_family; // THESIS LAYER: O1 family detect/split, O2 copy_assign, O3 ASJ, O4 absent-copy

pub use types::{
    cjunctions_to_junction_stats, junction_stats_to_cjunctions, AssemblyMode, Bundle, BundleData,
    CBundlenode, CExon, CJunction, CMaxIntv, CPred, CPrediction, GArray, GEdge, GPVec, GVec,
    Junction, JunctionStats, RunConfig,
};
