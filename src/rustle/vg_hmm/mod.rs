//! Variation-graph HMM rescue for novel gene-family copies.
//!
//! Replaces the flat k-mer rescue in `crate::vg::discover_novel_copies`
//! with a per-exon profile-HMM family model. See
//! `docs/superpowers/specs/2026-04-28-vg-novel-copy-hmm-design.md` for the
//! design, and `docs/UNMAPPED_FAMILY_RESCUE.md` for why the aligner misses
//! the reads this module rescues.

pub mod family_graph;
pub mod profile;
pub mod scorer;
pub mod rescue;
pub mod diagnostic;

pub use family_graph::{ExonClass, FamilyGraph, JunctionEdge};
pub use profile::ProfileHmm;                                     // Task 2.1
// pub use rescue::RescueResult;                                    // Task 5.x
pub use diagnostic::{RescueClass, classify_internal, classify_external, cigar_has_long_indel};  // Task 6.1
