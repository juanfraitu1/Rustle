//! Variation-graph family analysis for novel gene-family copies.
//!
//! Builds a per-exon family graph (`family_graph`) over paralog copies and
//! drives structural detectors (mosaic, segdup, hidden_copy, positional) plus
//! k-mer-based novel-copy rescue. See `docs/UNMAPPED_FAMILY_RESCUE.md` for why
//! the aligner misses the reads this module rescues.

pub mod family_graph;
pub mod rescue;
pub mod diagnostic;
pub mod positional;
pub mod tandem;
pub mod consensus; // cross-copy consensus error-correction (subtractive precision lever)
pub mod mosaic;
pub mod segdup;
pub mod hidden_copy;
pub mod phasing;

pub use family_graph::{ExonClass, FamilyGraph, JunctionEdge};
pub use diagnostic::{RescueClass, classify_internal, classify_external, cigar_has_long_indel};  // Task 6.1
