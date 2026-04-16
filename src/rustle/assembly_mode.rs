//! Assembly mode: long-read only.
//!
//! ## Long-read assembly (-L)
//! - All reads contribute to transfrag **abundance**; longstart/longend from read extent.
//! - Graph: **longtrim** (split at read boundaries lstart/lend); coverage trim optional.
//! - No srabund; no pair-aware localdist.
//! - BUNDLEDIST_LONGREAD = 0; localdist = 0 for long

use crate::types::AssemblyMode;

/// Bundle merge distance for long-read: 0 (BUNDLEDIST_LONGREAD).
pub const BUNDLE_MERGE_DIST_LONGREAD: u64 = 0;

/// Long intron anchor (header longintronanchor).
pub const LONGINTRONANCHOR: u64 = 25;

/// Whether to apply coverage-based node splitting (trimnode_all) in this mode.
/// Long-read uses longtrim/read boundaries, not short-read trimnode_all behavior.
pub fn use_coverage_trim(_mode: AssemblyMode) -> bool {
    false
}

/// Whether to apply longtrim (read-boundary splits).
/// longtrim when longreads && (lstart.Count() || lend.Count()).
pub fn use_longtrim(_mode: AssemblyMode) -> bool {
    true
}
