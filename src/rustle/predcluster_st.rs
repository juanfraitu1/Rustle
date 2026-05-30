//! StringTie-faithful prediction SELECTION (sub-project 1). Default OFF
//! (RUSTLE_PREDCLUSTER_ST=1). Runs ST's selection sub-stages in ST's order on the
//! candidate predictions Rustle's flow already produced. See
//! docs/superpowers/specs/2026-05-30-predcluster-selection-parity-design.md.
use crate::path_extract::Transcript;

/// Pass-through skeleton: returns all candidates kept (no selection yet).
/// Subsequent tasks add ST-faithful containment, isofrac, and collapse parity.
///
/// NOTE: the predcluster input/output transcript type is
/// `crate::path_extract::Transcript` (NOT `crate::types::Transcript`, which does
/// not exist — `types.rs` has no `Transcript` struct). Future containment/isofrac
/// ports should target `path_extract::Transcript`.
pub fn select_predictions_st(candidates: Vec<Transcript>) -> Vec<Transcript> {
    candidates
}
