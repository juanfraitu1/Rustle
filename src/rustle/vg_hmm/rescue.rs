//! Rescue pipeline: read prefilter → HMM scoring → clustering → synthetic bundle emission.

use crate::types::{Bundle, RunConfig};
use crate::vg::{FamilyGroup, NovelCandidate};
use anyhow::Result;

/// Top-level entry point. Phase 5 will replace this stub with the real
/// HMM-based rescue. Returning an empty Vec from the stub means
/// `--vg-discover-novel-mode hmm` is functionally equivalent to disabling
/// novel-copy discovery until Phase 5 lands.
pub fn run_rescue(
    _bam_path: &std::path::Path,
    _families: &[FamilyGroup],
    _bundles: &[Bundle],
    _config: &RunConfig,
) -> Result<Vec<NovelCandidate>> {
    Ok(Vec::new())
}
