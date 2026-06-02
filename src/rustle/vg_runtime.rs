//! Process-global "VG mode active" flag.
//!
//! Mirrors `debug_stage::STAGE_IO_ENABLED`: a few hot, config-less functions
//! (`map_reads::collect_read_nodes_exact` and the three `map_reads_to_graph*`
//! entrypoints) need to know whether `--vg` is active so they can drive the
//! parallel `anchored_coverage` channel WITHOUT touching the byte-identical
//! default de-novo path. Set once at pipeline start from `config.vg_mode`.

use std::sync::atomic::{AtomicBool, Ordering};

static VG_MODE_ENABLED: AtomicBool = AtomicBool::new(false);

/// Set at pipeline start from `config.vg_mode`. Idempotent.
#[inline]
pub fn set_vg_mode(enabled: bool) {
    VG_MODE_ENABLED.store(enabled, Ordering::Relaxed);
}

/// Fast path: `false` outside `--vg`, so the anchored-coverage channel is inert
/// and default de-novo output stays byte-identical.
#[inline]
pub fn vg_mode() -> bool {
    VG_MODE_ENABLED.load(Ordering::Relaxed)
}
