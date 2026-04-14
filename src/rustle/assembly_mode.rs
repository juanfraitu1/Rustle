//! Assembly mode: long-read, short-read, and mixed (C++ reference / reference script).
//!
//! ## Long-read assembly (-L)
//! - **C++ ref**: `longreads=true`; **reference assembler**: `-L`, `long_read_min_len=0`.
//! - All reads contribute to transfrag **abundance**; longstart/longend from read extent.
//! - Graph: when implemented, **longtrim** (split at read boundaries lstart/lend); else coverage trim optional.
//! - No srabund; no pair-aware localdist.
//! - BUNDLEDIST_LONGREAD = 0 (C++ reference reference script); localdist = 0 for long (C++ reference).
//!
//! ## Short-read assembly (no -L)
//! - **C++ ref**: `longreads=false`; **reference assembler**: omit `-L`.
//! - Graph: **trimnode_all** (coverage-based node splitting) only; no longtrim (C++ reference, 4071, 4142).
//! - Reads contribute to **srabund**; pair-aware; localdist = bundledist + longintronanchor (C++ reference).
//! - reference assembler.cpp: bundledist=50 when !longreads; mergeMode uses 250.
//!
//! ## Mixed mode (-L --long-read-min-len BP)
//! - **C++ ref**: `mixedMode=true` (--mix); **reference assembler**: `-L` and `--long-read-min-len BP`.
//! - Reads with query_length >= BP → long (abundance, longread=true); others → short (srabund).
//! - Graph: longtrim for long-read boundaries when present; trimnode_all for short (C++ reference uses longtrim when longreads && lstart/lend).
//! - process_srfrag: redistribute srabund to compatible transfrags (C++ reference; reference assembler process_srfrag when mixed_mode).
//! - localdrop more lenient in mixed (C++ reference localdrop=DROP*DROP).

use crate::types::AssemblyMode;

/// Bundle merge distance for long-read: 0 (BUNDLEDIST_LONGREAD).
pub const BUNDLE_MERGE_DIST_LONGREAD: u64 = 0;

/// Default bundle merge distance for short-read (reference assembler bundledist when !longreads).
pub const BUNDLE_MERGE_DIST_SHORTREAD_DEFAULT: u64 = 50;

/// Long intron anchor (C++ header longintronanchor); used for short-read localdist = bundledist + this.
pub const LONGINTRONANCHOR: u64 = 25;

/// Whether to apply coverage-based node splitting (trimnode_all) in this mode.
/// trimnode_all only when !longreads && !mergeMode; longtrim when longreads && lstart/lend.
pub fn use_coverage_trim(mode: AssemblyMode) -> bool {
    match mode {
        AssemblyMode::ShortRead => true,
        // C++ reference long-read uses longtrim/read boundaries, not short-read trimnode_all behavior.
        AssemblyMode::LongRead | AssemblyMode::Mixed => false,
    }
}

/// Whether to apply longtrim (read-boundary splits) when implemented.
/// longtrim when longreads && (lstart.Count() || lend.Count()).
pub fn use_longtrim(mode: AssemblyMode) -> bool {
    matches!(mode, AssemblyMode::LongRead | AssemblyMode::Mixed)
}
