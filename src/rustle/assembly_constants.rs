//! C++ header / C++ reference constants — single source of truth for exon extension, trimming, and tolerance.
//!
//! ## Constants (C++ header)
//!
//! | Constant           | Value | Meaning |
//! |--------------------|-------|---------|
//! | `longintronanchor` | 25    | Window around splice sites; used in three contexts below. |
//! | `DROP`             | 0.5   | Coverage drop ratio threshold. |
//! | `ERROR_PERC`       | 0.1   | 10% noise tolerance. |
//! | `CHI_WIN`          | 100   | Sliding window size (bp). |
//! | `CHI_THR`          | 50    | Threshold window / min trim spacing (bp). |
//!
//! ## longintronanchor = 25 — three contexts (C++ reference)
//!
//! 1. **find_trims() / find_all_trims()** (lines 1872, 2135): the coverage scan loop skips the first
//!    and last 25 bp of every exon: `for i in start+longintronanchor .. end-longintronanchor`.
//!    Rationale: reads spanning junctions have degraded alignment near splice sites; scanning too
//!    close produces false trim signals.
//!
//! 2. **longtrim()** (line 2603): a detected read boundary is only accepted as a node split if it's
//!    >25 bp away from the node's start/end (same-strand guard).
//!
//! 3. **included_pred()** (line 17548): when checking if a small transcript's exon extends into
//!    the intron of a big transcript, the bpcov drop window at the exon boundary is ±25 bp
//!    (exondrop / introndrop).
//!
//! ## Coverage-based trimming — find_all_trims() localdrop by region
//!
//! | Region                  | localdrop              | Meaning |
//! |--------------------------|------------------------|---------|
//! | Very short (len < CHI_WIN) | ERROR_PERC/(10×DROP) = 0.02 | Very strict |
//! | Normal main region      | ERROR_PERC/DROP = 0.2  | Accept if left/right ratio < 0.2 |
//! | Mixed mode (main)        | DROP×DROP = 0.25       | More tolerant — 25% ratio |
//! | End region              | ERROR_PERC = 0.1       | Strict 10% at exon ends |
//!
//! ## included_pred() intronextension
//!
//! - Short-read: `2 * longintronanchor` = 50 bp.
//! - Long-read: `CHI_WIN` = 100 bp (more aggressive inclusion rules).

/// C++ header: longintronanchor = 25 — "window around erroneous splice sites".
/// Used in: (1) find_trims/find_all_trims loop guard, (2) longtrim boundary guard, (3) included_pred bpcov ±25 bp.
pub const LONGINTRONANCHOR: u64 = 25;

/// C++ header: DROP = 0.5 — coverage drop ratio threshold.
pub const DROP: f64 = 0.5;

/// C++ header: ERROR_PERC = 0.1 — 10% noise tolerance.
pub const ERROR_PERC: f64 = 0.1;

/// C++ header: CHI_WIN = 100 — sliding window size (bp).
pub const CHI_WIN: u64 = 100;

/// C++ header: CHI_THR = 50 — threshold window / min trim spacing (bp).
pub const CHI_THR: u64 = 50;

/// trthr — minimum abundance for synthetic trim transfrags.
pub const TRTHR: f64 = 1.0;

// --- Derived (for documentation / consistency) ---

/// Very short exon localdrop: ERROR_PERC/(10*DROP) = 0.02.
pub const LOCALDROP_VERY_SHORT: f64 = ERROR_PERC / (10.0 * DROP);

/// Normal main region localdrop: ERROR_PERC/DROP = 0.2.
pub const LOCALDROP_NORMAL: f64 = ERROR_PERC / DROP;

/// Mixed mode main region localdrop: DROP*DROP = 0.25.
pub const LOCALDROP_MIXED: f64 = DROP * DROP;

/// End region localdrop: ERROR_PERC = 0.1.
pub const LOCALDROP_END: f64 = ERROR_PERC;

/// included_pred intronextension for short-read: 2 * longintronanchor = 50 bp.
pub const INTRON_EXTENSION_SHORT: u64 = 2 * LONGINTRONANCHOR;

/// included_pred intronextension for long-read: CHI_WIN = 100 bp.
pub const INTRON_EXTENSION_LONG: u64 = CHI_WIN;
