//! Parallel "StringTie-exact" mode scaffold.
//!
//! When `RUSTLE_STRINGTIE_EXACT=1` is set, Rustle disables every known
//! Rustle-specific relaxation at once and enables every identified
//! StringTie-parity gate. The goal is to see where Rustle's accumulated
//! divergence from StringTie's algorithm actually lives, instead of
//! chasing divergences one at a time.
//!
//! This is a BUILDING MODE: as new divergences are discovered via trace
//! comparison (e.g. `RUSTLE_PATHPAT_OR_TRACE=1 RUSTLE_TRACE_LOCUS=...`),
//! they should be gated on `stringtie_exact()` here and wired into the
//! relevant site. Over time the mode accumulates StringTie-parity until
//! it converges.
//!
//! # Usage
//!
//! ```bash
//! RUSTLE_STRINGTIE_EXACT=1 ./target/release/rustle -L input.bam -o out.gtf
//! ```
//!
//! # Currently gated divergences
//!
//! | Site | Default behavior | StringTie-exact behavior |
//! |------|------------------|--------------------------|
//! | `path_extract::fwd_to_sink_fast_long` `past_seed` gate | Skip non-sink children when `i_coord > maxpath_coord` (Rustle-specific) | Gate disabled; StringTie has no equivalent |
//! | `path_extract::onpath_long` reach relaxation (line ~3142) | Allow tf past maxp when `node_can_reach` returns true | Strict: `trnode.last() > maxp` → false |
//! | `transfrag_process::process_transfrags` `novel_splice_rescue` | Promote contained long-read tfs with novel junctions to keeptrf | Disabled; StringTie strictly sets weak=1 |
//! | `path_extract::extract_transcripts` strict_redist filter | Already off-by-default (2026-04-21 flip) | Same — loose best_trf_match like StringTie |
//!
//! # Adding a new divergence gate
//!
//! 1. Identify a Rustle-specific deviation from StringTie's rlink.cpp behavior.
//! 2. Add an opt-in env flag for the specific site (e.g. `RUSTLE_FOO_OFF=1`).
//! 3. In the site, check `stringtie_exact() || specific_flag_set()`.
//! 4. Add an entry to the table above.
//! 5. Re-run with `RUSTLE_STRINGTIE_EXACT=1` to see the compound effect.

/// Returns true when the caller has opted into the full StringTie-exact mode.
/// Individual sites can still have their own finer-grained opt-out flags —
/// this is the meta-flag that enables everything in sync.
#[inline]
pub fn stringtie_exact() -> bool {
    std::env::var_os("RUSTLE_STRINGTIE_EXACT").is_some()
}

/// Returns true if StringTie-exact mode is active OR the site's specific
/// opt-out flag is set. Use at each divergence site to check whether the
/// Rustle-specific relaxation should be disabled.
#[inline]
pub fn parity_requested(specific_env_var: &str) -> bool {
    stringtie_exact() || std::env::var_os(specific_env_var).is_some()
}

/// One-time startup banner. Call from the top of the main pipeline entry
/// when stringtie_exact is active, so runs are clearly labeled.
pub fn maybe_emit_banner() {
    if stringtie_exact() {
        eprintln!(
            "[STRINGTIE_EXACT] RUSTLE_STRINGTIE_EXACT=1 — disabling Rustle-specific relaxations. \
             Gated sites: fwd_to_sink past_seed gate, onpath_long reach relaxation, \
             novel_splice_rescue. See src/rustle/stringtie_parity.rs for details."
        );
    }
}
