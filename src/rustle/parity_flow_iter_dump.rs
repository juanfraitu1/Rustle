//! Per-iteration max-flow augmenting-path dump for cross-tool parity / Python visualizers.
//!
//! Enable with `RUSTLE_FLOW_ITER_TSV=/path/flow_iter.tsv`. Records the state of every
//! Edmonds-Karp augmenting-path iteration during max-flow assembly. Complements
//! `RUSTLE_PARITY_GRAPH_EDGES_TSV` (which dumps static graph topology); downstream
//! Python visualizers join these two TSVs to render per-iteration flow decompositions.
//!
//! Schema (header line):
//! `source	chrom	bdstart	bdend	strand	bundle_run	iter_idx	path_len	path_nodes	bottleneck	total_flow_after`
//!
//! Field semantics:
//! - `source` = literal `"rustle"`.
//! - `chrom`, `bdstart`, `bdend`, `strand` = bundle metadata (parsed from `bundle_id`
//!   shaped like `chrom:start-end`, populated via `set_bundle_context` thread-local).
//!   When the bundle context is unset (e.g. paths inside `max_flow.rs` reached without
//!   going through `extract_transcripts`), placeholders `"unknown"`/`'.'`/`0`/`0` are
//!   emitted; downstream visualizers join on `path_nodes` against the graph dump anyway.
//! - `bundle_run` = process-wide counter; `next_run_id()` increments at the start of
//!   each Edmonds-Karp invocation so multiple max-flow runs against the same bundle
//!   can be told apart.
//! - `iter_idx` = 1-based augmenting-path counter within a single `bundle_run`.
//! - `path_len` = number of nodes in the augmenting path (including source/sink endpoints).
//! - `path_nodes` = comma-separated string of numeric node indices, source-to-sink order.
//!   Source/sink are emitted as their numeric indices (no `S`/`T` markers); the downstream
//!   visualizer is designed to handle either convention but this module picks numeric.
//! - `bottleneck` = the residual capacity along the path (the value the algorithm
//!   subtracts from each edge this iteration), formatted `{:.4}`.
//! - `total_flow_after` = cumulative flow pushed including this iteration, formatted `{:.4}`.
//!
//! Mirrors the inline emitter pattern used by `parity_graph_edges_dump.rs`
//! (verbose by design — do not DRY).

use std::io::Write;
use std::sync::{Mutex, OnceLock};

static WR: OnceLock<Mutex<Option<std::fs::File>>> = OnceLock::new();
static HDR: OnceLock<Mutex<bool>> = OnceLock::new();
static RUN: OnceLock<Mutex<u64>> = OnceLock::new();

thread_local! {
    /// Thread-local bundle context. Set by `extract_transcripts` (and analogous
    /// callers) to carry chrom/strand/bundle_id into max_flow's EK loops without
    /// threading these through every function signature. Mirrors the existing
    /// `CURRENT_SEED_TF` pattern in `max_flow.rs`.
    static BUNDLE_CTX: std::cell::RefCell<Option<(String, char, String)>> =
        const { std::cell::RefCell::new(None) };
}

/// Set the active bundle context for the current thread. Returns the previous
/// context so the caller can restore it. Private so callers must use
/// `BundleCtxGuard`, which guarantees panic-safe restore.
fn set_bundle_context(
    bundle_chrom: &str,
    bundle_strand: char,
    bundle_id: &str,
) -> Option<(String, char, String)> {
    BUNDLE_CTX.with(|cell| {
        cell.borrow_mut().replace((
            bundle_chrom.to_string(),
            bundle_strand,
            bundle_id.to_string(),
        ))
    })
}

/// Restore (or clear) the bundle context for the current thread.
fn clear_bundle_context(prev: Option<(String, char, String)>) {
    BUNDLE_CTX.with(|cell| {
        *cell.borrow_mut() = prev;
    });
}

/// RAII guard that sets the bundle context on construction and restores the
/// previous value on drop. Use at the top of any function that calls into
/// `max_flow.rs` so the EK loops can dump bundle-tagged rows.
pub struct BundleCtxGuard {
    prev: Option<(String, char, String)>,
}

impl BundleCtxGuard {
    pub fn new(bundle_chrom: &str, bundle_strand: char, bundle_id: &str) -> Self {
        let prev = set_bundle_context(bundle_chrom, bundle_strand, bundle_id);
        Self { prev }
    }
}

impl Drop for BundleCtxGuard {
    fn drop(&mut self) {
        clear_bundle_context(self.prev.take());
    }
}

/// Process-wide monotonically increasing run id. Increments on each call.
pub fn next_run_id() -> u64 {
    let m = RUN.get_or_init(|| Mutex::new(0));
    let Ok(mut g) = m.lock() else { return 0 };
    *g += 1;
    *g
}

/// Emit one row for a single Edmonds-Karp augmenting-path iteration.
///
/// Bundle metadata (chrom/strand/bundle_id) is read from the thread-local
/// `BUNDLE_CTX` (typically set via `BundleCtxGuard` at the top of
/// `extract_transcripts`). Silently skips when the env var is unset or the
/// file can't be opened.
pub fn emit(
    bundle_run: u64,
    iter_idx: usize,
    path: &[usize],
    bottleneck: f64,
    total_flow_after: f64,
) {
    let gp = match std::env::var("RUSTLE_FLOW_ITER_TSV") {
        Ok(p) if !p.is_empty() => p,
        _ => return,
    };
    let wr = WR.get_or_init(|| {
        Mutex::new(
            std::fs::OpenOptions::new()
                .create(true).append(true).open(&gp).ok(),
        )
    });
    let hdr = HDR.get_or_init(|| Mutex::new(false));

    // Pull bundle metadata from thread-local context, with placeholders
    // when the context is unset.
    let (bundle_chrom, bundle_strand, bundle_id) = BUNDLE_CTX.with(|cell| {
        cell.borrow()
            .clone()
            .unwrap_or_else(|| ("unknown".to_string(), '.', String::new()))
    });
    let (bs, be) = bundle_id.split_once(':')
        .and_then(|(_, rng)| rng.split_once('-'))
        .map(|(s, e)| (
            s.trim().parse::<u64>().unwrap_or(0),
            e.trim().parse::<u64>().unwrap_or(0),
        ))
        .unwrap_or((0, 0));

    let mut path_nodes = String::new();
    for (i, n) in path.iter().enumerate() {
        if i > 0 {
            path_nodes.push(',');
        }
        path_nodes.push_str(&n.to_string());
    }

    if let (Ok(mut f_opt), Ok(mut hdr_w)) = (wr.lock(), hdr.lock()) {
        if let Some(f) = f_opt.as_mut() {
            if !*hdr_w {
                let _ = writeln!(f, "source\tchrom\tbdstart\tbdend\tstrand\tbundle_run\titer_idx\tpath_len\tpath_nodes\tbottleneck\ttotal_flow_after");
                *hdr_w = true;
            }
            let _ = writeln!(
                f,
                "rustle\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{:.4}\t{:.4}",
                bundle_chrom, bs, be, bundle_strand,
                bundle_run, iter_idx, path.len(), path_nodes,
                bottleneck, total_flow_after,
            );
        }
    }
}
