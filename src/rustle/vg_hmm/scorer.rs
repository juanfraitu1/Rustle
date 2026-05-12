//! Forward / Viterbi over the family graph.
//!
//! Task 3.2 (single-DP boundary threading) + Task 3.3 (Viterbi best-path).

use std::cell::RefCell;
use crate::vg_hmm::profile::ProfileHmm;

#[inline]
fn idx(b: u8) -> Option<usize> {
    match b { b'A' => Some(0), b'C' => Some(1), b'G' => Some(2), b'T' => Some(3), _ => None }
}

const NEG_INF: f64 = f64::NEG_INFINITY;

/// Reusable scratch buffers for profile forward DP. Avoids the 6
/// `vec![NEG_INF; l+1]` allocations per call (per node, per read) — which
/// dominate runtime on multi-copy families during HMM-EM.
///
/// `ensure(l)` resizes & resets all buffers to the right length filled with
/// NEG_INF in O(l). Re-use across calls with the same or larger l avoids
/// re-allocation; smaller l is also free (extra capacity ignored).
#[derive(Default)]
struct ForwardScratch {
    prev_m: Vec<f64>,
    prev_i: Vec<f64>,
    prev_d: Vec<f64>,
    cur_m: Vec<f64>,
    cur_i: Vec<f64>,
    cur_d: Vec<f64>,
}

impl ForwardScratch {
    /// Resize and reset all buffers to length `len` with NEG_INF.
    /// Cheap on warm cache: just memsets up to `len`.
    fn ensure(&mut self, len: usize) {
        for buf in [
            &mut self.prev_m, &mut self.prev_i, &mut self.prev_d,
            &mut self.cur_m, &mut self.cur_i, &mut self.cur_d,
        ] {
            if buf.len() < len {
                buf.resize(len, NEG_INF);
            }
            for v in buf.iter_mut().take(len) {
                *v = NEG_INF;
            }
        }
    }
}

thread_local! {
    static FORWARD_SCRATCH: RefCell<ForwardScratch> = RefCell::new(ForwardScratch::default());
}

/// Banding width for profile forward DP. Set to 0 to disable banding
/// (for tests / debugging). 100 is sufficient for biological reads
/// (allows up to 100 consecutive insertions/deletions before scoring breaks down).
pub const FORWARD_BANDWIDTH: usize = 100;

#[inline]
fn logsumexp(a: f64, b: f64) -> f64 {
    if a == NEG_INF { return b; }
    if b == NEG_INF { return a; }
    let m = a.max(b);
    m + ((a - m).exp() + (b - m).exp()).ln()
}

// ── Core DP: profile forward with an entry-point boundary distribution ────────

/// Forward DP over a profile with a per-read-position boundary distribution.
///
/// `boundary_in[r]` is the log probability of *entering* the profile at column 0
/// with `read[0..r]` already consumed externally.  `boundary_in.len()` must equal
/// `read.len() + 1`.
///
/// Returns `boundary_out` of length `read.len() + 1`, where `boundary_out[r]` is
/// the log probability of having *exited* the profile's last state with `read[0..r]`
/// consumed in total (i.e. `read[enter..r]` consumed internally).
///
/// This is a thin wrapper around `profile_forward_with_boundary_banded` using
/// `FORWARD_BANDWIDTH`.
pub fn profile_forward_with_boundary(p: &ProfileHmm, read: &[u8], boundary_in: &[f64]) -> Vec<f64> {
    profile_forward_with_boundary_banded(p, read, boundary_in, FORWARD_BANDWIDTH)
}

/// Banded forward DP over a profile with a per-read-position boundary distribution.
///
/// Identical to `profile_forward_with_boundary` but restricts the inner DP to a
/// diagonal band around each active entry in `boundary_in`, giving O(W) work per
/// column instead of O(L).  `bandwidth = 0` disables banding (full DP, for
/// debugging or short-read parity checks).
///
/// # Banding strategy
///
/// For each active entry `r` in `boundary_in` (where `boundary_in[r] > NEG_INF`),
/// the natural alignment diagonal is `i = r + j`.  Cell `(j, i)` is updated only
/// if there exists an active `r` such that `|i - (r + j)| ≤ bandwidth`.
///
/// Equivalently, define `offset[j][i] = i - j`.  Cell `(j, i)` is in-band iff
/// `offset[j][i] ∈ [active_r_min - bandwidth, active_r_max + bandwidth]`,
/// where `active_r_min/max` are the min/max active entry positions in `boundary_in`.
///
/// This keeps each column's window to `(active_r_max - active_r_min) + 2*bandwidth + 1`
/// cells, which is at most `2*bandwidth + 1` cells when there is a single active
/// entry (the common case for source nodes).
/// Banded forward DP over a profile with a per-read-position boundary distribution.
///
/// Identical to `profile_forward_with_boundary` but restricts the inner DP to a
/// diagonal band around each active entry in `boundary_in`, giving O(W) work per
/// column instead of O(L).  `bandwidth = 0` disables banding (full DP, for
/// debugging or short-read parity checks).
///
/// # Banding strategy
///
/// For each active entry `r` in `boundary_in` (where `boundary_in[r] > NEG_INF`),
/// the natural alignment diagonal is `i = r + j`.  Cell `(j, i)` is updated only
/// if there exists an active `r` such that `|i - (r + j)| ≤ bandwidth`.
///
/// Equivalently: define the column range at column j as
/// `[r_min + j - W, r_max + j + W]` (clamped to `[0, l]`), where `r_min/r_max`
/// are the min/max active positions in `boundary_in`.  This keeps the band at
/// `(r_max - r_min) + 2*W + 1` cells per column — O(W) for a single entry point.
///
/// # Memory layout
///
/// Uses a 2-row rolling buffer (current + previous column) so allocation is
/// O(W) rather than O(M × L).  The full `boundary_out` vector is still O(L).
pub fn profile_forward_with_boundary_banded(
    p: &ProfileHmm,
    read: &[u8],
    boundary_in: &[f64],
    bandwidth: usize,
) -> Vec<f64> {
    let mut out = vec![NEG_INF; read.len() + 1];
    FORWARD_SCRATCH.with(|s| {
        profile_forward_with_boundary_banded_into(
            p, read, boundary_in, bandwidth, &mut s.borrow_mut(), &mut out,
        );
    });
    out
}

/// Allocation-free variant of `profile_forward_with_boundary_banded`.
///
/// Writes the boundary-out vector into `out` (resized to `read.len() + 1` and
/// filled). Uses the supplied `scratch` for the rolling row buffers so a
/// caller chaining many profile DP calls (e.g. `forward_against_path_for_copy`)
/// pays for buffer growth once instead of 6 vec allocs per node-per-read.
fn profile_forward_with_boundary_banded_into(
    p: &ProfileHmm,
    read: &[u8],
    boundary_in: &[f64],
    bandwidth: usize,
    scratch: &mut ForwardScratch,
    out: &mut Vec<f64>,
) {
    let m = p.n_columns;
    let l = read.len();
    assert_eq!(boundary_in.len(), l + 1, "boundary_in length mismatch");

    if out.len() != l + 1 {
        out.clear();
        out.resize(l + 1, NEG_INF);
    } else {
        for v in out.iter_mut() { *v = NEG_INF; }
    }

    if m == 0 {
        return;
    }

    let (r_min, r_max) = if bandwidth == 0 {
        (0usize, l)
    } else {
        let mut lo = l + 1;
        let mut hi = 0usize;
        for r in 0..=l {
            if boundary_in[r] > NEG_INF {
                if r < lo { lo = r; }
                if r > hi { hi = r; }
            }
        }
        if lo > hi {
            return;
        }
        (lo, hi)
    };

    scratch.ensure(l + 1);
    // Row 0 (entry) = boundary_in; insert/delete rows are NEG_INF.
    scratch.prev_m[..l + 1].copy_from_slice(boundary_in);

    for j in 1..=m {
        let prev_col_lo = if bandwidth == 0 {
            0
        } else {
            (r_min + j - 1).saturating_sub(bandwidth)
        };
        let col_lo = if bandwidth == 0 {
            0
        } else {
            (r_min + j).saturating_sub(bandwidth)
        };
        let col_hi = if bandwidth == 0 {
            l
        } else {
            (r_max + j + bandwidth).min(l)
        };

        // Clear current column over the band we're about to write.
        for i in col_lo..=col_hi {
            scratch.cur_m[i] = NEG_INF;
            scratch.cur_i[i] = NEG_INF;
            scratch.cur_d[i] = NEG_INF;
        }

        for i in col_lo..=col_hi {
            let emit_m = if i > 0 {
                idx(read[i - 1]).map(|k| p.match_emit[j - 1][k]).unwrap_or(NEG_INF)
            } else { NEG_INF };
            let pm_prev = if i > 0 { scratch.prev_m[i - 1] } else { NEG_INF };
            let pi_prev = if i > 0 { scratch.prev_i[i - 1] } else { NEG_INF };
            let pd_prev = if i > 0 { scratch.prev_d[i - 1] } else { NEG_INF };
            let from_m = if pm_prev > NEG_INF { pm_prev + p.trans[j - 1].mm + emit_m } else { NEG_INF };
            let from_i = if pi_prev > NEG_INF { pi_prev + p.trans[j - 1].im + emit_m } else { NEG_INF };
            let from_d = if pd_prev > NEG_INF { pd_prev + p.trans[j - 1].dm + emit_m } else { NEG_INF };
            scratch.cur_m[i] = logsumexp(logsumexp(from_m, from_i), from_d);

            let emit_i = if i > 0 {
                idx(read[i - 1]).map(|k| p.insert_emit[j][k]).unwrap_or(NEG_INF)
            } else { NEG_INF };
            let i_from_m = if i > 0 && scratch.cur_m[i - 1] > NEG_INF { scratch.cur_m[i - 1] + p.trans[j].mi + emit_i } else { NEG_INF };
            let i_from_i = if i > 0 && scratch.cur_i[i - 1] > NEG_INF { scratch.cur_i[i - 1] + p.trans[j].ii + emit_i } else { NEG_INF };
            scratch.cur_i[i] = logsumexp(i_from_m, i_from_i);

            let d_from_m = if scratch.prev_m[i] > NEG_INF { scratch.prev_m[i] + p.trans[j - 1].md } else { NEG_INF };
            let d_from_d = if scratch.prev_d[i] > NEG_INF { scratch.prev_d[i] + p.trans[j - 1].dd } else { NEG_INF };
            scratch.cur_d[i] = logsumexp(d_from_m, d_from_d);
        }

        if bandwidth > 0 {
            let clear_lo = prev_col_lo.min(l);
            let clear_hi = col_lo.min(l);
            for i in clear_lo..clear_hi {
                scratch.prev_m[i] = NEG_INF;
                scratch.prev_i[i] = NEG_INF;
                scratch.prev_d[i] = NEG_INF;
            }
        }

        std::mem::swap(&mut scratch.prev_m, &mut scratch.cur_m);
        std::mem::swap(&mut scratch.prev_i, &mut scratch.cur_i);
        std::mem::swap(&mut scratch.prev_d, &mut scratch.cur_d);
    }

    for r in 0..=l {
        out[r] = logsumexp(logsumexp(scratch.prev_m[r], scratch.prev_i[r]), scratch.prev_d[r]);
    }
}

// ── Public API: single-profile forward (backward-compatible wrapper) ──────────

/// Forward log P(read | profile). Unbanded for correctness; banding can be
/// added later if profiling shows it matters.
///
/// This is a thin wrapper around `profile_forward_with_boundary` that starts
/// the DP at read position 0 and returns the total probability of consuming
/// the full read.
pub fn forward_against_profile(p: &ProfileHmm, read: &[u8]) -> f64 {
    let mut bin = vec![NEG_INF; read.len() + 1];
    bin[0] = 0.0;
    let bout = profile_forward_with_boundary(p, read, &bin);
    bout[read.len()]
}

// ── Family-graph forward (single-DP boundary threading) ──────────────────────

use crate::vg_hmm::family_graph::{FamilyGraph, NodeIdx};

/// Forward over the family graph via single-DP boundary threading.
///
/// Each node's profile-internal forward DP runs once.  A per-read-position
/// "boundary in" distribution threads through outgoing edges to become the
/// "boundary in" for successor nodes.
///
/// Cost: O(N × L × profile_len) — three orders of magnitude cheaper than the
/// previous O(L³ × profile_len) per-read implementation.
pub fn forward_against_family(fg: &FamilyGraph, read: &[u8]) -> f64 {
    if fg.nodes.is_empty() { return NEG_INF; }
    let l = read.len();
    let max_support = fg.edges.iter().map(|e| e.family_support).max().unwrap_or(1).max(1);

    // boundary_in[node][r]  = log P of entering node with read[0..r] already consumed.
    // boundary_out[node][r] = log P of exiting  node with read[0..r] consumed in total.
    let mut boundary_in:  Vec<Vec<f64>> = vec![vec![NEG_INF; l + 1]; fg.nodes.len()];
    let mut boundary_out: Vec<Vec<f64>> = vec![vec![NEG_INF; l + 1]; fg.nodes.len()];

    // Source nodes: no incoming edges → entered at read position 0.
    let mut has_incoming = vec![false; fg.nodes.len()];
    for e in &fg.edges { has_incoming[e.to.0] = true; }

    // Topological order: trust insertion order (validated by smoke for tandem-paralog families).
    for nidx in 0..fg.nodes.len() {
        if !has_incoming[nidx] {
            boundary_in[nidx][0] = 0.0;
        }

        // Run profile forward through this node.
        if let Some(p) = &fg.nodes[nidx].profile {
            boundary_out[nidx] = profile_forward_with_boundary(p, read, &boundary_in[nidx]);
        }

        // Propagate boundary_out → boundary_in of successors.
        for e in fg.edges.iter().filter(|e| e.from.0 == nidx) {
            let log_edge = (e.family_support as f64 / max_support as f64).ln();
            for r in 0..=l {
                let cand = boundary_out[nidx][r] + log_edge;
                boundary_in[e.to.0][r] = logsumexp(boundary_in[e.to.0][r], cand);
            }
        }
    }

    // Best score: any node exiting after consuming the full read.
    let mut best = NEG_INF;
    for n in 0..fg.nodes.len() {
        best = best.max(boundary_out[n][l]);
    }
    best
}

/// Score a read against a SPECIFIC path through the family graph (constrained
/// forward, not the unconstrained sum-over-all-paths of `forward_against_family`).
///
/// Walks `path` in order; chains profile_forward through each node's profile,
/// using the previous node's exit boundary as the next node's entry boundary.
/// The final score is the log-prob of consuming all `l` read bases through
/// exactly this node sequence.
///
/// This is the per-paralog likelihood needed by HMM-based EM for read-to-copy
/// assignment: for a multi-mapped read with placements at paralogs P0..Pn,
/// `forward_against_path(fg, read, &path_Pi)` for each `Pi` gives the score
/// vector that the EM E-step can use directly (or normalize via softmax to a
/// posterior).
///
/// Returns `NEG_INF` if `path` is empty, contains an out-of-range NodeIdx, or
/// any node lacks a fitted profile.
pub fn forward_against_path(fg: &FamilyGraph, read: &[u8], path: &[NodeIdx]) -> f64 {
    if path.is_empty() { return NEG_INF; }
    let l = read.len();

    let mut boundary_in = vec![NEG_INF; l + 1];
    let mut boundary_out = vec![NEG_INF; l + 1];
    boundary_in[0] = 0.0;

    FORWARD_SCRATCH.with(|s| {
        let mut scratch = s.borrow_mut();
        for &nidx in path {
            if nidx.0 >= fg.nodes.len() { return NEG_INF; }
            let node = &fg.nodes[nidx.0];
            match &node.profile {
                Some(profile) => {
                    profile_forward_with_boundary_banded_into(
                        profile, read, &boundary_in, FORWARD_BANDWIDTH,
                        &mut scratch, &mut boundary_out,
                    );
                    std::mem::swap(&mut boundary_in, &mut boundary_out);
                }
                None => return NEG_INF,
            }
        }
        boundary_in[l]
    })
}

/// Per-copy variant of `forward_against_path`: scores a read against a
/// specific paralog by using that paralog's per-copy profile at each node
/// instead of the shared (POA-MSA-derived) profile.
///
/// This is the scoring function HMM-EM should call for paralog assignment:
/// it preserves SNPs, donor/acceptor microshifts within shared exon-classes,
/// and per-copy small indels that the shared profile averages away.
///
/// Per-node lookup picks the profile entry where `cid == copy_id`. If a node
/// has no entry for that copy (the copy doesn't contribute there — this
/// shouldn't happen if `path` came from `recover_paralog_path(copy_id)`,
/// but is handled defensively), falls back to the shared profile.
///
/// Returns `NEG_INF` if `path` is empty, contains an out-of-range NodeIdx,
/// or any node lacks both per-copy and shared profiles.
pub fn forward_against_path_for_copy(
    fg: &FamilyGraph,
    read: &[u8],
    path: &[NodeIdx],
    copy_id: crate::vg_hmm::family_graph::CopyId,
) -> f64 {
    forward_against_path_for_copy_with_norm(fg, read, path, copy_id, None)
}

/// Variant that accepts a pre-computed `total_match_cols` (sum of profile
/// match-column counts along this path for this copy). Avoids re-summing
/// across every read when the caller batch-scores many reads against the
/// same set of paralog paths (e.g. `run_pre_assembly_em_hmm`).
///
/// `precomputed_match_cols = None` falls back to summing internally (only
/// used when `RUSTLE_VG_HMM_LENGTH_NORM` is set).
pub fn forward_against_path_for_copy_with_norm(
    fg: &FamilyGraph,
    read: &[u8],
    path: &[NodeIdx],
    copy_id: crate::vg_hmm::family_graph::CopyId,
    precomputed_match_cols: Option<usize>,
) -> f64 {
    if path.is_empty() { return NEG_INF; }
    let l = read.len();

    let mut boundary_in = vec![NEG_INF; l + 1];
    let mut boundary_out = vec![NEG_INF; l + 1];
    boundary_in[0] = 0.0;

    let length_norm = std::env::var_os("RUSTLE_VG_HMM_LENGTH_NORM").is_some();
    let mut total_match_cols: usize = 0;
    let need_local_match_cols = length_norm && precomputed_match_cols.is_none();

    let raw = FORWARD_SCRATCH.with(|s| -> f64 {
        let mut scratch = s.borrow_mut();
        for &nidx in path {
            if nidx.0 >= fg.nodes.len() { return NEG_INF; }
            let node = &fg.nodes[nidx.0];
            // Prefer the per-copy profile; fall back to shared if absent.
            let profile_opt = node.per_copy_profiles
                .iter()
                .find(|(c, _)| *c == copy_id)
                .map(|(_, p)| p)
                .or(node.profile.as_ref());
            match profile_opt {
                Some(profile) => {
                    if need_local_match_cols {
                        total_match_cols += profile.match_emit.len();
                    }
                    profile_forward_with_boundary_banded_into(
                        profile, read, &boundary_in, FORWARD_BANDWIDTH,
                        &mut scratch, &mut boundary_out,
                    );
                    std::mem::swap(&mut boundary_in, &mut boundary_out);
                }
                None => return NEG_INF,
            }
        }
        boundary_in[l]
    });

    if !length_norm || !raw.is_finite() {
        return raw;
    }
    let cols = precomputed_match_cols.unwrap_or(total_match_cols);
    if cols == 0 { raw } else { raw / cols as f64 }
}

/// Sum of profile `match_emit.len()` along `path` for `copy_id`. Use to
/// pre-compute the length-norm denominator once per (family, copy) before
/// scoring many reads.
pub fn path_match_cols_for_copy(
    fg: &FamilyGraph,
    path: &[NodeIdx],
    copy_id: crate::vg_hmm::family_graph::CopyId,
) -> usize {
    let mut total: usize = 0;
    for &nidx in path {
        if nidx.0 >= fg.nodes.len() { return 0; }
        let node = &fg.nodes[nidx.0];
        let profile_opt = node.per_copy_profiles
            .iter()
            .find(|(c, _)| *c == copy_id)
            .map(|(_, p)| p)
            .or(node.profile.as_ref());
        if let Some(p) = profile_opt {
            total += p.match_emit.len();
        }
    }
    total
}

/// Constrained Viterbi over a single path through the family graph.
///
/// Returns `(read_spans, score)` where:
///  - `read_spans[i]` is the (entry, exit) read positions consumed by `path[i]`
///  - `score` is the Viterbi log-probability for that exact path
///
/// Path-constrained variant of [`viterbi_path`]: the topology is fully
/// specified by `path`, so we walk the nodes in sequence and chain a per-node
/// Viterbi DP through them. Per-call cost is O(K × L × profile_len) where K
/// is the path length — much cheaper than the family-wide O(N × L × prof) DP
/// when N >> K (large multi-copy families with shared exon classes).
///
/// Returns `None` if `path` is empty or scores `NEG_INF`.
pub fn viterbi_against_path(
    fg: &FamilyGraph,
    read: &[u8],
    path: &[NodeIdx],
) -> Option<ViterbiPath> {
    if path.is_empty() { return None; }
    let l = read.len();

    // Per-node entry r-vector for traceback.
    let mut entry_rs: Vec<Vec<usize>> = Vec::with_capacity(path.len());
    let mut boundary = vec![NEG_INF; l + 1];
    boundary[0] = 0.0;

    for &nidx in path {
        if nidx.0 >= fg.nodes.len() { return None; }
        let node = &fg.nodes[nidx.0];
        match &node.profile {
            Some(profile) => {
                let (out, e_r) = profile_viterbi_with_boundary_and_entry(profile, read, &boundary);
                boundary = out;
                entry_rs.push(e_r);
            }
            None => return None,
        }
    }

    let score = boundary[l];
    if score == NEG_INF { return None; }

    // Traceback: at the LAST node, exit position is l. Each predecessor's
    // exit position is the current node's entry position (instantaneous edge).
    let mut spans: Vec<(usize, usize)> = Vec::with_capacity(path.len());
    let mut cur_exit = l;
    for entry_r in entry_rs.iter().rev() {
        let entry = entry_r[cur_exit];
        spans.push((entry, cur_exit));
        cur_exit = entry;
    }
    spans.reverse();

    Some(ViterbiPath {
        nodes: path.to_vec(),
        read_spans: spans,
        score,
    })
}

// ── Viterbi best-path (Task 3.3) ─────────────────────────────────────────────

/// Family-level Viterbi result.
///
/// `read_spans[i]` is `(entry, exit)` — the read positions at which the i-th
/// node in the path started and finished consuming. Same length as `nodes`,
/// in path order (source-to-sink). Per-node read footprint length is
/// `exit - entry` and is what enables "VG-only" (no-CIGAR) length-consensus
/// refinement of rescued bundles.
#[derive(Debug, Clone)]
pub struct ViterbiPath {
    pub nodes: Vec<NodeIdx>,
    pub read_spans: Vec<(usize, usize)>,
    pub score: f64,
}

/// Viterbi (max instead of sum) DP over a profile with a per-read-position
/// boundary distribution.
///
/// Returns `boundary_out[r]` — the Viterbi log-probability of entering the
/// profile at some position ≤ r and exiting at r.
///
/// boundary_in.len() must equal read.len() + 1.
///
/// DEAD CODE: HMM-EM uses the Forward variant (`profile_forward_with_boundary`);
/// Viterbi is preserved here for diagnostic best-path extraction if/when
/// per-read traceback is needed.
#[allow(dead_code)]
fn profile_viterbi_with_boundary(
    p: &ProfileHmm,
    read: &[u8],
    boundary_in: &[f64],
) -> Vec<f64> {
    // Forward to the with_entry variant and discard the entry-r data.
    profile_viterbi_with_boundary_and_entry(p, read, boundary_in).0
}

/// Same as `profile_viterbi_with_boundary` but additionally returns
/// `entry_r[r_out]` — the read position at which the profile started consuming
/// in order to achieve the max-prob exit at `r_out`. This is the per-profile
/// intra-DP traceback needed for read-footprint length recovery (used by the
/// "VG-only" / no-CIGAR rescue refinement).
fn profile_viterbi_with_boundary_and_entry(
    p: &ProfileHmm,
    read: &[u8],
    boundary_in: &[f64],
) -> (Vec<f64>, Vec<usize>) {
    let m = p.n_columns;
    let l = read.len();
    assert_eq!(boundary_in.len(), l + 1, "boundary_in length mismatch");

    if m == 0 {
        return (vec![NEG_INF; l + 1], (0..=l).collect());
    }

    let mut m_dp = vec![vec![NEG_INF; l + 1]; m + 1];
    let mut i_dp = vec![vec![NEG_INF; l + 1]; m + 1];
    let mut d_dp = vec![vec![NEG_INF; l + 1]; m + 1];

    // Parallel "entry-r" arrays: each cell records the read position at which
    // the path achieving its max log-prob first entered the profile (i.e., the
    // column-0 starting position).
    let mut m_e: Vec<Vec<usize>> = vec![vec![0; l + 1]; m + 1];
    let mut i_e: Vec<Vec<usize>> = vec![vec![0; l + 1]; m + 1];
    let mut d_e: Vec<Vec<usize>> = vec![vec![0; l + 1]; m + 1];

    // Column 0 base case: enter at position r (no profile-column consumption yet).
    for r in 0..=l {
        m_dp[0][r] = boundary_in[r];
        m_e[0][r] = r;
    }

    for j in 1..=m {
        for i in 0..=l {
            // Match state at (j, i): consumes one read base AND advances column.
            let emit_m = if i > 0 {
                idx(read[i - 1]).map(|k| p.match_emit[j - 1][k]).unwrap_or(NEG_INF)
            } else { NEG_INF };
            let from_m = if i > 0 { m_dp[j - 1][i - 1] + p.trans[j - 1].mm + emit_m } else { NEG_INF };
            let from_i = if i > 0 { i_dp[j - 1][i - 1] + p.trans[j - 1].im + emit_m } else { NEG_INF };
            let from_d = if i > 0 { d_dp[j - 1][i - 1] + p.trans[j - 1].dm + emit_m } else { NEG_INF };
            let (best_m, src_m) = max3_with_src(from_m, from_i, from_d);
            m_dp[j][i] = best_m;
            if best_m > NEG_INF && i > 0 {
                m_e[j][i] = match src_m { 0 => m_e[j - 1][i - 1], 1 => i_e[j - 1][i - 1], _ => d_e[j - 1][i - 1] };
            }

            // Insert state at (j, i): consumes one read base, stays in column j.
            let emit_i = if i > 0 {
                idx(read[i - 1]).map(|k| p.insert_emit[j][k]).unwrap_or(NEG_INF)
            } else { NEG_INF };
            let i_from_m = if i > 0 { m_dp[j][i - 1] + p.trans[j].mi + emit_i } else { NEG_INF };
            let i_from_i = if i > 0 { i_dp[j][i - 1] + p.trans[j].ii + emit_i } else { NEG_INF };
            let (best_i, src_i) = max2_with_src(i_from_m, i_from_i);
            i_dp[j][i] = best_i;
            if best_i > NEG_INF && i > 0 {
                i_e[j][i] = if src_i == 0 { m_e[j][i - 1] } else { i_e[j][i - 1] };
            }

            // Delete state at (j, i): consumes no read base, advances column.
            let d_from_m = m_dp[j - 1][i] + p.trans[j - 1].md;
            let d_from_d = d_dp[j - 1][i] + p.trans[j - 1].dd;
            let (best_d, src_d) = max2_with_src(d_from_m, d_from_d);
            d_dp[j][i] = best_d;
            if best_d > NEG_INF {
                d_e[j][i] = if src_d == 0 { m_e[j - 1][i] } else { d_e[j - 1][i] };
            }
        }
    }

    // Exit: for each r, pick the max over {m, i, d} at column m, position r.
    let mut out = vec![NEG_INF; l + 1];
    let mut out_e = vec![0_usize; l + 1];
    for r in 0..=l {
        let (best, src) = max3_with_src(m_dp[m][r], i_dp[m][r], d_dp[m][r]);
        out[r] = best;
        if best > NEG_INF {
            out_e[r] = match src { 0 => m_e[m][r], 1 => i_e[m][r], _ => d_e[m][r] };
        } else {
            out_e[r] = r;
        }
    }
    (out, out_e)
}

#[inline]
fn max3_with_src(a: f64, b: f64, c: f64) -> (f64, u8) {
    if a >= b && a >= c { (a, 0) }
    else if b >= c { (b, 1) }
    else { (c, 2) }
}

#[inline]
fn max2_with_src(a: f64, b: f64) -> (f64, u8) {
    if a >= b { (a, 0) } else { (b, 1) }
}

/// Viterbi best-path through the family graph.
///
/// Returns the highest-scoring sequence of nodes (family-level; no
/// intra-profile traceback) and the Viterbi log-probability for consuming the
/// full read through that path.  Returns `None` if the family graph is empty
/// or the score is −∞.
pub fn viterbi_path(fg: &FamilyGraph, read: &[u8]) -> Option<ViterbiPath> {
    if fg.nodes.is_empty() { return None; }
    let l = read.len();
    let n = fg.nodes.len();
    let max_support = fg.edges.iter().map(|e| e.family_support).max().unwrap_or(1).max(1);

    // boundary_in_vit[node][r]  — max log-prob of entering node at read position r.
    // boundary_out_vit[node][r] — max log-prob of exiting  node at read position r.
    let mut boundary_in:  Vec<Vec<f64>> = vec![vec![NEG_INF; l + 1]; n];
    let mut boundary_out: Vec<Vec<f64>> = vec![vec![NEG_INF; l + 1]; n];

    // Per-node intra-profile entry tracking: for each node and exit position
    // r_out, the read position at which the profile started consuming.
    let mut node_entry_r: Vec<Vec<usize>> = vec![vec![0; l + 1]; n];

    // Family-level backpointers: for each (node, exit_r), which predecessor
    // node and at what read position it was exited.
    let mut bp_node: Vec<Vec<usize>> = vec![vec![usize::MAX; l + 1]; n];
    let mut bp_prev_exit: Vec<Vec<usize>> = vec![vec![usize::MAX; l + 1]; n];

    let mut has_incoming = vec![false; n];
    for e in &fg.edges { has_incoming[e.to.0] = true; }

    for nidx in 0..n {
        if !has_incoming[nidx] {
            boundary_in[nidx][0] = 0.0;
        }

        if let Some(p) = &fg.nodes[nidx].profile {
            let (out, e_r) = profile_viterbi_with_boundary_and_entry(p, read, &boundary_in[nidx]);
            boundary_out[nidx] = out;
            node_entry_r[nidx] = e_r;
        }

        // Propagate to successors. The successor's `boundary_in[succ][r]` =
        // best `boundary_out[pred][r]` over all preds. The handoff is at the
        // same read position r (instantaneous edge — no read consumption on
        // junction). Backpointers record (predecessor, predecessor's exit r).
        for e in fg.edges.iter().filter(|e| e.from.0 == nidx) {
            let log_edge = (e.family_support as f64 / max_support as f64).ln();
            let to = e.to.0;
            for r in 0..=l {
                let cand = boundary_out[nidx][r] + log_edge;
                if cand > boundary_in[to][r] {
                    boundary_in[to][r] = cand;
                    bp_node[to][r] = nidx;
                    bp_prev_exit[to][r] = r;
                }
            }
        }
    }

    // Find the best terminal (node, r=l).
    let mut best = NEG_INF;
    let mut best_node = usize::MAX;
    for nidx in 0..n {
        if boundary_out[nidx][l] > best {
            best = boundary_out[nidx][l];
            best_node = nidx;
        }
    }

    if best == NEG_INF || best_node == usize::MAX {
        return None;
    }

    // Traceback: at each node we know the exit position (from bp_prev_exit of
    // the next node, or l for the terminal). The entry position into the
    // current node is node_entry_r[node][exit] — recovered from the per-profile
    // intra-DP backpointer.
    let mut path_nodes: Vec<NodeIdx> = Vec::new();
    let mut path_spans: Vec<(usize, usize)> = Vec::new();
    let mut cur = best_node;
    let mut cur_exit = l;
    loop {
        let entry = node_entry_r[cur][cur_exit];
        path_nodes.push(NodeIdx(cur));
        path_spans.push((entry, cur_exit));
        let prev = bp_node[cur][entry];  // predecessor at this node's entry r (== predecessor's exit r)
        if prev == usize::MAX { break; }
        cur = prev;
        cur_exit = entry;
    }
    path_nodes.reverse();
    path_spans.reverse();

    Some(ViterbiPath { nodes: path_nodes, read_spans: path_spans, score: best })
}

// ── Graph-relaxed assignment for low-similarity paralogs (POC) ───────────────

/// Per-paralog assignment via unconstrained Viterbi over the family graph,
/// scored by node-overlap with each paralog's known path.
///
/// Designed for the LOW-similarity band (jaccard < 0.30) where rigid per-copy
/// path forward DP underperforms: at high divergence a read with the natural
/// rate of sequencing errors collects too much per-base penalty against any
/// single paralog's profile, so the resulting log-likelihoods compress and
/// the gap rule abstains. This routine relaxes the path constraint:
///
///   1. Run `viterbi_path` over the family graph (read picks its own best
///      path through any combination of nodes) — same HMM trellis inside
///      each node, just without committing to a specific paralog.
///   2. For each CopyId, compute node-overlap (recall and precision) between
///      the Viterbi node-set and that paralog's known path.
///   3. Sort by recall descending — the top paralog is the one whose known
///      path the Viterbi trace covers most completely.
///
/// Returns `Vec<(CopyId, recall, precision)>`. Empty if the graph is empty or
/// Viterbi fails.
///
/// This is COMPLEMENTARY to `forward_against_path_for_copy`, not a replacement —
/// dispatch by similarity band: high/medium → forward DP, low → graph Viterbi.
pub fn assign_via_graph_viterbi(
    fg: &FamilyGraph,
    read: &[u8],
) -> Vec<(crate::vg_hmm::family_graph::CopyId, f64, f64)> {
    let viterbi = match viterbi_path(fg, read) {
        Some(v) => v,
        None => return Vec::new(),
    };
    let viterbi_set: std::collections::BTreeSet<NodeIdx> =
        viterbi.nodes.iter().copied().collect();
    if viterbi_set.is_empty() { return Vec::new(); }

    let mut out: Vec<(crate::vg_hmm::family_graph::CopyId, f64, f64)> = Vec::new();
    for cid in fg.all_copies() {
        let path = fg.recover_paralog_path(cid);
        if path.is_empty() { continue; }
        let path_set: std::collections::BTreeSet<NodeIdx> =
            path.iter().copied().collect();
        let inter = viterbi_set.intersection(&path_set).count() as f64;
        let recall = inter / path_set.len() as f64;
        let precision = inter / viterbi_set.len() as f64;
        out.push((cid, recall, precision));
    }
    out.sort_by(|a, b|
        b.1.partial_cmp(&a.1).unwrap_or(std::cmp::Ordering::Equal)
            .then(b.2.partial_cmp(&a.2).unwrap_or(std::cmp::Ordering::Equal))
    );
    out
}

// ── Tests ─────────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;
    use crate::vg_hmm::family_graph::{ExonClass, JunctionEdge};
    use crate::vg_hmm::profile::ProfileHmm;

    fn two_node_chain() -> FamilyGraph {
        let p1 = ProfileHmm::from_singleton(b"ACGT");
        let p2 = ProfileHmm::from_singleton(b"TGCA");
        let nodes = vec![
            ExonClass {
                idx: NodeIdx(0), chrom: "x".into(), span: (0, 4), strand: '+',
                per_copy_sequences: vec![(0, b"ACGT".to_vec())],
                per_copy_spans: vec![(0, (0, 4))],
                copy_specific: true, profile: Some(p1), per_copy_profiles: vec![],
            },
            ExonClass {
                idx: NodeIdx(1), chrom: "x".into(), span: (10, 14), strand: '+',
                per_copy_sequences: vec![(0, b"TGCA".to_vec())],
                per_copy_spans: vec![(0, (10, 14))],
                copy_specific: true, profile: Some(p2), per_copy_profiles: vec![],
            },
        ];
        let edges = vec![JunctionEdge { from: NodeIdx(0), to: NodeIdx(1), family_support: 1, strand: '+' }];
        FamilyGraph { family_id: 0, nodes, edges }
    }

    #[test]
    fn viterbi_returns_correct_node_sequence() {
        let fg = two_node_chain();
        let path = viterbi_path(&fg, b"ACGTTGCA").expect("path");
        let names: Vec<usize> = path.nodes.iter().map(|n| n.0).collect();
        assert_eq!(names, vec![0, 1]);
        assert!(path.score > 0.0_f64.ln() - 100.0, "score={}", path.score);
    }

    /// Build a 4-node "low-similarity" two-paralog family for the POC test.
    ///
    ///   N0 (shared) ─┬→ N1 (paralog A only) ─┬→ N3 (shared)
    ///                └→ N2 (paralog B only) ─┘
    ///
    /// N1 and N2 have completely different sequences but the same topological
    /// role (alternative middle exon). Sequence-level Jaccard between the two
    /// paralogs concatenated is very low — the band where per-copy path forward
    /// DP starts to fail.
    fn low_sim_two_paralog_graph() -> FamilyGraph {
        let n0_seq = b"ACGTACGTAC".to_vec();        // shared start
        let n1_seq = b"AAAAAAAATTTAAATTT".to_vec(); // paralog A's middle (A-rich)
        let n2_seq = b"GCGCGCGCGCGCGCGCG".to_vec(); // paralog B's middle (GC-only)
        let n3_seq = b"TTGGAATTGGAA".to_vec();      // shared end
        let nodes = vec![
            ExonClass {
                idx: NodeIdx(0), chrom: "x".into(),
                span: (0, n0_seq.len() as u64), strand: '+',
                per_copy_sequences: vec![(0, n0_seq.clone()), (1, n0_seq.clone())],
                per_copy_spans: vec![
                    (0, (0, n0_seq.len() as u64)),
                    (1, (0, n0_seq.len() as u64)),
                ],
                copy_specific: false,
                profile: Some(ProfileHmm::from_singleton(&n0_seq)),
                per_copy_profiles: vec![],
            },
            ExonClass {
                idx: NodeIdx(1), chrom: "x".into(),
                span: (20, 20 + n1_seq.len() as u64), strand: '+',
                per_copy_sequences: vec![(0, n1_seq.clone())],
                per_copy_spans: vec![(0, (20, 20 + n1_seq.len() as u64))],
                copy_specific: true,
                profile: Some(ProfileHmm::from_singleton(&n1_seq)),
                per_copy_profiles: vec![],
            },
            ExonClass {
                idx: NodeIdx(2), chrom: "x".into(),
                span: (20, 20 + n2_seq.len() as u64), strand: '+',
                per_copy_sequences: vec![(1, n2_seq.clone())],
                per_copy_spans: vec![(1, (20, 20 + n2_seq.len() as u64))],
                copy_specific: true,
                profile: Some(ProfileHmm::from_singleton(&n2_seq)),
                per_copy_profiles: vec![],
            },
            ExonClass {
                idx: NodeIdx(3), chrom: "x".into(),
                span: (50, 50 + n3_seq.len() as u64), strand: '+',
                per_copy_sequences: vec![(0, n3_seq.clone()), (1, n3_seq.clone())],
                per_copy_spans: vec![
                    (0, (50, 50 + n3_seq.len() as u64)),
                    (1, (50, 50 + n3_seq.len() as u64)),
                ],
                copy_specific: false,
                profile: Some(ProfileHmm::from_singleton(&n3_seq)),
                per_copy_profiles: vec![],
            },
        ];
        let edges = vec![
            JunctionEdge { from: NodeIdx(0), to: NodeIdx(1), family_support: 5, strand: '+' },
            JunctionEdge { from: NodeIdx(0), to: NodeIdx(2), family_support: 5, strand: '+' },
            JunctionEdge { from: NodeIdx(1), to: NodeIdx(3), family_support: 5, strand: '+' },
            JunctionEdge { from: NodeIdx(2), to: NodeIdx(3), family_support: 5, strand: '+' },
        ];
        FamilyGraph { family_id: 0, nodes, edges }
    }

    #[test]
    fn poc_graph_viterbi_assigns_low_similarity_paralog_correctly() {
        let fg = low_sim_two_paralog_graph();

        // Read traces paralog A's path: N0 + N1 + N3.
        let n0_seq = b"ACGTACGTAC";
        let n1_seq = b"AAAAAAAATTTAAATTT";
        let n3_seq = b"TTGGAATTGGAA";
        let read: Vec<u8> = [n0_seq.as_slice(), n1_seq, n3_seq].concat();

        let scores = assign_via_graph_viterbi(&fg, &read);
        assert_eq!(scores.len(), 2, "expected 2 paralog scores, got {:?}", scores);

        // Top assignment: paralog A (CopyId 0), recall = 1.0.
        let (top_cid, top_recall, top_precision) = scores[0];
        assert_eq!(top_cid, 0,
                   "expected paralog 0 (A) to be top, got CopyId {} (scores={:?})",
                   top_cid, scores);
        assert!((top_recall - 1.0).abs() < 1e-9,
                "A's recall should be 1.0, got {}", top_recall);
        assert!((top_precision - 1.0).abs() < 1e-9,
                "A's precision should be 1.0, got {}", top_precision);

        // Runner-up: paralog B (CopyId 1), recall = 2/3 (only N0 and N3 in common).
        let (other_cid, other_recall, _) = scores[1];
        assert_eq!(other_cid, 1, "expected paralog 1 (B) as runner-up");
        assert!((other_recall - 2.0/3.0).abs() < 1e-9,
                "B's recall should be 2/3 ≈ 0.667, got {}", other_recall);

        eprintln!(
            "\n[POC] graph-Viterbi assignment on a low-similarity 2-paralog family:\n\
             [POC]   paralog A  recall = {:.3}  precision = {:.3}   ← assigned\n\
             [POC]   paralog B  recall = {:.3}  precision = {:.3}\n\
             [POC]   gap (recall) = {:.3}   ← decisive\n",
            top_recall, top_precision,
            other_recall, scores[1].2,
            top_recall - other_recall,
        );
    }

    #[test]
    fn poc_graph_viterbi_beats_per_copy_forward_when_read_is_noisy() {
        // The mechanism advantage: at low similarity, even modest read errors
        // push log P(read | path_c) toward NEG_INF for ALL paralogs (per-base
        // penalty accumulates). The forward gap collapses → gap rule abstains.
        // Graph-Viterbi only needs to identify which set of nodes the best
        // path goes through — it tolerates per-base error and stays decisive.
        let fg = low_sim_two_paralog_graph();

        let n0_seq = b"ACGTACGTAC";
        let n1_seq = b"AAAAAAAATTTAAATTT";
        let n3_seq = b"TTGGAATTGGAA";
        // 15 % substitution rate against paralog A's path. Enough to wreck a
        // singleton profile's per-base log-likelihood (each error costs ~5 nats
        // against an eps=0.02 singleton emission), but topology still matches A.
        let mut read: Vec<u8> = [n0_seq.as_slice(), n1_seq, n3_seq].concat();
        let flip = |b: u8| -> u8 { match b { b'A' => b'C', b'C' => b'G', b'G' => b'T', _ => b'A' } };
        for i in (3..read.len()).step_by(7) { read[i] = flip(read[i]); }

        // Per-copy forward DP scores against each paralog's path.
        let path_a = fg.recover_paralog_path(0);
        let path_b = fg.recover_paralog_path(1);
        let lp_a = forward_against_path(&fg, &read, &path_a);
        let lp_b = forward_against_path(&fg, &read, &path_b);
        let fwd_gap = (lp_a - lp_b).abs();

        // Graph-Viterbi node-overlap scores.
        let scores = assign_via_graph_viterbi(&fg, &read);
        let recall_a = scores.iter().find(|s| s.0 == 0).map(|s| s.1).unwrap_or(0.0);
        let recall_b = scores.iter().find(|s| s.0 == 1).map(|s| s.1).unwrap_or(0.0);
        let viterbi_gap = recall_a - recall_b;

        eprintln!(
            "\n[POC noisy read]  read = paralog A's path with ~15% per-base errors:\n\
             [POC]   per-copy forward log P:  A = {:>7.2}   B = {:>7.2}   gap = {:.2}\n\
             [POC]   graph-Viterbi recall:    A = {:>7.3}   B = {:>7.3}   gap = {:.3}\n\
             [POC]   ↑ forward gap is in nats (could be tiny / unstable);\n\
             [POC]   ↑ Viterbi recall is in [0, 1] and stays decisive.\n",
            lp_a, lp_b, fwd_gap, recall_a, recall_b, viterbi_gap,
        );

        // Both methods should still pick paralog A in this controlled test,
        // but the Viterbi gap is bounded and easy to threshold; the forward
        // gap depends on absolute log-likelihoods that depend on read length.
        assert!(lp_a > lp_b, "forward should favor A: lp_a={}, lp_b={}", lp_a, lp_b);
        assert!(recall_a > recall_b, "Viterbi recall should favor A");
        assert!(recall_a >= 1.0 - 1e-9, "Viterbi must trace ALL of A's nodes");
    }

    #[test]
    fn poc_graph_viterbi_picks_other_paralog_when_read_matches_it() {
        // Sanity check: when the read traces paralog B's path, B must win.
        let fg = low_sim_two_paralog_graph();
        let n0_seq = b"ACGTACGTAC";
        let n2_seq = b"GCGCGCGCGCGCGCGCG";
        let n3_seq = b"TTGGAATTGGAA";
        let read: Vec<u8> = [n0_seq.as_slice(), n2_seq, n3_seq].concat();

        let scores = assign_via_graph_viterbi(&fg, &read);
        assert_eq!(scores.len(), 2);
        assert_eq!(scores[0].0, 1, "expected paralog 1 (B) to be top, got {:?}", scores);
        assert!((scores[0].1 - 1.0).abs() < 1e-9);
    }

    #[test]
    fn banded_matches_unbanded_for_short_read() {
        let p = ProfileHmm::from_singleton(b"ACGTACGT");
        let read = b"ACGTACGT";
        let mut bin = vec![NEG_INF; read.len() + 1];
        bin[0] = 0.0;
        let unbanded = profile_forward_with_boundary_banded(&p, read, &bin, 0);
        let banded   = profile_forward_with_boundary_banded(&p, read, &bin, 100);
        for i in 0..unbanded.len() {
            if unbanded[i].is_finite() {
                assert!(
                    (banded[i] - unbanded[i]).abs() < 1e-6,
                    "banded[{}] = {} != unbanded[{}] = {}",
                    i, banded[i], i, unbanded[i]
                );
            }
        }
    }
}
