//! Forward / Viterbi over the family graph.
//!
//! Task 3.2 (single-DP boundary threading) + Task 3.3 (Viterbi best-path).

use crate::vg_hmm::profile::ProfileHmm;

#[inline]
fn idx(b: u8) -> Option<usize> {
    match b { b'A' => Some(0), b'C' => Some(1), b'G' => Some(2), b'T' => Some(3), _ => None }
}

const NEG_INF: f64 = f64::NEG_INFINITY;

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
    let m = p.n_columns;
    let l = read.len();
    assert_eq!(boundary_in.len(), l + 1, "boundary_in length mismatch");

    if m == 0 {
        return vec![NEG_INF; l + 1];
    }

    // For banding: find the range of active read positions in boundary_in.
    // The diagonal for entry r at column j is i = r + j.
    // The allowed i-range at column j is [r_min + j - W, r_max + j + W].
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
            // No active entries — entire output is NEG_INF.
            return vec![NEG_INF; l + 1];
        }
        (lo, hi)
    };

    // Rolling 2-row buffers: prev_m/i/d = row j-1, cur_m/i/d = row j.
    // Size is l+1 (full read length) to allow direct indexing by i.
    // We fill only the band cells; the rest stay NEG_INF.
    let mut prev_m: Vec<f64> = boundary_in.to_vec(); // row 0 = boundary_in
    let mut prev_i: Vec<f64> = vec![NEG_INF; l + 1];
    let mut prev_d: Vec<f64> = vec![NEG_INF; l + 1];
    let mut cur_m:  Vec<f64> = vec![NEG_INF; l + 1];
    let mut cur_i:  Vec<f64> = vec![NEG_INF; l + 1];
    let mut cur_d:  Vec<f64> = vec![NEG_INF; l + 1];

    // We need the final row (j == m) in full to build boundary_out.
    // With rolling buffers we naturally have the last row in cur_{m,i,d} after the loop.

    for j in 1..=m {
        // Reset cur row to NEG_INF only over the band we're about to write.
        // (Previous iterations may have written to adjacent bands; clear them.)
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

        // Clear current column buffer for the range we'll write.
        for i in col_lo..=col_hi {
            cur_m[i] = NEG_INF;
            cur_i[i] = NEG_INF;
            cur_d[i] = NEG_INF;
        }

        for i in col_lo..=col_hi {
            // ── Match state ──────────────────────────────────────────────────
            let emit_m = if i > 0 {
                idx(read[i - 1]).map(|k| p.match_emit[j - 1][k]).unwrap_or(NEG_INF)
            } else { NEG_INF };
            let pm_prev = if i > 0 { prev_m[i - 1] } else { NEG_INF };
            let pi_prev = if i > 0 { prev_i[i - 1] } else { NEG_INF };
            let pd_prev = if i > 0 { prev_d[i - 1] } else { NEG_INF };
            let from_m = if pm_prev > NEG_INF { pm_prev + p.trans[j - 1].mm + emit_m } else { NEG_INF };
            let from_i = if pi_prev > NEG_INF { pi_prev + p.trans[j - 1].im + emit_m } else { NEG_INF };
            let from_d = if pd_prev > NEG_INF { pd_prev + p.trans[j - 1].dm + emit_m } else { NEG_INF };
            cur_m[i] = logsumexp(logsumexp(from_m, from_i), from_d);

            // ── Insert state ─────────────────────────────────────────────────
            let emit_i = if i > 0 {
                idx(read[i - 1]).map(|k| p.insert_emit[j][k]).unwrap_or(NEG_INF)
            } else { NEG_INF };
            let i_from_m = if i > 0 && cur_m[i - 1] > NEG_INF { cur_m[i - 1] + p.trans[j].mi + emit_i } else { NEG_INF };
            let i_from_i = if i > 0 && cur_i[i - 1] > NEG_INF { cur_i[i - 1] + p.trans[j].ii + emit_i } else { NEG_INF };
            cur_i[i] = logsumexp(i_from_m, i_from_i);

            // ── Delete state ─────────────────────────────────────────────────
            let d_from_m = if prev_m[i] > NEG_INF { prev_m[i] + p.trans[j - 1].md } else { NEG_INF };
            let d_from_d = if prev_d[i] > NEG_INF { prev_d[i] + p.trans[j - 1].dd } else { NEG_INF };
            cur_d[i] = logsumexp(d_from_m, d_from_d);
        }

        // Clear the previous band cells so stale data doesn't affect future columns.
        // (Only cells not covered by the new band need clearing, and only in [0, l].)
        if bandwidth > 0 {
            let clear_lo = prev_col_lo.min(l);
            let clear_hi = col_lo.min(l);
            for i in clear_lo..clear_hi {
                prev_m[i] = NEG_INF;
                prev_i[i] = NEG_INF;
                prev_d[i] = NEG_INF;
            }
        }

        // Swap rolling buffers: cur → prev for next column.
        std::mem::swap(&mut prev_m, &mut cur_m);
        std::mem::swap(&mut prev_i, &mut cur_i);
        std::mem::swap(&mut prev_d, &mut cur_d);
    }

    // After the loop, prev_* holds row m (the last column).
    // boundary_out[r] = log P exiting after consuming read[0..r] total.
    let mut out = vec![NEG_INF; l + 1];
    for r in 0..=l {
        out[r] = logsumexp(logsumexp(prev_m[r], prev_i[r]), prev_d[r]);
    }
    out
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

    let mut boundary = vec![NEG_INF; l + 1];
    boundary[0] = 0.0;

    for &nidx in path {
        if nidx.0 >= fg.nodes.len() { return NEG_INF; }
        let node = &fg.nodes[nidx.0];
        match &node.profile {
            Some(profile) => {
                boundary = profile_forward_with_boundary(profile, read, &boundary);
            }
            None => return NEG_INF,
        }
        // No edge log-probability is added: the path is given (deterministic
        // topology). Including edge weight here would double-count the path
        // identity, which the caller has already chosen.
    }

    boundary[l]
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
                copy_specific: true, profile: Some(p1),
            },
            ExonClass {
                idx: NodeIdx(1), chrom: "x".into(), span: (10, 14), strand: '+',
                per_copy_sequences: vec![(0, b"TGCA".to_vec())],
                per_copy_spans: vec![(0, (10, 14))],
                copy_specific: true, profile: Some(p2),
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
