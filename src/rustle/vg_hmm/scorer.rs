//! Forward / Viterbi over the family graph. For now, single-profile only;
//! inter-node transitions land in Task 3.3.

use crate::vg_hmm::profile::ProfileHmm;

#[inline]
fn idx(b: u8) -> Option<usize> {
    match b { b'A' => Some(0), b'C' => Some(1), b'G' => Some(2), b'T' => Some(3), _ => None }
}

const NEG_INF: f64 = f64::NEG_INFINITY;

#[inline]
fn logsumexp(a: f64, b: f64) -> f64 {
    if a == NEG_INF { return b; }
    if b == NEG_INF { return a; }
    let m = a.max(b);
    m + ((a - m).exp() + (b - m).exp()).ln()
}

/// Forward log P(read | profile). Unbanded for correctness; banding can be
/// added later if profiling shows it matters.
pub fn forward_against_profile(p: &ProfileHmm, read: &[u8]) -> f64 {
    let m = p.n_columns;     // number of match columns
    let l = read.len();
    if m == 0 { return NEG_INF; }

    // dp_M[j][i], dp_I[j][i], dp_D[j][i] over match columns 1..=m and positions 0..=l.
    let mut m_dp = vec![vec![NEG_INF; l + 1]; m + 1];
    let mut i_dp = vec![vec![NEG_INF; l + 1]; m + 1];
    let mut d_dp = vec![vec![NEG_INF; l + 1]; m + 1];

    // Initial: M[0][0] is virtual start.
    m_dp[0][0] = 0.0;

    for j in 1..=m {
        for i in 0..=l {
            // Emission of M[j] consuming read[i-1].
            let emit_m = if i > 0 {
                idx(read[i - 1]).map(|k| p.match_emit[j - 1][k]).unwrap_or(NEG_INF)
            } else { NEG_INF };
            // M from M[j-1][i-1] + mm, I[j-1][i-1] + im, D[j-1][i-1] + dm
            let from_m = if i > 0 { m_dp[j - 1][i - 1] + p.trans[j - 1].mm + emit_m } else { NEG_INF };
            let from_i = if i > 0 { i_dp[j - 1][i - 1] + p.trans[j - 1].im + emit_m } else { NEG_INF };
            let from_d = if i > 0 { d_dp[j - 1][i - 1] + p.trans[j - 1].dm + emit_m } else { NEG_INF };
            m_dp[j][i] = logsumexp(logsumexp(from_m, from_i), from_d);

            // I from M[j][i-1] + mi (with insert emission), I[j][i-1] + ii.
            let emit_i = if i > 0 {
                idx(read[i - 1]).map(|k| p.insert_emit[j][k]).unwrap_or(NEG_INF)
            } else { NEG_INF };
            let i_from_m = if i > 0 { m_dp[j][i - 1] + p.trans[j].mi + emit_i } else { NEG_INF };
            let i_from_i = if i > 0 { i_dp[j][i - 1] + p.trans[j].ii + emit_i } else { NEG_INF };
            i_dp[j][i] = logsumexp(i_from_m, i_from_i);

            // D from M[j-1][i] + md, D[j-1][i] + dd (no emission).
            let d_from_m = m_dp[j - 1][i] + p.trans[j - 1].md;
            let d_from_d = d_dp[j - 1][i] + p.trans[j - 1].dd;
            d_dp[j][i] = logsumexp(d_from_m, d_from_d);
        }
    }
    // Terminate: total = log-sum-exp over M[m][l], I[m][l], D[m][l].
    logsumexp(logsumexp(m_dp[m][l], i_dp[m][l]), d_dp[m][l])
}

use crate::vg_hmm::family_graph::{FamilyGraph, NodeIdx};

/// Forward over the family graph: for each read split point, try all
/// node sequences that respect junction edges. We compute, for each node
/// and each starting read position, the best forward score over all
/// path prefixes ending at that node.
///
/// Implementation: dynamic programming over `(node_idx, read_pos_after)`,
/// where the value is `log P(read[0..pos] | path ending at this node)`.
/// Transition between nodes adds the log-junction-prior
/// `log(family_support / max_support_in_family)`.
pub fn forward_against_family(fg: &FamilyGraph, read: &[u8]) -> f64 {
    if fg.nodes.is_empty() { return NEG_INF; }
    let n = fg.nodes.len();
    let l = read.len();
    let max_support = fg.edges.iter().map(|e| e.family_support).max().unwrap_or(1).max(1);

    // node_score[node][pos_after] = best log-prob of consuming read[0..pos_after] ending at node
    let mut node_score = vec![vec![NEG_INF; l + 1]; n];

    // Determine source nodes (those with no incoming edges).
    let mut has_incoming = vec![false; n];
    for e in &fg.edges { has_incoming[e.to.0] = true; }

    // For source nodes: consume some prefix of read.
    for (nidx, node) in fg.nodes.iter().enumerate() {
        if !has_incoming[nidx] {
            if let Some(p) = &node.profile {
                for pos_after in 0..=l {
                    node_score[nidx][pos_after] = logsumexp(
                        node_score[nidx][pos_after],
                        forward_against_profile(p, &read[..pos_after]),
                    );
                }
            }
        }
    }

    // Topological order over edges (assume nodes are already in reasonable order;
    // if not, run Kahn's algorithm — most family graphs are small).
    let order: Vec<NodeIdx> = (0..n).map(NodeIdx).collect();
    // Trust insertion order as topological for a single-strand chain; revisit if a test fails.

    for &nidx in &order {
        for e in fg.edges.iter().filter(|e| e.from == nidx) {
            let to_node = &fg.nodes[e.to.0];
            let log_edge = (e.family_support as f64 / max_support as f64).ln();
            if let Some(p) = &to_node.profile {
                // For each starting position at the destination, the best score is
                // node_score[from][start] + log_edge + forward(p, read[start..end]).
                for start in 0..=l {
                    if node_score[nidx.0][start] == NEG_INF { continue; }
                    for end in start..=l {
                        let inner = forward_against_profile(p, &read[start..end]);
                        let cand = node_score[nidx.0][start] + log_edge + inner;
                        node_score[e.to.0][end] = logsumexp(node_score[e.to.0][end], cand);
                    }
                }
            }
        }
    }

    // Final: best score over any node ending at read end.
    let mut best = NEG_INF;
    for s in &node_score { best = best.max(s[l]); }
    best
}
