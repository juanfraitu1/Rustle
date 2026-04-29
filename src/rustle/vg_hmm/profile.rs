//! Per-exon-class profile HMM (Krogh-style: Match / Insert / Delete states
//! with column-wise emissions). For multi-copy classes, columns come from
//! a partial-order alignment of the per-copy sequences. Singleton classes
//! get a degenerate 1-hot profile.

use anyhow::{anyhow, Result};

/// Log-probability for each base A/C/G/T (in that order).
pub type LogEmit = [f64; 4];

/// Transitions out of a single column. log P over (M->M, M->I, M->D, I->M, I->I, D->M, D->D).
#[derive(Debug, Clone, Copy, Default)]
pub struct TransRow {
    pub mm: f64, pub mi: f64, pub md: f64,
    pub im: f64, pub ii: f64,
    pub dm: f64, pub dd: f64,
}

#[derive(Debug, Clone)]
pub struct ProfileHmm {
    pub n_columns: usize,
    pub match_emit: Vec<LogEmit>,
    pub insert_emit: Vec<LogEmit>,
    pub trans: Vec<TransRow>,
    pub n_seed_copies: usize,
}

const SMOOTH_FLOOR: f64 = 0.001;

fn idx(b: u8) -> Option<usize> {
    match b { b'A' => Some(0), b'C' => Some(1), b'G' => Some(2), b'T' => Some(3), _ => None }
}

impl ProfileHmm {
    pub fn empty(n_columns: usize) -> Self {
        let bg = [(0.25_f64).ln(); 4];
        Self {
            n_columns,
            match_emit: vec![bg; n_columns],
            insert_emit: vec![bg; n_columns + 1],
            trans: vec![TransRow::default(); n_columns + 1],
            n_seed_copies: 0,
        }
    }

    /// Singleton: one copy → 1-hot match emissions with smoothing floor.
    pub fn from_singleton(seq: &[u8]) -> Self {
        let n = seq.len();
        let mut p = Self::empty(n);
        let log_main = (1.0 - 3.0 * SMOOTH_FLOOR).ln();
        let log_floor = SMOOTH_FLOOR.ln();
        for (i, &b) in seq.iter().enumerate() {
            let mut row = [log_floor; 4];
            if let Some(k) = idx(b) { row[k] = log_main; }
            p.match_emit[i] = row;
        }
        p.n_seed_copies = 1;
        p
    }

    /// Build a profile HMM from an equal-length MSA (rows = sequences,
    /// columns = aligned positions, b'-' for gaps).
    pub fn from_msa(rows: &[Vec<u8>]) -> Result<Self> {
        if rows.is_empty() { return Err(anyhow!("from_msa: no rows")); }
        let n_seq = rows.len();
        let n_col = rows[0].len();
        if !rows.iter().all(|r| r.len() == n_col) {
            return Err(anyhow!("from_msa: rows have unequal lengths"));
        }

        // Per-row gap fraction; columns where >= 50% of rows have a gap become
        // INSERT columns (no match state); others are match columns.
        let mut is_match: Vec<bool> = vec![false; n_col];
        for c in 0..n_col {
            let gaps = rows.iter().filter(|r| r[c] == b'-').count();
            is_match[c] = gaps * 2 < n_seq;
        }

        // Background composition (family-wide), used as the smoothing prior.
        let mut bg_counts = [0u32; 4];
        for r in rows {
            for &b in r { if let Some(k) = idx(b) { bg_counts[k] += 1; } }
        }
        let bg_total: u32 = bg_counts.iter().sum::<u32>().max(1);
        let bg: [f64; 4] = [
            bg_counts[0] as f64 / bg_total as f64,
            bg_counts[1] as f64 / bg_total as f64,
            bg_counts[2] as f64 / bg_total as f64,
            bg_counts[3] as f64 / bg_total as f64,
        ];

        // Smoothing weight: alpha = 1 / (1 + n_seq) — fewer copies → lean on prior.
        let alpha = 1.0 / (1.0 + n_seq as f64);

        // Build only match columns; record their original column index for trans fitting.
        let match_col_idx: Vec<usize> = (0..n_col).filter(|&c| is_match[c]).collect();
        let n_match = match_col_idx.len();

        let mut match_emit: Vec<LogEmit> = Vec::with_capacity(n_match);
        for &c in &match_col_idx {
            let mut counts = [0u32; 4];
            let mut observed = 0u32;
            for r in rows {
                if let Some(k) = idx(r[c]) { counts[k] += 1; observed += 1; }
            }
            let obs_total = observed.max(1) as f64;
            let mut row_p = [0.0_f64; 4];
            for k in 0..4 {
                let emp = counts[k] as f64 / obs_total;
                let p = alpha * bg[k] + (1.0 - alpha) * emp;
                row_p[k] = p.max(SMOOTH_FLOOR);
            }
            // Renormalize and log.
            let z: f64 = row_p.iter().sum();
            let log_row = [
                (row_p[0] / z).ln(), (row_p[1] / z).ln(),
                (row_p[2] / z).ln(), (row_p[3] / z).ln(),
            ];
            match_emit.push(log_row);
        }

        // Insert emissions: family background, log-space.
        let log_bg = [
            (bg[0].max(SMOOTH_FLOOR)).ln(),
            (bg[1].max(SMOOTH_FLOOR)).ln(),
            (bg[2].max(SMOOTH_FLOOR)).ln(),
            (bg[3].max(SMOOTH_FLOOR)).ln(),
        ];
        let insert_emit = vec![log_bg; n_match + 1];

        // Fit transitions by counting per-row state walks across columns.
        // State at column c (for row r):
        //   Match  if r[c] != '-' AND is_match[c]
        //   Delete if r[c] == '-' AND is_match[c]
        //   Insert if !is_match[c]   (we don't track per-insert sub-states)
        //
        // Transitions out of MATCH column m_idx are accumulated by walking the
        // row from that column to the next match column; intervening insert
        // columns count as M->I->...->I->M (one M->I out, then I->I within,
        // and I->M into the next match). Gaps at the next match column count
        // M->D or D->D depending on whether the *current* column is M or D.
        let mut counts: Vec<[u32; 7]> = vec![[0; 7]; n_match + 1]; // M->M, M->I, M->D, I->M, I->I, D->M, D->D

        for row in rows {
            // Walk left->right tracking previous match-column state.
            let mut prev_match_state: Option<bool> = None; // Some(true) = match, Some(false) = delete, None = before-first-match
            let mut had_insert_after_prev = false;
            let mut prev_match_idx_in_match: Option<usize> = None;
            for c in 0..n_col {
                if is_match[c] {
                    let cur_is_delete = row[c] == b'-';
                    let mi = match_col_idx.iter().position(|&x| x == c).unwrap();
                    if let Some(prev_was_match) = prev_match_state {
                        let prev_mi = prev_match_idx_in_match.unwrap();
                        if had_insert_after_prev {
                            // prev → I → ... → I → cur
                            if prev_was_match { counts[prev_mi][1] += 1; } // M->I
                            // I->I edges within would need to be counted; skip for simplicity (single insert run).
                            if cur_is_delete { counts[prev_mi + 1][3] += 0; /* I->D not modelled */ }
                            else { counts[prev_mi][3] += 1; } // I->M (counted at prev's row for simplicity)
                        } else {
                            if prev_was_match {
                                if cur_is_delete { counts[prev_mi][2] += 1; } // M->D
                                else { counts[prev_mi][0] += 1; }             // M->M
                            } else {
                                if cur_is_delete { counts[prev_mi][6] += 1; } // D->D
                                else { counts[prev_mi][5] += 1; }             // D->M
                            }
                        }
                    }
                    prev_match_state = Some(!cur_is_delete);
                    prev_match_idx_in_match = Some(mi);
                    had_insert_after_prev = false;
                } else if row[c] != b'-' {
                    had_insert_after_prev = true;
                }
            }
        }

        // Convert counts → smoothed log-probabilities per row.
        let mut trans: Vec<TransRow> = Vec::with_capacity(n_match + 1);
        for c in &counts {
            let m_total = (c[0] + c[1] + c[2]).max(1) as f64;
            let i_total = (c[3] + c[4]).max(1) as f64;
            let d_total = (c[5] + c[6]).max(1) as f64;
            let s = SMOOTH_FLOOR;
            let mm = ((c[0] as f64 + s) / (m_total + 3.0 * s)).ln();
            let mi = ((c[1] as f64 + s) / (m_total + 3.0 * s)).ln();
            let md = ((c[2] as f64 + s) / (m_total + 3.0 * s)).ln();
            let im = ((c[3] as f64 + s) / (i_total + 2.0 * s)).ln();
            let ii = ((c[4] as f64 + s) / (i_total + 2.0 * s)).ln();
            let dm = ((c[5] as f64 + s) / (d_total + 2.0 * s)).ln();
            let dd = ((c[6] as f64 + s) / (d_total + 2.0 * s)).ln();
            trans.push(TransRow { mm, mi, md, im, ii, dm, dd });
        }

        Ok(Self {
            n_columns: n_match,
            match_emit,
            insert_emit,
            trans,
            n_seed_copies: n_seq,
        })
    }
}

/// Multiple sequence alignment via progressive Needleman-Wunsch.
/// Returns rows of equal length, padded with b'-' for gaps.
/// One row per input sequence, in the same order.
/// Errors if fewer than 2 sequences provided.
///
/// Note: poasta's API is graph-based (POAGraph + PoastaAligner) and does not
/// expose a direct MSA-rows output. Per the plan's §12 risk mitigation, we
/// implement a pure-Rust banded progressive NW here. For the exon sizes in
/// scope (<300 bp), this is trivially fast.
pub fn poa_msa(seqs: &[Vec<u8>]) -> Result<Vec<Vec<u8>>> {
    if seqs.len() < 2 {
        return Err(anyhow!("poa_msa requires at least 2 sequences"));
    }
    // Progressive pairwise NW: align seqs[1..] against the growing consensus.
    // We build an MSA incrementally by aligning each new sequence to the
    // current profile-consensus (gapless projection) and inserting gaps.
    let mut msa: Vec<Vec<u8>> = vec![seqs[0].clone()];
    for next_seq in &seqs[1..] {
        // Consensus of current MSA: take most-common non-gap base at each column.
        let consensus = msa_consensus(&msa);
        // Align next_seq to consensus.
        let (aln_cons, aln_next) = nw_align(&consensus, next_seq);
        // Expand existing MSA rows to match aln_cons (insert gap columns where
        // aln_cons has a gap == insertion in next_seq).
        let mut new_msa: Vec<Vec<u8>> = msa.iter().map(|_| Vec::new()).collect();
        // We need to know which columns in the original consensus correspond to
        // which positions in aln_cons.
        let mut orig_pos = 0usize;
        for (&ac, &an) in aln_cons.iter().zip(aln_next.iter()) {
            if ac == b'-' {
                // Insertion in next_seq: all existing rows get a gap here.
                for r in &mut new_msa { r.push(b'-'); }
            } else {
                // Match or deletion in next_seq: copy column orig_pos from each old row.
                for (i, r) in new_msa.iter_mut().enumerate() {
                    r.push(msa[i][orig_pos]);
                }
                orig_pos += 1;
            }
            // next_seq column:
            let _ = an; // already encoded in aln_next
        }
        // Add next_seq row (using aln_next directly, but we need it indexed by
        // the new column layout — aln_next is already padded correctly).
        new_msa.push(aln_next);
        msa = new_msa;
    }
    // Verify equal lengths.
    let l = msa[0].len();
    assert!(msa.iter().all(|r| r.len() == l));
    Ok(msa)
}

/// Compute majority-vote consensus from an MSA (ignoring gaps).
fn msa_consensus(msa: &[Vec<u8>]) -> Vec<u8> {
    if msa.is_empty() { return Vec::new(); }
    let n_col = msa[0].len();
    let mut cons = Vec::with_capacity(n_col);
    for c in 0..n_col {
        let mut counts = [0u32; 256];
        for row in msa {
            let b = row[c];
            if b != b'-' { counts[b as usize] += 1; }
        }
        // Pick the most common non-gap base; default to 'N' if all gaps.
        let best = counts.iter().enumerate()
            .filter(|&(i, _)| i != b'-' as usize)
            .max_by_key(|&(_, &c)| c)
            .map(|(i, _)| i as u8)
            .unwrap_or(b'N');
        cons.push(best);
    }
    cons
}

/// Global Needleman-Wunsch alignment. Returns (aligned_a, aligned_b) with b'-' for gaps.
fn nw_align(a: &[u8], b: &[u8]) -> (Vec<u8>, Vec<u8>) {
    let m = a.len();
    let n = b.len();
    // Scoring: match=1, mismatch=-1, gap=-2.
    const MATCH: i32 = 1;
    const MISMATCH: i32 = -1;
    const GAP: i32 = -2;

    // DP matrix.
    let mut dp = vec![vec![0i32; n + 1]; m + 1];
    for i in 0..=m { dp[i][0] = i as i32 * GAP; }
    for j in 0..=n { dp[0][j] = j as i32 * GAP; }

    for i in 1..=m {
        for j in 1..=n {
            let s = if a[i-1] == b[j-1] { MATCH } else { MISMATCH };
            dp[i][j] = (dp[i-1][j-1] + s)
                .max(dp[i-1][j] + GAP)
                .max(dp[i][j-1] + GAP);
        }
    }

    // Traceback.
    let mut aa = Vec::new();
    let mut bb = Vec::new();
    let mut i = m;
    let mut j = n;
    while i > 0 || j > 0 {
        if i > 0 && j > 0 {
            let s = if a[i-1] == b[j-1] { MATCH } else { MISMATCH };
            if dp[i][j] == dp[i-1][j-1] + s {
                aa.push(a[i-1]); bb.push(b[j-1]);
                i -= 1; j -= 1;
                continue;
            }
        }
        if i > 0 && dp[i][j] == dp[i-1][j] + GAP {
            aa.push(a[i-1]); bb.push(b'-');
            i -= 1;
        } else {
            aa.push(b'-'); bb.push(b[j-1]);
            j -= 1;
        }
    }
    aa.reverse();
    bb.reverse();
    (aa, bb)
}
