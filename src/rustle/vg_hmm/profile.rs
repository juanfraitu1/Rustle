//! Per-exon-class profile HMM (Krogh-style: Match / Insert / Delete states
//! with column-wise emissions). For multi-copy classes, columns come from
//! a partial-order alignment of the per-copy sequences. Singleton classes
//! get a degenerate 1-hot profile.

use anyhow::{anyhow, Result};
use poasta::graphs::poa::{POAGraph, POANodeIndex};
use poasta::aligner::PoastaAligner;
use poasta::aligner::config::AffineDijkstra;
use poasta::aligner::scoring::{GapAffine, AlignmentType};
use poasta::io::fasta::poa_graph_to_fasta;

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

/// Sensible prior for unobserved transition rows: M→M dominant.
/// All values are log-probabilities; rows must sum (in linear space)
/// to ~1 within their out-group (M-out, I-out, D-out).
fn default_trans_row() -> TransRow {
    TransRow {
        mm: (0.95_f64).ln(), mi: (0.025_f64).ln(), md: (0.025_f64).ln(),
        im: (0.5_f64).ln(),  ii: (0.5_f64).ln(),
        dm: (0.95_f64).ln(), dd: (0.05_f64).ln(),
    }
}

impl ProfileHmm {
    pub fn empty(n_columns: usize) -> Self {
        let bg = [(0.25_f64).ln(); 4];
        Self {
            n_columns,
            match_emit: vec![bg; n_columns],
            insert_emit: vec![bg; n_columns + 1],
            trans: vec![default_trans_row(); n_columns + 1],
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

        // Match-column detection: a column is a match column if a strict majority
        // of rows have a non-gap base (i.e., more than half are non-gap).
        // Columns where ≥50% of rows have gaps are INSERT columns with no match state.
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

        // Smoothing weight: alpha = 1 / (1 + n_seq) as a mixing fraction.
        // Applied as p[k] = alpha * bg[k] + (1 - alpha) * emp[k]
        // where emp[k] = counts[k] / obs_total.
        // Fewer copies → lean more on prior.
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
            // Mixing-fraction smoothing: p[k] = alpha * bg[k] + (1 - alpha) * emp[k].
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
        let mut trans_counts: Vec<[u32; 7]> = vec![[0; 7]; n_match + 1]; // M->M, M->I, M->D, I->M, I->I, D->M, D->D

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
                            if prev_was_match { trans_counts[prev_mi][1] += 1; } // M->I
                            // I->I edges within would need to be counted; skip for simplicity (single insert run).
                            // I->D transitions are not modelled (insert ends in a delete is rare and design defers it).
                            if !cur_is_delete { trans_counts[prev_mi][3] += 1; } // I->M (counted at prev's row for simplicity)
                        } else {
                            if prev_was_match {
                                if cur_is_delete { trans_counts[prev_mi][2] += 1; } // M->D
                                else { trans_counts[prev_mi][0] += 1; }             // M->M
                            } else {
                                if cur_is_delete { trans_counts[prev_mi][6] += 1; } // D->D
                                else { trans_counts[prev_mi][5] += 1; }             // D->M
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
        // Prior: when no observed transitions, bias toward M->M by adding a
        // "virtual" M->M prior count (0.95) vs M->I/D (0.025 each). This ensures
        // terminal and unobserved columns still have M->M dominant.
        const PRIOR_MM: f64 = 0.95;
        const PRIOR_MI: f64 = 0.025;
        const PRIOR_MD: f64 = 0.025;
        let mut trans: Vec<TransRow> = Vec::with_capacity(n_match + 1);
        for c in &trans_counts {
            let observed_m = c[0] + c[1] + c[2];
            let s = SMOOTH_FLOOR;
            let (mm, mi, md) = if observed_m == 0 {
                // No data: use M->M-biased prior.
                ((PRIOR_MM).ln(), (PRIOR_MI).ln(), (PRIOR_MD).ln())
            } else {
                let m_total = observed_m as f64;
                (
                    ((c[0] as f64 + s) / (m_total + 3.0 * s)).ln(),
                    ((c[1] as f64 + s) / (m_total + 3.0 * s)).ln(),
                    ((c[2] as f64 + s) / (m_total + 3.0 * s)).ln(),
                )
            };
            let i_total = (c[3] + c[4]).max(1) as f64;
            let d_total = (c[5] + c[6]).max(1) as f64;
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

/// Multiple sequence alignment via poasta POA graph traversal.
/// Returns rows of equal length, padded with b'-' for gaps.
/// One row per input sequence, in the same order.
/// Errors if fewer than 2 sequences provided.
pub fn poa_msa(seqs: &[Vec<u8>]) -> Result<Vec<Vec<u8>>> {
    if seqs.len() < 2 {
        return Err(anyhow!("poa_msa requires at least 2 sequences"));
    }

    // Build a POA graph by aligning each sequence progressively.
    // Uses affine-gap Dijkstra (no A* heuristic) with mismatch=1, gap_extend=1, gap_open=2.
    let gap_costs = GapAffine::new(1, 1, 2);
    let aligner: PoastaAligner<'_, AffineDijkstra> =
        PoastaAligner::new(AffineDijkstra(gap_costs), AlignmentType::Global);

    let mut graph: POAGraph<u32> = POAGraph::new();
    let unit_weights: Vec<usize> = vec![1; seqs.iter().map(|s| s.len()).max().unwrap_or(0)];

    for (i, seq) in seqs.iter().enumerate() {
        let w = &unit_weights[..seq.len()];
        if graph.is_empty() {
            // First sequence: add unaligned (no existing graph to align to).
            graph
                .add_alignment_with_weights(&format!("seq{i}"), seq, None, w)
                .map_err(|e| anyhow!("poasta add_alignment_with_weights failed: {e}"))?;
        } else {
            // Subsequent sequences: align against the current graph.
            let aln_result = aligner.align::<u32, _>(&graph, seq);
            let alignment = if aln_result.alignment.is_empty() {
                None
            } else {
                Some(&aln_result.alignment as &poasta::aligner::Alignment<POANodeIndex<u32>>)
            };
            graph
                .add_alignment_with_weights(&format!("seq{i}"), seq, alignment, w)
                .map_err(|e| anyhow!("poasta add_alignment_with_weights failed: {e}"))?;
        }
    }

    // Extract the MSA rows by using poasta's public poa_graph_to_fasta function,
    // which walks the graph in topological order and writes aligned FASTA records.
    // We write into a Vec<u8> buffer and parse the resulting FASTA lines.
    let mut fasta_buf: Vec<u8> = Vec::new();
    poa_graph_to_fasta(&graph, &mut fasta_buf)
        .map_err(|e| anyhow!("poa_graph_to_fasta failed: {e}"))?;

    // Parse the FASTA output: each sequence record is a set of lines starting with '>'.
    // We collect the sequences in order and convert them to Vec<u8> rows.
    let mut msa: Vec<Vec<u8>> = Vec::with_capacity(seqs.len());
    let mut current_seq: Option<Vec<u8>> = None;
    for line in fasta_buf.split(|&b| b == b'\n') {
        if line.is_empty() {
            continue;
        }
        if line[0] == b'>' {
            if let Some(seq) = current_seq.take() {
                msa.push(seq);
            }
            current_seq = Some(Vec::new());
        } else if let Some(ref mut seq) = current_seq {
            seq.extend_from_slice(line);
        }
    }
    if let Some(seq) = current_seq {
        msa.push(seq);
    }

    if msa.len() != seqs.len() {
        return Err(anyhow!(
            "poa_msa: expected {} rows from graph but got {}",
            seqs.len(),
            msa.len()
        ));
    }

    // poasta's fasta_aln_for_seq has an off-by-one in the trailing-gap fill for
    // sequences whose path ends before max_col, so rows can occasionally come
    // back one column short. Pad with trailing gaps to the max row length —
    // this is content-neutral for an MSA (gaps are not part of the ungapped
    // sequence) and the round-trip ungap(row) ≡ original input is preserved.
    let n_cols = msa.iter().map(|r| r.len()).max().unwrap_or(0);
    for row in msa.iter_mut() {
        if row.len() < n_cols {
            row.extend(std::iter::repeat(b'-').take(n_cols - row.len()));
        }
    }
    Ok(msa)
}
