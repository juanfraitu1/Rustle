//! Copy ASSIGNMENT (resolution): assign ONE read to a specific paralog COPY, given the family's known
//! copies, using PSV bases + copy-specific JUNCTIONS, behind an IDENTIFIABILITY gate.
//!
//! This is distinct from `copy_split` (which DISCOVERS the copies). Here the copy set is known (from the
//! family layer) and we RESOLVE which copy a read came from -- the "pick one mapping" step for hard
//! multimappers minimap2 leaves at MAPQ 0. A read carries its true copy's alleles regardless of which
//! paralog it aligned to, so matching its observed features against every copy's feature vector resolves it.
//!
//! Identifiability theorem (operational): a read is RESOLVABLE iff it spans >= 1 feature (a PSV column or
//! a junction boundary) where the candidate copies DIFFER. Reads spanning 0 such features are genuinely
//! TIED ("N equally good places"). Among resolvable reads, a log-likelihood-ratio margin over the
//! runner-up calls ASSIGNED vs AMBIGUOUS. Mirrors `bench/copy_assign.py::assign_read`.

use super::copy_split::{allele_at, intron_chain_of, AlignedRead};

/// Per-copy feature profile over the family's PSV columns + intron-boundary set.
#[derive(Clone, Debug)]
pub struct CopyProfile {
    pub copy_id: usize,
    /// Allele base at each family PSV column (index = column). `None` = this copy has a gap there.
    /// A copy that MATCHES the reference allele must carry that base here, not `None` (else the
    /// likelihood would favour copies it shares no columns with).
    pub alleles: Vec<Option<u8>>,
    /// Intron-boundary offsets in transcription-spliced space (a copy's copy-specific junctions).
    pub junctions: Vec<i64>,
}

/// One PSV bubble of the family's ad-hoc-reference variation graph.
#[derive(Clone, Debug)]
pub struct Bubble {
    /// PSV column index (= bubble id). Bubbles are in ascending column order (== the matrix order).
    pub col: usize,
    /// Allele-node each copy PATH visits here (index = copy). `None` = the copy has a gap.
    pub copy_allele: Vec<Option<u8>>,
    /// The copies carry >= 2 distinct non-`None` alleles here (read-independent). Matches the old
    /// per-column `differ` test in `read_copy_evidence`.
    pub decisive: bool,
}

/// The per-family variation graph the O2 decision threads reads through: one bubble per PSV column,
/// each copy a path over the allele-nodes. Built once per family; the ad-hoc, auditable reference.
#[derive(Clone, Debug)]
pub struct BubbleGraph {
    pub bubbles: Vec<Bubble>,
    pub n_copies: usize,
}

impl BubbleGraph {
    /// Build the family's bubble graph from the copy profiles. Deterministic; `decisive` is computed
    /// exactly as `read_copy_evidence`'s inner `differ` loop (>= 2 distinct non-`None` alleles).
    pub fn from_copies(copies: &[CopyProfile]) -> BubbleGraph {
        let n_cols = copies.iter().map(|c| c.alleles.len()).max().unwrap_or(0);
        let mut bubbles = Vec::with_capacity(n_cols);
        for col in 0..n_cols {
            let copy_allele: Vec<Option<u8>> =
                copies.iter().map(|c| c.alleles.get(col).copied().flatten()).collect();
            let mut seen: Option<u8> = None;
            let mut decisive = false;
            for a in copy_allele.iter().flatten() {
                match seen {
                    None => seen = Some(*a),
                    Some(s) => {
                        if s != *a {
                            decisive = true;
                        }
                    }
                }
            }
            bubbles.push(Bubble { col, copy_allele, decisive });
        }
        BubbleGraph { bubbles, n_copies: copies.len() }
    }
}

/// One read's observed features in the family's column/boundary space.
#[derive(Clone, Debug)]
pub struct ReadFeatures {
    /// Observed base per PSV column (`None` = the read does not span that column).
    pub psv_obs: Vec<Option<u8>>,
    /// Phred quality of the observed base per PSV column, parallel to `psv_obs`. `None` (or a
    /// shorter/empty vector) -> the likelihood uses the flat `error_rate` for that column;
    /// `Some(q)` -> the column is weighted by its own per-base error `10^(-q/10)`.
    pub psv_qual: Vec<Option<u8>>,
    /// The read's intron-boundary offsets (transcription-spliced space).
    pub junctions: Vec<i64>,
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum AssignStatus {
    /// Resolvable and the best copy clears the likelihood-ratio margin.
    Assigned,
    /// Resolvable (spans a decisive feature) but the margin over the runner-up is too small.
    Ambiguous,
    /// Spans no feature where the copies differ -> genuinely unassignable.
    Tied,
}

#[derive(Clone, Debug)]
pub struct Assignment {
    pub best_copy: usize,
    pub log_lr_margin: f64,
    /// Number of features (PSV columns + junction boundaries) the read observes where copies differ.
    pub n_decisive: usize,
    pub resolvable: bool,
    pub status: AssignStatus,
    /// IsoCon certificate: P(this assignment by error if the read were the least-distinguishable competitor).
    pub p_value: f64,
    /// Identifiability bound: the best attainable `p_value` (read supports `best` at every distinguishing
    /// position) against the hardest competitor. `>= alpha` ⇒ the read is unresolvable (`Tied`).
    pub min_p_value: f64,
    /// Set when this read was assigned to a reference-ABSENT (collapsed) copy that was admitted only because
    /// the absent-copy discovery feature was on — i.e. the assignment is COUPLED to that discovery and would
    /// not exist in the frozen ref-only result. Always `false` on the default (absent-copy-OFF) path.
    pub discovery_coupled: bool,
    /// Per-copy POSTERIOR over the candidate copies, `softmax(logl)` (likelihood-normalized, i.e. a UNIFORM
    /// prior), indexed parallel to the `copies` slice. For an *assigned* read it is ~one-hot at `best_copy`;
    /// for a *Tied* read it spreads over the consistent ZONE (the copies the read cannot be distinguished
    /// from) — the soft/Bayesian complement to the hard assign/abstain. An informative prior (copy abundance,
    /// DNA parCN) is applied downstream by re-weighting and renormalizing. Empty for `n == 0`.
    pub posterior: Vec<f64>,
}

/// The decisive log-likelihood-ratio margin `τ` for a target per-read misassignment rate `p`.
///
/// `τ(p) = ln((1-p)/p)` is the Bayes-optimal decision threshold: accept the top copy iff its
/// log-LR over the runner-up exceeds `τ`, giving posterior misassignment ≤ `p` on the resolved reads.
/// This makes `margin` a PRINCIPLED operating point (the user chooses `p`), not an arbitrary number —
/// the two values seen in the codebase are simply two choices of `p`:
///   * `p = 1e-3 → τ ≈ 6.90` — conservative, the PSV-space analog of Eichler's AS ≥ 10 (precision mode).
///   * `p ≈ 0.119 → τ = 2.0` — the permissive default inherited from the Python prototype (recall mode).
/// At flat HiFi error the vote-count gate is algebraically equivalent (`τ ≈ votes · ln(3(1-e)/e)`),
/// which is why the production vote engine's `margin = 1` ≡ `τ ≈ 6.9` (kill-test, 16/16).
pub fn tau_from_p(p: f64) -> f64 {
    ((1.0 - p) / p).ln()
}

/// Inverse of [`tau_from_p`]: the target per-read misassignment rate implied by a margin `τ`.
pub fn p_from_tau(tau: f64) -> f64 {
    1.0 / (1.0 + tau.exp())
}

/// Exact upper tail of a Poisson-binomial: `P(Σ Bernoulli(probs[i]) >= k)`. O(n^2) DP convolution of the
/// per-trial success probabilities (the number of trials = distinguishing positions, typically < ~100, so
/// exact is cheap and avoids any normal-approximation error). Conventions: `k == 0 -> 1.0`; `k > probs.len()
/// -> 0.0` (also `k > 0` with no trials). Probabilities are clamped to `[0, 1]`.
pub(crate) fn poisson_binomial_upper_tail(k: usize, probs: &[f64]) -> f64 {
    if k == 0 {
        return 1.0;
    }
    if k > probs.len() {
        return 0.0;
    }
    // dp[s] = P(exactly s successes) after the processed trials. Iterate s downward so dp[s-1] is still the
    // pre-trial value when read.
    let mut dp = vec![0.0f64; probs.len() + 1];
    dp[0] = 1.0;
    for &p in probs {
        let p = p.clamp(0.0, 1.0);
        for s in (0..dp.len()).rev() {
            let from_prev = if s > 0 { dp[s - 1] * p } else { 0.0 };
            dp[s] = dp[s] * (1.0 - p) + from_prev;
        }
    }
    dp[k..].iter().sum()
}

/// Tunable assignment parameters. `margin` is the decisive log-LR threshold `τ` — see [`tau_from_p`];
/// set it via [`AssignParams::for_target_misassignment`] to choose the operating point by `p` directly.
#[derive(Clone, Copy, Debug)]
pub struct AssignParams {
    /// Per-base PSV error rate used in the likelihood (HiFi ~ 0.003).
    pub error_rate: f64,
    /// Per-junction log-odds: a matching copy-specific junction is +this, a non-match is -this.
    pub junction_weight: f64,
    /// Junction-boundary match tolerance (bp) for splice-site jitter.
    pub boundary_tol: i64,
    /// Decisive log-LR threshold `τ` over the runner-up — LEGACY, INERT by default. The production gate uses
    /// the IsoCon significance certificate (`min_p >= alpha`), not this margin; `margin`/`τ` are consulted
    /// ONLY when `use_margin_gate` is set (the A/B comparison). `tau_from_p(p)` stays as the principled
    /// operating-point map (`p → τ`), kept for that comparison — it is not part of the shipped decision.
    pub margin: f64,
    /// Significance level / target per-read false-assignment rate for the IsoCon gate. Default 1e-3.
    pub alpha: f64,
    /// Per-distinguishing-junction error probability ε used in the significance test (junctions are sharp).
    pub junction_err: f64,
    /// Revert to the legacy τ-margin gate (for reproducing legacy numbers / the A/B comparison). Default false.
    pub use_margin_gate: bool,
    /// RNA-editing filter: when on, a distinguishing PSV column flagged as an A-to-I editing site (passed as
    /// `editing_cols` to `assign_read_editing`) gets its εⱼ inflated to `edit_rate` in the significance test,
    /// so an edited base cannot fake copy-support. Default true (acts only on detected-heterogeneous A↔G
    /// columns — genuine A/G paralog SNVs are unaffected).
    pub rna_editing_filter: bool,
    /// εⱼ used for an editing-flagged distinguishing column (the rate at which a true base shows the other
    /// allele by editing rather than sequencing error). Default 0.2 — downweights the flagged column to
    /// near-uninformative for the certificate while leaving the likelihood ranking intact.
    pub edit_rate: f64,
    /// IsoCon-style iterative candidate pruning: after the per-read assignment, repeatedly merge copies that
    /// have no read with significant evidence distinguishing them from their nearest neighbor, reassigning
    /// reads until all surviving copies are defensible. Default false (byte-identical baseline).
    pub iterative_prune: bool,
}

impl Default for AssignParams {
    fn default() -> Self {
        // margin = 2.0 ≡ p ≈ 0.119 (recall mode, the Python-prototype operating point). Kept as the
        // default for behavioral continuity; use `for_target_misassignment(1e-3)` for the conservative
        // τ ≈ 6.9 / Eichler-AS≥10 precision operating point.
        AssignParams {
            error_rate: 0.003,
            junction_weight: 5.0,
            boundary_tol: 4,
            margin: 2.0,
            alpha: 1e-3,
            junction_err: 1e-4,
            use_margin_gate: false,
            rna_editing_filter: true,
            edit_rate: 0.2,
            iterative_prune: false,
        }
    }
}

impl AssignParams {
    /// Construct with the decisive margin set from a target per-read misassignment rate `p` (the
    /// principled knob): `margin = tau_from_p(p)`. e.g. `p = 1e-3` → the conservative τ ≈ 6.9.
    pub fn for_target_misassignment(p: f64) -> Self {
        AssignParams { margin: tau_from_p(p), ..Self::default() }
    }

    /// Construct with the significance level set directly (the IsoCon-gate knob): `alpha = p`.
    pub fn for_alpha(alpha: f64) -> Self {
        AssignParams { alpha, ..Self::default() }
    }
}

pub(crate) fn boundary_present(jb: i64, junctions: &[i64], tol: i64) -> bool {
    junctions.iter().any(|&x| (jb - x).abs() <= tol)
}

/// Pairwise IsoCon significance: probability that `read`'s evidence supporting `target` over
/// `competitor` arose by error under H0 "read came from competitor". Also returns the best
/// attainable p-value (product of per-trial error probabilities over the distinguishing observations).
/// This is the same computation used inside [`assign_read_editing`], exposed for iterative copy-pruning
/// and other per-pair analyses.
pub fn copy_pair_significance(
    read: &ReadFeatures,
    target: &CopyProfile,
    competitor: &CopyProfile,
    p: &AssignParams,
    editing_cols: &[bool],
) -> (f64, f64) {
    let spanned: Vec<usize> = (0..read.psv_obs.len()).filter(|&j| read.psv_obs[j].is_some()).collect();
    let mut eps: Vec<f64> = Vec::new();
    let mut k = 0usize;
    // distinguishing PSV columns the read spans
    for &j in &spanned {
        let obs = read.psv_obs[j].expect("spanned column carries an observation");
        let ba = target.alleles.get(j).copied().flatten();
        let ca = competitor.alleles.get(j).copied().flatten();
        if let (Some(ba), Some(ca)) = (ba, ca) {
            if ba != ca {
                let e = match read.psv_qual.get(j).copied().flatten() {
                    Some(q) => super::copy_split::phred_err(q),
                    None => p.error_rate,
                };
                let mut eps_j = (e / 3.0).clamp(0.0, 1.0);
                if p.rna_editing_filter && editing_cols.get(j).copied().unwrap_or(false) {
                    eps_j = eps_j.max(p.edit_rate).clamp(0.0, 1.0);
                }
                eps.push(eps_j);
                if obs == ba {
                    k += 1;
                }
            }
        }
    }
    // distinguishing junctions (from the read's own junction set)
    for &jb in &read.junctions {
        let in_target = boundary_present(jb, &target.junctions, p.boundary_tol);
        let in_comp = boundary_present(jb, &competitor.junctions, p.boundary_tol);
        if in_target != in_comp {
            eps.push(p.junction_err.clamp(0.0, 1.0));
            if in_target {
                k += 1; // read carries a junction `target` has and `competitor` lacks -> supports target
            }
        }
    }
    let p_value = poisson_binomial_upper_tail(k, &eps);
    let min_p = if eps.is_empty() { 1.0 } else { eps.iter().product::<f64>() };
    (p_value, min_p)
}

/// Per-read, per-copy evidence: the SAME log-likelihood + identifiability bound the one-shot gate
/// ([`assign_read_editing`]) computes internally, extracted so the EM copy-assignment engine can reuse
/// the exact same emission model (DRY -- no duplicated likelihood logic between the gate and EM).
#[derive(Clone, Debug)]
pub(crate) struct ReadEvidence {
    /// Per-copy log-likelihood (PSV term + junction term), parallel to `copies`.
    pub logl: Vec<f64>,
    /// Identifiability bound: the best attainable p-value the read could achieve against its hardest
    /// competitor (given `best = argmax(logl)`). `>= alpha` (family-wide-adjusted) ⇒ genuinely `Tied`.
    pub min_p: f64,
    /// Number of features (PSV columns + junction boundaries) the read observes where copies differ.
    pub n_decisive: usize,
}

/// Compute a read's per-copy log-likelihood, decisive-feature count, and identifiability bound against
/// `copies`. This is the single source of the emission model: [`assign_read_editing`] (the one-shot gate)
/// and the EM engine both go through this function, so neither can drift from the other.
///
/// `min_p` depends on `best = argmax(logl)`, so the three quantities are computed in the same order as
/// the gate always has: `logl` + `n_decisive` first, then `best`, then the `copy_pair_significance`
/// competitor loop for `min_p`.
pub(crate) fn read_copy_evidence(
    read: &ReadFeatures,
    copies: &[CopyProfile],
    p: &AssignParams,
    editing_cols: &[bool],
) -> ReadEvidence {
    let n = copies.len();
    let mut logl = vec![0.0f64; n];
    let mut n_decisive = 0usize;

    // --- PSV term (each column weighted by its own base quality when available) ---
    let n_cols = read.psv_obs.len();
    // Columns this read actually spans (psv_obs = Some). Built once and reused by the PSV term and the
    // significance gate below, so neither re-scans the full (often thousands-wide) column set per copy.
    // Reads span a small slice of the family's columns, so this turns the gate's per-competitor O(n_cols)
    // scan into O(spanned) — same observations, same arithmetic, just no walking over absent columns.
    let spanned: Vec<usize> = (0..n_cols).filter(|&j| read.psv_obs[j].is_some()).collect();
    for &j in &spanned {
        let obs = read.psv_obs[j].expect("spanned column carries an observation");
        // per-base error: the column's QV if the read carried one, else the flat default.
        let e = match read.psv_qual.get(j).copied().flatten() {
            Some(q) => super::copy_split::phred_err(q),
            None => p.error_rate,
        };
        let lp_match = (1.0 - e).ln();
        let lp_mis = (e / 3.0).ln();
        // decisive iff the copies disagree at this column (among those with a defined allele)
        let mut seen: Option<u8> = None;
        let mut differ = false;
        for c in copies {
            if let Some(a) = c.alleles.get(j).copied().flatten() {
                match seen {
                    None => seen = Some(a),
                    Some(s) => {
                        if s != a {
                            differ = true;
                        }
                    }
                }
            }
        }
        if differ {
            n_decisive += 1;
        }
        for (ci, c) in copies.iter().enumerate() {
            if let Some(a) = c.alleles.get(j).copied().flatten() {
                logl[ci] += if obs == a { lp_match } else { lp_mis };
            }
        }
    }

    // --- junction (copy-specific splicing) term ---
    for &jb in &read.junctions {
        let present: Vec<bool> = copies
            .iter()
            .map(|c| boundary_present(jb, &c.junctions, p.boundary_tol))
            .collect();
        let np = present.iter().filter(|&&x| x).count();
        if np > 0 && np < n {
            n_decisive += 1; // some copies have this junction, others lack it
        }
        for ci in 0..n {
            logl[ci] += if present[ci] { p.junction_weight } else { -p.junction_weight };
        }
    }

    // argmax (earliest on ties) -- needed before `min_p` since the competitor loop below is best-relative.
    let mut best = 0usize;
    for i in 1..n {
        if logl[i] > logl[best] {
            best = i;
        }
    }

    // --- IsoCon identifiability bound: `best` vs each competitor on their DISTINGUISHING obs ---
    let mut min_p = 0.0f64;
    for c in 0..n {
        if c == best {
            continue;
        }
        let (_pbc, attain) = copy_pair_significance(read, &copies[best], &copies[c], p, editing_cols);
        if attain > min_p {
            min_p = attain;
        }
    }

    ReadEvidence { logl, min_p, n_decisive }
}

/// Assign a read to its most likely copy. Returns `None` only if `copies` is empty.
/// Deterministic: copies are scored in slice order and ties resolve to the earliest (lowest index).
pub fn assign_read(read: &ReadFeatures, copies: &[CopyProfile], p: &AssignParams) -> Option<Assignment> {
    assign_read_editing(read, copies, p, &[])
}

/// Like [`assign_read`], but `editing_cols[j] == true` marks PSV column `j` as an A-to-I RNA-editing site
/// (detected by `copy_assign_pipeline::detect_editing_columns`). When `p.rna_editing_filter` is on, such a
/// column gets its εⱼ inflated to `p.edit_rate` in the significance certificate — an edited base then cannot
/// fake copy-support. The likelihood ranking is unchanged; only the certificate stops trusting the column.
/// `editing_cols` may be shorter than the PSV columns (missing = not an editing site).
pub fn assign_read_editing(
    read: &ReadFeatures,
    copies: &[CopyProfile],
    p: &AssignParams,
    editing_cols: &[bool],
) -> Option<Assignment> {
    let n = copies.len();
    if n == 0 {
        return None;
    }
    let ev = read_copy_evidence(read, copies, p, editing_cols);
    let ReadEvidence { logl, min_p, n_decisive } = ev;

    // argmax (earliest on ties) + runner-up
    let mut best = 0usize;
    for i in 1..n {
        if logl[i] > logl[best] {
            best = i;
        }
    }
    let mut second = f64::NEG_INFINITY;
    for i in 0..n {
        if i != best && logl[i] > second {
            second = logl[i];
        }
    }
    let margin = if n > 1 { logl[best] - second } else { f64::INFINITY };

    // --- IsoCon significance certificate: test `best` vs each competitor on their DISTINGUISHING obs ---
    // p_read = least-significant (max) pairwise p (min_p was already computed by `read_copy_evidence`,
    // relative to the SAME `best`, so it is reused as-is here).
    let mut p_read = 0.0f64;
    for c in 0..n {
        if c == best {
            continue;
        }
        let (pbc, _attain) = copy_pair_significance(read, &copies[best], &copies[c], p, editing_cols);
        if pbc > p_read {
            p_read = pbc;
        }
    }
    let (resolvable, status) = if p.use_margin_gate {
        let resolvable = n_decisive >= 1;
        let status = if !resolvable {
            AssignStatus::Tied
        } else if margin >= p.margin {
            AssignStatus::Assigned
        } else {
            AssignStatus::Ambiguous
        };
        (resolvable, status)
    } else {
        // Bonferroni: `best` is the argmax over n copies, so to bound the FAMILY-WIDE misassignment rate at
        // alpha the per-competitor certificate must clear alpha/(n-1) (union over the n-1 competitors).
        let thr = p.alpha / (n.saturating_sub(1).max(1) as f64);
        let resolvable = min_p < thr;
        let status = if !resolvable {
            AssignStatus::Tied
        } else if p_read < thr && margin > 0.0 {
            AssignStatus::Assigned
        } else {
            AssignStatus::Ambiguous
        };
        (resolvable, status)
    };

    // Per-copy posterior under a UNIFORM prior = softmax(logl), parallel to `copies`. (Likelihood-normalized;
    // an informative prior is applied downstream.) For a Tied read this spreads over the consistent zone.
    let posterior = {
        let m = logl.iter().cloned().fold(f64::NEG_INFINITY, f64::max);
        let mut e: Vec<f64> = logl.iter().map(|&l| (l - m).exp()).collect();
        let z: f64 = e.iter().sum();
        if z > 0.0 {
            for x in &mut e {
                *x /= z;
            }
        }
        e
    };

    Some(Assignment {
        best_copy: copies[best].copy_id,
        log_lr_margin: margin,
        n_decisive,
        resolvable,
        status,
        p_value: p_read,
        min_p_value: min_p,
        discovery_coupled: false,
        posterior,
    })
}

/// Build the PSV-observation vector for a read aligned in the frame of `psv_positions` (genomic, 0-based,
/// in the mapped copy's coordinate frame). Reuses the `copy_split::allele_at` CIGAR bridge. Bases are
/// forward-genome; the caller reverse-complements for a minus-strand copy before comparing (so allele
/// space matches the copies' transcription-strand alleles). Junction boundaries are computed by the
/// integration layer (they need per-copy exon context) and supplied separately.
pub fn read_psv_obs(read: &AlignedRead, psv_positions: &[u64]) -> Vec<Option<u8>> {
    psv_positions.iter().map(|&p| allele_at(read, p)).collect()
}

/// The read's intron donors/acceptors (genomic), straight from its CIGAR N ops -- the raw input the
/// integration layer maps into transcription-spliced boundary offsets.
pub fn read_introns(read: &AlignedRead) -> Vec<(u64, u64)> {
    intron_chain_of(read)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn cp(copy_id: usize, alleles: &[Option<u8>], junctions: &[i64]) -> CopyProfile {
        CopyProfile { copy_id, alleles: alleles.to_vec(), junctions: junctions.to_vec() }
    }
    fn rf(psv: &[Option<u8>], junctions: &[i64]) -> ReadFeatures {
        ReadFeatures { psv_obs: psv.to_vec(), psv_qual: vec![], junctions: junctions.to_vec() }
    }
    const A: u8 = b'A';
    const C: u8 = b'C';
    const G: u8 = b'G';
    const T: u8 = b'T';

    #[test]
    fn read_copy_evidence_matches_assignment_internals() {
        // two copies differ at column 0 (A vs C); read observes A -> favors copy 0.
        let copies = vec![
            CopyProfile { copy_id: 0, alleles: vec![Some(b'A')], junctions: vec![] },
            CopyProfile { copy_id: 1, alleles: vec![Some(b'C')], junctions: vec![] },
        ];
        let read = ReadFeatures { psv_obs: vec![Some(b'A')], psv_qual: vec![], junctions: vec![] };
        let p = AssignParams::for_alpha(1e-3);
        let ev = read_copy_evidence(&read, &copies, &p, &[false]);
        // logl favors copy 0; min_p equals the Assignment's identifiability bound; one decisive column.
        assert!(ev.logl[0] > ev.logl[1]);
        assert_eq!(ev.n_decisive, 1);
        let a = assign_read_editing(&read, &copies, &p, &[false]).unwrap();
        assert!((ev.min_p - a.min_p_value).abs() < 1e-12);
        // posterior consistency: softmax(logl) == Assignment.posterior
        let m = ev.logl.iter().cloned().fold(f64::NEG_INFINITY, f64::max);
        let mut e: Vec<f64> = ev.logl.iter().map(|&l| (l - m).exp()).collect();
        let z: f64 = e.iter().sum(); for x in &mut e { *x /= z; }
        for (i, &pv) in a.posterior.iter().enumerate() { assert!((e[i] - pv).abs() < 1e-9); }
    }

    // Golden fixture: freezes the EXACT current ReadEvidence so the Task-3 graph refactor is provably
    // bit-identical. If this ever changes, the decision changed — which this project forbids.
    fn golden_fixture() -> (ReadFeatures, Vec<CopyProfile>, AssignParams) {
        let cp = |a: Vec<Option<u8>>, j: Vec<i64>| CopyProfile { copy_id: 0, alleles: a, junctions: j };
        let copies = vec![
            cp(vec![Some(b'A'), Some(b'C'), Some(b'G'), Some(b'T')], vec![100]),
            cp(vec![Some(b'A'), Some(b'G'), Some(b'G'), None],       vec![]),
            cp(vec![Some(b'A'), Some(b'C'), Some(b'A'), Some(b'T')], vec![100, 250]),
        ];
        let read = ReadFeatures {
            psv_obs: vec![Some(b'A'), Some(b'C'), None, Some(b'T')],
            psv_qual: vec![Some(30), None, None, Some(20)],
            junctions: vec![100],
        };
        (read, copies, AssignParams::default())
    }

    const EXPECT_LOGL: [f64; 3] = [4.9859446547926165, -11.90875577931572, 4.9859446547926165];
    const EXPECT_MIN_P: f64 = 1.0;
    const EXPECT_N_DECISIVE: usize = 2;

    #[test]
    fn read_copy_evidence_golden_is_stable() {
        let (read, copies, p) = golden_fixture();
        let ev = read_copy_evidence(&read, &copies, &p, &[]);
        // Capture the CURRENT bit-exact values: run once, read the assert-failure message, paste them back,
        // then this test freezes them. (Use {:?} on f64 for the exact decimal that round-trips.)
        // EXPECTED (fill from the first run, then it must never change):
        let want_logl: Vec<f64> = EXPECT_LOGL.to_vec();
        let want_min_p: f64 = EXPECT_MIN_P;
        let want_n_decisive: usize = EXPECT_N_DECISIVE;
        assert_eq!(ev.logl, want_logl, "logl drifted -> decision changed");
        assert_eq!(ev.min_p, want_min_p, "min_p drifted");
        assert_eq!(ev.n_decisive, want_n_decisive, "n_decisive drifted");
    }

    #[test]
    fn two_copies_one_psv_resolves() {
        // Legacy τ-margin gate: one PSV at flat error resolves (margin >> τ=2).
        let copies = [cp(0, &[Some(A)], &[]), cp(1, &[Some(C)], &[])];
        let r = rf(&[Some(A)], &[]);
        let p = AssignParams { use_margin_gate: true, ..AssignParams::default() };
        let a = assign_read(&r, &copies, &p).unwrap();
        assert_eq!(a.best_copy, 0);
        assert!(a.resolvable);
        assert_eq!(a.n_decisive, 1);
        assert_eq!(a.status, AssignStatus::Assigned);
    }

    /// Per-base QV weights the PSV likelihood: the SAME matched allele is far more decisive at
    /// high quality (Q40, e~1e-4) than at low quality (Q3 -> capped e=0.25), i.e. a low-QV base
    /// at a decisive column is trusted less. An empty/None QV falls back to the flat error rate
    /// (legacy behaviour), so existing callers are unchanged.
    #[test]
    fn per_base_quality_weights_psv_likelihood() {
        let copies = [cp(0, &[Some(A)], &[]), cp(1, &[Some(C)], &[])];
        let p = AssignParams::default();
        let mk = |q: u8| ReadFeatures { psv_obs: vec![Some(A)], psv_qual: vec![Some(q)], junctions: vec![] };
        let hi = assign_read(&mk(40), &copies, &p).unwrap();
        let lo = assign_read(&mk(3), &copies, &p).unwrap();
        assert_eq!(hi.best_copy, 0);
        assert_eq!(lo.best_copy, 0);
        assert!(hi.log_lr_margin > lo.log_lr_margin + 5.0,
            "Q40 margin {} should dominate Q3 margin {}", hi.log_lr_margin, lo.log_lr_margin);
        // a missing QV (None / empty psv_qual) == the flat default error rate
        let flat = assign_read(&rf(&[Some(A)], &[]), &copies, &p).unwrap();
        let none_q = assign_read(
            &ReadFeatures { psv_obs: vec![Some(A)], psv_qual: vec![None], junctions: vec![] },
            &copies, &p).unwrap();
        assert!((flat.log_lr_margin - none_q.log_lr_margin).abs() < 1e-9);
    }

    #[test]
    fn no_features_is_tied() {
        let copies = [cp(0, &[Some(A)], &[100]), cp(1, &[Some(C)], &[200])];
        let r = rf(&[None], &[]); // spans no PSV column, no junctions
        let a = assign_read(&r, &copies, &AssignParams::default()).unwrap();
        assert!(!a.resolvable);
        assert_eq!(a.n_decisive, 0);
        assert_eq!(a.status, AssignStatus::Tied);
    }

    #[test]
    fn shared_allele_not_decisive() {
        // both copies carry A at the only spanned column -> not decisive -> tied
        let copies = [cp(0, &[Some(A)], &[]), cp(1, &[Some(A)], &[])];
        let r = rf(&[Some(A)], &[]);
        let a = assign_read(&r, &copies, &AssignParams::default()).unwrap();
        assert_eq!(a.n_decisive, 0);
        assert_eq!(a.status, AssignStatus::Tied);
    }

    #[test]
    fn junction_only_resolves() {
        // A single copy-specific junction resolves under the significance gate now that
        // junction_err=1e-4 < alpha=1e-3 (one clean junction => p_read=1e-4 < alpha, margin >> 0).
        let copies = [cp(0, &[None], &[100]), cp(1, &[None], &[500])];
        let r = rf(&[None], &[100]);
        let a = assign_read(&r, &copies, &AssignParams::default()).unwrap();
        assert_eq!(a.best_copy, 0);
        assert!(a.resolvable);
        assert_eq!(a.n_decisive, 1);
        assert_eq!(a.status, AssignStatus::Assigned);
    }

    #[test]
    fn psv_error_tolerated() {
        // read matches copy0 at 2 of 3 columns (1 sporadic error) -> still copy0, assigned
        let copies = [
            cp(0, &[Some(A), Some(A), Some(A)], &[]),
            cp(1, &[Some(C), Some(C), Some(C)], &[]),
        ];
        let r = rf(&[Some(A), Some(C), Some(A)], &[]); // middle base is an error
        let a = assign_read(&r, &copies, &AssignParams::default()).unwrap();
        assert_eq!(a.best_copy, 0);
        assert_eq!(a.status, AssignStatus::Assigned);
    }

    #[test]
    fn five_copies_two_psv_unique() {
        // sim5x K=2 analogue: base-4 allele vectors over 2 columns give 5 distinct copies.
        // Two Q40 PSVs: each competitor pair has at least one distinguishing column with eps~3.3e-5
        // so min_p <= 3.3e-5 << alpha=1e-3 (resolvable) and p_read <= 3.3e-5 < alpha (Assigned).
        let bases = [A, C, G, T];
        let allele = |c: usize, j: usize| bases[(c / 4usize.pow(j as u32)) % 4];
        let copies: Vec<CopyProfile> = (0..5)
            .map(|c| cp(c, &[Some(allele(c, 0)), Some(allele(c, 1))], &[]))
            .collect();
        // a read carrying copy 3's alleles at Q40 -> uniquely copy 3
        let r = ReadFeatures {
            psv_obs: vec![Some(allele(3, 0)), Some(allele(3, 1))],
            psv_qual: vec![Some(40), Some(40)],
            junctions: vec![],
        };
        let a = assign_read(&r, &copies, &AssignParams::default()).unwrap();
        assert_eq!(a.best_copy, 3);
        assert_eq!(a.status, AssignStatus::Assigned);
    }

    #[test]
    fn five_copies_one_column_ambiguous_for_collision() {
        // K=1 collision: copies 0 and 4 share allele A at column 0 (base-4: 0%4==0, 4%4==0).
        // Legacy τ-margin gate: the column is decisive (copies 1,2,3 differ) but margin=0 < τ=2 -> Ambiguous.
        // Under the significance gate, copy0 and copy4 are pairwise indistinguishable (no distinguishing
        // column between them) -> Tied, which is semantically more correct; use_margin_gate to preserve
        // the original Ambiguous assertion.
        let bases = [A, C, G, T];
        let allele = |c: usize, j: usize| bases[(c / 4usize.pow(j as u32)) % 4];
        let copies: Vec<CopyProfile> = (0..5)
            .map(|c| cp(c, &[Some(allele(c, 0))], &[])) // only column 0
            .collect();
        // a read from copy 4 spanning only column 0 -> allele A matches copy0 AND copy4 -> ambiguous
        let r = rf(&[Some(allele(4, 0))], &[]);
        let p = AssignParams { use_margin_gate: true, ..AssignParams::default() };
        let a = assign_read(&r, &copies, &p).unwrap();
        assert!(a.resolvable); // the column IS decisive (copies 1,2,3 differ)
        assert_eq!(a.status, AssignStatus::Ambiguous); // but copy0/copy4 tie -> margin 0 < 2
    }

    #[test]
    fn junction_augments_decisive_count() {
        // 1 decisive PSV column + 1 decisive junction -> n_decisive == 2
        let copies = [cp(0, &[Some(A)], &[100]), cp(1, &[Some(C)], &[500])];
        let r = rf(&[Some(A)], &[100]);
        let a = assign_read(&r, &copies, &AssignParams::default()).unwrap();
        assert_eq!(a.n_decisive, 2);
        assert_eq!(a.best_copy, 0);
    }

    #[test]
    fn boundary_tolerance_matches_jitter() {
        // Single junction at junction_err=1e-4 < alpha resolves under the default gate; this asserts the boundary tolerance (102 vs 100 within tol=4) is what makes the junction match.
        let copies = [cp(0, &[None], &[100]), cp(1, &[None], &[900])];
        let r = rf(&[None], &[102]);
        let a = assign_read(&r, &copies, &AssignParams::default()).unwrap();
        assert_eq!(a.best_copy, 0);
        assert_eq!(a.status, AssignStatus::Assigned);
    }

    #[test]
    fn copy_matching_ref_allele_inherits_it() {
        // the Python bug guard: a copy that matches the ref allele must carry that base, not None.
        // copy0=(A,A), copy1=(A,C): a read (A,C) must NOT be mis-assigned; only col1 is decisive.
        // Legacy τ-margin gate: one decisive PSV at flat error gives margin >> τ=2.
        let copies = [cp(0, &[Some(A), Some(A)], &[]), cp(1, &[Some(A), Some(C)], &[])];
        let r = rf(&[Some(A), Some(C)], &[]); // matches copy1 at col1 (the decisive one)
        let p = AssignParams { use_margin_gate: true, ..AssignParams::default() };
        let a = assign_read(&r, &copies, &p).unwrap();
        assert_eq!(a.best_copy, 1);
        assert_eq!(a.status, AssignStatus::Assigned);
    }

    #[test]
    fn deterministic_tiebreak_to_lowest_index() {
        // identical profiles -> not decisive -> tied, and best resolves to the earliest copy.
        let copies = [cp(7, &[Some(A)], &[]), cp(3, &[Some(A)], &[])];
        let r = rf(&[Some(A)], &[]);
        let a = assign_read(&r, &copies, &AssignParams::default()).unwrap();
        assert_eq!(a.best_copy, 7); // first in slice order
        assert_eq!(a.status, AssignStatus::Tied);
    }

    #[test]
    fn empty_copies_returns_none() {
        let r = rf(&[Some(A)], &[]);
        assert!(assign_read(&r, &[], &AssignParams::default()).is_none());
    }

    #[test]
    fn read_psv_obs_reuses_cigar_bridge() {
        // 5M2N5M over ref_start 10: positions 10..15 and 17..22 matched; 15,16 inside the intron.
        let read = AlignedRead { ref_start: 10, cigar: vec![('M', 5), ('N', 2), ('M', 5)], seq: b"AAAAACCCCC".to_vec(), qual: vec![] };
        let obs = read_psv_obs(&read, &[12, 16, 18]);
        assert_eq!(obs, vec![Some(b'A'), None, Some(b'C')]);
    }

    #[test]
    fn poisson_binomial_matches_binomial_for_equal_probs() {
        // n=4, p=0.5: P(X>=2) = (6+4+1)/16 = 11/16.
        let probs = [0.5_f64; 4];
        let got = poisson_binomial_upper_tail(2, &probs);
        assert!((got - 11.0 / 16.0).abs() < 1e-12, "got {got}");
    }

    #[test]
    fn poisson_binomial_edge_cases() {
        assert_eq!(poisson_binomial_upper_tail(0, &[]), 1.0); // P(>=0) = 1
        assert_eq!(poisson_binomial_upper_tail(1, &[]), 0.0); // 0 trials cannot reach 1
        assert_eq!(poisson_binomial_upper_tail(1, &[0.0, 0.0]), 0.0); // no possible success
        assert!((poisson_binomial_upper_tail(1, &[1.0]) - 1.0).abs() < 1e-12); // certain success
        // k == n: tail equals the product of probs (all must succeed)
        let probs = [0.1_f64, 0.2, 0.05];
        let prod: f64 = probs.iter().product();
        assert!((poisson_binomial_upper_tail(3, &probs) - prod).abs() < 1e-12);
    }

    #[test]
    fn poisson_binomial_monotone_in_k() {
        let probs = [0.3_f64, 0.4, 0.2, 0.6];
        let mut prev = poisson_binomial_upper_tail(0, &probs);
        for k in 1..=probs.len() {
            let cur = poisson_binomial_upper_tail(k, &probs);
            assert!(cur <= prev + 1e-12, "P(>={k}) must not exceed P(>={})", k - 1);
            prev = cur;
        }
    }

    #[test]
    fn assign_params_alpha_defaults_and_constructor() {
        let d = AssignParams::default();
        assert_eq!(d.alpha, 1e-3);
        assert_eq!(d.junction_err, 1e-4);
        assert!(!d.use_margin_gate);
        let a = AssignParams::for_alpha(0.05);
        assert_eq!(a.alpha, 0.05);
        assert_eq!(a.error_rate, d.error_rate); // other fields inherit the default
    }

    fn two_copies_2psv() -> Vec<CopyProfile> {
        vec![
            CopyProfile { copy_id: 0, alleles: vec![Some(b'A'), Some(b'C')], junctions: vec![] },
            CopyProfile { copy_id: 1, alleles: vec![Some(b'G'), Some(b'T')], junctions: vec![] },
        ]
    }

    #[test]
    fn sig_two_highq_psv_supporting_best_is_assigned() {
        let r = ReadFeatures { psv_obs: vec![Some(b'A'), Some(b'C')], psv_qual: vec![Some(40), Some(40)], junctions: vec![] };
        let a = assign_read(&r, &two_copies_2psv(), &AssignParams::default()).unwrap();
        assert_eq!(a.best_copy, 0);
        assert_eq!(a.status, AssignStatus::Assigned);
        assert!(a.p_value < 1e-3, "p_value {} should clear alpha", a.p_value);
    }

    #[test]
    fn sig_no_distinguishing_columns_is_tied() {
        let r = ReadFeatures { psv_obs: vec![None, None], psv_qual: vec![None, None], junctions: vec![] };
        let a = assign_read(&r, &two_copies_2psv(), &AssignParams::default()).unwrap();
        assert_eq!(a.status, AssignStatus::Tied);
        assert!(a.min_p_value >= 1e-3);
    }

    #[test]
    fn sig_conflicting_support_is_ambiguous() {
        let r = ReadFeatures { psv_obs: vec![Some(b'A'), Some(b'T')], psv_qual: vec![Some(40), Some(40)], junctions: vec![] };
        let a = assign_read(&r, &two_copies_2psv(), &AssignParams::default()).unwrap();
        assert_eq!(a.status, AssignStatus::Ambiguous);
        assert!(a.min_p_value < 1e-3, "two Q40 columns CAN resolve in principle");
        assert!(a.p_value < 1e-3, "best-support IS significant, but the LLR margin is 0 (balanced) -> abstain");
    }

    #[test]
    fn sig_legacy_margin_gate_still_selectable() {
        let copies = vec![
            CopyProfile { copy_id: 0, alleles: vec![Some(b'A')], junctions: vec![] },
            CopyProfile { copy_id: 1, alleles: vec![Some(b'G')], junctions: vec![] },
        ];
        let r = ReadFeatures { psv_obs: vec![Some(b'A')], psv_qual: vec![Some(40)], junctions: vec![] };
        let p = AssignParams { use_margin_gate: true, ..AssignParams::default() };
        let a = assign_read(&r, &copies, &p).unwrap();
        assert_eq!(a.status, AssignStatus::Assigned);
    }

    #[test]
    fn editing_params_defaults() {
        let d = AssignParams::default();
        assert!(d.rna_editing_filter);
        assert_eq!(d.edit_rate, 0.2);
    }

    #[test]
    fn sig_editing_col_downweights_to_abstain() {
        // two copies differing at ONE A<->G column. A high-Q read resolves it normally; flagging the column
        // as an editing site inflates εⱼ to edit_rate so the certificate abstains (the edited base can't be
        // trusted as copy-distinguishing). Disabling the filter restores the assignment.
        let copies = vec![
            CopyProfile { copy_id: 0, alleles: vec![Some(b'G')], junctions: vec![] },
            CopyProfile { copy_id: 1, alleles: vec![Some(b'A')], junctions: vec![] },
        ];
        let r = ReadFeatures { psv_obs: vec![Some(b'G')], psv_qual: vec![Some(40)], junctions: vec![] };
        let a = assign_read_editing(&r, &copies, &AssignParams::default(), &[]).unwrap();
        assert_eq!(a.status, AssignStatus::Assigned, "unflagged single PSV resolves");
        let a2 = assign_read_editing(&r, &copies, &AssignParams::default(), &[true]).unwrap();
        assert_eq!(a2.status, AssignStatus::Tied, "flagged editing column can't resolve -> Tied");
        assert!(a2.min_p_value >= 1e-3);
        let p_off = AssignParams { rna_editing_filter: false, ..AssignParams::default() };
        let a3 = assign_read_editing(&r, &copies, &p_off, &[true]).unwrap();
        assert_eq!(a3.status, AssignStatus::Assigned, "filter off ignores the flag");
        // assign_read wrapper == assign_read_editing with empty flags
        assert_eq!(assign_read(&r, &copies, &AssignParams::default()).unwrap().status, AssignStatus::Assigned);
    }

    #[test]
    fn tau_calibration_operating_points() {
        // the two operating points seen in the codebase are just two choices of p:
        assert!((tau_from_p(1e-3) - 6.9068).abs() < 1e-3, "p=1e-3 => tau~6.9 (Eichler AS>=10 analog)");
        assert!((tau_from_p(0.1192029) - 2.0).abs() < 1e-3, "tau=2.0 => p~0.119 (recall mode)");
        // round-trip
        for &p in &[1e-4, 1e-3, 0.05, 0.119, 0.3] {
            assert!((p_from_tau(tau_from_p(p)) - p).abs() < 1e-9);
        }
        // the constructor wires margin = tau(p)
        assert!((AssignParams::for_target_misassignment(1e-3).margin - 6.9068).abs() < 1e-3);
        // monotone: stricter p => larger margin
        assert!(tau_from_p(1e-4) > tau_from_p(1e-2));
    }

    #[test]
    fn sig_gate_bonferroni_controls_familywide_error() {
        // 11 copies; each competitor c=1..10 differs from the true copy (copy0) at ONE column -> a low-power,
        // many-competitor regime (the K>=3 / 1-PSV tail). Without the Bonferroni correction the argmax winner's
        // curse inflates realized misassignment to ~K*(e/3) > alpha; with it (threshold alpha/(n-1)) the gate
        // abstains on the marginal reads and keeps realized misassignment <= alpha.
        let n_comp = 10usize;
        let m = n_comp; // one distinguishing column per competitor
        let mk = |c: usize| -> Vec<Option<u8>> {
            (0..m).map(|j| if c >= 1 && j == c - 1 { Some(b'C') } else { Some(b'A') }).collect()
        };
        let copies: Vec<CopyProfile> =
            (0..=n_comp).map(|c| CopyProfile { copy_id: c, alleles: mk(c), junctions: vec![] }).collect();
        let bases = [b'A', b'C', b'G', b'T'];
        let q = 30u8;
        let e = 10f64.powf(-(q as f64) / 10.0);
        let run = |alpha: f64| -> (usize, usize) {
            let mut state = 0x9E3779B97F4A7C15u64;
            let mut next = || {
                state ^= state << 13;
                state ^= state >> 7;
                state ^= state << 17;
                state
            };
            let p = AssignParams::for_alpha(alpha);
            let (mut assigned, mut wrong) = (0usize, 0usize);
            for _ in 0..50_000 {
                let psv_obs: Vec<Option<u8>> = (0..m)
                    .map(|_| {
                        if (next() as f64 / u64::MAX as f64) < e {
                            Some(bases[(next() % 4) as usize])
                        } else {
                            Some(b'A')
                        }
                    })
                    .collect();
                let r = ReadFeatures { psv_obs, psv_qual: vec![Some(q); m], junctions: vec![] };
                let a = assign_read(&r, &copies, &p).unwrap();
                if a.status == AssignStatus::Assigned {
                    assigned += 1;
                    if a.best_copy != 0 {
                        wrong += 1; // all reads are truly copy0
                    }
                }
            }
            (assigned, wrong)
        };
        // Corrected default alpha=1e-3 -> internal threshold 1e-3/10 = 1e-4: the marginal reads abstain, error controlled.
        let (a_corr, w_corr) = run(1e-3);
        let rate_corr = if a_corr > 0 { w_corr as f64 / a_corr as f64 } else { 0.0 };
        assert!(rate_corr <= 1e-3, "Bonferroni gate realized misassignment {rate_corr} must respect alpha=1e-3");
        // alpha=1e-2 -> internal threshold 1e-2/10 = 1e-3 (= the UN-corrected per-pair level): reproduces the
        // winner's-curse inflation the correction removes.
        let (a_unc, w_unc) = run(1e-2);
        let rate_unc = if a_unc > 0 { w_unc as f64 / a_unc as f64 } else { 0.0 };
        assert!(a_unc > 1000, "under-corrected level should assign the marginal reads ({a_unc})");
        assert!(rate_unc > 1e-3, "without the correction the union-bound inflation exceeds alpha ({rate_unc})");
    }

    #[test]
    fn sig_gate_is_calibrated_realized_error_tracks_alpha() {
        // Two copies differing at 6 PSV columns. Simulate reads from a known true copy, inject errors at a
        // realistic HiFi rate, run the gate, and check the realized misassignment among ASSIGNED reads.
        let m = 6usize;
        let copy0: Vec<Option<u8>> = (0..m).map(|_| Some(b'A')).collect();
        let copy1: Vec<Option<u8>> = (0..m).map(|_| Some(b'C')).collect();
        let copies = vec![
            CopyProfile { copy_id: 0, alleles: copy0.clone(), junctions: vec![] },
            CopyProfile { copy_id: 1, alleles: copy1.clone(), junctions: vec![] },
        ];
        // deterministic xorshift RNG
        let mut state = 0x2545F4914F6CDD1Du64;
        let mut next = || {
            state ^= state << 13;
            state ^= state >> 7;
            state ^= state << 17;
            state
        };
        let bases = [b'A', b'C', b'G', b'T'];
        let q = 30u8; // e = 1e-3
        let e = 10f64.powf(-(q as f64) / 10.0); // inline (phred_err's path differs inside mod tests)
        let run = |alpha: f64, next: &mut dyn FnMut() -> u64| -> (usize, usize) {
            let p = AssignParams::for_alpha(alpha);
            let (mut assigned, mut wrong) = (0usize, 0usize);
            for _ in 0..20_000 {
                let truth = (next() % 2) as usize;
                let template = if truth == 0 { &copy0 } else { &copy1 };
                let psv_obs: Vec<Option<u8>> = template
                    .iter()
                    .map(|t| {
                        let tb = t.unwrap();
                        if (next() as f64 / u64::MAX as f64) < e {
                            Some(bases[(next() % 4) as usize])
                        } else {
                            Some(tb)
                        }
                    })
                    .collect();
                let r = ReadFeatures { psv_obs, psv_qual: vec![Some(q); m], junctions: vec![] };
                let a = assign_read(&r, &copies, &p).unwrap();
                if a.status == AssignStatus::Assigned {
                    assigned += 1;
                    if a.best_copy != truth {
                        wrong += 1;
                    }
                }
            }
            (assigned, wrong)
        };
        let (a_hi, w_hi) = run(1e-2, &mut next);
        let (a_lo, w_lo) = run(1e-4, &mut next);
        assert!(a_hi > 1000 && a_lo > 1000, "should assign a substantial fraction ({a_hi}, {a_lo})");
        let rate_hi = w_hi as f64 / a_hi as f64;
        let rate_lo = w_lo as f64 / a_lo as f64;
        assert!(rate_hi <= 3e-2, "alpha=1e-2 realized {rate_hi}");
        assert!(rate_lo <= 3e-4, "alpha=1e-4 realized {rate_lo}");
        assert!(rate_lo <= rate_hi + 1e-9, "stricter alpha must not increase error ({rate_lo} vs {rate_hi})");
    }

    #[test]
    fn bubble_graph_from_copies_structure() {
        // 3 copies over 3 columns: col0 all 'A' (not decisive); col1 A/C/A (decisive); col2 copy2 gap
        let mk = |a: Vec<Option<u8>>| CopyProfile { copy_id: 0, alleles: a, junctions: vec![] };
        let copies = vec![
            mk(vec![Some(b'A'), Some(b'A'), Some(b'G')]),
            mk(vec![Some(b'A'), Some(b'C'), Some(b'G')]),
            mk(vec![Some(b'A'), Some(b'A'), None]),
        ];
        let g = BubbleGraph::from_copies(&copies);
        assert_eq!(g.n_copies, 3);
        assert_eq!(g.bubbles.len(), 3);
        assert_eq!(g.bubbles[0].col, 0);
        assert!(!g.bubbles[0].decisive);                 // all 'A'
        assert!(g.bubbles[1].decisive);                  // A vs C
        assert!(!g.bubbles[2].decisive);                 // G, G, gap -> one distinct allele
        assert_eq!(g.bubbles[1].copy_allele, vec![Some(b'A'), Some(b'C'), Some(b'A')]);
        assert_eq!(g.bubbles[2].copy_allele, vec![Some(b'G'), Some(b'G'), None]);
    }
}
