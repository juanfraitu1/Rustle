//! EM copy-assignment core: maximum-likelihood soft relaxation of SDA's PSV
//! correlation-clustering (Vollger 2019) on the PSV-aware variation graph.
//!
//! Where SDA hard-clusters PSV-supporting reads into copy groups via
//! correlation clustering, this module treats copy assignment as a mixture
//! model over copy-paths through the family VG and finds a local-max
//! likelihood assignment by Expectation-Maximization:
//! - `e_step` = soft read -> copy-path assignment (posterior responsibility
//!   of each copy-path for each read, given per-read per-copy log-likelihoods
//!   and the current copy-path abundance prior).
//! - `m_step` = copy-path abundance re-estimation (mean responsibility =
//!   path-usage fraction across the VG's copy-paths).
//! - `loglik` = observed-data log-likelihood of the current (logl, pi)
//!   mixture; its monotone non-decrease across an E-step/M-step sweep is the
//!   EM convergence certificate.
//!
//! These three functions are pure numerics on plain slices (`Vec<Vec<f64>>`
//! log-likelihoods, `Vec<f64>` abundances) with no VG/SDA types at this
//! layer -- that grounding is documented here, not coded; the coupling to
//! Task 1's `ReadEvidence.logl` arrives in Task 3.

/// log-sum-exp over a slice, ignoring `-inf` terms so no NaN is produced when
/// some copies have zero prior mass (`ln 0 = -inf`).
fn logsumexp(xs: &[f64]) -> f64 {
    let m = xs
        .iter()
        .cloned()
        .fold(f64::NEG_INFINITY, f64::max);
    if m == f64::NEG_INFINITY {
        return f64::NEG_INFINITY;
    }
    let sum: f64 = xs
        .iter()
        .filter(|&&x| x != f64::NEG_INFINITY)
        .map(|&x| (x - m).exp())
        .sum();
    m + sum.ln()
}

/// `ln(p)`, treating `p == 0` as `-inf` (guards `ln(0)` -> NaN-free).
fn ln_pi(p: f64) -> f64 {
    if p > 0.0 { p.ln() } else { f64::NEG_INFINITY }
}

/// E-step: soft read -> copy-path assignment.
///
/// `gamma[r][k] = softmax_k(logl[r][k] + ln(pi[k]))`, i.e. the posterior
/// responsibility of copy `k` for read `r`. Each output row sums to 1.
pub fn e_step(logl: &[Vec<f64>], pi: &[f64]) -> Vec<Vec<f64>> {
    let ln_pi_k: Vec<f64> = pi.iter().map(|&p| ln_pi(p)).collect();
    logl.iter()
        .map(|row| {
            let joint: Vec<f64> = row
                .iter()
                .zip(ln_pi_k.iter())
                .map(|(&l, &lp)| l + lp)
                .collect();
            let denom = logsumexp(&joint);
            joint
                .iter()
                .map(|&j| {
                    if denom == f64::NEG_INFINITY {
                        0.0
                    } else {
                        (j - denom).exp()
                    }
                })
                .collect()
        })
        .collect()
}

/// M-step: copy-path abundance re-estimation. `pi[k] = mean_r(gamma[r][k])`.
pub fn m_step(gamma: &[Vec<f64>]) -> Vec<f64> {
    let n = gamma.len();
    if n == 0 {
        return Vec::new();
    }
    let k = gamma[0].len();
    let mut sums = vec![0.0_f64; k];
    for row in gamma {
        for (s, &g) in sums.iter_mut().zip(row.iter()) {
            *s += g;
        }
    }
    sums.into_iter().map(|s| s / n as f64).collect()
}

/// Observed-data log-likelihood of the mixture: `sum_r logsumexp_k(logl[r][k]
/// + ln(pi[k]))`. Its non-decrease across an E-step/M-step sweep is the EM
/// convergence certificate.
pub fn loglik(logl: &[Vec<f64>], pi: &[f64]) -> f64 {
    let ln_pi_k: Vec<f64> = pi.iter().map(|&p| ln_pi(p)).collect();
    logl.iter()
        .map(|row| {
            let joint: Vec<f64> = row
                .iter()
                .zip(ln_pi_k.iter())
                .map(|(&l, &lp)| l + lp)
                .collect();
            logsumexp(&joint)
        })
        .sum()
}

/// K-frontier identifiability label for a read, given its EM posterior context.
///
/// This is the SAME test the one-shot gate uses to decide `Tied` vs a confident
/// call ([`ReadEvidence::min_p`] against the family-wide Bonferroni threshold
/// `alpha/(k-1)`), just attached to the EM soft-relaxation's output instead of
/// the gate's hard decision: `Certified` reads sit outside the K-frontier
/// (their best-vs-runner-up p-value beats the correction), `SoftZone` reads
/// are inside it (no amount of EM iteration can identify them -- the mixture
/// legitimately keeps their mass spread across copies).
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum EmLabel {
    /// `min_p` beats the family-wide Bonferroni bound: the read's copy call is identifiable.
    Certified,
    /// `min_p` does not beat the bound: the read is in the K-frontier soft zone (genuinely ambiguous).
    SoftZone,
}

/// K-frontier test: `Certified` iff `min_p < alpha/((k-1).max(1))` (family-wide Bonferroni over the
/// `k-1` pairwise comparisons against the best copy), else `SoftZone`. Identical rule to the one-shot
/// significance gate ([`super::significance_gate`] et al.) -- no separate threshold is introduced here.
pub fn label_read(min_p: f64, alpha: f64, k: usize) -> EmLabel {
    let denom = (k.saturating_sub(1)).max(1) as f64;
    if min_p < alpha / denom {
        EmLabel::Certified
    } else {
        EmLabel::SoftZone
    }
}

/// Output of the EM soft-relaxation driver: final soft assignments, recovered copy abundances, a
/// per-read K-frontier identifiability label, and the convergence trace.
pub struct EmResult {
    /// `posteriors[r][k]` = final E-step responsibility of copy `k` for read `r` (rows sum to 1).
    pub posteriors: Vec<Vec<f64>>,
    /// Final `pi`: recovered per-copy abundance (mean responsibility across all reads).
    pub abundances: Vec<f64>,
    /// Per-read K-frontier identifiability label ([`label_read`] applied to each read's `min_p`).
    pub labels: Vec<EmLabel>,
    /// Number of E/M sweeps actually run before convergence or `max_iter`.
    pub n_iter: usize,
    /// Observed-data log-likelihood after each completed E/M sweep; monotone non-decreasing (the EM
    /// convergence certificate).
    pub loglik_trace: Vec<f64>,
}

/// EM driver: maximum-likelihood soft relaxation of SDA's PSV correlation-clustering, run to
/// convergence on a family's [`ReadEvidence`], then labeled per-read via the K-frontier test
/// ([`label_read`]).
///
/// Starts from a uniform `pi` (abundance-from-certified-reads is a documented follow-up, not
/// implemented here -- no extra threshold is smuggled in to warm-start it) and alternates
/// [`e_step`]/[`m_step`], tracking [`loglik`] after every sweep. Stops when the log-likelihood
/// improvement falls below `eps * |loglik|` (relative convergence, the only numeric tolerance this
/// driver introduces) or after `max_iter` sweeps, whichever comes first.
pub fn em_assign(
    evidence: &[super::copy_assign::ReadEvidence],
    k: usize,
    alpha: f64,
    eps: f64,
    max_iter: usize,
) -> EmResult {
    let logl: Vec<Vec<f64>> = evidence.iter().map(|ev| ev.logl.clone()).collect();
    let mut pi = vec![1.0 / k.max(1) as f64; k];

    let mut loglik_trace = Vec::new();
    let mut gamma = e_step(&logl, &pi);
    let mut n_iter = 0usize;
    let mut prev_ll = loglik(&logl, &pi);
    loglik_trace.push(prev_ll);

    while n_iter < max_iter {
        pi = m_step(&gamma);
        gamma = e_step(&logl, &pi);
        let ll = loglik(&logl, &pi);
        loglik_trace.push(ll);
        n_iter += 1;
        let converged = (ll - prev_ll).abs() < eps * prev_ll.abs();
        prev_ll = ll;
        if converged {
            break;
        }
    }

    let labels = evidence
        .iter()
        .map(|ev| label_read(ev.min_p, alpha, k))
        .collect();

    EmResult {
        posteriors: gamma,
        abundances: pi,
        labels,
        n_iter,
        loglik_trace,
    }
}

#[cfg(test)]
mod em_driver_tests {
    use super::{em_assign, label_read, EmLabel};

    #[test]
    fn em_recovers_planted_abundance_3copy() {
        // 3 copies, one decisive column each (one-hot logl), planted 5:3:2 reads -> pi ~ 0.5,0.3,0.2
        let ev = |best: usize| super::super::copy_assign::ReadEvidence {
            logl: (0..3).map(|c| if c==best {0.0} else {-12.0}).collect(), min_p: 1e-6, n_decisive: 1 };
        let mut reads = vec![]; for _ in 0..5 {reads.push(ev(0));} for _ in 0..3 {reads.push(ev(1));} for _ in 0..2 {reads.push(ev(2));}
        let r = em_assign(&reads, 3, 1e-3, 1e-9, 200);
        assert!((r.abundances[0]-0.5).abs()<0.02 && (r.abundances[1]-0.3).abs()<0.02 && (r.abundances[2]-0.2).abs()<0.02);
        for w in r.loglik_trace.windows(2) { assert!(w[1] >= w[0]-1e-9); }   // monotone
        assert!(r.labels.iter().all(|l| matches!(l, EmLabel::Certified)));   // all identifiable
    }
    #[test]
    fn em_k0_reads_stay_soft_and_abundance_proportional() {
        // no distinguishing feature: logl equal across copies, min_p >= alpha -> SoftZone, posterior == pi.
        let flat = super::super::copy_assign::ReadEvidence { logl: vec![0.0, 0.0], min_p: 1.0, n_decisive: 0 };
        let reads = vec![flat.clone(), flat.clone(), flat];
        let r = em_assign(&reads, 2, 1e-3, 1e-9, 50);
        assert!(r.labels.iter().all(|l| matches!(l, EmLabel::SoftZone)));
        for row in &r.posteriors { assert!((row[0]-r.abundances[0]).abs() < 1e-6); } // soft == prior, never a hard 1/k call
    }
    #[test]
    fn label_read_uses_bonferroni_min_p() {
        assert!(matches!(label_read(1e-6, 1e-3, 3), EmLabel::Certified));  // 1e-6 < 1e-3/2
        assert!(matches!(label_read(0.4,  1e-3, 3), EmLabel::SoftZone));
    }
}

#[cfg(test)]
mod tests {
    use super::{e_step, loglik, m_step};

    #[test]
    fn e_step_rows_are_normalized_posteriors() {
        let logl = vec![vec![0.0_f64, -10.0], vec![-10.0, 0.0]];
        let pi = vec![0.5, 0.5];
        let g = e_step(&logl, &pi);
        for row in &g {
            assert!((row.iter().sum::<f64>() - 1.0).abs() < 1e-12);
        }
        assert!(g[0][0] > 0.99 && g[1][1] > 0.99); // each read favors its copy
    }

    #[test]
    fn e_step_applies_the_abundance_prior() {
        // read equidistant in likelihood; prior 0.9/0.1 must tilt the posterior to 0.9/0.1.
        let g = e_step(&vec![vec![0.0, 0.0]], &vec![0.9, 0.1]);
        assert!((g[0][0] - 0.9).abs() < 1e-9 && (g[0][1] - 0.1).abs() < 1e-9);
    }

    #[test]
    fn m_step_is_mean_responsibility() {
        let g = vec![vec![1.0, 0.0], vec![0.0, 1.0], vec![1.0, 0.0]];
        let pi = m_step(&g);
        assert!((pi[0] - 2.0 / 3.0).abs() < 1e-12 && (pi[1] - 1.0 / 3.0).abs() < 1e-12);
    }

    #[test]
    fn loglik_is_nondecreasing_under_em() {
        // 3 reads, 2 copies, well-separated; one full E/M sweep must not decrease l.
        let logl = vec![vec![0.0, -8.0], vec![0.0, -8.0], vec![-8.0, 0.0]];
        let mut pi = vec![0.5, 0.5];
        let l0 = loglik(&logl, &pi);
        let g = e_step(&logl, &pi);
        pi = m_step(&g);
        let l1 = loglik(&logl, &pi);
        assert!(l1 >= l0 - 1e-12, "EM decreased loglik: {l0} -> {l1}");
    }
}
