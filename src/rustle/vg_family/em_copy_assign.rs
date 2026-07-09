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
