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
//!
//! **STATUS:** OPT-IN — `--em` (src/bin/copy_assign.rs:310-311, `#[arg(long, default_value_t = false)]`) OR `--vg-realign` (src/bin/copy_assign.rs:340-341, default false)  (docs/MODULE_STATUS.md; assigned by reachability, not by this header)

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
/// call (`ReadEvidence::min_p` against the family-wide Bonferroni threshold
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
/// significance gate (`copy_assign::assign_read_editing` et al.) -- no separate threshold is introduced here.
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
/// convergence on a family's `ReadEvidence`, then labeled per-read via the K-frontier test
/// ([`label_read`]).
///
/// Starts from a uniform `pi` (abundance-from-certified-reads is a documented follow-up, not
/// implemented here -- no extra threshold is smuggled in to warm-start it) and alternates
/// [`e_step`]/[`m_step`], tracking [`loglik`] after every sweep. Stops when the log-likelihood
/// improvement falls below `eps * (1 + |loglik|)` (absolute+relative convergence, so a fixed point at
/// exactly `loglik == 0.0` -- e.g. a K=0 flat-evidence family -- still converges immediately instead of
/// running to `max_iter`; the only numeric tolerance this driver introduces) or after `max_iter` sweeps,
/// whichever comes first.
///
/// ⚠ DEFECT ON RECORD (measured 2026-08-11, NOT fixed). THE TOLERANCE IS SET BY A QUANTITY THAT
/// CANNOT AFFECT THE ESTIMATE. `loglik` sums over records and carries the junction term, and when
/// the family's copies share their splice structure that term adds the SAME constant to every copy
/// of every read -- so it cancels in the E-step and cannot move `gamma` or `pi`, but it does not
/// cancel inside `|loglik|`, where it inflates the stopping tolerance. The driver therefore stops
/// EARLY relative to the same rule applied to the likelihood that actually drives the estimate:
/// measured 1 sweep against 3, 7 against 21, 7 against 24, 2 against 12 -- 3.0x to 6.0x -- on the
/// four sim5x families where it can be checked. (The check is exact: a PSV-only reconstruction that
/// omits the junction term entirely reproduces the binary's `gamma` on 26,774/26,774 records and
/// its `pi` to 5e-5, which is what proves the term is a per-read constant there.)
/// NOT FIXED HERE because it is not inert: `n_iter`, `posteriors` and `abundances` all move, and
/// `abundances` reaches `fa.copy_abundance` and therefore default output via
/// `denovo_pipeline::recompute_realign_abundance`. Running to true convergence changes NO decision
/// on the balanced sim5x panel (0/3081 co-committed reads move), but that panel's true `pi` is
/// uniform, so it cannot exercise the failure mode; the one skewed family available
/// (`unb_e003t7`) has copies of five different spans, so its junction term is NOT constant and the
/// reconstruction cannot score it. THE STOPPING RULE IS UNMEASURED ON SKEWED ABUNDANCE.
/// Evidence: `close_o1o2/EM_DEFECTS.txt` sections 1, 2, 4.
pub(crate) fn em_assign(
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
        let converged = (ll - prev_ll).abs() < eps * (1.0 + prev_ll.abs());
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

/// Public entry point for running the EM engine on a family straight off `denovo_pipeline`'s
/// already-extracted PSV data (`FamilyAssignment::read_psv_obs` / `copy_psv_alleles`), for callers
/// outside this crate (the `copy_assign` binary's `--em` mode). `em_assign`/`read_copy_evidence`/
/// `ReadEvidence` are `pub(crate)` and so unreachable from a binary crate; this wrapper takes only
/// plain/public types and goes through [`super::copy_assign::read_copy_evidence`] with the SAME
/// per-read/per-copy likelihood the one-shot significance gate uses (Task 1), now including the
/// Clair3-RNA-style A->I editing-column filter (`copy_assign_pipeline::detect_editing_columns`) --
/// so an editing column is downweighted here exactly as it is in the hard gate.
///
/// Task H2 (VG-harmony): threads O3's copy-specific JUNCTIONS through this wrapper alongside the
/// PSV alleles, so a copy distinguished only by a novel splice (not present in the reference, and
/// therefore not a PSV column) is resolvable by the EM exactly as it already is by the one-shot
/// gate. `read_junctions[i]` / `copy_junctions[k]` are parallel to `read_obs[i]` / `copy_alleles[k]`;
/// passing empty slices (`&[]`) for both reproduces the PREVIOUS PSV-only behavior byte-for-byte
/// (`junctions.get(..).cloned().unwrap_or_default()` yields `vec![]` for every read/copy, which is
/// exactly what this wrapper always built before). No per-base quality is threaded (`psv_qual` stays
/// empty), mirroring `copy_assign_pipeline::soft_quantify_em`'s PSV-only abundance model -- only the
/// junction gap this task closes is filled. Consequently a read whose hard-gate call depends on
/// per-base-quality evidence can still get a *different* per-read label here than `.assignments.tsv`
/// reports -- the abundance estimate is a light-weight complement to the hard PSV+junction
/// assignment, not a re-derivation of it, and is not claimed to reproduce it read-for-read.
///
/// ⚠ DEFECT ON RECORD (measured 2026-08-11, NOT fixed -- the fix is not local, see below).
/// The empty `psv_qual` is not a cosmetic simplification and its blast radius is NOT limited to
/// `--em`:
///   * `read_copy_evidence` falls through to the flat `p.error_rate` (default 3e-3) for EVERY
///     column, where the hard gate uses `phred_err(q)` clamped to 1e-4. On Q40+ reads that is a
///     30x inflation of the per-column error, hence a 30^m inflation of `min_p` over `m`
///     distinguishing columns, hence of the `label_read` Bonferroni test. On the sim5x panel every
///     copy pair differs at exactly ONE column, so `min_p` = e/3 = 1.0e-3 against a bound of
///     alpha/(K-1) = 2.5e-4: the EM certifies **0 of 26,774 read-records**, where the same reads
///     at 1e-4 certify 3,746. The gate's own `.assignments.tsv` prints `min_p_value = 3.333e-5`
///     for those reads, so the two halves of one likelihood model disagree 30-fold on the error
///     rate they assume.
///   * It also reaches DEFAULT output. `denovo_pipeline::recompute_realign_abundance` calls this
///     same wrapper and overwrites `fa.copy_abundance`, which `copy_assign.rs` writes out; so a
///     run with `--vg-realign` and no `--em` is already consuming this likelihood.
/// What it does NOT do is change any assignment: `e` enters every copy through the same two
/// constants `ln(1-e)` and `ln(e/3)`, so it stretches the gaps between copies without reordering
/// them. Measured: the EM argmax agrees with the gate on 3081/3081 co-committed reads at BOTH
/// error rates. It changes what is CERTIFIED, never what is CHOSEN.
/// NOT FIXED HERE because the fix is not one line: `FamilyAssignment` carries no `read_psv_qual`
/// parallel to `read_psv_obs`, two of the four sites that write `read_psv_obs` (the vg-realign
/// correction and novel-pool admission paths) have no quality to write, and threading it would
/// move `copy_abundance` on default runs. Evidence: `close_o1o2/EM_DEFECTS.txt` sections 1, 3, 4.
pub fn em_assign_family(
    read_obs: &[Vec<Option<u8>>],
    copy_alleles: &[Vec<Option<u8>>],
    read_junctions: &[Vec<i64>],
    copy_junctions: &[Vec<i64>],
    params: &super::copy_assign::AssignParams,
    eps: f64,
    max_iter: usize,
) -> EmResult {
    let copies: Vec<super::copy_assign::CopyProfile> = copy_alleles
        .iter()
        .enumerate()
        .map(|(k, alleles)| super::copy_assign::CopyProfile {
            copy_id: k,
            alleles: alleles.clone(),
            junctions: copy_junctions.get(k).cloned().unwrap_or_default(),
        })
        .collect();
    let editing = crate::vg_family::copy_assign_pipeline::detect_editing_columns(read_obs, &copies);
    let graph = super::copy_assign::BubbleGraph::from_copies(&copies);
    let evidence: Vec<super::copy_assign::ReadEvidence> = read_obs
        .iter()
        .enumerate()
        .map(|(i, obs)| {
            let rf = super::copy_assign::ReadFeatures {
                psv_obs: obs.clone(),
                psv_qual: vec![],
                junctions: read_junctions.get(i).cloned().unwrap_or_default(),
            };
            super::copy_assign::read_copy_evidence(&rf, &graph, &copies, params, &editing)
        })
        .collect();
    em_assign(&evidence, copy_alleles.len(), params.alpha, eps, max_iter)
}

#[cfg(test)]
mod em_assign_family_tests {
    use super::em_assign_family;
    use super::super::copy_assign::AssignParams;

    /// Planted 2-copy family, one distinguishing PSV column (A vs C): read 0 carries copy 0's
    /// allele at both columns, read 1 carries copy 1's -- the `--em` binary flag's exact input
    /// shape (`FamilyAssignment::read_psv_obs` / `copy_psv_alleles`, no junctions).
    #[test]
    fn em_assign_family_recovers_planted_copies() {
        let copy_alleles: Vec<Vec<Option<u8>>> = vec![
            vec![Some(b'A'), Some(b'A')],
            vec![Some(b'C'), Some(b'C')],
        ];
        let read_obs: Vec<Vec<Option<u8>>> = vec![
            vec![Some(b'A'), Some(b'A')],
            vec![Some(b'C'), Some(b'C')],
        ];
        let params = AssignParams::for_alpha(1e-3);
        let r = em_assign_family(&read_obs, &copy_alleles, &[], &[], &params, 1e-6, 200);

        assert_eq!(r.posteriors.len(), 2);
        let argmax = |row: &Vec<f64>| {
            row.iter().enumerate().max_by(|a, b| a.1.partial_cmp(b.1).unwrap()).map(|(k, _)| k).unwrap()
        };
        assert_eq!(argmax(&r.posteriors[0]), 0, "read 0 (all-A) must favor copy 0");
        assert_eq!(argmax(&r.posteriors[1]), 1, "read 1 (all-C) must favor copy 1");

        let abundance_sum: f64 = r.abundances.iter().sum();
        assert!((abundance_sum - 1.0).abs() < 1e-9, "abundances must sum to 1: {abundance_sum}");
    }

    /// Task H2 TDD pin: a copy pair distinguished ONLY by a novel junction (identical PSVs, so
    /// PSV-only evidence is a flat tie) is unresolvable when junctions are not threaded (empty
    /// slices -- the OLD `em_assign_family` behavior, byte-identical), and becomes `Certified` /
    /// correctly argmax'd once the real per-read/per-copy junctions are threaded through. This is
    /// the O2<->O3 harmony this task ships: a copy defined by a splice, not a PSV, is now
    /// resolvable by the EM exactly as the one-shot gate already resolves it
    /// (`copy_assign::junction_only_resolves`).
    #[test]
    fn junction_only_copy_is_resolved_by_em() {
        // Both copies carry the SAME allele at the one PSV column -- no PSV signal distinguishes
        // them at all.
        let copy_alleles: Vec<Vec<Option<u8>>> = vec![vec![Some(b'A')], vec![Some(b'A')]];
        // Copy-specific junctions: copy 0 splices at offset 100, copy 1 at offset 200 (well outside
        // the default boundary_tol=4, so they never cross-match).
        let copy_junctions: Vec<Vec<i64>> = vec![vec![100], vec![200]];

        // 4 reads, no PSV signal (all carry the shared allele); half carry copy 0's junction, half
        // copy 1's.
        let read_obs: Vec<Vec<Option<u8>>> = vec![vec![Some(b'A')]; 4];
        let read_junctions: Vec<Vec<i64>> =
            vec![vec![100], vec![100], vec![200], vec![200]];

        let params = AssignParams::for_alpha(1e-3);

        // (a) empty junction slices == the previous PSV-only behavior: PSV alleles are identical
        // across copies, so no read has any decisive feature -> every read stays SoftZone.
        let no_junctions = em_assign_family(&read_obs, &copy_alleles, &[], &[], &params, 1e-6, 200);
        assert!(
            no_junctions.labels.iter().all(|l| matches!(l, super::EmLabel::SoftZone)),
            "PSV-only (empty junction slices) must leave an identical-PSV copy pair unresolved: {:?}",
            no_junctions.labels
        );

        // (b) threading the real junctions resolves every read, argmax'ing the copy whose
        // junction it carries.
        let with_junctions = em_assign_family(
            &read_obs,
            &copy_alleles,
            &read_junctions,
            &copy_junctions,
            &params,
            1e-6,
            200,
        );
        let argmax = |row: &Vec<f64>| {
            row.iter().enumerate().max_by(|a, b| a.1.partial_cmp(b.1).unwrap()).map(|(k, _)| k).unwrap()
        };
        let expected_copy = [0usize, 0, 1, 1];
        for (i, &exp) in expected_copy.iter().enumerate() {
            assert!(
                matches!(with_junctions.labels[i], super::EmLabel::Certified),
                "read {i} must be Certified once its junction is threaded: {:?}",
                with_junctions.labels[i]
            );
            assert_eq!(
                argmax(&with_junctions.posteriors[i]),
                exp,
                "read {i} must argmax the copy whose junction it carries"
            );
        }
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

/// Coverage-sweep DEMONSTRATION (Task 6): shows the EM soft-relaxation converges to the planted
/// truth as coverage grows -- the "it works with enough data, not an opinion" evidence for the
/// advisor, and (since it exercises `read_copy_evidence` -> `em_assign` end-to-end on simulated
/// reads) it subsumes the integration test. Fully deterministic: a seeded LCG stands in for `rand`
/// (no crate, no system time), so re-running prints byte-identical numbers.
///
/// Two planted families:
/// - a Strong-Separated K=3 family (every pair of copies differs at >=2 of 8 PSV columns) with
///   known abundances `pi_star = [0.5, 0.3, 0.2]`. As coverage (reads/copy-share unit) grows from
///   2 to 100, the recovered abundance should converge to `pi_star` (L1 -> 0) and `Certified` reads
///   should be assigned to their true copy essentially always -- the `bench/em_consistency.md`
///   convergence theorem, confirmed empirically.
/// - a K=0 family (two copies with IDENTICAL allele vectors -- no distinguishing bubble at all):
///   no amount of coverage can identify these reads, so every read must stay `SoftZone` and the
///   posterior must track the abundance prior exactly (no forced 1/k hard call), at every coverage.
#[cfg(test)]
mod coverage_sweep {
    use super::super::copy_assign::{read_copy_evidence, AssignParams, BubbleGraph, CopyProfile, ReadFeatures};
    use super::{em_assign, EmLabel};

    /// Seeded 64-bit LCG (Knuth MMIX constants) -- deterministic stand-in for `rand`. Returns the
    /// top 31 bits of the updated state (avoids the low-bit periodicity LCGs are known for).
    fn lcg(s: &mut u64) -> u64 {
        *s = s.wrapping_mul(6364136223846793005).wrapping_add(1442695040888963407);
        *s >> 33
    }

    /// Uniform draw in `[0, 1)` from the LCG stream.
    fn lcg_unit(s: &mut u64) -> f64 {
        (lcg(s) as f64) / ((1u64 << 31) as f64)
    }

    const BASES: [u8; 4] = [b'A', b'C', b'G', b'T'];
    const N_COLS: usize = 8;

    /// Strong-Separated K=3 family: copy0 = all A, copy1 = all C, copy2 = G,G,G,G,A,A,A,A.
    /// Pairwise Hamming distances over the 8 columns are 8 (c0,c1), 4 (c0,c2), 8 (c1,c2) -- every
    /// pair clears the `>= 2` Strong-Separation bar with wide margin.
    fn planted_copies() -> Vec<CopyProfile> {
        let c0 = vec![Some(b'A'); N_COLS];
        let c1 = vec![Some(b'C'); N_COLS];
        let mut c2 = vec![Some(b'G'); N_COLS];
        for a in c2.iter_mut().skip(4) {
            *a = Some(b'A');
        }
        vec![
            CopyProfile { copy_id: 0, alleles: c0, junctions: vec![] },
            CopyProfile { copy_id: 1, alleles: c1, junctions: vec![] },
            CopyProfile { copy_id: 2, alleles: c2, junctions: vec![] },
        ]
    }

    /// K=0 family: two copies with byte-identical allele vectors -- no PSV column distinguishes
    /// them, so no read can ever carry evidence for one over the other.
    fn k0_copies() -> Vec<CopyProfile> {
        vec![
            CopyProfile { copy_id: 0, alleles: vec![Some(b'A'); N_COLS], junctions: vec![] },
            CopyProfile { copy_id: 1, alleles: vec![Some(b'A'); N_COLS], junctions: vec![] },
        ]
    }

    /// Draw a true copy index `z ~ pi_star` from the LCG stream.
    fn draw_copy(pi_star: &[f64], seed: &mut u64) -> usize {
        let r = lcg_unit(seed);
        let mut acc = 0.0;
        for (k, &pi) in pi_star.iter().enumerate() {
            acc += pi;
            if r < acc {
                return k;
            }
        }
        pi_star.len() - 1
    }

    /// Corrupt `base` at per-base error rate `err_rate`: with probability `err_rate`, replace it
    /// with a uniformly-drawn *different* base (LCG-driven); otherwise leave it untouched.
    fn corrupt(base: u8, err_rate: f64, seed: &mut u64) -> u8 {
        if lcg_unit(seed) < err_rate {
            loop {
                let cand = BASES[(lcg(seed) % 4) as usize];
                if cand != base {
                    return cand;
                }
            }
        } else {
            base
        }
    }

    /// Simulate `n_reads` reads from `copies` under `pi_star`, each column independently corrupted
    /// at `err_rate`. Returns `(true_copy, ReadFeatures)` pairs (truth kept alongside for scoring,
    /// never fed to the EM engine itself). `psv_qual = None` throughout so the likelihood model
    /// uses `p.error_rate` directly -- set that to the SAME `err_rate` so the emission model the EM
    /// engine assumes matches the process that generated the reads.
    fn simulate_reads(
        copies: &[CopyProfile],
        pi_star: &[f64],
        n_reads: usize,
        err_rate: f64,
        seed: &mut u64,
    ) -> Vec<(usize, ReadFeatures)> {
        (0..n_reads)
            .map(|_| {
                let z = draw_copy(pi_star, seed);
                let psv_obs: Vec<Option<u8>> = copies[z]
                    .alleles
                    .iter()
                    .map(|&a| Some(corrupt(a.expect("planted copy defines every column"), err_rate, seed)))
                    .collect();
                let psv_qual = vec![None; psv_obs.len()];
                (z, ReadFeatures { psv_obs, psv_qual, junctions: vec![] })
            })
            .collect()
    }

    #[test]
    fn coverage_sweep_demonstrates_em_convergence() {
        let pi_star = [0.5_f64, 0.3, 0.2];
        let err_rate = 0.01;
        let copies = planted_copies();
        // error_rate matches the simulator's err_rate exactly (psv_qual is None throughout, so this
        // is the only per-base error the likelihood model uses); alpha = 1e-3 per the task recap.
        let params = AssignParams { error_rate: err_rate, ..AssignParams::for_alpha(1e-3) };
        let editing_cols = vec![false; N_COLS];
        let coverages = [2usize, 5, 10, 20, 50, 100];
        let mut seed: u64 = 0xC0FFEE_2026_u64;
        let graph = BubbleGraph::from_copies(&copies);

        let mut l1_by_cov = std::collections::HashMap::new();
        let mut acc_by_cov = std::collections::HashMap::new();

        println!(
            "{:>8} {:>8} {:>15} {:>15} {:>14}",
            "coverage", "n_reads", "certified_frac", "certified_acc", "abundance_L1"
        );
        for &cov in &coverages {
            let n_reads = cov * 100;
            let sim = simulate_reads(&copies, &pi_star, n_reads, err_rate, &mut seed);
            let evidence: Vec<_> = sim
                .iter()
                .map(|(_, r)| read_copy_evidence(r, &graph, &copies, &params, &editing_cols))
                .collect();
            let truth: Vec<usize> = sim.iter().map(|(z, _)| *z).collect();

            let result = em_assign(&evidence, 3, 1e-3, 1e-6, 500);

            let mut n_cert = 0usize;
            let mut n_cert_correct = 0usize;
            for (i, lbl) in result.labels.iter().enumerate() {
                if matches!(lbl, EmLabel::Certified) {
                    n_cert += 1;
                    let argmax = result.posteriors[i]
                        .iter()
                        .enumerate()
                        .max_by(|a, b| a.1.partial_cmp(b.1).unwrap())
                        .map(|(k, _)| k)
                        .unwrap();
                    if argmax == truth[i] {
                        n_cert_correct += 1;
                    }
                }
            }
            let certified_frac = n_cert as f64 / n_reads as f64;
            let certified_acc = if n_cert > 0 { n_cert_correct as f64 / n_cert as f64 } else { f64::NAN };
            let l1: f64 = result.abundances.iter().zip(pi_star.iter()).map(|(h, s)| (h - s).abs()).sum();

            println!(
                "{:>8} {:>8} {:>15.4} {:>15.4} {:>14.4}",
                cov, n_reads, certified_frac, certified_acc, l1
            );

            l1_by_cov.insert(cov, l1);
            acc_by_cov.insert(cov, certified_acc);
        }

        let l1_2 = l1_by_cov[&2];
        let l1_100 = l1_by_cov[&100];
        let acc_100 = acc_by_cov[&100];

        assert!(
            l1_100 < l1_2 * 0.5,
            "abundance L1 must shrink substantially with coverage: cov2={l1_2} cov100={l1_100}"
        );
        assert!(l1_100 < 0.05, "abundance L1 at coverage 100 should be small: {l1_100}");
        assert!(
            acc_100 >= 0.95,
            "certified-read accuracy at high coverage should be near-perfect: {acc_100}"
        );
    }

    #[test]
    fn k0_family_stays_soft_zone_at_all_coverages() {
        let copies = k0_copies();
        let pi_star = [0.5_f64, 0.5];
        let err_rate = 0.01;
        let params = AssignParams { error_rate: err_rate, ..AssignParams::for_alpha(1e-3) };
        let editing_cols = vec![false; N_COLS];
        let mut seed: u64 = 0xDEADBEEF_u64;
        let graph = BubbleGraph::from_copies(&copies);

        for &cov in &[2usize, 10, 100] {
            let n_reads = cov * 10;
            let sim = simulate_reads(&copies, &pi_star, n_reads, err_rate, &mut seed);
            let evidence: Vec<_> = sim
                .iter()
                .map(|(_, r)| read_copy_evidence(r, &graph, &copies, &params, &editing_cols))
                .collect();
            let result = em_assign(&evidence, 2, 1e-3, 1e-6, 500);

            assert!(
                result.labels.iter().all(|l| matches!(l, EmLabel::SoftZone)),
                "K=0 family (identical alleles) must never certify a read, coverage={cov}"
            );
            for row in &result.posteriors {
                assert!(
                    (row[0] - result.abundances[0]).abs() < 1e-6,
                    "posterior must equal the abundance prior exactly (no forced 1/k call), coverage={cov}"
                );
            }
        }
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
