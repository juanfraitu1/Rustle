//! O3 flag-pass detector (`docs/superpowers/specs/2026-09-10-o3-flag-pass-integration-design.md`):
//! ports `bench/o3_flag_pass.py`'s missing-copy detector natively into `copy_assign`.
//!
//! **STATUS:** INFRASTRUCTURE  (pure core functions for O3 modules; not independently reachable)

use crate::vg_family::allele_specific_junctions::lgamma;

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Class {
    Divergent,
    Structural,
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum LocusClass {
    OtherFamily,
    AnnotatedNoUnit,
    Unannotated,
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Flag {
    MissingCopy,
    Untestable,
    NoFlag,
}

#[derive(Clone, Debug)]
pub struct RawPair {
    pub family_id: String,
    pub copy_idx: String,
    pub is_partner: bool,
    pub n_rejected: usize,
    pub n_aligned: usize,
    pub covered_kb: f64,
    pub n_sites: usize,
    pub ctl_n: usize,
    pub ctl_covered_kb: f64,
    pub ctl_n_sites: usize,
    pub p_uncorrected: Option<f64>,
    pub class: Class,
    pub med_mismatch: f64,
    pub med_unaligned: f64,
}

#[derive(Clone, Debug)]
pub struct OrphanLocus {
    pub chrom: String,
    pub start: u64,
    pub end: u64,
    pub n_reads: usize,
    pub n_orphans: usize,
    pub class: LocusClass,
    pub n_genes_overlapping: usize,
    pub other_family_units: Vec<String>,
}

#[derive(Clone, Debug)]
pub struct FlaggedPair {
    pub pair: RawPair,
    pub flag: Flag,
    pub p_corrected_threshold: f64,
}

#[derive(Clone, Debug)]
pub struct O3Params {
    pub alpha: f64,
    pub max_reads: usize,
    pub min_reads: usize,
}

impl Default for O3Params {
    fn default() -> Self {
        O3Params { alpha: 0.001, max_reads: 500, min_reads: 3 }
    }
}

/// P(X >= k) for X ~ Poisson(lam). Mirrors `bench/o3_flag_pass.py`'s `poisson_tail` exactly (same
/// closed-form sum, same edge cases): k==0 -> 1.0 always true; lam<=0 with k>=1 -> 0.0 (a Poisson(0) is
/// a point mass at 0).
pub fn poisson_tail(k: usize, lam: f64) -> f64 {
    if k == 0 {
        return 1.0;
    }
    if lam <= 0.0 {
        return 0.0;
    }
    let mut s = 0.0_f64;
    for i in 0..k {
        s += (-lam + (i as f64) * lam.ln() - lgamma((i + 1) as f64)).exp();
    }
    (1.0 - s).clamp(0.0, 1.0)
}

/// Genome-wide Bonferroni aggregation (Phase 2): n_pairs = every RawPair with a computed p-value, from
/// the WHOLE run, not just one family. threshold = alpha / max(n_pairs, 1). Pure, no I/O.
pub fn finalize_flags(all_pairs: &[RawPair], alpha: f64) -> Vec<FlaggedPair> {
    let n_pairs = all_pairs.iter().filter(|p| p.p_uncorrected.is_some()).count();
    let threshold = alpha / (n_pairs.max(1) as f64);
    all_pairs
        .iter()
        .map(|p| {
            let flag = match p.p_uncorrected {
                Some(pv) if pv < threshold => Flag::MissingCopy,
                Some(_) => Flag::NoFlag,
                None => Flag::Untestable,
            };
            FlaggedPair { pair: p.clone(), flag, p_corrected_threshold: threshold }
        })
        .collect()
}

#[cfg(test)]
mod tests {
    use super::*;

    fn pair(p: Option<f64>) -> RawPair {
        RawPair {
            family_id: "F".into(), copy_idx: "0".into(), is_partner: false,
            n_rejected: 5, n_aligned: 5, covered_kb: 1.0, n_sites: 1,
            ctl_n: 5, ctl_covered_kb: 1.0, ctl_n_sites: 1,
            p_uncorrected: p, class: Class::Divergent, med_mismatch: 1.0, med_unaligned: 0.0,
        }
    }

    #[test]
    fn poisson_tail_k_zero_is_always_one() {
        assert_eq!(poisson_tail(0, 5.0), 1.0);
        assert_eq!(poisson_tail(0, 0.0), 1.0);
    }

    #[test]
    fn poisson_tail_zero_lambda_with_positive_k_is_zero() {
        assert_eq!(poisson_tail(3, 0.0), 0.0);
    }

    #[test]
    fn poisson_tail_matches_hand_computed_value() {
        // P(X >= 1) for Poisson(1) = 1 - P(X=0) = 1 - e^-1 = 0.6321205588...
        let p = poisson_tail(1, 1.0);
        assert!((p - 0.6321205588).abs() < 1e-6, "got {p}");
        // P(X >= 5) for Poisson(1) is small; hand value from scipy.stats.poisson.sf(4, 1) = 0.003659846827...
        let p2 = poisson_tail(5, 1.0);
        assert!((p2 - 0.003659846827).abs() < 1e-6, "got {p2}");
    }

    #[test]
    fn finalize_flags_labels_below_threshold_as_missing_copy() {
        // 2 testable pairs -> threshold = 0.001 / 2 = 0.0005
        let pairs = vec![pair(Some(0.0001)), pair(Some(0.9))];
        let flagged = finalize_flags(&pairs, 0.001);
        assert_eq!(flagged[0].flag, Flag::MissingCopy);
        assert_eq!(flagged[1].flag, Flag::NoFlag);
        assert!((flagged[0].p_corrected_threshold - 0.0005).abs() < 1e-12);
    }

    #[test]
    fn finalize_flags_untestable_when_p_is_none() {
        let pairs = vec![pair(None), pair(Some(0.0001))];
        let flagged = finalize_flags(&pairs, 0.001);
        assert_eq!(flagged[0].flag, Flag::Untestable);
        // n_pairs counts only the testable one, so threshold = 0.001 / 1 = 0.001
        assert!((flagged[1].p_corrected_threshold - 0.001).abs() < 1e-12);
    }

    #[test]
    fn finalize_flags_on_empty_input_does_not_divide_by_zero() {
        let flagged = finalize_flags(&[], 0.001);
        assert!(flagged.is_empty());
    }
}
