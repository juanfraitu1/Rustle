//! Collapse gate: is a single-rep locus actually several copies the aligner could not separate?
//!
//! SDA (Vollger et al., Nat Methods 2019) detects a collapse by read-depth excess and only THEN defines PSVs,
//! "requiring sequence coverages consistent with a single-copy locus in order to distinguish PSVs from allelic
//! variants". We ran that second stage alone, and a single-copy gene (TSPYL1) reported 12 collapsed copies
//! against DAZ's 3 — see `bench/COLLAPSED_COPY_GATE.md`. **Collapse first, haplotypes second.**
//!
//! In DNA a collapse shows as excess depth. In RNA depth is copy number × expression (Clair3-RNA: "the coverage
//! is uneven across genomic regions in RNA-seq"), and allele-specific expression destroys the allelic balance
//! that would otherwise separate a het allele from a PSV ("zygosity flipping can happen"). So we keep SDA's
//! structure and change its instrument: a collapse shows instead as **reads the aligner cannot place uniquely**.
//!
//! Measured on GGO Iso-Seq (primary records only): **0 MAPQ-0 primaries across 9449 reads** at five single-copy
//! loci (TSPYL1, DERPC, ATXN7L3B, GSPT2, EEF1A1), against 19/20 at DAZ2 and 30/34 at TSPY. The statistic is
//! expression-invariant — TSPYL1 has 2151 reads and no ambiguity; DAZ2 has 20 reads and 95%.

use crate::vg_family::copy_assign::poisson_binomial_upper_tail;
use crate::vg_family::readonly_copy_number::chi_h;

/// Ambiguously-placed primary reads (`k`) out of primary reads (`n`) at a locus.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct Ambiguity {
    pub n: usize,
    pub k: usize,
}

/// Background per-read ambiguity rate, under a Jeffreys prior so it is never exactly zero.
///
/// `bg` is pooled over the region's uniquely-mappable reps (those that are not gate candidates). The controls
/// observe `k = 0`, whose MLE is 0 — and a zero rate would make a single stray MAPQ-0 read infinitely
/// significant. `None` when there is no background to estimate from: the caller must then ABSTAIN, never fire.
pub fn estimate_eps_amb(bg: Ambiguity) -> Option<f64> {
    if bg.n == 0 {
        return None;
    }
    Some((bg.k as f64 + 0.5) / (bg.n as f64 + 1.0))
}

/// Upper-tail probability of seeing `obs.k` or more ambiguous reads among `obs.n`, if the locus were unique.
///
/// Under the null `k ~ Binomial(n, eps_amb)`. Reuses the shipped Poisson-binomial tail with a constant
/// probability vector rather than adding a second distribution implementation. That routine is O(n²) in the
/// number of trials, and here a trial is a READ rather than a distinguishing position — but `k == 0` returns
/// immediately, and a clean locus is exactly the `k == 0` case, so uniquely-mapping loci cost nothing.
pub fn collapse_pvalue(obs: Ambiguity, eps_amb: f64) -> f64 {
    if obs.k == 0 || obs.n == 0 {
        return 1.0;
    }
    let probs = vec![eps_amb; obs.n];
    poisson_binomial_upper_tail(obs.k, &probs)
}

/// What the gate decided about a locus that has only one assembled rep.
#[derive(Clone, Debug, PartialEq)]
pub enum CollapseVerdict {
    /// Collapsed, and its reads resolve into `chi_h >= min_copies` conflicting haplotypes.
    Fire { chi_h: usize, p_value: f64 },
    /// The locus places its reads unambiguously, or its haplotypes do not reach `min_copies`.
    NotCollapsed { p_value: f64 },
    /// Cannot decide — never fire on an unbounded statistic.
    Abstain(&'static str),
}

/// Two legs, in SDA's order. Leg 2 is NOT consulted unless leg 1 fires: a single-copy gene reports plenty of
/// haplotypes (a het allele *is* a haplotype), and gating on them alone made TSPYL1 report 12 copies.
///
/// `haplotypes` are the `allele_vector`s of the **identifiable** copies at the locus.
///
/// χ(H) is a **lower bound** on copy number, never an estimate: two copies × two alleles also yields four
/// haplotypes. That is why the caller emits a copy NUMBER with reads certified tied, not an assignment.
pub fn collapse_verdict(
    obs: Ambiguity,
    bg: Ambiguity,
    haplotypes: &[Vec<Option<u8>>],
    alpha: f64,
    min_copies: usize,
) -> CollapseVerdict {
    // leg 1 — is this locus collapsed at all?
    let Some(eps_amb) = estimate_eps_amb(bg) else {
        return CollapseVerdict::Abstain("no uniquely-mappable rep to estimate the background ambiguity rate");
    };
    let p_value = collapse_pvalue(obs, eps_amb);
    if p_value >= alpha {
        return CollapseVerdict::NotCollapsed { p_value };
    }
    // leg 2 — and how many copies collapsed?
    let chi = chi_h(haplotypes);
    if chi < min_copies {
        return CollapseVerdict::NotCollapsed { p_value };
    }
    CollapseVerdict::Fire { chi_h: chi, p_value }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn eps_amb_is_never_zero_even_when_no_background_read_is_ambiguous() {
        // The five single-copy controls give 0 ambiguous reads in 9449. The MLE is 0, under which ONE stray
        // MAPQ-0 read would be infinitely significant. Jeffreys keeps it strictly positive.
        let eps = estimate_eps_amb(Ambiguity { n: 9449, k: 0 }).unwrap();
        assert!(eps > 0.0, "eps_amb must be strictly positive, got {eps}");
        assert!((eps - 0.5 / 9450.0).abs() < 1e-12, "Jeffreys: (k + 1/2) / (n + 1)");
    }

    #[test]
    fn eps_amb_abstains_without_background_reads() {
        assert_eq!(estimate_eps_amb(Ambiguity { n: 0, k: 0 }), None, "no background => cannot estimate => abstain");
    }

    #[test]
    fn collapse_pvalue_is_significant_for_daz2_and_not_for_a_clean_locus() {
        let eps = estimate_eps_amb(Ambiguity { n: 9449, k: 0 }).unwrap();
        let p_daz2 = collapse_pvalue(Ambiguity { n: 20, k: 19 }, eps); // DAZ2: 19 of 20 ambiguous
        assert!(p_daz2 < 1e-6, "DAZ2 must be overwhelmingly significant, got {p_daz2}");
        let p_clean = collapse_pvalue(Ambiguity { n: 2151, k: 0 }, eps); // TSPYL1
        assert!((p_clean - 1.0).abs() < 1e-12, "k = 0 => p = 1, got {p_clean}");
    }

    #[test]
    fn collapse_pvalue_of_a_single_stray_read_is_not_significant_at_alpha() {
        let eps = estimate_eps_amb(Ambiguity { n: 9449, k: 0 }).unwrap();
        let p = collapse_pvalue(Ambiguity { n: 500, k: 1 }, eps);
        assert!(p > 1e-3, "a single stray MAPQ-0 read must not fire the gate, got {p}");
    }

    /// Allele vectors: haplotypes differing at a shared column conflict, so `chi_h` counts them separately.
    fn haps(rows: &[&[u8]]) -> Vec<Vec<Option<u8>>> {
        rows.iter().map(|r| r.iter().map(|&b| if b == b'.' { None } else { Some(b) }).collect()).collect()
    }

    /// DAZ: 19/20 ambiguous, background clean, two distinguishable haplotypes.
    #[test]
    fn verdict_fires_on_a_collapsed_locus_with_two_haplotypes() {
        let v = collapse_verdict(
            Ambiguity { n: 20, k: 19 },
            Ambiguity { n: 9449, k: 0 },
            &haps(&[b"ACGT", b"ACGA"]),
            1e-3,
            2,
        );
        match v {
            CollapseVerdict::Fire { chi_h, .. } => assert_eq!(chi_h, 2),
            other => panic!("expected Fire, got {other:?}"),
        }
    }

    /// TSPYL1: a single-copy gene whose reads are ALL uniquely placed. Leg 2 would report many haplotypes
    /// (het alleles, editing, isoform noise) — leg 1 must stop it before leg 2 is ever consulted.
    #[test]
    fn verdict_rejects_a_unique_locus_however_many_haplotypes_it_reports() {
        let twelve: Vec<Vec<Option<u8>>> = (0..12u8).map(|i| vec![Some(b'A' + i), Some(b'C')]).collect();
        let v = collapse_verdict(Ambiguity { n: 2151, k: 0 }, Ambiguity { n: 9449, k: 0 }, &twelve, 1e-3, 2);
        assert!(matches!(v, CollapseVerdict::NotCollapsed { .. }), "unique locus must never fire, got {v:?}");
    }

    #[test]
    fn verdict_abstains_without_a_background_estimate() {
        let v =
            collapse_verdict(Ambiguity { n: 20, k: 19 }, Ambiguity { n: 0, k: 0 }, &haps(&[b"AC", b"AG"]), 1e-3, 2);
        assert!(matches!(v, CollapseVerdict::Abstain(_)), "no background => abstain, got {v:?}");
    }

    /// `min_copies` applies to χ(H), not to the rep count: one haplotype is not a family.
    #[test]
    fn verdict_rejects_a_collapse_that_resolves_to_one_haplotype() {
        let v = collapse_verdict(Ambiguity { n: 20, k: 19 }, Ambiguity { n: 9449, k: 0 }, &haps(&[b"ACGT"]), 1e-3, 2);
        assert!(
            matches!(v, CollapseVerdict::NotCollapsed { .. }),
            "chi_h = 1 < min_copies => no family, got {v:?}"
        );
    }
}
