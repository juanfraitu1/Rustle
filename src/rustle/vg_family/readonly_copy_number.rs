//! Reference-free per-family copy number (Task R1 of the reference-free copy-number plan).
//!
//! The advisor's key idea: even reads that can't be ASSIGNED to a specific copy still INDICATE
//! copy number. Two reference-free legs (no genome/assembly needed):
//!   * `chi_h`    -- distinguishable copies from the PSV conflict structure (Rust port of
//!                   `bench/family_copy_number.py`'s `copyonly_K`: the number of distinct
//!                   pairwise-conflicting hap-vectors). A LOWER bound -- it collapses identical/
//!                   compatible copies into one group.
//!                   families e.g. `chi_H=1` on a locus whose true copy number is ~11.
//!   * `depth_cn` -- read-depth estimate (E_fam / lambda_global) that recovers identical/
//!                   collapsed copies `chi_h` misses, because it counts ALL family reads
//!                   (including the unassignable ones) rather than distinct hap-vectors.
//!
//! Both legs are genome-free: `chi_h` only needs the per-copy PSV allele vectors already
//! computed for assignment, and `depth_cn` only needs a read count and an externally-supplied
//! single-copy expression floor (`lambda_global`, itself an RNA-only quantity -- see
//! `bench/rna_copy_number_depth.py::global_single_copy_anchor`).

/// Number of DISTINCT pairwise-conflicting copy hap-vectors (chi(H) / MCC, THEORY.md Lemma 1).
///
/// Two copies CONFLICT iff they differ at >= 1 PSV column where BOTH have an observed allele
/// (`Some(_)`); columns where either copy is `None` (no evidence) never contribute a conflict.
/// Copies are grouped so that within a group no two members conflict; `chi_h` is the resulting
/// group count. Implemented greedily (matches `copyonly_K`'s semantics: a family's colors are the
/// distinct pairwise-conflicting hap-vectors) -- place each copy in the first existing group none
/// of whose members it conflicts with, else start a new group.
pub fn chi_h(copy_alleles: &[Vec<Option<u8>>]) -> usize {
    fn conflicts(a: &[Option<u8>], b: &[Option<u8>]) -> bool {
        a.iter().zip(b.iter()).any(|(x, y)| matches!((x, y), (Some(xa), Some(yb)) if xa != yb))
    }

    let mut groups: Vec<Vec<&Vec<Option<u8>>>> = Vec::new();
    for alleles in copy_alleles {
        let home = groups.iter().position(|g| g.iter().all(|m| !conflicts(m, alleles)));
        match home {
            Some(gi) => groups[gi].push(alleles),
            None => groups.push(vec![alleles]),
        }
    }
    groups.len()
}

/// Reference-free read-depth copy-number estimate: `E_fam / lambda_global`.
///
/// `e_fam` = total reads over the family (E_fam, includes unassignable reads); `lambda_global` =
/// the genome-wide RNA single-copy expression floor (median n_reads over single-copy transcripts,
/// precomputed by `bench/rna_copy_number_depth.py::global_single_copy_anchor` -- a genome-wide RNA
/// quantity, NOT genomic). Returns `NaN` when `lambda_global <= 0` (undefined / not supplied).
pub fn depth_cn(e_fam: usize, lambda_global: f64) -> f64 {
    if lambda_global <= 0.0 {
        return f64::NAN;
    }
    e_fam as f64 / lambda_global
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::vg_family::copy_assign::AssignParams;
    use crate::vg_family::em_copy_assign::{em_assign_family, EmLabel};

    /// Task R2 (O1<->O2 harmony pin): `chi_h` (this module, O1's conflict-graph copy COUNT) and
    /// `em_assign_family` (`em_copy_assign`, O2's EM ASSIGNMENT) both consume the SAME per-copy
    /// hap-vector (`copy_alleles` / `CopyProfile.alleles`) -- there is exactly one copy object, not
    /// two independently-tunable ones. A better O1 (every copy pairwise-conflicts, no de-tie needed)
    /// must raise `chi_h` AND let the EM certify every read against that K; collapsing two copies to
    /// an identical hap-vector (a worse/over-merged O1) must drop `chi_h` by exactly 1 AND make the EM
    /// unable to separate reads from that pair -- both legs hit the K-frontier (`SoftZone`) together.
    /// This is a REGRESSION PIN: no new behavior, just nailing down that the two already compose.
    #[test]
    fn o1_o2_share_one_copy_object() {
        let argmax = |row: &Vec<f64>| {
            row.iter()
                .enumerate()
                .max_by(|a, b| a.1.partial_cmp(b.1).unwrap())
                .map(|(k, _)| k)
                .unwrap()
        };
        let params = AssignParams::for_alpha(1e-3);

        // --- Regime 1: distinguishable (good O1) -- 3 copies, each privately marked at its own
        // column against a shared background, so every pair conflicts at 2 of the 3 columns (well
        // clear of the K-frontier).
        let distinguishable: Vec<Vec<Option<u8>>> = vec![
            vec![Some(b'C'), Some(b'A'), Some(b'A')],
            vec![Some(b'A'), Some(b'C'), Some(b'A')],
            vec![Some(b'A'), Some(b'A'), Some(b'C')],
        ];
        assert_eq!(chi_h(&distinguishable), 3, "O1: 3 pairwise-conflicting copies -> chi_h = 3");

        // one read per true copy, carrying that copy's exact alleles.
        let reads = distinguishable.clone();
        let result = em_assign_family(&reads, &distinguishable, &params, 1e-6, 500);
        assert_eq!(
            result.abundances.len(),
            distinguishable.len(),
            "O2 (EM) must sum abundance over the SAME K copies O1's chi_h counts"
        );
        for (true_copy, row) in result.posteriors.iter().enumerate() {
            assert_eq!(argmax(row), true_copy, "read {true_copy} must be assigned to its true copy");
        }
        assert!(
            result.labels.iter().all(|l| matches!(l, EmLabel::Certified)),
            "good O1 (all copies pairwise-conflict) -> the EM certifies every read: {:?}",
            result.labels
        );

        // --- Regime 2: collapsed (worse O1) -- copy 1 and copy 2 forced to the IDENTICAL hap-vector
        // (a de-tie / over-merge failure upstream would produce exactly this).
        let mut collapsed = distinguishable.clone();
        collapsed[2] = collapsed[1].clone();
        assert_eq!(
            chi_h(&collapsed),
            chi_h(&distinguishable) - 1,
            "collapsing 2 copies to one hap-vector must drop chi_h by exactly 1 (3 -> 2)"
        );

        // reads carrying the POST-collapse hap-vectors: one from the untouched copy 0 (control),
        // and one from each of the two now-identical copies.
        let collapsed_reads: Vec<Vec<Option<u8>>> =
            vec![distinguishable[0].clone(), collapsed[1].clone(), collapsed[2].clone()];
        let collapsed_result = em_assign_family(&collapsed_reads, &collapsed, &params, 1e-6, 500);

        assert!(
            matches!(collapsed_result.labels[0], EmLabel::Certified),
            "the un-collapsed copy must remain identifiable: {:?}",
            collapsed_result.labels
        );
        // reads from EITHER of the two now-identical copies can no longer be separated: the EM hits
        // the K-frontier and abstains (SoftZone) on both, in lockstep with chi_h's drop.
        assert!(
            matches!(collapsed_result.labels[1], EmLabel::SoftZone),
            "read from collapsed copy 1 must go SoftZone: {:?}",
            collapsed_result.labels
        );
        assert!(
            matches!(collapsed_result.labels[2], EmLabel::SoftZone),
            "read from collapsed copy 2 must go SoftZone: {:?}",
            collapsed_result.labels
        );
    }

    #[test]
    fn chi_h_three_pairwise_distinct_private_alleles() {
        // each copy has a private allele at its own column -> all three pairwise conflict.
        let copies = vec![
            vec![Some(b'A'), Some(b'C'), Some(b'C')],
            vec![Some(b'C'), Some(b'A'), Some(b'C')],
            vec![Some(b'C'), Some(b'C'), Some(b'A')],
        ];
        assert_eq!(chi_h(&copies), 3);
    }

    #[test]
    fn chi_h_three_identical_copies_collapse_to_one() {
        let copies = vec![
            vec![Some(b'A'), Some(b'C')],
            vec![Some(b'A'), Some(b'C')],
            vec![Some(b'A'), Some(b'C')],
        ];
        assert_eq!(chi_h(&copies), 1);
    }

    #[test]
    fn chi_h_mixed_two_groups() {
        // [A,A] vs [A,C]: conflict at column 1 (A vs C). The two [A,C] copies agree everywhere ->
        // share a group. So groups = {[A,A]}, {[A,C], [A,C]} = 2.
        let copies = vec![
            vec![Some(b'A'), Some(b'A')],
            vec![Some(b'A'), Some(b'C')],
            vec![Some(b'A'), Some(b'C')],
        ];
        assert_eq!(chi_h(&copies), 2);
    }

    #[test]
    fn chi_h_empty_is_zero() {
        let copies: Vec<Vec<Option<u8>>> = Vec::new();
        assert_eq!(chi_h(&copies), 0);
    }

    #[test]
    fn chi_h_missing_evidence_never_conflicts() {
        // Columns with None on either side never contribute a conflict; these two copies never
        // disagree where both are observed, so they share a group.
        let copies = vec![
            vec![Some(b'A'), None, Some(b'C')],
            vec![None, Some(b'G'), Some(b'C')],
        ];
        assert_eq!(chi_h(&copies), 1);
    }

    #[test]
    fn depth_cn_basic_ratio() {
        assert!((depth_cn(100, 25.0) - 4.0).abs() < 1e-9);
    }

    #[test]
    fn depth_cn_zero_e_fam_is_zero() {
        assert_eq!(depth_cn(0, 10.0), 0.0);
    }

    #[test]
    fn depth_cn_nonpositive_lambda_is_nan() {
        assert!(depth_cn(50, 0.0).is_nan());
        assert!(depth_cn(50, -1.0).is_nan());
    }
}
