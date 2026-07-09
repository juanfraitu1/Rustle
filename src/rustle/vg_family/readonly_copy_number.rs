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
