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

/// `chi_h`, but two copies also conflict when BOTH carry a non-empty, DIFFERING copy-specific junction set.
///
/// The reference-free copy COUNT must see the same evidence as the copy ASSIGNMENT it certifies. The gate and
/// the EM assign on PSVs **and** copy-specific junctions; `chi_h` sees only PSVs. On the recovered DAZ family
/// that contradiction is visible in the output: `n_copies = 2`, O2 assigns 2213 of 2353 placements across the
/// two copies, and yet `chi_H = 1`, so `famcn_readonly` reports **1 copy** for a two-copy family. DAZ1 and DAZ2
/// are near-identical exonically (one PSV column) and are separated by their junction structure.
///
/// The `both non-empty` guard mirrors the PSV rule's `both Some`: absence of junction evidence is not evidence
/// of sameness, so `chi_H` remains a valid LOWER bound (`max(n_loci, chi_H) <= true`). With no junction
/// evidence at all this is exactly [`chi_h`].
pub fn chi_h_with_junctions(copy_alleles: &[Vec<Option<u8>>], copy_junctions: &[Vec<i64>]) -> usize {
    fn psv_conflict(a: &[Option<u8>], b: &[Option<u8>]) -> bool {
        a.iter().zip(b.iter()).any(|(x, y)| matches!((x, y), (Some(xa), Some(yb)) if xa != yb))
    }
    fn junction_conflict(a: &[i64], b: &[i64]) -> bool {
        if a.is_empty() || b.is_empty() {
            return false; // no evidence is not evidence of difference
        }
        let (sa, sb): (std::collections::BTreeSet<_>, std::collections::BTreeSet<_>) =
            (a.iter().collect(), b.iter().collect());
        sa != sb
    }
    let empty: Vec<i64> = Vec::new();
    let junc = |i: usize| copy_junctions.get(i).unwrap_or(&empty);

    let mut groups: Vec<Vec<usize>> = Vec::new();
    for i in 0..copy_alleles.len() {
        let home = groups.iter().position(|g| {
            g.iter().all(|&m| {
                !psv_conflict(&copy_alleles[m], &copy_alleles[i]) && !junction_conflict(junc(m), junc(i))
            })
        });
        match home {
            Some(gi) => groups[gi].push(i),
            None => groups.push(vec![i]),
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
mod chi_h_junction_tests {
    use super::*;

    /// The recovered DAZ family: two copies that are exonically near-identical (one PSV column, PSV-compatible)
    /// but carry different junction structures. chi_h alone reports 1 for a family O2 assigns across 2 copies.
    #[test]
    fn chi_h_with_junctions_separates_psv_compatible_junction_distinct_copies() {
        let alleles = vec![vec![None], vec![None]];
        assert_eq!(chi_h(&alleles), 1, "PSV-only evidence cannot see the second copy");
        assert_eq!(chi_h_with_junctions(&alleles, &[vec![0i64], vec![7i64]]), 2, "junctions separate them");
    }

    /// Absence of junction evidence is not evidence of sameness -- chi_H stays a LOWER bound.
    #[test]
    fn chi_h_with_junctions_without_junction_evidence_is_exactly_chi_h() {
        let alleles = vec![vec![Some(b'A')], vec![Some(b'C')], vec![None]];
        assert_eq!(chi_h_with_junctions(&alleles, &[vec![], vec![], vec![]]), chi_h(&alleles));
        // one copy has junctions, the other does not => no junction conflict is inferred
        assert_eq!(chi_h_with_junctions(&[vec![None], vec![None]], &[vec![3i64], vec![]]), 1);
    }

    /// Identical junction sets do not manufacture a conflict, whatever the ordering.
    #[test]
    fn chi_h_with_junctions_ignores_order_and_identical_sets() {
        let alleles = vec![vec![None], vec![None]];
        assert_eq!(chi_h_with_junctions(&alleles, &[vec![5i64, 1], vec![1i64, 5]]), 1);
    }

    /// A PSV conflict still separates copies even when their junctions agree.
    #[test]
    fn chi_h_with_junctions_still_honours_psv_conflicts() {
        let alleles = vec![vec![Some(b'A')], vec![Some(b'G')]];
        assert_eq!(chi_h_with_junctions(&alleles, &[vec![2i64], vec![2i64]]), 2);
    }

    /// It can never exceed the copy count, so it remains a lower bound on true copy number.
    #[test]
    fn chi_h_with_junctions_never_exceeds_n_copies() {
        let alleles = vec![vec![Some(b'A')], vec![Some(b'C')], vec![Some(b'G')]];
        let j = vec![vec![1i64], vec![2i64], vec![3i64]];
        assert!(chi_h_with_junctions(&alleles, &j) <= alleles.len());
    }
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
        let result = em_assign_family(&reads, &distinguishable, &[], &[], &params, 1e-6, 500);
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
        let collapsed_result = em_assign_family(&collapsed_reads, &collapsed, &[], &[], &params, 1e-6, 500);

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

    /// Task H1 (VG-harmony pin): a copy NOT in the reference genome — admitted by O4's
    /// `absent_copy::admit_candidate` gate as a synthetic `DenovoTranscript` — must be a
    /// first-class copy to BOTH downstream consumers, not just threaded through pipeline
    /// plumbing that happens to compile.
    ///
    /// Wiring trace (verified by inspection; this test simulates the RESULT of that trace rather
    /// than re-running it, since the genome/BAM/minimap2 machinery is exercised by
    /// `absent_copy.rs`'s own hermetic tests):
    /// `admit_candidate_with_remap` (`absent_copy.rs`) returns `Admission::Copy(t)` once a
    /// candidate clears all five gates (cluster count, `min_p_distinct` from host, strand-
    /// symmetric spectrum, placeable overlay, remap identity `< 98%`) — `t` is a synthetic
    /// `DenovoTranscript` whose sequence bakes in the copy's distinguishing alleles.
    /// `denovo_pipeline.rs` (~line 626-654, behind `--absent-copies`) collects every
    /// `Admission::Copy(t)` into `admitted`, then extends the family's copy list (`all_copies`)
    /// with it and re-runs `assign_family_detailed` (Stage-2) over the augmented set — so the
    /// admitted copy sits at a real index alongside the reference copies, indistinguishable in
    /// type from them. From there `assign_family_detailed` -> `build_family_profiles`
    /// (`copy_assign_pipeline.rs`) extracts each copy's per-PSV-column allele vector
    /// (`CopyProfile.alleles`, i.e. `copy_psv_alleles`) for every copy in the set — there is no
    /// branch for "reference" vs "admitted-absent" once a `DenovoTranscript` exists. That shared
    /// `copy_alleles` vector is exactly what both `chi_h` (this module — O1's reference-free
    /// copy-number COUNT) and `em_assign_family` (`em_copy_assign` — O2's read ASSIGNMENT)
    /// consume next. This test starts from that shared post-admission state directly: a
    /// `copy_alleles` with 2 reference copies + 1 admitted-absent copy, and reads carrying the
    /// absent copy's alleles.
    #[test]
    fn absent_copy_is_assigned_and_counted() {
        let argmax = |row: &Vec<f64>| {
            row.iter()
                .enumerate()
                .max_by(|a, b| a.1.partial_cmp(b.1).unwrap())
                .map(|(k, _)| k)
                .unwrap()
        };

        // copy 0, copy 1 = "reference" copies (present in the linear genome); copy 2 = the
        // O4-admitted absent copy (its alleles came from `admit_candidate`'s synthetic
        // DenovoTranscript, not from any genome coordinate). Same private-column-per-copy layout
        // as `o1_o2_share_one_copy_object`: every pair differs at 2 of the 3 columns, well clear
        // of the alpha=1e-3, K=3 Bonferroni bound.
        let copy_alleles: Vec<Vec<Option<u8>>> = vec![
            vec![Some(b'C'), Some(b'A'), Some(b'A')], // ref copy 0
            vec![Some(b'A'), Some(b'C'), Some(b'A')], // ref copy 1
            vec![Some(b'A'), Some(b'A'), Some(b'C')], // absent copy 2 (O4-admitted)
        ];
        let absent_idx = 2;

        // (a) chi_h counts the absent copy as its own color: the reference-free copy number
        // rises to 3, not 2 -- the admitted copy is not silently dropped or merged.
        assert_eq!(
            chi_h(&copy_alleles),
            3,
            "reference-free copy number must count the O4-admitted absent copy as a distinct color"
        );

        // A handful of reads carrying the absent copy's exact alleles (as if minimap2/discovery
        // had routed them to this locus and copy_psv_alleles had extracted this vector for them).
        let absent_reads: Vec<Vec<Option<u8>>> = vec![
            copy_alleles[absent_idx].clone(),
            copy_alleles[absent_idx].clone(),
            copy_alleles[absent_idx].clone(),
        ];
        let params = AssignParams::for_alpha(1e-3);
        let result = em_assign_family(&absent_reads, &copy_alleles, &[], &[], &params, 1e-6, 500);

        // (b) the EM assigns those reads to the absent copy, confidently.
        for (i, row) in result.posteriors.iter().enumerate() {
            assert_eq!(
                argmax(row),
                absent_idx,
                "read {i} carrying the absent copy's alleles must be assigned to the absent copy"
            );
            assert!(
                matches!(result.labels[i], EmLabel::Certified),
                "read {i} must be Certified, not stuck in the K-frontier soft zone: {:?}",
                result.labels[i]
            );
        }
        assert!(
            result.abundances[absent_idx] > 0.0,
            "the EM's recovered abundance for the absent copy must be > 0: {:?}",
            result.abundances
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
