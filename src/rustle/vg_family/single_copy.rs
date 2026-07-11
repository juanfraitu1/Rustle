//! Single-copy baseline: the χ(H)=1 loci that calibrate copy number.
//!
//! A TRANSCRIPT is the shipped `DenovoTranscript` object, defined by the `assemble_gate` predicate
//! (`denovo_assemble.rs`): an exact-intron-chain read cluster whose LOCUS (junction-incidence component,
//! `locus_support`) carries >= GATE_MIN_READS reads, whose junctions are all canonical and consistent-strand,
//! and whose spliced length is in [MIN_SPLICED, MAX_SPLICED]. A locus is a copy family of size >= 1:
//! single-copy is the degenerate χ(H)=1 case (0 PSVs); its transcripts are the λ_global baseline that
//! calibrates `depth_cn = E_fam / λ_global`. This module is the copy-number baseline, NOT an isoform catalog.

use super::denovo_pipeline::ColocatedFamily;
use super::family_detect::DenovoTranscript;

/// One single-copy locus (χ(H)=1). Carries NO `seq` — this is a lightweight baseline record, so a genome-wide
/// accumulator of ~100k of these is cheap.
#[derive(Clone, Debug, PartialEq)]
pub struct SingleCopyLocus {
    pub chrom: String,
    pub start: u64,
    pub end: u64,
    pub strand: char,
    pub n_reads: u32,
    pub n_exons: usize,
}

/// Reps that NO family claims, as lightweight `SingleCopyLocus` records. Membership is by `tid`.
///
/// `rep_totals[i]` is the LOCUS TOTAL read count for `reps[i]` — the summed reads over every isoform that
/// collapsed into that rep's locus (from `collapse_loci_span_aware_with_totals`), NOT the rep isoform's own
/// count. This is the expression basis λ_global needs: a gene's total, so it is on the same footing as the
/// family-total `E_fam` in `depth_cn = E_fam / λ_global`.
pub fn single_copy_loci(
    reps: &[DenovoTranscript],
    rep_totals: &[u32],
    families: &[ColocatedFamily],
) -> Vec<SingleCopyLocus> {
    assert_eq!(reps.len(), rep_totals.len(), "rep_totals must be parallel to reps");
    let claimed: std::collections::HashSet<&str> =
        families.iter().flat_map(|f| f.copies.iter().map(|c| c.tid.as_str())).collect();
    reps.iter()
        .zip(rep_totals.iter())
        .filter(|(r, _)| !claimed.contains(r.tid.as_str()))
        .map(|(r, &total)| SingleCopyLocus {
            chrom: r.chrom.clone(),
            start: r.start,
            end: r.end,
            strand: r.strand,
            n_reads: total,
            n_exons: r.introns.len() + 1,
        })
        .collect()
}

/// Median read count over the single-copy loci — the genome-wide single-copy expression floor `λ_global`.
/// `None` on an empty set.
pub fn lambda_global(loci: &[SingleCopyLocus]) -> Option<f64> {
    if loci.is_empty() {
        return None;
    }
    let mut v: Vec<u32> = loci.iter().map(|l| l.n_reads).collect();
    v.sort_unstable();
    let n = v.len();
    Some(if n % 2 == 1 {
        v[n / 2] as f64
    } else {
        (v[n / 2 - 1] as f64 + v[n / 2] as f64) / 2.0
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    fn tx(
        tid: &str,
        chrom: &str,
        start: u64,
        end: u64,
        n_reads: u32,
        strand: char,
        introns: Vec<(u64, u64)>,
    ) -> DenovoTranscript {
        DenovoTranscript { tid: tid.into(), chrom: chrom.into(), start, end, n_reads, strand, introns, seq: vec![] }
    }
    fn fam(id: &str, copies: Vec<DenovoTranscript>) -> ColocatedFamily {
        let chrom = copies[0].chrom.clone();
        ColocatedFamily { family_id: id.into(), chrom, start: 0, end: 0, copies }
    }

    #[test]
    fn single_copy_loci_are_the_reps_no_family_claims() {
        let a = tx("a", "c1", 100, 200, 10, '+', vec![(120, 150)]);
        let b = tx("b", "c1", 300, 400, 20, '-', vec![]);
        let c = tx("c", "c1", 500, 600, 30, '+', vec![(520, 540), (560, 580)]);
        let families = vec![fam("F0", vec![a.clone(), b.clone()])];
        let sc = single_copy_loci(&[a, b, c], &[10, 20, 30], &families);
        assert_eq!(sc.len(), 1);
        assert_eq!(sc[0].chrom, "c1");
        assert_eq!(sc[0].start, 500);
        assert_eq!(sc[0].end, 600);
        assert_eq!(sc[0].strand, '+');
        assert_eq!(sc[0].n_reads, 30);
        assert_eq!(sc[0].n_exons, 3, "n_exons = introns.len() + 1");
    }

    #[test]
    fn single_copy_locus_n_reads_is_the_locus_total_not_the_rep_isoform() {
        // A single-copy locus whose rep isoform holds 30 reads but whose locus (all collapsed isoforms) totals
        // 62. n_reads must be the total, matching the E_fam basis in depth_cn = E_fam / lambda_global.
        let c = tx("c", "c1", 500, 600, 30, '+', vec![(520, 540)]);
        let sc = single_copy_loci(&[c], &[62], &[]);
        assert_eq!(sc.len(), 1);
        assert_eq!(sc[0].n_reads, 62, "n_reads = locus total (summed isoforms), not the rep's 30");
    }

    #[test]
    fn single_copy_loci_membership_is_by_tid() {
        let a = tx("a", "c1", 100, 200, 10, '+', vec![]);
        let claimed = tx("a", "c1", 100, 200, 10, '+', vec![]);
        let other = tx("z", "c1", 100, 200, 10, '+', vec![]);
        let families = vec![fam("F0", vec![claimed, other])];
        assert!(
            single_copy_loci(&[a], &[10], &families).is_empty(),
            "a's tid is claimed -> not single-copy"
        );
    }

    #[test]
    fn single_copy_locus_carries_no_seq() {
        let a = tx("a", "c1", 1, 2, 5, '+', vec![]);
        let sc = single_copy_loci(&[a], &[5], &[]);
        let _ = SingleCopyLocus {
            chrom: sc[0].chrom.clone(),
            start: sc[0].start,
            end: sc[0].end,
            strand: sc[0].strand,
            n_reads: sc[0].n_reads,
            n_exons: sc[0].n_exons,
        };
    }

    #[test]
    fn lambda_global_is_the_median_n_reads() {
        let mk =
            |n: u32| SingleCopyLocus { chrom: "c".into(), start: 0, end: 1, strand: '+', n_reads: n, n_exons: 1 };
        assert_eq!(lambda_global(&[mk(10), mk(20), mk(30), mk(40)]), Some(25.0), "even -> mean of middle two");
        assert_eq!(lambda_global(&[mk(30), mk(10), mk(20)]), Some(20.0), "odd, unsorted -> middle after sort");
        assert_eq!(lambda_global(&[]), None, "empty -> None");
        assert_eq!(lambda_global(&[mk(7)]), Some(7.0));
    }
}
