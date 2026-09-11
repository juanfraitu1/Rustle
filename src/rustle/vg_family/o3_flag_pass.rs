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

/// Parses one `(\d+)([=XIDNS])` CIGAR run list, e.g. "88=3D510=1I436=" -> [(88,'='),(3,'D'),(510,'='),(1,'I'),(436,'=')].
fn parse_cigar_ops(cg: &str) -> Vec<(u64, char)> {
    let mut out = Vec::new();
    let mut num = 0u64;
    for ch in cg.chars() {
        if ch.is_ascii_digit() {
            num = num * 10 + (ch as u64 - '0' as u64);
        } else {
            out.push((num, ch));
            num = 0;
        }
    }
    out
}

pub(crate) struct AlignmentSummary {
    pub covered_kb: f64,
    pub n_sites: usize,
    pub per_read: std::collections::HashMap<String, (usize, i64, usize)>,
}

/// Parses `minimap2 -x splice -c --eqx -N 1` PAF output. Target-position coverage/mismatch tallies (no
/// query-sequence lookup needed — see the module doc comment on why allele identity is dropped).
/// PAF columns used: [0]=query name [1]=query len [2]=query start [3]=query end [7]=target start,
/// [12..]=tags (the `cg:Z:` CIGAR tag).
pub(crate) fn parse_paf_consistency(paf_text: &str) -> AlignmentSummary {
    let mut cov: std::collections::HashMap<u64, u32> = std::collections::HashMap::new();
    let mut mism: std::collections::HashMap<u64, u32> = std::collections::HashMap::new();
    let mut per_read: std::collections::HashMap<String, (usize, i64, usize)> = std::collections::HashMap::new();
    for line in paf_text.lines() {
        let f: Vec<&str> = line.split('\t').collect();
        if f.len() < 13 {
            continue;
        }
        let cg = match f[12..].iter().find(|t| t.starts_with("cg:Z:")) {
            Some(t) => &t[5..],
            None => continue,
        };
        let qlen: i64 = f[1].parse().unwrap_or(0);
        let qstart: i64 = f[2].parse().unwrap_or(0);
        let qend: i64 = f[3].parse().unwrap_or(0);
        let mut t: u64 = f[7].parse().unwrap_or(0);
        let mut nx: usize = 0;
        for (num, op) in parse_cigar_ops(cg) {
            match op {
                '=' => {
                    for k in 0..num {
                        *cov.entry(t + k).or_insert(0) += 1;
                    }
                    t += num;
                }
                'X' => {
                    nx += num as usize;
                    for k in 0..num {
                        *cov.entry(t + k).or_insert(0) += 1;
                        *mism.entry(t + k).or_insert(0) += 1;
                    }
                    t += num;
                }
                'D' | 'N' => t += num,
                _ => {} // 'I'/'S' advance only the query position, which this function never tracks
            }
        }
        let unal = qlen - (qend - qstart);
        let aligned_len = (qend - qstart) as usize;
        let name = f[0].to_string();
        let better = per_read.get(&name).map_or(true, |&(_, _, prev)| aligned_len > prev);
        if better {
            per_read.insert(name, (nx, unal, aligned_len));
        }
    }
    let n_sites = mism
        .iter()
        .filter(|(p, &c)| c >= 3 && (c as f64) >= 0.5 * (*cov.get(p).unwrap_or(&0) as f64))
        .count();
    let covered_kb = cov.values().filter(|&&c| c >= 3).count() as f64 / 1000.0;
    AlignmentSummary { covered_kb, n_sites, per_read }
}

#[cfg(test)]
mod paf_tests {
    use super::*;

    #[test]
    fn single_clean_read_no_mismatches() {
        // read "r1", qlen 100, aligned 0..100 on target starting at 1000, 100 matched bases
        let paf = "r1\t100\t0\t100\t+\tY\t2000\t1000\t1100\t100\t100\t60\tcg:Z:100=";
        let s = parse_paf_consistency(paf);
        assert_eq!(s.n_sites, 0);
        assert_eq!(s.per_read.get("r1"), Some(&(0, 0, 100)));
        assert!((s.covered_kb - 0.0).abs() < 1e-9); // 100 positions >= coverage 3? NO - coverage is 1 here
    }

    #[test]
    fn covered_kb_requires_coverage_at_least_three() {
        // three reads all covering the same 10bp target window -> those 10 positions reach coverage 3
        let paf = concat!(
            "r1\t10\t0\t10\t+\tY\t100\t0\t10\t10\t10\t60\tcg:Z:10=\n",
            "r2\t10\t0\t10\t+\tY\t100\t0\t10\t10\t10\t60\tcg:Z:10=\n",
            "r3\t10\t0\t10\t+\tY\t100\t0\t10\t10\t10\t60\tcg:Z:10=",
        );
        let s = parse_paf_consistency(paf);
        assert!((s.covered_kb - 0.010).abs() < 1e-9, "10 positions at cov>=3 => 0.010 kb, got {}", s.covered_kb);
    }

    #[test]
    fn consistent_mismatch_site_needs_majority_and_floor_of_three() {
        // 4 reads: 3 mismatch at target pos 5, 1 matches -> count 3 >= 3 AND 3 >= 0.5*4 (cov=4) -> consistent
        let paf = concat!(
            "a\t10\t0\t10\t+\tY\t100\t0\t10\t10\t10\t60\tcg:Z:5=1X4=\n",
            "b\t10\t0\t10\t+\tY\t100\t0\t10\t10\t10\t60\tcg:Z:5=1X4=\n",
            "c\t10\t0\t10\t+\tY\t100\t0\t10\t10\t10\t60\tcg:Z:5=1X4=\n",
            "d\t10\t0\t10\t+\tY\t100\t0\t10\t10\t10\t60\tcg:Z:10=",
        );
        let s = parse_paf_consistency(paf);
        assert_eq!(s.n_sites, 1);
    }

    #[test]
    fn mismatch_below_majority_is_not_consistent() {
        // 4 reads, only 1 mismatches at pos 5 (cov 4, mismatch count 1 < 0.5*4=2) -> not consistent
        let paf = concat!(
            "a\t10\t0\t10\t+\tY\t100\t0\t10\t10\t10\t60\tcg:Z:5=1X4=\n",
            "b\t10\t0\t10\t+\tY\t100\t0\t10\t10\t10\t60\tcg:Z:10=\n",
            "c\t10\t0\t10\t+\tY\t100\t0\t10\t10\t10\t60\tcg:Z:10=\n",
            "d\t10\t0\t10\t+\tY\t100\t0\t10\t10\t10\t60\tcg:Z:10=",
        );
        let s = parse_paf_consistency(paf);
        assert_eq!(s.n_sites, 0);
    }

    #[test]
    fn keeps_the_longer_alignment_when_a_read_has_two_paf_lines() {
        let paf = concat!(
            "r1\t100\t0\t40\t+\tY\t1000\t0\t40\t40\t40\t60\tcg:Z:40=\n",
            "r1\t100\t0\t90\t+\tY\t1000\t100\t190\t90\t90\t60\tcg:Z:90=",
        );
        let s = parse_paf_consistency(paf);
        assert_eq!(s.per_read.get("r1").map(|&(_, _, len)| len), Some(90));
    }

    #[test]
    fn unaligned_bases_is_query_len_minus_aligned_span() {
        // qlen 100, aligned query 10..90 -> unal = 100 - (90-10) = 20
        let paf = "r1\t100\t10\t90\t+\tY\t1000\t0\t80\t80\t80\t60\tcg:Z:80=";
        let s = parse_paf_consistency(paf);
        assert_eq!(s.per_read.get("r1"), Some(&(0, 20, 80)));
    }
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

/// Orphan-locus classification (§ "Scope note" in the design doc): OtherFamily > AnnotatedNoUnit >
/// Unannotated precedence. `all_units_by_chrom`/`genes_by_chrom` built ONCE per run, shared read-only.
pub(crate) fn classify_orphan_locus(
    chrom: &str,
    start: u64,
    end: u64,
    own_family_id: &str,
    all_units_by_chrom: &std::collections::BTreeMap<String, Vec<(u64, u64, String, String)>>,
    genes_by_chrom: &std::collections::BTreeMap<String, Vec<(u64, u64)>>,
) -> (LocusClass, usize, Vec<String>) {
    let overlaps = |a0: u64, a1: u64, b0: u64, b1: u64| a0 < b1 && b0 < a1;
    let other_family_units: Vec<String> = all_units_by_chrom
        .get(chrom)
        .map(|v| {
            v.iter()
                .filter(|(s, e, fid, _)| overlaps(start, end, *s, *e) && fid != own_family_id)
                .take(3)
                .map(|(_, _, fid, cidx)| format!("{fid}:{cidx}"))
                .collect()
        })
        .unwrap_or_default();
    let n_genes_overlapping = genes_by_chrom
        .get(chrom)
        .map(|v| v.iter().filter(|(s, e)| overlaps(start, end, *s, *e)).count())
        .unwrap_or(0);
    let class = if !other_family_units.is_empty() {
        LocusClass::OtherFamily
    } else if n_genes_overlapping > 0 {
        LocusClass::AnnotatedNoUnit
    } else {
        LocusClass::Unannotated
    };
    (class, n_genes_overlapping, other_family_units)
}

#[cfg(test)]
mod locus_tests {
    use super::*;
    use std::collections::BTreeMap;

    fn units() -> BTreeMap<String, Vec<(u64, u64, String, String)>> {
        let mut m = BTreeMap::new();
        m.insert("chr1".to_string(), vec![(100, 200, "OTHERFAM".to_string(), "3".to_string())]);
        m
    }

    fn genes() -> BTreeMap<String, Vec<(u64, u64)>> {
        let mut m = BTreeMap::new();
        m.insert("chr1".to_string(), vec![(500, 600)]);
        m
    }

    #[test]
    fn overlapping_another_familys_unit_wins_other_family() {
        let (class, _, units_hit) = classify_orphan_locus("chr1", 150, 250, "MYFAM", &units(), &genes());
        assert_eq!(class, LocusClass::OtherFamily);
        assert_eq!(units_hit, vec!["OTHERFAM:3".to_string()]);
    }

    #[test]
    fn own_family_unit_does_not_count_as_other_family() {
        let mut u = units();
        u.get_mut("chr1").unwrap().push((150, 250, "MYFAM".to_string(), "0".to_string()));
        let (class, _, units_hit) = classify_orphan_locus("chr1", 150, 250, "MYFAM", &u, &genes());
        // still hits OTHERFAM's unit at 100-200 too, so still OtherFamily -- verifies own-family rows
        // are excluded, not that the whole overlap set is
        assert_eq!(class, LocusClass::OtherFamily);
        assert!(!units_hit.iter().any(|s| s.starts_with("MYFAM:")));
    }

    #[test]
    fn no_other_family_unit_but_gene_overlap_is_annotated_no_unit() {
        let (class, n_genes, units_hit) = classify_orphan_locus("chr1", 550, 650, "MYFAM", &units(), &genes());
        assert_eq!(class, LocusClass::AnnotatedNoUnit);
        assert_eq!(n_genes, 1);
        assert!(units_hit.is_empty());
    }

    #[test]
    fn nothing_overlapping_is_unannotated() {
        let (class, n_genes, units_hit) = classify_orphan_locus("chr1", 9000, 9100, "MYFAM", &units(), &genes());
        assert_eq!(class, LocusClass::Unannotated);
        assert_eq!(n_genes, 0);
        assert!(units_hit.is_empty());
    }

    #[test]
    fn other_family_units_capped_at_three() {
        let mut u = BTreeMap::new();
        u.insert(
            "chr1".to_string(),
            (0..5).map(|i| (100u64, 200u64, format!("FAM{i}"), "0".to_string())).collect(),
        );
        let (_, _, units_hit) = classify_orphan_locus("chr1", 100, 200, "MYFAM", &u, &BTreeMap::new());
        assert_eq!(units_hit.len(), 3);
    }

    #[test]
    fn unknown_chrom_is_unannotated() {
        let (class, n_genes, units_hit) = classify_orphan_locus("chrZZZ", 0, 10, "MYFAM", &units(), &genes());
        assert_eq!(class, LocusClass::Unannotated);
        assert_eq!(n_genes, 0);
        assert!(units_hit.is_empty());
    }
}
