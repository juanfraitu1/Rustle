//! O3 flag-pass detector (`docs/superpowers/specs/2026-09-10-o3-flag-pass-integration-design.md`):
//! ports `bench/o3_flag_pass.py`'s missing-copy detector natively into `copy_assign`.
//!
//! **STATUS:** OPT-IN  (reachable via `copy_assign --flag-missing-copies`, src/bin/copy_assign.rs; default off)

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
/// `pub` (not `pub(crate)`): Task 5 (`src/bin/copy_assign.rs`) is a separate binary crate and calls this
/// directly per the integration plan's "Consumes" list -- `pub(crate)` here was invisible to it (E0603).
/// Pure visibility widening, no signature/behavior change; every other item Task 5 needs was already `pub`.
///
/// `own_family_ids` (final whole-branch-review fix round, follow-up): a SET, not a single id --
/// `bench/o3_flag_pass.py`'s own exclusion (`u[2] != cp[next(iter(cp))]['family_id'].split('_')[0]`) is
/// trivially always correct because each Python sweep is ISOLATED to exactly one catalog family. The
/// caller here is a locally co-located physical sweep group that can legitimately bundle MORE THAN ONE
/// true catalog family (the same fact Fix 1 addressed for the pair detector) -- passing a single id (the
/// group's own, possibly arbitrary, local label) would incorrectly classify an orphan cluster's OWN
/// bundled sibling family's real unit as "OtherFamily" whenever that sibling isn't the one label chosen to
/// represent the group. `own_family_ids` is every TRUE catalog family id actually present among the
/// calling group's own copies; a unit belongs to "this same group" (excluded from `OtherFamily`) iff its
/// `fid` is a member of this set.
pub fn classify_orphan_locus(
    chrom: &str,
    start: u64,
    end: u64,
    own_family_ids: &std::collections::HashSet<String>,
    all_units_by_chrom: &std::collections::BTreeMap<String, Vec<(u64, u64, String, String)>>,
    genes_by_chrom: &std::collections::BTreeMap<String, Vec<(u64, u64)>>,
) -> (LocusClass, usize, Vec<String>) {
    let overlaps = |a0: u64, a1: u64, b0: u64, b1: u64| a0 < b1 && b0 < a1;
    let other_family_units: Vec<String> = all_units_by_chrom
        .get(chrom)
        .map(|v| {
            v.iter()
                .filter(|(s, e, fid, _)| overlaps(start, end, *s, *e) && !own_family_ids.contains(fid))
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

    fn myfam() -> std::collections::HashSet<String> {
        std::collections::HashSet::from(["MYFAM".to_string()])
    }

    #[test]
    fn overlapping_another_familys_unit_wins_other_family() {
        let (class, _, units_hit) = classify_orphan_locus("chr1", 150, 250, &myfam(), &units(), &genes());
        assert_eq!(class, LocusClass::OtherFamily);
        assert_eq!(units_hit, vec!["OTHERFAM:3".to_string()]);
    }

    #[test]
    fn own_family_unit_does_not_count_as_other_family() {
        let mut u = units();
        u.get_mut("chr1").unwrap().push((150, 250, "MYFAM".to_string(), "0".to_string()));
        let (class, _, units_hit) = classify_orphan_locus("chr1", 150, 250, &myfam(), &u, &genes());
        // still hits OTHERFAM's unit at 100-200 too, so still OtherFamily -- verifies own-family rows
        // are excluded, not that the whole overlap set is
        assert_eq!(class, LocusClass::OtherFamily);
        assert!(!units_hit.iter().any(|s| s.starts_with("MYFAM:")));
    }

    #[test]
    fn no_other_family_unit_but_gene_overlap_is_annotated_no_unit() {
        let (class, n_genes, units_hit) = classify_orphan_locus("chr1", 550, 650, &myfam(), &units(), &genes());
        assert_eq!(class, LocusClass::AnnotatedNoUnit);
        assert_eq!(n_genes, 1);
        assert!(units_hit.is_empty());
    }

    #[test]
    fn nothing_overlapping_is_unannotated() {
        let (class, n_genes, units_hit) = classify_orphan_locus("chr1", 9000, 9100, &myfam(), &units(), &genes());
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
        let (_, _, units_hit) = classify_orphan_locus("chr1", 100, 200, &myfam(), &u, &BTreeMap::new());
        assert_eq!(units_hit.len(), 3);
    }

    #[test]
    fn unknown_chrom_is_unannotated() {
        let (class, n_genes, units_hit) = classify_orphan_locus("chrZZZ", 0, 10, &myfam(), &units(), &genes());
        assert_eq!(class, LocusClass::Unannotated);
        assert_eq!(n_genes, 0);
        assert!(units_hit.is_empty());
    }

    #[test]
    fn own_family_ids_set_excludes_every_member_not_just_the_first() {
        // Final whole-branch-review fix round, follow-up: a caller that bundles TWO true catalog
        // families (e.g. a locally co-located `fa` spanning more than one catalog family, the same fact
        // Fix 1 addressed for the pair detector) must have BOTH of its own family ids excluded from
        // `OtherFamily`, not just whichever single label happened to be chosen to represent the group.
        let mut u = BTreeMap::new();
        u.insert(
            "chr1".to_string(),
            vec![
                (100, 200, "SIBLING_B".to_string(), "0".to_string()), // bundled INTO the same caller group
                (100, 200, "TRULY_OTHER".to_string(), "0".to_string()), // NOT part of the caller's group
            ],
        );
        let own = std::collections::HashSet::from(["SIBLING_A".to_string(), "SIBLING_B".to_string()]);
        let (class, _, units_hit) = classify_orphan_locus("chr1", 100, 200, &own, &u, &BTreeMap::new());
        assert_eq!(class, LocusClass::OtherFamily, "TRULY_OTHER's unit is real external evidence");
        assert_eq!(
            units_hit,
            vec!["TRULY_OTHER:0".to_string()],
            "SIBLING_B must be excluded as self (it is in own_family_ids) even though it is not the \
             single id a bare-string comparison would have used"
        );
    }

    #[test]
    fn own_family_ids_set_with_no_external_unit_is_not_other_family() {
        // The bundle's OWN two sibling families both overlap the locus; with no truly external unit
        // present, the class must fall through to AnnotatedNoUnit/Unannotated, not OtherFamily -- the
        // single-id bug (comparing against only one arbitrary label) would have wrongly reported
        // OtherFamily here for whichever sibling wasn't the chosen label.
        let mut u = BTreeMap::new();
        u.insert(
            "chr1".to_string(),
            vec![
                (100, 200, "SIBLING_A".to_string(), "0".to_string()),
                (100, 200, "SIBLING_B".to_string(), "0".to_string()),
            ],
        );
        let own = std::collections::HashSet::from(["SIBLING_A".to_string(), "SIBLING_B".to_string()]);
        let (class, _, units_hit) = classify_orphan_locus("chr1", 100, 200, &own, &u, &genes());
        assert_eq!(class, LocusClass::Unannotated);
        assert!(units_hit.is_empty());
    }
}

/// One candidate copy Y's test batch (its own origin-rejected reads) and control batch (its own
/// certificate-accepted reads), both as `(name, sequence)` pairs ready to realign.
pub struct PairInput {
    pub copy_idx: String,
    pub is_partner: bool,
    pub rejected: Vec<(String, Vec<u8>)>,
    pub accepted: Vec<(String, Vec<u8>)>,
}

/// Realigns one batch of reads to a target window via `minimap2 -x splice -c --eqx -N 1`, mirroring
/// `copy_assign_pipeline.rs`'s `minimap2_msa_pair` temp-file convention (pid+atomic-nonce names,
/// `RUSTLE_MINIMAP2` env override, `Drop`-based cleanup) so concurrent region-parallel workers never
/// collide on the same path.
fn realign_batch(target_seq: &[u8], reads: &[(String, Vec<u8>)]) -> anyhow::Result<AlignmentSummary> {
    use std::io::Write;
    if reads.is_empty() {
        return Ok(AlignmentSummary { covered_kb: 0.0, n_sites: 0, per_read: std::collections::HashMap::new() });
    }
    let mm2 = std::env::var("RUSTLE_MINIMAP2").unwrap_or_else(|_| "minimap2".to_string());
    let dir = std::env::temp_dir();
    let pid = std::process::id();
    use std::sync::atomic::{AtomicUsize, Ordering};
    static O3_NONCE: AtomicUsize = AtomicUsize::new(0);
    let nonce = O3_NONCE.fetch_add(1, Ordering::Relaxed);
    let ypath = dir.join(format!("rustle_o3_y_{pid}_{nonce}.fa"));
    let rpath = dir.join(format!("rustle_o3_reads_{pid}_{nonce}.fa"));
    struct Cleanup(std::path::PathBuf, std::path::PathBuf);
    impl Drop for Cleanup {
        fn drop(&mut self) {
            let _ = std::fs::remove_file(&self.0);
            let _ = std::fs::remove_file(&self.1);
        }
    }
    let _cl = Cleanup(ypath.clone(), rpath.clone());
    {
        let mut y = std::fs::File::create(&ypath)?;
        y.write_all(b">Y\n")?;
        y.write_all(target_seq)?;
        y.write_all(b"\n")?;
        let mut r = std::fs::File::create(&rpath)?;
        for (name, seq) in reads {
            writeln!(r, ">{name}")?;
            r.write_all(seq)?;
            r.write_all(b"\n")?;
        }
    }
    // Minor (final whole-branch review): `-t 4` matches the Python reference's own hardcoded thread count,
    // but this codebase's existing `minimap2_msa_pair` pattern this function is modeled on uses `-t 1`.
    // Under `copy_assign --region-parallel N` this multiplies to up to 4xN minimap2 threads concurrently
    // -- on this project's 5-core WSL2 box (see the crash-rule notes), combining `--region-parallel` with
    // `--flag-missing-copies` could overload the machine. Left at `-t 4` (not changed here -- this is a
    // documentation-only note), but `--flag-missing-copies` runs on this machine should stay serial
    // (no `--region-parallel`) until/unless this is revisited.
    let out = std::process::Command::new(&mm2)
        .args(["-x", "splice", "-c", "--eqx", "-N", "1", "-t", "4"])
        .arg(&ypath)
        .arg(&rpath)
        .output()?;
    if !out.status.success() {
        anyhow::bail!("minimap2 failed: {}", String::from_utf8_lossy(&out.stderr));
    }
    Ok(parse_paf_consistency(&String::from_utf8_lossy(&out.stdout)))
}

fn median(mut v: Vec<f64>) -> Option<f64> {
    if v.is_empty() {
        return None;
    }
    v.sort_by(|a, b| a.partial_cmp(b).unwrap());
    Some(v[v.len() / 2])
}

/// Realignment target window for one Y-copy, extracted as a pure function so Task 7's window-selection
/// fix (using the catalog's L2 locus extent when present, matching `bench/o3_flag_pass.py`'s `detector()`
/// lines 45-48) is directly unit-testable without a real `minimap2`/`GenomeIndex`. `locus` is the
/// catalog's `(locus_start, locus_end)` for this copy when one is recorded; `longest_rejected` is the
/// length of the longest rejected read, used only in the `None` fallback to pad the bare copy span (so a
/// read realigned against the window cannot hang off its edge).
fn locus_or_padded_window(s: u64, e: u64, locus: Option<(u64, u64)>, longest_rejected: u64) -> (u64, u64) {
    match locus {
        Some((locus_s, locus_e)) => (locus_s.min(s), locus_e.max(e)),
        None => (s.saturating_sub(longest_rejected), e + longest_rejected),
    }
}

/// Phase 1: one family's raw (uncorrected) missing-copy pair statistics. Skips any Y with fewer than
/// `params.min_reads` rejected reads (matches the Python's `len(names) < 3: continue` -- not emitted as
/// Untestable, simply absent from the output). A `minimap2` failure for one Y is logged to stderr and
/// that pair is skipped -- see the design doc's Error Handling section for why this must not abort.
///
/// `copy_span_by_catalog_idx` maps to `(chrom, start, end, locus_extent)`: the realignment TARGET window
/// mirrors `bench/o3_flag_pass.py`'s `detector()` exactly (lines 45-48) -- `(min(locus_start, start),
/// max(locus_end, end))` when the catalog carries an L2 locus extent for this copy, else the copy's own
/// span padded by the longest rejected read on each side (see [`locus_or_padded_window`], which does this
/// computation and is unit-tested directly). Using the bare copy span unconditionally (the pre-fix
/// behaviour) mis-sizes the window whenever a copy's locus extent differs from its own span -- found
/// during Task 7's reproduction gate: `covered_kb`/`n_sites`/`rate_per_kb`/`p` all move even when
/// `n_rejected` matches exactly (e.g. MCL117_073244 copy 1: same 4 rejected reads, Python rate 4.27/kb
/// vs the unfixed Rust 59.81/kb).
pub fn detect_missing_copy_pairs(
    family_id: &str,
    copy_span_by_catalog_idx: &std::collections::HashMap<String, (String, u64, u64, Option<(u64, u64)>)>,
    genome: &crate::genome::GenomeIndex,
    inputs: &[PairInput],
    params: &O3Params,
) -> Vec<RawPair> {
    let mut out = Vec::new();
    for input in inputs {
        if input.rejected.len() < params.min_reads {
            continue;
        }
        let Some((chrom, s, e, locus)) = copy_span_by_catalog_idx.get(&input.copy_idx) else {
            continue;
        };
        // Task 7 fix (reproduction-gate finding): the Python caps each side at `max_reads` by NAME order
        // (`sorted(names)[:max_reads]` / `sorted(n for n in truth if ...)[:max_reads]`), not by collection
        // order -- irrelevant when a family stays under the cap (most pairs), but for a copy with more
        // than `max_reads` eligible reads (e.g. MCL121_073244 copy 0's 500-capped control pool) the two
        // orders sample a DIFFERENT subset, changing `covered_kb`/`n_sites`/`p` even when the uncapped
        // rejected side matches exactly. Sorting first makes the cap deterministic and Python-identical.
        let mut rejected_sorted = input.rejected.clone();
        rejected_sorted.sort_by(|a, b| a.0.cmp(&b.0));
        let rejected: Vec<_> = rejected_sorted.into_iter().take(params.max_reads).collect();
        let pad = rejected.iter().map(|(_, seq)| seq.len() as u64).max().unwrap_or(0);
        let (ls, le) = locus_or_padded_window(*s, *e, *locus, pad);
        let target = match genome.fetch_sequence(chrom, ls, le) {
            Some(t) => t,
            None => continue,
        };
        let test = match realign_batch(&target, &rejected) {
            Ok(s) => s,
            Err(err) => {
                eprintln!("[o3-flag-pass] realignment failed for {family_id}:{}: {err}", input.copy_idx);
                continue;
            }
        };
        let mut accepted_sorted = input.accepted.clone(); // same name-order cap as `rejected` above
        accepted_sorted.sort_by(|a, b| a.0.cmp(&b.0));
        let accepted: Vec<_> = accepted_sorted.into_iter().take(params.max_reads).collect();
        let ctl = if accepted.len() >= params.min_reads {
            match realign_batch(&target, &accepted) {
                Ok(s) => s,
                Err(err) => {
                    eprintln!("[o3-flag-pass] control realignment failed for {family_id}:{}: {err}", input.copy_idx);
                    AlignmentSummary { covered_kb: 0.0, n_sites: 0, per_read: std::collections::HashMap::new() }
                }
            }
        } else {
            AlignmentSummary { covered_kb: 0.0, n_sites: 0, per_read: std::collections::HashMap::new() }
        };
        // rate floored at 1 site over the control's covered kb, so a control with zero observed sites
        // never claims a zero rate (which would trivially pass every test) -- design doc, RawPair section.
        let p_uncorrected = if test.covered_kb > 0.0 && ctl.covered_kb > 0.0 {
            let rate = (ctl.n_sites.max(1) as f64) / ctl.covered_kb;
            Some(poisson_tail(test.n_sites, rate * test.covered_kb))
        } else {
            None
        };
        let mismatches: Vec<f64> = test.per_read.values().map(|&(nx, _, _)| nx as f64).collect();
        let unaligned: Vec<f64> = test.per_read.values().map(|&(_, unal, _)| unal as f64).collect();
        let med_mismatch = median(mismatches).unwrap_or(f64::NAN);
        let med_unaligned = median(unaligned).unwrap_or(f64::NAN);
        let class = if !med_mismatch.is_nan() && med_mismatch > med_unaligned {
            Class::Divergent
        } else {
            Class::Structural
        };
        out.push(RawPair {
            family_id: family_id.to_string(),
            copy_idx: input.copy_idx.clone(),
            is_partner: input.is_partner,
            n_rejected: input.rejected.len(),
            n_aligned: rejected.len(),
            covered_kb: test.covered_kb,
            n_sites: test.n_sites,
            ctl_n: accepted.len(),
            ctl_covered_kb: ctl.covered_kb,
            ctl_n_sites: ctl.n_sites,
            p_uncorrected,
            class,
            med_mismatch,
            med_unaligned,
        });
    }
    out
}

#[cfg(test)]
mod pair_detector_tests {
    use super::*;
    use crate::genome::GenomeIndex;

    #[test]
    fn fewer_than_min_reads_is_skipped_entirely() {
        let genome = GenomeIndex::from_seqs(&[("chrT", &[b'A'; 100])]);
        let mut spans = std::collections::HashMap::new();
        spans.insert("0".to_string(), ("chrT".to_string(), 0u64, 100u64, None));
        let inputs = vec![PairInput {
            copy_idx: "0".to_string(),
            is_partner: false,
            rejected: vec![("r1".to_string(), vec![b'A'; 50])], // only 1, min_reads default is 3
            accepted: vec![],
        }];
        let pairs = detect_missing_copy_pairs("F", &spans, &genome, &inputs, &O3Params::default());
        assert!(pairs.is_empty());
    }

    #[test]
    fn missing_copy_span_is_skipped_not_errored() {
        let genome = GenomeIndex::from_seqs(&[("chrT", &[b'A'; 100])]);
        let spans = std::collections::HashMap::new(); // no span for "0"
        let inputs = vec![PairInput {
            copy_idx: "0".to_string(),
            is_partner: false,
            rejected: vec![
                ("r1".to_string(), vec![b'A'; 50]),
                ("r2".to_string(), vec![b'A'; 50]),
                ("r3".to_string(), vec![b'A'; 50]),
            ],
            accepted: vec![],
        }];
        let pairs = detect_missing_copy_pairs("F", &spans, &genome, &inputs, &O3Params::default());
        assert!(pairs.is_empty());
    }

    #[test]
    fn locus_extent_wider_than_the_copy_span_wins() {
        // Locus extent (40,160) strictly contains the copy span (50,120): the window must widen to the
        // locus, matching `bench/o3_flag_pass.py`'s `(min(locus_start, start), max(locus_end, end))`.
        assert_eq!(locus_or_padded_window(50, 120, Some((40, 160)), 999), (40, 160));
    }

    #[test]
    fn locus_extent_narrower_than_the_copy_span_does_not_shrink_it() {
        // Locus extent (60,110) sits INSIDE the copy span (50,120): min/max never shrinks below the
        // copy's own span, so the window collapses to the bare span, not the narrower locus.
        assert_eq!(locus_or_padded_window(50, 120, Some((60, 110)), 999), (50, 120));
    }

    #[test]
    fn no_locus_falls_back_to_padding_by_the_longest_rejected_read() {
        // No catalog locus extent recorded (`None`): pad the bare copy span by the longest rejected
        // read's length on each side.
        assert_eq!(locus_or_padded_window(1000, 2000, None, 250), (750, 2250));
        // saturating_sub must not underflow when the pad exceeds the start coordinate.
        assert_eq!(locus_or_padded_window(100, 200, None, 500), (0, 700));
    }

    #[test]
    fn min_reads_guard_below_threshold_returns_empty() {
        // Smoke test for detect_missing_copy_pairs when input is below the min_reads threshold.
        // Verifies the function returns early without crashing when given insufficient rejected reads.
        // The Some(locus) destructure and locus_or_padded_window call logic are covered separately
        // by the three locus_or_padded_window unit tests directly above.
        let genome = GenomeIndex::from_seqs(&[("chrT", &[b'A'; 300])]);
        let mut spans = std::collections::HashMap::new();
        // locus extent (0,300) is wider than the bare copy span (100,200).
        spans.insert("0".to_string(), ("chrT".to_string(), 100u64, 200u64, Some((0u64, 300u64))));
        let inputs = vec![PairInput {
            copy_idx: "0".to_string(),
            is_partner: false,
            rejected: vec![("r1".to_string(), vec![b'A'; 50])], // 1 < default min_reads (3)
            accepted: vec![],
        }];
        let pairs = detect_missing_copy_pairs("F", &spans, &genome, &inputs, &O3Params::default());
        assert!(pairs.is_empty(), "below min_reads, skipped before realign_batch as usual");
    }
}
