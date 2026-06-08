//! Within-locus haplotype phasing for `--vg-phase`. Operates on one copy's
//! read set (a single VG-family bundle). Detects heterozygous sites from the
//! per-read mismatch pileup (B1), phases reads with exact Minimum-Error-
//! Correction (B2), and assigns hp/ps tags. Reference-base-free: alleles are
//! binary (Alt = read carries the site's dominant alt mismatch, Ref = spans
//! the site without it).

use crate::types::BundleRead;

/// A candidate heterozygous site within one copy's read set.
#[derive(Debug, Clone, PartialEq)]
pub struct HetSite {
    pub pos: u64,       // reference position (0-based; matches BundleRead.mismatches)
    pub alt_allele: u8, // dominant alternate base (ASCII A/C/G/T)
    pub n_ref: u32,
    pub n_alt: u32,
}

#[derive(Debug, Clone)]
pub struct PhasingConfig {
    pub min_maf: f64,          // minor-allele-fraction floor (default 0.25)
    pub max_maf: f64,          // ceiling — excludes fixed/ref-artifact sites (default 0.75)
    pub min_allele_reads: u32, // per-allele read floor (default 3)
    pub max_coverage: usize,   // MEC-DP active-read cap (default 15)
    pub ext_hp_frac: f64,      // external-HP-tag precedence threshold (default 0.5)
}

impl Default for PhasingConfig {
    fn default() -> Self {
        PhasingConfig {
            min_maf: 0.25,
            max_maf: 0.75,
            min_allele_reads: 3,
            max_coverage: 15,
            ext_hp_frac: 0.5,
        }
    }
}

/// True if `pos` falls inside any of the read's aligned exon intervals.
fn read_spans(read: &BundleRead, pos: u64) -> bool {
    read.exons.iter().any(|&(s, e)| pos >= s && pos < e)
}

/// Detect biallelic heterozygous sites from the mismatch pileup.
pub fn detect_het_sites(reads: &[BundleRead], cfg: &PhasingConfig) -> Vec<HetSite> {
    use std::collections::HashMap;
    // pos -> (alt_base -> count)
    let mut alt_counts: HashMap<u64, HashMap<u8, u32>> = HashMap::new();
    for r in reads {
        for &(pos, base) in &r.mismatches {
            *alt_counts.entry(pos).or_default().entry(base).or_default() += 1;
        }
    }
    let mut sites: Vec<HetSite> = Vec::new();
    for (&pos, bases) in &alt_counts {
        // Dominant alt base at this position (deterministic tie-break by base value).
        let (&alt_base, &n_alt) = match bases
            .iter()
            .max_by(|a, b| a.1.cmp(b.1).then(b.0.cmp(a.0)))
        {
            Some(v) => v,
            None => continue,
        };
        let coverage = reads.iter().filter(|r| read_spans(r, pos)).count() as u32;
        if coverage < n_alt {
            continue; // defensive
        }
        let n_ref = coverage - n_alt;
        let denom = (n_ref + n_alt) as f64;
        if denom == 0.0 {
            continue;
        }
        let maf = n_alt as f64 / denom;
        if n_alt >= cfg.min_allele_reads
            && n_ref >= cfg.min_allele_reads
            && maf >= cfg.min_maf
            && maf <= cfg.max_maf
        {
            sites.push(HetSite { pos, alt_allele: alt_base, n_ref, n_alt });
        }
    }
    sites.sort_by(|a, b| a.pos.cmp(&b.pos).then(a.alt_allele.cmp(&b.alt_allele)));
    sites
}

/// A read's allele at each het site: Some(true)=Alt, Some(false)=Ref, None=not covered.
pub type AlleleRow = Vec<Option<bool>>;

/// Build the read × het-site allele matrix (row order matches `reads`).
pub fn allele_matrix(reads: &[BundleRead], sites: &[HetSite]) -> Vec<AlleleRow> {
    reads
        .iter()
        .map(|r| {
            sites
                .iter()
                .map(|s| {
                    if !read_spans(r, s.pos) {
                        None
                    } else {
                        Some(r.mismatches.iter().any(|&(p, b)| p == s.pos && b == s.alt_allele))
                    }
                })
                .collect()
        })
        .collect()
}

/// Cost of a read row against a haplotype allele vector (Hamming over covered sites).
fn row_cost(row: &AlleleRow, hap: &[bool]) -> u32 {
    row.iter()
        .zip(hap)
        .filter_map(|(a, &h)| a.map(|av| if av != h { 1 } else { 0 }))
        .sum()
}

/// Exact MEC over all 2^M haplotype-A allele assignments (M = #sites).
/// Each read greedily joins the cheaper of {hapA, complement(hapA)}.
/// Returns (hapA alleles, side per read [false=A,true=B], total cost).
/// Caller must ensure M is small (<= ~20).
pub fn mec_brute(matrix: &[AlleleRow], n_sites: usize) -> (Vec<bool>, Vec<bool>, u32) {
    let mut best: Option<(Vec<bool>, Vec<bool>, u32)> = None;
    for mask in 0u32..(1u32 << n_sites) {
        let hap_a: Vec<bool> = (0..n_sites).map(|j| (mask >> j) & 1 == 1).collect();
        let hap_b: Vec<bool> = hap_a.iter().map(|&x| !x).collect();
        let mut sides = Vec::with_capacity(matrix.len());
        let mut total = 0u32;
        for row in matrix {
            let ca = row_cost(row, &hap_a);
            let cb = row_cost(row, &hap_b);
            if ca <= cb {
                sides.push(false);
                total += ca;
            } else {
                sides.push(true);
                total += cb;
            }
        }
        if best.as_ref().map_or(true, |b| total < b.2) {
            best = Some((hap_a, sides, total));
        }
    }
    best.unwrap_or((vec![false; n_sites], vec![false; matrix.len()], 0))
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::types::BundleRead;
    use std::sync::Arc;

    // Minimal BundleRead builder: one exon span [start,end), given mismatches.
    fn mk_read(name_hash: u64, start: u64, end: u64, mismatches: Vec<(u64, u8)>) -> BundleRead {
        BundleRead {
            read_uid: name_hash,
            read_name: Arc::from(""),
            read_name_hash: name_hash,
            ref_id: None,
            mate_ref_id: None,
            mate_start: None,
            hi: 0,
            ref_start: start,
            ref_end: end,
            exons: vec![(start, end)],
            junctions: Vec::new(),
            junction_valid: Vec::new(),
            junctions_raw: Vec::new(),
            junctions_del: Vec::new(),
            weight: 1.0,
            is_reverse: false,
            strand: '+',
            has_poly_start: false,
            has_poly_end: false,
            has_poly_start_aligned: false,
            has_poly_start_unaligned: false,
            has_poly_end_aligned: false,
            has_poly_end_unaligned: false,
            unaligned_poly_t: 0,
            unaligned_poly_a: 0,
            has_last_exon_polya: false,
            has_first_exon_polyt: false,
            query_length: None,
            clip_left: 0,
            clip_right: 0,
            nh: 1,
            nm: 0,
            de: None,
            md: None,
            insertion_sites: Vec::new(),
            unitig: false,
            unitig_cov: 0.0,
            read_count_yc: 1.0,
            countfrag_len: 0.0,
            countfrag_num: 0.0,
            junc_mismatch_weight: 0.0,
            pair_idx: Vec::new(),
            pair_count: Vec::new(),
            mapq: 60,
            mismatches,
            seq: Vec::new(),
            hp_tag: None,
            ps_tag: None,
            is_primary_alignment: true,
            em_weight_gap: -1.0,
            em_n_sites: 0,
            em_anchored: false,
            em_ev_decisive: false,
        }
    }

    #[test]
    fn detects_balanced_het() {
        let cfg = PhasingConfig::default();
        let mut reads = Vec::new();
        for i in 0..5 {
            reads.push(mk_read(i, 50, 200, vec![(100, b'A')]));
        }
        for i in 5..10 {
            reads.push(mk_read(i, 50, 200, vec![]));
        }
        let sites = detect_het_sites(&reads, &cfg);
        assert_eq!(sites.len(), 1);
        assert_eq!(sites[0].pos, 100);
        assert_eq!(sites[0].alt_allele, b'A');
        assert_eq!(sites[0].n_alt, 5);
        assert_eq!(sites[0].n_ref, 5);
    }

    #[test]
    fn rejects_sequencing_error_low_maf() {
        let cfg = PhasingConfig::default();
        let mut reads = Vec::new();
        reads.push(mk_read(0, 50, 200, vec![(100, b'A')]));
        for i in 1..20 {
            reads.push(mk_read(i, 50, 200, vec![]));
        }
        assert!(detect_het_sites(&reads, &cfg).is_empty());
    }

    #[test]
    fn rejects_fixed_difference_high_maf() {
        let cfg = PhasingConfig::default();
        let reads: Vec<_> = (0..10).map(|i| mk_read(i, 50, 200, vec![(100, b'A')])).collect();
        assert!(detect_het_sites(&reads, &cfg).is_empty());
    }

    #[test]
    fn mec_brute_two_clean_haplotypes() {
        let matrix = vec![
            vec![Some(true), Some(false)],
            vec![Some(true), Some(false)],
            vec![Some(false), Some(true)],
            vec![Some(false), Some(true)],
        ];
        let (_hap, sides, cost) = mec_brute(&matrix, 2);
        assert_eq!(cost, 0);
        assert_eq!(sides[0], sides[1]);
        assert_eq!(sides[2], sides[3]);
        assert_ne!(sides[0], sides[2]);
    }

    #[test]
    fn mec_brute_tolerates_one_error() {
        let matrix = vec![
            vec![Some(true), Some(true)],
            vec![Some(true), Some(false)],
            vec![Some(false), Some(true)],
            vec![Some(false), Some(true)],
        ];
        let (_hap, _sides, cost) = mec_brute(&matrix, 2);
        assert_eq!(cost, 1);
    }
}
