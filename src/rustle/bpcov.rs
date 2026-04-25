//! Per-base coverage from bundle reads (bpcov / get_cov_sign).
//! Used for coverage-based node splitting (find_trims_wsign, trimnode_all).
//!
//! **Strand indexing ( 16925):** sno = (int)strand+1 → bpcov[0]=minus, bpcov[1]=unstranded/all, bpcov[2]=plus.
//! get_cov_sign(sno, ...) uses bpcov[sno]. We build plus and minus separately; index 1 (all) = plus + minus when needed.

use crate::coord::len_half_open;
use crate::types::BundleRead;

/// Per-base coverage for one strand. Index 0 = bundle_start.
#[derive(Debug, Clone)]
pub struct Bpcov {
    /// Coverage at each position [0..len) relative to bundle_start.
    pub cov: Vec<f64>,
    pub bundle_start: u64,
    pub bundle_end: u64,
    /// Cached prefix sum for O(1) range queries. Computed lazily.
    prefix_sum: Option<Vec<f64>>,
}

impl Bpcov {
    /// Build directly from a coverage vector (primarily for tests and debugging).
    pub fn from_cov(cov: Vec<f64>, bundle_start: u64, bundle_end: u64) -> Self {
        Self {
            cov,
            bundle_start,
            bundle_end,
            prefix_sum: None,
        }
    }

    /// Build strand-specific per-base coverage from bundle reads.
    /// Uses inferred transcript strand (`BundleRead.strand`), not raw alignment orientation.
    pub fn from_reads(
        reads: &[BundleRead],
        bundle_start: u64,
        bundle_end: u64,
        use_plus_strand: bool,
    ) -> Self {
        let len = len_half_open(bundle_start, bundle_end) as usize;
        let len = len.saturating_add(1).max(1);
        let mut cov = vec![0.0; len];

        // ST-faithful: include unstranded ('.') reads in both strands'
        // bpcov, mirroring ST's `bpcov[1]=all` semantic (rlink.cpp:533-535)
        // where stranded reads also add to the all-strand layer. Empirical:
        // +2 matches / +0.1 Sn / -0.1 Pr on GGO_19 vs the prior strand-
        // exclusive default. Opt out via RUSTLE_BPCOV_EXCLUDE_NEUTRAL=1.
        let include_neutral =
            std::env::var_os("RUSTLE_BPCOV_EXCLUDE_NEUTRAL").is_none();
        for r in reads {
            // get_cov_sign(sno,...) expects strand-specific coverage:
            // - plus bundles use only '+' reads
            // - minus bundles use only '-' reads
            // Neutralized '.' reads should not be counted as minus (asymmetric) because that can
            // mask real boundary drops and prevent longtrim splits.
            if use_plus_strand {
                if r.strand != '+' && !(include_neutral && r.strand == '.') {
                    continue;
                }
            } else if r.strand != '-' && !(include_neutral && r.strand == '.') {
                continue;
            }
            let w = r.weight;
            for &(s, e) in &r.exons {
                let start = s.saturating_sub(bundle_start);
                let end = e.saturating_sub(bundle_start);
                let i0 = (start as usize).min(len);
                let i1 = (end as usize).min(len);
                // Delta encoding: O(1) per exon instead of O(exon_length)
                if i0 < len {
                    cov[i0] += w;
                }
                if i1 < len {
                    cov[i1] -= w;
                }
            }
        }

        // Convert delta-encoded array to actual per-base coverage via prefix sum scan
        // RUSTLE_BPCOV_NEG_CLAMP=1: clamp negative values to 0 after the
        // prefix-sum scan, mirroring ST's negative-clamp at rlink.cpp:18222
        // (`if(bpcov[s][i]<0) bpcov[s][i]=0`). Floating-point ordering of
        // diff-array updates can leave small negatives after summation.
        let neg_clamp = std::env::var_os("RUSTLE_BPCOV_NEG_CLAMP").is_some();
        let mut acc = 0.0f64;
        for v in cov.iter_mut() {
            acc += *v;
            *v = if neg_clamp && acc < 0.0 { 0.0 } else { acc };
        }

        Self {
            cov,
            bundle_start,
            bundle_end,
            prefix_sum: None,
        }
    }

    /// Compute prefix sum for O(1) range queries.
    fn compute_prefix_sum(&mut self) {
        if self.prefix_sum.is_some() {
            return;
        }
        let mut ps = vec![0.0; self.cov.len() + 1];
        for (i, &v) in self.cov.iter().enumerate() {
            ps[i + 1] = ps[i] + v;
        }
        self.prefix_sum = Some(ps);
    }

    /// Get prefix sum slice (computes if needed).
    fn get_prefix_sum(&mut self) -> &[f64] {
        self.compute_prefix_sum();
        self.prefix_sum.as_ref().unwrap()
    }

    /// Sum of coverage in [start_idx, end_idx] (array indices).
    #[inline]
    pub fn get_cov_range(&self, start_idx: usize, end_idx: usize) -> f64 {
        let s = start_idx.min(self.cov.len());
        let e = end_idx.min(self.cov.len());
        if e <= s {
            return 0.0;
        }
        // Use prefix sum if available for O(1) query
        if let Some(ref ps) = self.prefix_sum {
            return ps[e] - ps[s];
        }
        // Fall back to linear sum for small ranges or when prefix not built
        self.cov[s..e].iter().sum()
    }

    /// Sum of coverage in [start_idx, end_idx] with mutable self to enable prefix sum caching.
    #[inline]
    pub fn get_cov_range_cached(&mut self, start_idx: usize, end_idx: usize) -> f64 {
        let s = start_idx.min(self.cov.len());
        let e = end_idx.min(self.cov.len());
        if e <= s {
            return 0.0;
        }
        let ps = self.get_prefix_sum();
        ps[e] - ps[s]
    }

    /// Array index for a genomic position.
    #[inline]
    pub fn idx(&self, pos: u64) -> usize {
        pos.saturating_sub(self.bundle_start) as usize
    }

    /// Add coverage v to [start, end) (cov_edge_add for one strand).
    pub fn add_coverage(&mut self, start: u64, end: u64, v: f64) {
        let i0 = self.idx(start).min(self.cov.len());
        let i1 = self.idx(end).min(self.cov.len());
        if i0 < i1 {
            for i in i0..i1 {
                self.cov[i] += v;
            }
            // Invalidate prefix sum since data changed
            self.prefix_sum = None;
        }
    }

    /// Ensure prefix sum is computed for subsequent O(1) range queries.
    pub fn build_prefix_sum(&mut self) {
        self.compute_prefix_sum();
    }
}

/// Strand index for bpcov (ref: sno = strand+1). 0 = minus, 1 = unstranded/all, 2 = plus.
pub const BPCOV_STRAND_MINUS: usize = 0;
pub const BPCOV_STRAND_ALL: usize = 1;
pub const BPCOV_STRAND_PLUS: usize = 2;

/// Per-base coverage for both strands and optional "all" (bpcov[0], bpcov[1], bpcov[2]).
#[derive(Debug, Clone)]
pub struct BpcovStranded {
    pub minus: Bpcov,
    pub plus: Bpcov,
    /// All/unstranded = minus + plus (bpcov[1] when neutral reads added to both).
    pub all: Option<Vec<f64>>,
}

impl BpcovStranded {
    /// Create an empty BpcovStranded for a given coordinate range.
    /// Used for cross-strand coverage accumulation.
    pub fn empty(bundle_start: u64, bundle_end: u64) -> Self {
        let len = len_half_open(bundle_start, bundle_end) as usize + 1;
        let len = len.max(1);
        Self {
            minus: Bpcov {
                cov: vec![0.0; len],
                bundle_start,
                bundle_end,
                prefix_sum: None,
            },
            plus: Bpcov {
                cov: vec![0.0; len],
                bundle_start,
                bundle_end,
                prefix_sum: None,
            },
            all: Some(vec![0.0; len]),
        }
    }

    /// Build both-strand coverage from reads in a single pass (avoids scanning reads twice).
    pub fn from_reads(reads: &[BundleRead], bundle_start: u64, bundle_end: u64) -> Self {
        let len = len_half_open(bundle_start, bundle_end) as usize;
        let len = len.saturating_add(1).max(1);
        let mut cov_minus = vec![0.0f64; len];
        let mut cov_plus = vec![0.0f64; len];

        for r in reads {
            let w = r.weight;
            let cov = if r.strand == '+' {
                &mut cov_plus
            } else {
                &mut cov_minus
            };
            for &(s, e) in &r.exons {
                let i0 = (s.saturating_sub(bundle_start) as usize).min(len);
                let i1 = (e.saturating_sub(bundle_start) as usize).min(len);
                // Delta encoding: O(1) per exon instead of O(exon_length)
                if i0 < len {
                    cov[i0] += w;
                }
                if i1 < len {
                    cov[i1] -= w;
                }
            }
        }

        // Convert delta-encoded arrays to actual per-base coverage
        {
            let mut acc = 0.0f64;
            for v in cov_minus.iter_mut() {
                acc += *v;
                *v = acc;
            }
        }
        {
            let mut acc = 0.0f64;
            for v in cov_plus.iter_mut() {
                acc += *v;
                *v = acc;
            }
        }

        let mut all = vec![0.0f64; len];
        for i in 0..len {
            all[i] = cov_minus[i] + cov_plus[i];
        }

        Self {
            minus: Bpcov {
                cov: cov_minus,
                bundle_start,
                bundle_end,
                prefix_sum: None,
            },
            plus: Bpcov {
                cov: cov_plus,
                bundle_start,
                bundle_end,
                prefix_sum: None,
            },
            all: Some(all),
        }
    }

    /// get_cov_sign(sno, ...): coverage for strand index 0/1/2.
    pub fn get_cov_range(&self, strand_idx: usize, start_idx: usize, end_idx: usize) -> f64 {
        match strand_idx {
            BPCOV_STRAND_MINUS => self.minus.get_cov_range(start_idx, end_idx),
            BPCOV_STRAND_PLUS => self.plus.get_cov_range(start_idx, end_idx),
            BPCOV_STRAND_ALL => {
                if let Some(ref a) = self.all {
                    let s = start_idx.min(a.len());
                    let e = end_idx.min(a.len());
                    if e <= s {
                        return 0.0;
                    }
                    a[s..e].iter().sum()
                } else {
                    self.minus.get_cov_range(start_idx, end_idx)
                        + self.plus.get_cov_range(start_idx, end_idx)
                }
            }
            _ => 0.0,
        }
    }

    /// Add coverage v to strand strand_idx over [start, end); when strand_idx != BPCOV_STRAND_ALL also add to all (cov_edge_add).
    pub fn add_coverage(&mut self, strand_idx: usize, start: u64, end: u64, v: f64) {
        match strand_idx {
            BPCOV_STRAND_MINUS => {
                self.minus.add_coverage(start, end, v);
                if let Some(ref mut a) = self.all {
                    let i0 = self.minus.idx(start).min(a.len());
                    let i1 = self.minus.idx(end).min(a.len());
                    for i in i0..i1 {
                        a[i] += v;
                    }
                }
            }
            BPCOV_STRAND_PLUS => {
                self.plus.add_coverage(start, end, v);
                if let Some(ref mut a) = self.all {
                    let i0 = self.plus.idx(start).min(a.len());
                    let i1 = self.plus.idx(end).min(a.len());
                    for i in i0..i1 {
                        a[i] += v;
                    }
                }
            }
            BPCOV_STRAND_ALL => {
                if let Some(ref mut a) = self.all {
                    let i0 = self.minus.idx(start).min(a.len());
                    let i1 = self.minus.idx(end).min(a.len());
                    for i in i0..i1 {
                        a[i] += v;
                    }
                }
            }
            _ => {}
        }
    }

    /// Add another stranded coverage grid for the same `(bundle_start, bundle_end)` interval.
    /// Used when merging parallel partial accumulations for cross-strand `BpcovStranded`.
    pub fn accumulate_equal_extent(&mut self, other: &Self) {
        debug_assert_eq!(self.minus.bundle_start, other.minus.bundle_start);
        debug_assert_eq!(self.minus.bundle_end, other.minus.bundle_end);
        debug_assert_eq!(self.minus.cov.len(), other.minus.cov.len());

        self.minus.prefix_sum = None;
        self.plus.prefix_sum = None;
        for (a, b) in self.minus.cov.iter_mut().zip(&other.minus.cov) {
            *a += b;
        }
        for (a, b) in self.plus.cov.iter_mut().zip(&other.plus.cov) {
            *a += b;
        }
        if let (Some(ar), Some(br)) = (&mut self.all, other.all.as_ref()) {
            for (a, b) in ar.iter_mut().zip(br) {
                *a += b;
            }
        }
    }
}
