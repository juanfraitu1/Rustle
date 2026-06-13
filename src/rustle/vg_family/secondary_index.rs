//! C1 — Secondary side-index for Layer-2 family variation graph.
//!
//! Phase 1 (commit 664919c) stopped VG from ingesting secondary/supplementary
//! alignments into `bundle.reads` (they inflated per-gene-graph transfrags and
//! caused VG to DROP a whole baseline region — see
//! `project_vg_drops_baseline_region_rootcause`). Phase 1 only *dropped* them.
//!
//! Layer 2 needs those secondaries back as evidence — but NEVER in bundles.
//! This side-index is the *only* place secondaries live. It restores the
//! family-discovery signal (`build_multimap_index` measured 2125 → 313 reads
//! when secondaries left bundles) and feeds graph amendment, without touching
//! Layer-1 bundling at all.
//!
//! Determinism: all maps are `DetHashMap`/`DetHashSet` (FxHash, no seed);
//! every iteration that feeds output is sorted.

use crate::types::{DetHashMap, DetHashSet};

/// One secondary/supplementary alignment that Layer 1 dropped from `bundle.reads`.
/// All coordinates are 0-based half-open, identical to `BundleRead` conventions.
#[derive(Debug, Clone)]
pub struct SecondaryAlignment {
    /// Hash of the read QNAME (matches `BundleRead::read_name_hash`); links a
    /// secondary back to the primary it shadows. MUST be produced by
    /// `crate::vg::fnv1a64(name_bytes)` so it equals the primary's hash.
    pub read_name_hash: u64,
    /// Chromosome name this placement is on.
    pub chrom: String,
    /// Alignment span on this placement (0-based, half-open).
    pub ref_start: u64,
    pub ref_end: u64,
    /// Intron chain on this placement: (donor_site, acceptor_site) per junction.
    pub introns: Vec<(u64, u64)>,
    /// Edit distance (NM tag) for this placement — used for PSV / decisive evidence.
    pub nm: u32,
    /// `true` for supplementary alignments, `false` for secondary alignments.
    pub is_supplementary: bool,
    /// Layer-1 locus (bundle index) this placement overlaps, filled in after
    /// bundling by `assign_loci`. `None` until then (or if it overlaps no locus).
    pub locus: Option<usize>,
}

/// Chromosome-global store of dropped secondary/supplementary alignments.
///
/// Built once per chromosome by `collect_secondary_index_from_bam`. Cross-map
/// links (a primary in locus A whose secondary lands in locus B) are derived
/// from it for family discovery; per-locus views are derived for amendment.
#[derive(Debug, Default)]
pub struct SecondaryIndex {
    /// All captured secondary/supplementary alignments, in capture order.
    alignments: Vec<SecondaryAlignment>,
    /// read_name_hash → indices into `alignments` (one read can have many).
    by_read: DetHashMap<u64, Vec<usize>>,
}

impl SecondaryIndex {
    pub fn new() -> Self {
        SecondaryIndex {
            alignments: Vec::new(),
            by_read: DetHashMap::default(),
        }
    }

    /// Record one dropped secondary/supplementary alignment.
    pub fn push(&mut self, a: SecondaryAlignment) {
        let idx = self.alignments.len();
        self.by_read.entry(a.read_name_hash).or_default().push(idx);
        self.alignments.push(a);
    }

    /// Total number of stored alignments.
    pub fn len(&self) -> usize {
        self.alignments.len()
    }

    pub fn is_empty(&self) -> bool {
        self.alignments.is_empty()
    }

    /// Number of distinct reads represented.
    pub fn n_reads(&self) -> usize {
        self.by_read.len()
    }

    /// Read-only access to all alignments.
    pub fn alignments(&self) -> &[SecondaryAlignment] {
        &self.alignments
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn sa(hash: u64, start: u64, end: u64) -> SecondaryAlignment {
        SecondaryAlignment {
            read_name_hash: hash,
            chrom: "chrT".to_string(),
            ref_start: start,
            ref_end: end,
            introns: Vec::new(),
            nm: 0,
            is_supplementary: false,
            locus: None,
        }
    }

    #[test]
    fn empty_index_reports_empty() {
        let idx = SecondaryIndex::new();
        assert!(idx.is_empty(), "fresh index is empty");
        assert_eq!(idx.len(), 0);
        assert_eq!(idx.n_reads(), 0);
    }

    #[test]
    fn push_groups_alignments_by_read() {
        let mut idx = SecondaryIndex::new();
        idx.push(sa(7, 100, 200));
        idx.push(sa(7, 5000, 5100));
        idx.push(sa(9, 100, 200));
        assert_eq!(idx.len(), 3, "three alignments stored");
        assert_eq!(idx.n_reads(), 2, "two distinct reads");
        assert_eq!(idx.alignments()[0].read_name_hash, 7);
        // both placements of read 7 are retained (not collapsed to one)
        assert_eq!(
            idx.alignments()
                .iter()
                .filter(|a| a.read_name_hash == 7)
                .count(),
            2,
            "both alignments for read 7 stored"
        );
    }
}
