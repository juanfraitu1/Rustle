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

    /// Assign each stored alignment to the Layer-1 locus (bundle index) it overlaps.
    /// `locus_spans[i] = (chrom, start, end)` of bundle `i`. On overlap ties the
    /// lowest bundle index wins (deterministic).
    pub fn assign_loci(&mut self, locus_spans: &[(String, u64, u64)]) {
        // Build a per-chrom list of (start, end, idx) for overlap lookup.
        let mut by_chrom: DetHashMap<&str, Vec<(u64, u64, usize)>> = DetHashMap::default();
        for (i, (c, s, e)) in locus_spans.iter().enumerate() {
            by_chrom.entry(c.as_str()).or_default().push((*s, *e, i));
        }
        for v in by_chrom.values_mut() {
            v.sort_unstable();
        }
        for a in self.alignments.iter_mut() {
            a.locus = None;
            if let Some(spans) = by_chrom.get(a.chrom.as_str()) {
                // pick the locus with maximal overlap (deterministic: lowest idx on tie)
                let mut best: Option<(u64, usize)> = None;
                for (s, e, i) in spans {
                    // spans are sorted by start; once a locus starts at/after this
                    // alignment's end, no later locus can overlap → stop scanning.
                    if *s >= a.ref_end {
                        break;
                    }
                    if a.ref_start < *e {
                        let ov = a.ref_end.min(*e).saturating_sub(a.ref_start.max(*s));
                        match best {
                            Some((bov, bi)) if ov < bov || (ov == bov && *i >= bi) => {}
                            _ => best = Some((ov, *i)),
                        }
                    }
                }
                a.locus = best.map(|(_, i)| i);
            }
        }
    }

    /// Cross-map links for family discovery. `primary_locus[read_name_hash] = locus`
    /// is the Layer-1 locus the read's PRIMARY (in `bundle.reads`) lives in.
    /// A link (a,b) with a<b is emitted whenever a read's primary is in locus `a`
    /// and one of its secondaries is assigned to a DISTINCT locus `b` (or vice-versa).
    /// Returns sorted, deduplicated `((a,b), count)` pairs.
    pub fn cross_map_links(
        &self,
        primary_locus: &DetHashMap<u64, usize>,
    ) -> Vec<((usize, usize), u32)> {
        let mut counts: DetHashMap<(usize, usize), u32> = DetHashMap::default();
        for (&h, idxs) in self.by_read.iter() {
            let Some(&pl) = primary_locus.get(&h) else { continue };
            // distinct secondary loci for this read
            let mut sec_loci: DetHashSet<usize> = DetHashSet::default();
            for &ai in idxs {
                if let Some(sl) = self.alignments[ai].locus {
                    if sl != pl {
                        sec_loci.insert(sl);
                    }
                }
            }
            for sl in sec_loci {
                let key = if pl < sl { (pl, sl) } else { (sl, pl) };
                *counts.entry(key).or_default() += 1;
            }
        }
        let mut out: Vec<((usize, usize), u32)> = counts.into_iter().collect();
        out.sort_unstable();
        out
    }

    /// All secondaries assigned to `locus`, in capture order (deterministic).
    pub fn secondaries_for_locus(&self, locus: usize) -> Vec<&SecondaryAlignment> {
        self.alignments
            .iter()
            .filter(|a| a.locus == Some(locus))
            .collect()
    }

    /// Drop every alignment whose locus is not in `keep`. Open-decision (1):
    /// prune the index to family-candidate loci after discovery's first pass.
    /// Returns the number of alignments dropped (caller logs it; never silent).
    pub fn prune_to_loci(&mut self, keep: &DetHashSet<usize>) -> usize {
        let before = self.alignments.len();
        let kept: Vec<SecondaryAlignment> = self
            .alignments
            .drain(..)
            .filter(|a| a.locus.map(|l| keep.contains(&l)).unwrap_or(false))
            .collect();
        self.rebuild(kept);
        before - self.alignments.len()
    }

    /// Cap the number of secondaries kept per locus at `cap`, keeping the first
    /// `cap` in capture order. Open-decision (1): logged drop, no silent truncation.
    /// Alignments with no locus are retained (they feed C5 all-secondary detection).
    /// Returns total dropped across all loci.
    pub fn cap_per_locus(&mut self, cap: usize) -> usize {
        let before = self.alignments.len();
        let mut seen: DetHashMap<usize, usize> = DetHashMap::default();
        let kept: Vec<SecondaryAlignment> = self
            .alignments
            .drain(..)
            .filter(|a| match a.locus {
                Some(l) => {
                    let c = seen.entry(l).or_default();
                    if *c < cap {
                        *c += 1;
                        true
                    } else {
                        false
                    }
                }
                None => true,
            })
            .collect();
        self.rebuild(kept);
        before - self.alignments.len()
    }

    /// Rebuild `by_read` after a filtering pass.
    fn rebuild(&mut self, kept: Vec<SecondaryAlignment>) {
        self.alignments = kept;
        self.by_read.clear();
        for (i, a) in self.alignments.iter().enumerate() {
            self.by_read.entry(a.read_name_hash).or_default().push(i);
        }
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

    fn sa_loc(hash: u64, start: u64, end: u64, locus: Option<usize>) -> SecondaryAlignment {
        let mut a = sa(hash, start, end);
        a.locus = locus;
        a
    }

    #[test]
    fn cross_map_links_pair_distinct_loci() {
        // Read 7: primary-locus link encoded via two placements in loci 0 and 3.
        // Read 9: only one locus → no link.
        let mut idx = SecondaryIndex::new();
        idx.push(sa_loc(7, 100, 200, Some(0)));
        idx.push(sa_loc(7, 9000, 9100, Some(3)));
        idx.push(sa_loc(9, 100, 200, Some(0)));
        // primary loci for the two reads: read 7 primary in locus 0, read 9 in 0.
        let mut primary_locus: DetHashMap<u64, usize> = DetHashMap::default();
        primary_locus.insert(7, 0);
        primary_locus.insert(9, 0);
        let links = idx.cross_map_links(&primary_locus);
        // read 7: primary 0, secondary in 3 → link (0,3). read 9: secondary 0 == primary → none.
        assert_eq!(links, vec![((0usize, 3usize), 1u32)]);
    }

    #[test]
    fn per_locus_view_returns_overlapping_secondaries() {
        let mut idx = SecondaryIndex::new();
        idx.push(sa_loc(7, 100, 200, Some(5)));
        idx.push(sa_loc(8, 300, 400, Some(5)));
        idx.push(sa_loc(9, 100, 200, Some(2)));
        let v = idx.secondaries_for_locus(5);
        assert_eq!(v.len(), 2, "two secondaries assigned to locus 5");
        assert_eq!(v[0].read_name_hash, 7);
        assert_eq!(v[1].read_name_hash, 8);
    }

    #[test]
    fn prune_drops_alignments_outside_candidate_loci() {
        let mut idx = SecondaryIndex::new();
        idx.push(sa_loc(7, 100, 200, Some(0)));
        idx.push(sa_loc(8, 300, 400, Some(1)));
        idx.push(sa_loc(9, 100, 200, None)); // no locus → pruned
        let mut keep: DetHashSet<usize> = DetHashSet::default();
        keep.insert(0);
        let dropped = idx.prune_to_loci(&keep);
        assert_eq!(dropped, 2, "locus-1 and no-locus alignments dropped");
        assert_eq!(idx.len(), 1, "only the locus-0 alignment survives");
        assert_eq!(idx.alignments()[0].read_name_hash, 7);
    }

    #[test]
    fn cap_per_locus_logs_overflow_count() {
        let mut idx = SecondaryIndex::new();
        for h in 0..10u64 {
            idx.push(sa_loc(h, 100, 200, Some(0)));
        }
        let dropped = idx.cap_per_locus(4);
        assert_eq!(dropped, 6, "kept 4 of 10 at locus 0; 6 dropped (logged, not silent)");
        assert_eq!(idx.secondaries_for_locus(0).len(), 4);
    }

    #[test]
    fn assign_loci_picks_max_overlap_and_breaks_ties_low_index() {
        // Three loci on chrT; one alignment off-chrom; one with equal overlap to
        // two loci (lowest index wins); one with no overlap (→ None).
        let mut idx = SecondaryIndex::new();
        idx.push(sa(7, 100, 200)); // overlaps locus 0 (0..150) and locus 1 (150..300)? see spans
        idx.push(sa(8, 5000, 5100)); // overlaps nothing → None
        let mut off = sa(9, 100, 200);
        off.chrom = "chrOther".to_string();
        idx.push(off); // different chrom → None
        // locus 0: 90..160 (overlap with 7 = 60), locus 1: 140..210 (overlap with 7 = 60) → tie → idx 0
        let spans = vec![
            ("chrT".to_string(), 90u64, 160u64),
            ("chrT".to_string(), 140u64, 210u64),
        ];
        idx.assign_loci(&spans);
        assert_eq!(idx.alignments()[0].locus, Some(0), "equal overlap → lowest index wins");
        assert_eq!(idx.alignments()[1].locus, None, "no overlapping locus → None");
        assert_eq!(idx.alignments()[2].locus, None, "different chromosome → None");
    }
}
