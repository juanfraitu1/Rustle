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
    /// Placement strand (`'+'`/`'-'`), from the read's alignment strand. A recovered
    /// all-secondary new copy needs a real strand: hardcoding `'+'` is wrong for
    /// spliced transcripts and breaks gffcompare strand-matching downstream.
    pub strand: char,
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

    /// read_name_hash -> number of placements (alignments) this read has in the
    /// side-index. Used by M5's repeat-spread filter: a read placed many ways is a
    /// homology-shadow/repeat artifact, not a genuine new-copy molecule.
    pub fn read_placement_counts(&self) -> DetHashMap<u64, usize> {
        self.by_read.iter().map(|(h, idxs)| (*h, idxs.len())).collect()
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

    /// Cluster the unassigned (locus == None) alignments — those overlapping NO
    /// Layer-1 bundle, i.e. genomic regions Layer 1 dropped because every read there
    /// is a cross-map secondary — into candidate new-copy regions by genomic
    /// proximity. Single-linkage on [ref_start, ref_end] per chromosome: an alignment
    /// joins the running cluster when same chrom AND ref_start <= cluster_end + max_gap
    /// (cluster_end = max ref_end so far). Deterministic. Returns one Vec of
    /// alignment refs per cluster, in genomic order. NO support/min-reads filter here
    /// (the emission stage gates on distinct-read support); this is pure clustering.
    ///
    /// Must run BEFORE `prune_to_loci`, which deletes the locus == None alignments
    /// this method consumes.
    pub fn all_secondary_regions(&self, max_gap: u64) -> Vec<Vec<&SecondaryAlignment>> {
        // 1. Collect refs to the unassigned alignments, paired with capture index.
        //    The capture index is the final sort tiebreak: it gives a TOTAL order so
        //    the clustering is independent of map/iteration nondeterminism — this is
        //    load-bearing for byte-identical output.
        let mut unassigned: Vec<(usize, &SecondaryAlignment)> = self
            .alignments
            .iter()
            .enumerate()
            .filter(|(_, a)| a.locus.is_none())
            .collect();
        // 2. Total-order sort: (chrom, start, end, hash, capture_index).
        unassigned.sort_unstable_by(|(ia, a), (ib, b)| {
            a.chrom
                .cmp(&b.chrom)
                .then(a.ref_start.cmp(&b.ref_start))
                .then(a.ref_end.cmp(&b.ref_end))
                .then(a.read_name_hash.cmp(&b.read_name_hash))
                .then(ia.cmp(ib))
        });
        // 3. Single-linkage cluster: same chrom AND start within cluster_end + max_gap.
        let mut clusters: Vec<Vec<&SecondaryAlignment>> = Vec::new();
        let mut cluster_chrom: &str = "";
        let mut cluster_end: u64 = 0;
        for (_, a) in unassigned {
            let joins = !clusters.is_empty()
                && a.chrom == cluster_chrom
                && a.ref_start <= cluster_end.saturating_add(max_gap);
            if joins {
                cluster_end = cluster_end.max(a.ref_end);
                clusters.last_mut().unwrap().push(a);
            } else {
                cluster_chrom = a.chrom.as_str();
                cluster_end = a.ref_end;
                clusters.push(vec![a]);
            }
        }
        clusters
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

/// Re-export the bundle-side collector so tests + the pipeline have one path.
pub use crate::bundle::collect_secondary_index_from_bam;

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
            strand: '+',
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

    #[test]
    fn collector_captures_secondaries_not_primaries() {
        // Fixture: one primary + one secondary (FLAG 0x100), same QNAME on chrT.
        // The collector must store ONLY the secondary. config.long_reads true so
        // the gate-mirroring predicate admits the record.
        let mut config = crate::types::RunConfig::default();
        config.long_reads = true;
        let idx = super::collect_secondary_index_from_bam(
            "bench/fixtures/mini_secondary.bam",
            Some("chrT"),
            &config,
        )
        .expect("collect from fixture bam");
        assert_eq!(idx.len(), 1, "exactly one secondary captured");
        assert!(!idx.alignments()[0].is_supplementary, "it is a SECONDARY");
        assert_eq!(idx.n_reads(), 1, "one distinct read");
        assert_eq!(idx.alignments()[0].nm, 1, "NM tag parsed via the shared helper");
        assert_eq!(
            idx.alignments()[0].read_name_hash,
            crate::vg::fnv1a64(b"read_xmap")
        );
    }

    #[test]
    fn collector_intron_chain_matches_bundle_parser() {
        // Prove intron-chain derivation uses the shared exon parser
        // (bam::exons_from_cigar). "10M100N10M" at 0-based 100 → one intron (110,210).
        use noodles_sam::alignment::record::cigar::{op::Kind, Op};
        use noodles_sam::alignment::record_buf::Cigar;
        let cigar: Cigar = vec![Op::new(Kind::Match, 10), Op::new(Kind::Skip, 100), Op::new(Kind::Match, 10)]
            .into_iter()
            .collect();
        let exons = crate::bam::exons_from_cigar(100, &cigar).unwrap();
        let introns: Vec<(u64, u64)> = exons
            .windows(2)
            .map(|w| (w[0].1, w[1].0))
            .collect();
        assert_eq!(introns, vec![(110, 210)], "intron chain from shared exon parser");
    }

    // ── M5.0 / M5.1: strand field + all-secondary-region clustering ──────────

    /// Lock in that the `strand` field exists, is `pub`, and round-trips. A
    /// recovered all-secondary copy on the minus strand must keep `'-'`.
    #[test]
    fn secondary_alignment_carries_strand() {
        let mut a = sa(7, 100, 200);
        a.strand = '-';
        assert_eq!(a.strand, '-', "strand field is accessible and preserved");
        // the default-`'+'` helper path also works
        assert_eq!(sa(8, 0, 10).strand, '+');
    }

    #[test]
    fn all_secondary_regions_clusters_unassigned() {
        // Three overlapping locus=None spans + one far away + one ASSIGNED (Some(0))
        // sitting inside the overlapping cluster. The assigned one must be excluded.
        let mut idx = SecondaryIndex::new();
        idx.push(sa_loc(1, 1000, 1200, None));
        idx.push(sa_loc(2, 1100, 1300, None));
        idx.push(sa_loc(3, 1150, 1350, None));
        idx.push(sa_loc(4, 5000, 5200, None)); // far → own cluster
        idx.push(sa_loc(5, 1120, 1320, Some(0))); // assigned → excluded entirely
        let clusters = idx.all_secondary_regions(200);
        assert_eq!(clusters.len(), 2, "two clusters: the overlap group and the far one");
        assert_eq!(clusters[0].len(), 3, "the three overlapping spans cluster together");
        assert_eq!(clusters[1].len(), 1, "the far span is its own cluster");
        // the Some(0) alignment appears in NO cluster
        let total: usize = clusters.iter().map(|c| c.len()).sum();
        assert_eq!(total, 4, "the assigned (Some(0)) alignment is excluded");
        // the size-3 cluster carries exactly the three overlapping spans (genomic order)
        let spans: Vec<(u64, u64)> = clusters[0].iter().map(|a| (a.ref_start, a.ref_end)).collect();
        assert_eq!(spans, vec![(1000, 1200), (1100, 1300), (1150, 1350)]);
    }

    #[test]
    fn all_secondary_regions_respects_max_gap() {
        // Two non-overlapping locus=None spans with a 300bp gap between them.
        let mut idx = SecondaryIndex::new();
        idx.push(sa_loc(1, 1000, 1100, None));
        idx.push(sa_loc(2, 1400, 1500, None));
        // gap (1400 - 1100 = 300) exceeds max_gap=200 → separate clusters
        assert_eq!(idx.all_secondary_regions(200).len(), 2, "gap > max_gap splits");
        // max_gap=400 bridges the gap → single cluster
        assert_eq!(idx.all_secondary_regions(400).len(), 1, "gap <= max_gap merges");
    }

    #[test]
    fn all_secondary_regions_separates_chromosomes() {
        // Identical coordinates but different chromosomes must never merge.
        let mut idx = SecondaryIndex::new();
        idx.push(sa_loc(1, 1000, 1200, None));
        let mut other = sa_loc(2, 1000, 1200, None);
        other.chrom = "chrU".to_string();
        idx.push(other);
        assert_eq!(idx.all_secondary_regions(1000).len(), 2, "different chrom → 2 clusters");
    }

    #[test]
    fn all_secondary_regions_empty_when_all_assigned() {
        let mut idx = SecondaryIndex::new();
        idx.push(sa_loc(1, 1000, 1200, Some(0)));
        idx.push(sa_loc(2, 2000, 2200, Some(1)));
        assert!(idx.all_secondary_regions(200).is_empty(), "no None alignments → no clusters");
    }

    #[test]
    fn all_secondary_regions_deterministic() {
        // Same alignment set pushed in two DIFFERENT orders must cluster identically.
        // The total-order sort (chrom, start, end, hash, capture_index) makes the
        // result independent of insertion / map-iteration order — load-bearing for
        // byte-identical output.
        let items = [
            (1u64, 1000u64, 1200u64),
            (2, 1100, 1300),
            (3, 1150, 1350),
            (4, 5000, 5200),
            (5, 5100, 5300),
        ];
        let mut a = SecondaryIndex::new();
        for &(h, s, e) in items.iter() {
            a.push(sa_loc(h, s, e, None));
        }
        let mut b = SecondaryIndex::new();
        for &(h, s, e) in items.iter().rev() {
            b.push(sa_loc(h, s, e, None));
        }
        let ca = a.all_secondary_regions(200);
        let cb = b.all_secondary_regions(200);
        assert_eq!(ca.len(), cb.len(), "same cluster count regardless of push order");
        let key = |cs: &Vec<Vec<&SecondaryAlignment>>| -> Vec<Vec<(u64, u64, u64)>> {
            cs.iter()
                .map(|c| c.iter().map(|a| (a.ref_start, a.ref_end, a.read_name_hash)).collect())
                .collect()
        };
        assert_eq!(key(&ca), key(&cb), "each cluster's sorted spans are identical");
    }

    #[test]
    fn read_placement_counts_reports_per_read() {
        // M5 repeat-spread filter input: a read with 3 placements + a read with 1.
        // The map must report the per-read placement count so a read placed many
        // ways (a homology-shadow/repeat artifact) can be excluded downstream.
        let mut idx = SecondaryIndex::new();
        idx.push(sa(101, 100, 200));
        idx.push(sa(101, 5000, 5100));
        idx.push(sa(101, 9000, 9100));
        idx.push(sa(202, 100, 200));
        let counts = idx.read_placement_counts();
        assert_eq!(counts.get(&101).copied(), Some(3), "read 101 placed 3 ways");
        assert_eq!(counts.get(&202).copied(), Some(1), "read 202 placed once");
        assert_eq!(counts.len(), 2, "exactly two distinct reads counted");
    }

    #[test]
    fn collector_derives_intron_chain_for_spliced_secondary() {
        // End-to-end: a SPLICED secondary (10M100N10M @ 1-based 101) flows through
        // the collector and yields the intron chain [(110, 210)] — proving the
        // collector's intron derivation works on a real record, not just the
        // exon-parser unit. Fixture built out-of-band with samtools.
        let mut config = crate::types::RunConfig::default();
        config.long_reads = true;
        let idx = super::collect_secondary_index_from_bam(
            "bench/fixtures/mini_secondary_spliced.bam",
            Some("chrT"),
            &config,
        )
        .expect("collect from spliced fixture bam");
        assert_eq!(idx.len(), 1, "one spliced secondary captured");
        assert_eq!(
            idx.alignments()[0].introns,
            vec![(110, 210)],
            "spliced secondary's intron chain derived through the collector"
        );
    }
}
