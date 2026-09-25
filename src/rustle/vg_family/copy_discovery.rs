//! Discovery of candidate gene-family copies from read alignment ties.
//!
//! **STATUS:** OPT-IN
//!
//! Reached from `copy_assign` behind `--discover-copies` (`src/bin/copy_assign.rs`, `default_value_t =
//! false`): `discover_copies_for_family` calls [`cluster_tie_partners`] once per family inside that
//! flag's own `if args.discover_copies {` block. REPORT ONLY -- the result is written to
//! `<out>.discovered_copies.tsv` and never mutates the input catalog or this run's own assignments.
//!
//! ⚠ A [`DiscoveredCopy`] row is a candidate ALIGNED-BLOCK CLUSTER, not necessarily a whole candidate
//! copy: clustering is per-block (see [`aligned_blocks`]), so a genuine multi-exon copy supported by only
//! a few reads can fragment into several rows, one per exon, each independently only needing its own
//! `TIE_PARTNER_MIN_SUPPORT` reads -- more chances for a coincidental cluster than a per-placement scheme
//! would give. Not detected or merged here (docs/o1_ledger.md §6l6 addendum, "Named limitation"); a
//! follow-up would regroup blocks sharing a supporting placement into one candidate with its own
//! `exon_blocks` column, mirroring `catalog_input::exon_blocks_str`.

use std::collections::{HashMap, HashSet};
use crate::vg_family::copy_split::AlignedRead;
use crate::vg_family::denovo_assemble::BamRead;

pub const TIE_PARTNER_MERGE_DISTANCE_BP: u64 = 500;
pub const TIE_PARTNER_MIN_SUPPORT: usize = 2;

#[derive(Clone, Debug, PartialEq)]
pub struct DiscoveredCopy {
    pub family_id: String,
    pub chrom: String,
    pub start: u64,
    pub end: u64,
    /// MAJORITY strand across the placements supporting this cluster (`false` -> `+`, `true` -> `-` from
    /// each record's own SAM FLAG 0x10), an exact tie defaulting to `+` -- the rule the design doc's own
    /// "Open Question" section resolved, matching `majority_read_strand`'s existing convention
    /// (`denovo_pipeline.rs`) and `build_footprint_seq`'s documented `+`-on-tie placeholder caveat.
    pub strand: char,
    pub n_supporting_reads: usize,
    pub read_names: Vec<String>,
    pub nearest_copy_tid: String,
    /// `None` when this family has NO copy on the candidate's own chromosome (reachable for a
    /// cross-chromosome family) -- written as `NA`, never as a `u64::MAX` sentinel.
    pub nearest_copy_distance: Option<u64>,
}

/// One max-AS placement of one AS-tied read, carried as ALIGNED BLOCKS rather than a bounding span.
///
/// ⚠ The blocks, not a `(ref_start, ref_end)` pair, are the load-bearing part. A bounding span counts a
/// spliced-out intron (`N`) as if the read covered it, which made `inside_any_copy` call a placement
/// "inside" a catalog copy that only an intron spans (no aligned base ever landing there) and inflated
/// reported cluster widths by whole intron lengths. This is the same bug class `block_overlap`
/// (`src/bin/copy_assign.rs`) was introduced to fix at the O3 truth gate.
#[derive(Clone, Debug, PartialEq)]
pub struct TiePlacement {
    pub chrom: String,
    /// One `(start, end)` per `M`/`=`/`X` run, in reference order. `D`/`N` advance the reference position
    /// without producing a block (exactly `block_overlap`'s own accounting).
    pub blocks: Vec<(u64, u64)>,
    /// SAM FLAG 0x10 of the record this placement came from -- the strand majority vote's input.
    pub reverse: bool,
}

/// Every `M`/`=`/`X` run of an alignment as its own reference interval `[pos, pos + n)`. The copy_assign
/// binary's `block_overlap` (pysam `get_blocks()` semantics) is computed from these blocks.
/// `D` and `N` advance the reference cursor without emitting a block, so a deletion splits a run here too
/// (matching `block_overlap`) -- usually immaterial since alignment `D` runs are short, but a `D` longer
/// than `TIE_PARTNER_MERGE_DISTANCE_BP` would split one placement into two reported clusters, same as a
/// real gap between two placements. Not observed on the real substrate this module was built against.
pub fn aligned_blocks(read: &AlignedRead) -> Vec<(u64, u64)> {
    let mut pos = read.ref_start;
    let mut out = Vec::new();
    for &(op, n) in &read.cigar {
        match op {
            'M' | '=' | 'X' => {
                out.push((pos, pos + n));
                pos += n;
            }
            'D' | 'N' => pos += n,
            _ => {}
        }
    }
    out
}

/// Does ANY aligned block of this placement overlap a copy already in the catalog?
///
/// Per-block, not a bounding-box test: a read whose intron merely SPANS a catalog copy has no aligned
/// base inside it and is NOT "inside" that copy.
fn inside_any_copy(existing_copies: &[(String, u64, u64, String)], p: &TiePlacement) -> bool {
    existing_copies.iter().any(|(c_chrom, c_start, c_end, _)| {
        *c_chrom == p.chrom && p.blocks.iter().any(|(b_start, b_end)| b_start < c_end && b_end > c_start)
    })
}

/// Nearest catalog copy of this family ON THE SAME CHROMOSOME, and its distance in bp (`0` when the
/// candidate overlaps it). `("NA", None)` when the family has no copy on that chromosome at all.
fn nearest_copy(existing_copies: &[(String, u64, u64, String)], chrom: &str, start: u64, end: u64) -> (String, Option<u64>) {
    existing_copies
        .iter()
        .filter(|(c_chrom, ..)| c_chrom == chrom)
        .map(|(_, c_start, c_end, tid)| {
            let d = if end <= *c_start { c_start - end } else if start >= *c_end { start - c_end } else { 0 };
            (tid.clone(), Some(d))
        })
        .min_by_key(|(_, d)| *d)
        .unwrap_or(("NA".to_string(), None))
}

/// Cluster out-of-catalog AS-tie-partner positions into candidate copies for one family.
///
/// `tied_reads`: `(read_name, that read's max-AS placements)` -- ALREADY RESTRICTED BY THE CALLER to the
/// reads this family actually considered (`discover_copies_for_family`, `src/bin/copy_assign.rs`). Passing
/// a whole region's tied reads to every family is the cross-family pooling bug the final whole-branch
/// review caught: it attributed one identical site, with an identical read list, to 3-8 different
/// `family_id`s.
/// `existing_copies`: `(chrom, start, end, tid)` for every copy already catalogued in this family (the
/// caller zips `FamilyAssignment::copy_spans` with `copy_tids`).
/// Positions already inside an existing copy span are defensively re-excluded here (never surfaced), even
/// though the caller is expected to have filtered them out already.
///
/// Pure: no BAM, no catalog types, no I/O.
pub fn cluster_tie_partners(
    tied_reads: &[(String, Vec<TiePlacement>)],
    family_id: &str,
    existing_copies: &[(String, u64, u64, String)],
    merge_distance: u64,
    min_support: usize,
) -> Vec<DiscoveredCopy> {
    // Flatten to one site per ALIGNED BLOCK of every placement that lands outside every existing copy.
    // `pid` identifies the contributing placement so the strand vote counts each placement once, not once
    // per exon block. Exclusion is per-PLACEMENT (a placement with any block inside a catalog copy is
    // dropped whole); clustering is per-BLOCK, so a spliced read's intron never chains two real loci into
    // one candidate.
    struct Site {
        pid: usize,
        name: String,
        chrom: String,
        start: u64,
        end: u64,
        reverse: bool,
    }
    let mut sites: Vec<Site> = Vec::new();
    let mut pid = 0usize;
    for (name, placements) in tied_reads {
        for p in placements {
            if !inside_any_copy(existing_copies, p) {
                for &(start, end) in &p.blocks {
                    sites.push(Site { pid, name: name.clone(), chrom: p.chrom.clone(), start, end, reverse: p.reverse });
                }
            }
            pid += 1;
        }
    }
    // Sort by (chrom, start, end) so overlap/proximity clustering is a single linear pass. `sort_by` is
    // stable, so equal keys keep `tied_reads`' own (already deterministic) order.
    sites.sort_by(|a, b| (a.chrom.as_str(), a.start, a.end).cmp(&(b.chrom.as_str(), b.start, b.end)));

    struct Cluster {
        chrom: String,
        start: u64,
        end: u64,
        read_names: Vec<String>,
        seen_names: HashSet<String>,
        voted: HashSet<usize>,
        n_fwd: usize,
        n_rev: usize,
    }
    let mut clusters: Vec<Cluster> = Vec::new();
    for s in sites {
        let merge = match clusters.last() {
            Some(last) => last.chrom == s.chrom && s.start <= last.end + merge_distance,
            None => false,
        };
        if !merge {
            clusters.push(Cluster {
                chrom: s.chrom.clone(),
                start: s.start,
                end: s.end,
                read_names: Vec::new(),
                seen_names: HashSet::new(),
                voted: HashSet::new(),
                n_fwd: 0,
                n_rev: 0,
            });
        }
        let c = clusters.last_mut().expect("just pushed or matched");
        c.end = c.end.max(s.end);
        // `seen_names` keeps membership O(1) while `read_names` keeps insertion order for the report.
        if c.seen_names.insert(s.name.clone()) {
            c.read_names.push(s.name);
        }
        if c.voted.insert(s.pid) {
            if s.reverse {
                c.n_rev += 1;
            } else {
                c.n_fwd += 1;
            }
        }
    }

    clusters
        .into_iter()
        .filter(|c| c.read_names.len() >= min_support)
        .map(|c| {
            let (nearest_copy_tid, nearest_copy_distance) = nearest_copy(existing_copies, &c.chrom, c.start, c.end);
            DiscoveredCopy {
                family_id: family_id.to_string(),
                chrom: c.chrom,
                start: c.start,
                end: c.end,
                // Majority vote; an exact tie defaults to `+` (design doc's Open Question, resolved).
                strand: if c.n_rev > c.n_fwd { '-' } else { '+' },
                n_supporting_reads: c.read_names.len(),
                read_names: c.read_names,
                nearest_copy_tid,
                nearest_copy_distance,
            }
        })
        .collect()
}

/// Every AS-tied read in the region, with all of its own max-AS placements as aligned blocks.
///
/// A read is AS-tied here iff at least two of its non-supplementary records share its maximum `AS`.
pub fn tie_partner_placements(bam_reads: &[BamRead]) -> Vec<(String, Vec<TiePlacement>)> {
    let mut by_name: HashMap<&str, Vec<&BamRead>> = HashMap::new();
    for br in bam_reads.iter().filter(|b| !b.is_supplementary) {
        by_name.entry(br.name.as_str()).or_default().push(br);
    }
    let mut out = Vec::new();
    for (name, placements) in by_name {
        if placements.len() < 2 {
            continue;
        }
        let max_as = placements.iter().map(|b| b.as_score).max().unwrap();
        let tied: Vec<&BamRead> = placements.iter().copied().filter(|b| b.as_score == max_as).collect();
        if tied.len() >= 2 {
            out.push((
                name.to_string(),
                tied.into_iter()
                    .map(|b| TiePlacement {
                        chrom: b.chrom.clone(),
                        blocks: aligned_blocks(&b.read),
                        reverse: b.reverse,
                    })
                    .collect(),
            ));
        }
    }
    // Sort by read name for deterministic output order (HashMap iteration is randomized).
    out.sort_by(|a, b| a.0.cmp(&b.0));
    out
}

#[cfg(test)]
mod tests {
    use super::*;

    fn one_copy(tid: &str, chrom: &str, start: u64, end: u64) -> Vec<(String, u64, u64, String)> {
        vec![(chrom.to_string(), start, end, tid.to_string())]
    }

    /// A single-block (unspliced) placement on the forward strand.
    fn pl(chrom: &str, start: u64, end: u64) -> TiePlacement {
        TiePlacement { chrom: chrom.to_string(), blocks: vec![(start, end)], reverse: false }
    }

    fn pl_rev(chrom: &str, start: u64, end: u64) -> TiePlacement {
        TiePlacement { chrom: chrom.to_string(), blocks: vec![(start, end)], reverse: true }
    }

    #[test]
    fn defensively_excludes_positions_inside_a_catalog_copy() {
        let existing = one_copy("c0", "chr1", 1000, 2000);
        // this "tied" position sits INSIDE c0's span -- must never surface as a discovery
        let tied = vec![("read1".to_string(), vec![pl("chr1", 1200, 1300)])];
        let out = cluster_tie_partners(&tied, "FAM0", &existing, 500, 2);
        assert!(out.is_empty(), "a position already inside a catalog copy must never be reported");
    }

    #[test]
    fn merges_positions_within_merge_distance_and_respects_min_support() {
        let existing = one_copy("c0", "chr1", 1000, 2000);
        let tied = vec![
            ("read1".to_string(), vec![pl("chr1", 5000, 5100)]),
            ("read2".to_string(), vec![pl("chr1", 5050, 5150)]), // within 500bp of read1's site
        ];
        let out = cluster_tie_partners(&tied, "FAM0", &existing, 500, 2);
        assert_eq!(out.len(), 1, "two nearby out-of-catalog positions with 2 supporting reads = 1 cluster");
        assert_eq!(out[0].n_supporting_reads, 2);
        assert_eq!(out[0].nearest_copy_tid, "c0");
        assert_eq!(out[0].nearest_copy_distance, Some(3000)); // 5000 - 2000
    }

    #[test]
    fn keeps_clusters_separate_beyond_merge_distance() {
        let existing = one_copy("c0", "chr1", 1000, 2000);
        let tied = vec![
            ("read1".to_string(), vec![pl("chr1", 5000, 5100)]),
            ("read2".to_string(), vec![pl("chr1", 5100, 5100)]),
            ("read3".to_string(), vec![pl("chr1", 9000, 9100)]),
            ("read4".to_string(), vec![pl("chr1", 9050, 9150)]),
        ];
        let out = cluster_tie_partners(&tied, "FAM0", &existing, 500, 2);
        assert_eq!(out.len(), 2, "two far-apart pairs must stay two separate clusters");
    }

    #[test]
    fn drops_clusters_below_min_support() {
        let existing = one_copy("c0", "chr1", 1000, 2000);
        let tied = vec![("read1".to_string(), vec![pl("chr1", 5000, 5100)])];
        let out = cluster_tie_partners(&tied, "FAM0", &existing, 500, 2);
        assert!(out.is_empty(), "a single supporting read must not clear min_support=2");
    }

    #[test]
    fn a_copy_inside_a_spliced_out_intron_does_not_exclude_the_placement() {
        // Mirrors `block_overlap_ignores_an_intron_that_merely_spans_the_window`
        // (src/bin/copy_assign.rs) at the containment layer. ref_start=1000, "10M5000N10M": aligned blocks
        // are [1000,1010) and [6010,6020); the intron covers [1010,6010) with no aligned base in it. A
        // catalog copy sitting entirely inside that intron must NOT swallow the placement -- with the old
        // bounding-span test (`ref_start..ref_end` = 1000..6020) it did, silently killing the candidate.
        let spliced = TiePlacement {
            chrom: "chr1".to_string(),
            blocks: aligned_blocks(&AlignedRead {
                ref_start: 1000,
                cigar: vec![('M', 10), ('N', 5000), ('M', 10)],
                seq: vec![],
                qual: vec![],
            }),
            reverse: false,
        };
        assert_eq!(spliced.blocks, vec![(1000, 1010), (6010, 6020)]);
        let copy_in_the_intron = one_copy("c0", "chr1", 3000, 3100);
        let tied = vec![
            ("read1".to_string(), vec![spliced.clone()]),
            ("read2".to_string(), vec![spliced.clone()]),
        ];
        let out = cluster_tie_partners(&tied, "FAM0", &copy_in_the_intron, 500, 2);
        assert_eq!(out.len(), 2, "both aligned blocks survive as their own candidates; neither is 'inside' c0");
        assert_eq!((out[0].start, out[0].end), (1000, 1010));
        assert_eq!((out[1].start, out[1].end), (6010, 6020), "the 5kb intron never chains the two blocks");

        // Contrast: a copy that genuinely overlaps one of the ALIGNED blocks does exclude the placement.
        let copy_on_the_block = one_copy("c0", "chr1", 1005, 1100);
        let out2 = cluster_tie_partners(&tied, "FAM0", &copy_on_the_block, 500, 2);
        assert!(out2.is_empty(), "a real aligned-base overlap must still exclude the whole placement");
    }

    #[test]
    fn strand_is_the_majority_vote_and_a_tie_defaults_to_plus() {
        let existing = one_copy("c0", "chr1", 1000, 2000);
        // 2 forward + 1 reverse -> '+'
        let fwd_majority = vec![
            ("read1".to_string(), vec![pl("chr1", 5000, 5100)]),
            ("read2".to_string(), vec![pl("chr1", 5050, 5150)]),
            ("read3".to_string(), vec![pl_rev("chr1", 5060, 5160)]),
        ];
        let out = cluster_tie_partners(&fwd_majority, "FAM0", &existing, 500, 2);
        assert_eq!(out.len(), 1);
        assert_eq!(out[0].strand, '+', "2 forward vs 1 reverse");

        // 2 reverse + 1 forward -> '-'
        let rev_majority = vec![
            ("read1".to_string(), vec![pl_rev("chr1", 5000, 5100)]),
            ("read2".to_string(), vec![pl_rev("chr1", 5050, 5150)]),
            ("read3".to_string(), vec![pl("chr1", 5060, 5160)]),
        ];
        let out = cluster_tie_partners(&rev_majority, "FAM0", &existing, 500, 2);
        assert_eq!(out.len(), 1);
        assert_eq!(out[0].strand, '-', "2 reverse vs 1 forward");

        // exact tie -> '+' by the documented default
        let tie = vec![
            ("read1".to_string(), vec![pl("chr1", 5000, 5100)]),
            ("read2".to_string(), vec![pl_rev("chr1", 5050, 5150)]),
        ];
        let out = cluster_tie_partners(&tie, "FAM0", &existing, 500, 2);
        assert_eq!(out.len(), 1);
        assert_eq!(out[0].strand, '+', "an exact strand tie defaults to '+'");
    }

    #[test]
    fn a_multi_exon_placement_votes_once_for_strand() {
        // One reverse placement with 3 exon blocks close enough to land in ONE cluster must not outvote
        // two forward single-block reads: the vote is per PLACEMENT, not per block.
        let existing = one_copy("c0", "chr1", 1000, 2000);
        let three_exons = TiePlacement {
            chrom: "chr1".to_string(),
            blocks: vec![(5000, 5050), (5100, 5150), (5200, 5250)],
            reverse: true,
        };
        let tied = vec![
            ("read1".to_string(), vec![pl("chr1", 5000, 5100)]),
            ("read2".to_string(), vec![pl("chr1", 5050, 5150)]),
            ("read3".to_string(), vec![three_exons]),
        ];
        let out = cluster_tie_partners(&tied, "FAM0", &existing, 500, 2);
        assert_eq!(out.len(), 1);
        assert_eq!(out[0].n_supporting_reads, 3, "distinct read NAMES, not blocks");
        assert_eq!(out[0].strand, '+', "2 forward placements outvote 1 reverse placement's 3 blocks");
    }

    #[test]
    fn no_copy_on_this_chromosome_yields_na_not_a_sentinel() {
        // A cross-chromosome family: the candidate lands on chr2, every catalog copy is on chr1.
        let existing = one_copy("c0", "chr1", 1000, 2000);
        let tied = vec![
            ("read1".to_string(), vec![pl("chr2", 5000, 5100)]),
            ("read2".to_string(), vec![pl("chr2", 5050, 5150)]),
        ];
        let out = cluster_tie_partners(&tied, "FAM0", &existing, 500, 2);
        assert_eq!(out.len(), 1);
        assert_eq!(out[0].nearest_copy_tid, "NA");
        assert_eq!(out[0].nearest_copy_distance, None, "no copy on chr2 -> NA, never u64::MAX");
    }

    #[test]
    fn tie_partner_placements_finds_reads_tied_at_their_own_max_as() {
        use crate::vg_family::denovo_assemble::BamRead;
        use crate::vg_family::copy_split::AlignedRead;
        let mk = |name: &str, chrom: &str, start: u64, as_score: i32| BamRead {
            chrom: chrom.into(),
            read: AlignedRead { ref_start: start, cigar: vec![('M', 100)], seq: vec![], qual: vec![] },
            mapq: 0, name: name.into(), as_score, de: 0.0,
            is_supplementary: false, is_secondary: as_score != 200, reverse: false, ts: None,
        };
        let reads = vec![
            mk("tied_read", "chr1", 1000, 200),   // best
            mk("tied_read", "chr1", 5000, 200),   // tied with the above
            mk("tied_read", "chr1", 9000, 150),   // worse, not part of the tie
            mk("unique_read", "chr1", 2000, 300), // only one placement, never tied
        ];
        let out = tie_partner_placements(&reads);
        assert_eq!(out.len(), 1);
        assert_eq!(out[0].0, "tied_read");
        assert_eq!(out[0].1.len(), 2, "only the 2 max-scoring placements, not the 150-scoring one");
        assert_eq!(out[0].1[0].blocks, vec![(1000, 1100)], "placements carry aligned blocks, not a span");
    }

    #[test]
    fn tie_partner_placements_carries_the_strand_flag_and_splits_on_introns() {
        use crate::vg_family::denovo_assemble::BamRead;
        use crate::vg_family::copy_split::AlignedRead;
        let mk = |start: u64, reverse: bool| BamRead {
            chrom: "chr1".into(),
            read: AlignedRead { ref_start: start, cigar: vec![('M', 10), ('N', 500), ('M', 10)], seq: vec![], qual: vec![] },
            mapq: 0, name: "r".into(), as_score: 200, de: 0.0,
            is_supplementary: false, is_secondary: false, reverse, ts: None,
        };
        let out = tie_partner_placements(&[mk(100, false), mk(9000, true)]);
        assert_eq!(out.len(), 1);
        let p = &out[0].1;
        assert_eq!(p.len(), 2);
        assert_eq!(p[0].blocks, vec![(100, 110), (610, 620)], "the intron is not part of any block");
        assert!(!p[0].reverse);
        assert!(p[1].reverse, "FLAG 0x10 is carried through per placement");
    }
}
