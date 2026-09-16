//! Discovery of candidate gene-family copies from read alignment ties.
//!
//! **STATUS:** INFRASTRUCTURE

use std::collections::HashMap;
use crate::vg_family::denovo_assemble::BamRead;

pub const TIE_PARTNER_MERGE_DISTANCE_BP: u64 = 500;
pub const TIE_PARTNER_MIN_SUPPORT: usize = 2;

#[derive(Clone, Debug, PartialEq)]
pub struct DiscoveredCopy {
    pub family_id: String,
    pub chrom: String,
    pub start: u64,
    pub end: u64,
    pub n_supporting_reads: usize,
    pub read_names: Vec<String>,
    pub nearest_copy_tid: String,
    pub nearest_copy_distance: u64,
}

fn inside_any_copy(existing_copies: &[(String, u64, u64, String)], chrom: &str, start: u64, end: u64) -> bool {
    existing_copies.iter().any(|(c_chrom, c_start, c_end, _)| c_chrom == chrom && start < *c_end && end > *c_start)
}

fn nearest_copy(existing_copies: &[(String, u64, u64, String)], chrom: &str, start: u64, end: u64) -> (String, u64) {
    existing_copies
        .iter()
        .filter(|(c_chrom, ..)| c_chrom == chrom)
        .map(|(_, c_start, c_end, tid)| {
            let d = if end <= *c_start { c_start - end } else if start >= *c_end { start - c_end } else { 0 };
            (tid.clone(), d)
        })
        .min_by_key(|(_, d)| *d)
        .unwrap_or(("NA".to_string(), u64::MAX))
}

/// Cluster out-of-catalog AS-tie-partner positions into candidate copies for one family.
/// `tied_reads`: (read_name, all of that read's max-AS placements as (chrom, start, end)).
/// `existing_copies`: (chrom, start, end, tid) for every copy already catalogued in this family (caller
/// converts from whatever its own copy-list type is -- `ColocatedFamily.copies: Vec<DenovoTranscript>` at
/// the real call site, Task 3).
/// Positions already inside an existing copy span are defensively re-excluded here (never surfaced),
/// even though the caller (Task 3) is expected to have filtered them out already.
pub fn cluster_tie_partners(
    tied_reads: &[(String, Vec<(String, u64, u64)>)],
    family_id: &str,
    existing_copies: &[(String, u64, u64, String)],
    merge_distance: u64,
    min_support: usize,
) -> Vec<DiscoveredCopy> {
    // Flatten to (read_name, chrom, start, end), keeping only positions outside every existing copy.
    let mut sites: Vec<(String, String, u64, u64)> = Vec::new();
    for (name, placements) in tied_reads {
        for (chrom, start, end) in placements {
            if !inside_any_copy(existing_copies, chrom, *start, *end) {
                sites.push((name.clone(), chrom.clone(), *start, *end));
            }
        }
    }
    // Sort by (chrom, start) so overlap/proximity clustering is a single linear pass.
    sites.sort_by(|a, b| (a.1.clone(), a.2).cmp(&(b.1.clone(), b.2)));

    let mut clusters: Vec<(String, u64, u64, Vec<String>)> = Vec::new(); // (chrom, start, end, read_names)
    for (name, chrom, start, end) in sites {
        if let Some(last) = clusters.last_mut() {
            let (lchrom, _lstart, lend, names) = last;
            if *lchrom == chrom && start <= *lend + merge_distance {
                *lend = (*lend).max(end);
                if !names.contains(&name) {
                    names.push(name);
                }
                continue;
            }
        }
        clusters.push((chrom, start, end, vec![name]));
    }

    clusters
        .into_iter()
        .filter(|(_, _, _, names)| names.len() >= min_support)
        .map(|(chrom, start, end, read_names)| {
            let (nearest_copy_tid, nearest_copy_distance) = nearest_copy(existing_copies, &chrom, start, end);
            DiscoveredCopy {
                family_id: family_id.to_string(),
                chrom,
                start,
                end,
                n_supporting_reads: read_names.len(),
                read_names,
                nearest_copy_tid,
                nearest_copy_distance,
            }
        })
        .collect()
}

fn ref_end(br: &BamRead) -> u64 {
    br.read.ref_start
        + br.read.cigar.iter().filter(|(op, _)| matches!(op, 'M' | '=' | 'X' | 'D' | 'N')).map(|(_, n)| *n).sum::<u64>()
}

pub fn tie_partner_placements(bam_reads: &[BamRead]) -> Vec<(String, Vec<(String, u64, u64)>)> {
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
                tied.into_iter().map(|b| (b.chrom.clone(), b.read.ref_start, ref_end(b))).collect(),
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

    #[test]
    fn defensively_excludes_positions_inside_a_catalog_copy() {
        let existing = one_copy("c0", "chr1", 1000, 2000);
        // this "tied" position sits INSIDE c0's span -- must never surface as a discovery
        let tied = vec![("read1".to_string(), vec![("chr1".to_string(), 1200, 1300)])];
        let out = cluster_tie_partners(&tied, "FAM0", &existing, 500, 2);
        assert!(out.is_empty(), "a position already inside a catalog copy must never be reported");
    }

    #[test]
    fn merges_positions_within_merge_distance_and_respects_min_support() {
        let existing = one_copy("c0", "chr1", 1000, 2000);
        let tied = vec![
            ("read1".to_string(), vec![("chr1".to_string(), 5000, 5100)]),
            ("read2".to_string(), vec![("chr1".to_string(), 5050, 5150)]), // within 500bp of read1's site
        ];
        let out = cluster_tie_partners(&tied, "FAM0", &existing, 500, 2);
        assert_eq!(out.len(), 1, "two nearby out-of-catalog positions with 2 supporting reads = 1 cluster");
        assert_eq!(out[0].n_supporting_reads, 2);
        assert_eq!(out[0].nearest_copy_tid, "c0");
        assert_eq!(out[0].nearest_copy_distance, 3000); // 5000 - 2000
    }

    #[test]
    fn keeps_clusters_separate_beyond_merge_distance() {
        let existing = one_copy("c0", "chr1", 1000, 2000);
        let tied = vec![
            ("read1".to_string(), vec![("chr1".to_string(), 5000, 5100)]),
            ("read2".to_string(), vec![("chr1".to_string(), 5100, 5100)]),
            ("read3".to_string(), vec![("chr1".to_string(), 9000, 9100)]),
            ("read4".to_string(), vec![("chr1".to_string(), 9050, 9150)]),
        ];
        let out = cluster_tie_partners(&tied, "FAM0", &existing, 500, 2);
        assert_eq!(out.len(), 2, "two far-apart pairs must stay two separate clusters");
    }

    #[test]
    fn drops_clusters_below_min_support() {
        let existing = one_copy("c0", "chr1", 1000, 2000);
        let tied = vec![("read1".to_string(), vec![("chr1".to_string(), 5000, 5100)])];
        let out = cluster_tie_partners(&tied, "FAM0", &existing, 500, 2);
        assert!(out.is_empty(), "a single supporting read must not clear min_support=2");
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
    }
}
