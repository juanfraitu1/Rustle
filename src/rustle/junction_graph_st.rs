/// StringTie-style coverage-based junction inference
///
/// Attempts to port StringTie's bundlenode graph algorithm:
/// 1. Identify coverage gaps (potential intron regions)
/// 2. Collect exon endpoints from reads spanning gaps
/// 3. Infer junction candidates from endpoint statistics
///
/// Current implementation uses simple endpoint aggregation.
/// Full StringTie parity would require:
/// - Bundlenode graph construction with graph topology
/// - Multi-pass junction inference from graph structure
/// - EM-based path abundance estimation
/// (See STRINGTIE_PARALLEL_PORT_SUMMARY.md for details)

use crate::types::{BundleRead, Junction, JunctionStat, JunctionStats};
use std::collections::{BTreeMap, HashMap};

/// StringTie-style junction inference from graph topology
pub struct BundlenodeGraphInferer {
    /// Candidate junctions inferred from graph paths
    candidates: BTreeMap<(u64, u64), f64>,
}

impl BundlenodeGraphInferer {
    /// Build junction candidates from bundlenode graph inference.
    ///
    /// Algorithm:
    /// 1. Identify coverage islands (maximal contiguous regions with >0 coverage)
    /// 2. For reads connecting distant islands, infer junctions at island boundaries
    /// 3. Aggregate junction support across all reads
    pub fn from_reads(
        reads: &[BundleRead],
        bundle_start: u64,
        bundle_end: u64,
        min_intron: u64,
        max_intron: u64,
    ) -> Self {
        let mut candidates = BTreeMap::new();

        // NOTE: This is a placeholder implementation matching the coverage-gap approach.
        // True StringTie parity would require full bundlenode graph algorithm (1000+ lines).
        // Current implementation just uses coverage gaps and bridge-read detection,
        // same as the original approach.

        // Build per-base coverage to identify gaps
        let coverage = Self::build_coverage(reads, bundle_start, bundle_end);
        let gaps = Self::find_gaps(&coverage, bundle_start, min_intron, max_intron);

        // For each gap, find bridge reads and aggregate their junction endpoints
        for (gap_left, gap_right) in gaps {
            for read in reads {
                if read.exons.is_empty() {
                    continue;
                }

                // Find rightmost exon before gap and leftmost exon after gap
                let left_exon = read.exons.iter()
                    .filter(|exon| exon.1 <= gap_left)
                    .max_by_key(|exon| exon.1);
                let right_exon = read.exons.iter()
                    .filter(|exon| exon.0 >= gap_right)
                    .min_by_key(|exon| exon.0);

                // If read bridges the gap, record its junction
                if let (Some(&(_, left_end)), Some(&(right_start, _))) = (left_exon, right_exon) {
                    let key = (left_end, right_start);
                    *candidates.entry(key).or_insert(0.0) += read.weight;
                }
            }
        }

        BundlenodeGraphInferer { candidates }
    }

    /// Find zero-coverage gaps
    fn find_gaps(coverage: &[f64], bundle_start: u64, min_intron: u64, max_intron: u64) -> Vec<(u64, u64)> {
        let mut gaps = Vec::new();
        let mut gap_start: Option<usize> = None;

        for (i, &cov) in coverage.iter().enumerate() {
            if cov == 0.0 {
                if gap_start.is_none() {
                    gap_start = Some(i);
                }
            } else if let Some(start) = gap_start {
                let gap_len = i - start;
                if gap_len >= min_intron as usize {
                    let gap_left = bundle_start + start as u64;
                    let gap_right = bundle_start + i as u64;
                    if gap_right - gap_left <= max_intron {
                        gaps.push((gap_left, gap_right));
                    }
                }
                gap_start = None;
            }
        }

        // Handle gap that extends to end
        if let Some(start) = gap_start {
            let gap_len = coverage.len() - start;
            if gap_len >= min_intron as usize {
                let gap_left = bundle_start + start as u64;
                let gap_right = bundle_start + coverage.len() as u64;
                if gap_right - gap_left <= max_intron {
                    gaps.push((gap_left, gap_right));
                }
            }
        }

        gaps
    }

    /// Build per-base coverage from exon intervals
    fn build_coverage(reads: &[BundleRead], bundle_start: u64, bundle_end: u64) -> Vec<f64> {
        let len = (bundle_end - bundle_start) as usize;
        let mut coverage = vec![0.0; len];

        for read in reads {
            for &(exon_start, exon_end) in &read.exons {
                let start_idx = (exon_start.saturating_sub(bundle_start)) as usize;
                let end_idx = ((exon_end.saturating_sub(bundle_start)).min(bundle_end - bundle_start))
                    as usize;

                if start_idx < len {
                    for i in start_idx..end_idx.min(len) {
                        coverage[i] += read.weight;
                    }
                }
            }
        }

        coverage
    }

    /// Identify coverage islands: maximal contiguous regions with >0 coverage
    fn identify_coverage_islands(coverage: &[f64], bundle_start: u64) -> Vec<(u64, u64)> {
        let mut islands = Vec::new();
        let mut island_start: Option<usize> = None;

        for (i, &cov) in coverage.iter().enumerate() {
            if cov > 0.0 {
                if island_start.is_none() {
                    island_start = Some(i);
                }
            } else if let Some(start) = island_start {
                islands.push((
                    bundle_start + start as u64,
                    bundle_start + i as u64,
                ));
                island_start = None;
            }
        }

        // Handle island that extends to end
        if let Some(start) = island_start {
            islands.push((
                bundle_start + start as u64,
                bundle_start + coverage.len() as u64,
            ));
        }

        islands
    }

    /// Inject inferred junctions into existing JunctionStats.
    /// Only adds junctions NOT already present.
    pub fn inject_into(&self, stats: &mut JunctionStats, min_bridge_reads: f64) {
        for ((donor, acceptor), weight) in &self.candidates {
            if *weight >= min_bridge_reads {
                let junc = Junction::new(*donor, *acceptor);
                // Only add if not already present
                if !stats.contains_key(&junc) {
                    let mut junc_stat = JunctionStat::default();
                    junc_stat.mrcount = *weight;
                    junc_stat.nreads_good = *weight;
                    junc_stat.leftsupport = *weight;
                    junc_stat.rightsupport = *weight;
                    junc_stat.mm = *weight;
                    junc_stat.strand = None;
                    junc_stat.nm = 0.0;
                    junc_stat.consleft = -1;
                    junc_stat.consright = -1;
                    stats.insert(junc, junc_stat);
                }
            }
        }
    }
}
