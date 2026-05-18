/// Junction Graph Construction: Exon-endpoint clustering + bridge inference
///
/// Two approaches:
/// 1. Endpoint clustering (RUSTLE_JG_CLUSTER_RADIUS > 0): Normalizes slight position jitter
///    from indel-containing reads by clustering and remapping to canonical positions.
/// 2. Bridge inference (always on when RUSTLE_JUNCTION_GRAPH=1): Finds coverage gaps
///    between exon islands and detects reads bridging the gaps, inferring new junctions.

use crate::types::{BundleRead, Junction, JunctionStat, JunctionStats};
use std::collections::{BTreeMap, HashMap};

/// Maps original exon endpoint positions to canonical (clustered) positions.
pub struct EndpointClusterer {
    /// donor (exon end) position mapping: original_pos -> canonical_pos
    donor_map: BTreeMap<u64, u64>,
    /// acceptor (exon start) position mapping: original_pos -> canonical_pos
    acceptor_map: BTreeMap<u64, u64>,
}

impl EndpointClusterer {
    /// Build an endpoint clusterer from all reads in a bundle.
    ///
    /// Collects all exon endpoints, clusters nearby positions (within cluster_radius),
    /// and computes weighted-mode canonical position per cluster.
    pub fn from_bundle(
        reads: &[BundleRead],
        cluster_radius: u64,
    ) -> Self {
        let mut donor_endpoints: Vec<(u64, f64)> = Vec::new();
        let mut acceptor_endpoints: Vec<(u64, f64)> = Vec::new();

        // Collect all exon endpoints from reads, with weights
        for read in reads {
            for exon in &read.exons {
                let donor = exon.1;
                let acceptor = exon.0;
                donor_endpoints.push((donor, read.weight));
                acceptor_endpoints.push((acceptor, read.weight));
            }
        }

        let donor_map = Self::cluster_and_map(&donor_endpoints, cluster_radius);
        let acceptor_map = Self::cluster_and_map(&acceptor_endpoints, cluster_radius);

        EndpointClusterer {
            donor_map,
            acceptor_map,
        }
    }

    /// Cluster similar positions and return a mapping to canonical positions.
    ///
    /// Groups positions within cluster_radius bp of each other into clusters,
    /// then selects the weighted-mode position (highest total weight) as canonical.
    fn cluster_and_map(
        endpoints: &[(u64, f64)],
        cluster_radius: u64,
    ) -> BTreeMap<u64, u64> {
        if endpoints.is_empty() {
            return BTreeMap::new();
        }

        // Group positions by value
        let mut position_weights: BTreeMap<u64, f64> = BTreeMap::new();
        for &(pos, weight) in endpoints {
            *position_weights.entry(pos).or_insert(0.0) += weight;
        }

        let mut canonical_map = BTreeMap::new();
        let mut sorted_positions: Vec<u64> = position_weights.keys().copied().collect();

        let mut cluster_start = 0;
        while cluster_start < sorted_positions.len() {
            // Find all positions within cluster_radius of the cluster start
            let start_pos = sorted_positions[cluster_start];
            let mut cluster_end = cluster_start;

            while cluster_end < sorted_positions.len()
                && sorted_positions[cluster_end] <= start_pos + cluster_radius
            {
                cluster_end += 1;
            }

            // Find the position with highest weight in this cluster
            let mut max_weight = 0.0;
            let mut canonical = start_pos;

            for i in cluster_start..cluster_end {
                let pos = sorted_positions[i];
                let weight = position_weights[&pos];
                if weight > max_weight {
                    max_weight = weight;
                    canonical = pos;
                }
            }

            // Map all positions in cluster to canonical
            for i in cluster_start..cluster_end {
                let pos = sorted_positions[i];
                canonical_map.insert(pos, canonical);
            }

            cluster_start = cluster_end;
        }

        canonical_map
    }

    /// Remap junctions from a read to use canonical endpoint positions.
    pub fn remap_junctions(&self, junctions: &[Junction]) -> Vec<Junction> {
        junctions
            .iter()
            .map(|j| {
                let canonical_donor =
                    self.donor_map.get(&j.donor).copied().unwrap_or(j.donor);
                let canonical_acceptor =
                    self.acceptor_map.get(&j.acceptor).copied().unwrap_or(j.acceptor);
                Junction::new(canonical_donor, canonical_acceptor)
            })
            .collect()
    }

    /// Remap exon boundaries to canonical endpoint positions.
    pub fn remap_exons(&self, exons: &[(u64, u64)]) -> Vec<(u64, u64)> {
        exons
            .iter()
            .map(|&(start, end)| {
                let canonical_start =
                    self.acceptor_map.get(&start).copied().unwrap_or(start);
                let canonical_end = self.donor_map.get(&end).copied().unwrap_or(end);
                (canonical_start, canonical_end)
            })
            .collect()
    }

    /// Get the mapping of donor positions
    pub fn donor_mapping(&self) -> &BTreeMap<u64, u64> {
        &self.donor_map
    }

    /// Get the mapping of acceptor positions
    pub fn acceptor_mapping(&self) -> &BTreeMap<u64, u64> {
        &self.acceptor_map
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_cluster_endpoints() {
        // Create test endpoints with position jitter
        let endpoints = vec![
            (100u64, 10.0),
            (101u64, 5.0),
            (102u64, 8.0),
            // Cluster 1: 100-102, should canonical to 100 (weight 10.0 is max)
            (200u64, 15.0),
            (202u64, 3.0),
            // Cluster 2: 200-202, should canonical to 200 (weight 15.0 is max)
        ];

        let map = EndpointClusterer::cluster_and_map(&endpoints, 3);

        assert_eq!(map.get(&100), Some(&100)); // canonical
        assert_eq!(map.get(&101), Some(&100)); // mapped
        assert_eq!(map.get(&102), Some(&100)); // mapped
        assert_eq!(map.get(&200), Some(&200)); // canonical
        assert_eq!(map.get(&202), Some(&200)); // mapped
    }

    #[test]
    fn test_remap_junctions() {
        let endpoints = vec![
            (100u64, 10.0),
            (101u64, 5.0),
            (102u64, 8.0),
            (200u64, 15.0),
            (202u64, 3.0),
        ];

        let map = EndpointClusterer::cluster_and_map(&endpoints, 3);
        let donor_map = map;
        let acceptor_map = donor_map.clone(); // In real use, would be separate

        let clusterer = EndpointClusterer {
            donor_map,
            acceptor_map,
        };

        let junction = Junction::new(101, 200); // (donor 101 -> 100, acceptor 200 -> 200)
        let remapped = clusterer.remap_junctions(&[junction]);

        assert_eq!(remapped[0].donor, 100);
        assert_eq!(remapped[0].acceptor, 200);
    }
}

/// Bridge-read junction inference: finds coverage gaps between exon islands
/// and detects reads that bridge the gaps, inferring new junction candidates.
pub struct BridgeJunctionInferer {
    /// Candidate junctions: (donor, acceptor) -> total_bridge_weight
    candidates: BTreeMap<(u64, u64), f64>,
}

impl BridgeJunctionInferer {
    /// Build bridge-junction candidates from all reads in a bundle.
    ///
    /// Algorithm:
    /// 1. Build per-base coverage from reads' exon intervals
    /// 2. Find zero-coverage gaps >= min_intron in length
    /// 3. For each gap, find bridge reads (with exons touching both sides)
    /// 4. Record (left_exon.end, right_exon.start) with read weight
    pub fn from_reads(
        reads: &[BundleRead],
        bundle_start: u64,
        bundle_end: u64,
        min_intron: u64,
        max_intron: u64,
        min_anchor: u64,
    ) -> Self {
        let mut candidates = BTreeMap::new();

        // Build per-base coverage from exon intervals
        let coverage = Self::build_coverage(reads, bundle_start, bundle_end);

        // Find coverage gaps (zero-coverage runs of length >= min_intron)
        let gaps = Self::find_gaps(&coverage, bundle_start, min_intron, max_intron);

        // For each gap, find reads that bridge it and infer junction candidates.
        for (gap_left, gap_right) in gaps {
            for read in reads {
                if read.exons.is_empty() {
                    continue;
                }

                // Find the rightmost exon before the gap
                let left_exon = read.exons.iter()
                    .filter(|exon| exon.1 <= gap_left)
                    .max_by_key(|exon| exon.1);

                // Find the leftmost exon after the gap
                let right_exon = read.exons.iter()
                    .filter(|exon| exon.0 >= gap_right)
                    .min_by_key(|exon| exon.0);

                // If this read bridges the gap, record its junction
                if let (Some(&(left_start, left_end)), Some(&(right_start, right_end))) = (left_exon, right_exon) {
                    let left_len = left_end.saturating_sub(left_start);
                    let right_len = right_end.saturating_sub(right_start);

                    if left_len >= min_anchor && right_len >= min_anchor {
                        let key = (left_end, right_start);
                        *candidates.entry(key).or_insert(0.0) += read.weight;
                    }
                }
            }
        }

        BridgeJunctionInferer { candidates }
    }

    /// Build per-base coverage array from reads' exon intervals.
    /// Returns a Vec where index i represents position (bundle_start + i).
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

    /// Find zero-coverage gaps of length >= min_intron.
    /// Returns list of (gap_left, gap_right) positions.
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

        // Handle gap that extends to end of bundle
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
