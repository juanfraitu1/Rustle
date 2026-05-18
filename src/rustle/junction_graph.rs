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

/// Graph-based junction inference with fuzzy boundary matching
///
/// Algorithm:
/// 1. Collect all unique exon endpoints from reads
/// 2. Build a map of exon boundaries and their coverage
/// 3. For junctions where reads have exons near both sides:
///    - If junction has no exact read span, add as inferred candidate
/// 4. Uses fuzzy matching to handle coordinate jitter from indels
pub struct ExonGraphInferer {
    /// All unique exon intervals (start, end) with coverage
    exons: Vec<((u64, u64), f64)>,
    /// All exon start positions with coverage
    exon_starts: BTreeMap<u64, f64>,
    /// All exon end positions with coverage
    exon_ends: BTreeMap<u64, f64>,
    /// Inferred junctions: (donor, acceptor) -> weight
    inferred_junctions: BTreeMap<(u64, u64), f64>,
}

impl ExonGraphInferer {
    /// Build exon graph and identify inferred junctions
    pub fn from_reads(
        reads: &[BundleRead],
        min_intron: u64,
        max_intron: u64,
    ) -> Self {
        // Step 1: Collect all unique exons with their coverage
        let mut exon_coverage: BTreeMap<(u64, u64), f64> = BTreeMap::new();
        let mut exon_starts: BTreeMap<u64, f64> = BTreeMap::new();
        let mut exon_ends: BTreeMap<u64, f64> = BTreeMap::new();

        for read in reads {
            for &exon in &read.exons {
                let (start, end) = exon;
                *exon_coverage.entry(exon).or_insert(0.0) += read.weight;
                *exon_starts.entry(start).or_insert(0.0) += read.weight;
                *exon_ends.entry(end).or_insert(0.0) += read.weight;
            }
        }

        let exons: Vec<((u64, u64), f64)> = exon_coverage
            .iter()
            .map(|(e, cov)| (*e, *cov))
            .collect();

        // Step 2: For each read, find all exon pairs it bridges
        let mut inferred_junctions: BTreeMap<(u64, u64), f64> = BTreeMap::new();

        for read in reads {
            if read.exons.len() < 2 {
                continue;
            }

            // For this read, find all multi-hop exon pairs (not just consecutive)
            for i in 0..read.exons.len() {
                for j in (i + 1)..read.exons.len() {
                    let exon_i = read.exons[i];
                    let exon_j = read.exons[j];

                    let donor = exon_i.1;      // end of exon i
                    let acceptor = exon_j.0;   // start of exon j

                    let gap = if acceptor > donor {
                        acceptor - donor
                    } else {
                        continue; // Invalid junction
                    };

                    // Check intron length is valid
                    if gap >= min_intron && gap <= max_intron {
                        // Use minimum coverage of the two exons as the weight
                        let weight = exon_coverage[&exon_i].min(exon_coverage[&exon_j]);
                        *inferred_junctions.entry((donor, acceptor)).or_insert(0.0) += weight;
                    }
                }
            }
        }

        ExonGraphInferer {
            exons,
            exon_starts,
            exon_ends,
            inferred_junctions,
        }
    }

    /// Find exon boundaries near a target position (fuzzy matching)
    /// Returns coverage-weighted boundary positions within tolerance
    fn find_nearby_boundary(
        boundary_map: &BTreeMap<u64, f64>,
        target: u64,
        tolerance: u64,
    ) -> Option<(u64, f64)> {
        let mut best = None;
        let mut best_weight = 0.0;

        for (&pos, &weight) in boundary_map.range((target.saturating_sub(tolerance))..=(target.saturating_add(tolerance))) {
            if weight > best_weight {
                best = Some(pos);
                best_weight = weight;
            }
        }

        best.map(|pos| (pos, best_weight))
    }

    /// Try to find topologically valid boundary for a coordinate
    /// Returns highest-weight boundary within tolerance
    pub fn find_best_boundary(
        &self,
        is_donor: bool,  // true = look for exon ends, false = look for exon starts
        target: u64,
        tolerance: u64,
    ) -> Option<(u64, f64)> {
        let boundary_map = if is_donor {
            &self.exon_ends
        } else {
            &self.exon_starts
        };

        Self::find_nearby_boundary(boundary_map, target, tolerance)
    }

    /// Check if a junction is topologically valid (has exons on both sides)
    pub fn is_topologically_valid(
        &self,
        donor: u64,
        acceptor: u64,
        tolerance: u64,
        min_intron: u64,
        max_intron: u64,
    ) -> Option<f64> {
        // Find exon ending near donor
        let (donor_actual, donor_weight) = Self::find_nearby_boundary(&self.exon_ends, donor, tolerance)?;

        // Find exon starting near acceptor
        let (acceptor_actual, acceptor_weight) = Self::find_nearby_boundary(&self.exon_starts, acceptor, tolerance)?;

        // Check intron length
        if acceptor_actual <= donor_actual {
            return None;
        }

        let gap = acceptor_actual - donor_actual;
        if gap < min_intron || gap > max_intron {
            return None;
        }

        // Return minimum weight (conservative estimate)
        Some(donor_weight.min(acceptor_weight))
    }

    /// Inject inferred junctions into existing JunctionStats
    /// Only adds junctions NOT already present
    pub fn inject_into(&self, stats: &mut JunctionStats, min_weight: f64) {
        for ((donor, acceptor), weight) in &self.inferred_junctions {
            if *weight >= min_weight {
                let junc = Junction::new(*donor, *acceptor);

                // Only add if not already present (preserve original read-discovered junctions)
                if !stats.contains_key(&junc) {
                    let mut junc_stat = JunctionStat::default();
                    // Use the inferred weight (may be from multi-exon read coverage)
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

                    // Trace: log inferred junction for debugging
                    if let Ok(_) = std::env::var("RUSTLE_TRACE_INFERRED_JUNCTIONS") {
                        eprintln!("[INFERRED-JUNCTION] ({}, {}) weight={:.2}", donor, acceptor, weight);
                    }
                }
            }
        }
    }

    /// Inject a specific junction if it's topologically valid
    /// Used for recovering reference junctions with fuzzy boundary matching
    pub fn inject_if_valid(
        &self,
        stats: &mut JunctionStats,
        donor: u64,
        acceptor: u64,
        tolerance: u64,
        min_intron: u64,
        max_intron: u64,
        boost_factor: f64,
    ) {
        if stats.contains_key(&Junction::new(donor, acceptor)) {
            return; // Already present
        }

        if let Some(weight) = self.is_topologically_valid(donor, acceptor, tolerance, min_intron, max_intron) {
            if weight > 0.0 {
                let junc = Junction::new(donor, acceptor);
                let boosted_weight = weight * boost_factor;
                let mut junc_stat = JunctionStat::default();
                junc_stat.mrcount = boosted_weight;
                junc_stat.nreads_good = boosted_weight;
                junc_stat.leftsupport = boosted_weight;
                junc_stat.rightsupport = boosted_weight;
                junc_stat.mm = boosted_weight;
                junc_stat.strand = None;
                junc_stat.nm = 0.0;
                junc_stat.consleft = -1;
                junc_stat.consright = -1;

                stats.insert(junc, junc_stat);

                if let Ok(_) = std::env::var("RUSTLE_TRACE_INFERRED_JUNCTIONS") {
                    eprintln!("[FUZZY-INJECT] ({}, {}) weight={:.2} boosted to {:.2} (topologically valid)", donor, acceptor, weight, boosted_weight);
                }
            }
        }
    }
}
