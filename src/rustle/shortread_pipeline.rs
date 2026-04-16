//! Short-read pipeline: dedicated processing for short-read RNA-seq assembly.
//!
//! differences for short-read mode:
//! - trimnode_all: coverage-based node splitting (always on for short-read)
//! - srabund: short-read abundance tracking
//! - pair-aware: paired-end fragment bridging
//! - localdist = bundledist + longintronanchor = 50 + 25 = 75

use crate::types::{BundleRead, RunConfig};

/// Short-read specific constants (header / original algorithm.cpp)
pub mod constants {
    /// Default bundle merge distance for short-read (original algorithm bundledist=50)
    pub const BUNDLE_MERGE_DIST_SHORTREAD: u64 = 50;

    /// localdist for short-read: bundledist + longintronanchor = 50 + 25 = 75
    pub const LOCALDIST_SHORTREAD: u64 = 75;

    /// Coverage drop threshold for trimnode_all (DROP = 0.5)
    pub const COVERAGE_DROP_THRESHOLD: f64 = 0.5;

    /// Minimum coverage for node to be kept in trimnode_all
    pub const MIN_NODE_COVERAGE: f64 = 1.0;
}

/// Classify reads for mixed mode (long vs short)
pub struct ReadClassifier {
    /// Minimum query length to be considered long-read
    long_read_min_len: u64,
}

impl ReadClassifier {
    pub fn new(long_read_min_len: u64) -> Self {
        Self { long_read_min_len }
    }

    /// Classify a read as long or short based on query length
    /// mixedMode classifies reads with query_length >= long_read_min_len as long
    pub fn classify(&self, read: &BundleRead) -> ReadType {
        match read.query_length {
            Some(len) if len >= self.long_read_min_len => ReadType::Long,
            _ => ReadType::Short,
        }
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ReadType {
    Long,
    Short,
}

/// Short-read graph builder with coverage-based node splitting
pub struct ShortReadGraphBuilder {
    _config: RunConfig,
}

impl ShortReadGraphBuilder {
    pub fn new(config: &RunConfig) -> Self {
        Self {
            _config: config.clone(),
        }
    }

    /// Build graph for short-read mode
    /// create_graph with trimnode_all, no longtrim
    pub fn build_graph(&self, reads: &[BundleRead]) -> anyhow::Result<crate::graph::Graph> {
        // TODO: Implement short-read specific graph building
        // - Use bundledist = 50 for grouping
        // - Apply trimnode_all for coverage-based splitting
        // - No longtrim (only for long-read mode)
        // - Pair-aware localdist = 75

        // Placeholder: would call existing graph building with short-read params
        let _ = reads;
        todo!("Short-read graph building not yet implemented")
    }
}

/// Short-read coverage trimmer (trimnode_all)
pub struct CoverageTrimmer {
    _drop_threshold: f64,
}

impl CoverageTrimmer {
    pub fn new(drop_threshold: f64) -> Self {
        Self {
            _drop_threshold: drop_threshold,
        }
    }

    /// Trim nodes based on coverage drops
    /// ref: trimnode_all - split nodes at significant coverage drops
    pub fn trim_nodes(
        &self,
        graph: &mut crate::graph::Graph,
        bpcov: &crate::bpcov::Bpcov,
    ) -> Vec<usize> {
        // TODO: Implement coverage-based node trimming
        // - Find coverage drops > threshold
        // - Split nodes at drop points
        // - Return list of new node IDs

        let _ = (graph, bpcov);
        vec![]
    }
}

/// Short-read abundance tracker (srabund)
pub struct ShortReadAbundance {
    /// Total short-read abundance per node
    node_srabund: Vec<f64>,
}

impl ShortReadAbundance {
    pub fn new(n_nodes: usize) -> Self {
        Self {
            node_srabund: vec![0.0; n_nodes],
        }
    }

    /// Add short-read abundance to a node
    pub fn add_abundance(&mut self, node_id: usize, abundance: f64) {
        if let Some(entry) = self.node_srabund.get_mut(node_id) {
            *entry += abundance;
        }
    }

    /// Get srabund for a node
    pub fn get(&self, node_id: usize) -> f64 {
        self.node_srabund.get(node_id).copied().unwrap_or(0.0)
    }

    /// Compute srabund for all nodes from reads
    pub fn compute_from_reads(&mut self, reads: &[BundleRead], graph: &crate::graph::Graph) {
        // TODO: Accumulate srabund from short reads
        // - For each short read, add its weight to covered nodes
        // - Handle paired-end reads (split weight appropriately)

        let _ = (reads, graph);
    }
}

/// Pre-process reads for short-read mode
/// - Trim polyA/T
/// - Filter by quality
/// - Sort by position
pub fn preprocess_short_reads(mut reads: Vec<BundleRead>) -> Vec<BundleRead> {
    // Sort by start position (required for bundle building)
    reads.sort_by_key(|r| (r.ref_start, r.ref_end));

    // TODO: Additional preprocessing
    // - Trim low-quality ends
    // - Filter by NM/MD tags
    // - Mark duplicates

    reads
}

/// Build bundles from short reads
/// bundledist=50, localdist=75 for short-read
pub fn build_shortread_bundles(reads: &[BundleRead], bundledist: u64) -> Vec<crate::types::Bundle> {
    // TODO: Implement short-read bundle building
    // - Group reads within bundledist
    // - Use pair info for bundle boundaries
    // - Return bundles ready for graph building

    let _ = (reads, bundledist);
    vec![]
}
