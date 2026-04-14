//! Mixed mode coordinator: handle both long and short reads together.
//!
//! C++ reference mixedMode: process_srfrag redistributes srabund to compatible transfrags

use crate::shortread_pipeline::{ReadClassifier, ReadType};
use crate::types::{BundleRead, RunConfig};

/// Mixed mode read set with separated long and short reads
pub struct MixedReadSet {
    /// Long reads (processed for abundance, longstart/longend)
    pub long_reads: Vec<BundleRead>,
    /// Short reads (processed for srabund, pair-aware)
    pub short_reads: Vec<BundleRead>,
}

impl MixedReadSet {
    /// Separate reads into long and short based on query length threshold
    pub fn classify_reads(reads: Vec<BundleRead>, long_read_min_len: u64) -> Self {
        let classifier = ReadClassifier::new(long_read_min_len);

        let (long_reads, short_reads): (Vec<_>, Vec<_>) = reads
            .into_iter()
            .partition(|r| classifier.classify(r) == ReadType::Long);

        Self {
            long_reads,
            short_reads,
        }
    }

    /// Total read count
    pub fn total_reads(&self) -> usize {
        self.long_reads.len() + self.short_reads.len()
    }

    /// Check if we have both types (true mixed mode)
    pub fn is_truly_mixed(&self) -> bool {
        !self.long_reads.is_empty() && !self.short_reads.is_empty()
    }
}

/// Mixed mode graph coordinator
/// Handles the interaction between long and short read processing
pub struct MixedModeCoordinator {
    config: RunConfig,
}

impl MixedModeCoordinator {
    pub fn new(config: &RunConfig) -> Self {
        Self {
            config: config.clone(),
        }
    }

    /// Process mixed reads
    /// 1. Separate into long and short
    /// 2. Build graph with long-read boundaries (if any)
    /// 3. Apply coverage trim for short-read coverage patterns
    /// 4. Redistribute srabund via process_srfrag
    pub fn process_mixed_reads(
        &self,
        reads: Vec<BundleRead>,
    ) -> anyhow::Result<MixedProcessResult> {
        // Separate reads
        let read_set = MixedReadSet::classify_reads(reads, self.config.long_read_min_len);

        // TODO: Implement mixed mode processing
        // 1. Build base graph from long reads (longtrim)
        // 2. Add short-read coverage and srabund
        // 3. Apply trimnode_all for short-read patterns
        // 4. Redistribute srabund via process_srfrag

        let _ = read_set;
        todo!("Mixed mode processing not yet fully implemented")
    }
}

/// Result of mixed mode processing
pub struct MixedProcessResult {
    /// The unified graph containing both long and short read info
    pub graph: crate::graph::Graph,
    /// Long-read transfrags (for abundance)
    pub long_transfrags: Vec<crate::graph::GraphTransfrag>,
    /// Short-read transfrags (for srabund)
    pub short_transfrags: Vec<crate::graph::GraphTransfrag>,
    /// Statistics for reporting
    pub stats: MixedModeStats,
}

/// Statistics for mixed mode processing
#[derive(Debug, Default)]
pub struct MixedModeStats {
    pub n_long_reads: usize,
    pub n_short_reads: usize,
    pub n_long_transfrags: usize,
    pub n_short_transfrags: usize,
    pub srabund_redistributed: f64,
}

/// Process short-read abundance (srabund) redistribution
/// C++ ref: process_srfrag redistributes srabund to compatible long-read transfrags
pub fn redistribute_srabund(
    short_reads: &[BundleRead],
    transfrags: &mut [crate::graph::GraphTransfrag],
    graph: &crate::graph::Graph,
) -> f64 {
    // TODO: Implement srabund redistribution
    // - For each short read, find compatible long-read transfrags
    // - Distribute srabund among compatible transfrags
    // - Update transfrag.srabund field

    let _ = (short_reads, transfrags, graph);
    0.0
}

/// Mixed mode extraction order
/// Long reads first (by abundance), then short reads (by srabund)
pub fn mixed_extraction_order(transfrags: &[crate::graph::GraphTransfrag]) -> Vec<usize> {
    // TODO: Implement mixed mode ordering
    // - Long transfrags sorted by (abundance, node_count) desc
    // - Short transfrags sorted by (srabund + abundance, node_count) desc
    // - Combine: long first, then short

    let _ = transfrags;
    vec![]
}
