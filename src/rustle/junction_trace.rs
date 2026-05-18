//! Per-junction tracing: follow a junction through each pipeline stage.
//!
//! This module helps diagnose where missing junctions get lost by tracking
//! them from raw reads all the way through graph construction and filtering.

use std::collections::HashSet;
use crate::types::{Junction, BundleRead};

/// Trace a single junction through the pipeline stages.
#[derive(Debug, Clone)]
pub struct JunctionTrace {
    pub donor: u64,
    pub acceptor: u64,

    // Stage 1: In reads
    pub in_read_junctions: bool,
    pub spanning_reads: usize,
    pub read_coverage: f64,
    pub read_details: Vec<String>,  // List of read names/indices that span this junction

    // Stage 2: In bundle
    pub in_bundle: bool,
    pub bundle_start: Option<u64>,
    pub bundle_end: Option<u64>,

    // Stage 3: In computed junction_stats
    pub in_junction_stats: bool,
    pub stats_weight: Option<f64>,
    pub stats_mrcount: Option<f64>,
    pub stats_nreads: Option<f64>,

    // Stage 4: Exon nodes exist
    pub left_exon_found: bool,
    pub right_exon_found: bool,

    // Stage 5: Edge in graph
    pub edge_exists: bool,

    // Stage 6: Filtering
    pub passed_all_filters: bool,
    pub filters_applied: Vec<(String, bool)>,  // (filter_name, passed)

    // Summary
    pub bottleneck_stage: Option<String>,
}

impl JunctionTrace {
    /// Create a new trace for a junction.
    pub fn new(donor: u64, acceptor: u64) -> Self {
        JunctionTrace {
            donor,
            acceptor,

            in_read_junctions: false,
            spanning_reads: 0,
            read_coverage: 0.0,
            read_details: Vec::new(),

            in_bundle: false,
            bundle_start: None,
            bundle_end: None,

            in_junction_stats: false,
            stats_weight: None,
            stats_mrcount: None,
            stats_nreads: None,

            left_exon_found: false,
            right_exon_found: false,

            edge_exists: false,

            passed_all_filters: true,
            filters_applied: Vec::new(),

            bottleneck_stage: None,
        }
    }

    /// Check if this junction is in the read's junctions.
    pub fn check_read_support(&mut self, reads: &[BundleRead]) {
        let junction = Junction::new(self.donor, self.acceptor);

        for (read_idx, read) in reads.iter().enumerate() {
            if read.junctions.contains(&junction) {
                self.spanning_reads += 1;
                self.read_coverage += read.weight;
                self.read_details.push(format!("read[{}] weight={:.2}", read_idx, read.weight));
            }
        }

        self.in_read_junctions = self.spanning_reads > 0;

        // Record bottleneck if no read support
        if !self.in_read_junctions && self.bottleneck_stage.is_none() {
            self.bottleneck_stage = Some("NO_READ_SUPPORT".to_string());
        }
    }

    /// Record that this junction exists in junction_stats.
    pub fn mark_in_stats(&mut self, weight: f64, mrcount: f64, nreads: f64) {
        self.in_junction_stats = true;
        self.stats_weight = Some(weight);
        self.stats_mrcount = Some(mrcount);
        self.stats_nreads = Some(nreads);
    }

    /// Record that this junction was NOT in junction_stats.
    pub fn mark_not_in_stats(&mut self) {
        self.in_junction_stats = false;
        if self.bottleneck_stage.is_none() && self.in_read_junctions {
            self.bottleneck_stage = Some("NOT_IN_JUNCTION_STATS".to_string());
        }
    }

    /// Record node existence.
    pub fn mark_exons_found(&mut self, left: bool, right: bool) {
        self.left_exon_found = left;
        self.right_exon_found = right;

        if !left || !right {
            if self.bottleneck_stage.is_none() {
                self.bottleneck_stage = Some("EXON_NOT_FOUND".to_string());
            }
        }
    }

    /// Record edge existence.
    pub fn mark_edge(&mut self, exists: bool) {
        self.edge_exists = exists;

        if !exists && self.bottleneck_stage.is_none() && self.in_junction_stats {
            self.bottleneck_stage = Some("EDGE_NOT_CREATED".to_string());
        }
    }

    /// Record filter result.
    pub fn record_filter(&mut self, filter_name: &str, passed: bool) {
        self.filters_applied.push((filter_name.to_string(), passed));

        if !passed {
            self.passed_all_filters = false;
            if self.bottleneck_stage.is_none() {
                self.bottleneck_stage = Some(format!("FILTERED_BY_{}", filter_name));
            }
        }
    }

    /// Generate human-readable report.
    pub fn print_report(&self) {
        eprintln!("\n[Trace] Junction ({}, {}):", self.donor, self.acceptor);

        // Stage 1: Read support
        eprintln!("  Stage 1 - Read Support:");
        if self.in_read_junctions {
            eprintln!("    ✓ Found in {} reads, coverage={:.2}", self.spanning_reads, self.read_coverage);
            if !self.read_details.is_empty() && self.read_details.len() <= 3 {
                for detail in &self.read_details {
                    eprintln!("      {}", detail);
                }
            }
        } else {
            eprintln!("    ✗ NOT found in any read ← BOTTLENECK");
        }

        // Stage 2: Bundle
        eprintln!("  Stage 2 - Bundle Coverage:");
        if self.in_bundle {
            eprintln!("    ✓ Within bundle {:?}-{:?}", self.bundle_start, self.bundle_end);
        } else {
            eprintln!("    ? Bundle info not recorded");
        }

        // Stage 3: Junction stats
        eprintln!("  Stage 3 - Junction Stats:");
        if self.in_junction_stats {
            eprintln!("    ✓ In stats: weight={:.2}, mrcount={:.2}, nreads={:.2}",
                     self.stats_weight.unwrap_or(0.0),
                     self.stats_mrcount.unwrap_or(0.0),
                     self.stats_nreads.unwrap_or(0.0));
        } else if self.in_read_junctions {
            eprintln!("    ✗ NOT in stats despite read support ← BOTTLENECK");
        }

        // Stage 4: Exons/nodes
        eprintln!("  Stage 4 - Exon Nodes:");
        eprintln!("    Left ({}):  {}", self.donor, if self.left_exon_found { "✓" } else { "✗" });
        eprintln!("    Right ({}):{}", self.acceptor, if self.right_exon_found { "✓" } else { "✗" });

        // Stage 5: Edge
        eprintln!("  Stage 5 - Graph Edge:");
        if self.edge_exists {
            eprintln!("    ✓ Edge created");
        } else if self.in_junction_stats {
            eprintln!("    ✗ Edge NOT created ← BOTTLENECK");
        }

        // Stage 6: Filters
        eprintln!("  Stage 6 - Filtering:");
        if self.filters_applied.is_empty() {
            eprintln!("    ? No filters recorded");
        } else {
            let mut any_failed = false;
            for (filter, passed) in &self.filters_applied {
                let status = if *passed { "✓" } else { "✗" };
                eprintln!("    {} {}", status, filter);
                if !*passed {
                    any_failed = true;
                }
            }
            if any_failed {
                eprintln!("    ← BOTTLENECK");
            }
        }

        // Summary
        eprintln!("\n  ROOT CAUSE: {}",
                 self.bottleneck_stage.as_ref()
                     .unwrap_or(&"UNKNOWN".to_string()));
    }
}

/// Collect traces for a set of target junctions.
pub fn parse_target_junctions(env_var: &str) -> Vec<(u64, u64)> {
    if let Ok(targets_str) = std::env::var(env_var) {
        let mut targets = Vec::new();
        // Format: "(donor1,acceptor1);(donor2,acceptor2);..."
        for pair in targets_str.split(';') {
            let pair = pair.trim();
            if pair.starts_with('(') && pair.ends_with(')') {
                let inner = &pair[1..pair.len()-1];
                let parts: Vec<&str> = inner.split(',').collect();
                if parts.len() == 2 {
                    if let (Ok(d), Ok(a)) = (parts[0].parse::<u64>(), parts[1].parse::<u64>()) {
                        targets.push((d, a));
                    }
                }
            }
        }
        targets
    } else {
        Vec::new()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_parse_targets() {
        std::env::set_var("TEST_TARGETS", "(100,200);(300,400)");
        let targets = parse_target_junctions("TEST_TARGETS");
        assert_eq!(targets.len(), 2);
        assert_eq!(targets[0], (100, 200));
        assert_eq!(targets[1], (300, 400));
    }
}
