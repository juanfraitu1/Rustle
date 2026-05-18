//! Graph comparison: diagnose missing junctions by comparing against reference GTF.
//!
//! This module enables systematic diagnosis of missing junctions by comparing
//! de novo assembly output against a reference GTF. For each missing junction,
//! we trace through the pipeline to identify where it gets lost.

use std::collections::{HashMap, HashSet};
use std::path::Path;
use anyhow::{anyhow, Result};
use crate::types::{Junction, JunctionStats, DetHashMap};

/// Reference graph structure extracted from GTF.
#[derive(Debug, Clone)]
pub struct ReferenceGraphStructure {
    /// All junctions found in reference GTF: (donor, acceptor) -> transcript IDs using it
    pub junctions: HashMap<(u64, u64), Vec<String>>,
    /// Metadata about each junction
    pub junction_info: HashMap<(u64, u64), JunctionInfo>,
}

#[derive(Debug, Clone)]
pub struct JunctionInfo {
    pub donor: u64,
    pub acceptor: u64,
    pub transcript_ids: Vec<String>,
    pub gene_ids: Vec<String>,
    pub count: usize,
}

/// Comparison results between reference and discovered junctions.
#[derive(Debug)]
pub struct ComparisonResult {
    pub found: HashSet<(u64, u64)>,
    pub missing: HashSet<(u64, u64)>,
    pub extra: HashSet<(u64, u64)>,
    pub metrics: ComparisonMetrics,
    pub missing_details: Vec<MissingJunctionDetail>,
}

#[derive(Debug, Clone)]
pub struct ComparisonMetrics {
    pub total_reference: usize,
    pub found_count: usize,
    pub missing_count: usize,
    pub extra_count: usize,
    pub sensitivity: f64,  // found / (found + missing)
    pub specificity: f64,  // found / (found + extra)
}

#[derive(Debug, Clone)]
pub struct MissingJunctionDetail {
    pub junction: (u64, u64),
    pub donor: u64,
    pub acceptor: u64,
    pub reference_transcripts: Vec<String>,
    pub reference_genes: Vec<String>,
    pub count_in_reference: usize,
}

impl ReferenceGraphStructure {
    /// Load reference junctions from GTF file.
    pub fn from_gtf(gtf_path: &Path) -> Result<Self> {
        eprintln!("[Graph-Comparison] Loading reference from {:?}...", gtf_path);

        let content = std::fs::read_to_string(gtf_path)
            .map_err(|e| anyhow!("Failed to read GTF: {}", e))?;

        let mut by_transcript: HashMap<String, Vec<(u64, u64)>> = HashMap::new();
        let mut transcript_genes: HashMap<String, String> = HashMap::new();

        for line in content.lines() {
            if line.starts_with('#') || line.trim().is_empty() {
                continue;
            }

            let cols: Vec<&str> = line.split('\t').collect();
            if cols.len() < 9 {
                continue;
            }

            let feature = cols[2];
            if feature != "exon" {
                continue;
            }

            let start_1b: u64 = cols[3].parse().ok().unwrap_or(0);
            let end_1b: u64 = cols[4].parse().ok().unwrap_or(0);

            if start_1b == 0 || end_1b == 0 || start_1b >= end_1b {
                continue;
            }

            // Convert 1-based to 0-based
            let start_0 = start_1b.saturating_sub(1);
            let end_0_excl = end_1b;

            // Extract transcript_id and gene_id
            let attrs = cols[8];
            let transcript_id = extract_gtf_attr(attrs, "transcript_id")
                .unwrap_or_else(|| format!("tx_{}", line.len()));
            let gene_id = extract_gtf_attr(attrs, "gene_id")
                .unwrap_or_else(|| format!("gene_{}", line.len()));

            transcript_genes.insert(transcript_id.clone(), gene_id);
            by_transcript
                .entry(transcript_id)
                .or_insert_with(Vec::new)
                .push((start_0, end_0_excl));
        }

        // Extract junctions from exons
        let mut junctions: HashMap<(u64, u64), Vec<String>> = HashMap::new();
        let mut junction_info: HashMap<(u64, u64), JunctionInfo> = HashMap::new();

        for (transcript_id, mut exons) in by_transcript {
            exons.sort_by_key(|e| e.0);
            exons.dedup();

            // Create junctions between consecutive exons
            for i in 0..exons.len() - 1 {
                let donor = exons[i].1;     // exon end
                let acceptor = exons[i + 1].0;  // next exon start

                let junction = (donor, acceptor);
                junctions
                    .entry(junction)
                    .or_insert_with(Vec::new)
                    .push(transcript_id.clone());

                let gene_id = transcript_genes
                    .get(&transcript_id)
                    .cloned()
                    .unwrap_or_else(|| "unknown".to_string());

                junction_info
                    .entry(junction)
                    .and_modify(|info| {
                        info.transcript_ids.push(transcript_id.clone());
                        if !info.gene_ids.contains(&gene_id) {
                            info.gene_ids.push(gene_id.clone());
                        }
                        info.count += 1;
                    })
                    .or_insert_with(|| JunctionInfo {
                        donor,
                        acceptor,
                        transcript_ids: vec![transcript_id.clone()],
                        gene_ids: vec![gene_id],
                        count: 1,
                    });
            }
        }

        eprintln!("[Graph-Comparison] Found {} junctions in reference", junctions.len());

        Ok(ReferenceGraphStructure {
            junctions,
            junction_info,
        })
    }

    /// Compare reference junctions against discovered junctions.
    pub fn compare_junctions(
        &self,
        discovered: &DetHashMap<Junction, crate::types::JunctionStat>,
    ) -> ComparisonResult {
        let mut found = HashSet::new();
        let mut missing = HashSet::new();
        let mut missing_details = Vec::new();

        // Check each reference junction
        for (&(donor, acceptor), txs) in &self.junctions {
            let junction = Junction::new(donor, acceptor);

            if discovered.contains_key(&junction) {
                found.insert((donor, acceptor));
            } else {
                missing.insert((donor, acceptor));

                let info = self.junction_info.get(&(donor, acceptor)).cloned();
                if let Some(info) = info {
                    missing_details.push(MissingJunctionDetail {
                        junction: (donor, acceptor),
                        donor,
                        acceptor,
                        reference_transcripts: info.transcript_ids,
                        reference_genes: info.gene_ids,
                        count_in_reference: info.count,
                    });
                }
            }
        }

        // Check for extra junctions (in discovered but not in reference)
        let mut extra = HashSet::new();
        for junction in discovered.keys() {
            let key = (junction.donor, junction.acceptor);
            if !self.junctions.contains_key(&key) {
                extra.insert(key);
            }
        }

        let found_count = found.len();
        let missing_count = missing.len();
        let extra_count = extra.len();
        let total_reference = self.junctions.len();

        let sensitivity = if found_count + missing_count > 0 {
            found_count as f64 / (found_count + missing_count) as f64
        } else {
            0.0
        };

        let specificity = if found_count + extra_count > 0 {
            found_count as f64 / (found_count + extra_count) as f64
        } else {
            0.0
        };

        // Sort missing details by count (most important first)
        missing_details.sort_by(|a, b| b.count_in_reference.cmp(&a.count_in_reference));

        ComparisonResult {
            found,
            missing,
            extra,
            metrics: ComparisonMetrics {
                total_reference,
                found_count,
                missing_count,
                extra_count,
                sensitivity,
                specificity,
            },
            missing_details,
        }
    }

    /// Print detailed comparison report.
    pub fn print_report(&self, comparison: &ComparisonResult) {
        eprintln!("\n[Graph-Comparison] ════════════════════════════════════════");
        eprintln!("[Graph-Comparison] COMPARISON RESULTS");
        eprintln!("[Graph-Comparison] ════════════════════════════════════════");
        eprintln!(
            "[Graph-Comparison] Total reference junctions: {}",
            comparison.metrics.total_reference
        );
        eprintln!(
            "[Graph-Comparison] Found in discovered: {} ({:.1}%)",
            comparison.metrics.found_count,
            comparison.metrics.sensitivity * 100.0
        );
        eprintln!(
            "[Graph-Comparison] Missing from discovered: {}",
            comparison.metrics.missing_count
        );
        eprintln!(
            "[Graph-Comparison] Extra in discovered: {} ({:.1}%)",
            comparison.metrics.extra_count,
            (comparison.metrics.extra_count as f64 / comparison.metrics.found_count as f64) * 100.0
        );

        eprintln!(
            "[Graph-Comparison] Sensitivity (Sn): {:.1}%",
            comparison.metrics.sensitivity * 100.0
        );
        eprintln!(
            "[Graph-Comparison] Specificity (Pr): {:.1}%",
            comparison.metrics.specificity * 100.0
        );

        // Show first 10 missing junctions
        eprintln!("\n[Graph-Comparison] Top {} missing junctions:", comparison.missing_details.len().min(15));
        for (i, detail) in comparison.missing_details.iter().take(15).enumerate() {
            eprintln!(
                "  {}. ({}, {}) — used by {} transcripts: {}",
                i + 1,
                detail.donor,
                detail.acceptor,
                detail.reference_transcripts.len(),
                detail.reference_transcripts.join(", ")
            );
            if !detail.reference_genes.is_empty() {
                eprintln!(
                    "     Genes: {}",
                    detail.reference_genes.join(", ")
                );
            }
        }

        eprintln!("\n[Graph-Comparison] ════════════════════════════════════════\n");
    }
}

/// Extract attribute value from GTF attributes string.
fn extract_gtf_attr(attrs: &str, key: &str) -> Option<String> {
    for part in attrs.split(';') {
        let part = part.trim();
        if part.starts_with(key) {
            let rest = part.strip_prefix(key)?;
            let rest = rest.trim_start();
            let rest = rest.strip_prefix(' ')?;
            let rest = rest.strip_prefix('"')?;
            let end = rest.find('"')?;
            return Some(rest[..end].to_string());
        }
    }
    None
}
