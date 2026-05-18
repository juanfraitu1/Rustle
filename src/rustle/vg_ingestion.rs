//! GTF/GFF ingestion mode: Parse transcripts and create template families for VG analysis.
//!
//! This mode allows Rustle to ingest a reference GTF/GFF file (e.g., StringTie output)
//! and use those transcripts as templates for discovering new isoforms or family members
//! via multi-mapped read reassignment with HMM-EM.
//!
//! Flow:
//! 1. Parse GTF → extract transcripts with exons
//! 2. Group transcripts into families (by gene, overlap, or k-mer similarity)
//! 3. Find bundles overlapping each family
//! 4. Link bundles into FamilyGroups
//! 5. Feed into standard HMM-EM + assembly pipeline

use crate::types::Bundle;
use crate::vg::FamilyGroup;
use crate::genome::GenomeIndex;
use anyhow::{anyhow, Result};
use std::collections::HashMap;
use std::path::Path;

// ── Data Structures ────────────────────────────────────────────────────────

/// A transcript template parsed from GTF.
#[derive(Debug, Clone)]
pub struct TemplateTranscript {
    pub transcript_id: String,
    pub gene_id: Option<String>,
    pub chrom: String,
    pub strand: char,
    /// Exons (0-based, end-exclusive) sorted by position
    pub exons: Vec<(u64, u64)>,
}

/// A family of template transcripts (grouped by gene or similarity).
#[derive(Debug, Clone)]
pub struct TemplateFamily {
    pub family_id: usize,
    pub gene_id: Option<String>,
    pub chrom: String,
    pub strand: char,
    pub transcripts: Vec<TemplateTranscript>,
    pub start: u64,
    pub end: u64,
}

/// Strategy for grouping transcripts into families.
#[derive(Debug, Clone, Copy)]
pub enum FamilyGroupingStrategy {
    /// Group by gene_id attribute (simplest, most conservative)
    ByGene,
    /// Group overlapping transcripts (by genomic range)
    ByOverlap,
}

// ── GTF Parsing ────────────────────────────────────────────────────────────

/// Parse a GTF/GFF file and extract transcripts.
///
/// Returns a list of transcripts with:
/// - transcript_id and gene_id (if available)
/// - exons sorted by start position
/// - 0-based coordinates (start inclusive, end exclusive)
pub fn parse_gtf_transcripts<P: AsRef<Path>>(path: P) -> Result<Vec<TemplateTranscript>> {
    let content = std::fs::read_to_string(path.as_ref())
        .map_err(|e| anyhow!("Failed to read GTF: {}", e))?;

    let mut by_transcript: HashMap<String, TemplateTranscript> = HashMap::new();

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

        let chrom = cols[0].to_string();
        let start_1b: u64 = cols[3].parse().ok().unwrap_or(0);
        let end_1b: u64 = cols[4].parse().ok().unwrap_or(0);
        let strand = cols[6].chars().next().unwrap_or('.');

        if start_1b == 0 || end_1b == 0 || start_1b >= end_1b {
            continue;
        }

        // Convert 1-based (inclusive) to 0-based (exclusive)
        let start_0 = start_1b.saturating_sub(1);
        let end_0_excl = end_1b;

        // Extract transcript_id and gene_id from attributes
        let attrs = cols[8];
        let transcript_id = extract_gtf_attr(attrs, "transcript_id")
            .unwrap_or_else(|| format!("transcript_{}", line.len()));
        let gene_id = extract_gtf_attr(attrs, "gene_id");

        let entry = by_transcript.entry(transcript_id.clone())
            .or_insert_with(|| TemplateTranscript {
                transcript_id: transcript_id.clone(),
                gene_id: gene_id.clone(),
                chrom: chrom.clone(),
                strand,
                exons: Vec::new(),
            });

        entry.exons.push((start_0, end_0_excl));
    }

    let mut transcripts: Vec<TemplateTranscript> = by_transcript.into_values().collect();

    // Sort exons within each transcript
    for tx in &mut transcripts {
        tx.exons.sort_by_key(|e| e.0);
        // Remove duplicates
        tx.exons.dedup();
    }

    // Filter out single-position or malformed transcripts
    transcripts.retain(|tx| !tx.exons.is_empty() && tx.exons.iter().all(|(s, e)| s < e));

    Ok(transcripts)
}

/// Extract a value from GTF attributes (key="value" format).
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

// ── Family Grouping ────────────────────────────────────────────────────────

/// Group transcripts into families using the specified strategy.
pub fn group_transcripts_into_families(
    transcripts: Vec<TemplateTranscript>,
    strategy: FamilyGroupingStrategy,
) -> Vec<TemplateFamily> {
    match strategy {
        FamilyGroupingStrategy::ByGene => group_by_gene(transcripts),
        FamilyGroupingStrategy::ByOverlap => group_by_overlap(transcripts),
    }
}

/// Group transcripts by gene_id. Fall back to overlapping regions if gene_id missing.
fn group_by_gene(transcripts: Vec<TemplateTranscript>) -> Vec<TemplateFamily> {
    let mut by_gene: HashMap<(String, Option<String>, char), Vec<TemplateTranscript>> = HashMap::new();

    for tx in transcripts {
        let key = (tx.chrom.clone(), tx.gene_id.clone(), tx.strand);
        by_gene.entry(key).or_insert_with(Vec::new).push(tx);
    }

    let mut families = Vec::new();
    let mut family_id = 0;

    for ((chrom, gene_id, strand), mut txs) in by_gene {
        if txs.is_empty() {
            continue;
        }

        // Compute range
        let start = txs.iter().flat_map(|tx| tx.exons.iter().map(|e| e.0)).min().unwrap_or(0);
        let end = txs.iter().flat_map(|tx| tx.exons.iter().map(|e| e.1)).max().unwrap_or(0);

        txs.sort_by_key(|tx| tx.transcript_id.clone());

        families.push(TemplateFamily {
            family_id,
            gene_id,
            chrom,
            strand,
            transcripts: txs,
            start,
            end,
        });

        family_id += 1;
    }

    families.sort_by_key(|f| (f.chrom.clone(), f.start));
    families
}

/// Group transcripts by genomic overlap.
fn group_by_overlap(transcripts: Vec<TemplateTranscript>) -> Vec<TemplateFamily> {
    // Simple approach: sort by position, group contiguous/overlapping regions
    let mut sorted = transcripts;
    sorted.sort_by_key(|tx| {
        let first_exon = tx.exons.first().map(|e| e.0).unwrap_or(0);
        (tx.chrom.clone(), first_exon)
    });

    let mut families = Vec::new();
    let mut family_id = 0;
    let mut current_group: Vec<TemplateTranscript> = Vec::new();
    let mut current_max_end: u64 = 0;
    const OVERLAP_THRESHOLD: u64 = 100; // Min gap for new group

    for tx in sorted {
        let tx_start = tx.exons.first().map(|e| e.0).unwrap_or(0);
        let tx_end = tx.exons.last().map(|e| e.1).unwrap_or(0);

        if current_group.is_empty() {
            current_group.push(tx);
            current_max_end = tx_end;
        } else {
            let prev_chrom = &current_group[0].chrom;
            let prev_strand = current_group[0].strand;

            // If different chrom/strand or gap too large, start new family
            if &tx.chrom != prev_chrom || tx.strand != prev_strand || tx_start > current_max_end + OVERLAP_THRESHOLD {
                if !current_group.is_empty() {
                    let start = current_group.iter()
                        .flat_map(|t| t.exons.iter().map(|e| e.0))
                        .min()
                        .unwrap_or(0);
                    let end = current_group.iter()
                        .flat_map(|t| t.exons.iter().map(|e| e.1))
                        .max()
                        .unwrap_or(0);

                    families.push(TemplateFamily {
                        family_id,
                        gene_id: current_group[0].gene_id.clone(),
                        chrom: current_group[0].chrom.clone(),
                        strand: current_group[0].strand,
                        transcripts: current_group.clone(),
                        start,
                        end,
                    });
                    family_id += 1;
                }
                current_group.clear();
                current_max_end = 0;
            }

            current_group.push(tx);
            current_max_end = current_max_end.max(tx_end);
        }
    }

    // Finalize last group
    if !current_group.is_empty() {
        let start = current_group.iter()
            .flat_map(|t| t.exons.iter().map(|e| e.0))
            .min()
            .unwrap_or(0);
        let end = current_group.iter()
            .flat_map(|t| t.exons.iter().map(|e| e.1))
            .max()
            .unwrap_or(0);

        families.push(TemplateFamily {
            family_id,
            gene_id: current_group[0].gene_id.clone(),
            chrom: current_group[0].chrom.clone(),
            strand: current_group[0].strand,
            transcripts: current_group,
            start,
            end,
        });
    }

    families
}

// ── Synthetic Bundle Creation ──────────────────────────────────────────────

/// Create family groups from template families by linking overlapping bundles.
///
/// For each template family:
/// - Find all bundles that overlap the family's genomic range
/// - Link those bundles into a FamilyGroup
/// - The FamilyGroup will be enriched by real multi-mapped reads in the pipeline
pub fn create_family_groups_from_templates(
    families: Vec<TemplateFamily>,
    bundles: &[Bundle],
    _genome: Option<&GenomeIndex>,
) -> Result<Vec<FamilyGroup>> {
    let mut family_groups = Vec::new();

    for template_fam in families {
        let mut bundle_indices = Vec::new();

        // Find all bundles overlapping this family's genomic range
        for (bundle_idx, bundle) in bundles.iter().enumerate() {
            if bundle.chrom != template_fam.chrom || bundle.strand != template_fam.strand {
                continue;
            }

            // Check if bundle overlaps family range (half-open: [start, end))
            if bundle.end > template_fam.start && bundle.start < template_fam.end {
                bundle_indices.push(bundle_idx);
            }
        }

        // Only create a family group if we found overlapping bundles
        if !bundle_indices.is_empty() {
            let family_group = FamilyGroup {
                family_id: template_fam.family_id,
                bundle_indices,
                multimap_reads: HashMap::new(), // Will be populated by multi-mappers in pipeline
            };
            family_groups.push(family_group);
        }
    }

    Ok(family_groups)
}

// ── Public API ────────────────────────────────────────────────────────────

/// Full pipeline: parse GTF → group → create families → link to bundles.
pub fn ingest_gtf_families(
    gtf_path: &Path,
    grouping_strategy: FamilyGroupingStrategy,
    bundles: &[Bundle],
    genome: Option<&GenomeIndex>,
) -> Result<Vec<FamilyGroup>> {
    eprintln!("[GTF-Ingestion] Parsing GTF: {:?}", gtf_path);
    let transcripts = parse_gtf_transcripts(gtf_path)?;
    eprintln!("[GTF-Ingestion] Parsed {} transcripts", transcripts.len());

    eprintln!("[GTF-Ingestion] Grouping into families ({:?})...", grouping_strategy);
    let families = group_transcripts_into_families(transcripts, grouping_strategy);
    eprintln!("[GTF-Ingestion] Created {} families", families.len());

    for (i, fam) in families.iter().enumerate() {
        eprintln!(
            "  Family {}: {} transcripts, {}:{}-{} ({})",
            i, fam.transcripts.len(), fam.chrom, fam.start, fam.end, fam.strand
        );
    }

    eprintln!("[GTF-Ingestion] Linking families to bundles...");
    let family_groups = create_family_groups_from_templates(families, bundles, genome)?;
    eprintln!(
        "[GTF-Ingestion] Created {} family groups from {} bundles",
        family_groups.len(),
        bundles.len()
    );

    Ok(family_groups)
}

// ── Family-Based Naming for Discovered Transcripts ────────────────────────

/// Map of annotated families keyed by (chrom, strand) for fast lookup.
pub type FamilyIndex = HashMap<(String, char), Vec<TemplateFamily>>;

/// Build an index of families from template transcripts for fast lookup.
pub fn build_family_index(families: Vec<TemplateFamily>) -> FamilyIndex {
    let mut idx: FamilyIndex = HashMap::new();
    for fam in families {
        let key = (fam.chrom.clone(), fam.strand);
        idx.entry(key).or_insert_with(Vec::new).push(fam);
    }
    idx
}

/// Check if a discovered transcript's junctions match a family template.
/// Returns (family_name, match_quality) if found.
fn match_transcript_to_family(
    tx_exons: &[(u64, u64)],
    families: &[TemplateFamily],
) -> Option<(String, f64)> {
    if tx_exons.is_empty() {
        return None;
    }

    let mut best_match: Option<(String, f64)> = None;

    for fam in families {
        // Score based on junction overlap with family transcripts
        for template_tx in &fam.transcripts {
            // Build junction set from discovered transcript
            let tx_junctions: std::collections::HashSet<_> = (0..tx_exons.len() - 1)
                .map(|i| (tx_exons[i].1, tx_exons[i + 1].0))
                .collect();

            // Build junction set from template
            let template_junctions: std::collections::HashSet<_> = (0..template_tx.exons.len() - 1)
                .map(|i| (template_tx.exons[i].1, template_tx.exons[i + 1].0))
                .collect();

            // Count matches
            let matching_junctions = tx_junctions.intersection(&template_junctions).count();
            let total_junctions = tx_junctions.union(&template_junctions).count();

            if total_junctions > 0 {
                let jaccard = matching_junctions as f64 / total_junctions as f64;
                // Lower threshold to catch partial isoforms and variants
                if jaccard > 0.3 {
                    // Match found: use gene_id if available, else family_id
                    let family_name = fam
                        .gene_id
                        .as_ref()
                        .cloned()
                        .unwrap_or_else(|| format!("FAMILY_{}", fam.family_id));

                    if best_match.is_none() || jaccard > best_match.as_ref().unwrap().1 {
                        best_match = Some((family_name, jaccard));
                    }
                }
            }
        }
    }

    best_match
}

/// Apply family-based naming to discovered transcripts.
/// Renames transcripts that match annotated families with a "-like" suffix.
pub fn apply_family_based_naming(
    transcripts: &mut [crate::path_extract::Transcript],
    gtf_path: &Path,
    grouping_strategy: FamilyGroupingStrategy,
) -> Result<()> {
    eprintln!("[Family-Naming] Loading reference families from: {:?}", gtf_path);

    // Parse reference transcripts
    let ref_transcripts = parse_gtf_transcripts(gtf_path)?;
    eprintln!("[Family-Naming] Parsed {} reference transcripts", ref_transcripts.len());
    let families = group_transcripts_into_families(ref_transcripts, grouping_strategy);
    let family_index = build_family_index(families.clone());

    eprintln!("[Family-Naming] Loaded {} families", family_index.len());
    for ((chrom, strand), fams) in family_index.iter() {
        eprintln!("  - {}:{}: {} families", chrom, strand, fams.len());
    }

    let mut renamed = 0;
    let mut isoform_counters: HashMap<String, usize> = HashMap::new();
    let mut checked = 0;
    let mut matched_families = HashMap::new();

    for tx in transcripts.iter_mut() {
        // Skip transcripts that already have forced IDs
        if tx.transcript_id.is_some() || tx.gene_id.is_some() {
            continue;
        }

        checked += 1;

        // Look for families on the same chromosome and strand
        let key = (tx.chrom.clone(), tx.strand);
        if let Some(families_on_locus) = family_index.get(&key) {
            if let Some((family_name, quality)) = match_transcript_to_family(&tx.exons, families_on_locus)
            {
                // Assign family-based name with isoform number
                let isoform_num = isoform_counters
                    .entry(family_name.clone())
                    .and_modify(|c| *c += 1)
                    .or_insert(1);

                let new_gene_id = format!("{}-like", family_name);
                let new_tx_id = format!("{}-like.{}", family_name, isoform_num);

                eprintln!(
                    "[Family-Naming] Renamed transcript to {} (Jaccard = {:.3})",
                    new_tx_id, quality
                );

                tx.gene_id = Some(new_gene_id);
                tx.transcript_id = Some(new_tx_id);
                renamed += 1;

                *matched_families.entry(family_name).or_insert(0) += 1;
            }
        } else if checked <= 10 {
            // Log first few chromosomes with no families for debugging
            eprintln!("[Family-Naming] No families found on {}:{}", tx.chrom, tx.strand);
        }
    }

    eprintln!("[Family-Naming] Checked {} discovered transcripts, renamed {} based on family membership", checked, renamed);
    for (fam_name, count) in matched_families.iter() {
        eprintln!("  - {}: {} transcripts", fam_name, count);
    }
    Ok(())
}
