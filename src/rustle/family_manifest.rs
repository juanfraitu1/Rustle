//! Family manifest ingestion: parse R-exported TSV of multi-copy gene family loci
//! and create FamilyGroups that the existing EM pipeline can consume.

use crate::vg::FamilyGroup;
use crate::types::Bundle;
use crate::util::coord::overlaps_half_open;
use anyhow::{anyhow, Context, Result};
use std::collections::HashMap;
use std::path::Path;

#[derive(Debug, Clone)]
pub struct FamilyLocus {
    pub family_id: String,
    pub gene_id:   String,
    pub chrom:     String,
    pub start:     u64,
    pub end:       u64,
    pub strand:    char,
}

/// Parse a tab-separated family manifest produced by `export_family_manifest()` in R.
///
/// Format (tab-separated, first row is a header):
/// ```text
/// family_id<TAB>gene_id<TAB>chrom<TAB>start<TAB>end<TAB>strand
/// RBMY<TAB>LOC129530259<TAB>chrY<TAB>1000000<TAB>1010000<TAB>+
/// ```
/// Blank lines and lines starting with `#` are ignored.
pub fn parse_family_manifest<P: AsRef<Path>>(path: P) -> Result<Vec<FamilyLocus>> {
    let content = std::fs::read_to_string(path.as_ref())
        .with_context(|| format!("reading family manifest {:?}", path.as_ref()))?;

    let mut loci = Vec::new();
    let mut header_seen = false;

    for (lineno, line) in content.lines().enumerate() {
        let line = line.trim();
        if line.is_empty() || line.starts_with('#') { continue; }

        if !header_seen {
            header_seen = true;
            continue; // skip header row
        }

        let cols: Vec<&str> = line.split('\t').collect();
        if cols.len() != 6 {
            return Err(anyhow!(
                "family manifest line {}: expected exactly 6 tab-separated columns, got {}: {:?}",
                lineno + 1, cols.len(), line
            ));
        }

        let start: u64 = cols[3].parse()
            .with_context(|| format!("manifest line {}: bad start {:?}", lineno + 1, cols[3]))?;
        let end: u64 = cols[4].parse()
            .with_context(|| format!("manifest line {}: bad end {:?}", lineno + 1, cols[4]))?;
        let strand_str = cols[5].trim();
        let strand_char = match strand_str {
            "+" | "-" | "." => strand_str.chars().next().unwrap(),
            other => return Err(anyhow!(
                "manifest line {}: invalid strand {:?}", lineno + 1, other
            )),
        };

        loci.push(FamilyLocus {
            family_id: cols[0].to_string(),
            gene_id:   cols[1].to_string(),
            chrom:     cols[2].to_string(),
            start,
            end,
            strand:    strand_char,
        });
    }

    Ok(loci)
}

/// Build `FamilyGroup`s from the manifest loci by matching overlapping bundles.
///
/// For each distinct `family_id`, collects every bundle whose `[start, end)` half-open
/// range overlaps any locus belonging to that family. The resulting `FamilyGroup`
/// has `multimap_reads` left empty — the EM pipeline populates it from the BAM.
///
/// A bundle may appear in more than one family group if it overlaps loci from multiple
/// distinct families. Each family's EM runs independently, so this does not cause
/// double-counting within a family.
pub fn create_family_groups_from_manifest(
    loci: &[FamilyLocus],
    bundles: &[Bundle],
) -> Vec<FamilyGroup> {
    // Group loci by family_id (preserve insertion order for determinism).
    let mut family_order: Vec<String> = Vec::new();
    let mut family_loci: HashMap<String, Vec<&FamilyLocus>> = HashMap::new();
    for loc in loci {
        if !family_loci.contains_key(&loc.family_id) {
            family_order.push(loc.family_id.clone());
        }
        family_loci.entry(loc.family_id.clone()).or_default().push(loc);
    }

    let mut groups = Vec::new();
    for (family_numeric_id, fid) in family_order.iter().enumerate() {
        let fam_loci = &family_loci[fid];
        let mut bundle_indices: Vec<usize> = Vec::new();

        for (bi, bundle) in bundles.iter().enumerate() {
            // Strand is intentionally not checked: bundles are already strand-separated
            // by the bundler, so chrom + coordinate overlap is sufficient.
            let overlaps_any = fam_loci.iter().any(|loc| {
                bundle.chrom == loc.chrom
                    && overlaps_half_open(bundle.start, bundle.end, loc.start, loc.end)
            });
            if overlaps_any {
                bundle_indices.push(bi);
            }
        }

        if !bundle_indices.is_empty() {
            groups.push(FamilyGroup {
                family_id: family_numeric_id,
                bundle_indices,
                multimap_reads: HashMap::new(),
            });
        }
    }

    groups
}
