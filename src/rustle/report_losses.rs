//! Parse gffcompare annotation to report where reference transcripts are lost.
//! Reads ref GTF (transcript IDs) and gffcompare tracking file; reports missed and best class per ref.

use crate::types::DetHashMap as HashMap;
use anyhow::{Context, Result};
use std::path::Path;

use crate::reference_gtf::parse_reference_gtf;

/// Best class code for a reference transcript from tracking (= exact, j alt splice, x/c/m/etc partial, - no ref).
fn class_rank(c: u8) -> u8 {
    match c {
        b'=' => 0,  // exact match
        b'j' => 1,  // alternative splice (count as recovered)
        b'c' => 2,  // contained
        b'm' => 3,  // retained intron / merge
        b'x' => 4,  // exonic overlap
        b'i' => 5,  // intron
        b'y' => 6,  // same strand
        b'p' => 7,  // polymerase run-on
        b'n' => 8,  // generic
        b's' => 9,  // opposite strand
        b'u' => 10, // unknown (no ref)
        _ => 11,
    }
}

/// Parse gffcompare tracking file: columns are tab-sep; col 3 = ref (gene|transcript or -), col 4 = class_code.
/// Returns map: ref_transcript_id -> best class code (single char).
pub fn parse_tracking<P: AsRef<Path>>(path: P) -> Result<HashMap<String, u8>> {
    let s = std::fs::read_to_string(path.as_ref())
        .with_context(|| format!("read tracking {}", path.as_ref().display()))?;
    let mut best: HashMap<String, u8> = Default::default();
    for line in s.lines() {
        let mut cols = line.split('\t');
        let _q = cols.next();
        let _xloc = cols.next();
        let ref_part = cols.next().unwrap_or("-");
        let class = cols.next().and_then(|c| c.bytes().next()).unwrap_or(b'-');
        let ref_tx_id = if ref_part == "-" || ref_part.is_empty() {
            continue;
        } else {
            ref_part.rsplit('|').next().unwrap_or(ref_part).to_string()
        };
        let entry = best.entry(ref_tx_id).or_insert(class);
        if class_rank(class) < class_rank(*entry) {
            *entry = class;
        }
    }
    Ok(best)
}

/// Status of a reference transcript relative to query assembly.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum LossStatus {
    /// Exact or alternative splice match (= or j)
    Matched,
    /// Partial (x, c, m, etc.)
    Partial(u8),
    /// Not found in tracking (no query overlapped this ref)
    Missed,
}

/// Run report: load ref GTF and tracking, print counts and list of missed/partial refs.
pub fn report_losses<P: AsRef<Path>>(
    ref_gtf_path: P,
    tracking_path: P,
    verbose: bool,
) -> Result<()> {
    let ref_txs = parse_reference_gtf(ref_gtf_path)?;
    let ref_ids: Vec<String> = ref_txs.iter().map(|t| t.id.clone()).collect();
    let tracking = parse_tracking(tracking_path)?;

    let mut matched = 0_usize;
    let mut partial = 0_usize;
    let mut missed = 0_usize;
    let mut partial_by_class: HashMap<u8, usize> = Default::default();
    let mut missed_ids: Vec<String> = Vec::new();
    let mut partial_list: Vec<(String, u8)> = Vec::new();

    for id in &ref_ids {
        let status = match tracking.get(id) {
            Some(&c) if c == b'=' || c == b'j' => {
                matched += 1;
                LossStatus::Matched
            }
            Some(&c) => {
                partial += 1;
                *partial_by_class.entry(c).or_insert(0) += 1;
                if verbose {
                    partial_list.push((id.clone(), c));
                }
                LossStatus::Partial(c)
            }
            None => {
                missed += 1;
                missed_ids.push(id.clone());
                LossStatus::Missed
            }
        };
        if verbose && matches!(status, LossStatus::Missed) {
            // missed_ids already pushed
        }
    }

    eprintln!("# Reference transcripts: {} total", ref_ids.len());
    eprintln!(
        "# Matched (= or j): {} ({:.1}%)",
        matched,
        100.0 * matched as f64 / ref_ids.len().max(1) as f64
    );
    eprintln!(
        "# Partial (x/c/m/...): {} ({:.1}%)",
        partial,
        100.0 * partial as f64 / ref_ids.len().max(1) as f64
    );
    eprintln!(
        "# Missed (no query): {} ({:.1}%)",
        missed,
        100.0 * missed as f64 / ref_ids.len().max(1) as f64
    );
    if !partial_by_class.is_empty() {
        let class_summary: Vec<String> = partial_by_class
            .iter()
            .map(|(&c, &n)| format!("{}:{}", c as char, n))
            .collect();
        eprintln!("# Partial by class: {}", class_summary.join(", "));
    }
    if !missed_ids.is_empty() {
        eprintln!("# Missed transcript IDs (first 100):");
        for id in missed_ids.iter().take(100) {
            eprintln!("  {}", id);
        }
        if missed_ids.len() > 100 {
            eprintln!("  ... and {} more", missed_ids.len() - 100);
        }
    }
    if verbose && !partial_list.is_empty() {
        eprintln!("# Partial transcript IDs (first 50):");
        for (id, c) in partial_list.iter().take(50) {
            eprintln!("  {}  class={}", id, *c as char);
        }
    }
    if missed > 0 {
        eprintln!("# Next: run with --trace-reference <ref.gtf> to see why missed (NotExtracted vs filter stage)");
    }

    Ok(())
}
