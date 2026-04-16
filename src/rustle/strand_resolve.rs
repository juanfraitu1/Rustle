//! Unstranded strand resolution.
//! Assigns strand to unstranded multi-exon reads using junction evidence and coverage.

use crate::bpcov::{BpcovStranded, BPCOV_STRAND_MINUS, BPCOV_STRAND_PLUS};
use crate::types::{BundleRead, Junction, JunctionStats};

/// Resolve strand for unstranded multi-exon reads.
///
/// Algorithm:
/// 1. For each unstranded multi-exon read, check junction strand evidence
/// 2. If junctions don't resolve: use plus vs minus coverage ratio
/// 3. Assign strand to read, propagate to junctions
///
/// Returns count of resolved reads.
pub fn resolve_unstranded_reads(
    reads: &mut [BundleRead],
    junction_stats: &mut JunctionStats,
    bpcov_stranded: &BpcovStranded,
    bundle_start: u64,
    long_reads: bool,
) -> usize {
    let mut resolved = 0;

    for read in reads.iter_mut() {
        // Only process unstranded multi-exon reads
        if read.strand != '.' || read.exons.len() < 2 {
            continue;
        }

        // Check bounds
        if read.ref_start < bundle_start {
            continue;
        }

        let mut pos_cov: f64 = 0.0;
        let mut neg_cov: f64 = 0.0;
        let mut posjunc = false;
        let mut negjunc = false;

        // Accumulate coverage and check junction strands
        for (ei, &(exon_start, exon_end)) in read.exons.iter().enumerate() {
            if long_reads {
                let s = exon_start.saturating_sub(bundle_start) as usize;
                let e = exon_end.saturating_sub(bundle_start) as usize;
                pos_cov += bpcov_stranded.get_cov_range(BPCOV_STRAND_PLUS, s, e);
                neg_cov += bpcov_stranded.get_cov_range(BPCOV_STRAND_MINUS, s, e);
            }

            // Check junction strand between this exon and the next
            if ei > 0 {
                let junc = if ei - 1 < read.junctions.len() {
                    read.junctions[ei - 1]
                } else {
                    Junction::new(read.exons[ei - 1].1, exon_start)
                };

                if let Some(stat) = junction_stats.get(&junc) {
                    match stat.strand {
                        Some(s) if s > 0 => posjunc = true,
                        Some(s) if s < 0 => negjunc = true,
                        _ => {
                            // Junction has no strand or killed — check if a stranded
                            // version exists in the stats (uses binary search by
                            // creating CJunction with strand=1 then strand=-1).
                            // In Rust, junctions are keyed by (donor, acceptor) only,
                            // so strand info is already in the stat. If strand is None/0,
                            // there's no stranded evidence for this junction.
                        }
                    }
                }
            }
        }

        // Resolve strand
        let new_strand = if posjunc && !negjunc {
            '+'
        } else if negjunc && !posjunc {
            '-'
        } else {
            // Coverage tiebreak
            if neg_cov < 1.0 {
                neg_cov = 0.0;
            }
            if pos_cov < 1.0 {
                pos_cov = 0.0;
            }
            if neg_cov > pos_cov {
                '-'
            } else if pos_cov > neg_cov {
                '+'
            } else {
                continue; // Can't resolve
            }
        };

        read.strand = new_strand;
        resolved += 1;

        // Propagate strand to unstranded junctions
        let strand_val: i8 = if new_strand == '+' { 1 } else { -1 };
        for j in &read.junctions {
            if let Some(stat) = junction_stats.get_mut(j) {
                if stat.strand.is_none() || stat.strand == Some(0) {
                    stat.strand = Some(strand_val);
                }
            }
        }
    }

    if resolved > 0 && std::env::var_os("RUSTLE_STRAND_DIAG").is_some() {
        eprintln!("STRAND_RESOLVE: resolved {} unstranded reads", resolved);
    }

    resolved
}
