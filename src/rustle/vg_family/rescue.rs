//! Family helpers shared by the VG novel-copy paths.
//!
//! The HMM-based rescue driver (profile fitting + forward/Viterbi scoring +
//! synthetic-bundle emission) was removed; the surviving piece is the
//! per-strand family partitioning used by both the k-mer rescue and the
//! family-graph identity diagnostics.

use crate::types::Bundle;
use crate::vg::FamilyGroup;

/// Split a mixed-strand FamilyGroup into per-strand sub-families.
///
/// Many discovered families contain bundles on both `+` and `-` strands —
/// shared multi-mappers across strands are common (sequence palindromy,
/// antisense paralogs, or cross-strand seed matches). `build_family_graph`
/// requires single-strand input, so we partition by strand here. Sub-families
/// with fewer than 2 bundles are dropped (a 1-bundle "family" produces no
/// useful signal).
///
/// Returns the original family if it's already single-strand; otherwise one
/// new FamilyGroup per strand. Each split family inherits a derived
/// `family_id` (`original * 10 + 0/1/2`) so the source family is recoverable
/// from the report.
pub fn partition_family_by_strand(
    family: &FamilyGroup,
    bundles: &[Bundle],
) -> Vec<FamilyGroup> {
    use std::collections::BTreeMap;
    let mut by_strand: BTreeMap<char, Vec<usize>> = BTreeMap::new();
    for &bi in &family.bundle_indices {
        by_strand
            .entry(bundles[bi].strand)
            .or_default()
            .push(bi);
    }
    if by_strand.len() <= 1 {
        return vec![FamilyGroup {
            family_id: family.family_id,
            bundle_indices: family.bundle_indices.clone(),
            multimap_reads: family.multimap_reads.clone(),
        }];
    }
    let mut out = Vec::new();
    for (i, (_strand, bis)) in by_strand.into_iter().enumerate() {
        if bis.len() < 2 {
            continue;
        }
        out.push(FamilyGroup {
            family_id: family.family_id * 10 + i,
            bundle_indices: bis,
            // Multimap reads are kept as-is; downstream consumers filter by bundle_indices.
            multimap_reads: family.multimap_reads.clone(),
        });
    }
    out
}
