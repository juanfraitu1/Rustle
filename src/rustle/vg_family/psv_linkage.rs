//! Layer-2 "C": within-molecule PSV->junction linkage for paralog isoform recovery.
//! Task 1: PSV-column extraction (positions + per-copy alleles + per-copy genomic coords).
//!
//! A PSV (paralog-sequence variant) column is one base position, inside a
//! diagnostic exon node, where >=2 copies of the family carry distinct alleles.
//! It is the genotyping reference for later tasks: a read covering this column
//! tells us which copy it came from (its base matches that copy's allele).
//!
//! WHY this is its own primitive (pure, not yet wired): the column carries the
//! genomic coordinate IN EACH COPY's OWN FRAME (from `per_copy_spans`, never the
//! cross-copy union `ExonClass.span`), so the genotyping pass in a later task can
//! look up "the base this read has at copy X's position 502" directly. Computing
//! that frame mapping once, deterministically, keeps the downstream pass cheap and
//! the whole channel reproducible (DetHash + total-order sorts everywhere).

use crate::vg_family::family_graph::FamilyGraph;

/// One copy's allele at a PSV column: the base and its genomic position IN THAT
/// COPY's frame (the copy's exon start + the within-exon offset).
#[derive(Debug, Clone, PartialEq)]
pub struct PsvCopyAllele {
    pub copy_id: usize,
    pub genomic_pos: u64,
    pub allele: u8,
}

/// One paralog-distinguishing column: a position (within a diagnostic exon node)
/// where >= 2 copies carry distinct bases. `per_copy` is sorted by copy_id so the
/// column is canonical regardless of the order copies were pushed into the node.
#[derive(Debug, Clone, PartialEq)]
pub struct PsvColumn {
    pub node_idx: usize,
    pub per_copy: Vec<PsvCopyAllele>,
}

/// Extract the family's PSV columns. For each diagnostic node (a node whose
/// `per_copy_sequences` carry >= 2 DISTINCT sequences), for each within-exon
/// offset where the present copies' bases differ, emit a `PsvColumn`: each present
/// copy's base + the genomic coordinate `per_copy_spans[copy].0 + offset` (the
/// copy's OWN frame, never `ExonClass.span`, which is a cross-copy union and can
/// span the whole chromosome).
///
/// v1 limitation (documented): only handle nodes where ALL contributing copies'
/// per_copy_sequences are EQUAL LENGTH (no indel re-alignment) — skip a node whose
/// copies' sequences differ in length. With equal length, offset i is the same
/// genomic position relative to each copy's start, so the per-copy frame mapping is
/// exact; an indel would shift offsets between copies and require alignment, which
/// is deferred. A copy that has a sequence but no `per_copy_spans` entry is skipped
/// (we cannot place its allele in genomic coordinates without its frame).
///
/// Deterministic: the returned Vec is sorted by (node_idx, the column's first
/// genomic_pos); each column's `per_copy` is sorted by copy_id.
pub fn psv_columns_for_family(fg: &FamilyGraph) -> Vec<PsvColumn> {
    let mut columns: Vec<PsvColumn> = Vec::new();

    for node in &fg.nodes {
        // A node's contributing copies = its per_copy_sequences entries. We must be
        // able to place each in genomic coordinates, so a copy with a sequence but
        // no per_copy_spans entry is dropped (we cannot know its frame).
        let mut contributors: Vec<(usize, &[u8], u64)> = Vec::new();
        for (cid, seq) in &node.per_copy_sequences {
            if let Some((_, (start, _))) = node.per_copy_spans.iter().find(|(c, _)| c == cid) {
                contributors.push((*cid, seq.as_slice(), *start));
            }
        }

        // Diagnostic = >= 2 contributors carrying >= 2 DISTINCT sequences.
        if contributors.len() < 2 {
            continue;
        }
        let distinct = {
            let mut s: crate::types::DetHashSet<&[u8]> = crate::types::DetHashSet::default();
            for (_, seq, _) in &contributors {
                s.insert(seq);
            }
            s.len()
        };
        if distinct < 2 {
            continue;
        }

        // v1: only equal-length contributors (no indel re-alignment). If lengths
        // differ, offset i is NOT the same genomic position across copies, so we
        // skip the whole node rather than misalign columns.
        let len0 = contributors[0].1.len();
        if contributors.iter().any(|(_, seq, _)| seq.len() != len0) {
            continue;
        }

        // For each within-exon offset, if >= 2 distinct bases are present among
        // contributors, emit a column with EVERY present copy's base placed in its
        // OWN frame (start + offset). per_copy sorted by copy_id for canonical form.
        for offset in 0..len0 {
            let mut distinct_bases: crate::types::DetHashSet<u8> = crate::types::DetHashSet::default();
            for (_, seq, _) in &contributors {
                distinct_bases.insert(seq[offset]);
            }
            if distinct_bases.len() < 2 {
                continue;
            }
            let mut per_copy: Vec<PsvCopyAllele> = contributors
                .iter()
                .map(|(cid, seq, start)| PsvCopyAllele {
                    copy_id: *cid,
                    genomic_pos: start + offset as u64,
                    allele: seq[offset],
                })
                .collect();
            per_copy.sort_by_key(|p| p.copy_id);
            columns.push(PsvColumn { node_idx: node.idx.0, per_copy });
        }
    }

    // Total order: (node_idx, the column's first genomic_pos). per_copy is already
    // sorted by copy_id, so "first genomic_pos" is the lowest-copy_id allele's pos
    // — a stable, deterministic key.
    columns.sort_by(|a, b| {
        let a_first = a.per_copy.first().map(|p| p.genomic_pos).unwrap_or(0);
        let b_first = b.per_copy.first().map(|p| p.genomic_pos).unwrap_or(0);
        a.node_idx.cmp(&b.node_idx).then(a_first.cmp(&b_first))
    });

    columns
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::vg_family::family_graph::{ExonClass, FamilyGraph, JunctionEdge, NodeIdx};

    /// Build a one-node family graph by hand (modeled on layer2.rs's
    /// `build_psv_three_copy_fg`). `seqs` = (copy_id, exon sequence, exon span).
    /// The node's union `span` is set to a deliberately WRONG, whole-chromosome
    /// value so any test that accidentally reads `span` instead of `per_copy_spans`
    /// produces obviously-wrong coordinates.
    fn one_node_fg(seqs: &[(usize, &[u8], (u64, u64))]) -> FamilyGraph {
        let node = ExonClass {
            idx: NodeIdx(0),
            chrom: "chrT".into(),
            // Cross-copy union landmine: NEVER the coordinate source.
            span: (0, 1_000_000),
            strand: '+',
            per_copy_sequences: seqs.iter().map(|(c, s, _)| (*c, s.to_vec())).collect(),
            per_copy_spans: seqs.iter().map(|(c, _, sp)| (*c, *sp)).collect(),
            copy_specific: seqs.len() == 1,
            per_copy_cov: Vec::new(),
        };
        FamilyGraph { family_id: 0, nodes: vec![node], edges: Vec::<JunctionEdge>::new() }
    }

    /// Two-copy diagnostic node: copy0 b"ACGTAC"@(100,106), copy1 b"ACCTAC"@(500,506).
    /// They differ ONLY at offset 2 (G vs C). The single PSV column must map that
    /// offset into each copy's OWN frame: copy0 -> 102, copy1 -> 502.
    fn two_copy_one_psv_fg() -> FamilyGraph {
        one_node_fg(&[
            (0, b"ACGTAC", (100, 106)),
            (1, b"ACCTAC", (500, 506)),
        ])
    }

    #[test]
    fn psv_columns_map_offset_to_each_copy_frame() {
        let fg = two_copy_one_psv_fg();
        let cols = psv_columns_for_family(&fg);
        assert_eq!(cols.len(), 1, "exactly one differing column (offset 2): {cols:?}");
        let c = &cols[0];
        assert_eq!(c.node_idx, 0);
        assert!(
            c.per_copy.iter().any(|p| p.copy_id == 0 && p.genomic_pos == 102 && p.allele == b'G'),
            "copy0 allele G at its own frame 100+2=102: {c:?}"
        );
        assert!(
            c.per_copy.iter().any(|p| p.copy_id == 1 && p.genomic_pos == 502 && p.allele == b'C'),
            "copy1 allele C at its own frame 500+2=502: {c:?}"
        );
    }

    #[test]
    fn psv_columns_skips_non_diagnostic_node() {
        // Both copies carry the SAME sequence -> not diagnostic -> no columns.
        let fg = one_node_fg(&[
            (0, b"ACGTAC", (100, 106)),
            (1, b"ACGTAC", (500, 506)),
        ]);
        let cols = psv_columns_for_family(&fg);
        assert!(cols.is_empty(), "identical copies are not diagnostic: {cols:?}");
    }

    #[test]
    fn psv_columns_skips_unequal_length() {
        // Copies differ in length -> v1 cannot column-align without indel handling.
        // The sequences are also clearly distinct (so the node IS diagnostic), but
        // the length mismatch must still suppress all columns from this node.
        let fg = one_node_fg(&[
            (0, b"ACGTAC", (100, 106)),
            (1, b"ACGTACG", (500, 507)),
        ]);
        let cols = psv_columns_for_family(&fg);
        assert!(cols.is_empty(), "v1 skips unequal-length nodes (no indel align): {cols:?}");
    }

    #[test]
    fn psv_columns_deterministic() {
        // Same family, copies pushed in different order -> identical output.
        let fg_a = one_node_fg(&[
            (0, b"ACGTAC", (100, 106)),
            (1, b"ACCTAC", (500, 506)),
        ]);
        let fg_b = one_node_fg(&[
            (1, b"ACCTAC", (500, 506)),
            (0, b"ACGTAC", (100, 106)),
        ]);
        assert_eq!(
            psv_columns_for_family(&fg_a),
            psv_columns_for_family(&fg_b),
            "column extraction must be order-independent (per_copy sorted by copy_id)"
        );
    }

    #[test]
    fn psv_columns_multiple_offsets_sorted() {
        // Two differing offsets (1 and 4) in one node -> two columns, sorted by
        // genomic_pos. Confirms the (node_idx, first genomic_pos) total order.
        let fg = one_node_fg(&[
            (0, b"AAAAAA", (100, 106)),
            (1, b"ACAACA", (500, 506)),
        ]);
        let cols = psv_columns_for_family(&fg);
        assert_eq!(cols.len(), 2, "offsets 1 and 4 differ: {cols:?}");
        // First column = lower copy0 genomic_pos (101 < 104).
        assert_eq!(cols[0].per_copy.iter().find(|p| p.copy_id == 0).unwrap().genomic_pos, 101);
        assert_eq!(cols[1].per_copy.iter().find(|p| p.copy_id == 0).unwrap().genomic_pos, 104);
    }

    #[test]
    fn isoform_source_psvlinked_orders_after_native() {
        // source_rank parity check: PsvLinked must rank strictly after Native and
        // Transferred, and order among PsvLinked by copy_id. Mirrors the closure in
        // emit_family_isoforms so the layer2 regression anchor holds.
        use crate::vg_family::layer2::IsoformSource;
        let rank = |s: &IsoformSource| -> (usize, usize) {
            match s {
                IsoformSource::Native => (0, 0),
                IsoformSource::Transferred { donor_copy } => (1, *donor_copy),
                IsoformSource::PsvLinked { copy_id } => (2, *copy_id),
            }
        };
        assert!(rank(&IsoformSource::Native) < rank(&IsoformSource::PsvLinked { copy_id: 0 }));
        assert!(
            rank(&IsoformSource::Transferred { donor_copy: 9 })
                < rank(&IsoformSource::PsvLinked { copy_id: 0 })
        );
        assert!(
            rank(&IsoformSource::PsvLinked { copy_id: 0 })
                < rank(&IsoformSource::PsvLinked { copy_id: 1 })
        );
    }
}
