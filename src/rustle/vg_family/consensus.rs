//! Cross-copy consensus error-correction (the one SUBTRACTIVE precision lever).
//!
//! Every other cross-copy mechanism is ADDITIVE (borrow/completion/topo inject structure into
//! sparse copies → over-enumeration). This module is the opposite: it VETOES a low-coverage copy's
//! off-consensus junction when the family's well-covered siblings agree on a different backbone.
//!
//! Consensus is computed in the family graph's HOMOLOGOUS-edge space (a node-pair `(from,to)`),
//! NOT absolute coordinates — paralog copies sit at different loci, so a junction shared between
//! copies has different absolute donor/acceptor in each but maps to the SAME ExonClass node pair.
//!
//! The veto rule is deliberately conservative to avoid the indistinguishability wall (deleting a
//! REAL copy-private junction = biological divergence between paralogs): a junction is vetoed ONLY
//! when its edge is supported by `< min_consensus_copies` copies AND its own read support is
//! `< support_floor`. A well-supported copy-private junction (support ≥ floor) is preserved even if
//! no sibling shares it.

use super::family_graph::{CopyId, NodeIdx};

/// One copy's junction within a family, mapped to the homologous family-graph edge it spans.
#[derive(Clone, Debug, PartialEq)]
pub struct FamilyJunction {
    pub copy: CopyId,
    /// The homologous edge (family-graph node pair) this junction maps to.
    pub edge: (NodeIdx, NodeIdx),
    /// The copy's own absolute junction coordinates (donor, acceptor).
    pub junc: (u64, u64),
    /// Read support for this junction at this copy.
    pub support: f64,
}

/// Return the (copy, junction) pairs to VETO: junctions whose homologous edge is present in fewer
/// than `min_consensus_copies` copies AND whose own read support is below `support_floor`.
/// Consensus edges (≥ min copies) and well-supported copy-private junctions (support ≥ floor) are
/// preserved. Pure function — no graph mutation.
pub fn consensus_vetoes(
    juncs: &[FamilyJunction],
    min_consensus_copies: usize,
    support_floor: f64,
) -> Vec<(CopyId, (u64, u64))> {
    use std::collections::{HashMap, HashSet};
    // Distinct copies supporting each homologous edge.
    let mut edge_copies: HashMap<(NodeIdx, NodeIdx), HashSet<CopyId>> = HashMap::new();
    for j in juncs {
        edge_copies.entry(j.edge).or_default().insert(j.copy);
    }
    juncs
        .iter()
        .filter(|j| {
            let n_copies = edge_copies.get(&j.edge).map_or(0, |s| s.len());
            n_copies < min_consensus_copies && j.support < support_floor
        })
        .map(|j| (j.copy, j.junc))
        .collect()
}

/// A family-graph node's span for one copy, flattened from `ExonClass.per_copy_spans`:
/// `(node, copy, (start, end))`.
pub type NodeCopySpan = (NodeIdx, CopyId, (u64, u64));

/// Map a copy's junctions to homologous family-graph edges via per-copy node spans.
/// A junction `(donor, acceptor)` at copy `C` maps to edge `(F, T)` where node `F`'s span for `C`
/// ends at ~`donor` and node `T`'s span for `C` starts at ~`acceptor` (within `tol`). Junctions
/// that don't map to a homologous edge pair (e.g. the family graph never unified that exon across
/// copies — the barrier regime) are dropped, so they are never vetoed: consensus correction only
/// acts where the homologous structure is actually established.
pub fn map_junctions_to_edges(
    node_spans: &[NodeCopySpan],
    copy_junctions: &[(CopyId, (u64, u64), f64)],
    tol: u64,
) -> Vec<FamilyJunction> {
    let mut out = Vec::new();
    for &(cid, (donor, acceptor), support) in copy_junctions {
        let f = node_spans
            .iter()
            .find(|(_, c, (_, e))| *c == cid && e.abs_diff(donor) <= tol)
            .map(|(n, _, _)| *n);
        let t = node_spans
            .iter()
            .find(|(_, c, (s, _))| *c == cid && s.abs_diff(acceptor) <= tol)
            .map(|(n, _, _)| *n);
        if let (Some(f), Some(t)) = (f, t) {
            out.push(FamilyJunction { copy: cid, edge: (f, t), junc: (donor, acceptor), support });
        }
    }
    out
}

#[cfg(test)]
mod tests {
    use super::*;

    fn fj(copy: CopyId, edge: (usize, usize), junc: (u64, u64), support: f64) -> FamilyJunction {
        FamilyJunction { copy, edge: (NodeIdx(edge.0), NodeIdx(edge.1)), junc, support }
    }

    #[test]
    fn vetoes_thin_off_consensus_keeps_consensus_and_well_supported_private() {
        // Family of 3 copies. e1,e2 are the consensus backbone (present in >=2 copies).
        // Copy 2 (low-cov) additionally has: a THIN spurious off-consensus junction, and a
        // WELL-SUPPORTED copy-private junction (real biological divergence).
        let juncs = vec![
            fj(0, (0, 1), (100, 200), 30.0),
            fj(1, (0, 1), (1100, 1200), 25.0),
            fj(2, (0, 1), (2100, 2200), 3.0),   // consensus edge, copy2 thin but PROTECTED by consensus
            fj(0, (1, 2), (300, 400), 28.0),
            fj(1, (1, 2), (1300, 1400), 24.0),
            fj(2, (1, 5), (2100, 2900), 1.0),   // SPURIOUS: 1-copy edge, support 1 -> VETO
            fj(2, (2, 6), (2300, 2950), 22.0),  // copy-private REAL: 1-copy edge but support 22 -> KEEP
        ];
        let vetoes = consensus_vetoes(&juncs, 2, 5.0);
        assert_eq!(
            vetoes,
            vec![(2usize, (2100u64, 2900u64))],
            "only the thin off-consensus junction is vetoed; consensus + well-supported private survive"
        );
    }

    #[test]
    fn never_vetoes_when_all_junctions_are_consensus() {
        let juncs = vec![
            fj(0, (0, 1), (100, 200), 1.0),
            fj(1, (0, 1), (1100, 1200), 1.0),
        ];
        assert!(consensus_vetoes(&juncs, 2, 5.0).is_empty());
    }

    #[test]
    fn end_to_end_maps_then_vetoes_spurious_via_homologous_edges() {
        // 3-copy family. Shared backbone = nodes 0,1,2 (each copy at its own locus). Copy 2 (low-cov)
        // also has node 5 (a spurious off-consensus target) and node 6 (a REAL copy-private exon).
        let node_spans: Vec<NodeCopySpan> = vec![
            (NodeIdx(0), 0, (100, 200)), (NodeIdx(1), 0, (300, 400)), (NodeIdx(2), 0, (500, 600)),
            (NodeIdx(0), 1, (1100, 1200)), (NodeIdx(1), 1, (1300, 1400)), (NodeIdx(2), 1, (1500, 1600)),
            (NodeIdx(0), 2, (2100, 2200)), (NodeIdx(1), 2, (2300, 2400)), (NodeIdx(2), 2, (2500, 2600)),
            (NodeIdx(5), 2, (2900, 3000)), (NodeIdx(6), 2, (3100, 3200)),
        ];
        let copy_junctions = vec![
            (0usize, (200u64, 300u64), 30.0), (0, (400, 500), 28.0),   // copyA backbone (0->1, 1->2)
            (1usize, (1200u64, 1300u64), 25.0), (1, (1400, 1500), 24.0), // copyB backbone
            (2usize, (2200u64, 2300u64), 3.0),  // copyC backbone (0->1): consensus, thin but PROTECTED
            (2, (2400, 2900), 1.0),             // copyC SPURIOUS (1->5), support 1 -> VETO
            (2, (2600, 3100), 22.0),            // copyC copy-private REAL (2->6), support 22 -> KEEP
        ];
        let mapped = map_junctions_to_edges(&node_spans, &copy_junctions, 2);
        // every junction maps to an edge here
        assert_eq!(mapped.len(), 7, "all junctions map to homologous edges");
        let vetoes = consensus_vetoes(&mapped, 2, 5.0);
        assert_eq!(
            vetoes,
            vec![(2usize, (2400u64, 2900u64))],
            "only the thin off-consensus junction is vetoed; consensus + copy-private survive"
        );
    }

    #[test]
    fn unmappable_junction_is_never_vetoed() {
        // A junction whose acceptor matches no node span (barrier regime: exon not unified) is
        // dropped by the mapper, so it can never be vetoed.
        let node_spans: Vec<NodeCopySpan> = vec![
            (NodeIdx(0), 0, (100, 200)), (NodeIdx(0), 1, (1100, 1200)),
        ];
        let copy_junctions = vec![(0usize, (200u64, 99999u64), 1.0)]; // acceptor 99999 maps to nothing
        let mapped = map_junctions_to_edges(&node_spans, &copy_junctions, 2);
        assert!(mapped.is_empty(), "unmappable junction dropped -> never vetoed");
        assert!(consensus_vetoes(&mapped, 2, 5.0).is_empty());
    }

    #[test]
    fn high_floor_does_not_veto_well_supported_off_consensus() {
        // A thin AND a strong off-consensus junction; floor=5 vetoes only the thin one.
        let juncs = vec![
            fj(0, (0, 1), (100, 200), 30.0),
            fj(1, (0, 1), (1100, 1200), 30.0),
            fj(0, (1, 9), (300, 800), 2.0),    // thin off-consensus -> veto
            fj(1, (1, 8), (1300, 1800), 40.0), // strong off-consensus -> keep
        ];
        assert_eq!(consensus_vetoes(&juncs, 2, 5.0), vec![(0usize, (300u64, 800u64))]);
    }
}
