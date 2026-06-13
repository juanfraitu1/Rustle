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
