//! Per-bundlenode junction-expanded graph processing.
//!
//! Splits a merged bundle into junction-connected bundlenode components,
//! then builds a separate graph per component. This matches StringTie's
//! architecture where each graph contains the bundlenodes reachable via
//! junctions from a seed bundlenode, not the entire bundle.
//!
//! The key difference from the BUNDLE_GRAPH approach: instead of ONE graph
//! for the entire merged bundle, this creates multiple smaller graphs
//! (one per junction-connected component). Each graph is independently
//! flow-decomposed, preventing large bundles from overwhelming the flow.

use crate::types::{CBundlenode, Junction};

/// Find junction-connected bundlenode components.
///
/// Two bundlenodes are connected if a junction has its donor in one and
/// acceptor in the other. Returns groups of bundlenode indices (into the
/// bundlenode list) that form connected components.
pub fn find_junction_connected_components(
    bundlenodes: &[(u64, u64, f64, usize)], // (start, end, cov, bid)
    junctions: &[Junction],
) -> Vec<Vec<usize>> {
    let n = bundlenodes.len();
    if n == 0 {
        return Vec::new();
    }
    if n == 1 {
        return vec![vec![0]];
    }

    // Union-Find
    let mut parent: Vec<usize> = (0..n).collect();
    let mut rank: Vec<usize> = vec![0; n];

    fn find(parent: &mut [usize], x: usize) -> usize {
        if parent[x] != x {
            parent[x] = find(parent, parent[x]);
        }
        parent[x]
    }
    fn union(parent: &mut [usize], rank: &mut [usize], a: usize, b: usize) {
        let ra = find(parent, a);
        let rb = find(parent, b);
        if ra == rb { return; }
        if rank[ra] < rank[rb] { parent[ra] = rb; }
        else if rank[ra] > rank[rb] { parent[rb] = ra; }
        else { parent[rb] = ra; rank[ra] += 1; }
    }

    // For each junction, find which bundlenodes contain the donor and acceptor
    for j in junctions {
        let donor_bn = bundlenodes.iter().position(|&(s, e, _, _)| j.donor >= s && j.donor <= e);
        let acceptor_bn = bundlenodes.iter().position(|&(s, e, _, _)| j.acceptor >= s && j.acceptor <= e);
        if let (Some(d), Some(a)) = (donor_bn, acceptor_bn) {
            if d != a {
                union(&mut parent, &mut rank, d, a);
            }
        }
    }

    // Also connect adjacent/overlapping bundlenodes (contiguous exonic regions)
    for i in 0..n.saturating_sub(1) {
        let (_, e1, _, _) = bundlenodes[i];
        let (s2, _, _, _) = bundlenodes[i + 1];
        // If bundlenodes are contiguous or overlapping, they're in the same component
        if s2 <= e1 + 1 {
            union(&mut parent, &mut rank, i, i + 1);
        }
    }

    // Collect components
    let mut components: std::collections::HashMap<usize, Vec<usize>> = Default::default();
    for i in 0..n {
        let root = find(&mut parent, i);
        components.entry(root).or_default().push(i);
    }

    components.into_values().collect()
}

/// Build a bundlenode linked list from a subset of bundlenodes.
pub fn build_bnode_chain(
    bundlenodes: &[(u64, u64, f64, usize)],
    indices: &[usize],
) -> Option<CBundlenode> {
    if indices.is_empty() {
        return None;
    }
    let mut sorted: Vec<(u64, u64, f64, usize)> = indices
        .iter()
        .map(|&i| bundlenodes[i])
        .collect();
    sorted.sort_by_key(|&(s, _, _, _)| s);

    let mut head = CBundlenode {
        start: sorted[0].0,
        end: sorted[0].1,
        cov: sorted[0].2,
        bid: sorted[0].3,
        next: None,
        hardstart: false,
        hardend: false,
    };
    let mut tail = &mut head;
    for &(s, e, c, bid) in &sorted[1..] {
        tail.next = Some(Box::new(CBundlenode {
            start: s, end: e, cov: c, bid,
            next: None, hardstart: false, hardend: false,
        }));
        tail = tail.next.as_mut().unwrap();
    }
    Some(head)
}

/// Get the coordinate range of a component.
pub fn component_range(
    bundlenodes: &[(u64, u64, f64, usize)],
    indices: &[usize],
) -> (u64, u64) {
    let start = indices.iter().map(|&i| bundlenodes[i].0).min().unwrap_or(0);
    let end = indices.iter().map(|&i| bundlenodes[i].1).max().unwrap_or(0);
    (start, end)
}

/// Filter junctions to only those within a component's coordinate range.
pub fn filter_junctions_for_component(
    junctions: &[Junction],
    range_start: u64,
    range_end: u64,
) -> Vec<Junction> {
    junctions
        .iter()
        .filter(|j| j.donor >= range_start && j.acceptor <= range_end)
        .copied()
        .collect()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_single_bundlenode() {
        let bnodes = vec![(100, 200, 10.0, 0)];
        let juncs = vec![];
        let comps = find_junction_connected_components(&bnodes, &juncs);
        assert_eq!(comps.len(), 1);
        assert_eq!(comps[0], vec![0]);
    }

    #[test]
    fn test_junction_connects_two_bundlenodes() {
        let bnodes = vec![(100, 200, 10.0, 0), (500, 600, 5.0, 1)];
        let juncs = vec![Junction::new(200, 500)]; // donor at 200 (in bn0), acceptor at 500 (in bn1)
        let comps = find_junction_connected_components(&bnodes, &juncs);
        assert_eq!(comps.len(), 1); // Connected by junction
    }

    #[test]
    fn test_disconnected_bundlenodes() {
        let bnodes = vec![(100, 200, 10.0, 0), (500, 600, 5.0, 1)];
        let juncs = vec![]; // No junctions connecting them
        let comps = find_junction_connected_components(&bnodes, &juncs);
        assert_eq!(comps.len(), 2); // Not connected
    }

    #[test]
    fn test_contiguous_bundlenodes() {
        let bnodes = vec![(100, 200, 10.0, 0), (201, 300, 5.0, 1)];
        let juncs = vec![];
        let comps = find_junction_connected_components(&bnodes, &juncs);
        assert_eq!(comps.len(), 1); // Contiguous = connected
    }
}
