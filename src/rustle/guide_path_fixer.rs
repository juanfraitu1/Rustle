//! Guide path fixer: ensures all guide transcripts are discoverable in the graph
//!
//! In guided mode, after graph creation, this verifies that each guide transcript
//! can be found as a path through the graph. If a guide can't be found, it
//! force-creates the missing nodes and edges needed to make it traversable.

use crate::graph::Graph;
use crate::reference_gtf::{find_guide_pat, RefTranscript};

/// Verify that all guides can be pathed through the graph.
/// For guides that can't be found, force-create missing nodes/edges.
pub fn verify_and_fix_guide_paths(
    graph: &mut Graph,
    guides: &[RefTranscript],
    bundle_start: u64,
    bundle_end: u64,
    ssdist: u64,
    verbose: bool,
) -> (usize, usize, Vec<(u64, u64)>) {
    let mut found = 0usize;
    let mut fixed = 0usize;
    let mut missing_junctions = Vec::new();

    for guide in guides {
        // Try to find the guide in the current graph
        if find_guide_pat(guide, graph, bundle_start, bundle_end, '.', ssdist).is_some() {
            found += 1;
            continue;
        }

        // Guide not found — reconstruct its path
        if let Some(junctions) = reconstruct_guide_path(graph, guide, ssdist, verbose) {
            missing_junctions.extend(junctions);
            fixed += 1;
            found += 1;
        }
    }

    if verbose && fixed > 0 {
        eprintln!(
            "      [Guide path fixer] Fixed {} guides (found {}/{} total, {} new junctions)",
            fixed,
            found,
            guides.len(),
            missing_junctions.len()
        );
    }

    (found, fixed, missing_junctions)
}

/// Reconstruct a guide transcript's path by force-creating missing nodes/edges.
///
/// Algorithm:
/// 1. For each guide exon, ensure a node exists covering it
/// 2. For each pair of consecutive exons, ensure an edge (junction) connects them
/// 3. Return the junctions created
fn reconstruct_guide_path(
    graph: &mut Graph,
    guide: &RefTranscript,
    ssdist: u64,
    _verbose: bool,
) -> Option<Vec<(u64, u64)>> {
    if guide.exons.is_empty() {
        return None;
    }

    let source_id = graph.source_id;
    let sink_id = graph.sink_id;

    // Create or find nodes for each guide exon
    let mut node_ids: Vec<usize> = Vec::new();
    for (exon_idx, (exon_start, exon_end)) in guide.exons.iter().enumerate() {
        // Find or create a node covering this exon
        let node_id = find_or_create_exon_node(graph, *exon_start, *exon_end, ssdist);
        node_ids.push(node_id);

        // Ensure this node is connected to source if it's the first exon
        if exon_idx == 0 && !graph.nodes[node_id].parents.contains(source_id) {
            graph.add_edge(source_id, node_id);
        }

        // Ensure this node is connected to sink if it's the last exon
        if exon_idx == guide.exons.len() - 1 && !graph.nodes[node_id].children.contains(sink_id) {
            graph.add_edge(node_id, sink_id);
        }
    }

    // Create junctions (edges) between consecutive exons
    let mut new_junctions = Vec::new();
    for i in 0..node_ids.len().saturating_sub(1) {
        let from_node = &graph.nodes[node_ids[i]];
        let to_node = &graph.nodes[node_ids[i + 1]];

        let donor = from_node.end;
        let acceptor = to_node.start;

        // Ensure edge exists
        if !graph.nodes[node_ids[i]].children.contains(node_ids[i + 1]) {
            graph.add_edge(node_ids[i], node_ids[i + 1]);
            new_junctions.push((donor, acceptor));
        }
    }

    Some(new_junctions)
}

/// Find an existing node covering the exon, or create a new one.
fn find_or_create_exon_node(graph: &mut Graph, exon_start: u64, exon_end: u64, ssdist: u64) -> usize {
    // First, try to find an existing node that covers this exon
    for (nid, node) in graph.nodes.iter().enumerate() {
        if nid == graph.source_id || nid == graph.sink_id {
            continue;
        }
        // Node covers exon if it overlaps and boundaries are close
        if node.start <= exon_start + ssdist && node.end >= exon_end - ssdist {
            return nid;
        }
    }

    // No covering node found — create one
    let _node = graph.add_node(exon_start, exon_end);
    // The new node ID is the index of the last node in the vector
    graph.nodes.len() - 1
}
