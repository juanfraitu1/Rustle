//! Deferred synthetic transfrag links (`futuretr`-style lifecycle).
//!
//! stores deferred triples `(n1, n2, abundance)` during graph construction,
//! updates/deletes them as graph evolves, then materializes transfrags later.
//! This module provides the same deferred data flow for Rust pipeline stages.

use crate::graph::{Graph, GraphTransfrag};
use crate::graph_build::PruneRedirect;
use crate::types::DetHashMap as HashMap;

/// Deferred synthetic link:
/// - `from == graph.source_id` means source->to
/// - `to == graph.sink_id` means from->sink
#[derive(Debug, Clone)]
pub struct FutureLink {
    pub from: usize,
    pub to: usize,
    pub abundance: f64,
    pub longread: bool,
    /// Apply futuretr sink suppression heuristic
    /// for node->sink synthetic links.
    pub suppress_sink_heuristic: bool,
}

impl FutureLink {
    pub fn new(
        from: usize,
        to: usize,
        abundance: f64,
        longread: bool,
        suppress_sink_heuristic: bool,
    ) -> Self {
        Self {
            from,
            to,
            abundance,
            longread,
            suppress_sink_heuristic,
        }
    }
}

/// Convert an already-built synthetic transfrag list into deferred future links.
pub fn collect_from_transfrags(
    out: &mut Vec<FutureLink>,
    transfrags: &[GraphTransfrag],
    suppress_sink_heuristic: bool,
    _source_id: usize,
    _sink_id: usize,
) {
    for tf in transfrags {
        if tf.node_ids.len() != 2 {
            continue;
        }
        let from = tf.node_ids[0];
        let to = tf.node_ids[1];
        out.push(FutureLink::new(
            from,
            to,
            tf.abundance,
            tf.longread,
            suppress_sink_heuristic,
        ));
    }
}

/// Remove obviously invalid deferred links and collapse duplicates.
pub fn normalize_links(graph: &Graph, links: &mut Vec<FutureLink>) {
    links.retain(|l| {
        let from_active = l.from == graph.source_id
            || l.from == graph.sink_id
            || graph
                .nodes
                .get(l.from)
                .map(|n| !n.parents.is_empty() || !n.children.is_empty())
                .unwrap_or(false);
        let to_active = l.to == graph.source_id
            || l.to == graph.sink_id
            || graph
                .nodes
                .get(l.to)
                .map(|n| !n.parents.is_empty() || !n.children.is_empty())
                .unwrap_or(false);
        l.from < graph.n_nodes
            && l.to < graph.n_nodes
            && l.from != l.to
            && l.abundance > 0.0
            && l.from != graph.sink_id
            && l.to != graph.source_id
            && from_active
            && to_active
    });
    links.sort_by(|a, b| {
        a.from
            .cmp(&b.from)
            .then(a.to.cmp(&b.to))
            .then(a.longread.cmp(&b.longread))
            .then(a.suppress_sink_heuristic.cmp(&b.suppress_sink_heuristic))
            .then_with(|| {
                b.abundance
                    .partial_cmp(&a.abundance)
                    .unwrap_or(std::cmp::Ordering::Equal)
            })
    });
    let mut i = 0usize;
    while i + 1 < links.len() {
        if links[i].from == links[i + 1].from
            && links[i].to == links[i + 1].to
            && links[i].longread == links[i + 1].longread
            && links[i].suppress_sink_heuristic == links[i + 1].suppress_sink_heuristic
        {
            links[i].abundance += links[i + 1].abundance;
            links.remove(i + 1);
        } else {
            i += 1;
        }
    }
}

/// Apply prune-time node redirects to deferred links (`futuretr` update/delete intent).
pub fn apply_prune_redirects(
    links: &mut [FutureLink],
    redirects: &HashMap<usize, PruneRedirect>,
    source_id: usize,
    sink_id: usize,
) {
    for l in links.iter_mut() {
        if let Some(redir) = redirects.get(&l.from) {
            l.from = redir.upstream.unwrap_or(source_id);
        }
        if let Some(redir) = redirects.get(&l.to) {
            l.to = redir.downstream.unwrap_or(sink_id);
        }
    }
}

/// suppression:
/// if immediate downstream contiguous chain (within `sink_anchor`) already contains
/// a node with direct sink child, skip adding this extra future node->sink link.
fn suppress_sink_future_link(graph: &Graph, from: usize, sink_anchor: u64) -> bool {
    let sink = graph.sink_id;
    if from + 1 >= sink || sink >= graph.n_nodes {
        return false;
    }

    let mut prev_id = from;
    let mut next_id = from + 1;
    let Some(next0) = graph.nodes.get(next_id) else {
        return false;
    };
    let mut dist = next0.length().saturating_sub(1);

    while next_id < sink {
        let Some(prev) = graph.nodes.get(prev_id) else {
            break;
        };
        let Some(next) = graph.nodes.get(next_id) else {
            break;
        };
        if next.start != prev.end || dist >= sink_anchor {
            break;
        }
        if next.children.contains(sink) {
            return true;
        }

        next_id = next_id.saturating_add(1);
        if next_id == sink {
            break;
        }
        if let Some(nn) = graph.nodes.get(next_id) {
            dist = dist.saturating_add(nn.length());
        } else {
            break;
        }
        prev_id = next_id.saturating_sub(1);
    }
    false
}

/// Materialize deferred links into synthetic graph transfrags.
pub fn materialize_links(
    graph: &mut Graph,
    links: &[FutureLink],
    _trthr: f64,
    _sink_anchor: u64,
) -> Vec<GraphTransfrag> {
    let mut out = Vec::new();
    for l in links {
        if l.from >= graph.n_nodes || l.to >= graph.n_nodes {
            continue;
        }
        if l.to == graph.sink_id
            && l.suppress_sink_heuristic
            && suppress_sink_future_link(graph, l.from, _sink_anchor)
        {
            continue;
        }
        graph.add_edge(l.from, l.to);
        let mut tf = GraphTransfrag::new(vec![l.from, l.to], graph.pattern_size());
        graph.set_pattern_edges_for_path(&mut tf.pattern, &tf.node_ids);
        // futuretr abundance is used directly,
        // not clamped to trthr. The abundance comes from coverage-drop calculation.
        tf.abundance = l.abundance.max(0.0);
        tf.read_count = tf.abundance;
        tf.longread = l.longread;
        out.push(tf);
    }
    if !out.is_empty() {
        graph.compute_reachability();
    }
    out
}
