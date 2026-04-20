//! Build splice graph from junctions and bundlenodes (create_graph).

use crate::bitset::NodeSet;
use crate::bpcov::{Bpcov, BpcovStranded, BPCOV_STRAND_ALL, BPCOV_STRAND_MINUS, BPCOV_STRAND_PLUS};
use crate::graph::{Graph, GraphTransfrag, NodeRole};
use crate::longtrim::{
    apply_longtrim_direct, LongtrimBundleSchedule, LongtrimNodeCall, LongtrimStats,
};
use crate::read_boundaries::collect_longtrim_boundary_map;
use crate::trace_events::{
    create_graph_new_node, create_graph_node_final_end, create_graph_node_shrink_jend,
    create_graph_node_shrink_jstart,
};
use crate::types::{
    Bundle2Graph, BundleRead, CBundlenode, DetHashMap as HashMap, DetHashSet as HashSet, Junction,
    JunctionStats, ReadBoundary,
};

/// Constants for coverage-based source/sink edge addition (header).
const ERROR_PERC: f64 = 0.1;
const DROP: f64 = 0.5;

#[inline]
fn trace_strand_index(bundle_strand: char) -> usize {
    if bundle_strand == '+' {
        1
    } else {
        0
    }
}

fn reachable_from_source(graph: &Graph) -> crate::bitset::SmallBitset {
    let mut visited = crate::bitset::SmallBitset::with_capacity(graph.n_nodes.min(64));
    if graph.n_nodes == 0 {
        return visited;
    }
    let src = graph.source_id.min(graph.n_nodes - 1);
    let mut stack = vec![src];
    while let Some(nid) = stack.pop() {
        if visited.contains(nid) {
            continue;
        }
        visited.insert_grow(nid);
        for child in graph.nodes[nid].children.ones().collect::<Vec<_>>() {
            if !visited.contains(child) {
                stack.push(child);
            }
        }
    }
    visited
}

#[derive(Debug, Clone, Copy, Default)]
pub struct PruneRedirect {
    pub upstream: Option<usize>,
    pub downstream: Option<usize>,
}

/// Add source/sink edges where coverage drops sharply (4145-4209).
///
/// For each real node: if total parent coverage << node coverage, add source→node edge.
/// If total child coverage << node coverage, add node→sink edge.
/// Uses strand-specific per-base coverage (`get_cov` equivalent via `BpcovStranded`).
///
/// Returns coverage-proportional transfrags matching futuretr behavior:
/// source edges get abundance `(icov - parcov) / DROP`,
/// sink edges get abundance `(icov - chcov) / DROP`.
pub fn add_coverage_source_sink_edges(
    graph: &mut Graph,
    bpcov: &BpcovStranded,
    strand_idx: usize,
) -> Vec<GraphTransfrag> {
    let sno = match strand_idx {
        BPCOV_STRAND_MINUS | BPCOV_STRAND_PLUS | BPCOV_STRAND_ALL => strand_idx,
        _ => BPCOV_STRAND_ALL,
    };
    let source_id = graph.source_id;
    let sink_id = graph.sink_id;
    let n_nodes = graph.n_nodes;
    let threshold_frac = ERROR_PERC * DROP; // 0.05

    // Precompute per-node average coverage
    let mut node_avg_cov: Vec<f64> = Vec::with_capacity(n_nodes);
    for i in 0..n_nodes {
        let node = &graph.nodes[i];
        if i == source_id || i == sink_id || node.end <= node.start {
            node_avg_cov.push(0.0);
            continue;
        }
        let len = node.end.saturating_sub(node.start);
        if len == 0 {
            node_avg_cov.push(0.0);
            continue;
        }
        let si = bpcov.plus.idx(node.start);
        let ei = bpcov.plus.idx(node.end);
        let cov_sum = bpcov.get_cov_range(sno, si, ei);
        node_avg_cov.push(cov_sum / len as f64);
    }

    // Collect edges and their coverage values (avoid borrow issues)
    let mut add_source_edges: Vec<(usize, f64)> = Vec::new(); // (nid, abundance)
    let mut add_sink_edges: Vec<(usize, f64)> = Vec::new(); // (nid, abundance)

    for i in 0..n_nodes {
        if i == source_id || i == sink_id {
            continue;
        }
        // Skip nodes whose role opts out of coverage aggregation.
        if !graph.nodes[i].role.accrues_coverage() {
            continue;
        }
        // source coverage-drop check starts at i>1 (skips node 1,
        // which is the first real node after source). This prevents spurious source
        // edges to nodes immediately after the source.
        let skip_source_check = i <= 1;
        let icov = node_avg_cov[i];
        if icov <= 0.0 {
            continue;
        }

        // Check if node needs source edge (4147-4172)
        let parents: Vec<usize> = graph.nodes[i].parents.ones().collect();
        if !skip_source_check && !parents.is_empty() && !parents.contains(&source_id) {
            let mut parcov = 0.0;
            for &p in &parents {
                if p == source_id {
                    parcov = f64::MAX;
                    break;
                }
                parcov += node_avg_cov.get(p).copied().unwrap_or(0.0);
            }
            // Recurse through zero-cov parent chain to find real ancestor coverage.
            // Rustle can create phantom single-base zero-cov nodes between exon-end and
            // intergenic nodes that hide the true parent coverage from this check.
            // StringTie doesn't create such nodes, so its parcov directly sees the
            // high-coverage exon node. Recurse if enabled.
            // Default ON: recurse through zero-cov phantom nodes to find real parent/child cov.
            // Disable via RUSTLE_COVLINK_RECURSE_ZERO_OFF=1.
            let recurse_zero = std::env::var_os("RUSTLE_COVLINK_RECURSE_ZERO_OFF").is_none();
            let mut effective_parcov = parcov;
            if recurse_zero && effective_parcov < icov * threshold_frac {
                let mut seen: std::collections::HashSet<usize> = std::collections::HashSet::new();
                let mut stack: Vec<usize> = parents.iter().copied().filter(|&p| node_avg_cov.get(p).copied().unwrap_or(0.0) <= 0.01).collect();
                while let Some(p) = stack.pop() {
                    if !seen.insert(p) { continue; }
                    let grand: Vec<usize> = graph.nodes[p].parents.ones().collect();
                    for &gp in &grand {
                        if gp == source_id { continue; }
                        let c = node_avg_cov.get(gp).copied().unwrap_or(0.0);
                        if c > 0.01 {
                            effective_parcov += c;
                        } else {
                            stack.push(gp);
                        }
                    }
                }
            }
            if effective_parcov < icov * threshold_frac {
                // ref:4168: abundance = (icov - parcov) / DROP
                let abundance = (icov - parcov) / DROP;
                add_source_edges.push((i, abundance));
            }
        }

        // Check if node needs sink edge (4175-4206)
        let children: Vec<usize> = graph.nodes[i].children.ones().collect();
        if !children.is_empty() && !children.contains(&sink_id) {
            let mut chcov = 0.0;
            for &c in &children {
                if c == sink_id {
                    chcov = f64::MAX;
                    break;
                }
                chcov += node_avg_cov.get(c).copied().unwrap_or(0.0);
            }
            // Default ON: recurse through zero-cov phantom nodes to find real parent/child cov.
            // Disable via RUSTLE_COVLINK_RECURSE_ZERO_OFF=1.
            let recurse_zero = std::env::var_os("RUSTLE_COVLINK_RECURSE_ZERO_OFF").is_none();
            let mut effective_chcov = chcov;
            if recurse_zero && effective_chcov < icov * threshold_frac {
                let mut seen: std::collections::HashSet<usize> = std::collections::HashSet::new();
                let mut stack: Vec<usize> = children.iter().copied().filter(|&c| node_avg_cov.get(c).copied().unwrap_or(0.0) <= 0.01).collect();
                while let Some(c) = stack.pop() {
                    if !seen.insert(c) { continue; }
                    let grand: Vec<usize> = graph.nodes[c].children.ones().collect();
                    for &gc in &grand {
                        if gc == sink_id { continue; }
                        let cc = node_avg_cov.get(gc).copied().unwrap_or(0.0);
                        if cc > 0.01 {
                            effective_chcov += cc;
                        } else {
                            stack.push(gc);
                        }
                    }
                }
            }
            if effective_chcov < icov * threshold_frac {
                // ref:4200: abundance = (icov - chcov) / DROP
                let abundance = (icov - chcov) / DROP;
                add_sink_edges.push((i, abundance));
            }
        }
    }

    let pattern_size = graph.pattern_size();
    let mut synth_transfrags: Vec<GraphTransfrag> = Vec::new();

    let trace = std::env::var_os("RUSTLE_TRACE_COVLINKS").is_some();
    for (nid, abundance) in add_source_edges {
        if trace {
            let n = &graph.nodes[nid];
            eprintln!("COVLINK_SRC nid={} range={}..{} abund={:.2}", nid, n.start, n.end, abundance);
        }
        graph.add_edge(source_id, nid);
        let mut tf = GraphTransfrag::new(vec![source_id, nid], pattern_size);
        tf.abundance = abundance;
        synth_transfrags.push(tf);
    }
    for (nid, abundance) in add_sink_edges {
        if trace {
            let n = &graph.nodes[nid];
            eprintln!("COVLINK_SNK nid={} range={}..{} abund={:.2}", nid, n.start, n.end, abundance);
        }
        graph.add_edge(nid, sink_id);
        let mut tf = GraphTransfrag::new(vec![nid, sink_id], pattern_size);
        tf.abundance = abundance;
        synth_transfrags.push(tf);
    }

    synth_transfrags
}

fn collect_bundlenodes(bn: Option<&CBundlenode>) -> Vec<(usize, u64, u64, f64)> {
    let mut out = Vec::new();
    let mut cur = bn;
    while let Some(n) = cur {
        if n.end > n.start {
            out.push((n.bid, n.start, n.end, n.cov));
        }
        cur = n.next.as_deref();
    }
    out
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord)]
enum JunctionEventType {
    End = 0,
    Start = 1,
}

#[inline]
fn target_bundle_strand(bundle_strand: char) -> Option<i8> {
    match bundle_strand {
        '-' => Some(-1i8),
        '+' => Some(1i8),
        _ => None,
    }
}

fn filter_junctions_for_bundle<'a>(
    junctions: &'a [Junction],
    bundle_strand: char,
    junction_stats: Option<&JunctionStats>,
) -> Vec<&'a Junction> {
    let target_strand = target_bundle_strand(bundle_strand);
    junctions
        .iter()
        .filter(|j| {
            if let Some(js) = junction_stats {
                if let Some(stat) = js.get(j) {
                    if let Some(ts) = target_strand {
                        if stat.strand != Some(ts) {
                            return false;
                        }
                    }
                    if stat.mm < 0.0 {
                        return false;
                    }
                    if stat.strand == Some(0) {
                        return false;
                    }
                }
            }
            true
        })
        .collect()
}

fn collect_bundlenode_events(
    filtered_juncs: &[&Junction],
    currentstart: u64,
    endbundle: u64,
    bundle_end: u64,
) -> Vec<(u64, JunctionEventType)> {
    let mut events: Vec<(u64, JunctionEventType)> = Vec::new();
    for j in filtered_juncs {
        if j.donor > currentstart && j.donor <= endbundle {
            if j.acceptor <= bundle_end {
                events.push((j.donor, JunctionEventType::Start));
            }
        }
        if j.acceptor >= currentstart && j.acceptor <= endbundle {
            events.push((j.acceptor, JunctionEventType::End));
        }
    }
    events.sort_unstable();
    events.dedup_by_key(|e| (e.0, e.1));
    events
}

fn build_longtrim_bundle_schedules(
    graph: &Graph,
    junctions: &[Junction],
    bundlenodes: Option<&CBundlenode>,
    bundle_end: u64,
    bundle_strand: char,
    junction_stats: Option<&JunctionStats>,
) -> Vec<LongtrimBundleSchedule> {
    let filtered_juncs = filter_junctions_for_bundle(junctions, bundle_strand, junction_stats);
    let mut nodes_by_bid: HashMap<usize, Vec<usize>> = Default::default();
    for (nid, node) in graph.nodes.iter().enumerate() {
        if nid == graph.source_id || nid == graph.sink_id {
            continue;
        }
        if !node.role.longtrim_schedule() {
            continue;
        }
        if let Some(bid) = node.source_bnode {
            nodes_by_bid.entry(bid).or_default().push(nid);
        }
    }
    for node_ids in nodes_by_bid.values_mut() {
        node_ids.sort_unstable_by_key(|nid| graph.nodes[*nid].start);
    }

    let mut schedules = Vec::new();
    let mut bnode_opt = bundlenodes;
    while let Some(bn) = bnode_opt {
        let node_ids = nodes_by_bid.get(&bn.bid).cloned().unwrap_or_default();
        if node_ids.is_empty() {
            bnode_opt = bn.next.as_deref();
            continue;
        }

        let events = collect_bundlenode_events(&filtered_juncs, bn.start, bn.end, bundle_end);
        let mut calls: Vec<LongtrimNodeCall> = Vec::new();
        let mut node_idx = 0usize;
        let mut current_start = bn.start;
        let mut ei = 0usize;

        while ei < events.len() && node_idx < node_ids.len() {
            let (pos, ev_type) = events[ei];
            match ev_type {
                JunctionEventType::Start => {
                    calls.push(LongtrimNodeCall {
                        source_bid: bn.bid,
                        node_id: node_ids[node_idx],
                        endcov: true,
                        post_startcov: Some(true),
                    });
                    node_idx += 1;
                    while ei < events.len()
                        && events[ei].0 == pos
                        && events[ei].1 == JunctionEventType::Start
                    {
                        ei += 1;
                    }
                    if pos < bn.end {
                        current_start = pos;
                    } else {
                        break;
                    }
                }
                JunctionEventType::End => {
                    while ei < events.len()
                        && events[ei].0 == pos
                        && events[ei].1 == JunctionEventType::End
                    {
                        ei += 1;
                    }
                    if current_start < pos {
                        calls.push(LongtrimNodeCall {
                            source_bid: bn.bid,
                            node_id: node_ids[node_idx],
                            endcov: false,
                            post_startcov: Some(false),
                        });
                        node_idx += 1;
                        current_start = pos;
                    }
                }
            }
        }

        for (i, &node_id) in node_ids.iter().enumerate().skip(node_idx) {
            calls.push(LongtrimNodeCall {
                source_bid: bn.bid,
                node_id,
                endcov: true,
                post_startcov: if i + 1 == node_ids.len() {
                    None
                } else {
                    Some(false)
                },
            });
        }

        schedules.push(LongtrimBundleSchedule {
            source_bid: bn.bid,
            calls,
        });
        bnode_opt = bn.next.as_deref();
    }

    schedules
}

/// Build bundle2graph mapping (CGraphinfo / bundle2graph).
/// For each bundlenode (bid), record graph node ids that were created from it.
pub fn build_bundle2graph(graph: &Graph, bundlenodes: Option<&CBundlenode>) -> Bundle2Graph {
    let mut bids: Vec<usize> = Vec::new();
    let mut cur = bundlenodes;
    while let Some(n) = cur {
        bids.push(n.bid);
        cur = n.next.as_deref();
    }
    if bids.is_empty() {
        return Vec::new();
    }
    let max_bid = bids.into_iter().max().unwrap_or(0);
    let mut map: Bundle2Graph = vec![Vec::new(); max_bid + 1];
    for (nid, node) in graph.nodes.iter().enumerate() {
        if nid == graph.source_id || nid == graph.sink_id {
            continue;
        }
        // Skip nodes whose role opts out of read routing.
        if !node.role.accepts_reads() {
            continue;
        }
        if let Some(bid) = node.source_bnode {
            if bid < map.len() {
                map[bid].push((0, nid));
            }
        }
    }
    for bucket in &mut map {
        if !bucket.is_empty() {
            bucket.sort_unstable_by_key(|&(_, nid)| graph.nodes[nid].start);
        }
    }
    map
}

/// Build graph from junctions and bundlenodes. Source=0, sink=last.
/// CRITICAL: Filters junctions by strand, matching behavior.
/// uses formula `(junction.strand+1) == 2*s` where s is bundle strand (0,1,2).
/// Only junctions matching the bundle strand are used for graph edges.
///
/// Uses Original-style dynamic junction processing
/// processes junction start/end events in coordinate order, splitting the
/// current graphnode at each event. This matches node structure including
/// the short-tail skip optimization.
pub fn create_graph(
    junctions: &[Junction],
    _bundle_start: u64,
    bundle_end: u64,
    bundlenodes: Option<&CBundlenode>,
    junction_support: u64,
    _reads: Option<&[BundleRead]>,
    bundle_strand: char,
    junction_stats: Option<&JunctionStats>,
    bpcov: Option<&Bpcov>,
) -> Graph {
    create_graph_inner(junctions, _bundle_start, bundle_end, bundlenodes,
        junction_support, _reads, bundle_strand, junction_stats, bpcov, None, &[], &[])
}

fn create_graph_inner(
    junctions: &[Junction],
    _bundle_start: u64,
    bundle_end: u64,
    bundlenodes: Option<&CBundlenode>,
    junction_support: u64,
    _reads: Option<&[BundleRead]>,
    bundle_strand: char,
    junction_stats: Option<&JunctionStats>,
    bpcov: Option<&Bpcov>,
    bpcov_stranded: Option<&BpcovStranded>,
    lstart: &[ReadBoundary],
    lend: &[ReadBoundary],
) -> Graph {
    let mut graph = Graph::new();
    let trace_s = trace_strand_index(bundle_strand);
    let trace_g = 0usize;

    graph.add_node(0, 0);
    graph.source_id = 0;

    let good_junctions: HashSet<Junction> = junctions.iter().copied().collect();
    if good_junctions.is_empty() {
        let mut segs = collect_bundlenodes(bundlenodes);
        segs.sort_unstable_by_key(|(_, s, _, _)| *s);
        let mut node_ids = Vec::new();
        for (bid, s, e, _cov) in segs {
            let node = graph.add_node(s, e);
            node.source_bnode = Some(bid);
            // `CGraphnode::cov` starts at zero and is accumulated later from mapped reads.
            // Do not seed it from bundlenode coverage.
            node.coverage = 0.0;
            create_graph_new_node(trace_s, trace_g, node.node_id, s, e, "initial bundlenode");
            create_graph_node_final_end(trace_s, trace_g, node.node_id, s, e);
            node_ids.push(node.node_id);
        }
        graph.add_node(0, 0);
        let sink_id = graph.n_nodes - 1;
        for (i, &nid) in node_ids.iter().enumerate() {
            graph.add_edge(graph.source_id, nid);
            graph.add_edge(nid, sink_id);
            if i + 1 < node_ids.len() {
                let a = &graph.nodes[nid];
                let b = &graph.nodes[node_ids[i + 1]];
                if a.end >= b.start {
                    graph.add_edge(nid, node_ids[i + 1]);
                }
            }
        }
        graph.sink_id = sink_id;
        graph.compute_reachability();
        return graph;
    }

    // Filter junctions by strand and killed status (once, reuse for all bundlenodes).
    let filtered_juncs = filter_junctions_for_bundle(junctions, bundle_strand, junction_stats);

    // ends[] hashmap: tracks junction endpoints for cross-junction linking.
    // Persists across bundlenodes ( 4080-4088).
    let mut ends: HashMap<u64, Vec<usize>> = Default::default();

    let mut sink_parents: Vec<usize> = Vec::new();
    let mut bnode_opt = bundlenodes;
    while let Some(bn) = bnode_opt {
        let currentstart = bn.start;
        let endbundle = bn.end;
        let source_bid = bn.bid;

        // Check if any previously-recorded junction ends at currentstart
        let has_end_at_start = ends.contains_key(&currentstart);

        // Collect junction events within this bundlenode.
        // START event at donor: junction starts here (exon ends, intron begins).
        // END event at acceptor: junction ends here (intron ends, exon begins).
        let events =
            collect_bundlenode_events(&filtered_juncs, currentstart, endbundle, bundle_end);

        // Create initial graphnode [currentstart, endbundle)
        let node = graph.add_node(currentstart, endbundle);
        node.source_bnode = Some(source_bid);
        node.coverage = 0.0;
        let mut graphnode_id = node.node_id;
        create_graph_new_node(
            trace_s,
            trace_g,
            graphnode_id,
            currentstart,
            endbundle,
            "initial bundlenode",
        );

        // Link initial node: either from ends[] or from source
        let mut linked_initial = false;
        if has_end_at_start {
            if let Some(parent_ids) = ends.get(&currentstart) {
                for &pid in parent_ids {
                    graph.add_edge(pid, graphnode_id);
                    linked_initial = true;
                }
            }
        }
        if !linked_initial {
            graph.add_edge(graph.source_id, graphnode_id);
        }

        let mut completed = false;
        // Collect overlap-alias node IDs created within this bundlenode,
        // so we can patch their children (= final split node's children)
        // after the event loop finishes.
        let mut bundlenode_aliases: Vec<usize> = Vec::new();
        // (-3880): compute longtrim boundaries PER
        // BUNDLENODE using the coverage derivative sliding window. the original implementation
        // runs this detection inside the create_graph per-bundlenode loop,
        // not externally per-bundle. When external lstart/lend are empty but
        // bpcov is available, generate per-bundlenode boundaries here.
        let bnode_lstart_owned;
        let bnode_lend_owned;
        // compute boundaries per-bundlenode when longtrim is active
        // and no external boundaries were provided. External boundaries come from
        // CPAS/poly-A evidence; when absent, use coverage-derivative detection.
        let bnode_span = endbundle.saturating_sub(currentstart);
        let min_lt_span = 2 * (100 + 50) + 50; // 2*(CHI_WIN+CHI_THR) + margin = 350
        let (effective_lstart, effective_lend): (&[ReadBoundary], &[ReadBoundary]) =
            if lstart.is_empty() && lend.is_empty()
                && std::env::var_os("RUSTLE_DISABLE_LONGTRIM").is_none()
                && bnode_span >= min_lt_span
            {
                if let Some(bpc) = bpcov {
                    // Build strand-specific bpcov for this bundlenode
                    let strand_bpc = if let Some(bps) = bpcov_stranded {
                        let n = bpc.cov.len();
                        let mut cov = vec![0.0f64; n];
                        let opp = match bundle_strand {
                            '-' => BPCOV_STRAND_PLUS,
                            '+' => BPCOV_STRAND_MINUS,
                            _ => BPCOV_STRAND_ALL,
                        };
                        for k in 0..n {
                            let total = bps.get_cov_range(BPCOV_STRAND_ALL, k, k + 1);
                            let o = if opp != BPCOV_STRAND_ALL {
                                bps.get_cov_range(opp, k, k + 1)
                            } else { 0.0 };
                            cov[k] = (total - o).max(0.0);
                        }
                        crate::bpcov::Bpcov::from_cov(cov, bpc.bundle_start, bpc.bundle_end)
                    } else {
                        bpc.clone()
                    };
                    let (ls, le) = crate::read_boundaries::collect_longtrim_boundaries_in_span(
                        &strand_bpc, currentstart, endbundle, &[], &[],
                    );
                    bnode_lstart_owned = ls;
                    bnode_lend_owned = le;
                    (bnode_lstart_owned.as_slice(), bnode_lend_owned.as_slice())
                } else {
                    bnode_lstart_owned = Vec::new();
                    bnode_lend_owned = Vec::new();
                    (lstart, lend)
                }
            } else {
                bnode_lstart_owned = Vec::new();
                bnode_lend_owned = Vec::new();
                (lstart, lend)
            };
        // Longtrim state: pointers into lstart/lend arrays (never backtrack).
        let mut nls = 0usize;
        let mut nle = 0usize;
        let has_longtrim = !effective_lstart.is_empty() || !effective_lend.is_empty();
        // Skip events before this bundlenode.
        while nls < effective_lstart.len() && effective_lstart[nls].pos < currentstart { nls += 1; }
        while nle < effective_lend.len() && effective_lend[nle].pos < currentstart { nle += 1; }

        // Process events in coordinate order (do-while loop).
        let mut ei = 0;
        while ei < events.len() && !completed {
            let (pos, ev_type) = events[ei];

            // call longtrim BEFORE each junction event to process
            // any lstart/lend events between the current node and this junction.
            if has_longtrim {
                if let Some(bpc) = bpcov {
                    let nodeend = pos; // Process up to this junction position.
                    let has_junction_start = graph.nodes[graphnode_id].parents.count_ones() > 0
                        && !graph.nodes[graphnode_id].parents.contains(graph.source_id);
                    let has_junction_end = ev_type == JunctionEventType::Start; // junction starts here = exon boundary
                    longtrim_inline(
                        &mut graph, &mut graphnode_id, &mut nls, &mut nle,
                        effective_lstart, effective_lend, bpc, bpcov_stranded, bundle_strand,
                        _bundle_start, nodeend,
                        has_junction_start, has_junction_end, source_bid,
                        &mut sink_parents,
                    );
                }
            }

            match ev_type {
                JunctionEventType::Start => {
                    // Junction start at donor position
                    // Set current graphnode end to donor.
                    let old_end = graph.nodes[graphnode_id].end;
                    graph.nodes[graphnode_id].end = pos;
                    create_graph_node_shrink_jstart(
                        trace_s,
                        trace_g,
                        graphnode_id,
                        graph.nodes[graphnode_id].start,
                        old_end,
                        pos,
                        pos,
                    );

                    // Record all junctions starting at this position in ends[]
                    for j in &filtered_juncs {
                        if j.donor == pos {
                            ends.entry(j.acceptor).or_default().push(graphnode_id);
                        }
                    }

                    // Advance past all START events at this position.
                    while ei < events.len()
                        && events[ei].0 == pos
                        && events[ei].1 == JunctionEventType::Start
                    {
                        ei += 1;
                    }

                    if pos < endbundle {
                        // Short-tail skip
                        // If remaining bundlenode is shorter than junction_support and
                        // no more junctions exist in this bundlenode, check coverage.
                        if endbundle - pos < junction_support {
                            let has_more_events = ei < events.len();
                            if !has_more_events {
                                if let Some(bpc) = bpcov {
                                    // Compare coverage: left window vs right window.
                                    // Left: [2*pos - endbundle, pos), Right: [pos, endbundle)
                                    let tail_len = endbundle - pos;
                                    let left_start = if pos >= endbundle + currentstart {
                                        pos - tail_len
                                    } else {
                                        currentstart // clamp to bundlenode start
                                    };
                                    let si_left = bpc.idx(left_start);
                                    let ei_left = bpc.idx(pos);
                                    let si_right = bpc.idx(pos);
                                    let ei_right = bpc.idx(endbundle);
                                    let len_left = ei_left.saturating_sub(si_left).max(1);
                                    let len_right = ei_right.saturating_sub(si_right).max(1);
                                    let covleft =
                                        bpc.get_cov_range(si_left, ei_left) / len_left as f64;
                                    let covright =
                                        bpc.get_cov_range(si_right, ei_right) / len_right as f64;
                                    if covright < covleft * (1.0 - ERROR_PERC) {
                                        completed = true;
                                    }
                                }
                            }
                        }

                        if !completed {
                            // Create next node [donor, endbundle)
                            let next = graph.add_node(pos, endbundle);
                            next.source_bnode = Some(source_bid);
                            next.coverage = 0.0;
                            let next_id = next.node_id;
                            create_graph_new_node(
                                trace_s,
                                trace_g,
                                next_id,
                                pos,
                                endbundle,
                                "nextnode after junction",
                            );
                            graph.add_edge(graphnode_id, next_id);
                            graphnode_id = next_id;
                        }
                    } else {
                        completed = true;
                    }
                }
                JunctionEventType::End => {
                    // Junction end at acceptor position
                    let acceptor = pos;

                    // Advance past all END events at this position
                    while ei < events.len()
                        && events[ei].0 == acceptor
                        && events[ei].1 == JunctionEventType::End
                    {
                        ei += 1;
                    }

                    // Split current node if it started before this acceptor
                    if graph.nodes[graphnode_id].start < acceptor {
                        // Capture the bundlenode-level pre-shrink state for
                        // overlap-alias scaffold (StringTie keeps a full-span
                        // "anchor" node 33 alongside the narrower junction-
                        // acceptor node 34 — see trace PASS_ID=260 STRG.294).
                        let pre_shrink_start = graph.nodes[graphnode_id].start;
                        let old_end = graph.nodes[graphnode_id].end;
                        let shrunk_id = graphnode_id;
                        graph.nodes[graphnode_id].end = acceptor;
                        create_graph_node_shrink_jend(
                            trace_s,
                            trace_g,
                            graphnode_id,
                            graph.nodes[graphnode_id].start,
                            old_end,
                            acceptor,
                            acceptor,
                        );
                        // Create next node [acceptor, endbundle)
                        let next = graph.add_node(acceptor, endbundle);
                        next.source_bnode = Some(source_bid);
                        next.coverage = 0.0;
                        let next_id = next.node_id;
                        create_graph_new_node(
                            trace_s,
                            trace_g,
                            next_id,
                            acceptor,
                            endbundle,
                            "nextnode after junction",
                        );
                        graph.add_edge(graphnode_id, next_id);
                        graphnode_id = next_id;

                        // Overlap-alias scaffold. Three modes:
                        //
                        // RUSTLE_OVERLAP_NODE_MEASURE=1 — log-only counter
                        //   of where an alias WOULD be created.
                        //
                        // RUSTLE_OVERLAP_NODE_SCAFFOLD=1 — create inert alias
                        //   nodes (no edges). CURRENTLY REGRESSES F1.
                        //
                        // RUSTLE_OVERLAP_NODE_EDGES=1 — create alias nodes AND
                        //   wire parents (from shrunk node, i.e., original
                        //   pre-shrink parents) + children (patched at end of
                        //   bundlenode to match final split node's children).
                        //   No transfrags mapped through alias yet, so flow can't
                        //   use it (no capacity). Purely structural parallel edge.
                        let emit_measure = std::env::var_os("RUSTLE_OVERLAP_NODE_MEASURE").is_some();
                        let emit_alias = std::env::var_os("RUSTLE_OVERLAP_NODE_SCAFFOLD").is_some()
                            || std::env::var_os("RUSTLE_OVERLAP_NODE_EDGES").is_some();
                        let wire_edges = std::env::var_os("RUSTLE_OVERLAP_NODE_EDGES").is_some();
                        if (emit_measure || emit_alias)
                            && pre_shrink_start < acceptor
                            && acceptor < endbundle
                        {
                            if emit_measure {
                                eprintln!(
                                    "OVERLAP_ALIAS_CANDIDATE bid={} span={}-{} acceptor={} endbundle={}",
                                    source_bid, pre_shrink_start, endbundle, acceptor, endbundle
                                );
                            }
                            if emit_alias {
                                let alias = graph.add_node(pre_shrink_start, endbundle);
                                alias.source_bnode = Some(source_bid);
                                alias.coverage = 0.0;
                                alias.role = NodeRole::OverlapAnchor;
                                let alias_id = alias.node_id;
                                create_graph_new_node(
                                    trace_s,
                                    trace_g,
                                    alias_id,
                                    pre_shrink_start,
                                    endbundle,
                                    "overlap_alias(jend)",
                                );
                                if wire_edges {
                                    // Inherit parents from the shrunk node
                                    // (shrink doesn't change parents). These are
                                    // the original pre-shrink parents.
                                    let shrunk_parents: Vec<usize> =
                                        graph.nodes[shrunk_id].parents.ones().collect();
                                    for p in shrunk_parents {
                                        if p != alias_id {
                                            graph.add_edge(p, alias_id);
                                        }
                                    }
                                    bundlenode_aliases.push(alias_id);
                                }
                            }
                        }
                    }

                    // Link nodes from ends[] to current graphnode
                    if let Some(parent_ids) = ends.get(&acceptor).cloned() {
                        for pid in parent_ids {
                            graph.add_edge(pid, graphnode_id);
                        }
                    }
                }
            }
        }

        // call longtrim for remaining events up to endbundle.
        if has_longtrim && !completed {
            if let Some(bpc) = bpcov {
                longtrim_inline(
                    &mut graph, &mut graphnode_id, &mut nls, &mut nle,
                    lstart, lend, bpc, bpcov_stranded, bundle_strand,
                    _bundle_start, endbundle,
                    true, true, source_bid,
                    &mut sink_parents,
                );
            }
        }

        if !completed {
            // Set final graphnode end to endbundle
            graph.nodes[graphnode_id].end = endbundle;
            create_graph_node_final_end(
                trace_s,
                trace_g,
                graphnode_id,
                graph.nodes[graphnode_id].start,
                endbundle,
            );
            sink_parents.push(graphnode_id);
        }

        // Patch overlap-alias children to match the final split node
        // (graphnode_id at this point) of this bundlenode. An alias
        // spans [pre_shrink_start, endbundle) and shares the same
        // outgoing junctions/continuations as the last split node.
        // Executed only when RUSTLE_OVERLAP_NODE_EDGES=1 (aliases were
        // already emitted with parents).
        if !bundlenode_aliases.is_empty() {
            let final_children: Vec<usize> =
                graph.nodes[graphnode_id].children.ones().collect();
            for &alias_id in &bundlenode_aliases {
                for &c in &final_children {
                    if c != alias_id {
                        graph.add_edge(alias_id, c);
                    }
                }
            }
        }

        bnode_opt = bn.next.as_deref();
    }

    // Add sink node.
    graph.add_node(0, 0);
    graph.sink_id = graph.n_nodes - 1;
    let sink_id = graph.sink_id;

    // sink links are attached to the tail graphnode of each non-completed
    // bundlenode segment during construction.
    for nid in sink_parents {
        graph.add_edge(nid, sink_id);
    }

    graph.compute_reachability();
    graph
}

#[derive(Debug, Default, Clone, Copy)]
pub struct CreateGraphLongtrimStats {
    pub lstart_events: usize,
    pub lend_events: usize,
    pub applied: bool,
    pub longtrim: LongtrimStats,
}

/// Validate longtrim boundary events using bpcov contrast logic.
///
/// the original implementation only splits at a read start/end position when the coverage contrast
/// across a CHI_THR (50bp) window is positive:
/// - For starts: coverage to the RIGHT > coverage to the LEFT (reads begin here)
/// - For ends: coverage to the LEFT > coverage to the RIGHT (reads end here)
///
/// This pre-filters the raw lstart/lend arrays (which can have 40K+ events)
/// down to ~200-400 validated split points.
#[allow(dead_code)]
#[allow(dead_code)]
fn validate_longtrim_boundaries(
    lstart: &[crate::types::ReadBoundary],
    lend: &[crate::types::ReadBoundary],
    bpcov: &Bpcov,
    bundle_start: u64,
    bundle_strand: char,
) -> (Vec<crate::types::ReadBoundary>, Vec<crate::types::ReadBoundary>) {
    use crate::types::ReadBoundary;
    const CHI_THR: u64 = 50;
    const DROP: f64 = 0.5;
    const LONGINTRONANCHOR: u64 = 25;

    let _strand_idx = match bundle_strand {
        '+' => crate::bpcov::BPCOV_STRAND_PLUS,
        '-' => crate::bpcov::BPCOV_STRAND_MINUS,
        _ => crate::bpcov::BPCOV_STRAND_ALL,
    };

    let get_cov = |start: u64, end: u64| -> f64 {
        if end < start || start < bundle_start {
            return 0.0;
        }
        let s = (start - bundle_start) as usize;
        let e = (end - bundle_start) as usize;
        bpcov.get_cov_range(s, e)
    };

    let mut valid_starts: Vec<ReadBoundary> = Vec::new();
    for b in lstart {
        let pos = b.pos;
        // Skip if too close to node boundaries (longintronanchor check)
        if pos < bundle_start + LONGINTRONANCHOR {
            continue;
        }
        // Contrast: coverage RIGHT of position vs LEFT
        let right_start = pos;
        let right_end = pos + CHI_THR - 1;
        let left_start = pos.saturating_sub(CHI_THR);
        let left_end = pos - 1;
        let right_cov = get_cov(right_start, right_end);
        let left_cov = get_cov(left_start, left_end);
        let tmpcov = (right_cov - left_cov) / (DROP * CHI_THR as f64);
        if tmpcov > 0.0 {
            valid_starts.push(*b);
        }
    }

    let mut valid_ends: Vec<ReadBoundary> = Vec::new();
    for b in lend {
        let pos = b.pos;
        if pos < bundle_start + LONGINTRONANCHOR {
            continue;
        }
        // Contrast: coverage LEFT of position vs RIGHT
        let left_start = pos.saturating_sub(CHI_THR - 1);
        let left_end = pos;
        let right_start = pos + 1;
        let right_end = pos + CHI_THR;
        let left_cov = get_cov(left_start, left_end);
        let right_cov = get_cov(right_start, right_end);
        let tmpcov = (left_cov - right_cov) / (DROP * CHI_THR as f64);
        if tmpcov > 0.0 {
            valid_ends.push(*b);
        }
    }

    (valid_starts, valid_ends)
}

/// Iterative node-split pass matching longtrim().
///
/// For each existing graph node, processes sorted lstart/lend events that fall
/// within the node's range. At each event position, computes bpcov contrast in
/// a CHI_THR=50bp window. If the contrast is positive, the node is split.
///
/// Start splits: creates two nodes [node.start, pos-1] and [pos, node.end].
///   The new node gets a source edge (hardstart).
/// End splits: creates two nodes [node.start, pos] and [pos+1, node.end].
///   The left node gets a sink edge (hardend).
#[allow(dead_code)]
#[allow(dead_code)]
fn apply_iterative_longtrim_splits(
    graph: &mut Graph,
    lstart: &[ReadBoundary],
    lend: &[ReadBoundary],
    bpcov: &Bpcov,
    bundle_start: u64,
    _bundle_strand: char,
) -> usize {
    const CHI_THR: i64 = 50;
    const DROP: f64 = 0.5;
    const LONGINTRONANCHOR: u64 = 25;
    const ERROR_PERC: f64 = 0.1;

    if lstart.is_empty() && lend.is_empty() {
        return 0;
    }

    let source_id = graph.source_id;
    let sink_id = graph.sink_id;
    let bpcov_len = bpcov.cov.len();

    let get_cov_range = |s: i64, e: i64| -> f64 {
        if s < 0 || e < s || s as usize >= bpcov_len {
            return 0.0;
        }
        let si = s as usize;
        let ei = (e as usize).min(bpcov_len.saturating_sub(1));
        bpcov.get_cov_range(si, ei)
    };

    // Collect non-source/sink nodes sorted by start position.
    let mut node_ids: Vec<usize> = (0..graph.nodes.len())
        .filter(|&i| i != source_id && i != sink_id && graph.nodes[i].end > graph.nodes[i].start)
        .collect();
    node_ids.sort_unstable_by_key(|&i| graph.nodes[i].start);

    let mut nls: usize = 0; // Pointer into lstart (never backtracks)
    let mut nle: usize = 0; // Pointer into lend

    let mut splits = 0usize;

    for &orig_nid in &node_ids {
        let node_start = graph.nodes[orig_nid].start;
        let node_end = graph.nodes[orig_nid].end;

        // Is this node bounded by junctions?
        let has_junction_at_start = graph.nodes[orig_nid].parents.count_ones() > 1
            || (graph.nodes[orig_nid].parents.count_ones() == 1
                && !graph.nodes[orig_nid].parents.contains(source_id));
        let has_junction_at_end = graph.nodes[orig_nid].children.count_ones() > 1
            || (graph.nodes[orig_nid].children.count_ones() == 1
                && !graph.nodes[orig_nid].children.contains(sink_id));

        // Skip lstart/lend events before this node.
        while nls < lstart.len() && lstart[nls].pos < node_start {
            nls += 1;
        }
        while nle < lend.len() && lend[nle].pos < node_start {
            nle += 1;
        }

        let mut cur_nid = orig_nid;
        let mut ls = nls;
        let mut le = nle;

        // Process events within this node.
        while (ls < lstart.len() && lstart[ls].pos < node_end)
            || (le < lend.len() && lend[le].pos < node_end)
        {
            let use_start = if le >= lend.len() {
                true
            } else if ls >= lstart.len() {
                false
            } else {
                lstart[ls].pos <= lend[le].pos
            };

            if use_start {
                let pos = lstart[ls].pos;
                let cur_start = graph.nodes[cur_nid].start;
                let cur_end = graph.nodes[cur_nid].end;

                // Proximity check: not too close to junction boundaries
                let start_ok = has_junction_at_start || pos > cur_start + LONGINTRONANCHOR;
                let end_ok = has_junction_at_end || pos < node_end + LONGINTRONANCHOR;

                if start_ok && end_ok && pos > cur_start && pos < cur_end {
                    // Bpcov contrast: coverage RIGHT vs LEFT of pos
                    let startpos = (pos - bundle_start) as i64;
                    let winstart = (startpos - CHI_THR).max(0);
                    let winend = (startpos + CHI_THR - 1).min(bpcov_len as i64 - 1);
                    let right_cov = get_cov_range(startpos, winend);
                    let left_cov = get_cov_range(winstart, startpos - 1);
                    let mut tmpcov = (right_cov - left_cov) / (DROP * CHI_THR as f64);

                    // Negative longcov: force a small positive cov for re-estimation
                    if tmpcov <= 0.0 && lstart[ls].cov < 0.0 {
                        tmpcov = ERROR_PERC;
                    }

                    if tmpcov > 0.0 {
                        // SPLIT: [cur_start, pos-1] and [pos, cur_end]
                        graph.nodes[cur_nid].end = pos; // Rustle uses half-open [start, end)

                        let new_nid = graph.nodes.len();
                        let mut new_node = crate::graph::GraphNode::new(0, 0, 0);
                        new_node.node_id = new_nid;
                        new_node.start = pos;
                        new_node.end = cur_end;
                        new_node.hardstart = true;
                        graph.nodes.push(new_node);
                        graph.n_nodes = graph.nodes.len();

                        // Edge: source → new node (hardstart)
                        graph.nodes[source_id].children.insert_grow(new_nid);
                        graph.nodes[new_nid].parents.insert_grow(source_id);
                        // Edge: prev → new (contiguous)
                        graph.nodes[cur_nid].children.insert_grow(new_nid);
                        graph.nodes[new_nid].parents.insert_grow(cur_nid);

                        cur_nid = new_nid;
                        splits += 1;
                    }
                }
                ls += 1;
            } else {
                let pos = lend[le].pos;
                let cur_start = graph.nodes[cur_nid].start;
                let cur_end = graph.nodes[cur_nid].end;

                let start_ok = !has_junction_at_start || pos > cur_start + LONGINTRONANCHOR;
                let end_ok = !has_junction_at_end || pos < node_end + LONGINTRONANCHOR;

                if start_ok && end_ok && pos > cur_start && pos < cur_end {
                    // Bpcov contrast: coverage LEFT vs RIGHT of pos
                    let endpos = (pos - bundle_start) as i64;
                    let winstart = (endpos - CHI_THR + 1).max(0);
                    let winend = (endpos + CHI_THR).min(bpcov_len as i64 - 1);
                    let left_cov = get_cov_range(winstart, endpos);
                    let right_cov = get_cov_range(endpos + 1, winend);
                    let mut tmpcov = (left_cov - right_cov) / (DROP * CHI_THR as f64);

                    if tmpcov <= 0.0 && lend[le].cov < 0.0 {
                        tmpcov = ERROR_PERC;
                    }

                    if tmpcov > 0.0 {
                        // SPLIT: [cur_start, pos] and [pos+1, cur_end]
                        // In half-open: [cur_start, pos+1) and [pos+1, cur_end)
                        let split_pos = pos + 1; // half-open boundary
                        graph.nodes[cur_nid].end = split_pos;
                        graph.nodes[cur_nid].hardend = true;

                        let new_nid = graph.nodes.len();
                        let mut new_node = crate::graph::GraphNode::new(0, 0, 0);
                        new_node.node_id = new_nid;
                        new_node.start = split_pos;
                        new_node.end = cur_end;
                        graph.nodes.push(new_node);
                        graph.n_nodes = graph.nodes.len();

                        // Edge: prev → sink (hardend)
                        graph.nodes[sink_id].parents.insert_grow(cur_nid);
                        // Edge: prev → new (contiguous)
                        graph.nodes[cur_nid].children.insert_grow(new_nid);
                        graph.nodes[new_nid].parents.insert_grow(cur_nid);

                        cur_nid = new_nid;
                        splits += 1;
                    }
                }
                le += 1;
            }
        }
    }

    splits
}

/// Inline longtrim: process lstart/lend events within the current graphnode
/// up to `nodeend`. Matches longtrim() (-2740).
///
/// Splits the current graphnode at validated read boundary positions.
/// - lstart events: split at pos, new node gets source edge + hardstart.
/// - lend events: split at pos+1, left half gets sink edge + hardend.
///
/// `nls`/`nle` are advancing pointers into the sorted lstart/lend arrays.
#[allow(clippy::too_many_arguments)]
fn longtrim_inline(
    graph: &mut Graph,
    graphnode_id: &mut usize,
    nls: &mut usize,
    nle: &mut usize,
    lstart: &[ReadBoundary],
    lend: &[ReadBoundary],
    bpcov: &Bpcov,
    bpcov_stranded: Option<&BpcovStranded>,
    bundle_strand: char,
    bundle_start: u64,
    nodeend: u64,
    startcov: bool, // true if junction exists at current node start
    endcov: bool,   // true if junction exists at nodeend
    source_bid: usize,
    sink_parents: &mut Vec<usize>,
) {
    const CHI_THR: i64 = 50;
    const DROP: f64 = 0.5;
    const LONGINTRONANCHOR: u64 = 25;
    const ERROR_PERC: f64 = 0.1;

    let bpcov_len = bpcov.cov.len();
    let source_id = graph.source_id;

    // : use strand-specific coverage for longtrim.
    // the original implementation uses get_cov_sign(2*s, ...) which computes total - opposite_strand,
    // giving the coverage on the current strand only. This prevents false splits
    // where total coverage is constant but strand-specific coverage drops.
    let get_cov = |s: i64, e: i64| -> f64 {
        if s < 0 || e < s || s as usize >= bpcov_len { return 0.0; }
        let su = s as usize;
        let eu = (e as usize).min(bpcov_len - 1);
        if let Some(bps) = bpcov_stranded {
            let total = bps.get_cov_range(BPCOV_STRAND_ALL, su, eu);
            let opposite = match bundle_strand {
                '-' => bps.get_cov_range(BPCOV_STRAND_PLUS, su, eu),
                '+' => bps.get_cov_range(BPCOV_STRAND_MINUS, su, eu),
                _ => 0.0, // unstranded: use total
            };
            (total - opposite).max(0.0)
        } else {
            bpcov.get_cov_range(su, eu)
        }
    };

    // Skip events before current node.
    let gstart = graph.nodes[*graphnode_id].start;
    while *nls < lstart.len() && lstart[*nls].pos < gstart { *nls += 1; }
    while *nle < lend.len() && lend[*nle].pos < gstart { *nle += 1; }

    // Process events within [gstart, nodeend).
    while (*nls < lstart.len() && lstart[*nls].pos < nodeend)
        || (*nle < lend.len() && lend[*nle].pos < nodeend)
    {
        let use_start = if *nle >= lend.len() {
            true
        } else if *nls >= lstart.len() {
            false
        } else {
            lstart[*nls].pos <= lend[*nle].pos
        };

        // Minimum accumulated reads at this boundary to be considered for splitting.
        // Positions with only 1-2 reads are alignment jitter, not real TSS/TES.
        const MIN_BOUNDARY_READS: f64 = 3.0;
        // Minimum tmpcov (window coverage contrast) to accept a split.
        // StringTie has no threshold (splits on any positive tmpcov); keeping 25
        // here avoids false splits at ambiguous boundaries. Override via
        // RUSTLE_LONGTRIM_MIN_TMPCOV.
        let min_tmpcov: f64 = std::env::var("RUSTLE_LONGTRIM_MIN_TMPCOV")
            .ok()
            .and_then(|v| v.parse::<f64>().ok())
            .unwrap_or(25.0);

        if use_start {
            let pos = lstart[*nls].pos;
            let boundary_cov = lstart[*nls].cov;
            let cur_start = graph.nodes[*graphnode_id].start;

            // Skip boundaries with insufficient read support.
            if boundary_cov.abs() < MIN_BOUNDARY_READS {
                *nls += 1;
                continue;
            }

            // Proximity check ( startcov/endcov + longintronanchor).
            let start_ok = startcov || pos > cur_start + LONGINTRONANCHOR;
            let end_ok = endcov || pos < nodeend + LONGINTRONANCHOR;

            let mut tmpcov = 0.0;
            if start_ok && end_ok && pos > cur_start {
                let startpos = (pos - bundle_start) as i64;
                let winstart = (startpos - CHI_THR).max(0);
                let winend = (startpos + CHI_THR - 1).min(bpcov_len as i64 - 1);
                tmpcov = (get_cov(startpos, winend) - get_cov(winstart, startpos - 1))
                    / (DROP * CHI_THR as f64);
            }
            if tmpcov <= 0.0 && lstart[*nls].cov < 0.0 {
                tmpcov = ERROR_PERC;
            }
            if tmpcov > min_tmpcov {
                let cur_end = graph.nodes[*graphnode_id].end;
                if pos > graph.nodes[*graphnode_id].start && pos < cur_end {
                    // Split: [cur_start, pos) and [pos, cur_end)
                    let prev_id = *graphnode_id;
                    graph.nodes[prev_id].end = pos;
                    let new_node = graph.add_node(pos, cur_end);
                    new_node.source_bnode = Some(source_bid);
                    new_node.hardstart = true;
                    let new_id = new_node.node_id;
                    // Source → new (hardstart boundary)
                    graph.add_edge(source_id, new_id);
                    // Prev → new (contiguous)
                    graph.add_edge(prev_id, new_id);
                    // Prev → sink: when reads start at pos (new TSS), the preceding
                    // transcript should be able to terminate at pos-1. Without this
                    // edge, flow must continue through the new hardstart node, which
                    // produces gene-chimeric paths (see STRG.125 case: TSS at
                    // 22459860 with ~18 supporting read starts). Disable via
                    // RUSTLE_NO_LSTART_SINK=1.
                    if std::env::var_os("RUSTLE_NO_LSTART_SINK").is_none() {
                        graph.nodes[prev_id].hardend = true;
                        sink_parents.push(prev_id);
                    }
                    *graphnode_id = new_id;
                }
            }
            *nls += 1;
        } else {
            let pos = lend[*nle].pos;
            let boundary_cov = lend[*nle].cov;
            let cur_start = graph.nodes[*graphnode_id].start;

            if boundary_cov.abs() < MIN_BOUNDARY_READS {
                *nle += 1;
                continue;
            }

            let start_ok = !startcov || pos > cur_start + LONGINTRONANCHOR;
            let end_ok = !endcov || pos < nodeend + LONGINTRONANCHOR;

            let mut tmpcov = 0.0;
            if start_ok && end_ok && pos > cur_start {
                let endpos = (pos - bundle_start) as i64;
                let winstart = (endpos - CHI_THR + 1).max(0);
                let winend = (endpos + CHI_THR).min(bpcov_len as i64 - 1);
                tmpcov = (get_cov(winstart, endpos) - get_cov(endpos + 1, winend))
                    / (DROP * CHI_THR as f64);
            }
            if tmpcov <= 0.0 && lend[*nle].cov < 0.0 {
                tmpcov = ERROR_PERC;
            }
            if tmpcov > min_tmpcov {
                let split = pos + 1; // half-open boundary
                let cur_end = graph.nodes[*graphnode_id].end;
                if split > graph.nodes[*graphnode_id].start && split < cur_end {
                    // Split: [cur_start, split) and [split, cur_end)
                    graph.nodes[*graphnode_id].end = split;
                    graph.nodes[*graphnode_id].hardend = true;
                    let new_node = graph.add_node(split, cur_end);
                    new_node.source_bnode = Some(source_bid);
                    let new_id = new_node.node_id;
                    // Prev → new (contiguous)
                    graph.add_edge(*graphnode_id, new_id);
                    // Prev → sink (hardend)
                    sink_parents.push(*graphnode_id);
                    *graphnode_id = new_id;
                }
            }
            *nle += 1;
        }
    }
}

/// Build graph and, when enabled, apply longtrim at graph-build stage.
/// This keeps long-read boundary splitting in the create_graph flow instead of
/// as a later pipeline-only post-step.
pub fn create_graph_with_longtrim(
    junctions: &[Junction],
    bundle_start: u64,
    bundle_end: u64,
    bundlenodes: Option<&CBundlenode>,
    junction_support: u64,
    reads: Option<&[BundleRead]>,
    bundle_strand: char,
    junction_stats: Option<&JunctionStats>,
    bpcov: &Bpcov,
    bpcov_stranded: Option<&BpcovStranded>,
    lstart: &[ReadBoundary],
    lend: &[ReadBoundary],
    enable_longtrim: bool,
    longtrim_min_boundary_cov: f64,
) -> (Graph, Vec<GraphTransfrag>, CreateGraphLongtrimStats) {
    // Oracle splits: read exact split positions from a file (for debugging/validation).
    // Format: one line per split, "position type" where type is "start" or "end".
    // Positions are 0-based. Set RUSTLE_ORACLE_SPLITS=/path/to/splits.txt
    let oracle_starts_owned;
    let oracle_ends_owned;
    let (lt_starts, lt_ends): (&[ReadBoundary], &[ReadBoundary]) = if let Some(path) = std::env::var_os("RUSTLE_ORACLE_SPLITS") {
        let mut starts = Vec::new();
        let mut ends = Vec::new();
        if let Ok(content) = std::fs::read_to_string(path.to_str().unwrap_or("")) {
            for line in content.lines() {
                let parts: Vec<&str> = line.split_whitespace().collect();
                if parts.len() >= 2 {
                    if let Ok(pos) = parts[0].parse::<u64>() {
                        if pos >= bundle_start && pos <= bundle_end {
                            let b = ReadBoundary { pos, cov: -1.0 }; // cov=-1 forces split via ERROR_PERC
                            if parts[1] == "start" { starts.push(b); }
                            else { ends.push(b); }
                        }
                    }
                }
            }
        }
        starts.sort_unstable_by_key(|b| b.pos);
        ends.sort_unstable_by_key(|b| b.pos);
        oracle_starts_owned = starts;
        oracle_ends_owned = ends;
        (&oracle_starts_owned, &oracle_ends_owned)
    } else if enable_longtrim && std::env::var_os("RUSTLE_DISABLE_LONGTRIM").is_none() {
        // Pass empty externals so per-bundlenode coverage-derivative detection
        // runs inside create_graph_inner ( the original implementation generates
        // boundaries per-bundlenode, not externally per-bundle).
        (&[] as &[ReadBoundary], &[] as &[ReadBoundary])
    } else {
        (&[] as &[ReadBoundary], &[] as &[ReadBoundary])
    };
    let mut graph = create_graph_inner(
        junctions,
        bundle_start,
        bundle_end,
        bundlenodes,
        junction_support,
        reads,
        bundle_strand,
        junction_stats,
        Some(bpcov),
        bpcov_stranded,
        lt_starts,
        lt_ends,
    );

    let mut stats = CreateGraphLongtrimStats {
        lstart_events: lstart.len(),
        lend_events: lend.len(),
        ..Default::default()
    };
    let mut synthetic: Vec<GraphTransfrag> = Vec::new();

    if enable_longtrim {
        // the original algorithm's longtrim uses only chi-square coverage-based boundary detection
        // (the sliding diffval window in create_graph).  It does NOT inject
        // raw read start/end positions as feature points.  Passing lstart/lend here
        // was adding hundreds of candidate splits per locus (vs ~5-10 in the original algorithm),
        // creates node boundaries at read start/end positions
        // (lstart/lend) when coverage evidence supports a split (bpcov contrast
        // window test). Rustle's boundary collection passes ALL read boundaries
        // without the inline contrast check, causing 4-5x over-segmentation
        // . Pass empty arrays until the exact
        // contrast logic is reimplemented in apply_longtrim_direct.
        // NOTE: Iterative longtrim splits (apply_iterative_longtrim_splits) must
        // happen INSIDE create_graph, before transfrags are built. Post-creation
        // splits break transfrag patterns and edge consistency.
        // TODO: Integrate longtrim splitting into create_graph's node iteration loop.
        //
        // Pass the caller's per-bundle lstart/lend, optionally clustered
        // to avoid 4-5x over-segmentation vs the original algorithm's
        // inline-contrast check. Raw read-ends produce a split candidate per
        // unique position; clustering collapses nearby ones into one peak
        // keeping the MAX cov seen. Gate behind RUSTLE_LONGTRIM_APPLY_ON
        // so the default remains the historical empty-map (baseline safe).
        let apply_on = std::env::var_os("RUSTLE_LONGTRIM_APPLY_ON").is_some();
        let cluster_win: u64 = std::env::var("RUSTLE_LONGTRIM_APPLY_CLUSTER_WIN")
            .ok().and_then(|v| v.parse().ok()).unwrap_or(50);
        let min_cluster_cov: f64 = std::env::var("RUSTLE_LONGTRIM_APPLY_MIN_COV")
            .ok().and_then(|v| v.parse().ok()).unwrap_or(longtrim_min_boundary_cov.max(5.0));
        let cluster_max = |src: &[ReadBoundary]| -> Vec<ReadBoundary> {
            let mut out: Vec<ReadBoundary> = Vec::new();
            for b in src {
                let keep_feat = b.cov < 0.0;
                if !keep_feat && b.cov < min_cluster_cov { continue; }
                if let Some(last) = out.last_mut() {
                    if b.pos.abs_diff(last.pos) <= cluster_win {
                        if b.cov > last.cov {
                            last.pos = b.pos;
                            last.cov = b.cov;
                        }
                        continue;
                    }
                }
                out.push(*b);
            }
            out
        };
        let (lt_lstart, lt_lend): (Vec<ReadBoundary>, Vec<ReadBoundary>) = if apply_on {
            (cluster_max(lstart), cluster_max(lend))
        } else {
            (Vec::new(), Vec::new())
        };
        let boundary_map = collect_longtrim_boundary_map(
            bpcov,
            bundlenodes,
            &lt_lstart,
            &lt_lend,
            longtrim_min_boundary_cov,
        );
        let schedule = build_longtrim_bundle_schedules(
            &graph,
            junctions,
            bundlenodes,
            bundle_end,
            bundle_strand,
            junction_stats,
        );
        let (lt_synth, lt_stats) = apply_longtrim_direct(
            &mut graph,
            &boundary_map,
            &schedule,
            bpcov,
            longtrim_min_boundary_cov,
            trace_strand_index(bundle_strand),
            bundle_start,
            bundle_end,
        );
        stats.applied = true;
        stats.longtrim = lt_stats;
        stats.lstart_events = lt_stats.start_boundary_events;
        stats.lend_events = lt_stats.end_boundary_events;
        synthetic = lt_synth;
    }

    (graph, synthetic, stats)
}

/// Rebuild all transfrag patterns after graph modification (pruning, trimming).
/// transfrag patterns must reflect current graph topology for capacity network.
pub fn rebuild_transfrag_patterns(graph: &Graph, transfrags: &mut [GraphTransfrag]) {
    for tf in transfrags.iter_mut() {
        tf.rebuild_pattern(graph);
    }
}

/// Original-style graph complexity reduction pass (prune_graph_nodes intent):
/// iteratively remove lowest-coverage internal nodes and reconnect parents->children
/// until active internal node count <= allowed_nodes.
/// Also rebuilds transfrag patterns to match new graph topology.
pub fn prune_graph_nodes_with_transfrags(
    graph: &mut Graph,
    transfrags: &mut [GraphTransfrag],
    allowed_nodes: usize,
    verbose: bool,
) -> usize {
    let (removed, old_to_new) = prune_graph_nodes(graph, allowed_nodes, verbose);
    if removed > 0 && !transfrags.is_empty() {
        // Remap transfrag node_ids to match new graph indices
        for tf in transfrags.iter_mut() {
            tf.node_ids = tf
                .node_ids
                .iter()
                .filter_map(|&old_id| {
                    if old_id < old_to_new.len() {
                        let new_id = old_to_new[old_id];
                        // new_id == 0 means the node was disconnected/removed
                        if new_id != 0 || old_id == 0 {
                            Some(new_id)
                        } else {
                            None
                        }
                    } else {
                        None
                    }
                })
                .collect();
        }
        rebuild_transfrag_patterns(graph, transfrags);
    }
    removed
}

/// Original-style graph complexity reduction pass (prune_graph_nodes intent):
/// Iteratively remove lowest-coverage internal nodes and reconnect parents->children
/// until active internal node count <= allowed_nodes.
///
/// After disconnecting nodes, physically remove them from the vector
/// and remap all indices to match behavior (~3929).
///
/// Returns (removed_count, old_to_new_mapping) where mapping can be used to update
/// external data structures that reference node indices (e.g., transfrags).
pub fn prune_graph_nodes(graph: &mut Graph, allowed_nodes: usize, verbose: bool) -> (usize, Vec<usize>) {
    let source_id = graph.source_id;
    let sink_id = graph.sink_id;
    if graph.n_nodes <= 2 {
        return (0, (0..graph.n_nodes).collect());
    }

    let mut removed = 0usize;
    if allowed_nodes > 0 {
        loop {
            let active: Vec<usize> = (0..graph.n_nodes)
                .filter(|&i| i != source_id && i != sink_id)
                .filter(|&i| {
                    let n = &graph.nodes[i];
                    !n.parents.is_empty() || !n.children.is_empty()
                })
                .collect();
            if active.len() <= allowed_nodes {
                break;
            }

            let Some(nid) = active.into_iter().min_by(|&a, &b| {
                graph.nodes[a]
                    .coverage
                    .partial_cmp(&graph.nodes[b].coverage)
                    .unwrap_or(std::cmp::Ordering::Equal)
            }) else {
                break;
            };

            let parents: Vec<usize> = graph.nodes[nid].parents.ones().collect();
            let children: Vec<usize> = graph.nodes[nid].children.ones().collect();
            for &p in &parents {
                graph.remove_edge(p, nid);
            }
            for &c in &children {
                graph.remove_edge(nid, c);
            }
            for &p in &parents {
                for &c in &children {
                    if p == c {
                        continue;
                    }
                    graph.add_edge(p, c);
                }
            }
            removed += 1;
        }
    }

    for id in 1..sink_id {
        // Skip nodes whose role opts out of prune auto-attach. Overlap
        // nodes get edges programmatically; auto-attach would create
        // spurious single-node paths.
        if !graph.nodes[id].role.prune_autoattach() {
            continue;
        }
        if graph.nodes[id].parents.is_empty() {
            graph.add_edge(source_id, id);
        }
    }
    // only sink-link nodes reachable from source traversal.
    let source_reach = reachable_from_source(graph);
    for id in 1..sink_id {
        if !graph.nodes[id].role.prune_autoattach() {
            continue;
        }
        if source_reach.contains(id) && graph.nodes[id].children.is_empty() {
            graph.add_edge(id, sink_id);
        }
    }

    // Remove disconnected nodes and remap indices (~3929)
    // After pruning, some nodes may have no parents and no children (disconnected).
    // physically deletes these nodes from the vector; we must do the same.
    // Always run remapping to ensure gno matches connected node count.
    let n_nodes_before = graph.n_nodes;
    let old_to_new = remap_graph_nodes_with_mapping(graph, verbose);
    let n_nodes_after = graph.n_nodes;
    let actually_removed = n_nodes_before.saturating_sub(n_nodes_after);

    if verbose && actually_removed > 0 {
        eprintln!(
            "      [Rustle] prune_graph_nodes: removed {} low-coverage internal nodes, remapped {} disconnected nodes",
            removed, actually_removed
        );
    }
    (removed, old_to_new)
}

/// Remap graph nodes and return the old->new mapping.
/// This is used by `prune_graph_nodes_with_redirects` to update redirects.
fn remap_graph_nodes_with_mapping(
    graph: &mut Graph,
    verbose: bool,
) -> Vec<usize> {
    let source_id = graph.source_id;
    let sink_id = graph.sink_id;

    // Identify connected nodes (have at least one parent or child, or are source/sink)
    let mut connected = NodeSet::with_capacity(graph.n_nodes);
    connected.insert(source_id);
    connected.insert(sink_id);

    for i in 0..graph.n_nodes {
        if i == source_id || i == sink_id {
            continue;
        }
        let node = &graph.nodes[i];
        // Node is connected if it has edges
        if !node.parents.is_empty() || !node.children.is_empty() {
            connected.insert(i);
        }
    }

    // Count connected nodes
    let connected_count = connected.count_ones();
    if connected_count == graph.n_nodes {
        // No disconnected nodes to remove - return identity mapping
        return (0..graph.n_nodes).collect();
    }

    // Create old -> new index mapping
    let mut old_to_new: Vec<usize> = vec![0; graph.n_nodes];
    let mut new_nodes: Vec<crate::graph::GraphNode> = Vec::with_capacity(connected_count);
    let mut new_source_id = 0usize;
    let mut new_sink_id = 0usize;

    for old_id in 0..graph.n_nodes {
        if connected.contains(old_id) {
            let new_id = new_nodes.len();
            old_to_new[old_id] = new_id;

            if old_id == source_id {
                new_source_id = new_id;
            }
            if old_id == sink_id {
                new_sink_id = new_id;
            }

            // Clone the node and update its node_id
            let mut new_node = graph.nodes[old_id].clone();
            new_node.node_id = new_id;
            new_nodes.push(new_node);
        }
    }

    // Update all parent/child references in the new nodes
    for new_node in &mut new_nodes {
        // Remap parents
        let old_parents: Vec<usize> = new_node.parents.ones().collect();
        new_node.parents.clear();
        for old_parent in old_parents {
            if old_parent < old_to_new.len() && connected.contains(old_parent) {
                new_node.parents.insert_grow(old_to_new[old_parent]);
            }
        }

        // Remap children
        let old_children: Vec<usize> = new_node.children.ones().collect();
        new_node.children.clear();
        for old_child in old_children {
            if old_child < old_to_new.len() && connected.contains(old_child) {
                new_node.children.insert_grow(old_to_new[old_child]);
            }
        }
    }

    // Update gpos (edge indices) - edge IDs don't change, but node references in keys do
    let mut new_gpos: HashMap<(usize, usize), usize> = HashMap::default();
    for ((from, to), edge_id) in &graph.gpos {
        if *from < old_to_new.len() && *to < old_to_new.len()
            && connected.contains(*from) && connected.contains(*to) {
            let new_from = old_to_new[*from];
            let new_to = old_to_new[*to];
            new_gpos.insert((new_from.min(new_to), new_from.max(new_to)), *edge_id);
        }
    }

    // Replace graph nodes
    graph.nodes = new_nodes;
    graph.n_nodes = graph.nodes.len();
    graph.source_id = new_source_id;
    graph.sink_id = new_sink_id;
    graph.gpos = new_gpos;

    // Recompute reachability with new indices
    graph.compute_reachability();

    if verbose {
        eprintln!(
            "      [Rustle] remap_graph_nodes: {} -> {} nodes (source={}, sink={})",
            old_to_new.len(),
            graph.n_nodes,
            graph.source_id,
            graph.sink_id
        );
    }

    old_to_new
}

/// Same pruning pass as `prune_graph_nodes`, but returns endpoint redirects for removed nodes.
pub fn prune_graph_nodes_with_redirects(
    graph: &mut Graph,
    allowed_nodes: usize,
    verbose: bool,
) -> (usize, HashMap<usize, PruneRedirect>) {
    let source_id = graph.source_id;
    let sink_id = graph.sink_id;
    let mut redirects: HashMap<usize, PruneRedirect> = Default::default();
    if allowed_nodes == 0 || graph.n_nodes <= 2 {
        return (0, redirects);
    }

    let mut removed = 0usize;
    loop {
        let active: Vec<usize> = (0..graph.n_nodes)
            .filter(|&i| i != source_id && i != sink_id)
            .filter(|&i| {
                let n = &graph.nodes[i];
                !n.parents.is_empty() || !n.children.is_empty()
            })
            .collect();
        if active.len() <= allowed_nodes {
            break;
        }

        let Some(nid) = active.into_iter().min_by(|&a, &b| {
            graph.nodes[a]
                .coverage
                .partial_cmp(&graph.nodes[b].coverage)
                .unwrap_or(std::cmp::Ordering::Equal)
        }) else {
            break;
        };

        let parents: Vec<usize> = graph.nodes[nid].parents.ones().collect();
        let children: Vec<usize> = graph.nodes[nid].children.ones().collect();
        let upstream = parents
            .iter()
            .copied()
            .filter(|&p| p != nid)
            .max_by(|&a, &b| {
                graph.nodes[a]
                    .coverage
                    .partial_cmp(&graph.nodes[b].coverage)
                    .unwrap_or(std::cmp::Ordering::Equal)
            });
        let downstream = children
            .iter()
            .copied()
            .filter(|&c| c != nid)
            .max_by(|&a, &b| {
                graph.nodes[a]
                    .coverage
                    .partial_cmp(&graph.nodes[b].coverage)
                    .unwrap_or(std::cmp::Ordering::Equal)
            });
        redirects.insert(
            nid,
            PruneRedirect {
                upstream,
                downstream,
            },
        );

        for &p in &parents {
            graph.remove_edge(p, nid);
        }
        for &c in &children {
            graph.remove_edge(nid, c);
        }
        for &p in &parents {
            for &c in &children {
                if p == c {
                    continue;
                }
                graph.add_edge(p, c);
            }
        }
        removed += 1;
    }

    for id in 1..sink_id {
        if id >= graph.nodes.len() || !graph.nodes[id].role.prune_autoattach() {
            continue;
        }
        if graph.nodes[id].parents.is_empty() {
            graph.add_edge(source_id, id);
        }
    }
    // only sink-link nodes reachable from source traversal.
    let source_reach = reachable_from_source(graph);
    for id in 1..sink_id {
        if id >= graph.nodes.len() || !graph.nodes[id].role.prune_autoattach() {
            continue;
        }
        if source_reach.contains(id) && graph.nodes[id].children.is_empty() {
            graph.add_edge(id, sink_id);
        }
    }

    // Remove disconnected nodes and remap indices (~3929)
    // Update redirects to use new indices after remapping.
    let n_nodes_before = graph.n_nodes;
    let old_to_new = remap_graph_nodes_with_mapping(graph, verbose);
    let n_nodes_after = graph.n_nodes;

    // Update redirects to use new indices
    if n_nodes_after < n_nodes_before {
        let mut new_redirects: HashMap<usize, PruneRedirect> = Default::default();
        for (old_nid, redirect) in &redirects {
            if let Some(&new_nid) = old_to_new.get(*old_nid) {
                let new_upstream = redirect.upstream.and_then(|u| old_to_new.get(u).copied());
                let new_downstream = redirect.downstream.and_then(|d| old_to_new.get(d).copied());
                new_redirects.insert(
                    new_nid,
                    PruneRedirect {
                        upstream: new_upstream,
                        downstream: new_downstream,
                    },
                );
            }
        }
        redirects = new_redirects;
    }

    graph.reindex_edge_bits_dense();
    graph.compute_reachability();
    if verbose && removed > 0 {
        eprintln!(
            "      [Rustle] prune_graph_nodes(post-split): removed {} low-coverage internal nodes",
            removed
        );
    }
    (removed, redirects)
}
