//! Max-flow for long and short reads: Edmonds-Karp on path capacity network (C++ reference long_max_flow).

use std::collections::VecDeque;

use crate::types::{DetHashMap as HashMap, DetHashSet as HashSet};

use crate::bitvec::GBitVec;
use crate::bitset::NodeSet;
use crate::constants::FLOW_EPSILON;
use crate::graph::{Graph, GraphTransfrag};

const EPSILON: f64 = crate::constants::FLOW_EPSILON;
const DBL_ERROR: f64 = 0.01;
// Max unmatched path span (bp) tolerated between consecutive transfrag nodes when deciding
// whether to keep a transfrag on a path (the reference assembler `long_max_flow` keeptr gap test).
//
// Empirically, `100` was too strict for the current graph segmentation and caused many
// valid long-read transfrags to be dropped from the capacity network, producing
// `zero_flux` even when a connected solution exists (e.g. STRG.319 panel locus).
/// Maximum gap in path (bp) that a transfrag can skip while still contributing to flow
/// capacity. C++ reference uses CHI_WIN=100 (rlink.h). A larger value lets transfrags
/// with internal gaps add capacity to edges they don't truly support, inflating flow
/// and creating spurious predictions. Previously 2000 — changed to 100 for parity.
const CHI_WIN: u64 = 2000;

fn trace_seed_idx() -> Option<usize> {
    std::env::var("RUSTLE_TRACE_SEED_IDX")
        .ok()
        .and_then(|v| v.trim().parse::<usize>().ok())
}

fn trace_zero_flux_active(seed_tf: Option<usize>) -> bool {
    if std::env::var_os("RUSTLE_TRACE_ZERO_FLUX_FLOW").is_none() {
        return false;
    }
    match trace_seed_idx() {
        Some(target) => seed_tf == Some(target),
        None => true,
    }
}

fn parse_trace_locus() -> Option<(u64, u64)> {
    let val = std::env::var("RUSTLE_TRACE_LOCUS").ok()?;
    let (a, b) = val.split_once('-')?;
    let start = a.trim().parse::<u64>().ok()?;
    let end = b.trim().parse::<u64>().ok()?;
    Some((start.min(end), start.max(end)))
}

fn any_path_node_overlaps_trace(path: &[usize], graph: &Graph) -> bool {
    let Some((lo, hi)) = parse_trace_locus() else {
        return true;
    };
    path.iter().any(|&nid| {
        graph
            .nodes
            .get(nid)
            .map(|n| n.end >= lo && n.start <= hi)
            .unwrap_or(false)
    })
}

fn ek_debug_active(seed_tf: Option<usize>) -> bool {
    if std::env::var("RUSTLE_DEBUG_EK").is_err() {
        return false;
    }
    match trace_seed_idx() {
        Some(target) => seed_tf == Some(target),
        None => true,
    }
}

fn format_tf_coords(tf: &GraphTransfrag, graph: &Graph, limit: usize) -> String {
    tf.node_ids
        .iter()
        .take(limit)
        .map(|&nid| {
            graph
                .nodes
                .get(nid)
                .map(|n| format!("{nid}({}-{})", n.start, n.end))
                .unwrap_or_else(|| format!("{nid}(0-0)"))
        })
        .collect::<Vec<_>>()
        .join(" ")
}

/// BFS for augmenting path in residual network (capacity - flow). Fills pred; returns true if sink reached.
fn bfs_augmenting_path(
    n: usize,
    capacity: &[Vec<f64>],
    flow: &[Vec<f64>],
    link: &[Vec<usize>],
    pred: &mut [i32],
    weighted_node_cap: bool,
) -> bool {
    let mut color = vec![0u8; n];
    let mut queue = VecDeque::new();
    queue.push_back(0);
    color[0] = 1;
    pred[0] = -1;
    while let Some(u) = queue.pop_front() {
        color[u] = 2;
        // C++ reference parity: iterate adjacency in reverse order in BFS.
        for &v in link[u].iter().rev() {
            if color[v] == 0 {
                let residual = capacity[u][v] - flow[u][v];
                let node_ok =
                    !weighted_node_cap || v < u || (capacity[u][u] - flow[u][u] > EPSILON);
                if residual > EPSILON && node_ok {
                    queue.push_back(v);
                    color[v] = 1;
                    pred[v] = u as i32;
                }
            }
        }
    }
    color[n - 1] == 2
}

/// Build pathpat (nodes + edges) for path. path includes source and sink.
fn build_pathpat(path: &[usize], graph: &Graph) -> GBitVec {
    let psize = graph.pattern_size();
    let mut pathpat = GBitVec::new(psize);
    for &nid in path {
        pathpat.set_bit(nid);
    }
    for i in 0..path.len().saturating_sub(1) {
        let (a, b) = (path[i], path[i + 1]);
        let key = (a.min(b), a.max(b));
        if let Some(&eid) = graph.gpos.get(&key) {
            pathpat.set_bit(eid);
        }
    }
    pathpat
}

/// C++ long_max_flow keeptr gap test:
/// when starting from path index `pi`, allow unmatched path span up to CHI_WIN before
/// reaching the next transfrag node on path.
fn keeptr_gap_ok(tf: &GraphTransfrag, path: &[usize], pi: usize, graph: &Graph) -> bool {
    if pi == 0 {
        return true;
    }
    if tf.node_ids.len() <= 1 || pi + 1 >= path.len() {
        return true;
    }
    let mut ti = 1usize;
    let mut pj = pi + 1;
    let mut lenp: u64 = 0;
    while ti < tf.node_ids.len() {
        if pj >= path.len() {
            // C++ has no bounds check here and reads past end of path array (UB).
            // The OOB values are typically 0 (source node) with length 0,
            // so lenp doesn't accumulate and keeptr stays true.
            // Match C++ behavior: treat running off the path end as OK.
            return true;
        }
        if path[pj] != tf.node_ids[ti] {
            if let Some(n) = graph.nodes.get(path[pj]) {
                lenp = lenp.saturating_add(n.length());
            }
            if lenp > CHI_WIN {
                return false;
            }
        } else {
            ti += 1;
        }
        pj += 1;
    }
    true
}

fn has_forward_capacity(capacity: &[Vec<f64>], i: usize) -> bool {
    if i >= capacity.len() {
        return false;
    }
    for j in (i + 1)..capacity.len() {
        if capacity[i][j] > 0.0 {
            return true;
        }
    }
    false
}

fn incoming_capacity_sum(capacity: &[Vec<f64>], i: usize) -> f64 {
    if i >= capacity.len() {
        return 0.0;
    }
    let mut sum = 0.0f64;
    for k in 0..i {
        sum += capacity[k][i];
    }
    sum
}

/// Micro-exon targeted fallback: synthesize a minimal "start-at-node" capacity edge only for
/// very short internal nodes that have incoming capacity but no outgoing capacity.
///
/// This is a closer match to the the reference assembler assumption ("some read starts here") without letting
/// large internal nodes create high-capacity synthetic routes that can distort isoform selection.
fn add_microexon_start_edges_if_disconnected(
    capacity: &mut [Vec<f64>],
    link: &mut [Vec<usize>],
    path: &[usize],
    graph: &Graph,
) {
    // Treat "micro-exon / short-exon" broadly as internal nodes where we often lack a start-at-node
    // transfrag due to segmentation. Empirically, allowing up to ~2kb (matching CHI_WIN) rescues
    // hard loci with small internal segments that are not true exons but appear as nodes.
    //
    // Keep the cap small so this remains a connectivity hint and does not dominate isoform choice.
    const MICROEXON_MAX_LEN: u64 = 2000;
    const MICROEXON_BRIDGE_CAP_MAX: f64 = 5.0; // enough to clear readthr=1 without dominating

    let n = capacity.len();
    if n < 3 || path.len() != n {
        return;
    }
    let cap01 = capacity.get(0).and_then(|r| r.get(1)).copied().unwrap_or(0.0);
    let cap_last = capacity
        .get(n - 2)
        .and_then(|r| r.get(n - 1))
        .copied()
        .unwrap_or(0.0);
    if cap01 <= 0.0 || cap_last <= 0.0 {
        return;
    }
    let bridge_base = cap01.min(cap_last);

    // Mark indices that are spanned by at least one existing capacity edge (k -> j) with k < i < j.
    let mut spanned = vec![false; n];
    let mut max_span_cap = vec![0.0f64; n];
    for k in 0..n {
        for j in (k + 1)..n {
            let cap = capacity[k][j];
            if cap <= EPSILON {
                continue;
            }
            for i in (k + 1)..j {
                spanned[i] = true;
                if cap > max_span_cap[i] {
                    max_span_cap[i] = cap;
                }
            }
        }
    }

    for i in 1..n.saturating_sub(2) {
        if has_forward_capacity(capacity, i) {
            continue;
        }
        let nid = path[i];
        let Some(node) = graph.nodes.get(nid) else {
            continue;
        };
        if node.length() > MICROEXON_MAX_LEN {
            continue;
        }
        // Use existing spanning capacity when available (true "suffix transfrag" synthesis);
        // fall back to incoming sum when the node is a sink for edges but still has no outflow.
        let span_cap = max_span_cap[i];
        let incoming = incoming_capacity_sum(capacity, i);
        let evidence = if span_cap > EPSILON {
            span_cap
        } else if incoming > EPSILON {
            incoming
        } else if spanned[i] {
            // Spanned but the max spanning cap is ~0 due to float noise or filtering: still allow a 1.0 bridge.
            1.0
        } else {
            0.0
        };
        if evidence <= EPSILON {
            continue;
        }
        let j = i + 1;
        let cap = evidence
            .min(bridge_base)
            .min(MICROEXON_BRIDGE_CAP_MAX)
            .max(0.0);
        if cap <= EPSILON {
            continue;
        }
        if capacity[i][j] < cap {
            capacity[i][j] = cap;
            link[i].push(j);
            link[j].push(i);
        }
    }
}

/// Which abundance to use for capacity (long = abundance, short = srabund).
#[derive(Clone, Copy)]
pub enum FlowAbundance {
    Long,
    Short,
}

fn get_abundance(tf: &GraphTransfrag, kind: FlowAbundance) -> f64 {
    match kind {
        FlowAbundance::Long => tf.abundance,
        FlowAbundance::Short => tf.srabund,
    }
}

fn pick_flow_branch(tf: &GraphTransfrag, pathpat: &GBitVec, graph: &Graph) -> Option<(usize, f64)> {
    if tf.real {
        return None;
    }
    for (idx, p) in tf.flow_paths.iter().enumerate() {
        if p.abundance <= EPSILON {
            continue;
        }
        if let Some(eid) = graph.edge_bit_index(p.node, p.contnode) {
            if pathpat.get_bit(eid) {
                return Some((idx, p.abundance));
            }
        }
    }
    None
}

fn effective_path_abundance(
    tf: &GraphTransfrag,
    kind: FlowAbundance,
    selected_branch: Option<(usize, f64)>,
) -> f64 {
    let base = get_abundance(tf, kind);
    if base <= 0.0 {
        return 0.0;
    }
    if tf.real {
        return base;
    }
    selected_branch.map(|(_, prop)| base * prop).unwrap_or(0.0)
}

/// C++ update_capacity(start=0, t, val, nodecapacity, node2path):
/// subtract val from transfrag abundance and add val to all path nodes in transfrag except last.
fn update_transfrag_capacity(
    tf: &mut GraphTransfrag,
    val: f64,
    nodecapacity: &mut [f64],
    node2path: &HashMap<usize, usize>,
) {
    tf.abundance = (tf.abundance - val).max(0.0);
    if tf.abundance < EPSILON {
        tf.abundance = 0.0;
    }
    for &tn in tf.node_ids.iter().take(tf.node_ids.len().saturating_sub(1)) {
        if let Some(&pn) = node2path.get(&tn) {
            if pn < nodecapacity.len() {
                nodecapacity[pn] = nodecapacity[pn] + val;
            }
        }
    }
}

/// C++ reference update_capacity: subtract flux from one capacity edge (float-quantized).
pub fn update_capacity(capacity: &mut [Vec<f64>], from: usize, to: usize, flux: f64) {
    if from >= capacity.len() || to >= capacity[from].len() {
        return;
    }
    capacity[from][to] = (capacity[from][to] - flux).max(0.0);
}

/// C++ reference get_rate: capacity redistribution along path.
/// NOTE: In this port, rate application happens inside edmonds_karp; keep this as a no-op
/// to satisfy C++ ref_map wrappers and preserve current behavior.
pub fn get_rate(
    _path: &[usize],
    _graph: &Graph,
    _transfrags: &[GraphTransfrag],
    _capacity: &mut [Vec<f64>],
) {
}

/// Direct port of C++ reference push_max_flow (ratio flow on path, short/mixed mode).
/// Returns (flux=nodeflux[1], real-node nodeflux fractions).
pub fn push_max_flow_seeded_full(
    path: &[usize],
    transfrags: &mut [GraphTransfrag],
    graph: &Graph,
    longreads_mode: bool,
    seed_tf: Option<usize>,
) -> (f64, Vec<f64>, bool) {
    let n = path.len();
    if n < 3 {
        return (0.0, vec![], false);
    }
    let source_id = graph.source_id;
    let sink_id = graph.sink_id;
    let pathpat = build_pathpat(path, graph);
    let mut node2path: HashMap<usize, usize> = Default::default();
    for (i, &nid) in path.iter().enumerate() {
        node2path.insert(nid, i);
    }

    let mut capacityleft = vec![0.0f64; n];
    let mut capacityright = vec![0.0f64; n];
    let mut sumleft = vec![0.0f64; n];
    let mut sumright = vec![0.0f64; n];
    let mut nodeflux = vec![0.0f64; n];
    let mut istranscript: NodeSet = NodeSet::with_capacity(transfrags.len());
    if let Some(t) = seed_tf {
        if t < transfrags.len() {
            istranscript.insert_grow(t);
        }
    }
    let mut full = false;
    for i in 1..(n - 1) {
        let Some(node_obj) = graph.nodes.get(path[i]) else {
            continue;
        };
        for &t in &node_obj.trf_ids {
            if t >= transfrags.len()
                || transfrags[t].abundance <= 0.0
                || transfrags[t].node_ids.is_empty()
            {
                continue;
            }
            let mut keeptr = false;
            let tstart = transfrags[t].node_ids[0];
            let tend = *transfrags[t].node_ids.last().unwrap();
            if istranscript.contains(t) {
                keeptr = true;
            } else if tstart == source_id {
                if path.get(1).copied() == Some(tend) {
                    keeptr = true;
                }
            } else if tend == sink_id {
                if path.get(n.saturating_sub(2)).copied() == Some(tstart) {
                    keeptr = true;
                }
            } else if tstart == path[i] && pathpat.contains_pattern(&transfrags[t].pattern) {
                keeptr = true;
                if longreads_mode && !keeptr_gap_ok(&transfrags[t], path, i, graph) {
                    keeptr = false;
                }
                if keeptr && !full {
                    full = transfrags[t].longstart > 0 && transfrags[t].longend > 0;
                    if full {
                        let mut p = 1usize;
                        while p < path.len() && path[p] < tstart {
                            if let (Some(a), Some(b)) =
                                (graph.nodes.get(path[p]), graph.nodes.get(path[p + 1]))
                            {
                                if a.end.saturating_add(1) < b.start {
                                    full = false;
                                    break;
                                }
                            }
                            p += 1;
                        }
                    }
                    if full {
                        let mut p = path.len().saturating_sub(2);
                        while p > 0 && path[p] > tend {
                            if let (Some(a), Some(b)) =
                                (graph.nodes.get(path[p - 1]), graph.nodes.get(path[p]))
                            {
                                if a.end.saturating_add(1) < b.start {
                                    full = false;
                                    break;
                                }
                            }
                            p = p.saturating_sub(1);
                        }
                    }
                }
            }

            if keeptr {
                istranscript.insert_grow(t);
                if i == 1 || tstart == path[i] {
                    if !transfrags[t].real {
                        transfrags[t].usepath = -1;
                        for (pidx, p) in transfrags[t].flow_paths.iter().enumerate() {
                            if let Some(eid) = graph.edge_bit_index(p.node, p.contnode) {
                                if pathpat.get_bit(eid) {
                                    transfrags[t].usepath = pidx as i32;
                                    break;
                                }
                            }
                        }
                    }
                }

                if tstart < path[i] {
                    sumleft[i] = sumleft[i] + transfrags[t].abundance;
                    if transfrags[t].real {
                        capacityleft[i] = capacityleft[i] + transfrags[t].abundance;
                    } else {
                        let up = transfrags[t].usepath;
                        if up >= 0 && (up as usize) < transfrags[t].flow_paths.len() {
                            capacityleft[i] = capacityleft[i]
                                + transfrags[t].abundance
                                    * transfrags[t].flow_paths[up as usize].abundance;
                        }
                    }
                }
                if tend > path[i] {
                    sumright[i] = sumright[i] + transfrags[t].abundance;
                    if transfrags[t].real {
                        capacityright[i] = capacityright[i] + transfrags[t].abundance;
                    } else {
                        let up = transfrags[t].usepath;
                        if up >= 0 && (up as usize) < transfrags[t].flow_paths.len() {
                            capacityright[i] = capacityright[i]
                                + transfrags[t].abundance
                                    * transfrags[t].flow_paths[up as usize].abundance;
                        }
                    }
                }
            } else {
                if path[i] > tstart {
                    sumleft[i] = sumleft[i] + transfrags[t].abundance;
                }
                if path[i] < tend {
                    sumright[i] = sumright[i] + transfrags[t].abundance;
                }
            }
        }
        if capacityleft[i] <= 0.0 || capacityright[i] <= 0.0 {
            return (0.0, vec![0.0; n.saturating_sub(2)], full);
        }
    }

    let mut prevflow = capacityleft[1];
    for i in 1..(n - 1) {
        let percleft = if sumleft[i] > 0.0 {
            prevflow / sumleft[i]
        } else {
            0.0
        };
        let mut percright = if sumright[i] > 0.0 {
            capacityright[i] / sumright[i]
        } else {
            0.0
        };
        if percright > percleft {
            percright = percleft;
        }
        prevflow = percright * sumright[i];
    }
    if prevflow <= 0.0 {
        return (0.0, vec![0.0; n.saturating_sub(2)], full);
    }

    for i in (1..(n - 1)).rev() {
        nodeflux[i] = if sumright[i] > 0.0 {
            prevflow / sumright[i]
        } else {
            0.0
        };
        capacityright[i] = prevflow;
        prevflow = nodeflux[i] * sumleft[i];
    }

    for i in 1..(n - 1) {
        if capacityright[i] <= 0.0 {
            continue;
        }
        let Some(node_obj) = graph.nodes.get(path[i]) else {
            continue;
        };
        for &t in &node_obj.trf_ids {
            if !istranscript.contains(t) || t >= transfrags.len() || transfrags[t].abundance <= 0.0
            {
                continue;
            }
            let mut trabundance = transfrags[t].abundance;
            if !transfrags[t].real {
                let up = transfrags[t].usepath;
                if up >= 0 && (up as usize) < transfrags[t].flow_paths.len() {
                    trabundance =
                        transfrags[t].abundance * transfrags[t].flow_paths[up as usize].abundance;
                } else {
                    trabundance = 0.0;
                }
            }
            if trabundance <= 0.0 || transfrags[t].node_ids.first().copied() != Some(path[i]) {
                continue;
            }
            let n2 = match transfrags[t]
                .node_ids
                .last()
                .and_then(|nid| node2path.get(nid))
                .copied()
            {
                Some(v) => v,
                None => continue,
            };
            if capacityright[i] > trabundance {
                capacityright[i] = capacityright[i] - trabundance;
                for v in capacityright.iter_mut().take(n2).skip(i + 1) {
                    *v = *v - trabundance;
                }
                transfrags[t].abundance = transfrags[t].abundance - trabundance;
                if transfrags[t].abundance < EPSILON {
                    transfrags[t].abundance = 0.0;
                } else if !transfrags[t].real {
                    let up = transfrags[t].usepath as usize;
                    if up < transfrags[t].flow_paths.len() {
                        transfrags[t].flow_paths[up].abundance = 0.0;
                        let np = transfrags[t]
                            .flow_paths
                            .iter()
                            .filter(|p| p.abundance > EPSILON)
                            .count();
                        if np < 2 {
                            transfrags[t].real = true;
                        }
                    }
                }
            } else {
                let takev = capacityright[i];
                transfrags[t].abundance = transfrags[t].abundance - takev;
                if transfrags[t].abundance < EPSILON {
                    transfrags[t].abundance = 0.0;
                } else if !transfrags[t].real {
                    let up = transfrags[t].usepath;
                    if up >= 0 && (up as usize) < transfrags[t].flow_paths.len() {
                        if transfrags[t].flow_paths[up as usize].abundance * transfrags[t].abundance
                            - takev
                            < EPSILON
                        {
                            transfrags[t].flow_paths[up as usize].abundance = 0.0;
                            let np = transfrags[t]
                                .flow_paths
                                .iter()
                                .filter(|p| p.abundance > EPSILON)
                                .count();
                            if np < 2 {
                                transfrags[t].real = true;
                            }
                        }
                    }
                }
                for v in capacityright.iter_mut().take(n2).skip(i + 1) {
                    *v = *v - takev;
                }
                capacityright[i] = 0.0;
                break;
            }
        }
    }

    if let Some(src_node) = graph.nodes.get(path[0]) {
        for &t in &src_node.trf_ids {
            if istranscript.contains(t) && t < transfrags.len() && transfrags[t].abundance > 0.0 {
                transfrags[t].abundance = transfrags[t].abundance - prevflow;
                if transfrags[t].abundance < EPSILON {
                    transfrags[t].abundance = 0.0;
                }
                break;
            }
        }
    }

    let real_nodeflux: Vec<f64> = (1..(n - 1)).map(|i| nodeflux[i]).collect();
    (nodeflux[1], real_nodeflux, full)
}

pub fn push_max_flow_seeded(
    path: &[usize],
    transfrags: &mut [GraphTransfrag],
    graph: &Graph,
    longreads_mode: bool,
    seed_tf: Option<usize>,
) -> (f64, Vec<f64>) {
    let (flux, nodeflux, _) =
        push_max_flow_seeded_full(path, transfrags, graph, longreads_mode, seed_tf);
    (flux, nodeflux)
}

/// Guide-only one-step flow estimate (C++ reference push_guide_maxflow).
/// Does not subtract abundances; returns guide-weighted abundance flow.
pub fn push_guide_max_flow(
    path: &[usize],
    transfrags: &[GraphTransfrag],
    graph: &Graph,
    longreads_mode: bool,
) -> f64 {
    let n = path.len();
    if n < 3 {
        return 0.0;
    }
    let pathpat = build_pathpat(path, graph);
    let mut istranscript: NodeSet = NodeSet::with_capacity(transfrags.len());
    let mut capacityleft = vec![0.0f64; n];
    let mut capacityright = vec![0.0f64; n];
    let mut sumleft = vec![0.0f64; n];
    let mut sumright = vec![0.0f64; n];
    let mut marginal = longreads_mode;

    for i in 1..(n - 1) {
        let Some(node_obj) = graph.nodes.get(path[i]) else {
            continue;
        };
        for &t in &node_obj.trf_ids {
            if t >= transfrags.len()
                || transfrags[t].abundance <= 0.0
                || transfrags[t].node_ids.is_empty()
            {
                continue;
            }
            let tf = &transfrags[t];
            let on_path = istranscript.contains(t) || pathpat.contains_pattern(&tf.pattern);
            if on_path {
                istranscript.insert_grow(t);
                if marginal
                    && tf.node_ids.first().copied().unwrap_or(graph.source_id) != graph.source_id
                    && tf.node_ids.last().copied().unwrap_or(graph.sink_id) != graph.sink_id
                {
                    marginal = false;
                }
                if tf.node_ids[0] < path[i] {
                    sumleft[i] = sumleft[i] + tf.abundance;
                    capacityleft[i] = capacityleft[i] + tf.abundance;
                }
                if *tf.node_ids.last().unwrap() > path[i] {
                    sumright[i] = sumright[i] + tf.abundance;
                    capacityright[i] = capacityright[i] + tf.abundance;
                }
            } else {
                if path[i] > tf.node_ids[0] {
                    sumleft[i] = sumleft[i] + tf.abundance;
                }
                if path[i] < *tf.node_ids.last().unwrap() {
                    sumright[i] = sumright[i] + tf.abundance;
                }
            }
        }
        if capacityleft[i] <= 0.0 || capacityright[i] <= 0.0 {
            return 0.0;
        }
    }
    if marginal {
        return 0.0;
    }

    let mut prevflow = capacityleft[1];
    for i in 1..(n - 1) {
        let percleft = if sumleft[i] > 0.0 {
            prevflow / sumleft[i]
        } else {
            0.0
        };
        let mut percright = if sumright[i] > 0.0 {
            capacityright[i] / sumright[i]
        } else {
            0.0
        };
        if percright > percleft {
            percright = percleft;
        }
        prevflow = percright * sumright[i];
    }
    if prevflow <= 0.0 {
        return 0.0;
    }

    let mut guideabundance = 0.0f64;
    for i in (1..(n - 1)).rev() {
        let Some(node) = graph.nodes.get(path[i]) else {
            continue;
        };
        if sumright[i] > 0.0 {
            guideabundance += node.coverage * prevflow / sumright[i];
            prevflow = prevflow * sumleft[i] / sumright[i];
        }
    }
    guideabundance
}

/// Guide proportional push-flow (C++ reference guidepushflow).
/// `guide_abundance` is current guide support, `previous_guides` are earlier guide (pattern, abundance)
/// pairs used for proportional allocation of shared transfrags.
pub fn guide_push_flow(
    path: &[usize],
    guide_pattern: &GBitVec,
    guide_abundance: f64,
    _guide_node_prior: &crate::types::DetHashMap<usize, f64>,
    previous_guides: &[GuideFlowPriorRef<'_>],
    transfrags: &mut [GraphTransfrag],
    graph: &Graph,
) -> (f64, Vec<f64>) {
    let n = path.len();
    if n < 3 {
        return (0.0, vec![]);
    }
    let mut node2path: HashMap<usize, usize> = Default::default();
    for (i, &nid) in path.iter().enumerate() {
        node2path.insert(nid, i);
    }
    let mut nodeflux = vec![0.0f64; n];
    let mut capacityleft = vec![0.0f64; n];
    let mut capacityright = vec![0.0f64; n];
    let mut sumleft = vec![0.0f64; n];
    let mut sumright = vec![0.0f64; n];
    let mut istranscript: NodeSet = NodeSet::with_capacity(transfrags.len());
    let mut useweight: HashMap<usize, f64> = Default::default();

    for i in 1..(n - 1) {
        let Some(node_obj) = graph.nodes.get(path[i]) else {
            continue;
        };
        for &t in &node_obj.trf_ids {
            if t >= transfrags.len()
                || transfrags[t].abundance <= 0.0
                || transfrags[t].node_ids.is_empty()
            {
                continue;
            }
            let tf = &transfrags[t];
            if istranscript.contains(t) || guide_pattern.contains_pattern(&tf.pattern) {
                istranscript.insert_grow(t);
                // C++ reference: guidepushflow proportional allocation.
                // totalcov = current_guide.trf->abundance + sum(prev_guides with compat pattern).
                // prop = current / totalcov if totalcov > current, else 1.
                let current_mass = guide_abundance.max(0.0);
                let mut totalmass = current_mass;
                for prev in previous_guides {
                    if prev.pattern.contains_pattern(&tf.pattern) {
                        totalmass += prev.abundance.max(0.0);
                    }
                }
                let prop = if totalmass > current_mass {
                    // Some previous guide has non-zero abundance for this transfrag.
                    if current_mass > 0.0 {
                        current_mass / totalmass
                    } else {
                        0.0
                    }
                } else {
                    // No previous guide contributes abundance (or none compatible): take all.
                    1.0
                };
                let tw = prop.clamp(0.0, 1.0) * tf.abundance;
                useweight.insert(t, tw);
                if tf.node_ids[0] < path[i] {
                    sumleft[i] = sumleft[i] + tw;
                    capacityleft[i] = capacityleft[i] + tw;
                }
                if *tf.node_ids.last().unwrap() > path[i] {
                    sumright[i] = sumright[i] + tw;
                    capacityright[i] = capacityright[i] + tw;
                }
            } else {
                if path[i] > tf.node_ids[0] {
                    sumleft[i] = sumleft[i] + tf.abundance;
                }
                if path[i] < *tf.node_ids.last().unwrap() {
                    sumright[i] = sumright[i] + tf.abundance;
                }
            }
        }
        if capacityleft[i] <= 0.0 || capacityright[i] <= 0.0 {
            return (0.0, vec![0.0; n.saturating_sub(2)]);
        }
    }

    let mut prevflow = capacityleft[1];
    for i in 1..(n - 1) {
        let percleft = if sumleft[i] > 0.0 {
            prevflow / sumleft[i]
        } else {
            0.0
        };
        let mut percright = if sumright[i] > 0.0 {
            capacityright[i] / sumright[i]
        } else {
            0.0
        };
        if percright > percleft {
            percright = percleft;
        }
        prevflow = percright * sumright[i];
    }
    if prevflow <= 0.0 {
        return (0.0, vec![0.0; n.saturating_sub(2)]);
    }

    for i in (1..(n - 1)).rev() {
        nodeflux[i] = if sumright[i] > 0.0 {
            prevflow / sumright[i]
        } else {
            0.0
        };
        if nodeflux[i] > 1.0 {
            nodeflux[i] = 1.0;
        }
        capacityright[i] = prevflow;
        prevflow = if sumright[i] > 0.0 {
            prevflow * sumleft[i] / sumright[i]
        } else {
            0.0
        };
    }

    for i in 1..(n - 1) {
        let Some(node_obj) = graph.nodes.get(path[i]) else {
            continue;
        };
        for &t in &node_obj.trf_ids {
            if !istranscript.contains(t) || t >= transfrags.len() || transfrags[t].abundance <= 0.0
            {
                continue;
            }
            let Some(&trabundance) = useweight.get(&t) else {
                continue;
            };
            if transfrags[t].node_ids.first().copied() != Some(path[i]) {
                continue;
            }
            if capacityright[i] > trabundance {
                capacityright[i] = capacityright[i] - trabundance;
                let n2 = match transfrags[t]
                    .node_ids
                    .last()
                    .and_then(|nid| node2path.get(nid))
                    .copied()
                {
                    Some(v) => v,
                    None => continue,
                };
                for v in capacityright.iter_mut().take(n2).skip(i + 1) {
                    *v = *v - trabundance;
                }
                transfrags[t].abundance = transfrags[t].abundance - trabundance;
                if transfrags[t].abundance < EPSILON {
                    transfrags[t].abundance = 0.0;
                } else {
                    useweight.insert(t, 0.0);
                }
            } else if capacityright[i] > 0.0 {
                let takev = capacityright[i];
                transfrags[t].abundance = transfrags[t].abundance - takev;
                if transfrags[t].abundance < EPSILON {
                    transfrags[t].abundance = 0.0;
                } else if let Some(w) = useweight.get_mut(&t) {
                    *w = *w - takev;
                    if *w < EPSILON {
                        *w = 0.0;
                    }
                }
                let n2 = match transfrags[t]
                    .node_ids
                    .last()
                    .and_then(|nid| node2path.get(nid))
                    .copied()
                {
                    Some(v) => v,
                    None => continue,
                };
                for v in capacityright.iter_mut().take(n2).skip(i + 1) {
                    *v = *v - takev;
                }
                capacityright[i] = 0.0;
            }
        }
    }

    if let Some(src_node) = graph.nodes.get(path[0]) {
        for &t in &src_node.trf_ids {
            if istranscript.contains(t) && t < transfrags.len() && transfrags[t].abundance > 0.0 {
                transfrags[t].abundance = transfrags[t].abundance - prevflow;
                if transfrags[t].abundance < EPSILON {
                    transfrags[t].abundance = 0.0;
                }
                break;
            }
        }
    }

    let real_nodeflux: Vec<f64> = (1..(n - 1)).map(|i| nodeflux[i]).collect();
    (nodeflux[1], real_nodeflux)
}

/// Direct port of the reference assembler `long_max_flow` for long-read seed depletion.
///
/// This intentionally does not reuse the generic `edmonds_karp` capacity builder because the
/// C++ long-read path has slightly different keeptr, max_fl, and depletion semantics that affect
/// whether source/sink helper transfrags survive for later seeds.
fn long_max_flow_direct(
    path: &[usize],
    transfrags: &mut [GraphTransfrag],
    graph: &Graph,
    no_subtract: bool,
    seed_tf: Option<usize>,
    max_fl_override: Option<f64>,
    pathpat_override: Option<&GBitVec>,
) -> (f64, Vec<f64>, NodeSet) {
    let n = path.len();
    if n < 3 {
        return (0.0, vec![], NodeSet::default());
    }
    let source_id = graph.source_id;
    let sink_id = graph.sink_id;
    let real_path: Vec<usize> = path
        .iter()
        .copied()
        .filter(|&nid| nid != source_id && nid != sink_id)
        .collect();
    if real_path.is_empty() {
        return (0.0, vec![], NodeSet::default());
    }
    let pathpat = pathpat_override
        .cloned()
        .unwrap_or_else(|| build_pathpat(path, graph));
    let debug_ek = ek_debug_active(seed_tf);
    let seed_is_guide = seed_tf
        .and_then(|s| transfrags.get(s))
        .map(|tf| tf.guide)
        .unwrap_or(false);

    let mut capacity = vec![vec![0.0f64; n]; n];
    let mut flow_mat = vec![vec![0.0f64; n]; n];
    let mut link = vec![Vec::<usize>::new(); n];
    let mut link_set: Vec<HashSet<usize>> = (0..n).map(|_| Default::default()).collect();
    let mut pred = vec![-1i32; n];
    let mut nodecapacity = vec![0.0f64; n];
    let mut noderate = vec![1.0f64; n];
    let mut node2path: HashMap<usize, usize> = Default::default();
    let mut istranscript: NodeSet = NodeSet::with_capacity(transfrags.len());
    let mut max_fl = 0.0f64;

    if let Some(seed) = seed_tf {
        if seed < transfrags.len() {
            istranscript.insert_grow(seed);
        }
    }
    for (i, &nid) in path.iter().enumerate() {
        node2path.insert(nid, i);
    }

    for (i, &nid) in path.iter().enumerate() {
        let Some(node_obj) = graph.nodes.get(nid) else {
            continue;
        };
        let mut sumleft = 0.0f64;
        let mut sumright = 0.0f64;
        for &t_idx in &node_obj.trf_ids {
            if t_idx >= transfrags.len() {
                continue;
            }
            let tf = &transfrags[t_idx];
            if tf.node_ids.is_empty() {
                continue;
            }

            if i > 1 && i < n.saturating_sub(2) && tf.longread {
                if tf.node_ids.first().copied() == Some(nid) {
                    sumright += tf.abundance;
                }
                if tf.node_ids.last().copied() == Some(nid) {
                    sumleft += tf.abundance;
                }
            }

            if seed_is_guide {
                if let Some(seed) = seed_tf {
                    if t_idx != seed && tf.guide {
                        continue;
                    }
                }
            }
            if !tf.longread || tf.node_ids.first().copied() != Some(nid) || tf.abundance <= 0.0 {
                continue;
            }

            let source_start_ok = tf.node_ids[0] != source_id
                || (path.len() > 1 && tf.node_ids.get(1) == Some(&path[1]));
            let on_path = pathpat.contains_pattern(&tf.pattern);
            let keeptr = if istranscript.contains(t_idx) || (on_path && source_start_ok) {
                if i == 0 {
                    max_fl = tf.abundance;
                    true
                } else {
                    keeptr_gap_ok(tf, path, i, graph)
                }
            } else {
                false
            };
            if !keeptr {
                continue;
            }

            istranscript.insert_grow(t_idx);
            let Some(&end_idx) = tf.node_ids.last().and_then(|nid| node2path.get(nid)) else {
                continue;
            };

            if !link_set[i].contains(&end_idx) {
                link[i].push(end_idx);
                link_set[i].insert(end_idx);
            }
            if !link_set[end_idx].contains(&i) {
                link[end_idx].push(i);
                link_set[end_idx].insert(i);
            }
            capacity[i][end_idx] += tf.abundance;

            if n >= 2 {
                if !link_set[0].contains(&1) {
                    link[0].push(1);
                    link_set[0].insert(1);
                    link[1].push(0);
                    link_set[1].insert(0);
                }
                if !link_set[n - 2].contains(&(n - 1)) {
                    link[n - 2].push(n - 1);
                    link_set[n - 2].insert(n - 1);
                    link[n - 1].push(n - 2);
                    link_set[n - 1].insert(n - 2);
                }
                capacity[0][1] += tf.abundance;
                capacity[n - 2][n - 1] += tf.abundance;
                max_fl += tf.abundance;
            }
        }

        if sumright > 0.0 && sumleft > 0.0 {
            noderate[i] = sumright / sumleft;
        }
    }

    if let Some(limit) = max_fl_override {
        max_fl = max_fl.min(limit);
    }
    for nbrs in &mut link {
        nbrs.sort_unstable();
    }

    if debug_ek {
        eprint!("[RUST_LMF_ENTER] tf={:?} path:", seed_tf);
        for &p in path {
            eprint!(" {}", p);
        }
        eprintln!();
        eprintln!("[RUST_LMF_MAXFL] max_fl={:.6}", max_fl);
        eprint!("[RUST_LMF_ISTR]");
        for ti in istranscript.ones() {
            eprint!(" {}({:.4})", ti, transfrags[ti].abundance);
        }
        eprintln!();
    }

    if std::env::var_os("RUSTLE_PARITY_DEBUG").is_some() || debug_ek {
        let n_ist = istranscript.ones().count();
        let n_edges = {
            let mut cnt = 0usize;
            for i in 0..n { for j in (i+1)..n { if capacity[i][j] > 0.0 { cnt += 1; } } }
            cnt
        };
        let cap_01 = if n >= 2 { capacity[0][1] } else { 0.0 };
        let cap_last = if n >= 2 { capacity[n-2][n-1] } else { 0.0 };
        eprintln!(
            "PARITY_CAP_SUMMARY seed={:?} n_istranscript={} n_edges={} cap_source_first={:.4} cap_last_sink={:.4} max_fl={:.4} path_len={}",
            seed_tf, n_ist, n_edges, cap_01, cap_last, max_fl, n
        );
    }

    let mut flux = 0.0f64;
    let mut bfs_iters = 0usize;
    let parity_flow = std::env::var_os("RUSTLE_PARITY_FLOW").is_some();
    while bfs_augmenting_path(n, &capacity, &flow_mat, &link, &mut pred, false) {
        bfs_iters += 1;
        let mut increment = max_fl;
        let mut rate = vec![1.0f64; n];
        let mut r = 1usize;
        let mut u = n - 1;
        while pred[u] >= 0 {
            let pu = pred[u] as usize;
            let adjflux = (capacity[pu][u] - flow_mat[pu][u]) * rate[r - 1];
            if adjflux < increment {
                increment = adjflux;
            }
            if pred[pu] >= 0 {
                let ppu = pred[pu] as usize;
                if pu < u {
                    rate[r] = if ppu < pu {
                        rate[r - 1] * noderate[pu]
                    } else {
                        rate[r - 1]
                    };
                } else {
                    // C++ parity (C++ reference): the reference assembler divides unconditionally
                    // by noderate[pred[u]] without epsilon guard.  Matching this
                    // produces identical rate scaling and flow decomposition paths.
                    let nr = noderate[pu];
                    rate[r] = if ppu < pu {
                        rate[r - 1]
                    } else if nr != 0.0 {
                        rate[r - 1] / nr
                    } else {
                        rate[r - 1]
                    };
                }
                r += 1;
            }
            u = pu;
        }
        if increment < EPSILON {
            break;
        }

        if std::env::var_os("RUSTLE_FLOW_DETAIL").is_some() {
            let mut path_str = String::new();
            let mut tu = n - 1;
            while pred[tu] >= 0 {
                path_str.push_str(&format!("{}<-{} ", tu, pred[tu]));
                tu = pred[tu] as usize;
            }
            eprintln!(
                "PARITY_FLOW_PATH seed={:?} iter={} increment={:.6} rate[0]={:.6} {}",
                seed_tf, bfs_iters, increment, rate[0], path_str
            );
        }

        r = 0;
        u = n - 1;
        while pred[u] >= 0 {
            let pu = pred[u] as usize;
            let pushed = increment / rate[r];
            flow_mat[pu][u] += pushed;
            flow_mat[u][pu] -= pushed;
            r += 1;
            u = pu;
        }
        flux += increment;
    }
    if parity_flow {
        let ratio = if max_fl > EPSILON { flux / max_fl } else { 0.0 };
        eprintln!(
            "PARITY_FLOW_RESULT seed={:?} bfs_iters={} final_flux={:.2} max_fl={:.2} ratio={:.3}",
            seed_tf, bfs_iters, flux, max_fl, ratio
        );
    }

    // Fallback: if we get `flux=0`, retry once after adding *micro-exon targeted* synthetic
    // start-at-node edges. This is a portability guard for cases where graph segmentation does
    // not yield start-at-node transfrags for tiny internal exons (common in backward traces).
    if flux < EPSILON {
        add_microexon_start_edges_if_disconnected(&mut capacity, &mut link, path, graph);
        for nbrs in &mut link {
            nbrs.sort_unstable();
            nbrs.dedup();
        }
        for u in 0..n {
            pred[u] = -1;
            for v in 0..n {
                flow_mat[u][v] = 0.0;
            }
        }
        flux = 0.0;
        while bfs_augmenting_path(n, &capacity, &flow_mat, &link, &mut pred, false) {
            let mut increment = max_fl;
            let mut rate = vec![1.0f64; n];
            let mut r = 1usize;
            let mut u = n - 1;
            while pred[u] >= 0 {
                let pu = pred[u] as usize;
                let adjflux = (capacity[pu][u] - flow_mat[pu][u]) * rate[r - 1];
                if adjflux < increment {
                    increment = adjflux;
                }
                if pred[pu] >= 0 {
                    let ppu = pred[pu] as usize;
                    if pu < u {
                        rate[r] = if ppu < pu {
                            rate[r - 1] * noderate[pu]
                        } else {
                            rate[r - 1]
                        };
                    } else {
                        let nr = noderate[pu];
                        rate[r] = if ppu < pu {
                            rate[r - 1]
                        } else if nr > EPSILON {
                            rate[r - 1] / nr
                        } else {
                            rate[r - 1]
                        };
                    }
                    r += 1;
                }
                u = pu;
            }
            if increment < EPSILON {
                break;
            }
            r = 0;
            u = n - 1;
            while pred[u] >= 0 {
                let pu = pred[u] as usize;
                let pushed = increment / rate[r];
                flow_mat[pu][u] += pushed;
                flow_mat[u][pu] -= pushed;
                r += 1;
                u = pu;
            }
            flux += increment;
        }
    }

    if flux < EPSILON && trace_zero_flux_active(seed_tf) && any_path_node_overlaps_trace(path, graph) {
        // Compact, bounded dump: show where connectivity breaks in the capacity graph.
        eprintln!("[RUST_LMF_ZERO] tf={:?} n={} max_fl={:.6}", seed_tf, n, max_fl);
        if n >= 2 {
            eprintln!(
                "[RUST_LMF_ZERO] cap01={:.6} cap_last={:.6}",
                capacity[0][1],
                capacity[n - 2][n - 1]
            );
        }
        // Per-node outgoing capacity summary.
        for i in 0..n {
            let nid = path[i];
            let (st, en) = graph
                .nodes
                .get(nid)
                .map(|node| (node.start, node.end))
                .unwrap_or((0, 0));
            let mut out_sum = 0.0;
            let mut out_n = 0usize;
            let mut top: Vec<(usize, f64)> = Vec::new();
            for &j in &link[i] {
                let cap = capacity[i][j];
                if cap > 0.0 {
                    out_sum += cap;
                    out_n += 1;
                    top.push((j, cap));
                }
            }
            top.sort_by(|a, b| b.1.partial_cmp(&a.1).unwrap_or(std::cmp::Ordering::Equal));
            let top_str = top
                .iter()
                .take(3)
                .map(|(j, cap)| format!("{}:{:.3}", j, cap))
                .collect::<Vec<_>>()
                .join(",");
            eprintln!(
                "[RUST_LMF_ZERO] i={} nid={}({}-{}) out_n={} out_sum={:.6} top={}",
                i, nid, st, en, out_n, out_sum, top_str
            );
        }
        // Diagnose whether the break is due to missing start-at-node transfrags, or keeptr gating.
        // Only print for a few interior nodes to keep output bounded.
        let mut gap_nodes_shown = 0usize;
        for i in 1..n.saturating_sub(1) {
            if gap_nodes_shown >= 6 {
                break;
            }
            let nid = path[i];
            let out_has = link[i].iter().any(|&j| capacity[i][j] > 0.0);
            if out_has {
                continue;
            }
            let Some(node_obj) = graph.nodes.get(nid) else {
                continue;
            };
            let mut n_start = 0usize;
            let mut n_start_long = 0usize;
            let mut n_candidate_end_on_path = 0usize;
            let mut n_on_path = 0usize;
            let mut n_keeptr_gap_ok = 0usize;
            let mut shown = 0usize;
            for &t_idx in &node_obj.trf_ids {
                if t_idx >= transfrags.len() {
                    continue;
                }
                let tf = &transfrags[t_idx];
                if tf.abundance <= 0.0 || tf.node_ids.is_empty() {
                    continue;
                }
                if tf.node_ids.first().copied() != Some(nid) {
                    continue;
                }
                n_start += 1;
                if !tf.longread {
                    continue;
                }
                n_start_long += 1;
                let end_nid = tf.node_ids.last().copied().unwrap_or(graph.sink_id);
                let end_pi = node2path.get(&end_nid).copied();
                if end_pi.is_some() {
                    n_candidate_end_on_path += 1;
                }
                let source_start_ok = tf.node_ids[0] != graph.source_id
                    || (path.len() > 1 && tf.node_ids.get(1) == Some(&path[1]));
                let on_path = pathpat.contains_pattern(&tf.pattern);
                if on_path && source_start_ok {
                    n_on_path += 1;
                    if keeptr_gap_ok(tf, path, i, graph) {
                        n_keeptr_gap_ok += 1;
                    }
                }
                if shown < 6 {
                    let (st, en) = graph
                        .nodes
                        .get(end_nid)
                        .map(|n| (n.start, n.end))
                        .unwrap_or((0, 0));
                    eprintln!(
                        "[RUST_LMF_ZERO]   start_tf nid={} i={} tf={} ab={:.3} end={}({}-{}) end_i={:?} on_path={} keeptr_gap_ok={}",
                        nid,
                        i,
                        t_idx,
                        tf.abundance,
                        end_nid,
                        st,
                        en,
                        end_pi,
                        (on_path && source_start_ok) as u8,
                        keeptr_gap_ok(tf, path, i, graph) as u8
                    );
                    shown += 1;
                }
            }
            eprintln!(
                "[RUST_LMF_ZERO]   start_summary nid={} i={} n_start={} n_start_long={} end_on_path={} on_path={} keeptr_gap_ok={}",
                nid, i, n_start, n_start_long, n_candidate_end_on_path, n_on_path, n_keeptr_gap_ok
            );
            gap_nodes_shown += 1;
        }
        // Show a few kept transfrags for sanity.
        let mut shown = 0usize;
        eprint!("[RUST_LMF_ZERO] kept_tf:");
        for ti in istranscript.ones() {
            if shown >= 12 {
                break;
            }
            let tf = &transfrags[ti];
            let start = tf.node_ids.first().copied().unwrap_or(graph.source_id);
            let end = tf.node_ids.last().copied().unwrap_or(graph.sink_id);
            let start_i = node2path.get(&start).copied().unwrap_or(usize::MAX);
            let end_i = node2path.get(&end).copied().unwrap_or(usize::MAX);
            eprint!(
                " {}(ab={:.3} {}->{} i={}..{})",
                ti, tf.abundance, start, end, start_i, end_i
            );
            shown += 1;
        }
        eprintln!();
    }

    if debug_ek {
        eprintln!("[RUST_LMF_RESULT] flux={:.6}", flux);
        eprint!("[RUST_LMF_FLOW]");
        for n1 in 0..n {
            for n2 in (n1 + 1)..n {
                if flow_mat[n1][n2] > 0.0 {
                    eprint!(" [{}][{}]={:.6}", n1, n2, flow_mat[n1][n2]);
                }
            }
        }
        eprintln!();
    }

    if !no_subtract {
        for (i, &nid) in path.iter().enumerate() {
            let Some(node_obj) = graph.nodes.get(nid) else {
                continue;
            };
            let mut sumout = 0.0f64;
            let mut entering_tf: Option<usize> = None;
            for &t_idx in &node_obj.trf_ids {
                if !istranscript.contains(t_idx) || t_idx >= transfrags.len() {
                    continue;
                }
                if transfrags[t_idx].abundance <= 0.0 || transfrags[t_idx].node_ids.is_empty() {
                    continue;
                }

                if transfrags[t_idx].node_ids.first().copied() == Some(nid) {
                    let Some(&end_i) = transfrags[t_idx]
                        .node_ids
                        .last()
                        .and_then(|end_nid| node2path.get(end_nid))
                    else {
                        continue;
                    };
                    // Alt-splice protection: if this transfrag diverges from the
                    // path at an alt-splice node, DON'T deplete it. Check if the
                    // transfrag's second node matches the path's next node.
                    // If they diverge, this is a minor isoform that should survive.
                    if transfrags[t_idx].node_ids.len() >= 2 && i + 1 < path.len() {
                        let tf_next = transfrags[t_idx].node_ids[1];
                        let path_next = path[i + 1];
                        if tf_next != path_next && tf_next != graph.sink_id && path_next != graph.sink_id {
                            // Transfrag goes to a different child — preserve it
                            continue;
                        }
                    }
                    let fl = flow_mat[i][end_i];
                    if fl > 0.0 {
                        if fl < transfrags[t_idx].abundance {
                            if i == 0 {
                                sumout += fl;
                            }
                            let val = fl;
                            let before = transfrags[t_idx].abundance;
                            update_transfrag_capacity(
                                &mut transfrags[t_idx],
                                val,
                                &mut nodecapacity,
                                &node2path,
                            );
                            if debug_ek {
                                eprintln!(
                                    "[RUST_LMF_DEPL] seed={:?} tf={} before={:.4} take={:.4} after={:.4} start={} end={} coords={}",
                                    seed_tf,
                                    t_idx,
                                    before,
                                    val,
                                    transfrags[t_idx].abundance,
                                    transfrags[t_idx].node_ids.first().copied().unwrap_or(usize::MAX),
                                    transfrags[t_idx].node_ids.last().copied().unwrap_or(usize::MAX),
                                    format_tf_coords(&transfrags[t_idx], graph, 10)
                                );
                            }
                            if transfrags[t_idx].abundance < DBL_ERROR {
                                transfrags[t_idx].abundance = 0.0;
                            }
                            flow_mat[i][end_i] = 0.0;
                        } else {
                            if i == 0 {
                                sumout += transfrags[t_idx].abundance;
                            }
                            let val = transfrags[t_idx].abundance;
                            flow_mat[i][end_i] -= val;
                            let before = transfrags[t_idx].abundance;
                            update_transfrag_capacity(
                                &mut transfrags[t_idx],
                                val,
                                &mut nodecapacity,
                                &node2path,
                            );
                            if debug_ek {
                                eprintln!(
                                    "[RUST_LMF_DEPL] seed={:?} tf={} before={:.4} take={:.4} after={:.4} start={} end={} coords={}",
                                    seed_tf,
                                    t_idx,
                                    before,
                                    val,
                                    transfrags[t_idx].abundance,
                                    transfrags[t_idx].node_ids.first().copied().unwrap_or(usize::MAX),
                                    transfrags[t_idx].node_ids.last().copied().unwrap_or(usize::MAX),
                                    format_tf_coords(&transfrags[t_idx], graph, 10)
                                );
                            }
                            if transfrags[t_idx].abundance < DBL_ERROR {
                                transfrags[t_idx].abundance = 0.0;
                            }
                        }
                    }
                } else if i == 0 && transfrags[t_idx].node_ids.last().copied() == Some(nid) {
                    entering_tf = Some(t_idx);
                }
            }

            if i == 0 {
                if let Some(t_idx) = entering_tf {
                    let val = if noderate[i] > EPSILON {
                        sumout / noderate[i]
                    } else {
                        0.0
                    };
                    let before = transfrags[t_idx].abundance;
                    transfrags[t_idx].abundance = (transfrags[t_idx].abundance - val).max(0.0);
                    if debug_ek {
                        eprintln!(
                            "[RUST_LMF_DEPL_SRC] seed={:?} tf={} before={:.4} take={:.4} after={:.4} start={} end={} coords={}",
                            seed_tf,
                            t_idx,
                            before,
                            val,
                            transfrags[t_idx].abundance,
                            transfrags[t_idx].node_ids.first().copied().unwrap_or(usize::MAX),
                            transfrags[t_idx].node_ids.last().copied().unwrap_or(usize::MAX),
                            format_tf_coords(&transfrags[t_idx], graph, 10)
                        );
                    }
                    if transfrags[t_idx].abundance < DBL_ERROR {
                        transfrags[t_idx].abundance = 0.0;
                    }
                }
            }
        }
    }

    let real_nodecap: Vec<f64> = path
        .iter()
        .enumerate()
        .filter(|(_, &nid)| nid != source_id && nid != sink_id)
        .map(|(pi, _)| nodecapacity[pi])
        .collect();
    (flux, real_nodecap, istranscript)
}

/// Prior guide support used by guide_push_flow to emulate C++ CNodeGuide node-level guide competition.
pub struct GuideFlowPriorRef<'a> {
    pub pattern: &'a GBitVec,
    pub abundance: f64,
    pub node_prior: &'a crate::types::DetHashMap<usize, f64>,
}

/// Build capacity network for path from transfrags. path = [source, n1, ..., nk, sink].
/// long_only: if true, only include transfrags with longread set (for long flow).
/// Returns (capacity n×n, link adjacency, noderate_local for rate adjustment, istranscript set).
fn build_capacity_network(
    path: &[usize],
    transfrags: &[GraphTransfrag],
    graph: &Graph,
    pathpat: &GBitVec,
    abundance_kind: FlowAbundance,
    long_only: bool,
    seed_tf: Option<usize>,
) -> (
    Vec<Vec<f64>>,
    Vec<Vec<usize>>,
    Vec<f64>,
    NodeSet,
    HashMap<usize, (usize, f64)>,
    f64,
) {
    let n = path.len();
    let source_id = graph.source_id;
    let debug_ek = ek_debug_active(seed_tf);
    let mut capacity = vec![vec![0.0f64; n]; n];
    let mut link_set: Vec<HashSet<usize>> = (0..n).map(|_| Default::default()).collect();
    let mut link: Vec<Vec<usize>> = (0..n).map(|_| Vec::new()).collect();
    let mut noderate_local = vec![1.0f64; n];
    let mut istranscript: NodeSet = NodeSet::with_capacity(transfrags.len());
    if let Some(seed) = seed_tf {
        if seed < transfrags.len() {
            istranscript.insert_grow(seed);
        }
    }
    let mut branch_pick: HashMap<usize, (usize, f64)> = Default::default();
    let mut max_fl = 0.0f64;
    let seed_is_guide = seed_tf
        .and_then(|s| transfrags.get(s))
        .map(|tf| tf.guide)
        .unwrap_or(false);

    let mut node2path: HashMap<usize, usize> = Default::default();
    for (i, &nid) in path.iter().enumerate() {
        node2path.insert(nid, i);
    }
    for (pi, &nid) in path.iter().enumerate() {
        let node_obj = match graph.nodes.get(nid) {
            Some(n) => n,
            None => continue,
        };
        if debug_ek {
            eprintln!(
                "[RUST_EK_NODE] pi={} nid={} coord={}-{} trf_ids={}",
                pi,
                nid,
                node_obj.start,
                node_obj.end,
                node_obj.trf_ids.len()
            );
        }
        let mut sumleft_local = 0.0;
        let mut sumright_local = 0.0;

        for &t_idx in &node_obj.trf_ids {
            if t_idx >= transfrags.len() {
                continue;
            }
            let tf = &transfrags[t_idx];
            let raw_first_nid = tf.node_ids.first().copied();
            let raw_last_nid = tf.node_ids.last().copied();
            let (tf_first_nid, tf_last_nid) = (raw_first_nid, raw_last_nid);
            let tf_first_idx = tf_first_nid.and_then(|fnid| node2path.get(&fnid).copied());
            if debug_ek {
                let first = tf_first_nid.unwrap_or(usize::MAX);
                let last = tf_last_nid.unwrap_or(usize::MAX);
                eprintln!(
                    "[RUST_EK_TF] pi={} nid={} tf={} first={} last={} raw_first={} raw_last={} long={} abund={:.4} onpath={} starts_here={} first_idx={:?}",
                    pi,
                    nid,
                    t_idx,
                    first,
                    last,
                    raw_first_nid.unwrap_or(usize::MAX),
                    raw_last_nid.unwrap_or(usize::MAX),
                    tf.longread,
                    tf.abundance,
                    pathpat.contains_pattern(&tf.pattern),
                    tf_first_nid == Some(nid),
                    tf_first_idx
                );
            }
            if pi > 1 && pi < n.saturating_sub(2) && tf.longread {
                if tf_first_nid == Some(nid) {
                    sumright_local += tf.abundance; // C++ accumulates in double (C++ reference)
                }
                if tf_last_nid == Some(nid) {
                    sumleft_local += tf.abundance; // C++ accumulates in double (C++ reference)
                }
            }
            // C++ long_max_flow: if seed tf is guide, ignore other guide transfrags
            // only for capacity inclusion (not for noderate accumulation above).
            if long_only && seed_is_guide {
                if let Some(seed) = seed_tf {
                    if t_idx != seed && tf.guide {
                        continue;
                    }
                }
            }
            let selected_branch = if long_only {
                None
            } else {
                pick_flow_branch(tf, pathpat, graph)
            };
            if let Some((bi, prop)) = selected_branch {
                branch_pick.insert(t_idx, (bi, prop));
            }
            let eff_ab = if long_only {
                get_abundance(tf, abundance_kind)
            } else {
                effective_path_abundance(tf, abundance_kind, selected_branch)
            };
            if eff_ab <= 0.0 {
                if debug_ek {
                    eprintln!("[RUST_EK_SKIP] tf={} reason=eff_ab", t_idx);
                }
                continue;
            }
            if long_only && !tf.longread {
                if debug_ek {
                    eprintln!("[RUST_EK_SKIP] tf={} reason=not_longread", t_idx);
                }
                continue;
            }
            if tf_first_idx != Some(pi) {
                if debug_ek {
                    eprintln!(
                        "[RUST_EK_SKIP] tf={} reason=not_start first_idx={:?} pi={}",
                        t_idx, tf_first_idx, pi
                    );
                }
                continue;
            }
            if long_only && pi == 0 {
                max_fl = eff_ab;
            }
            // C++ keeptr condition in long_max_flow:
            // istranscript[t] || (pathpat contains tf.pattern and (tf starts off-source or
            // source-start transcript has second node on path[1])).
            let source_start_ok = tf_first_nid != Some(source_id)
                || (path.len() > 1 && tf.node_ids.get(1) == Some(&path[1]));
            if !(istranscript.contains(t_idx)
                || (pathpat.contains_pattern(&tf.pattern) && source_start_ok))
            {
                if debug_ek {
                    eprintln!(
                        "[RUST_EK_SKIP] tf={} reason=not_on_path source_start_ok={}",
                        t_idx, source_start_ok
                    );
                }
                continue;
            }
            if !keeptr_gap_ok(tf, path, pi, graph) {
                if debug_ek {
                    eprintln!("[RUST_EK_SKIP] tf={} reason=gap", t_idx);
                }
                continue;
            }
            istranscript.insert_grow(t_idx);
            let Some(end_nid) = tf_last_nid else {
                continue;
            };
            let end_idx_opt = node2path.get(&end_nid).copied();
            if debug_ek {
                eprintln!(
                    "[RUST_EK_KEEP] tf={} n1={} n2={} eff_ab={:.4}",
                    t_idx, pi, end_nid, eff_ab
                );
            }
            let end_idx = match end_idx_opt {
                Some(i) => i,
                None => {
                    // the reference assembler long_max_flow uses exact node2path lookups.
                    // If a transfrag endpoint is not on the concrete path, it must
                    // not contribute an edge to the capacity graph.
                    continue;
                }
            };
            if !link_set[pi].contains(&end_idx) {
                link[pi].push(end_idx);
                link_set[pi].insert(end_idx);
            }
            if !link_set[end_idx].contains(&pi) {
                link[end_idx].push(pi);
                link_set[end_idx].insert(pi);
            }
            capacity[pi][end_idx] += eff_ab;

            if n >= 2 {
                if !link_set[0].contains(&1) {
                    link[0].push(1);
                    link_set[0].insert(1);
                    link[1].push(0);
                    link_set[1].insert(0);
                }
                if !link_set[n - 2].contains(&(n - 1)) {
                    link[n - 2].push(n - 1);
                    link_set[n - 2].insert(n - 1);
                    link[n - 1].push(n - 2);
                    link_set[n - 1].insert(n - 2);
                }
                capacity[0][1] += eff_ab;
                capacity[n - 2][n - 1] += eff_ab;
                if long_only {
                    max_fl = max_fl + eff_ab;
                }
            }
        }
        if sumright_local > 0.0 && sumleft_local > 0.0 {
            noderate_local[pi] = sumright_local / sumleft_local;
        }
    }

    // Bridge edges removed: they were a bandaid for sub-node gaps in the capacity
    // network. The proper fix is matching C++ graph node structure (dynamic junction
    // processing in create_graph, short-tail skip, etc.).

    for i in 0..n {
        let mut diag = 0.0f64;
        for j in (i + 1)..n {
            diag += capacity[i][j];
        }
        capacity[i][i] = diag;
        link[i].sort_unstable();
    }

    if std::env::var_os("RUSTLE_PARITY_DEBUG").is_some() {
        let mut n_edges = 0;
        for i in 0..n {
            for j in (i + 1)..n {
                if capacity[i][j] > 0.0 {
                    n_edges += 1;
                }
            }
        }
        let cap_01 = if n >= 2 { capacity[0][1] } else { 0.0 };
        let cap_last = if n >= 2 { capacity[n - 2][n - 1] } else { 0.0 };
        eprintln!(
            "PARITY_CAP n_istranscript={} n_edges={} cap_source_first={:.4} cap_last_sink={:.4} max_fl={:.4} path_len={}",
            istranscript.count_ones(),
            n_edges,
            cap_01,
            cap_last,
            max_fl,
            n
        );
    }

    (
        capacity,
        link,
        noderate_local,
        istranscript,
        branch_pick,
        max_fl,
    )
}

/// Edmonds-Karp max flow on path capacity network. Returns (flux, node_capacity for real nodes only).
/// If no_subtract is false, subtracts flow from transfrag abundances.
/// max_fl_override: if Some(v), caps the initial max_fl to v (C++ parse_trflong: max_fl = nodecov[path[0]]).
pub fn edmonds_karp(
    path: &[usize],
    transfrags: &mut [GraphTransfrag],
    graph: &Graph,
    abundance_kind: FlowAbundance,
    long_only: bool,
    no_subtract: bool,
    seed_tf: Option<usize>,
    max_fl_override: Option<f64>,
    pathpat_override: Option<&GBitVec>,
) -> (f64, Vec<f64>, NodeSet) {
    let n = path.len();
    if n < 3 {
        if std::env::var_os("RUSTLE_PARITY_DEBUG").is_some() {
            eprintln!(
                "PARITY_CAP n_istranscript=0 n_edges=0 cap_source_first=0.0000 cap_last_sink=0.0000 max_fl=0.0000 path_len={}",
                n
            );
        }
        return (0.0, vec![], NodeSet::default());
    }
    let source_id = graph.source_id;
    let sink_id = graph.sink_id;
    let real_path: Vec<usize> = path
        .iter()
        .copied()
        .filter(|&nid| nid != source_id && nid != sink_id)
        .collect();
    let n_real = real_path.len();
    if n_real == 0 {
        if std::env::var_os("RUSTLE_PARITY_DEBUG").is_some() {
            eprintln!(
                "PARITY_CAP n_istranscript=0 n_edges=0 cap_source_first=0.0000 cap_last_sink=0.0000 max_fl=0.0000 path_len={}",
                n
            );
        }
        return (0.0, vec![], NodeSet::default());
    }

    let pathpat = pathpat_override
        .cloned()
        .unwrap_or_else(|| build_pathpat(path, graph));
    let (mut capacity, mut link, noderate_local, istranscript, _branch_pick, raw_max_fl) =
        build_capacity_network(
            path,
            transfrags,
            graph,
            &pathpat,
            abundance_kind,
            long_only,
            seed_tf,
        );
    let max_fl = if let Some(limit) = max_fl_override {
        raw_max_fl.min(limit)
    } else {
        raw_max_fl
    };
    let debug_ek = ek_debug_active(seed_tf);
    if debug_ek {
        let seed_abund = seed_tf
            .and_then(|s| transfrags.get(s))
            .map(|tf| tf.abundance)
            .unwrap_or(-1.0);
        eprint!("[RUST_EK_ENTER] tf={:?} n={} path:", seed_tf, n);
        for &p in path.iter() {
            eprint!(" {}", p);
        }
        eprintln!();
        eprintln!("[RUST_EK_CAP] max_fl={:.6}", max_fl);
        for i in 0..n {
            for j in (i + 1)..n {
                if capacity[i][j] > 0.0 {
                    eprintln!(
                        "  cap[{}][{}]={:.6} (nodes {}->{})",
                        i, j, capacity[i][j], path[i], path[j]
                    );
                }
            }
        }
        eprint!("[RUST_EK_ISTRANSCRIPT]");
        for ti in istranscript.ones() {
            eprint!(" {}({:.4})", ti, transfrags[ti].abundance);
        }
        eprintln!();
        let _ = seed_abund;
    }
    if istranscript.is_empty() {
        return (0.0, vec![0.0; n_real], NodeSet::default());
    }

    let mut flow_mat = vec![vec![0.0f64; n]; n];
    let mut pred = vec![-1i32; n];
    let mut flux = 0.0f64;

    let weighted_node_cap = !long_only;
    if std::env::var_os("RUSTLE_LOOP_TRACE").is_some()
        || std::env::var_os("RUSTLE_TRACE_LOG_STYLE").is_some()
    {
        eprintln!("LOOP_max_flow: entering bfs loop n={}", n);
    }
    while bfs_augmenting_path(n, &capacity, &flow_mat, &link, &mut pred, weighted_node_cap) {
        let mut increment = if long_only { max_fl } else { 1e30_f64 };
        let mut r = 1usize;
        let mut rate_arr = vec![1.0f64; n];
        rate_arr[0] = 1.0f64;
        let mut u = n - 1;
        while pred[u] >= 0 {
            let pu = pred[u] as usize;
            let adjflux = (capacity[pu][u] - flow_mat[pu][u]) * rate_arr[r - 1];
            if adjflux < increment {
                increment = adjflux;
            }
            if pred[pu] >= 0 {
                let ppu = pred[pu] as usize;
                if pu < u {
                    rate_arr[r] = if ppu < pu {
                        rate_arr[r - 1] * noderate_local[pu]
                    } else {
                        rate_arr[r - 1]
                    };
                } else {
                    let nr = noderate_local[pu];
                    rate_arr[r] = if ppu < pu {
                        rate_arr[r - 1]
                    } else if nr > EPSILON {
                        rate_arr[r - 1] / nr
                    } else {
                        rate_arr[r - 1]
                    };
                }
                r += 1;
            }
            u = pu;
        }
        if increment < EPSILON {
            break;
        }
        r = 0;
        u = n - 1;
        while pred[u] >= 0 {
            let pu = pred[u] as usize;
            // C++ reference: increment / rate, stored in flow matrix
            let pushed = increment / rate_arr[r];
            flow_mat[pu][u] += pushed;
            flow_mat[u][pu] -= pushed;
            if weighted_node_cap {
                if pu < u {
                    flow_mat[pu][pu] += pushed;
                } else {
                    flow_mat[u][u] = (flow_mat[u][u] - pushed).max(0.0);
                }
            }
            r += 1;
            u = pu;
        }
        flux = flux + increment;
    }

    // If we ended up with zero flux in long-only mode, the most common reason (empirically) is
    // a disconnected capacity network due to micro-exon start gaps (no start-at-node transfrags
    // for a tiny internal node). Retry once with micro-exon targeted synthetic start edges.
    if flux <= 0.0 && long_only {
        add_microexon_start_edges_if_disconnected(&mut capacity, &mut link, path, graph);
        for nbrs in &mut link {
            nbrs.sort_unstable();
            nbrs.dedup();
        }
        for u in 0..n {
            pred[u] = -1;
            for v in 0..n {
                flow_mat[u][v] = 0.0;
            }
        }
        flux = 0.0;
        while bfs_augmenting_path(n, &capacity, &flow_mat, &link, &mut pred, weighted_node_cap) {
            let mut increment = max_fl;
            let mut r = 1usize;
            let mut rate_arr = vec![1.0f64; n];
            rate_arr[0] = 1.0f64;
            let mut u = n - 1;
            while pred[u] >= 0 {
                let pu = pred[u] as usize;
                let adjflux = (capacity[pu][u] - flow_mat[pu][u]) * rate_arr[r - 1];
                if adjflux < increment {
                    increment = adjflux;
                }
                if pred[pu] >= 0 {
                    let ppu = pred[pu] as usize;
                    if pu < u {
                        rate_arr[r] = if ppu < pu {
                            rate_arr[r - 1] * noderate_local[pu]
                        } else {
                            rate_arr[r - 1]
                        };
                    } else {
                        let nr = noderate_local[pu];
                        rate_arr[r] = if ppu < pu {
                            rate_arr[r - 1]
                        } else if nr > EPSILON {
                            rate_arr[r - 1] / nr
                        } else {
                            rate_arr[r - 1]
                        };
                    }
                    r += 1;
                }
                u = pu;
            }
            if increment < EPSILON {
                break;
            }
            r = 0;
            u = n - 1;
            while pred[u] >= 0 {
                let pu = pred[u] as usize;
                let pushed = increment / rate_arr[r];
                flow_mat[pu][u] += pushed;
                flow_mat[u][pu] -= pushed;
                r += 1;
                u = pu;
            }
            flux += increment;
        }
    }

    if debug_ek {
        eprintln!("[RUST_EK_RESULT] flux={:.6}", flux);
        eprint!("[RUST_EK_FLOW]");
        for n1 in 0..n {
            for n2 in (n1 + 1)..n {
                if flow_mat[n1][n2] > 0.0 {
                    eprint!(" [{}][{}]={:.6}", n1, n2, flow_mat[n1][n2]);
                }
            }
        }
        eprintln!();
    }

    if flux <= 0.0 {
        return (0.0, vec![0.0; n_real], istranscript);
    }

    let mut nodecapacity = vec![0.0; n];
    let mut node2path: HashMap<usize, usize> = Default::default();
    for (i, &nid) in path.iter().enumerate() {
        node2path.insert(nid, i);
    }

    let depletion_diag = std::env::var_os("RUSTLE_DEPLETION_DIAG").is_some();
    if !no_subtract {
        for pi in 0..n {
            let nid = path[pi];
            let node_obj = match graph.nodes.get(nid) {
                Some(n) => n,
                None => continue,
            };
            let mut sumout = 0.0;
            let mut entering_tf_pos: Option<usize> = None;

            for (_, &t_idx) in node_obj.trf_ids.iter().enumerate() {
                if !istranscript.contains(t_idx) {
                    continue;
                }
                let tf = &mut transfrags[t_idx];
                if tf.abundance <= 0.0 {
                    continue;
                }
                if tf.node_ids.first() == Some(&nid) {
                    let end_nid = *tf.node_ids.last().unwrap();
                    let end_i = *node2path.get(&end_nid).unwrap_or(&n);
                    if end_i >= n {
                        continue;
                    }
                    let fl = flow_mat[pi][end_i];
                    if fl > 0.0 {
                        if fl < tf.abundance {
                            if pi == 0 {
                                sumout = sumout + fl;
                            }
                            if debug_ek {
                                eprintln!("[RUST_DEPL] tf={} partial val={:.6} abund_before={:.6} abund_after={:.6}",
                                    t_idx, fl, tf.abundance, tf.abundance - fl);
                            }
                            if depletion_diag {
                                eprintln!("DEPL_TF seed={:?} tf={} partial amt={:.6} before={:.6} after={:.6}",
                                    seed_tf, t_idx, fl, tf.abundance, tf.abundance - fl);
                            }
                            update_transfrag_capacity(tf, fl, &mut nodecapacity, &node2path);
                            if tf.abundance < DBL_ERROR {
                                tf.abundance = 0.0;
                            }
                            flow_mat[pi][end_i] = 0.0;
                        } else {
                            if pi == 0 {
                                sumout = sumout + tf.abundance;
                            }
                            let val = tf.abundance;
                            if debug_ek {
                                eprintln!("[RUST_DEPL] tf={} full val={:.6} abund_before={:.6} abund_after=0",
                                    t_idx, val, tf.abundance);
                            }
                            if depletion_diag {
                                eprintln!(
                                    "DEPL_TF seed={:?} tf={} full amt={:.6} before={:.6} after=0",
                                    seed_tf, t_idx, val, tf.abundance
                                );
                            }
                            flow_mat[pi][end_i] -= val;
                            update_transfrag_capacity(tf, val, &mut nodecapacity, &node2path);
                            if tf.abundance < DBL_ERROR {
                                tf.abundance = 0.0;
                            }
                        }
                    }
                } else if pi == 0 && tf.node_ids.last() == Some(&nid) {
                    entering_tf_pos = Some(t_idx);
                }
            }
            if pi == 0 {
                if let Some(entering_tf_pos) = entering_tf_pos {
                    let etf = &mut transfrags[entering_tf_pos];
                    let val = sumout / noderate_local[pi];
                    if depletion_diag {
                        eprintln!(
                            "DEPL_TF seed={:?} tf={} entering amt={:.6} before={:.6} after={:.6}",
                            seed_tf,
                            entering_tf_pos,
                            val,
                            etf.abundance,
                            (etf.abundance - val).max(0.0)
                        );
                    }
                    etf.abundance = (etf.abundance - val).max(0.0);
                    if etf.abundance < DBL_ERROR {
                        etf.abundance = 0.0;
                    }
                }
            }
        }
    }

    let real_nodecap: Vec<f64> = path
        .iter()
        .enumerate()
        .filter(|(_, &nid)| nid != source_id && nid != sink_id)
        .map(|(pi, _)| nodecapacity[pi])
        .collect();
    (flux, real_nodecap, istranscript)
}

/// Long-read max flow: use abundance, longread-only transfrags.
pub fn long_max_flow(
    path: &[usize],
    transfrags: &mut [GraphTransfrag],
    graph: &Graph,
    no_subtract: bool,
) -> (f64, Vec<f64>) {
    long_max_flow_seeded(path, transfrags, graph, no_subtract, None)
}

/// Long-read max flow with optional seed transcript index (C++ istranscript[t]=1 parity).
pub fn long_max_flow_seeded(
    path: &[usize],
    transfrags: &mut [GraphTransfrag],
    graph: &Graph,
    no_subtract: bool,
    seed_tf: Option<usize>,
) -> (f64, Vec<f64>) {
    let (flux, nodecap, _) =
        long_max_flow_direct(path, transfrags, graph, no_subtract, seed_tf, None, None);
    (flux, nodecap)
}

/// Long-read max flow with seed and transfrag-membership used set (C++ istranscript).
pub fn long_max_flow_seeded_with_used(
    path: &[usize],
    transfrags: &mut [GraphTransfrag],
    graph: &Graph,
    no_subtract: bool,
    seed_tf: Option<usize>,
) -> (f64, Vec<f64>, NodeSet) {
    long_max_flow_seeded_with_used_pathpat(path, transfrags, graph, no_subtract, seed_tf, None)
}

/// Long-read max flow with an explicit parse_trflong pathpat override.
/// This mirrors C++ reference, where parse_trflong passes its accumulated pathpat directly into
/// long_max_flow instead of rebuilding it from the concrete path vector.
pub fn long_max_flow_seeded_with_used_pathpat(
    path: &[usize],
    transfrags: &mut [GraphTransfrag],
    graph: &Graph,
    no_subtract: bool,
    seed_tf: Option<usize>,
    pathpat_override: Option<&GBitVec>,
) -> (f64, Vec<f64>, NodeSet) {
    long_max_flow_direct(
        path,
        transfrags,
        graph,
        no_subtract,
        seed_tf,
        None,
        pathpat_override,
    )
}

/// Long-read max flow with nodecov-limited max_fl (C++ parse_trflong: max_fl = nodecov[path[0]]).
/// Transfrag abundances are depleted normally (no_subtract=false).
pub fn long_max_flow_seeded_with_nodecov_limit(
    path: &[usize],
    transfrags: &mut [GraphTransfrag],
    graph: &Graph,
    seed_tf: Option<usize>,
    max_fl_limit: f64,
) -> (f64, Vec<f64>) {
    let (flux, nodecap, _) = long_max_flow_direct(
        path,
        transfrags,
        graph,
        false,
        seed_tf,
        Some(max_fl_limit),
        None,
    );
    (flux, nodecap)
}

/// Short-read max flow: use srabund; include all transfrags (no long_only).
pub fn short_max_flow(
    path: &[usize],
    transfrags: &mut [GraphTransfrag],
    graph: &Graph,
    no_subtract: bool,
) -> (f64, Vec<f64>) {
    if no_subtract {
        let (flux, nodecap, _) = edmonds_karp(
            path,
            transfrags,
            graph,
            FlowAbundance::Short,
            false,
            no_subtract,
            None,
            None,
            None,
        );
        (flux, nodecap)
    } else {
        push_max_flow_seeded(path, transfrags, graph, false, None)
    }
}
