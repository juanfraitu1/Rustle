//! Longtrim direct-style node splitting from read-boundary points (lstart/lend).
//! This mirrors longtrim behavior more closely than simple peak splitting:
//! scan boundary events in genomic order, split when local coverage contrast supports
//! a start/end boundary, and add source/sink + continuity edges.

use crate::bpcov::Bpcov;
use crate::coord::inclusive_to_half_open;
use crate::graph::{Graph, GraphTransfrag};
use crate::read_boundaries::LongtrimBoundaryMap;
use crate::reference_gtf::GuideInfo;
use crate::reference_gtf::RefTranscript;
use crate::trace_events::{
    longtrim_bnode, longtrim_bound_end, longtrim_bound_start, longtrim_split,
};
use crate::types::ReadBoundary;

const CHI_THR: i64 = 50;
const DROP: f64 = 0.5;
const ERROR_PERC: f64 = 0.1;
const TRTHR: f64 = 1.0;
const LONGINTRONANCHOR: u64 = 25;

fn longtrim_trace_active() -> bool {
    std::env::var_os("RUSTLE_TRACE_LONGTRIM").is_some()
}

#[derive(Debug, Default, Clone, Copy)]
pub struct LongtrimStats {
    pub start_boundary_events: usize,
    pub end_boundary_events: usize,
    pub split_nodes: usize,
    pub new_nodes: usize,
    pub source_edges_added: usize,
    pub sink_edges_added: usize,
    pub synthetic_transfrags: usize,
}

#[derive(Debug, Clone, Copy)]
pub struct LongtrimNodeCall {
    pub source_bid: usize,
    pub node_id: usize,
    pub endcov: bool,
    pub post_startcov: Option<bool>,
}

#[derive(Debug, Default, Clone)]
pub struct LongtrimBundleSchedule {
    pub source_bid: usize,
    pub calls: Vec<LongtrimNodeCall>,
}

#[derive(Debug, Default, Clone, Copy)]
pub struct GuideSplitStats {
    pub split_nodes: usize,
    pub new_nodes: usize,
    pub source_edges_added: usize,
    pub sink_edges_added: usize,
    pub synthetic_transfrags: usize,
}

#[inline]
fn cov_range(bpcov: &Bpcov, start_pos: i64, end_pos_incl: i64) -> f64 {
    if end_pos_incl < start_pos {
        return 0.0;
    }
    let s = if start_pos < 0 { 0 } else { start_pos as u64 };
    let e = if end_pos_incl < 0 {
        0
    } else {
        end_pos_incl as u64
    };
    let (hs, he) = inclusive_to_half_open(s, e);
    bpcov.get_cov_range(hs as usize, he as usize)
}

#[inline]
fn clamp_i64(v: i64, lo: i64, hi: i64) -> i64 {
    if v < lo {
        lo
    } else if v > hi {
        hi
    } else {
        v
    }
}

fn split_node_keep_children(graph: &mut Graph, nid: usize, cut: u64) -> Option<usize> {
    if nid >= graph.nodes.len() {
        return None;
    }
    let node = graph.nodes[nid].clone();
    if cut <= node.start || cut >= node.end {
        return None;
    }
    let old_children: Vec<usize> = node.children.ones().collect();
    for &c in &old_children {
        graph.remove_edge(nid, c);
    }
    graph.nodes[nid].end = cut;
    let new_nid = graph.add_node(cut, node.end).node_id;
    graph.nodes[new_nid].source_bnode = node.source_bnode;
    graph.add_edge(nid, new_nid);
    for &c in &old_children {
        graph.add_edge(new_nid, c);
    }
    Some(new_nid)
}

/// Apply longtrim-like splitting for original graph nodes.
/// Returns synthetic transfrags to append before process_transfrags, plus stats.
///
/// Detects read-boundary events from weighted read start/end counts (lstart/lend),
/// validates each event with a 50-bp bpcov contrast window (matching longtrim()
/// re-validation: `tmpcov = (right_sum - left_sum) / (DROP * CHI_THR)`), and splits
/// nodes where coverage contrast confirms a real transcript start or end.
pub fn apply_longtrim_direct(
    graph: &mut Graph,
    boundary_map: &LongtrimBoundaryMap,
    schedule: &[LongtrimBundleSchedule],
    bpcov: &Bpcov,
    min_boundary_cov: f64,
    trace_s: usize,
    trace_bundle_start: u64,
    trace_bundle_end: u64,
) -> (Vec<GraphTransfrag>, LongtrimStats) {
    let mut synthetic: Vec<GraphTransfrag> = Vec::new();
    let mut stats = LongtrimStats::default();

    if boundary_map.is_empty() {
        return (synthetic, stats);
    }

    let source = graph.source_id;
    let sink = graph.sink_id;
    let trace = longtrim_trace_active();
    let paired_scan = std::env::var_os("RUSTLE_PAIRED_BOUND_SCAN").is_some();
    let paired_prune = std::env::var_os("RUSTLE_PAIRED_PRUNE").is_some();

    for bundle in schedule {
        let Some((bundlenode_lstart, bundlenode_lend)) = boundary_map.get(&bundle.source_bid)
        else {
            continue;
        };

        if paired_scan || paired_prune {
            if paired_scan {
                eprintln!(
                    "PAIRED_SCAN_BUNDLE bid={} calls={} lstart_boundaries={} lend_boundaries={}",
                    bundle.source_bid,
                    bundle.calls.len(),
                    bundlenode_lstart.len(),
                    bundlenode_lend.len()
                );
            }
            let mut prune_targets: Vec<usize> = Vec::new();
            for call in &bundle.calls {
                let nid = call.node_id;
                if nid == source || nid == sink || nid >= graph.nodes.len() {
                    continue;
                }
                let n = &graph.nodes[nid];
                let ns = n.start;
                let ne = n.end;
                if ne <= ns.saturating_add(50) {
                    continue;
                }
                let len = ne - ns;
                let s0 = ns.saturating_sub(bpcov.bundle_start) as i64;
                let max_idx = bpcov.cov.len().saturating_sub(1) as i64;
                let s1 = clamp_i64((ne.saturating_sub(1)).saturating_sub(bpcov.bundle_start) as i64, 0, max_idx);
                let total = cov_range(bpcov, clamp_i64(s0, 0, max_idx), s1);
                let avg = if len > 0 { total / (len as f64) } else { 0.0 };
                if avg >= 3.0 {
                    continue;
                }
                // Flanking signature: LEND just upstream of this node and LSTART just downstream.
                // Search a small window (± 100bp) outside the node.
                let win: u64 = 100;
                let up_lo = ns.saturating_sub(win);
                let up_hi = ns; // exclusive upper
                let dn_lo = ne;
                let dn_hi = ne.saturating_add(win);
                let max_lend_up = bundlenode_lend
                    .iter()
                    .filter(|b| b.pos >= up_lo && b.pos < up_hi)
                    .map(|b| b.cov)
                    .fold(0.0f64, f64::max);
                let max_lstart_dn = bundlenode_lstart
                    .iter()
                    .filter(|b| b.pos >= dn_lo && b.pos < dn_hi)
                    .map(|b| b.cov)
                    .fold(0.0f64, f64::max);
                let max_lend_in = bundlenode_lend
                    .iter()
                    .filter(|b| b.pos >= ns && b.pos < ne)
                    .map(|b| b.cov)
                    .fold(0.0f64, f64::max);
                let max_lstart_in = bundlenode_lstart
                    .iter()
                    .filter(|b| b.pos >= ns && b.pos < ne)
                    .map(|b| b.cov)
                    .fold(0.0f64, f64::max);
                let paired_bound_thr: f64 = std::env::var("RUSTLE_PAIRED_BOUND_THR")
                    .ok()
                    .and_then(|v| v.parse().ok())
                    .unwrap_or(10.0);
                if max_lend_up >= paired_bound_thr && max_lstart_dn >= paired_bound_thr {
                    // Flanking avg cov: 100bp upstream and downstream of the node.
                    let up_s = clamp_i64(
                        (ns.saturating_sub(100)).saturating_sub(bpcov.bundle_start) as i64,
                        0, max_idx);
                    let up_e = clamp_i64(
                        (ns.saturating_sub(1)).saturating_sub(bpcov.bundle_start) as i64,
                        0, max_idx);
                    let dn_s = clamp_i64(
                        ne.saturating_sub(bpcov.bundle_start) as i64, 0, max_idx);
                    let dn_e = clamp_i64(
                        (ne.saturating_add(99)).saturating_sub(bpcov.bundle_start) as i64,
                        0, max_idx);
                    let up_cov = cov_range(bpcov, up_s, up_e) / ((up_e - up_s + 1).max(1) as f64);
                    let dn_cov = cov_range(bpcov, dn_s, dn_e) / ((dn_e - dn_s + 1).max(1) as f64);
                    let ratio_up = if avg > 0.0 { up_cov / avg } else { up_cov };
                    let ratio_dn = if avg > 0.0 { dn_cov / avg } else { dn_cov };
                    let strict = ratio_up >= 10.0 && ratio_dn >= 10.0
                        && up_cov >= 10.0 && dn_cov >= 10.0;
                    if paired_scan {
                        eprintln!(
                            "PAIRED_BOUND bid={} nid={} node={}..{} len={} avg_cov={:.2} up_lend={:.2} dn_lstart={:.2} up_cov={:.2} dn_cov={:.2} ratio_up={:.1} ratio_dn={:.1} strict={}",
                            bundle.source_bid, nid, ns, ne, len, avg,
                            max_lend_up, max_lstart_dn, up_cov, dn_cov, ratio_up, ratio_dn,
                            if strict { 1 } else { 0 }
                        );
                    }
                    if paired_prune && strict {
                        prune_targets.push(nid);
                    }
                }
            }
            if paired_prune {
                for &nid in &prune_targets {
                    if nid >= graph.nodes.len() { continue; }
                    let parents: Vec<usize> = graph.nodes[nid].parents.ones().collect();
                    let children: Vec<usize> = graph.nodes[nid].children.ones().collect();
                    let ns = graph.nodes[nid].start;
                    let ne = graph.nodes[nid].end;
                    let mut removed_skip = 0usize;
                    for &p in &parents { graph.remove_edge(p, nid); }
                    for &c in &children { graph.remove_edge(nid, c); }
                    for &p in &parents {
                        if p == source { continue; }
                        let p_kids: Vec<usize> = graph.nodes[p].children.ones().collect();
                        for &c in &p_kids {
                            if c == nid || c == sink { continue; }
                            if graph.nodes[c].start >= ne {
                                graph.remove_edge(p, c);
                                removed_skip += 1;
                            }
                        }
                    }
                    // Stitch source/sink so both sides become complete assembly subgraphs.
                    if std::env::var_os("RUSTLE_PAIRED_PRUNE_STITCH").is_some() {
                        for &p in &parents {
                            if p != source && p != sink {
                                graph.add_edge(p, sink);
                            }
                        }
                        for &c in &children {
                            if c != source && c != sink {
                                graph.add_edge(source, c);
                                if let Some(n) = graph.nodes.get_mut(c) {
                                    n.hardstart = true;
                                }
                            }
                        }
                        for &p in &parents {
                            if let Some(n) = graph.nodes.get_mut(p) {
                                n.hardend = true;
                            }
                        }
                    }
                    eprintln!(
                        "PAIRED_PRUNE_FIRE bid={} nid={} node={}..{} parents={} children={} skip_removed={}",
                        bundle.source_bid, nid, ns, ne, parents.len(), children.len(), removed_skip
                    );
                }
            }
        }

        let mut nls = 0usize;
        let mut nle = 0usize;
        let mut startcov = false;

        for call in &bundle.calls {
            let nid = call.node_id;
            if nid == source || nid == sink || nid >= graph.nodes.len() {
                continue;
            }

            let node0 = graph.nodes[nid].clone();
            if node0.end <= node0.start.saturating_add(1) {
                if let Some(next_state) = call.post_startcov {
                    startcov = next_state;
                }
                continue;
            }
            let node_start0 = node0.start;
            let node_end0 = node0.end;

            while nls < bundlenode_lstart.len() && bundlenode_lstart[nls].pos < node_start0 {
                nls += 1;
            }
            while nle < bundlenode_lend.len() && bundlenode_lend[nle].pos < node_start0 {
                nle += 1;
            }

            let mut s_end = nls;
            while s_end < bundlenode_lstart.len() && bundlenode_lstart[s_end].pos < node_end0 {
                s_end += 1;
            }
            let mut e_end = nle;
            while e_end < bundlenode_lend.len() && bundlenode_lend[e_end].pos < node_end0 {
                e_end += 1;
            }

            let trace_starts: Vec<ReadBoundary> = bundlenode_lstart[nls..s_end]
                .iter()
                .copied()
                .filter(|b| {
                    b.cov < 0.0
                        || (b.cov >= min_boundary_cov && b.pos > node_start0 && b.pos < node_end0)
                })
                .collect();
            let trace_ends: Vec<ReadBoundary> = bundlenode_lend[nle..e_end]
                .iter()
                .copied()
                .filter(|b| {
                    b.cov < 0.0
                        || (b.cov >= min_boundary_cov
                            && b.pos >= node_start0
                            && b.pos < node_end0.saturating_sub(1))
                })
                .collect();

            if trace_starts.is_empty() && trace_ends.is_empty() {
                if let Some(next_state) = call.post_startcov {
                    startcov = next_state;
                }
                continue;
            }

            stats.start_boundary_events += trace_starts.len();
            stats.end_boundary_events += trace_ends.len();

            if trace {
                longtrim_bnode(
                    trace_s,
                    trace_bundle_start,
                    trace_bundle_end,
                    node_start0,
                    node_end0,
                    bpcov.bundle_start,
                );
                for b in &trace_starts {
                    longtrim_bound_start(trace_s, node_start0, node_end0, b.pos as i32, b.cov);
                }
                for b in &trace_ends {
                    longtrim_bound_end(trace_s, node_start0, node_end0, b.pos as i32, b.cov);
                }
            }

            let old_children: Vec<usize> = node0.children.ones().collect();
            for &c in &old_children {
                graph.remove_edge(nid, c);
            }

            let mut cur_nid = nid;
            let mut did_split_node = false;

            while ((nls < bundlenode_lstart.len() && bundlenode_lstart[nls].pos < node_end0)
                || (nle < bundlenode_lend.len() && bundlenode_lend[nle].pos < node_end0))
                && graph.nodes[cur_nid].start < node_end0
            {
                let next_is_start = if nle >= bundlenode_lend.len()
                    || bundlenode_lend[nle].pos >= node_end0
                {
                    true
                } else if nls >= bundlenode_lstart.len() || bundlenode_lstart[nls].pos >= node_end0
                {
                    false
                } else {
                    bundlenode_lstart[nls].pos <= bundlenode_lend[nle].pos
                };

                if next_is_start {
                    let b = bundlenode_lstart[nls];
                    nls += 1;
                    let pos = b.pos;
                    if !(b.cov < 0.0
                        || (b.cov >= min_boundary_cov && pos > node_start0 && pos < node_end0))
                    {
                        continue;
                    }
                    if pos <= graph.nodes[cur_nid].start || pos >= node_end0 {
                        continue;
                    }

                    // note:
                    // LONGINTRONANCHOR spacing avoids creating tiny nodes from noisy boundary spikes.
                    // However, in real loci we can have a low-support early end and a strong true end
                    // within <25 bp; the original algorithm still needs to represent the strong end.
                    // Treat high-support boundary events as "strong" and allow them to bypass spacing.
                    let strong_boundary = b.cov >= (min_boundary_cov.max(1.0) * 4.0);
                    let allow = ((startcov
                        || pos > graph.nodes[cur_nid].start.saturating_add(LONGINTRONANCHOR))
                        || strong_boundary)
                        && (call.endcov || pos < node_end0.saturating_add(LONGINTRONANCHOR));
                    if !allow {
                        continue;
                    }

                    let startpos = (pos.saturating_sub(bpcov.bundle_start)) as i64;
                    let max_idx = bpcov.cov.len().saturating_sub(1) as i64;
                    let winstart = clamp_i64(startpos - CHI_THR, 0, max_idx);
                    let winend = clamp_i64(startpos + CHI_THR - 1, 0, max_idx);
                    let right = cov_range(bpcov, startpos, winend);
                    let left = cov_range(bpcov, winstart, startpos - 1);
                    let mut tmpcov = (right - left) / (DROP * CHI_THR as f64);
                    if tmpcov <= 0.0 && b.cov < 0.0 {
                        tmpcov = ERROR_PERC;
                    }
                    if tmpcov <= 0.0 {
                        continue;
                    }

                    tmpcov += TRTHR;

                    let prev_end = graph.nodes[cur_nid].end;
                    if pos <= graph.nodes[cur_nid].start || pos >= prev_end {
                        continue;
                    }
                    graph.nodes[cur_nid].end = pos;
                    let new_nid = graph.add_node(pos, prev_end).node_id;
                    graph.nodes[new_nid].source_bnode = graph.nodes[cur_nid].source_bnode;
                    graph.nodes[new_nid].hardstart = true;
                    graph.add_edge(source, new_nid);
                    graph.add_edge(cur_nid, new_nid);

                    let mut tf_source =
                        GraphTransfrag::new(vec![source, new_nid], graph.pattern_size());
                    tf_source.abundance = tmpcov;
                    tf_source.longread = true;
                    tf_source.longstart = pos;
                    tf_source.longend = pos + 1;
                    synthetic.push(tf_source);

                    let mut tf_link =
                        GraphTransfrag::new(vec![cur_nid, new_nid], graph.pattern_size());
                    tf_link.abundance = TRTHR;
                    tf_link.longread = true;
                    synthetic.push(tf_link);

                    if trace {
                        longtrim_split(
                            trace_s,
                            "start",
                            pos as i32,
                            tmpcov,
                            node_start0,
                            node_end0,
                        );
                    }
                    did_split_node = true;
                    stats.new_nodes += 1;
                    stats.source_edges_added += 1;
                    stats.synthetic_transfrags += 2;
                    startcov = false;
                    cur_nid = new_nid;
                } else {
                    let b = bundlenode_lend[nle];
                    nle += 1;
                    let pos = b.pos;
                    if !(b.cov < 0.0
                        || (b.cov >= min_boundary_cov
                            && pos >= node_start0
                            && pos < node_end0.saturating_sub(1)))
                    {
                        continue;
                    }
                    if pos < graph.nodes[cur_nid].start || pos >= node_end0 {
                        continue;
                    }

                    let strong_boundary = b.cov >= (min_boundary_cov.max(1.0) * 4.0);
                    let allow = ((!startcov
                        || pos > graph.nodes[cur_nid].start.saturating_add(LONGINTRONANCHOR))
                        || strong_boundary)
                        && (!call.endcov || pos < node_end0.saturating_add(LONGINTRONANCHOR));
                    if !allow {
                        continue;
                    }

                    let endpos = (pos.saturating_sub(bpcov.bundle_start)) as i64;
                    let max_idx = bpcov.cov.len().saturating_sub(1) as i64;
                    let winstart = clamp_i64(endpos - CHI_THR + 1, 0, max_idx);
                    let winend = clamp_i64(endpos + CHI_THR, 0, max_idx);
                    let left = cov_range(bpcov, winstart, endpos);
                    let right = cov_range(bpcov, endpos + 1, winend);
                    let mut tmpcov = (left - right) / (DROP * CHI_THR as f64);

                    if let Ok(v) = std::env::var("RUSTLE_TRACE_LONGTRIM_POS") {
                        if let Ok(tp) = v.parse::<u64>() {
                            // Accept both "pos" in 0-based and "cut" in 1-based-ish (pos+1) for convenience.
                            if tp == pos || tp == pos.saturating_add(1) {
                                eprintln!(
                                    "[LONGTRIM_DBG] end_boundary pos={} cut={} node={}({}-{}) left={:.3} right={:.3} tmpcov={:.6} b_cov={:.3}",
                                    pos,
                                    pos.saturating_add(1),
                                    cur_nid,
                                    graph.nodes[cur_nid].start,
                                    node_end0,
                                    left,
                                    right,
                                    tmpcov,
                                    b.cov
                                );
                            }
                        }
                    }

                    if tmpcov <= 0.0 && b.cov < 0.0 {
                        tmpcov = ERROR_PERC;
                    }
                    if tmpcov <= 0.0 {
                        continue;
                    }

                    tmpcov += TRTHR;

                    let prev_end = graph.nodes[cur_nid].end;
                    let cut = pos.saturating_add(1);
                    if cut <= graph.nodes[cur_nid].start || cut >= prev_end {
                        continue;
                    }
                    graph.nodes[cur_nid].end = cut;
                    graph.nodes[cur_nid].hardend = true;
                    let new_nid = graph.add_node(cut, prev_end).node_id;
                    graph.nodes[new_nid].source_bnode = graph.nodes[cur_nid].source_bnode;
                    graph.add_edge(cur_nid, new_nid);
                    graph.add_edge(cur_nid, sink);

                    let mut tf_sink =
                        GraphTransfrag::new(vec![cur_nid, sink], graph.pattern_size());
                    tf_sink.abundance = tmpcov;
                    tf_sink.longread = true;
                    tf_sink.longstart = graph.nodes[cur_nid].start;
                    tf_sink.longend = cut;
                    synthetic.push(tf_sink);

                    let mut tf_link =
                        GraphTransfrag::new(vec![cur_nid, new_nid], graph.pattern_size());
                    tf_link.abundance = TRTHR;
                    tf_link.longread = true;
                    synthetic.push(tf_link);

                    if trace {
                        longtrim_split(trace_s, "end", pos as i32, tmpcov, node_start0, node_end0);
                    }
                    did_split_node = true;
                    stats.new_nodes += 1;
                    stats.sink_edges_added += 1;
                    stats.synthetic_transfrags += 2;
                    startcov = true;
                    cur_nid = new_nid;
                }
            }

            for &c in &old_children {
                graph.add_edge(cur_nid, c);
            }

            if did_split_node {
                stats.split_nodes += 1;
            }

            if let Some(next_state) = call.post_startcov {
                startcov = next_state;
            }
        }
    }

    if stats.split_nodes > 0 {
        graph.compute_reachability();
    }

    (synthetic, stats)
}

fn source2guide_abundance(bpcov: &Bpcov, node_start: u64, newstart: u64, newend: u64) -> f64 {
    const CHI_WIN: u64 = 100;
    let mut leftcov = 0.0;
    let mut rightcov = 0.0;
    if newstart > node_start {
        let gstart = newstart.saturating_sub(CHI_WIN).max(node_start);
        let len = newstart.saturating_sub(gstart).max(1);
        leftcov = bpcov.get_cov_range(bpcov.idx(gstart), bpcov.idx(newstart)) / len as f64;
    }
    if newstart < newend {
        let gend_excl = (newstart + CHI_WIN).min(newend);
        let len = gend_excl.saturating_sub(newstart).max(1);
        rightcov = bpcov.get_cov_range(bpcov.idx(newstart), bpcov.idx(gend_excl)) / len as f64;
    }
    (rightcov - leftcov).max(TRTHR)
}

fn guide2sink_abundance(bpcov: &Bpcov, node_start: u64, split_excl: u64, node_end: u64) -> f64 {
    const CHI_WIN: u64 = 100;
    let mut leftcov = 0.0;
    let mut rightcov = 0.0;
    if split_excl >= node_start {
        let gstart = split_excl.saturating_sub(CHI_WIN).max(node_start);
        let len = split_excl.saturating_sub(gstart).max(1);
        leftcov = bpcov.get_cov_range(bpcov.idx(gstart), bpcov.idx(split_excl)) / len as f64;
    }
    if split_excl < node_end {
        let gend_excl = (split_excl + CHI_WIN).min(node_end);
        let len = gend_excl.saturating_sub(split_excl).max(1);
        rightcov = bpcov.get_cov_range(bpcov.idx(split_excl), bpcov.idx(gend_excl)) / len as f64;
    }
    (leftcov - rightcov).max(TRTHR)
}

/// Guide-boundary splitting akin to source2guide/guide2sink in create_graph flow.
/// Splits first/last guide nodes at guide transcript boundaries and creates source/sink
/// support transfrags so guide boundaries participate in transfrag processing.
pub fn apply_guide_boundary_splits(
    graph: &mut Graph,
    guides: &[GuideInfo],
    bpcov: &Bpcov,
) -> (Vec<GraphTransfrag>, GuideSplitStats) {
    let mut synth = Vec::new();
    let mut st = GuideSplitStats::default();
    let source = graph.source_id;
    let sink = graph.sink_id;

    for g in guides {
        if g.node_ids.is_empty() {
            continue;
        }
        let first = g.node_ids[0];
        let last = *g.node_ids.last().unwrap_or(&first);
        if first == source || first == sink || last == source || last == sink {
            continue;
        }

        // source2guide-like split at guide tx_start.
        if first < graph.nodes.len() {
            let n = graph.nodes[first].clone();
            if g.tx_start > n.start && g.tx_start < n.end {
                if let Some(new_first) = split_node_keep_children(graph, first, g.tx_start) {
                    let had = graph.nodes[source].children.contains(new_first);
                    graph.add_edge(source, new_first);
                    graph.nodes[new_first].hardstart = true;
                    if !had {
                        st.source_edges_added += 1;
                    }
                    let mut tf = GraphTransfrag::new(vec![source, new_first], graph.pattern_size());
                    tf.abundance = source2guide_abundance(bpcov, n.start, g.tx_start, n.end);
                    tf.longread = true;
                    tf.guide = true;
                    tf.longstart = g.tx_start;
                    tf.longend = g.tx_start.saturating_add(1);
                    synth.push(tf);
                    st.split_nodes += 1;
                    st.new_nodes += 1;
                    st.synthetic_transfrags += 1;
                }
            }
        }

        // guide2sink-like split at guide tx_end (half-open end).
        if last < graph.nodes.len() {
            let n = graph.nodes[last].clone();
            if g.tx_end > n.start && g.tx_end < n.end {
                if split_node_keep_children(graph, last, g.tx_end).is_some() {
                    let had = graph.nodes[last].children.contains(sink);
                    graph.add_edge(last, sink);
                    graph.nodes[last].hardend = true;
                    if !had {
                        st.sink_edges_added += 1;
                    }
                    let mut tf = GraphTransfrag::new(vec![last, sink], graph.pattern_size());
                    tf.abundance = guide2sink_abundance(bpcov, n.start, g.tx_end, n.end);
                    tf.longread = true;
                    tf.guide = true;
                    tf.longstart = n.start;
                    tf.longend = g.tx_end;
                    synth.push(tf);
                    st.split_nodes += 1;
                    st.new_nodes += 1;
                    st.synthetic_transfrags += 1;
                }
            } else if g.tx_end == n.end {
                let had = graph.nodes[last].children.contains(sink);
                graph.add_edge(last, sink);
                graph.nodes[last].hardend = true;
                if !had {
                    st.sink_edges_added += 1;
                }
            }
        }
    }

    if st.split_nodes > 0 {
        graph.compute_reachability();
    }
    (synth, st)
}

/// Apply guide boundary splits directly from reference transcript boundaries, prior to guide mapping.
/// This is closer to create_graph-time source2guide/guide2sink handling.
pub fn apply_guide_boundary_splits_from_refs(
    graph: &mut Graph,
    guides: &[RefTranscript],
    bundle_start: u64,
    bundle_end: u64,
    strand: char,
    bpcov: &Bpcov,
    min_boundary_cov: f64,
) -> (Vec<GraphTransfrag>, GuideSplitStats) {
    let mut synth = Vec::new();
    let mut st = GuideSplitStats::default();
    let source = graph.source_id;
    let sink = graph.sink_id;

    for g in guides {
        if g.chrom.is_empty() || g.exons.is_empty() {
            continue;
        }
        if strand != '.' && g.strand != '.' && g.strand != strand {
            continue;
        }
        let tx_start = g.exons.first().map(|e| e.0).unwrap_or(0);
        let tx_end = g.exons.last().map(|e| e.1).unwrap_or(0);
        if tx_end <= bundle_start || tx_start >= bundle_end {
            continue;
        }

        let mut first_nid: Option<usize> = None;
        let mut last_nid: Option<usize> = None;
        for (nid, n) in graph.nodes.iter().enumerate() {
            if nid == source || nid == sink {
                continue;
            }
            if tx_start >= n.start && tx_start < n.end {
                first_nid = Some(nid);
            }
            let tx_end_last = tx_end.saturating_sub(1);
            if tx_end_last >= n.start && tx_end_last < n.end {
                last_nid = Some(nid);
            }
        }

        if let Some(first) = first_nid {
            let n = graph.nodes[first].clone();
            if tx_start > n.start.saturating_add(LONGINTRONANCHOR)
                && tx_start < n.end.saturating_sub(LONGINTRONANCHOR)
            {
                if let Some(new_first) = split_node_keep_children(graph, first, tx_start) {
                    let ab = source2guide_abundance(bpcov, n.start, tx_start, n.end);
                    if ab < min_boundary_cov {
                        continue;
                    }
                    let had = graph.nodes[source].children.contains(new_first);
                    graph.add_edge(source, new_first);
                    graph.nodes[new_first].hardstart = true;
                    if !had {
                        st.source_edges_added += 1;
                    }
                    let mut tf = GraphTransfrag::new(vec![source, new_first], graph.pattern_size());
                    tf.abundance = ab;
                    tf.longread = true;
                    tf.guide = true;
                    tf.longstart = tx_start;
                    tf.longend = tx_start.saturating_add(1);
                    synth.push(tf);
                    st.split_nodes += 1;
                    st.new_nodes += 1;
                    st.synthetic_transfrags += 1;
                }
            }
        }

        if let Some(last) = last_nid {
            let n = graph.nodes[last].clone();
            if tx_end > n.start.saturating_add(LONGINTRONANCHOR)
                && tx_end < n.end.saturating_sub(LONGINTRONANCHOR)
            {
                if split_node_keep_children(graph, last, tx_end).is_some() {
                    let ab = guide2sink_abundance(bpcov, n.start, tx_end, n.end);
                    if ab < min_boundary_cov {
                        continue;
                    }
                    let had = graph.nodes[last].children.contains(sink);
                    graph.add_edge(last, sink);
                    graph.nodes[last].hardend = true;
                    if !had {
                        st.sink_edges_added += 1;
                    }
                    let mut tf = GraphTransfrag::new(vec![last, sink], graph.pattern_size());
                    tf.abundance = ab;
                    tf.longread = true;
                    tf.guide = true;
                    tf.longstart = n.start;
                    tf.longend = tx_end;
                    synth.push(tf);
                    st.split_nodes += 1;
                    st.new_nodes += 1;
                    st.synthetic_transfrags += 1;
                }
            } else if tx_end == n.end {
                let had = graph.nodes[last].children.contains(sink);
                graph.add_edge(last, sink);
                graph.nodes[last].hardend = true;
                if !had {
                    st.sink_edges_added += 1;
                }
            }
        }
    }

    if st.split_nodes > 0 {
        graph.compute_reachability();
    }
    (synth, st)
}
