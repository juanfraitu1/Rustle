//! Global flow extraction layer — sole long-read extractor using proven edmonds_karp from max_flow.rs.
//!
//! This module implements the flow pipeline:
//! 1. `compute_global_nodecov`: get_trf_long nodecov precomputation
//! 2. `extract_transcripts_global_flow`: parse_trflong seed iteration
//!    with edmonds_karp flow + retry logic matching path_extract.rs

use crate::bitvec::GBitVec;
use crate::constants::FLOW_EPSILON;
use crate::coord::len_half_open;
use crate::graph::{Graph, GraphTransfrag};
use crate::max_flow::{long_max_flow_seeded_with_used, long_max_flow_seeded_with_nodecov_limit};
use crate::path_extract::{
    LongRecDiag, Transcript, back_to_source_fast_long, best_trf_match,
    fwd_to_sink_fast_long,
    redistribute_transfrag_to_matches,
};
use crate::types::{DetHashSet as HashSet, RunConfig};

const EPSILON: f64 = crate::constants::FLOW_EPSILON;
const TRTHR: f64 = 1.0;
const ERROR_PERC: f64 = 0.1;
const POLY_TAIL_STOP_COUNT: u16 = 8;
const CHI_WIN_ABUNDANCE_THR: f64 = 100.0;
const CHECKTRF_REDISTRIBUTE_INTRON_TOL: u64 = 5;

/// Debug flag for global flow tracing.
fn global_flow_debug() -> bool {
    std::env::var_os("RUSTLE_GLOBAL_FLOW_DEBUG").is_some()
}

// ─── Helper functions ───────────────────────────────────────────────────────

/// Set/clear an edge bit in pathpat.
fn edge_set(pathpat: &mut GBitVec, graph: &Graph, a: usize, b: usize, val: bool) {
    if let Some(eid) = graph.edge_bit_index(a, b) {
        if val {
            pathpat.set_bit(eid);
        } else {
            pathpat.clear_bit(eid);
        }
    }
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
        edge_set(&mut pathpat, graph, a, b, true);
    }
    pathpat
}

/// Are two consecutive nodes contiguous (no intron gap)?
fn nodes_are_contiguous(graph: &Graph, a: usize, b: usize) -> bool {
    if let (Some(na), Some(nb)) = (graph.nodes.get(a), graph.nodes.get(b)) {
        na.end == nb.start || na.end + 1 == nb.start
    } else {
        false
    }
}

/// Is there a splice (non-contiguous gap) between two path nodes?
fn is_splice_between(graph: &Graph, u: usize, v: usize) -> bool {
    if let (Some(nu), Some(nv)) = (graph.nodes.get(u), graph.nodes.get(v)) {
        nv.start > nu.end + 1
    } else {
        false
    }
}

/// Check if any transfrag contains both splice edges (la->ra, lb->rb) as consecutive node pairs in order.
/// Matches path_extract.rs has_lr_witness_two_splices_nodes — does NOT filter on abundance or longread,
/// since witnessing is about whether the splice was observed in the data, not remaining flow.
fn has_lr_witness_two_splices(
    la: usize,
    ra: usize,
    lb: usize,
    rb: usize,
    transfrags: &[GraphTransfrag],
) -> bool {
    for tf in transfrags {
        let ns = &tf.node_ids;
        let mut first_edge_found = false;
        for i in 1..ns.len() {
            if !first_edge_found {
                if ns[i - 1] == la && ns[i] == ra {
                    first_edge_found = true;
                }
            } else if ns[i - 1] == lb && ns[i] == rb {
                return true;
            }
        }
    }
    false
}

fn intron_chain_from_nodes(graph: &Graph, nodes: &[usize]) -> Vec<(u64, u64)> {
    let mut out: Vec<(u64, u64)> = Vec::new();
    for w in nodes.windows(2) {
        let (a, b) = (w[0], w[1]);
        let Some(na) = graph.nodes.get(a) else { continue };
        let Some(nb) = graph.nodes.get(b) else { continue };
        if nb.start > na.end + 1 {
            out.push((na.end, nb.start));
        }
    }
    out
}

fn intron_chains_equal_tol(a: &[(u64, u64)], b: &[(u64, u64)], tol: u64) -> bool {
    if a.len() != b.len() {
        return false;
    }
    a.iter()
        .zip(b.iter())
        .all(|((d1, a1), (d2, a2))| d1.abs_diff(*d2) <= tol && a1.abs_diff(*a2) <= tol)
}

// ─── Step 1: Flow engine (using proven edmonds_karp from max_flow.rs) ───────


// ─── Step 2: get_trf_long nodecov precomputation ─

/// Compute per-node coverage from transfrags: nodecov[i] = max(abundin, abundout).
/// Also computes noderate[i] = node.coverage / nodecov[i].
///
/// This matches get_trf_long exactly:
/// - Skip guide transfrags with abundance < trthr*ERROR_PERC+epsilon
/// - Skip source-to-sink transfrags unless single-exon gene
/// - abundin = sum(tf.abundance where tf.nodes[0] < i)
/// - abundout = sum(tf.abundance where tf.nodes.Last() > i)
pub fn compute_global_nodecov(
    graph: &Graph,
    transfrags: &[GraphTransfrag],
) -> (Vec<f64>, Vec<f64>) {
    let gno = graph.n_nodes;
    let source_id = graph.source_id;
    let sink_id = graph.sink_id;
    let mut nodecov = vec![0.0f64; gno];
    let mut noderate = vec![1.0f64; gno];

    let single_exon_poly_start = if gno == 3 {
        graph
            .nodes
            .get(1)
            .map(|n| n.poly_start_unaligned_total)
            .unwrap_or(0) as f64
    } else {
        0.0
    };
    let single_exon_poly_end = if gno == 3 {
        graph
            .nodes
            .get(1)
            .map(|n| n.poly_end_unaligned_total)
            .unwrap_or(0) as f64
    } else {
        0.0
    };

    for i in 1..gno.saturating_sub(1) {
        let Some(inode) = graph.nodes.get(i) else {
            continue;
        };
        let mut abundin = 0.0f64;
        let mut abundout = 0.0f64;

        for &t in &inode.trf_ids {
            let Some(tf) = transfrags.get(t) else {
                continue;
            };
            if tf.node_ids.is_empty() || tf.abundance <= 0.0 {
                continue;
            }

            // filter: skip guides with low abundance
            if tf.guide && tf.abundance < TRTHR * ERROR_PERC + 1e-9 {
                continue;
            }

            let first_node = tf.node_ids[0];
            let last_node = *tf.node_ids.last().unwrap();

            // filter: skip source-to-sink unless single-exon gene
            let single_exon_gene = gno == 3
                && (single_exon_poly_start > 10.0
                    || single_exon_poly_end > 10.0
                    || tf.abundance >= 50.0);

            if !single_exon_gene && (first_node == source_id && last_node == sink_id) {
                continue;
            }
            // also requires: first_node != source AND last_node != sink
            // (for non-single-exon-gene transfrags)
            if !single_exon_gene && (first_node == source_id || last_node == sink_id) {
                continue;
            }

            // uses node ID comparison (not coordinate comparison)
            if first_node < i {
                abundin += tf.abundance;
            }
            if last_node > i {
                abundout += tf.abundance;
            }
        }

        nodecov[i] = if abundin > abundout {
            abundin
        } else {
            abundout
        };

        let mut rate = nodecov[i];
        if rate <= 0.0 {
            rate = 1.0;
        }
        // inode->cov is per-base (total/len), so noderate = per_base_cov / nodecov
        // Rust: inode.coverage is total, so divide by length to match the original
        rate = inode.coverage / rate / inode.length() as f64;
        noderate[i] = rate;
    }

    (nodecov, noderate)
}

// ─── Step 3: parse_trflong ─────────────────────

/// Extract transcripts using the global flow pipeline.
///
/// This is a faithful port of `parse_trflong`
/// 1. Iterate trflong seeds in reverse order (highest abundance first)
/// 2. For each seed: construct path via back_to_source + fwd_to_sink
/// 3. Run edmonds_karp flow + retry to get flux and nodeflux
/// 4. Build transcript from exons with coverage = sum(nodecapacity[j] * noderate[path[j]]) / len
/// 5. Deplete nodecov by nodecapacity
pub fn extract_transcripts_global_flow(
    graph: &Graph,
    transfrags: &mut [GraphTransfrag],
    bundle_chrom: &str,
    bundle_strand: char,
    config: &RunConfig,
) -> Vec<Transcript> {
    let debug = global_flow_debug();
    let debug_detail = std::env::var_os("DEBUG_DEBUG").is_some();
    let source_id = graph.source_id;
    let sink_id = graph.sink_id;

    // Step 2: compute global nodecov + noderate
    let (mut nodecov, noderate) = compute_global_nodecov(graph, transfrags);

    // Build seed order: same as parse_trflong — iterate trflong in reverse
    let mut seeded: Vec<usize> = transfrags
        .iter()
        .enumerate()
        .filter(|(_, tf)| {
            tf.trflong_seed && tf.weak == 0 && !tf.node_ids.is_empty() && tf.usepath >= 0
        })
        .map(|(i, _)| i)
        .collect();
    seeded.sort_by(|&i, &j| {
        transfrags[i]
            .usepath
            .cmp(&transfrags[j].usepath)
            .then_with(|| transfrags[j].real.cmp(&transfrags[i].real))
            .then_with(|| transfrags[j].guide.cmp(&transfrags[i].guide))
            .then_with(|| {
                transfrags[j]
                    .abundance
                    .partial_cmp(&transfrags[i].abundance)
                    .unwrap_or(std::cmp::Ordering::Equal)
            })
    });

    if debug {
        eprintln!(
            "[GLOBAL_FLOW] extract: gno={} seeds={} strand={}",
            graph.n_nodes,
            seeded.len(),
            bundle_strand
        );
    }

    let mut out: Vec<Transcript> = Vec::new();
    let mut kept_paths: Vec<(Vec<usize>, f64, bool, usize)> = Vec::new(); // (inner nodes, support, guide, out_idx)
    let mut checktrf: Vec<usize> = Vec::new();

    for &idx in &seeded {
        if transfrags[idx].node_ids.is_empty() || transfrags[idx].weak != 0 {
            continue;
        }

        let real_nodes: Vec<usize> = transfrags[idx]
            .node_ids
            .iter()
            .copied()
            .filter(|&n| n != source_id && n != sink_id)
            .collect();
        if real_nodes.is_empty() {
            continue;
        }

        let longcov = transfrags[idx].abundance;

        if debug {
            let coord_start = real_nodes.first().and_then(|&n| graph.nodes.get(n)).map(|n| n.start).unwrap_or(0);
            let coord_end = real_nodes.last().and_then(|&n| graph.nodes.get(n)).map(|n| n.end).unwrap_or(0);
            eprintln!("[GLOBAL_FLOW] seed idx={} abund={:.2} nodes={} span={}-{} usepath={} guide={}",
                idx, longcov, real_nodes.len(), coord_start, coord_end,
                transfrags[idx].usepath, transfrags[idx].guide);
        }

        // Construct pathpat from seed transfrag pattern
        let mut pathpat = transfrags[idx].pattern.clone();
        let mut minp = real_nodes[0];
        let mut maxp = *real_nodes.last().unwrap();

        // add edge bits for consecutive node pairs in seed
        for w in transfrags[idx].node_ids.windows(2) {
            edge_set(&mut pathpat, graph, w[0], w[1], true);
        }
        // set source→minp and maxp→sink edge bits
        edge_set(&mut pathpat, graph, source_id, minp, true);
        edge_set(&mut pathpat, graph, maxp, sink_id, true);

        let maxi = minp;
        let mut path: Vec<usize> = vec![maxi];
        pathpat.set_bit(maxi);

        let mut diag = LongRecDiag::default();

        // back_to_source — use path_extract's proven implementation
        let mut visited_back: HashSet<usize> = Default::default();
        let mut weak_cache: Vec<Option<bool>> = vec![None; transfrags.len()];
        let back_ok = back_to_source_fast_long(
            idx,
            maxi,
            &mut path,
            &mut minp,
            &mut maxp,
            &mut pathpat,
            transfrags,
            graph,
            &nodecov,
            true, // require_longread
            false, // mixed_mode
            &mut diag,
            &mut visited_back,
            &mut weak_cache,
        );
        if !back_ok {
            if debug {
                eprintln!(
                    "[GLOBAL_FLOW] seed {} back_to_source FAILED → checktrf (unreachable={} noreach={} nochoice={} exclude={})",
                    idx, diag.back_unreachable_minpath, diag.back_no_reach, diag.back_no_choice, diag.back_exclude_no_support
                );
            }
            checktrf.push(idx);
            continue;
        }
        path.push(source_id);
        path.reverse();

        // fwd_to_sink — use path_extract's proven implementation
        let mut visited_fwd: HashSet<usize> = Default::default();
        let fwd_ok = fwd_to_sink_fast_long(
            idx,
            maxi,
            &mut path,
            &mut minp,
            &mut maxp,
            &mut pathpat,
            transfrags,
            graph,
            &nodecov,
            true, // require_longread
            false, // mixed_mode
            &mut diag,
            &mut visited_fwd,
            &mut weak_cache,
        );
        if !fwd_ok {
            if debug {
                eprintln!(
                    "[GLOBAL_FLOW] seed {} fwd_to_sink FAILED → checktrf (unreachable={} noreach={} nochoice={} exclude={})",
                    idx, diag.fwd_unreachable_maxpath, diag.fwd_no_reach, diag.fwd_no_choice, diag.fwd_exclude_no_support
                );
            }
            checktrf.push(idx);
            continue;
        }

        if debug {
            let path_coords: Vec<String> = path.iter().filter(|&&n| n != source_id && n != sink_id)
                .map(|&n| graph.nodes.get(n).map(|nd| format!("{}-{}", nd.start, nd.end)).unwrap_or_else(|| format!("?{}", n)))
                .collect();
            eprintln!("[GLOBAL_FLOW] seed {} path_len={} nodes=[{}]",
                idx, path.len(), path_coords.join(", "));
        }

        // Validate path has source at start and sink at end
        if path.first() != Some(&source_id) || path.last() != Some(&sink_id) {
            if debug {
                eprintln!("[GLOBAL_FLOW] seed {} INVALID path (no source/sink) → checktrf", idx);
            }
            checktrf.push(idx);
            continue;
        }
        let startnode = 1usize;
        let lastnode = path.len().saturating_sub(2);
        if lastnode < startnode || lastnode >= path.len() {
            if debug {
                eprintln!("[GLOBAL_FLOW] seed {} INVALID path (startnode/lastnode) → checktrf", idx);
            }
            checktrf.push(idx);
            continue;
        }

        // Poly-tail path trimming (10197-10275)
        let poly_start = transfrags[idx].poly_start_unaligned;
        let poly_end = transfrags[idx].poly_end_unaligned;
        let mut effective_start = startnode;
        let mut effective_last = lastnode;

        if bundle_strand == '+' && poly_end >= POLY_TAIL_STOP_COUNT {
            let tf_last = *real_nodes.last().unwrap();
            if path[effective_last] > tf_last {
                if let Some(cutpos) =
                    (effective_start..=effective_last).find(|&pi| path[pi] == tf_last)
                {
                    effective_last = cutpos;
                }
            }
        }
        if bundle_strand == '-' && poly_start >= POLY_TAIL_STOP_COUNT {
            let tf_first = real_nodes[0];
            if path[effective_start] < tf_first {
                if let Some(keep_from) =
                    (effective_start..=effective_last).find(|&pi| path[pi] == tf_first)
                {
                    effective_start = keep_from;
                }
            }
        }

        // Hardstart/hardend enforcement
        let tf_first_node = real_nodes[0];
        let tf_last_node = *real_nodes.last().unwrap();
        let thardstart = graph
            .nodes
            .get(tf_first_node)
            .map(|n| n.hardstart)
            .unwrap_or(false);
        let thardend = graph
            .nodes
            .get(tf_last_node)
            .map(|n| n.hardend)
            .unwrap_or(false);

        let mut checkpath = true;
        let (use_path, use_start, use_last) = if thardstart
            && thardend
            && (path[effective_start] != tf_first_node || path[effective_last] != tf_last_node)
        {
            if transfrags[idx].abundance > CHI_WIN_ABUNDANCE_THR {
                let mut newpath = vec![source_id];
                newpath.extend_from_slice(&real_nodes);
                newpath.push(sink_id);
                // Rebuild pathpat
                pathpat = build_pathpat(&newpath, graph);
                pathpat.or_assign(&transfrags[idx].pattern);
                checkpath = false; // Using transfrag's own nodes — no witness needed
                (newpath, 1, real_nodes.len())
            } else {
                checktrf.push(idx);
                continue; // deferred to checktrf
            }
        } else {
            (path, effective_start, effective_last)
        };

        // Splice witness check (skip for paths using transfrag's own nodes)
        if checkpath && use_last >= use_start {
            let mut splice_pos: Vec<usize> = Vec::new();
            for p in use_start..=use_last {
                if p == use_start {
                    continue;
                }
                let u = use_path[p - 1];
                let v = use_path[p];
                if is_splice_between(graph, u, v) {
                    splice_pos.push(p);
                }
            }
            let mut unwitnessed = false;
            for k in 1..splice_pos.len() {
                let pa = splice_pos[k - 1];
                let pb = splice_pos[k];
                let la = use_path[pa - 1];
                let ra = use_path[pa];
                let lb = use_path[pb - 1];
                let rb = use_path[pb];
                if !has_lr_witness_two_splices(la, ra, lb, rb, transfrags) {
                    unwitnessed = true;
                    break;
                }
            }
            if unwitnessed {
                if debug {
                    eprintln!("[GLOBAL_FLOW] seed {} UNWITNESSED splice → checktrf", idx);
                }
                checktrf.push(idx);
                continue;
            }
        }

        // DEBUG_DEBUG: seed entry + path
        if debug_detail {
            let minp = use_path.get(use_start).copied().unwrap_or(0);
            let maxp = use_path.get(use_last).copied().unwrap_or(0);
            eprintln!("DEBUG_SEED idx={} maxi={} minp={} maxp={} guide={} abund={}",
                idx, minp, minp, maxp, if transfrags[idx].guide { 1 } else { 0 }, longcov);
            eprint!("DEBUG_SEED_NODES");
            for &p in use_path.iter() { eprint!(" {}", p); }
            eprintln!();
        }

        // Run edmonds_karp flow with retry logic (matching path_extract.rs:3599-3679)
        let first_real_nid = use_path[use_start..=use_last]
            .iter()
            .find(|&&n| n != source_id && n != sink_id)
            .copied();
        let max_fl_limit = first_real_nid
            .and_then(|n| nodecov.get(n).copied())
            .unwrap_or(0.0);

        let (orig_flux, orig_nodeflux, _flow_used) = long_max_flow_seeded_with_used(
            &use_path,
            transfrags,
            graph,
            false,
            Some(idx),
        );

        let (flux, nodeflux) = if orig_flux >= config.readthr {
            (orig_flux, orig_nodeflux)
        } else if max_fl_limit > EPSILON {
            long_max_flow_seeded_with_nodecov_limit(
                &use_path,
                transfrags,
                graph,
                Some(idx),
                max_fl_limit,
            )
        } else {
            (orig_flux, orig_nodeflux)
        };

        // DEBUG_DEBUG: flow result
        if debug_detail {
            eprintln!("DEBUG_FLOW idx={} flux={}", idx, flux);
            eprint!("DEBUG_NODEFLUX");
            for v in &nodeflux { eprint!(" {}", v); }
            eprintln!();
            eprint!("DEBUG_NODECOV");
            for &p in use_path.iter() {
                eprint!(" {}", nodecov.get(p).copied().unwrap_or(0.0));
            }
            eprintln!();
        }

        if flux <= 0.0 && !transfrags[idx].guide {
            if debug_detail {
                eprintln!("DEBUG_SKIP idx={} reason=zero_flux", idx);
            }
            if debug {
                eprintln!("[GLOBAL_FLOW] seed {} flux=0 → checktrf", idx);
            }
            checktrf.push(idx);
            continue;
        }

        // Build exons from path[use_start..=use_last], compute coverage using nodeflux*noderate
        // nodeflux[j] capped to nodecov[path[j]], then
        // ecov = nodeflux[j] * noderate[path[j]]; nodecov[path[j]] -= nodeflux[j]
        let mut exons: Vec<(u64, u64)> = Vec::new();
        let mut cov_total = 0.0f64;
        let mut exon_cov_accum: Vec<f64> = Vec::new();

        // First pass: build exons from merged contiguous nodes
        let mut exon_ranges: Vec<(usize, usize)> = Vec::new(); // (path_start, path_end) for each exon
        {
            let mut ej = use_start;
            while ej <= use_last {
                let nid = use_path[ej];
                if nid >= graph.nodes.len() || nid == source_id || nid == sink_id {
                    ej += 1;
                    continue;
                }
                let node = &graph.nodes[nid];
                let nodestart = node.start;
                let mut nodeend = node.end;
                let exon_path_start = ej;
                while ej + 1 <= use_last {
                    let next_nid = use_path[ej + 1];
                    if next_nid < graph.nodes.len()
                        && nodes_are_contiguous(graph, use_path[ej], next_nid)
                    {
                        ej += 1;
                        nodeend = graph.nodes[next_nid].end;
                    } else {
                        break;
                    }
                }
                if nodeend > nodestart {
                    exons.push((nodestart, nodeend));
                    exon_ranges.push((exon_path_start, ej));
                }
                ej += 1;
            }
        }

        if exons.is_empty() {
            checktrf.push(idx);
            continue;
        }

        // Second pass: compute coverage from nodeflux * noderate, deplete nodecov
        exon_cov_accum.resize(exons.len(), 0.0);
        let mut k = 0usize;
        for (ei, &(ep_start, ep_end)) in exon_ranges.iter().enumerate() {
            for pp in ep_start..=ep_end {
                let nid = use_path[pp];
                if nid == source_id || nid == sink_id || nid >= nodecov.len() {
                    continue;
                }
                if k >= nodeflux.len() {
                    break;
                }
                let raw_nflux = nodeflux[k];
                // Cap to nodecov, then subtract
                let nflux = raw_nflux.min(nodecov[nid].max(0.0));
                nodecov[nid] = (nodecov[nid] - nflux).max(0.0);
                if nodecov[nid] < EPSILON {
                    nodecov[nid] = 0.0;
                }
                let ecov = nflux * noderate[nid];
                cov_total += ecov;
                exon_cov_accum[ei] += ecov;
                k += 1;
            }
        }

        // Normalize per-exon coverage by exon length
        let mut exon_cov: Vec<f64> = Vec::new();
        for (ei, &(es, ee)) in exons.iter().enumerate() {
            let elen = len_half_open(es, ee) as f64;
            if elen > 0.0 {
                exon_cov.push(exon_cov_accum[ei] / elen);
            } else {
                exon_cov.push(0.0);
            }
        }

        let length: u64 = exons.iter().map(|(s, e)| e - s).sum();
        // parse_trflong stores cov directly without dividing by length
        let coverage = if cov_total > 0.0 {
            cov_total
        } else if flux > 0.0 {
            flux
        } else {
            longcov
        };

        if debug {
            eprintln!(
                "[GLOBAL_FLOW] seed {} flux={:.4} cov={:.4} exons={} span={}-{}",
                idx,
                flux,
                coverage,
                exons.len(),
                exons.first().map(|e| e.0).unwrap_or(0),
                exons.last().map(|e| e.1).unwrap_or(0)
            );
        }

        // DEBUG_DEBUG: store prediction
        if debug_detail {
            eprint!("DEBUG_STORE idx={} cov={} len={} nexon={} exons=",
                idx, coverage, length, exons.len());
            for (ei, &(s, e)) in exons.iter().enumerate() {
                if ei > 0 { eprint!(","); }
                eprint!("{}-{}", s, e);
            }
            eprintln!();
        }

        // Track inner nodes for checktrf matching
        let inner_nodes: Vec<usize> = use_path[use_start..=use_last]
            .iter()
            .copied()
            .filter(|&n| n != source_id && n != sink_id)
            .collect();
        let out_idx = out.len();

        let first_node = use_path[use_start];
        let last_node = use_path[use_last];
        out.push(Transcript {
            chrom: bundle_chrom.to_string(),
            strand: bundle_strand,
            exons,
            coverage,
            exon_cov,
            tpm: 0.0,
            fpkm: 0.0,
            source: Some("global_flow".to_string()),
            is_longread: true,
            longcov,
            bpcov_cov: 0.0,
            transcript_id: None,
            gene_id: None,
            ref_transcript_id: None,
            ref_gene_id: None,
            hardstart: graph.nodes.get(first_node).map(|n| n.hardstart).unwrap_or(false),
            hardend: graph.nodes.get(last_node).map(|n| n.hardend).unwrap_or(false),
        });

        kept_paths.push((inner_nodes, coverage, transfrags[idx].guide, out_idx));
    }

    // ─── checktrf rescue pass ───
    // For failed seeds: try to redistribute to best-matched kept paths,
    // or create independent predictions from complete transfrags.
    if !checktrf.is_empty() {
        // Dedup while preserving insertion order.
        let mut seen: HashSet<usize> = Default::default();
        let mut uniq_check = Vec::new();
        for t in checktrf {
            if seen.insert(t) {
                uniq_check.push(t);
            }
        }
        let has_guide_kept = kept_paths.iter().any(|(_, _, g, _)| *g);
        if debug {
            eprintln!("[GLOBAL_FLOW_CHECKTRF] count={} kept_paths={}", uniq_check.len(), kept_paths.len());
        }
        // Compute local nodecov/noderate for independent rescue coverage
        let mut local_nodecov = nodecov.clone();
        let local_noderate = noderate.clone();
        for t in uniq_check {
            if t >= transfrags.len() || transfrags[t].node_ids.is_empty() {
                continue;
            }
            if !(transfrags[t].guide || transfrags[t].abundance >= config.readthr) {
                continue;
            }
            let tf_nodes: Vec<usize> = transfrags[t]
                .node_ids
                .iter()
                .copied()
                .filter(|&n| n != source_id && n != sink_id)
                .collect();
            if tf_nodes.is_empty() {
                continue;
            }

            // For multi-node longread: first try redistribution to best-matched kept paths
            if transfrags[t].longread && tf_nodes.len() > 1 {
                let mut tmatch: Vec<usize> = Vec::new();
                let mut abundancesum = 0.0;
                if has_guide_kept {
                    let (m, s) = best_trf_match(&tf_nodes, &kept_paths, graph, true);
                    tmatch = m;
                    abundancesum = s;
                }
                if tmatch.is_empty() {
                    let (m, s) = best_trf_match(&tf_nodes, &kept_paths, graph, false);
                    tmatch = m;
                    abundancesum = s;
                }
                if !tmatch.is_empty() {
                    // P4 isoform-aware check: skip redistribution when intron counts differ
                    let tf_intron_count = tf_nodes.windows(2).filter(|w| {
                        graph.nodes.get(w[0]).zip(graph.nodes.get(w[1]))
                            .map_or(false, |(na, nb)| nb.start > na.end)
                    }).count();
                    let tf_chain = intron_chain_from_nodes(graph, &tf_nodes);
                    let has_same_intron_count = tmatch.iter().any(|&kidx| {
                        let k_intron_count = kept_paths[kidx].0.windows(2).filter(|w| {
                            graph.nodes.get(w[0]).zip(graph.nodes.get(w[1]))
                                .map_or(false, |(na, nb)| nb.start > na.end)
                        }).count();
                        k_intron_count == tf_intron_count
                    });
                    let has_same_intron_chain = tmatch.iter().any(|&kidx| {
                        let kp_chain = intron_chain_from_nodes(graph, &kept_paths[kidx].0);
                        intron_chains_equal_tol(
                            &tf_chain,
                            &kp_chain,
                            CHECKTRF_REDISTRIBUTE_INTRON_TOL,
                        )
                    });
                    if !has_same_intron_count && tf_intron_count >= 2 {
                        if debug {
                            eprintln!("[GLOBAL_FLOW_CHECKTRF] t={} DIFFERENT_ISOFORM introns={} → independent rescue",
                                t, tf_intron_count);
                        }
                        // Fall through to independent rescue below
                    } else if !has_same_intron_chain && !tf_chain.is_empty() {
                        if debug {
                            eprintln!(
                                "[GLOBAL_FLOW_CHECKTRF] t={} DIFFERENT_INTRON_CHAIN (tol={}) → independent rescue",
                                t, CHECKTRF_REDISTRIBUTE_INTRON_TOL
                            );
                        }
                    } else {
                        if debug {
                            eprintln!("[GLOBAL_FLOW_CHECKTRF] t={} MATCHED kept_paths={:?} abundsum={:.4}",
                                t, &tmatch, abundancesum);
                        }
                        let updated = redistribute_transfrag_to_matches(
                            &transfrags[t],
                            &tmatch,
                            abundancesum,
                            &kept_paths,
                            graph,
                            &mut out,
                        );
                        if updated {
                            transfrags[t].abundance = 0.0;
                        }
                        continue; // matched: done
                    }
                }
            }

            // Independent rescue: store complete transfrag as its own prediction.
            // `else if(!eonly || guide)` — in non-eonly mode always rescues.
            let mut complete = true;
            for w in tf_nodes.windows(2) {
                let a = w[0];
                let b = w[1];
                let contiguous = nodes_are_contiguous(graph, a, b);
                if contiguous {
                    continue;
                }
                let has_edge = graph
                    .nodes
                    .get(a)
                    .map(|n| n.children.contains(b))
                    .unwrap_or(false);
                let has_pattern_edge = graph
                    .edge_bit_index(a, b)
                    .map(|eid| transfrags[t].pattern.get_bit(eid))
                    .unwrap_or(false);
                if !has_edge || !has_pattern_edge {
                    complete = false;
                    break;
                }
            }
            if !complete {
                transfrags[t].abundance = 0.0;
                continue;
            }
            // Build exons from tf_nodes, compute coverage using local_nodecov
            let mut exons: Vec<(u64, u64)> = Vec::new();
            let mut exoncov: Vec<f64> = Vec::new();
            let mut cov_bp_total = 0.0f64;
            let mut j = 0usize;
            while j < tf_nodes.len() {
                let nid = tf_nodes[j];
                let Some(node) = graph.nodes.get(nid) else {
                    j += 1;
                    continue;
                };
                let start = node.start;
                let mut end = node.end;
                let mut exon_bp_cov = 0.0f64;
                if nid < local_nodecov.len() {
                    let nlen = node.length() as f64;
                    if local_nodecov[nid] > EPSILON && nlen > 0.0 {
                        let mut addcov = transfrags[t].abundance * nlen;
                        let rate = local_noderate
                            .get(nid)
                            .copied()
                            .filter(|r| *r > 0.0)
                            .unwrap_or(1.0);
                        let mut newnodecov = local_nodecov[nid] - addcov / rate;
                        if newnodecov < 0.0 {
                            addcov = local_nodecov[nid] * rate;
                            newnodecov = 0.0;
                        }
                        local_nodecov[nid] = newnodecov;
                        exon_bp_cov += addcov;
                    }
                }
                while j + 1 < tf_nodes.len()
                    && nodes_are_contiguous(graph, tf_nodes[j], tf_nodes[j + 1])
                {
                    j += 1;
                    if let Some(nn) = graph.nodes.get(tf_nodes[j]) {
                        end = nn.end;
                        let nid2 = tf_nodes[j];
                        if nid2 < local_nodecov.len() {
                            let nlen = nn.length() as f64;
                            if local_nodecov[nid2] > EPSILON && nlen > 0.0 {
                                let mut addcov = transfrags[t].abundance * nlen;
                                let rate = local_noderate
                                    .get(nid2)
                                    .copied()
                                    .filter(|r| *r > 0.0)
                                    .unwrap_or(1.0);
                                let mut newnodecov = local_nodecov[nid2] - addcov / rate;
                                if newnodecov < 0.0 {
                                    addcov = local_nodecov[nid2] * rate;
                                    newnodecov = 0.0;
                                }
                                local_nodecov[nid2] = newnodecov;
                                exon_bp_cov += addcov;
                            }
                        }
                    }
                }
                if end > start {
                    exons.push((start, end));
                    let elen = len_half_open(start, end).max(1) as f64;
                    exoncov.push(exon_bp_cov / elen);
                    cov_bp_total += exon_bp_cov;
                }
                j += 1;
            }
            if exons.is_empty() {
                transfrags[t].abundance = 0.0;
                continue;
            }
            // Apply longstart/longend trimming
            if transfrags[t].longstart > 0 {
                let (s, e) = exons[0];
                exons[0].0 = s.max(transfrags[t].longstart).min(e);
            }
            if transfrags[t].longend > 0 {
                let last = exons.len() - 1;
                let (s, e) = exons[last];
                exons[last].1 = e.min(transfrags[t].longend).max(s);
            }
            let length: u64 = exons.iter().map(|(s, e)| len_half_open(*s, *e)).sum();
            if length < config.min_transcript_length {
                transfrags[t].abundance = 0.0;
                continue;
            }
            let coverage = if length > 0 {
                cov_bp_total / (length as f64)
            } else {
                0.0
            };
            if cov_bp_total <= EPSILON {
                transfrags[t].abundance = 0.0;
                continue;
            }
            if debug {
                eprintln!("[GLOBAL_FLOW_CHECKTRF] t={} INDEPENDENT rescue exons={} cov={:.4} span={}-{}",
                    t, exons.len(), coverage,
                    exons.first().map(|e| e.0).unwrap_or(0),
                    exons.last().map(|e| e.1).unwrap_or(0));
            }

            let inner: Vec<usize> = tf_nodes.clone();
            let oi = out.len();
            let first_node = tf_nodes[0];
            let last_node = *tf_nodes.last().unwrap_or(&first_node);
            out.push(Transcript {
                chrom: bundle_chrom.to_string(),
                strand: bundle_strand,
                exons: exons.clone(),
                coverage,
                exon_cov: if exoncov.len() == exons.len() && exoncov.iter().any(|v| *v > 0.0) {
                    exoncov
                } else {
                    vec![coverage; exons.len()]
                },
                tpm: 0.0,
                fpkm: 0.0,
                source: transfrags[t]
                    .guide_tid
                    .as_ref()
                    .map(|gid| format!("guide:{gid}")),
                is_longread: true,
                longcov: transfrags[t].abundance,
                bpcov_cov: 0.0,
                transcript_id: None,
                gene_id: None,
                ref_transcript_id: None,
                ref_gene_id: None,
                hardstart: graph.nodes.get(first_node).map(|n| n.hardstart).unwrap_or(false),
                hardend: graph.nodes.get(last_node).map(|n| n.hardend).unwrap_or(false),
            });
            kept_paths.push((inner, coverage, transfrags[t].guide, oi));
            transfrags[t].abundance = 0.0;
        }
    }

    // Post-pass: redistribute remaining long-read transfrags to best kept predictions
    // (second-pass redistribution)
    if !out.is_empty() && !kept_paths.is_empty() {
        let has_guide_kept = kept_paths.iter().any(|(_, _, g, _)| *g);
        for t in 0..transfrags.len() {
            if !transfrags[t].longread || transfrags[t].abundance <= EPSILON {
                continue;
            }
            if transfrags[t].node_ids.len() <= 1 {
                continue;
            }
            let first_raw = transfrags[t].node_ids.first().copied().unwrap_or(source_id);
            let last_raw = transfrags[t].node_ids.last().copied().unwrap_or(sink_id);
            if first_raw == source_id || last_raw == sink_id {
                continue;
            }
            let tf_nodes: Vec<usize> = transfrags[t]
                .node_ids
                .iter()
                .copied()
                .filter(|&n| n != source_id && n != sink_id)
                .collect();
            if tf_nodes.len() <= 1 {
                continue;
            }
            let mut tmatch: Vec<usize> = Vec::new();
            let mut abundancesum = 0.0;
            if has_guide_kept {
                let (m, s) = best_trf_match(&tf_nodes, &kept_paths, graph, true);
                tmatch = m;
                abundancesum = s;
            }
            if tmatch.is_empty() {
                let (m, s) = best_trf_match(&tf_nodes, &kept_paths, graph, false);
                tmatch = m;
                abundancesum = s;
            }
            if tmatch.is_empty() {
                continue;
            }
            let updated = redistribute_transfrag_to_matches(
                &transfrags[t],
                &tmatch,
                abundancesum,
                &kept_paths,
                graph,
                &mut out,
            );
            if updated {
                transfrags[t].abundance = 0.0;
            }
        }
    }

    // after ALL seeds processed, zero ALL longread transfrag abundances
    for tf in transfrags.iter_mut() {
        if tf.longread {
            tf.abundance = 0.0;
        }
    }

    if debug {
        eprintln!(
            "[GLOBAL_FLOW] extract done: {} transcripts (checktrf rescued some)",
            out.len()
        );
    }

    out
}
