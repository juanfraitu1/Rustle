//! Per-node coverage from transfrags (compute_nodecov).

use crate::graph::{Graph, GraphTransfrag};

const TRTHR: f64 = 1.0;
const ERROR_PERC: f64 = 0.1;
const CHI_WIN: f64 = 100.0;
/// Guide-transfrag abundance threshold: trthr * ERROR_PERC + epsilon (get_trf_long:11086).
const GUIDE_ABUND_THR: f64 = TRTHR * ERROR_PERC + 1e-9;

/// Compute per-node coverage: for each node, max(abundin, abundout).
/// - abundin: sum of abundance of transfrags that pass through and have first_node < nid.
/// - abundout: sum of abundance of transfrags that pass through and have last_node > nid.
/// Excludes source/sink-only transfrags unless single-exon gene (n_nodes==3 with poly evidence).
///
/// eliminate_transfrags_under_thr removes fully-internal non-guide transfrags with
/// abundance < trthr=1.0 before inode->trf[] is built. Guide transfrags are kept but filtered in
/// get_trf_long with abundance >= trthr*ERROR_PERC=0.1+eps. Source/sink-anchored transfrags are
/// not removed by eliminate_transfrags_under_thr but are filtered out in get_trf_long.
pub fn compute_nodecov(
    graph: &mut Graph,
    transfrags: &[GraphTransfrag],
    verbose: bool,
) -> Vec<f64> {
    let n_nodes = graph.n_nodes;
    let source_id = graph.source_id;
    let sink_id = graph.sink_id;

    let mut abundin = vec![0.0; n_nodes];
    let mut abundout = vec![0.0; n_nodes];
    let mut longcov = vec![0.0; n_nodes];

    let single_exon_poly_start = if n_nodes == 3 {
        graph
            .nodes
            .get(1)
            .map(|n| n.poly_start_unaligned_total as f64)
            .unwrap_or(0.0)
    } else {
        0.0
    };
    let single_exon_poly_end = if n_nodes == 3 {
        graph
            .nodes
            .get(1)
            .map(|n| n.poly_end_unaligned_total as f64)
            .unwrap_or(0.0)
    } else {
        0.0
    };

    for tf in transfrags {
        if tf.abundance <= 0.0 {
            continue;
        }
        let nodes = &tf.node_ids;
        if nodes.is_empty() {
            continue;
        }
        let first_node = nodes[0];
        let last_node = *nodes.last().unwrap();
        let full_length = first_node == source_id && last_node == sink_id;
        if full_length {
            let single_exon_gene = n_nodes == 3
                && (single_exon_poly_start > 10.0
                    || single_exon_poly_end > 10.0
                    || tf.abundance >= CHI_WIN / 2.0);
            if !single_exon_gene {
                continue;
            }
            // Single-exon gene source-to-sink: no abundance threshold (doesn't eliminate these)
        } else if first_node == source_id || last_node == sink_id {
            continue; // source/sink-anchored: filtered in get_trf_long, not in eliminate_transfrags_under_thr
        } else {
            // Fully internal transfrag: apply elimination threshold
            // Non-guides need abundance >= trthr=1.0; guides need >= trthr*ERROR_PERC+eps=0.1+eps
            let threshold = if tf.guide || tf.guide_tid.is_some() {
                GUIDE_ABUND_THR
            } else {
                TRTHR
            };
            if tf.abundance < threshold {
                continue;
            }
        }

        let real_nodes: Vec<usize> = nodes
            .iter()
            .copied()
            .filter(|&n| n != source_id && n != sink_id)
            .collect();
        if real_nodes.is_empty() {
            continue;
        }

        // Use genomic coordinates for entering/exiting transfrag classification.
        // StringTie uses node IDs (which are coordinate-ordered in C++), but after
        // longtrim, Rustle's node IDs are NOT coordinate-ordered. Using coordinates
        // ensures correct nodecov computation regardless of ID ordering.
        let first_coord = graph.nodes.get(real_nodes[0]).map(|n| n.start).unwrap_or(0);
        let last_coord = graph.nodes.get(*real_nodes.last().unwrap()).map(|n| n.start).unwrap_or(0);
        let abund = tf.abundance;
        for &nid in &real_nodes {
            let nid_coord = graph.nodes.get(nid).map(|n| n.start).unwrap_or(0);
            if first_coord < nid_coord {
                abundin[nid] += abund;
            }
            if last_coord > nid_coord {
                abundout[nid] += abund;
            }
            if tf.longread {
                longcov[nid] += abund;
            }
        }
    }

    let nodecov: Vec<f64> = abundin
        .iter()
        .zip(abundout.iter())
        .map(|(a, b)| a.max(*b))
        .collect();

    // Debug: print nc summary for large bundles when NODECOV_DEBUG env var is set.
    if std::env::var("NODECOV_DEBUG").is_ok() && n_nodes > 80 {
        // Print node coordinate range to identify the bundle
        let first_real = graph
            .nodes
            .get(1)
            .map(|n| (n.start, n.end))
            .unwrap_or((0, 0));
        let last_real = graph
            .nodes
            .get(n_nodes.saturating_sub(2))
            .map(|n| (n.start, n.end))
            .unwrap_or((0, 0));
        eprintln!(
            "[NODECOV_DEBUG] gno={} coord_range={}..{}",
            n_nodes, first_real.0, last_real.1
        );
        let active: Vec<_> = transfrags
            .iter()
            .filter(|tf| {
                if tf.abundance <= 0.0 {
                    return false;
                }
                let nodes = &tf.node_ids;
                if nodes.is_empty() {
                    return false;
                }
                let fn_ = nodes[0];
                let ln_ = *nodes.last().unwrap();
                let full = fn_ == source_id && ln_ == sink_id;
                if full {
                    return false;
                }
                if fn_ == source_id || ln_ == sink_id {
                    return false;
                }
                true
            })
            .collect();
        eprintln!(
            "[NODECOV_DEBUG] gno={} active_internal_tfs={}",
            n_nodes,
            active.len()
        );
        for (i, &tf) in active.iter().enumerate().take(30) {
            eprintln!(
                "  tf[{}] nodes=[{}..{}] len={} abund={:.1} guide={} longread={}",
                i,
                tf.node_ids[0],
                tf.node_ids.last().unwrap(),
                tf.node_ids.len(),
                tf.abundance,
                tf.guide,
                tf.longread
            );
        }
        if active.len() > 30 {
            eprintln!("  ... ({} more)", active.len() - 30);
        }
        // Print nc for all nodes with nc > 0
        eprintln!("[NODECOV_DEBUG] per-node nc (non-zero):");
        for (nid, &nc) in nodecov.iter().enumerate() {
            if nc > 0.0 {
                let bncov = graph.nodes.get(nid).map(|n| n.coverage).unwrap_or(0.0);
                let (ns, ne) = graph
                    .nodes
                    .get(nid)
                    .map(|n| (n.start, n.end))
                    .unwrap_or((0, 0));
                eprintln!(
                    "  node[{}] {}..{} nc={:.1} bn_cov={:.1} rate={:.4}",
                    nid,
                    ns,
                    ne,
                    nc,
                    bncov,
                    if nc > 0.0 { bncov / nc } else { 1.0 }
                );
            }
        }

        // For gno=264 specifically, dump which transfrags hit specific nodes of interest
        if n_nodes == 264 {
            let nodes_of_interest = [37usize, 111, 112, 116, 37]; // node 37 and STRG.15.1 range
            let unique_noi: Vec<usize> = {
                let mut v = nodes_of_interest.to_vec();
                v.sort_unstable();
                v.dedup();
                v
            };
            for &target in &unique_noi {
                let (ts, te) = graph
                    .nodes
                    .get(target)
                    .map(|n| (n.start, n.end))
                    .unwrap_or((0, 0));
                eprintln!(
                    "[NODECOV_DEBUG] contributors to node[{}] {}..{}:",
                    target, ts, te
                );
                for (i, tf) in transfrags.iter().enumerate() {
                    if tf.node_ids.contains(&target) {
                        let fn_ = tf.node_ids[0];
                        let ln_ = *tf.node_ids.last().unwrap();
                        let contrib_in = target > fn_;
                        let contrib_out = target < ln_;
                        if contrib_in || contrib_out {
                            eprintln!(
                                "  tf[{}] first={} last={} abund={:.1} contrib={{in:{},out:{}}}",
                                i, fn_, ln_, tf.abundance, contrib_in, contrib_out
                            );
                        }
                    }
                }
            }
        }
    }

    // Persist per-node coverage state with ref's node accounting.
    for nid in 0..n_nodes {
        if let Some(node) = graph.nodes.get_mut(nid) {
            let ai = abundin.get(nid).copied().unwrap_or(0.0);
            let ao = abundout.get(nid).copied().unwrap_or(0.0);
            let nc = nodecov.get(nid).copied().unwrap_or(0.0);
            let lc = longcov.get(nid).copied().unwrap_or(0.0);
            node.abundin = ai;
            node.abundout = ao;
            node.nodecov = nc;
            node.longcov = lc;
            node.noderate = if nc > 0.0 {
                let nr = node.coverage / nc;
                if nr.is_finite() && nr > 0.0 {
                    nr
                } else {
                    1.0
                }
            } else {
                1.0
            };
        }
    }

    if verbose {
        let nonzero = nodecov.iter().filter(|&&c| c > 0.0).count();
        let total: f64 = nodecov.iter().sum();
        let sum_in: f64 = abundin.iter().sum();
        let sum_out: f64 = abundout.iter().sum();
        let sum_long: f64 = longcov.iter().sum();
        let rate_nonunit = graph
            .nodes
            .iter()
            .filter(|n| (n.noderate - 1.0).abs() > 1e-9)
            .count();
        eprintln!(
            "      [Rustle] compute_nodecov: gno={} covered={}/{} total={:.1} abundin={:.1} abundout={:.1} longcov={:.1} noderate_nonunit={}",
            n_nodes,
            nonzero,
            n_nodes,
            total,
            sum_in,
            sum_out,
            sum_long,
            rate_nonunit
        );
    }
    nodecov
}

/// Mark transfrags that have a weak link (coverage drop between consecutive contiguous nodes) (compute_weak).
/// threshold = drop + error_perc (e.g. 0.6). Returns count of transfrags marked weak.
pub fn compute_weak(
    transfrags: &mut [GraphTransfrag],
    indices: &[usize],
    graph: &Graph,
    nodecov: &[f64],
    drop: f64,
    error_perc: f64,
    verbose: bool,
) -> usize {
    let threshold = drop + error_perc;
    let mut weak_count = 0usize;
    for &tf_idx in indices {
        if tf_idx >= transfrags.len() {
            continue;
        }
        let tf = &mut transfrags[tf_idx];
        tf.coverage_weak = false;
        let nodes = &tf.node_ids;
        if nodes.len() < 2 {
            continue;
        }
        for i in 0..nodes.len() - 1 {
            let nid_a = nodes[i];
            let nid_b = nodes[i + 1];
            let (n1, n2) = match (graph.nodes.get(nid_a), graph.nodes.get(nid_b)) {
                (Some(a), Some(b)) => (a, b),
                _ => continue,
            };
            // StringTie checks `nodes[n] == nodes[n-1]+1 && end+1 == start`. Node-ID
            // adjacency is broken after longtrim, so use coordinate adjacency only.
            // Rustle uses half-open coords: contiguous means `end == start` (not `end+1 == start`).
            if n1.end != n2.start {
                continue;
            }
            let cov_a = nodecov.get(nid_a).copied().unwrap_or(0.0);
            let cov_b = nodecov.get(nid_b).copied().unwrap_or(0.0);
            if (cov_a * threshold > cov_b) || (cov_b * threshold > cov_a) {
                tf.coverage_weak = true;
                weak_count += 1;
                break;
            }
        }
    }
    if verbose && weak_count > 0 {
        eprintln!(
            "      [Rustle] compute_weak: {}/{} transfrags with weak links (threshold={})",
            weak_count,
            indices.len(),
            threshold
        );
    }
    weak_count
}
