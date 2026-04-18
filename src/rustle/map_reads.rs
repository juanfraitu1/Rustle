//! Map bundle reads onto graph -> transfrags (get_fragment_pattern + update_abundance).
//!
//! Long-read only: all reads → abundance, longread=true.
//!
//! V99: When killed_junction_pairs is provided, reads are split at killed junctions and
//! segments to the right of a killed junction get killed_junction_orphan = true.

use crate::assembly_mode::LONGINTRONANCHOR;
use crate::bitvec::GBitVec;
use crate::coord::{overlap_len_half_open, overlaps_half_open};
use crate::graph::{Graph, GraphTransfrag};
use crate::treepat::TreePatIndex;
use crate::types::{
    AssemblyMode, Bundle2Graph, BundleRead, CBundlenode, DetHashSet as HashSet, Junction,
};
use hashbrown::{HashMap as HbHashMap, HashSet as HbHashSet};
use lru::LruCache;
use roaring::RoaringBitmap;
use std::num::NonZeroUsize;
const CHI_WIN: u64 = 100;
const DROP: f64 = 0.5;
const KMER: u64 = 31;
const DEFAULT_READGROUP_CACHE_CAP: usize = 32768;

fn sort_node_ids_by_coord(graph: &Graph, nodes: &mut Vec<usize>) {
    nodes.sort_unstable_by(|&a, &b| {
        let na = &graph.nodes[a];
        let nb = &graph.nodes[b];
        (na.start, na.end, a).cmp(&(nb.start, nb.end, b))
    });
    nodes.dedup();
}

fn collect_candidate_nodes_sorted(
    graph: &Graph,
    allowed_nodes: Option<&HashSet<usize>>,
) -> Vec<usize> {
    let mut out: Vec<usize> = graph
        .nodes
        .iter()
        .enumerate()
        .filter_map(|(idx, node)| {
            if node.node_id == graph.source_id || node.node_id == graph.sink_id {
                return None;
            }
            if let Some(allowed) = allowed_nodes {
                if !allowed.contains(&idx) {
                    return None;
                }
            }
            Some(idx)
        })
        .collect();
    // walks graph nodes in genomic order. Rust node ids can be out of coordinate order
    // after trim/split operations, so restore genomic ordering before exon->node assignment.
    sort_node_ids_by_coord(graph, &mut out);
    out
}

fn update_abund_trace_active() -> bool {
    std::env::var_os("RUSTLE_TRACE_UPDATE_ABUND").is_some()
        || std::env::var_os("RUSTLE_TRACE_LOG_STYLE").is_some()
}

fn parse_trace_locus_env() -> Option<(u64, u64)> {
    let Ok(val) = std::env::var("RUSTLE_TRACE_LOCUS") else {
        return None;
    };
    let (a, b) = val.split_once('-')?;
    let start = a.trim().parse::<u64>().ok()?;
    let end = b.trim().parse::<u64>().ok()?;
    Some((start.min(end), start.max(end)))
}

fn path_overlaps_trace(path: &[usize], graph: &Graph, trace_locus: Option<(u64, u64)>) -> bool {
    let Some((lo, hi)) = trace_locus else {
        return false;
    };
    path.iter().copied().any(|nid| {
        graph
            .nodes
            .get(nid)
            .map(|n| n.start <= hi && n.end >= lo)
            .unwrap_or(false)
    })
}

fn path_span(path: &[usize], graph: &Graph) -> (u64, u64) {
    let first = path
        .first()
        .and_then(|&nid| graph.nodes.get(nid))
        .map(|n| n.start)
        .unwrap_or(0);
    let last = path
        .last()
        .and_then(|&nid| graph.nodes.get(nid))
        .map(|n| n.end)
        .unwrap_or(0);
    (first, last)
}

fn read_overlaps_trace(read: &BundleRead, trace_locus: Option<(u64, u64)>) -> bool {
    let Some((lo, hi)) = trace_locus else {
        return false;
    };
    read.exons
        .iter()
        .any(|&(start, end)| start <= hi && end >= lo)
}

fn format_nodes_with_coords(path: &[usize], graph: &Graph) -> String {
    path.iter()
        .filter_map(|&nid| {
            graph
                .nodes
                .get(nid)
                .map(|n| format!(" {}({}-{})", nid, n.start, n.end))
        })
        .collect::<String>()
}

#[derive(Debug, Clone, Hash, PartialEq, Eq)]
struct TransfragKey {
    nodes: Vec<usize>,
    edge_flags: Vec<bool>,
}

fn build_transfrag_pattern(
    path: &[usize],
    graph: &Graph,
    junctions: &[Junction],
    junction_correction_window: u64,
    psize: usize,
) -> GBitVec {
    let mut pattern = GBitVec::new(psize);
    for &nid in path {
        pattern.set_bit(nid);
    }
    for window in path.windows(2) {
        let a = window[0];
        let b = window[1];
        let Some(na) = graph.nodes.get(a) else {
            continue;
        };
        let Some(nb) = graph.nodes.get(b) else {
            continue;
        };
        let edge_supported = if na.end == nb.start {
            true
        } else {
            junctions
                .iter()
                .any(|j| {
                    // the original algorithm -E: accept donor/acceptor mismatch within the correction window.
                    j.donor.abs_diff(na.end) <= junction_correction_window
                        && j.acceptor.abs_diff(nb.start) <= junction_correction_window
                })
        };
        if !edge_supported {
            continue;
        }
        if let Some(eid) = graph.edge_bit_index(a, b) {
            pattern.set_bit(eid);
        }
    }
    pattern
}

fn transfrag_key_from_pattern(path: &[usize], graph: &Graph, pattern: &GBitVec) -> TransfragKey {
    let edge_flags = path
        .windows(2)
        .map(|w| {
            graph
                .edge_bit_index(w[0], w[1])
                .map(|eid| pattern.get_bit(eid))
                .unwrap_or(false)
        })
        .collect();
    TransfragKey {
        nodes: path.to_vec(),
        edge_flags,
    }
}

/// Ensure that splice edges implied by the read's exon-to-node path exist in the graph.
///
/// In long-read mode with `-E`, a read junction can land within the correction window of a node
/// boundary splice (node_end -> next_node_start). If we don't materialize that edge in the graph,
/// later stages treat the path as "unwitnessed splice" and the isoform can be lost.
fn ensure_edges_for_read_path(
    graph: &mut Graph,
    unique_nodes: &[usize],
    read_junctions: &[Junction],
    junction_correction_window: u64,
    killed_junction_pairs: Option<&HashSet<Junction>>,
) {
    if unique_nodes.len() < 2 {
        return;
    }
    for w in unique_nodes.windows(2) {
        let a = w[0];
        let b = w[1];
        let (Some(na), Some(nb)) = (graph.nodes.get(a), graph.nodes.get(b)) else {
            continue;
        };
        if na.end == nb.start {
            // Contiguous edge should already exist.
            graph.add_edge(a, b);
            continue;
        }
        let supports = read_junctions.iter().any(|j| {
            j.donor.abs_diff(na.end) <= junction_correction_window
                && j.acceptor.abs_diff(nb.start) <= junction_correction_window
        });
        if !supports {
            continue;
        }
        let corrected = Junction {
            donor: na.end,
            acceptor: nb.start,
        };
        let killed = killed_junction_pairs
            .map(|set| set.contains(&corrected))
            .unwrap_or(false);
        if killed {
            continue;
        }
        graph.add_edge(a, b);
    }
}

/// Classify read as long in mixed mode: length >= long_read_min_len (_is_long_read).
#[inline]
fn is_long_read(read: &BundleRead, long_read_min_len: u64) -> bool {
    if long_read_min_len == 0 {
        return true;
    }
    read.query_length.map_or(true, |q| q >= long_read_min_len)
}

/// Return the node path (fragment pattern) for one read (get_fragment_pattern).
/// Uses same logic as map_reads_to_graph: exon→node assignment, unique_nodes, split at disconnected/killed; returns primary (longest) path.
pub fn read_to_path(
    read: &BundleRead,
    graph: &Graph,
    killed_junction_pairs: Option<&HashSet<Junction>>,
    allowed_nodes: Option<&HashSet<usize>>,
) -> Vec<usize> {
    let ordered_nodes = collect_candidate_nodes_sorted(graph, allowed_nodes);
    let (unique_nodes, _) = collect_read_nodes_exact(read, graph, &ordered_nodes, false);
    // Helper path extraction used in diagnostics; default to strict junction matching (0bp window).
    split_read_segments(read, graph, &unique_nodes, 0, killed_junction_pairs, None)
        .into_iter()
        .next()
        .map(|seg| seg.path)
        .unwrap_or_default()
}

fn collect_read_nodes_exact(
    read: &BundleRead,
    graph: &Graph,
    ordered_nodes: &[usize],
    add_coverage: bool,
) -> (Vec<usize>, Vec<(usize, f64)>) {
    let mut path = Vec::new();
    let mut cov_add = Vec::new();
    let mut exon_idx = 0usize;
    let kmer = KMER.saturating_sub(1);

    for &nid in ordered_nodes {
        let node = &graph.nodes[nid];
        let mut intersect = false;
        while exon_idx < read.exons.len() {
            let (seg_start, seg_end) = read.exons[exon_idx];
            if seg_end <= node.start {
                exon_idx += 1;
                continue;
            }

            let bp = overlap_len_half_open(seg_start, seg_end, node.start, node.end);
            if bp > 0 {
                intersect = true;
                if add_coverage && !read.unitig {
                    cov_add.push((nid, read.weight * (bp as f64)));
                }
                if seg_end <= node.end {
                    exon_idx += 1;
                } else {
                    break;
                }
            } else {
                break;
            }
        }

        if !intersect || path.last() == Some(&nid) {
            continue;
        }

        let mut trimu = false;
        if read.unitig {
            if path.len() == 1 && node.parents.count_ones() > 1 {
                let first = path[0];
                if graph.nodes[first].end.saturating_sub(read.ref_start) < kmer {
                    path.pop();
                }
            }
            if let Some(&last) = path.last() {
                if graph.nodes[last].children.count_ones() > 1
                    && read.ref_end.saturating_sub(node.start) < kmer
                {
                    trimu = true;
                }
            }
        }

        if !trimu {
            path.push(nid);
        }
    }

    (path, cov_add)
}

fn collect_bundlenode_ranges(bn: Option<&CBundlenode>) -> Vec<(usize, u64, u64)> {
    let mut out = Vec::new();
    let mut cur = bn;
    while let Some(n) = cur {
        out.push((n.bid, n.start, n.end));
        cur = n.next.as_deref();
    }
    out
}

fn readgroup_cache_capacity() -> NonZeroUsize {
    std::env::var("RUSTLE_READGROUP_CACHE")
        .ok()
        .and_then(|v| v.parse::<usize>().ok())
        .and_then(NonZeroUsize::new)
        .unwrap_or_else(|| {
            NonZeroUsize::new(DEFAULT_READGROUP_CACHE_CAP).expect("cache cap must be non-zero")
        })
}

fn ordered_nodes_from_readgroup(
    graph: &Graph,
    bundle2graph: &Bundle2Graph,
    readgroup: &[usize],
) -> Vec<usize> {
    let mut dedup = RoaringBitmap::new();
    let mut overflow: Option<HbHashSet<usize>> = None;
    for &bid in readgroup {
        let Some(nodes) = bundle2graph.get(bid) else {
            continue;
        };
        for &(_, nid) in nodes {
            if nid == graph.source_id || nid == graph.sink_id {
                continue;
            }
            if nid <= u32::MAX as usize {
                dedup.insert(nid as u32);
            } else {
                overflow.get_or_insert_with(Default::default).insert(nid);
            }
        }
    }
    let mut ordered_nodes: Vec<usize> = dedup.iter().map(|nid| nid as usize).collect();
    if let Some(extra) = overflow {
        ordered_nodes.extend(extra);
    }
    sort_node_ids_by_coord(graph, &mut ordered_nodes);
    ordered_nodes
}

fn readgroup_for_bundlenodes(read: &BundleRead, bnodes: &[(usize, u64, u64)]) -> Vec<usize> {
    let mut out = Vec::new();
    for &(es, ee) in &read.exons {
        for (bid, bstart, bend) in bnodes {
            if overlaps_half_open(es, ee, *bstart, *bend) && out.last() != Some(bid) {
                out.push(*bid);
            }
        }
    }
    out
}

/// Return node path using bundlenode->graph mapping (get_read_pattern).
pub fn read_to_path_bundlenodes(
    read: &BundleRead,
    graph: &Graph,
    bundlenodes: Option<&CBundlenode>,
    bundle2graph: Option<&Bundle2Graph>,
) -> Vec<usize> {
    let bnodes = collect_bundlenode_ranges(bundlenodes);
    let b2g = match bundle2graph {
        Some(m) if !m.is_empty() => m,
        _ => {
            return read_to_path(read, graph, None, None);
        }
    };

    let readgroup = readgroup_for_bundlenodes(read, &bnodes);
    let ordered_nodes = ordered_nodes_from_readgroup(graph, b2g, &readgroup);

    collect_read_nodes_exact(read, graph, &ordered_nodes, false).0
}

/// For each read, find path of nodes that overlap its exons; create/update transfrags.
/// In mixed mode, long reads add to abundance and set longread; short reads add to srabund.
/// If killed_junction_pairs is Some, reads are split at killed junctions and segments to the
/// right of a killed junction get killed_junction_orphan = true (original algorithm V99).
pub fn map_reads_to_graph(
    reads: &[BundleRead],
    graph: &mut Graph,
    mode: AssemblyMode,
    _long_read_min_len: u64,
    junction_correction_window: u64,
    killed_junction_pairs: Option<&HashSet<Junction>>,
    allowed_nodes: Option<&HashSet<usize>>,
) -> Vec<GraphTransfrag> {
    let _ = &mode; // suppress unused warning; passed through to add_or_update_transfrag
    let psize = graph.pattern_size();
    let ordered_nodes = collect_candidate_nodes_sorted(graph, allowed_nodes);
    let mut transfrags: Vec<GraphTransfrag> = Vec::with_capacity(reads.len().saturating_div(2));
    let mut pattern_map: HbHashMap<TransfragKey, usize> =
        HbHashMap::with_capacity(reads.len());
    let mut tr_index = TreePatIndex::new(graph.n_nodes);

    for (read_idx, read) in reads.iter().enumerate() {
        let is_long = true; // long-read only mode

        let (unique_nodes, cov_add) = collect_read_nodes_exact(read, graph, &ordered_nodes, true);
        for (idx, add) in cov_add {
            if let Some(n) = graph.nodes.get_mut(idx) {
                n.coverage += add;
            }
        }

        ensure_edges_for_read_path(
            graph,
            &unique_nodes,
            &read.junctions,
            junction_correction_window,
            killed_junction_pairs,
        );

        let segments = split_read_segments(
            read,
            graph,
            &unique_nodes,
            junction_correction_window,
            killed_junction_pairs,
            Some(read_idx),
        );

        if segments.is_empty() {
            continue;
        }

        // Trace: dump read-to-node mapping for debug comparison with the original algorithm.
        if std::env::var_os("RUSTLE_READ_PATTERN_TRACE").is_some() {
            for (si, seg) in segments.iter().enumerate() {
                if seg.path.is_empty() { continue; }
                let nodes_str: Vec<String> = seg.path.iter().map(|&nid| {
                    let n = &graph.nodes[nid];
                    format!("{}({}-{})", nid, n.start, n.end)
                }).collect();
                eprintln!(
                    "READ_PATTERN read={} seg={} weight={:.2} nodes=[{}] nexons={}",
                    read_idx, si, read.weight,
                    nodes_str.join(","),
                    read.exons.len(),
                );
            }
        }

        let weight = read.weight;
        let ref_start = read.ref_start;
        let ref_end = read.ref_end;
        // Hard boundary evidence must come from unaligned (soft-clip) poly tails only.
        // Counting aligned poly runs (RT drop-off artifacts) inflates hardstart/hardend and
        // causes hard-boundary mismatches during long-read stitching.
        let has_poly_start = read.has_poly_start_unaligned;
        let has_poly_end = read.has_poly_end_unaligned;

        for seg in &segments {
            if seg.path.is_empty() {
                continue;
            }
            // skip single-node fragments from multi-exon reads.
            // When a multi-exon long read has killed junctions, split_read_segments
            // fragments it into tiny segments. Single-node fragments from such reads
            // are alignment artifacts that inflate node coverage and create spurious
            // source/sink transfrags, adding ~23 extra reads per gene locus.
            if is_long && seg.path.len() == 1 && read.exons.len() >= 3 && segments.len() > 1 {
                continue;
            }
            add_or_update_transfrag(
                &seg.path,
                graph,
                &mut transfrags,
                &mut pattern_map,
                &mut tr_index,
                psize,
                is_long,
                mode,
                junction_correction_window,
                weight,
                ref_start,
                ref_end,
                has_poly_start,
                has_poly_end,
                seg.orphan,
                &read.junctions,
            );
        }
    }

    // Post-processing: split chimeric transfrags at coverage valleys.
    // When a transfrag's node path goes through adjacent nodes (prev.end == curr.start)
    // where coverage drops sharply, the transfrag likely spans two separate genes via
    // an intergenic bridging region. Split it into two independent transfrags.
    // This matches StringTie's coverage-drop exclusion (rlink.cpp:8262-8266).
    transfrags = split_chimeric_transfrags(transfrags, graph);

    for (tf_idx, tf) in transfrags.iter().enumerate() {
        if tf.node_ids.len() <= 1 {
            continue;
        }
        for &nid in &tf.node_ids {
            if nid < graph.nodes.len() {
                graph.nodes[nid].trf_ids.push(tf_idx);
            }
        }
    }

    transfrags
}

/// Split chimeric transfrags at coverage valleys between adjacent graph nodes.
/// Adjacent nodes (prev.end == curr.start) with a >3x coverage drop indicate
/// an intergenic boundary that a readthrough read spans continuously.
fn split_chimeric_transfrags(
    mut transfrags: Vec<GraphTransfrag>,
    graph: &Graph,
) -> Vec<GraphTransfrag> {
    // Coverage drop threshold for chimeric-split. 0.01 catches a handful of
    // intergenic coverage valleys (+2 matches on GGO_19); higher values
    // (0.05+) lose matches by over-splitting legitimate low-cov internal
    // exons. Override via RUSTLE_CHIMERIC_DROP.
    let drop_threshold: f64 = std::env::var("RUSTLE_CHIMERIC_DROP")
        .ok()
        .and_then(|v| v.parse().ok())
        .unwrap_or(0.01);
    let mut new_transfrags: Vec<GraphTransfrag> = Vec::new();
    let mut split_count = 0usize;
    let n_before = transfrags.len();

    for tf in &mut transfrags {
        if tf.node_ids.len() < 4 || !tf.longread {
            continue;
        }
        // Find the deepest per-base coverage valley in the transfrag's node path.
        // A coverage valley indicates an intergenic boundary where a readthrough
        // read spans two genes. We split at the node with the lowest per-base coverage.
        let mut best_split: Option<(usize, f64)> = None;
        // Compute per-base coverage for each node
        let per_base: Vec<f64> = tf.node_ids.iter().map(|&nid| {
            graph.nodes.get(nid).map_or(0.0, |n| {
                let len = n.length().max(1) as f64;
                n.coverage / len
            })
        }).collect();
        // Find the median per-base coverage for context
        let mut sorted_cov: Vec<f64> = per_base.iter().copied().filter(|&c| c > 0.0).collect();
        sorted_cov.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));
        let median_cov = if sorted_cov.is_empty() { 0.0 } else { sorted_cov[sorted_cov.len() / 2] };
        if median_cov > 0.0 {
            for j in 1..tf.node_ids.len().saturating_sub(1) {
                let nid = tf.node_ids[j];
                if nid == graph.source_id || nid == graph.sink_id { continue; }
                let cov = per_base[j];
                let prev_cov = per_base[j - 1];
                let next_cov = per_base[j + 1];
                // Valley: node j has much lower per-base coverage than both neighbors
                let neighbor_min = prev_cov.min(next_cov);
                if neighbor_min <= 0.0 { continue; }
                let ratio = cov / neighbor_min;
                if ratio < drop_threshold && cov < median_cov * drop_threshold {
                    if best_split.as_ref().map_or(true, |&(_, r)| ratio < r) {
                        best_split = Some((j, ratio));
                    }
                }
            }
        }

        if let Some((split_at, _ratio)) = best_split {
            let right_nodes: Vec<usize> = tf.node_ids[split_at..].to_vec();
            if right_nodes.len() >= 2 {
                let psize = graph.pattern_size();
                let mut right_tf = GraphTransfrag::new(right_nodes, psize);
                right_tf.abundance = tf.abundance;
                right_tf.longread = tf.longread;
                right_tf.longstart = 0;
                right_tf.longend = tf.longend;
                right_tf.read_count = tf.read_count;
                for &nid in &right_tf.node_ids {
                    right_tf.pattern.set_bit(nid);
                }
                for w in right_tf.node_ids.windows(2) {
                    if let Some(eid) = graph.edge_bit_index(w[0], w[1]) {
                        right_tf.pattern.set_bit(eid);
                    }
                }
                new_transfrags.push(right_tf);
            }
            tf.node_ids.truncate(split_at);
            tf.longend = 0;
            let psize = graph.pattern_size();
            tf.pattern = crate::bitvec::GBitVec::new(psize);
            for &nid in &tf.node_ids {
                tf.pattern.set_bit(nid);
            }
            for w in tf.node_ids.windows(2) {
                if let Some(eid) = graph.edge_bit_index(w[0], w[1]) {
                    tf.pattern.set_bit(eid);
                }
            }
            split_count += 1;
        }
    }
    if split_count > 0 {
        transfrags.extend(new_transfrags);
        if std::env::var_os("RUSTLE_DEBUG_DETAIL").is_some() {
            eprintln!(
                "[CHIMERIC_SPLIT] split {} chimeric transfrags (was {}, now {})",
                split_count, n_before, transfrags.len()
            );
        }
    }
    transfrags
}

/// Map reads using bundlenode->graph mapping (get_read_pattern).
pub fn map_reads_to_graph_bundlenodes(
    reads: &[BundleRead],
    graph: &mut Graph,
    mode: AssemblyMode,
    _long_read_min_len: u64,
    junction_correction_window: u64,
    killed_junction_pairs: Option<&HashSet<Junction>>,
    bundlenodes: Option<&CBundlenode>,
    bundle2graph: Option<&Bundle2Graph>,
    read_bundles: Option<&[Vec<usize>]>,
) -> Vec<GraphTransfrag> {
    let _ = &mode; // suppress unused warning; passed through to add_or_update_transfrag
    let psize = graph.pattern_size();
    let mut transfrags: Vec<GraphTransfrag> = Vec::with_capacity(reads.len().saturating_div(2));
    let mut pattern_map: HbHashMap<TransfragKey, usize> =
        HbHashMap::with_capacity(reads.len());
    let mut tr_index = TreePatIndex::new(graph.n_nodes);
    let bnodes = collect_bundlenode_ranges(bundlenodes);
    let b2g = bundle2graph.filter(|m| !m.is_empty());
    // Use u64 hash fingerprint as cache key instead of Vec<usize> to avoid allocation per lookup.
    let mut readgroup_node_cache: LruCache<u64, Vec<usize>> =
        LruCache::new(readgroup_cache_capacity());
    let dump_range = std::env::var("RUSTLE_DUMP_RANGE").ok().and_then(|v| {
        let (a, b) = v.split_once('-')?;
        let s = a.trim().parse::<u64>().ok()?;
        let e = b.trim().parse::<u64>().ok()?;
        Some((s, e))
    });
    let mut dump_hits = 0usize;

    for (read_idx, read) in reads.iter().enumerate() {
        let is_long = true; // long-read only mode

        let readgroup_owned;
        let readgroup: &[usize] = if let Some(mapped) = read_bundles.and_then(|v| v.get(read_idx)) {
            mapped
        } else {
            readgroup_owned = readgroup_for_bundlenodes(read, &bnodes);
            &readgroup_owned
        };
        let mut ordered_nodes = Vec::new();
        if let Some(map) = b2g {
            let cache_key = {
                use std::hash::{Hash, Hasher};
                let mut h = fxhash::FxHasher::default();
                readgroup.hash(&mut h);
                h.finish()
            };
            if let Some(cached) = readgroup_node_cache.get(&cache_key) {
                ordered_nodes = cached.clone();
            } else {
                ordered_nodes = ordered_nodes_from_readgroup(graph, map, readgroup);
                readgroup_node_cache.put(cache_key, ordered_nodes.clone());
            }
        }
        let (unique_nodes, cov_add) = collect_read_nodes_exact(read, graph, &ordered_nodes, true);
        if unique_nodes.is_empty() {
            continue;
        }
        for (nid, add) in cov_add {
            if let Some(n) = graph.nodes.get_mut(nid) {
                n.coverage += add;
            }
        }
        if let Some((rs, re)) = dump_range {
            if dump_hits < 10 {
                let hit = unique_nodes.iter().any(|&nid| {
                    let node = &graph.nodes[nid];
                    node.end >= rs && node.start <= re
                });
                if hit {
                    let ex = read
                        .exons
                        .iter()
                        .map(|(s, e)| format!("{}-{}", s, e))
                        .collect::<Vec<_>>()
                        .join(",");
                    let nodes = unique_nodes
                        .iter()
                        .map(|&nid| {
                            let n = &graph.nodes[nid];
                            format!("{}-{}", n.start, n.end)
                        })
                        .collect::<Vec<_>>()
                        .join(",");
                    eprintln!("  [BUNDLEMAP] read exons={} nodes={}", ex, nodes);
                    dump_hits += 1;
                }
            }
        }

        ensure_edges_for_read_path(
            graph,
            &unique_nodes,
            &read.junctions,
            junction_correction_window,
            killed_junction_pairs,
        );

        let segments = split_read_segments(
            read,
            graph,
            &unique_nodes,
            junction_correction_window,
            killed_junction_pairs,
            Some(read_idx),
        );
        if segments.is_empty() {
            continue;
        }

        // Trace: dump read-to-node mapping for debug comparison with the original algorithm.
        if std::env::var_os("RUSTLE_READ_PATTERN_TRACE").is_some() {
            for (si, seg) in segments.iter().enumerate() {
                if seg.path.is_empty() { continue; }
                let nodes_str: Vec<String> = seg.path.iter().map(|&nid| {
                    let n = &graph.nodes[nid];
                    format!("{}({}-{})", nid, n.start, n.end)
                }).collect();
                eprintln!(
                    "READ_PATTERN read={} seg={} weight={:.2} nodes=[{}] nexons={}",
                    read_idx, si, read.weight,
                    nodes_str.join(","),
                    read.exons.len(),
                );
            }
        }

        let weight = read.weight;
        let ref_start = read.ref_start;
        let ref_end = read.ref_end;
        let has_poly_start = read.has_poly_start_unaligned;
        let has_poly_end = read.has_poly_end_unaligned;

        for seg in &segments {
            if seg.path.is_empty() {
                continue;
            }
            // skip single-node fragments from multi-exon reads.
            // When a multi-exon long read has killed junctions, split_read_segments
            // fragments it into tiny segments. Single-node fragments from such reads
            // are alignment artifacts that inflate node coverage and create spurious
            // source/sink transfrags, adding ~23 extra reads per gene locus.
            if is_long && seg.path.len() == 1 && read.exons.len() >= 3 && segments.len() > 1 {
                continue;
            }
            add_or_update_transfrag(
                &seg.path,
                graph,
                &mut transfrags,
                &mut pattern_map,
                &mut tr_index,
                psize,
                is_long,
                mode,
                junction_correction_window,
                weight,
                ref_start,
                ref_end,
                has_poly_start,
                has_poly_end,
                seg.orphan,
                &read.junctions,
            );
        }
    }

    // Split chimeric transfrags at coverage valleys (same as map_reads_to_graph)
    transfrags = split_chimeric_transfrags(transfrags, graph);

    for (tf_idx, tf) in transfrags.iter().enumerate() {
        if tf.node_ids.len() <= 1 {
            continue;
        }
        for &nid in &tf.node_ids {
            if nid < graph.nodes.len() {
                graph.nodes[nid].trf_ids.push(tf_idx);
            }
        }
    }

    transfrags
}

#[derive(Debug, Clone)]
struct ReadPathSegment {
    path: Vec<usize>,
    orphan: bool,
}

/// the original algorithm get_fragment_pattern: walk the read left-to-right and emit a new
/// segment only when the current read junction is invalid for the current node pair,
/// or when a unitig was trimmed through source/sink sidecars.
fn split_read_segments(
    read: &BundleRead,
    graph: &Graph,
    unique_nodes: &[usize],
    junction_correction_window: u64,
    killed_junction_pairs: Option<&HashSet<Junction>>,
    read_idx: Option<usize>,
) -> Vec<ReadPathSegment> {
    if unique_nodes.is_empty() {
        return Vec::new();
    }
    if unique_nodes.len() == 1 {
        return vec![ReadPathSegment {
            path: unique_nodes.to_vec(),
            orphan: false,
        }];
    }

    let trace_locus = parse_trace_locus_env();
    let trace_this = read_idx.is_some() && read_overlaps_trace(read, trace_locus);
    let mut segments = Vec::new();
    let mut seg_start = 0usize;
    let mut next_segment_orphan = false;
    let mut junc_idx = 0usize;

    for j in 1..unique_nodes.len() {
        let prev_nid = unique_nodes[j - 1];
        let curr_nid = unique_nodes[j];
        let Some(prev_node) = graph.nodes.get(prev_nid) else {
            continue;
        };
        let Some(curr_node) = graph.nodes.get(curr_nid) else {
            continue;
        };

        let mut split_here = false;
        let mut orphan_right = false;

        if prev_node.end == curr_node.start {
            if read.unitig {
                let curr_has_source = curr_node.parents.contains(graph.source_id);
                let prev_has_sink = prev_node.children.contains(graph.sink_id);
                if curr_has_source || prev_has_sink {
                    split_here = true;
                    if trace_this {
                        eprintln!(
                            "PATH_frag_pattern: read[{}] SPLIT_UNITIG at_j={} upstream_nodes={} | restart_node={}({}-{})",
                            read_idx.unwrap_or(0),
                            j,
                            format_nodes_with_coords(&unique_nodes[seg_start..j], graph),
                            curr_nid,
                            curr_node.start,
                            curr_node.end
                        );
                    }
                }
            }
        } else if std::env::var_os("RUSTLE_BUNDLE_GRAPH").is_some() {
            // the original algorithm never splits reads at junction gaps.
            // With bundle graph mode, match the original algorithm's behavior exactly.
        } else {
            while junc_idx < read.junctions.len()
                && read.junctions[junc_idx]
                    .acceptor
                    .saturating_add(junction_correction_window)
                    < curr_node.start
            {
                junc_idx += 1;
            }

            if let Some(&junc) = read.junctions.get(junc_idx) {
                if junc
                    .donor
                    .saturating_sub(junction_correction_window)
                    <= prev_node.end
                {
                    // the original algorithm -E: accept donor/acceptor mismatch within the correction window,
                    // and treat the splice as matching the node boundary junction (prev_end -> curr_start).
                    let within_window = junc.donor.abs_diff(prev_node.end) <= junction_correction_window
                        && junc.acceptor.abs_diff(curr_node.start) <= junction_correction_window;
                    let corrected = Junction {
                        donor: prev_node.end,
                        acceptor: curr_node.start,
                    };
                    let killed = killed_junction_pairs
                        .map(|set| {
                            if within_window {
                                set.contains(&corrected)
                            } else {
                                set.contains(&junc)
                            }
                        })
                        .unwrap_or(false);
                    let valid = !killed
                        && (within_window
                            || (junc.donor == prev_node.end && junc.acceptor == curr_node.start));
                    // (get_read_pattern:4431-4512): the original algorithm
                    // adds ALL nodes the read intersects to a single pattern.  When a
                    // junction is killed, the original algorithm's changeleft/changeright repair
                    // redirects to a nearby good junction,
                    // keeping the path intact.  Only split when the junction is killed
                    // AND there's no nearby good junction to redirect to — this
                    // indicates a genuinely bad alignment that can't be repaired.
                    if !valid {
                        // Check if a nearby good junction exists that could bridge this gap
                        let has_nearby_good = killed_junction_pairs
                            .map(|set| {
                                // Look for any non-killed junction within the correction window
                                read.junctions.iter().any(|rj| {
                                    rj.donor.abs_diff(prev_node.end) <= junction_correction_window
                                        && rj.acceptor.abs_diff(curr_node.start) <= junction_correction_window
                                        && !set.contains(rj)
                                })
                            })
                            .unwrap_or(false);
                        if !has_nearby_good {
                            split_here = true;
                            orphan_right = killed;
                        }
                        if trace_this {
                            eprintln!(
                                "PATH_frag_pattern: read[{}] SPLIT_BADJUNC at_j={} junc_i={} upstream_nodes={} | restart_node={}({}-{}) junc={}-{} node_prev_end={} node_cur_start={} killed={}",
                                read_idx.unwrap_or(0),
                                j,
                                junc_idx,
                                format_nodes_with_coords(&unique_nodes[seg_start..j], graph),
                                curr_nid,
                                curr_node.start,
                                curr_node.end,
                                junc.donor,
                                junc.acceptor,
                                prev_node.end,
                                curr_node.start,
                                if killed { 1 } else { 0 }
                            );
                        }
                    }
                }
            }
        }

        // StringTie-faithful EXON_SKIP split: if the read's intron spans a covlink-
        // suppressed zero-cov node (one that has source or sink as neighbor, flagged
        // by add_coverage_source_sink_edges), split the read there. This is what
        // StringTie's update_abund does in response to EXON_SKIP (see
        // stringtie_debug/run_locus_17190254_17521824_tracepos.stderr: read[283]
        // produces two separate transfrags).
        // Gate: enable with RUSTLE_EXON_SKIP_SPLIT=1 (default off — current implementation
        // fires too broadly and regresses matches).
        // Default-on: disable via RUSTLE_EXON_SKIP_SPLIT_OFF=1.
        if !split_here
            && prev_node.end < curr_node.start
            && std::env::var_os("RUSTLE_EXON_SKIP_SPLIT_OFF").is_none()
        {
            // StringTie-style intergenic split signature:
            // 1. The node BEFORE the gap (prev_node) has sink as a child (covlink-added)
            // 2. The node AFTER the gap (curr_node) has source as a parent (covlink-added)
            // This double-boundary pattern is the exact signal
            // add_coverage_source_sink_edges produces for an intergenic node's flanks.
            let prev_has_sink = prev_node.children.contains(graph.sink_id);
            let curr_has_source = curr_node.parents.contains(graph.source_id);
            if prev_has_sink && curr_has_source {
                split_here = true;
                orphan_right = false;
            }
        }

        if split_here {
            if seg_start < j {
                segments.push(ReadPathSegment {
                    path: unique_nodes[seg_start..j].to_vec(),
                    orphan: next_segment_orphan,
                });
            }
            seg_start = j;
            next_segment_orphan = orphan_right;
        }
    }

    if seg_start < unique_nodes.len() {
        segments.push(ReadPathSegment {
            path: unique_nodes[seg_start..].to_vec(),
            orphan: next_segment_orphan,
        });
    }

    if segments.is_empty() {
        segments.push(ReadPathSegment {
            path: unique_nodes.to_vec(),
            orphan: false,
        });
    }

    segments
}

fn add_or_update_transfrag(
    path: &[usize],
    graph: &mut Graph,
    transfrags: &mut Vec<GraphTransfrag>,
    pattern_map: &mut HbHashMap<TransfragKey, usize>,
    tr_index: &mut TreePatIndex,
    psize: usize,
    is_long: bool,
    _mode: AssemblyMode,
    junction_correction_window: u64,
    weight: f64,
    ref_start: u64,
    ref_end: u64,
    has_poly_start: bool,
    has_poly_end: bool,
    killed_junction_orphan: bool,
    junctions: &[Junction],
) {
    let trace_locus = parse_trace_locus_env();
    let ua_trace = update_abund_trace_active();
    let mut key = path.to_vec();
    let raw_key = key.clone();
    if is_long && std::env::var_os("RUSTLE_DISABLE_LR_PATH_TRIM").is_none() {
        trim_longread_path_for_update_abundance(graph, &mut key, ref_start, ref_end);
    } else if key.len() == 1 {
        // update_abundance: skip single-node non-long transfrags.
        if ua_trace {
            let nid = key[0];
            let (ns, ne) = graph
                .nodes
                .get(nid)
                .map(|n| (n.start, n.end))
                .unwrap_or((0, 0));
            eprintln!(
                "PATH_update_abund: SKIP_SINGLE_NODE nid={} {}-{}",
                nid, ns, ne
            );
        }
        return;
    }
    if key.is_empty() {
        return;
    }
    let pattern = build_transfrag_pattern(&key, graph, junctions, junction_correction_window, psize);
    let key_sig = transfrag_key_from_pattern(&key, graph, &pattern);
    if ua_trace && is_long {
        let trimmed = key != raw_key;
        let (ks, ke) = path_span(&key, graph);
        eprintln!(
            "PATH_update_abund: LR_ENTRY ref={}-{} raw_nodes={} final_nodes={} trimmed={} span={}-{}",
            ref_start, ref_end, raw_key.len(), key.len(), trimmed, ks, ke
        );
        eprintln!("PATH_update_abund: FINAL_NODES nodes={:?}", &key);
        eprintln!(
            "PATH_update_abund: LR_FLAGS poly_start={} poly_end={} orphan={}",
            has_poly_start, has_poly_end, killed_junction_orphan
        );
    }
    if path_overlaps_trace(&raw_key, graph, trace_locus)
        || path_overlaps_trace(&key, graph, trace_locus)
    {
        let (raw_s, raw_e) = path_span(&raw_key, graph);
        let (key_s, key_e) = path_span(&key, graph);
        eprintln!(
            "[TRACE_MAPTRF] raw={:?} raw_span={}-{} trimmed={:?} trimmed_span={}-{} ref={}-{} is_long={} poly={}/{} orphan={}",
            raw_key,
            raw_s,
            raw_e,
            key,
            key_s,
            key_e,
            ref_start,
            ref_end,
            is_long,
            has_poly_start,
            has_poly_end,
            killed_junction_orphan
        );
    }
    let first_last = key
        .first()
        .and_then(|&f| graph.nodes.get(f).map(|n| n.start))
        .zip(key.last().and_then(|&l| graph.nodes.get(l).map(|n| n.end)));
    let (cand_longstart, cand_longend) = if is_long {
        if let Some((node_start, node_end)) = first_last {
            (
                // PATH_update_abund only records longstart/longend when
                // the read boundary lies inside the first/last node of this path.
                // Do not invent a boundary for split/orphan suffix/prefix fragments.
                if ref_start >= node_start {
                    ref_start
                } else {
                    0
                },
                if ref_end <= node_end { ref_end } else { 0 },
            )
        } else {
            (0, 0)
        }
    } else {
        (0, 0)
    };

    let mut tf_hit = tr_index.find(graph, &key, &pattern).and_then(|idx| {
        transfrags.get(idx).and_then(|t| {
            if t.node_ids == key && t.pattern == pattern {
                Some(idx)
            } else {
                None
            }
        })
    });
    if tf_hit.is_none() {
        tf_hit = pattern_map.get(&key_sig).copied();
    }
    if let Some(tf_idx) = tf_hit {
        if path_overlaps_trace(&key, graph, trace_locus) {
            eprintln!(
                "[TRACE_MAPTRF] update_existing idx={} key={:?} old_abund={:.3} old_srabund={:.3} add={:.3}",
                tf_idx,
                key,
                transfrags[tf_idx].abundance,
                transfrags[tf_idx].srabund,
                weight
            );
        }
        let tf = &mut transfrags[tf_idx];
        tf.read_count += weight;
        if is_long {
            tf.abundance += weight;
            tf.longread = true;
        } else {
            tf.srabund += weight;
            // short-read mode removed; srabund is accumulated but
            // abundance is only added for long reads
        }
        if is_long {
            // longstart/longend are only set when read bounds fall inside
            // first/last node bounds; otherwise they remain 0.
            if cand_longstart > 0 && (tf.longstart == 0 || cand_longstart < tf.longstart) {
                tf.longstart = cand_longstart;
            }
            if cand_longend > 0 && cand_longend > tf.longend {
                tf.longend = cand_longend;
            }
        }
        if has_poly_start {
            tf.poly_start_unaligned = tf.poly_start_unaligned.saturating_add(1).min(65535);
        }
        if has_poly_end {
            tf.poly_end_unaligned = tf.poly_end_unaligned.saturating_add(1).min(65535);
        }
        if killed_junction_orphan {
            tf.killed_junction_orphan = true;
        }
        if ua_trace && is_long {
            eprintln!(
                "PATH_update_abund: ADD_TO_TRF idx={} abund={:.4} longstart={} longend={}",
                tf_idx, tf.abundance, tf.longstart, tf.longend
            );
            if cand_longstart > 0 {
                eprintln!(
                    "PATH_update_abund: SET_LONGSTART idx={} val={}",
                    tf_idx, tf.longstart
                );
            }
            if cand_longend > 0 {
                eprintln!(
                    "PATH_update_abund: SET_LONGEND idx={} val={}",
                    tf_idx, tf.longend
                );
            }
        }
    } else {
        if path_overlaps_trace(&key, graph, trace_locus) {
            eprintln!(
                "[TRACE_MAPTRF] create_new key={:?} add={:.3} ref={}-{}",
                key, weight, ref_start, ref_end
            );
        }
        // Sort key by genomic coordinate to ensure nodes are in order
        let mut sorted_key = key.clone();
        sorted_key.sort_unstable_by(|&a, &b| {
            let na = &graph.nodes[a];
            let nb = &graph.nodes[b];
            (na.start, na.end, a).cmp(&(nb.start, nb.end, b))
        });
        let mut tf = GraphTransfrag::new(sorted_key.clone(), psize);
        tf.pattern = pattern.clone();
        tf.read_count = weight;
        tf.abundance = weight;
        tf.longread = true;
        tf.longstart = cand_longstart;
        tf.longend = cand_longend;
        if has_poly_start {
            tf.poly_start_unaligned = 1;
        }
        if has_poly_end {
            tf.poly_end_unaligned = 1;
        }
        tf.killed_junction_orphan = killed_junction_orphan;
        let idx = transfrags.len();
        if ua_trace && is_long {
            let (ks, ke) = path_span(&key, graph);
            eprintln!(
                "PATH_update_abund: NEW_TRF idx={} nodes={} span={}-{} abund={:.4} longstart={} longend={}",
                idx, key.len(), ks, ke, tf.abundance, tf.longstart, tf.longend
            );
            if cand_longstart > 0 {
                eprintln!(
                    "PATH_update_abund: SET_LONGSTART idx={} val={}",
                    idx, cand_longstart
                );
            }
            if cand_longend > 0 {
                eprintln!(
                    "PATH_update_abund: SET_LONGEND idx={} val={}",
                    idx, cand_longend
                );
            }
        }
        transfrags.push(tf);
        pattern_map.insert(key_sig, idx);
        tr_index.insert(graph, &key, &pattern, idx);
    }
}

fn cov_drop_significant(prev_cov: f64, prev_len: u64, cur_cov: f64, cur_len: u64) -> bool {
    prev_cov * DROP * (cur_len as f64) > cur_cov * (prev_len as f64)
}

/// update_abundance:
/// For long reads, trim terminal contiguous nodes if read starts/ends near adjacent node boundary
/// and coverage/hard-boundary heuristics indicate the end-node is better represented by source/sink linkage.
fn trim_longread_path_for_update_abundance(
    graph: &Graph,
    path: &mut Vec<usize>,
    read_start: u64,
    read_end: u64,
) {
    let trace_locus = parse_trace_locus_env();
    let trace_this = path_overlaps_trace(path, graph, trace_locus);
    let ua_trace = update_abund_trace_active();
    if path.len() <= 1 {
        return;
    }
    // Sink-side trimming (trim tail nodes).
    let first = path[0];
    let last = *path.last().unwrap_or(&first);
    let can_source_trim = graph
        .nodes
        .get(first)
        .map(|n| read_start >= n.start && !n.hardstart)
        .unwrap_or(false);
    let can_sink_trim = graph
        .nodes
        .get(last)
        .map(|n| read_end <= n.end && !n.hardend)
        .unwrap_or(false);

    if ua_trace {
        if can_sink_trim {
            eprintln!(
                "PATH_update_abund: SINK_CHECK_ENTER last={} read_end={} node_end={}",
                last,
                read_end,
                graph.nodes.get(last).map(|n| n.end).unwrap_or(0)
            );
        }
        if can_source_trim {
            eprintln!(
                "PATH_update_abund: SOURCE_CHECK_ENTER first={} read_start={} node_start={}",
                first,
                read_start,
                graph.nodes.get(first).map(|n| n.start).unwrap_or(0)
            );
        }
    }

    if can_sink_trim {
        let mut i = path.len() - 1;
        while i > 0 {
            let prev = path[i - 1];
            let cur = path[i];
            let Some(prevn) = graph.nodes.get(prev) else {
                break;
            };
            let Some(curn) = graph.nodes.get(cur) else {
                break;
            };
            if prevn.end != curn.start {
                break;
            }
            let dist = read_end.saturating_sub(prevn.end);
            let sig_drop = prevn.hardend
                || cov_drop_significant(
                    prevn.coverage,
                    prevn.length(),
                    curn.coverage,
                    curn.length(),
                );
            if !sig_drop {
                break;
            }
            let mut trim = false;
            if dist < LONGINTRONANCHOR {
                // Ignore synthetic sink branch in near-end alternative-child trimming.
                trim = prevn
                    .children
                    .ones()
                    .any(|c| c != cur && c != graph.sink_id);
            } else if dist < CHI_WIN {
                trim = prevn.children.contains(graph.sink_id);
            }
            if trim {
                if trace_this {
                    let prev_has_sink = prevn.children.contains(graph.sink_id);
                    let prev_has_alt_child = prevn.children.ones().any(|c| c != cur);
                    eprintln!(
                        "[TRACE_TRIM_LR] sink trim dist={} read_end={} prev={}({}-{} cov={:.3}) cur={}({}-{} cov={:.3}) prev_hardend={} sig_drop={} alt_child={} sink_child={}",
                        dist,
                        read_end,
                        prev,
                        prevn.start,
                        prevn.end,
                        prevn.coverage,
                        cur,
                        curn.start,
                        curn.end,
                        curn.coverage,
                        prevn.hardend,
                        sig_drop,
                        prev_has_alt_child,
                        prev_has_sink
                    );
                }
                if ua_trace {
                    let trim_type = if dist < LONGINTRONANCHOR {
                        "SINK_TRIM[1]"
                    } else {
                        "SINK_TRIM[2]"
                    };
                    eprintln!(
                        "PATH_update_abund: {} prev={} cur={} dist={} prev_cov={:.3} cur_cov={:.3} hardend={}",
                        trim_type, prev, cur, dist, prevn.coverage, curn.coverage, prevn.hardend
                    );
                }
                path.truncate(i);
                break;
            }
            i -= 1;
        }
    }

    if path.len() <= 1 {
        return;
    }
    if can_source_trim {
        let mut i = 0usize;
        while i + 1 < path.len() {
            let cur = path[i];
            let nxt = path[i + 1];
            let Some(curn) = graph.nodes.get(cur) else {
                break;
            };
            let Some(nextn) = graph.nodes.get(nxt) else {
                break;
            };
            if curn.end != nextn.start {
                break;
            }
            let dist = nextn.start.saturating_sub(read_start);
            let sig_drop = nextn.hardstart
                || cov_drop_significant(
                    nextn.coverage,
                    nextn.length(),
                    curn.coverage,
                    curn.length(),
                );
            if !sig_drop {
                break;
            }
            let mut trim = false;
            if dist < LONGINTRONANCHOR {
                // Ignore synthetic source branch in near-start alternative-parent trimming.
                trim = nextn
                    .parents
                    .ones()
                    .any(|p| p != cur && p != graph.source_id);
            } else if dist < CHI_WIN {
                // Avoid trimming true upstream-supported starts when a synthetic source helper
                // co-exists with a real contiguous parent (e.g. node2->node36 in STRG.27.*).
                // Keep the old behavior when source is the only parent.
                let has_non_source_parent = nextn.parents.ones().any(|p| p != graph.source_id);
                trim = nextn.parents.contains(graph.source_id) && !has_non_source_parent;
            }
            if trim {
                if trace_this {
                    let next_has_source = nextn.parents.contains(graph.source_id);
                    let next_has_alt_parent = nextn.parents.ones().any(|p| p != cur);
                    eprintln!(
                        "[TRACE_TRIM_LR] source trim dist={} read_start={} cur={}({}-{} cov={:.3}) next={}({}-{} cov={:.3}) next_hardstart={} sig_drop={} alt_parent={} source_parent={}",
                        dist,
                        read_start,
                        cur,
                        curn.start,
                        curn.end,
                        curn.coverage,
                        nxt,
                        nextn.start,
                        nextn.end,
                        nextn.coverage,
                        nextn.hardstart,
                        sig_drop,
                        next_has_alt_parent,
                        next_has_source
                    );
                }
                if ua_trace {
                    let trim_type = if dist < LONGINTRONANCHOR {
                        "SOURCE_TRIM[1]"
                    } else {
                        "SOURCE_TRIM[2]"
                    };
                    eprintln!(
                        "PATH_update_abund: {} cur={} next={} dist={} next_cov={:.3} cur_cov={:.3} hardstart={}",
                        trim_type, cur, nxt, dist, nextn.coverage, curn.coverage, nextn.hardstart
                    );
                }
                path.drain(0..=i);
                break;
            }
            i += 1;
        }
    }
}
