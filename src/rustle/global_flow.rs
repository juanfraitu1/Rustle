//! Global flow extraction layer.
//!
//! Modules:
//! 1. `compute_global_nodecov`: get_trf_long nodecov precomputation
//! 2. `extract_transcripts_greedy_decompose`: flow-decomposition seeds
//!    (RUSTLE_GREEDY_DECOMPOSE=1) — finds all source-to-sink paths via
//!    iterative bottleneck extraction, bypassing competing-junction misses.

use crate::util::coord::len_half_open;
use crate::graph::{Graph, GraphTransfrag};
use crate::path_extract::Transcript;
use crate::types::RunConfig;

const EPSILON: f64 = crate::util::constants::FLOW_EPSILON;
const CHECKTRF_REDISTRIBUTE_INTRON_TOL: u64 = 5;

// ─── Helper functions ───────────────────────────────────────────────────────

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
///
/// DEAD CODE: the two-splice witness check now lives in path_extract.rs; this
/// local copy is preserved for a potential global-flow-stage rewiring.
#[allow(dead_code)]
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

/// Compute intron chain from exon coordinate pairs (exon.end, next_exon.start).
fn intron_chain_from_exons(exons: &[(u64, u64)]) -> Vec<(u64, u64)> {
    exons.windows(2).map(|w| (w[0].1, w[1].0)).collect()
}

fn intron_chains_equal_tol(a: &[(u64, u64)], b: &[(u64, u64)], tol: u64) -> bool {
    if a.len() != b.len() {
        return false;
    }
    a.iter()
        .zip(b.iter())
        .all(|((d1, a1), (d2, a2))| d1.abs_diff(*d2) <= tol && a1.abs_diff(*a2) <= tol)
}

// ─── Nodecov precomputation ──────────────────────────────────────────────────

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
            if tf.guide && tf.abundance < 1.0 * 0.1 + 1e-9 {
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

        let nc = nodecov[i];
        // Match nodecov.rs line 252: noderate = node.coverage / nodecov.
        // Do NOT divide by length — that would give units of (cov/reads/bp) which
        // are too small to serve as a meaningful per-base coverage multiplier.
        // Default to 1.0 when nc=0 (no non-source-to-sink transfrags seen).
        noderate[i] = if nc > EPSILON {
            let nr = inode.coverage / nc;
            if nr.is_finite() && nr > 0.0 { nr } else { 1.0 }
        } else {
            1.0
        };
    }

    (nodecov, noderate)
}

// ─── Greedy flow decomposition (RUSTLE_GREEDY_DECOMPOSE=1) ──────────────────

/// Build per-edge capacity from long-read transfrag abundances.
fn build_lr_edge_capacities(
    graph: &Graph,
    transfrags: &[GraphTransfrag],
) -> std::collections::HashMap<(usize, usize), f64> {
    let source_id = graph.source_id;
    let sink_id = graph.sink_id;
    let mut cap: std::collections::HashMap<(usize, usize), f64> =
        std::collections::HashMap::new();
    for tf in transfrags {
        if !tf.longread || tf.abundance < EPSILON || tf.node_ids.is_empty() {
            continue;
        }
        // Add edges from the transfrag's internal node sequence.
        for w in tf.node_ids.windows(2) {
            *cap.entry((w[0], w[1])).or_insert(0.0) += tf.abundance;
        }
        // Synthesize source→first and last→sink boundary edges for partial reads
        // that don't themselves span the full source-to-sink path. Without these
        // the BFS can never reach sink from source through partial-read transfrags.
        let first = tf.node_ids[0];
        let last = *tf.node_ids.last().unwrap();
        if first != source_id {
            *cap.entry((source_id, first)).or_insert(0.0) += tf.abundance;
        }
        if last != sink_id {
            *cap.entry((last, sink_id)).or_insert(0.0) += tf.abundance;
        }
    }
    cap
}

/// BFS from source to sink following edges with positive remaining capacity.
/// At each node, children are tried in descending capacity order so the
/// dominant path (highest bottleneck) is found first.
fn bfs_capacity_path(
    source: usize,
    sink: usize,
    edge_cap: &std::collections::HashMap<(usize, usize), f64>,
    graph: &Graph,
) -> Option<Vec<usize>> {
    let n = graph.n_nodes;
    if source >= n || sink >= n {
        return None;
    }
    let mut parent: Vec<Option<usize>> = vec![None; n];
    let mut visited = vec![false; n];
    let mut queue = std::collections::VecDeque::new();

    visited[source] = true;
    queue.push_back(source);

    while let Some(node) = queue.pop_front() {
        if node == sink {
            let mut path = vec![sink];
            let mut cur = sink;
            while let Some(p) = parent[cur] {
                path.push(p);
                cur = p;
            }
            path.reverse();
            return Some(path);
        }
        let Some(gnode) = graph.nodes.get(node) else {
            continue;
        };
        let mut children: Vec<(usize, f64)> = gnode
            .children
            .ones()
            .filter_map(|c| {
                if visited[c] {
                    return None;
                }
                let cap = edge_cap.get(&(node, c)).copied().unwrap_or(0.0);
                if cap > EPSILON { Some((c, cap)) } else { None }
            })
            .collect();
        // Descend by capacity: most-supported path first.
        children.sort_by(|a, b| {
            b.1.partial_cmp(&a.1).unwrap_or(std::cmp::Ordering::Equal)
        });
        for (child, _) in children {
            visited[child] = true;
            parent[child] = Some(node);
            queue.push_back(child);
        }
    }
    None
}

/// Iteratively extract dominant source-to-sink paths by bottleneck.
/// Each extraction subtracts the bottleneck flow from all path edges.
/// Stops when the bottleneck drops below `min_flow` or no path remains.
fn greedy_flow_decompose_paths(
    graph: &Graph,
    transfrags: &[GraphTransfrag],
    min_flow: f64,
) -> Vec<(Vec<usize>, f64)> {
    let mut edge_cap = build_lr_edge_capacities(graph, transfrags);
    let source = graph.source_id;
    let sink = graph.sink_id;
    let mut results: Vec<(Vec<usize>, f64)> = Vec::new();

    loop {
        let path = match bfs_capacity_path(source, sink, &edge_cap, graph) {
            Some(p) => p,
            None => break,
        };

        let bottleneck = path
            .windows(2)
            .map(|w| edge_cap.get(&(w[0], w[1])).copied().unwrap_or(0.0))
            .fold(f64::INFINITY, f64::min);

        if !bottleneck.is_finite() || bottleneck < min_flow {
            break;
        }

        for w in path.windows(2) {
            if let Some(cap) = edge_cap.get_mut(&(w[0], w[1])) {
                *cap -= bottleneck;
            }
        }
        edge_cap.retain(|_, v| *v > EPSILON);

        results.push((path, bottleneck));
        if results.len() > 10_000 {
            break; // safety cap against runaway
        }
    }
    results
}

/// Extract transcripts using greedy flow decomposition.
///
/// Enable via `RUSTLE_GREEDY_DECOMPOSE=1`. Each source-to-sink path from the
/// decomposition becomes a transcript, bypassing the greedy back/fwd trace
/// that misdirects competing-junction minority paths (e.g. STRG.513.3).
pub fn extract_transcripts_greedy_decompose(
    graph: &Graph,
    transfrags: &mut [GraphTransfrag],
    bundle_chrom: &str,
    bundle_strand: char,
    config: &RunConfig,
) -> Vec<Transcript> {
    let debug = std::env::var_os("RUSTLE_GREEDY_DECOMPOSE_DEBUG").is_some();
    let source_id = graph.source_id;
    let sink_id = graph.sink_id;

    let (nodecov, noderate) = compute_global_nodecov(graph, transfrags);
    let mut remaining_nodecov = nodecov;

    let min_flow = config.readthr.max(EPSILON);
    let decomposed = greedy_flow_decompose_paths(graph, transfrags, min_flow);

    if debug {
        eprintln!(
            "[GREEDY_DECOMPOSE] {} paths from {} transfrags (min_flow={:.2})",
            decomposed.len(),
            transfrags.len(),
            min_flow
        );
    }

    let mut out: Vec<Transcript> = Vec::new();
    let mut kept_paths: Vec<(Vec<usize>, f64, bool, usize)> = Vec::new();

    for (path, flow) in &decomposed {
        // Slice to inner nodes (strip source/sink)
        let inner_start = if path.first() == Some(&source_id) { 1 } else { 0 };
        let inner_end = if path.last() == Some(&sink_id) {
            path.len().saturating_sub(1)
        } else {
            path.len()
        };
        if inner_start >= inner_end {
            if debug {
                eprintln!("[GREEDY_DECOMPOSE] flow={:.2} EMPTY_INNER path_len={} skipped", flow, path.len());
            }
            continue;
        }
        let inner_path = &path[inner_start..inner_end];
        if inner_path.is_empty() {
            continue;
        }

        // Per-junction witness check: each splice edge in the path must be
        // individually witnessed by at least one long-read transfrag.
        // We use per-junction (not consecutive-pair) witnessing because flow-
        // decomposed paths are built from aggregated edge capacities across
        // multiple transfrags — a consecutive-pair check would incorrectly
        // reject valid composite paths.
        let mut unwitnessed = false;
        for p in 1..inner_path.len() {
            let u = inner_path[p - 1];
            let v = inner_path[p];
            if !is_splice_between(graph, u, v) {
                continue; // contiguous exon boundary; no witness needed
            }
            let witnessed = transfrags.iter().any(|tf| {
                tf.node_ids.windows(2).any(|w| w[0] == u && w[1] == v)
            });
            if !witnessed {
                unwitnessed = true;
                break;
            }
        }
        if unwitnessed {
            if debug {
                eprintln!("[GREEDY_DECOMPOSE] flow={:.2} UNWITNESSED_JUNCTION skipped", flow);
            }
            continue;
        }

        // Build exons from inner_path. Coverage is derived from nodecov/noderate
        // where available, falling back to the actual node coverage when nodecov=0
        // (which happens for single-exon loci where all reads are source-to-sink
        // and compute_global_nodecov excludes them).
        let mut exons: Vec<(u64, u64)> = Vec::new();
        let mut cov_total = 0.0f64;

        let mut j = 0;
        while j < inner_path.len() {
            let nid = inner_path[j];
            let Some(node) = graph.nodes.get(nid) else {
                j += 1;
                continue;
            };
            let start = node.start;
            let mut end = node.end;

            // Per-node coverage contribution (deplete nodecov if available).
            let ecov0 = {
                let nc = remaining_nodecov.get(nid).copied().unwrap_or(0.0);
                if nc > EPSILON {
                    let nf = flow.min(nc);
                    remaining_nodecov[nid] = (nc - nf).max(0.0);
                    nf * noderate.get(nid).copied().unwrap_or(1.0)
                } else {
                    // nc=0 (e.g. single-exon bundle: all transfrags span source→sink,
                    // excluded from nodecov). Use raw n.coverage = Σ(weight×bp_overlap)
                    // so that coverage = n.coverage / exon_len = reads/bp, matching the
                    // baseline fallback behaviour in path_extract.rs.
                    graph.nodes.get(nid).map(|n| n.coverage).unwrap_or(0.0)
                }
            };
            let mut ecov_accum = ecov0;

            // Merge contiguous nodes into a single exon.
            while j + 1 < inner_path.len()
                && nodes_are_contiguous(graph, inner_path[j], inner_path[j + 1])
            {
                j += 1;
                let nid2 = inner_path[j];
                if let Some(nn) = graph.nodes.get(nid2) {
                    end = nn.end;
                    let nc2 = remaining_nodecov.get(nid2).copied().unwrap_or(0.0);
                    ecov_accum += if nc2 > EPSILON {
                        let nf2 = flow.min(nc2);
                        remaining_nodecov[nid2] = (nc2 - nf2).max(0.0);
                        nf2 * noderate.get(nid2).copied().unwrap_or(1.0)
                    } else {
                        nn.coverage
                    };
                }
            }

            if end > start {
                exons.push((start, end));
                cov_total += ecov_accum;
            }
            j += 1;
        }

        if exons.is_empty() {
            if debug {
                eprintln!("[GREEDY_DECOMPOSE] flow={:.2} NO_EXONS skipped (inner_path_len={})", flow, inner_path.len());
            }
            continue;
        }

        let length: u64 = exons.iter().map(|(s, e)| len_half_open(*s, *e)).sum();
        if length < config.min_transcript_length {
            if debug {
                eprintln!(
                    "[GREEDY_DECOMPOSE] flow={:.2} TOO_SHORT length={} < {} skipped exons={} span={}-{}",
                    flow, length, config.min_transcript_length,
                    exons.len(),
                    exons.first().map(|e| e.0).unwrap_or(0),
                    exons.last().map(|e| e.1).unwrap_or(0),
                );
            }
            continue;
        }

        // Use flow as coverage if cov_total is zero (coverage estimation failed).
        let coverage = if cov_total > EPSILON && length > 0 {
            cov_total / length as f64
        } else {
            *flow
        };
        if coverage < EPSILON {
            if debug {
                eprintln!("[GREEDY_DECOMPOSE] flow={:.2} ZERO_COV cov_total={:.4} length={} skipped", flow, cov_total, length);
            }
            continue;
        }

        if debug {
            eprintln!(
                "[GREEDY_DECOMPOSE] flow={:.2} exons={} span={}-{} cov={:.4}",
                flow,
                exons.len(),
                exons.first().map(|e| e.0).unwrap_or(0),
                exons.last().map(|e| e.1).unwrap_or(0),
                coverage
            );
        }

        let first_node = inner_path[0];
        let last_node = *inner_path.last().unwrap();
        let out_idx = out.len();

        out.push(Transcript {
            chrom: bundle_chrom.to_string(),
            strand: bundle_strand,
            exon_cov: vec![coverage; exons.len()],
            exons: exons.clone(),
            coverage,
            source: Some("greedy_decompose".to_string()),
            is_longread: true,
            longcov: *flow,
            hardstart: graph.nodes.get(first_node).map(|n| n.hardstart).unwrap_or(false),
            hardend: graph.nodes.get(last_node).map(|n| n.hardend).unwrap_or(false),
            ..Transcript::default()
        });

        kept_paths.push((inner_path.to_vec(), coverage, false, out_idx));
    }

    // Zero all long-read transfrag abundances (consumed by decomposition).
    for tf in transfrags.iter_mut() {
        if tf.longread {
            tf.abundance = 0.0;
        }
    }

    if debug {
        eprintln!("[GREEDY_DECOMPOSE] done: {} transcripts", out.len());
    }

    out
}

// ─── Read-chain enumeration (RUSTLE_READCHAIN=1) ────────────────────────────

/// Extract transcripts by grouping long-read transfrags by their exact intron chain.
///
/// Enable via `RUSTLE_READCHAIN=1`. Each unique (donor, acceptor) chain observed in
/// long-read transfrags becomes one transcript, provided total abundance >= readthr.
/// This avoids flow-composition artifacts: every emitted transcript is directly
/// supported by one or more actual reads with matching splice sites.
///
/// Single-exon transfrags (empty intron chain) are skipped; set RUSTLE_READCHAIN_SINGLE
/// to also group single-exon reads by their locus mid-point (coarse clustering).
pub fn extract_transcripts_readchain(
    graph: &Graph,
    transfrags: &mut [GraphTransfrag],
    bundle_chrom: &str,
    bundle_strand: char,
    config: &RunConfig,
) -> Vec<Transcript> {
    let debug = std::env::var_os("RUSTLE_READCHAIN_DEBUG").is_some();
    let source_id = graph.source_id;
    let sink_id = graph.sink_id;
    let min_flow = config.readthr.max(EPSILON);

    // Group long-read transfrags by intron chain.
    // Value: (total_abundance, min_start, max_end, representative_node_ids)
    let mut chain_groups: std::collections::HashMap<
        Vec<(u64, u64)>,
        (f64, u64, u64, Vec<usize>),
    > = std::collections::HashMap::new();

    for tf in transfrags.iter() {
        if !tf.longread || tf.abundance < EPSILON || tf.node_ids.is_empty() {
            continue;
        }
        // Strip virtual source/sink nodes.
        let inner_start = if tf.node_ids.first() == Some(&source_id) { 1 } else { 0 };
        let inner_end = if tf.node_ids.last() == Some(&sink_id) {
            tf.node_ids.len().saturating_sub(1)
        } else {
            tf.node_ids.len()
        };
        if inner_start >= inner_end {
            continue;
        }
        let inner = &tf.node_ids[inner_start..inner_end];
        if inner.is_empty() {
            continue;
        }

        let chain = intron_chain_from_nodes(graph, inner);

        // Skip single-exon reads — baseline handles those.
        if chain.is_empty() {
            continue;
        }

        let first_node = inner[0];
        let last_node = *inner.last().unwrap();
        let tf_start = graph.nodes.get(first_node).map(|n| n.start).unwrap_or(0);
        let tf_end = graph.nodes.get(last_node).map(|n| n.end).unwrap_or(0);

        let entry = chain_groups
            .entry(chain)
            .or_insert((0.0, tf_start, tf_end, inner.to_vec()));
        entry.0 += tf.abundance;
        if tf_start < entry.1 {
            entry.1 = tf_start;
        }
        if tf_end > entry.2 {
            entry.2 = tf_end;
        }
        // Keep highest-abundance transfrag's nodes as representative.
        // (entry is initialised with the first seen; no update needed for coverage
        //  because we use nodecov rather than raw node coords for coverage.)
    }

    // Isofrac filter: within this bundle, drop chains below isofrac fraction of
    // the most-abundant chain. This mirrors how extract_transcripts filters minority
    // isoforms per locus and is the main guard against FPs from rare mis-spliced reads.
    let max_abund_in_bundle = chain_groups
        .values()
        .map(|(a, ..)| *a)
        .fold(0.0_f64, f64::max);
    let isofrac_threshold = if max_abund_in_bundle > EPSILON && config.transcript_isofrac > 0.0 {
        (max_abund_in_bundle * config.transcript_isofrac).max(min_flow)
    } else {
        min_flow
    };

    let (nodecov, noderate) = compute_global_nodecov(graph, transfrags);
    let mut remaining_nodecov = nodecov;

    let mut out: Vec<Transcript> = Vec::new();

    // Sort by descending abundance for deterministic nodecov depletion order.
    let mut sorted_chains: Vec<(&Vec<(u64, u64)>, &(f64, u64, u64, Vec<usize>))> =
        chain_groups.iter().collect();
    sorted_chains.sort_by(|a, b| {
        b.1 .0
            .partial_cmp(&a.1 .0)
            .unwrap_or(std::cmp::Ordering::Equal)
    });

    for (chain, (total_abund, min_start, max_end, rep_nodes)) in sorted_chains {
        if *total_abund < isofrac_threshold {
            if debug {
                eprintln!(
                    "[READCHAIN] chain_len={} abund={:.2} < threshold={:.2} FILTERED",
                    chain.len(), total_abund, isofrac_threshold
                );
            }
            continue;
        }

        // Build exon list from intron chain + span.
        // Internal exon boundaries are fully determined by (acceptor_{i-1}, donor_i);
        // only the first-exon start and last-exon end vary across reads.
        let mut exons: Vec<(u64, u64)> = Vec::new();
        exons.push((*min_start, chain[0].0)); // first exon
        for i in 1..chain.len() {
            exons.push((chain[i - 1].1, chain[i].0)); // internal exons
        }
        exons.push((chain.last().unwrap().1, *max_end)); // last exon

        // Validate: all exons must have positive length.
        if exons.iter().any(|(s, e)| e <= s) {
            if debug {
                eprintln!("[READCHAIN] INVALID_EXONS chain_len={} skipped", chain.len());
            }
            continue;
        }

        let length: u64 = exons.iter().map(|(s, e)| len_half_open(*s, *e)).sum();
        if length < config.min_transcript_length {
            if debug {
                eprintln!(
                    "[READCHAIN] TOO_SHORT length={} < {} skipped",
                    length, config.min_transcript_length
                );
            }
            continue;
        }

        // Coverage from nodecov/noderate of representative nodes.
        let mut cov_total = 0.0f64;
        let mut j = 0usize;
        while j < rep_nodes.len() {
            let nid = rep_nodes[j];
            let nc = remaining_nodecov.get(nid).copied().unwrap_or(0.0);
            let mut ecov_accum = if nc > EPSILON {
                let nf = total_abund.min(nc);
                remaining_nodecov[nid] = (nc - nf).max(0.0);
                nf * noderate.get(nid).copied().unwrap_or(1.0)
            } else {
                graph.nodes.get(nid).map(|n| n.coverage).unwrap_or(0.0)
            };
            while j + 1 < rep_nodes.len()
                && nodes_are_contiguous(graph, rep_nodes[j], rep_nodes[j + 1])
            {
                j += 1;
                let nid2 = rep_nodes[j];
                let nc2 = remaining_nodecov.get(nid2).copied().unwrap_or(0.0);
                ecov_accum += if nc2 > EPSILON {
                    let nf2 = total_abund.min(nc2);
                    remaining_nodecov[nid2] = (nc2 - nf2).max(0.0);
                    nf2 * noderate.get(nid2).copied().unwrap_or(1.0)
                } else {
                    graph.nodes.get(nid2).map(|n| n.coverage).unwrap_or(0.0)
                };
            }
            cov_total += ecov_accum;
            j += 1;
        }

        let coverage = if cov_total > EPSILON && length > 0 {
            cov_total / length as f64
        } else {
            *total_abund
        };
        if coverage < EPSILON {
            continue;
        }

        if debug {
            eprintln!(
                "[READCHAIN] chain_len={} abund={:.2} exons={} span={}-{} cov={:.4}",
                chain.len(), total_abund, exons.len(),
                exons.first().map(|e| e.0).unwrap_or(0),
                exons.last().map(|e| e.1).unwrap_or(0),
                coverage
            );
        }

        let first_inner = rep_nodes[0];
        let last_inner = *rep_nodes.last().unwrap();

        out.push(Transcript {
            chrom: bundle_chrom.to_string(),
            strand: bundle_strand,
            exon_cov: vec![coverage; exons.len()],
            exons,
            coverage,
            source: Some("readchain".to_string()),
            is_longread: true,
            longcov: *total_abund,
            hardstart: graph.nodes.get(first_inner).map(|n| n.hardstart).unwrap_or(false),
            hardend: graph.nodes.get(last_inner).map(|n| n.hardend).unwrap_or(false),
            ..Transcript::default()
        });
    }

    // Zero all long-read transfrag abundances (consumed by readchain extraction).
    for tf in transfrags.iter_mut() {
        if tf.longread {
            tf.abundance = 0.0;
        }
    }

    if debug {
        eprintln!(
            "[READCHAIN] done: {} transcripts from {} chains",
            out.len(),
            chain_groups.len()
        );
    }

    out
}

// ─── Read-chain supplement (RUSTLE_READCHAIN_SUPPLEMENT=1) ──────────────────

/// Supplement baseline extraction with read-chain paths that have novel intron chains.
///
/// Enable via `RUSTLE_READCHAIN_SUPPLEMENT=1`. Runs after baseline extraction on the
/// RESTORED transfrag abundances. Groups long-read transfrags by intron chain, then emits
/// only chains with at least one intron absent from any baseline transcript.
/// Single-exon chains are always skipped.
pub fn readchain_supplement(
    graph: &Graph,
    transfrags: &[GraphTransfrag],
    baseline_txs: &[Transcript],
    bundle_chrom: &str,
    bundle_strand: char,
    config: &RunConfig,
) -> Vec<Transcript> {
    let debug = std::env::var_os("RUSTLE_READCHAIN_SUPPLEMENT_DEBUG").is_some();
    let source_id = graph.source_id;
    let sink_id = graph.sink_id;
    // Default min_flow for supplement is higher than readthr to avoid FPs from single reads.
    let min_flow = std::env::var("RUSTLE_READCHAIN_SUPPLEMENT_MIN_FLOW")
        .ok()
        .and_then(|v| v.parse::<f64>().ok())
        .unwrap_or(2.0_f64.max(config.readthr))
        .max(EPSILON);
    let min_introns = std::env::var("RUSTLE_READCHAIN_SUPPLEMENT_MIN_INTRONS")
        .ok()
        .and_then(|v| v.parse::<usize>().ok())
        .unwrap_or(2); // require ≥ 3 exons (2 introns) by default

    // Build per-transcript intron chains for exact-chain deduplication.
    let baseline_chains: Vec<Vec<(u64, u64)>> = baseline_txs
        .iter()
        .map(|tx| intron_chain_from_exons(&tx.exons))
        .collect();

    // Group long-read transfrags by intron chain.
    let mut chain_groups: std::collections::HashMap<
        Vec<(u64, u64)>,
        (f64, u64, u64, Vec<usize>),
    > = std::collections::HashMap::new();

    for tf in transfrags.iter() {
        if !tf.longread || tf.abundance < EPSILON || tf.node_ids.is_empty() {
            continue;
        }
        let inner_start = if tf.node_ids.first() == Some(&source_id) { 1 } else { 0 };
        let inner_end = if tf.node_ids.last() == Some(&sink_id) {
            tf.node_ids.len().saturating_sub(1)
        } else {
            tf.node_ids.len()
        };
        if inner_start >= inner_end {
            continue;
        }
        let inner = &tf.node_ids[inner_start..inner_end];
        if inner.is_empty() {
            continue;
        }
        let chain = intron_chain_from_nodes(graph, inner);
        if chain.is_empty() {
            continue; // skip single-exon
        }
        let first_node = inner[0];
        let last_node = *inner.last().unwrap();
        let tf_start = graph.nodes.get(first_node).map(|n| n.start).unwrap_or(0);
        let tf_end = graph.nodes.get(last_node).map(|n| n.end).unwrap_or(0);
        let entry = chain_groups
            .entry(chain)
            .or_insert((0.0, tf_start, tf_end, inner.to_vec()));
        entry.0 += tf.abundance;
        if tf_start < entry.1 { entry.1 = tf_start; }
        if tf_end > entry.2 { entry.2 = tf_end; }
    }

    let (nodecov, noderate) = compute_global_nodecov(graph, transfrags);
    let mut remaining_nodecov = nodecov;
    let mut out: Vec<Transcript> = Vec::new();

    // Sort by descending abundance for deterministic nodecov depletion order.
    let mut sorted_chains: Vec<(&Vec<(u64, u64)>, &(f64, u64, u64, Vec<usize>))> =
        chain_groups.iter().collect();
    sorted_chains.sort_by(|a, b| {
        b.1 .0
            .partial_cmp(&a.1 .0)
            .unwrap_or(std::cmp::Ordering::Equal)
    });

    for (chain, (total_abund, min_start, max_end, rep_nodes)) in sorted_chains {
        if *total_abund < min_flow {
            continue;
        }

        // Require minimum intron count to avoid trivially short chains.
        if chain.len() < min_introns {
            continue;
        }

        // Skip if baseline already has a transcript with this exact intron chain.
        // We check for exact-chain match (within tolerance), NOT sub-chain containment.
        // This correctly emits "shorter" isoforms whose complete chain isn't in any
        // baseline transcript, even if their individual introns appear in longer chains.
        let chain_already_in_baseline = baseline_chains.iter().any(|bc| {
            intron_chains_equal_tol(bc, chain, CHECKTRF_REDISTRIBUTE_INTRON_TOL)
        });
        if chain_already_in_baseline {
            if debug {
                eprintln!("[RC_SUPPL] chain_len={} abund={:.2} CHAIN_IN_BASELINE skipped", chain.len(), total_abund);
            }
            continue;
        }

        // Build exon coordinates from intron chain + span.
        let mut exons: Vec<(u64, u64)> = Vec::new();
        exons.push((*min_start, chain[0].0));
        for i in 1..chain.len() {
            exons.push((chain[i - 1].1, chain[i].0));
        }
        exons.push((chain.last().unwrap().1, *max_end));

        if exons.iter().any(|(s, e)| e <= s) {
            continue;
        }
        let length: u64 = exons.iter().map(|(s, e)| len_half_open(*s, *e)).sum();
        if length < config.min_transcript_length {
            continue;
        }

        // Coverage from nodecov.
        let mut cov_total = 0.0f64;
        let mut j = 0usize;
        while j < rep_nodes.len() {
            let nid = rep_nodes[j];
            let nc = remaining_nodecov.get(nid).copied().unwrap_or(0.0);
            let mut ecov_accum = if nc > EPSILON {
                let nf = total_abund.min(nc);
                remaining_nodecov[nid] = (nc - nf).max(0.0);
                nf * noderate.get(nid).copied().unwrap_or(1.0)
            } else {
                graph.nodes.get(nid).map(|n| n.coverage).unwrap_or(0.0)
            };
            while j + 1 < rep_nodes.len()
                && nodes_are_contiguous(graph, rep_nodes[j], rep_nodes[j + 1])
            {
                j += 1;
                let nid2 = rep_nodes[j];
                let nc2 = remaining_nodecov.get(nid2).copied().unwrap_or(0.0);
                ecov_accum += if nc2 > EPSILON {
                    let nf2 = total_abund.min(nc2);
                    remaining_nodecov[nid2] = (nc2 - nf2).max(0.0);
                    nf2 * noderate.get(nid2).copied().unwrap_or(1.0)
                } else {
                    graph.nodes.get(nid2).map(|n| n.coverage).unwrap_or(0.0)
                };
            }
            cov_total += ecov_accum;
            j += 1;
        }
        let coverage = if cov_total > EPSILON && length > 0 {
            cov_total / length as f64
        } else {
            *total_abund
        };
        if coverage < EPSILON {
            continue;
        }

        if debug {
            eprintln!(
                "[RC_SUPPL] chain_len={} abund={:.2} exons={} span={}-{}",
                chain.len(), total_abund, exons.len(),
                exons.first().map(|e| e.0).unwrap_or(0),
                exons.last().map(|e| e.1).unwrap_or(0),
            );
        }

        let first_inner = rep_nodes[0];
        let last_inner = *rep_nodes.last().unwrap();
        out.push(Transcript {
            chrom: bundle_chrom.to_string(),
            strand: bundle_strand,
            exon_cov: vec![coverage; exons.len()],
            exons,
            coverage,
            source: Some("readchain_supplement".to_string()),
            is_longread: true,
            longcov: *total_abund,
            hardstart: graph.nodes.get(first_inner).map(|n| n.hardstart).unwrap_or(false),
            hardend: graph.nodes.get(last_inner).map(|n| n.hardend).unwrap_or(false),
            ..Transcript::default()
        });
    }

    if debug {
        eprintln!("[RC_SUPPL] done: {} novel chains emitted", out.len());
    }
    out
}

// ─── Greedy supplement (RUSTLE_GREEDY_SUPPLEMENT=1) ─────────────────────────

/// Supplement baseline extraction with greedy-flow paths that have novel intron chains.
///
/// Enable via `RUSTLE_GREEDY_SUPPLEMENT=1`. Runs greedy flow decomposition on the
/// *restored* transfrags (caller must restore abundances before calling), then emits
/// only paths whose intron chain is not already present in `baseline_txs`. Single-exon
/// paths are always skipped — the baseline handles those. All splice edges must be
/// individually witnessed by at least one long-read transfrag.
pub fn greedy_supplement(
    graph: &Graph,
    transfrags: &[GraphTransfrag],
    baseline_txs: &[Transcript],
    bundle_chrom: &str,
    bundle_strand: char,
    config: &RunConfig,
) -> Vec<Transcript> {
    let debug = std::env::var_os("RUSTLE_GREEDY_SUPPLEMENT_DEBUG").is_some();
    let source_id = graph.source_id;
    let sink_id = graph.sink_id;

    // Build (a) flat set of all baseline introns and (b) per-transcript chains.
    // Two complementary guards protect baseline transcripts from supplement:
    //   - "novel intron" check: supplement must have ≥1 intron absent from baseline
    //   - "superset" check: supplement must NOT contain any baseline chain as a
    //     contiguous sub-sequence (prevents "extended" versions of baseline transcripts
    //     from displacing them in pairwise filtering without improving reference coverage)
    let baseline_introns: Vec<(u64, u64)> = baseline_txs
        .iter()
        .flat_map(|tx| intron_chain_from_exons(&tx.exons))
        .collect();
    let baseline_chains: Vec<Vec<(u64, u64)>> = baseline_txs
        .iter()
        .map(|tx| intron_chain_from_exons(&tx.exons))
        .collect();

    // RUSTLE_GREEDY_SUPPLEMENT_MIN_FLOW overrides the minimum flow for supplement paths.
    // Default: config.readthr (same as baseline). Raising this (e.g. to 3–5) filters
    // low-confidence paths and improves precision at the cost of fewer novel discoveries.
    let min_flow = std::env::var("RUSTLE_GREEDY_SUPPLEMENT_MIN_FLOW")
        .ok()
        .and_then(|v| v.parse::<f64>().ok())
        .unwrap_or(config.readthr)
        .max(EPSILON);
    let decomposed = greedy_flow_decompose_paths(graph, transfrags, min_flow);

    if debug {
        eprintln!(
            "[GREEDY_SUPPLEMENT] {} greedy paths, {} baseline introns, min_flow={:.2}",
            decomposed.len(),
            baseline_introns.len(),
            min_flow,
        );
    }

    let (nodecov, noderate) = compute_global_nodecov(graph, transfrags);
    let mut remaining_nodecov = nodecov;

    let mut out: Vec<Transcript> = Vec::new();

    for (path, flow) in &decomposed {
        // Strip source/sink virtual nodes.
        let inner_start = if path.first() == Some(&source_id) { 1 } else { 0 };
        let inner_end = if path.last() == Some(&sink_id) {
            path.len().saturating_sub(1)
        } else {
            path.len()
        };
        if inner_start >= inner_end {
            continue;
        }
        let inner_path = &path[inner_start..inner_end];
        if inner_path.is_empty() {
            continue;
        }

        // Compute intron chain for this greedy path.
        let chain = intron_chain_from_nodes(graph, inner_path);

        // Require at least one intron not seen in any baseline transcript.
        // This filters out: (a) single-exon paths (chain is empty → no novel intron),
        // (b) exact-chain duplicates, and (c) sub-path artifacts whose junctions are a
        // subset of baseline junctions — while admitting competing-junction paths that
        // have a genuinely new (donor, acceptor) pair.
        let has_novel_intron = chain.iter().any(|&(d, a)| {
            !baseline_introns.iter().any(|&(bd, ba)| {
                d.abs_diff(bd) <= CHECKTRF_REDISTRIBUTE_INTRON_TOL
                    && a.abs_diff(ba) <= CHECKTRF_REDISTRIBUTE_INTRON_TOL
            })
        });
        if !has_novel_intron {
            if debug {
                eprintln!(
                    "[GREEDY_SUPPLEMENT] flow={:.2} chain_len={} NO_NOVEL_INTRON skipped",
                    flow,
                    chain.len()
                );
            }
            continue;
        }

        // Skip if this chain is a strict superset of any baseline chain: the supplement
        // transcript is just a longer version of an existing baseline transcript and would
        // kill it in pairwise filtering without providing a better reference match.
        let is_superset_of_baseline = baseline_chains.iter().any(|bc| {
            if bc.is_empty() || chain.len() <= bc.len() {
                return false;
            }
            (0..=(chain.len() - bc.len())).any(|start| {
                bc.iter().enumerate().all(|(i, &(bd, ba))| {
                    let (sd, sa) = chain[start + i];
                    sd.abs_diff(bd) <= CHECKTRF_REDISTRIBUTE_INTRON_TOL
                        && sa.abs_diff(ba) <= CHECKTRF_REDISTRIBUTE_INTRON_TOL
                })
            })
        });
        if is_superset_of_baseline {
            if debug {
                eprintln!(
                    "[GREEDY_SUPPLEMENT] flow={:.2} chain_len={} SUPERSET_OF_BASELINE skipped",
                    flow,
                    chain.len()
                );
            }
            continue;
        }

        // Per-junction witness check: every splice edge must appear in at least one transfrag.
        let mut unwitnessed = false;
        for p in 1..inner_path.len() {
            let u = inner_path[p - 1];
            let v = inner_path[p];
            if !is_splice_between(graph, u, v) {
                continue;
            }
            let witnessed = transfrags
                .iter()
                .any(|tf| tf.node_ids.windows(2).any(|w| w[0] == u && w[1] == v));
            if !witnessed {
                unwitnessed = true;
                break;
            }
        }
        if unwitnessed {
            if debug {
                eprintln!(
                    "[GREEDY_SUPPLEMENT] flow={:.2} UNWITNESSED_JUNCTION skipped",
                    flow
                );
            }
            continue;
        }

        // Build exons and accumulate coverage using same nodecov/noderate approach.
        let mut exons: Vec<(u64, u64)> = Vec::new();
        let mut cov_total = 0.0f64;

        let mut j = 0;
        while j < inner_path.len() {
            let nid = inner_path[j];
            let Some(node) = graph.nodes.get(nid) else {
                j += 1;
                continue;
            };
            let start = node.start;
            let mut end = node.end;

            let ecov0 = {
                let nc = remaining_nodecov.get(nid).copied().unwrap_or(0.0);
                if nc > EPSILON {
                    let nf = flow.min(nc);
                    remaining_nodecov[nid] = (nc - nf).max(0.0);
                    nf * noderate.get(nid).copied().unwrap_or(1.0)
                } else {
                    graph.nodes.get(nid).map(|n| n.coverage).unwrap_or(0.0)
                }
            };
            let mut ecov_accum = ecov0;

            // Merge contiguous nodes into a single exon.
            while j + 1 < inner_path.len()
                && nodes_are_contiguous(graph, inner_path[j], inner_path[j + 1])
            {
                j += 1;
                let nid2 = inner_path[j];
                if let Some(nn) = graph.nodes.get(nid2) {
                    end = nn.end;
                    let nc2 = remaining_nodecov.get(nid2).copied().unwrap_or(0.0);
                    ecov_accum += if nc2 > EPSILON {
                        let nf2 = flow.min(nc2);
                        remaining_nodecov[nid2] = (nc2 - nf2).max(0.0);
                        nf2 * noderate.get(nid2).copied().unwrap_or(1.0)
                    } else {
                        nn.coverage
                    };
                }
            }

            if end > start {
                exons.push((start, end));
                cov_total += ecov_accum;
            }
            j += 1;
        }

        if exons.is_empty() {
            continue;
        }

        let length: u64 = exons.iter().map(|(s, e)| len_half_open(*s, *e)).sum();
        if length < config.min_transcript_length {
            if debug {
                eprintln!(
                    "[GREEDY_SUPPLEMENT] flow={:.2} TOO_SHORT length={} < {}",
                    flow, length, config.min_transcript_length
                );
            }
            continue;
        }

        let coverage = if cov_total > EPSILON && length > 0 {
            cov_total / length as f64
        } else {
            *flow
        };
        if coverage < EPSILON {
            continue;
        }

        if debug {
            eprintln!(
                "[GREEDY_SUPPLEMENT] NEW chain_len={} flow={:.2} exons={} span={}-{} cov={:.4}",
                chain.len(),
                flow,
                exons.len(),
                exons.first().map(|e| e.0).unwrap_or(0),
                exons.last().map(|e| e.1).unwrap_or(0),
                coverage,
            );
        }

        let first_node = inner_path[0];
        let last_node = *inner_path.last().unwrap();
        out.push(Transcript {
            chrom: bundle_chrom.to_string(),
            strand: bundle_strand,
            exon_cov: vec![coverage; exons.len()],
            exons,
            coverage,
            source: Some("greedy_supplement".to_string()),
            is_longread: true,
            longcov: *flow,
            hardstart: graph.nodes.get(first_node).map(|n| n.hardstart).unwrap_or(false),
            hardend: graph.nodes.get(last_node).map(|n| n.hardend).unwrap_or(false),
            ..Transcript::default()
        });
    }

    if debug {
        eprintln!("[GREEDY_SUPPLEMENT] done: {} new transcripts", out.len());
    }

    out
}

// ─── Donor-level competing-junction supplement (RUSTLE_COMPETING_JX_SUPPLEMENT) ──

/// For each splice donor with ≥2 adequately-witnessed acceptors, emit a transcript
/// for every minority acceptor whose intron is not already covered by baseline.
///
/// Unlike the global greedy supplement, this operates at individual junction-donor
/// nodes rather than flow-decomposing the full graph.  It targets exactly the case
/// where the greedy back-trace always picks the high-flow acceptor, leaving the
/// low-flow (but witness-supported) acceptor unextracted.
///
/// Enable via `RUSTLE_COMPETING_JX_SUPPLEMENT=1`.
/// Minimum witness reads for a minority acceptor: `RUSTLE_COMPETING_JX_MIN_WITNESSES`
/// (default 2).  The majority acceptor must have strictly more witness reads than
/// the minority acceptor for the donor to count as "competing".
pub fn competing_junction_supplement(
    graph: &Graph,
    transfrags: &[GraphTransfrag],
    baseline_txs: &[Transcript],
    bundle_chrom: &str,
    bundle_strand: char,
    config: &RunConfig,
) -> Vec<Transcript> {
    let debug = std::env::var_os("RUSTLE_COMPETING_JX_DEBUG").is_some();
    let source_id = graph.source_id;
    let sink_id = graph.sink_id;

    let min_witnesses: usize = std::env::var("RUSTLE_COMPETING_JX_MIN_WITNESSES")
        .ok()
        .and_then(|v| v.parse().ok())
        .unwrap_or(3);

    // Flat baseline intron set (individual introns) for per-junction novelty check
    let baseline_introns: Vec<(u64, u64)> = baseline_txs
        .iter()
        .flat_map(|tx| intron_chain_from_exons(&tx.exons))
        .collect();

    // Per-transcript intron chains for full-chain dedup (upstream extension check)
    let baseline_chains: Vec<Vec<(u64, u64)>> = baseline_txs
        .iter()
        .map(|tx| intron_chain_from_exons(&tx.exons))
        .collect();

    // Count per-splice-edge witness transfrags and total abundance
    let mut edge_tfs: std::collections::HashMap<(usize, usize), Vec<usize>> =
        std::collections::HashMap::new();
    for (ti, tf) in transfrags.iter().enumerate() {
        if !tf.longread || tf.abundance < EPSILON {
            continue;
        }
        for w in tf.node_ids.windows(2) {
            let (u, v) = (w[0], w[1]);
            if is_splice_between(graph, u, v) {
                edge_tfs.entry((u, v)).or_default().push(ti);
            }
        }
    }

    // Build donor → [(acceptor, witness_count)] map; only edges with ≥1 witness
    let mut donor_acceptors: std::collections::HashMap<usize, Vec<(usize, usize)>> =
        std::collections::HashMap::new();
    for (&(u, v), tfs) in &edge_tfs {
        donor_acceptors.entry(u).or_default().push((v, tfs.len()));
    }

    let (nodecov, noderate) = compute_global_nodecov(graph, transfrags);
    let mut remaining_nodecov = nodecov;
    let mut out: Vec<Transcript> = Vec::new();

    for (donor, mut acceptors) in donor_acceptors {
        // Need ≥2 acceptors each with ≥min_witnesses reads to be a competing donor
        acceptors.retain(|(_, cnt)| *cnt >= min_witnesses);
        if acceptors.len() < 2 {
            continue;
        }

        // Sort descending by witness count so acceptors[0] is dominant
        acceptors.sort_by(|a, b| b.1.cmp(&a.1));

        let dom_count = acceptors[0].1;

        // Process every minority acceptor (rank ≥ 1)
        for &(acceptor, min_count) in &acceptors[1..] {
            // Require dominant has strictly more witnesses (true competing junction)
            if min_count >= dom_count {
                continue;
            }

            let donor_coord = match graph.nodes.get(donor) {
                Some(n) => n.end,
                None => continue,
            };
            let acceptor_coord = match graph.nodes.get(acceptor) {
                Some(n) => n.start,
                None => continue,
            };

            // Skip if this specific intron is already in baseline
            let already_in_baseline = baseline_introns.iter().any(|&(bd, ba)| {
                donor_coord.abs_diff(bd) <= CHECKTRF_REDISTRIBUTE_INTRON_TOL
                    && acceptor_coord.abs_diff(ba) <= CHECKTRF_REDISTRIBUTE_INTRON_TOL
            });
            if already_in_baseline {
                if debug {
                    eprintln!(
                        "[COMPETING_JX] donor={} acceptor={} intron={}-{} ALREADY_IN_BASELINE",
                        donor, acceptor, donor_coord, acceptor_coord
                    );
                }
                continue;
            }

            // Find the most-abundant transfrag through this junction
            let best_tf_idx = edge_tfs
                .get(&(donor, acceptor))
                .and_then(|tfs| {
                    tfs.iter().copied().max_by(|&a, &b| {
                        transfrags[a]
                            .abundance
                            .partial_cmp(&transfrags[b].abundance)
                            .unwrap_or(std::cmp::Ordering::Equal)
                    })
                });
            let Some(tf_idx) = best_tf_idx else { continue };
            let tf = &transfrags[tf_idx];

            // Extract inner path (strip source/sink bookends)
            let inner_path: Vec<usize> = tf
                .node_ids
                .iter()
                .copied()
                .skip_while(|&n| n == source_id)
                .take_while(|&n| n != sink_id)
                .collect();

            if inner_path.is_empty() {
                continue;
            }

            let chain = intron_chain_from_nodes(graph, &inner_path);

            // Full-chain novelty check: skip if already in baseline
            let chain_already_covered = chain.iter().all(|&(d, a)| {
                baseline_introns.iter().any(|&(bd, ba)| {
                    d.abs_diff(bd) <= CHECKTRF_REDISTRIBUTE_INTRON_TOL
                        && a.abs_diff(ba) <= CHECKTRF_REDISTRIBUTE_INTRON_TOL
                })
            });
            if !chain.is_empty() && chain_already_covered {
                continue;
            }

            // Build exons by merging contiguous nodes
            let flow = tf.abundance;
            let mut exons: Vec<(u64, u64)> = Vec::new();
            let mut cov_total = 0.0f64;
            let mut j = 0;
            while j < inner_path.len() {
                let nid = inner_path[j];
                let Some(node) = graph.nodes.get(nid) else {
                    j += 1;
                    continue;
                };
                let start = node.start;
                let mut end = node.end;

                let ecov0 = {
                    let nc = remaining_nodecov.get(nid).copied().unwrap_or(0.0);
                    if nc > EPSILON {
                        let nf = flow.min(nc);
                        remaining_nodecov[nid] = (nc - nf).max(0.0);
                        nf * noderate.get(nid).copied().unwrap_or(1.0)
                    } else {
                        node.coverage
                    }
                };
                let mut ecov_accum = ecov0;

                while j + 1 < inner_path.len()
                    && nodes_are_contiguous(graph, inner_path[j], inner_path[j + 1])
                {
                    j += 1;
                    let nid2 = inner_path[j];
                    if let Some(nn) = graph.nodes.get(nid2) {
                        end = nn.end;
                        let nc2 = remaining_nodecov.get(nid2).copied().unwrap_or(0.0);
                        ecov_accum += if nc2 > EPSILON {
                            let nf2 = flow.min(nc2);
                            remaining_nodecov[nid2] = (nc2 - nf2).max(0.0);
                            nf2 * noderate.get(nid2).copied().unwrap_or(1.0)
                        } else {
                            nn.coverage
                        };
                    }
                }

                if end > start {
                    exons.push((start, end));
                    cov_total += ecov_accum;
                }
                j += 1;
            }

            if exons.is_empty() {
                continue;
            }

            let length: u64 = exons.iter().map(|(s, e)| len_half_open(*s, *e)).sum();
            if length < config.min_transcript_length {
                continue;
            }

            let coverage = if cov_total > EPSILON && length > 0 {
                cov_total / length as f64
            } else {
                flow
            };
            if coverage < EPSILON {
                continue;
            }

            if debug {
                eprintln!(
                    "[COMPETING_JX] NEW donor={} acceptor={} intron={}-{} witnesses={} flow={:.1} exons={} span={}-{}",
                    donor, acceptor, donor_coord, acceptor_coord,
                    min_count, flow, exons.len(),
                    exons.first().map(|e| e.0).unwrap_or(0),
                    exons.last().map(|e| e.1).unwrap_or(0),
                );
            }

            let first_node = inner_path[0];
            let last_node = *inner_path.last().unwrap();
            out.push(Transcript {
                chrom: bundle_chrom.to_string(),
                strand: bundle_strand,
                exon_cov: vec![coverage; exons.len()],
                exons,
                coverage,
                source: Some("competing_jx".to_string()),
                is_longread: true,
                longcov: flow,
                hardstart: graph.nodes.get(first_node).map(|n| n.hardstart).unwrap_or(false),
                hardend: graph.nodes.get(last_node).map(|n| n.hardend).unwrap_or(false),
                ..Transcript::default()
            });
        }
    }

    if debug {
        eprintln!("[COMPETING_JX] done: {} new transcripts", out.len());
    }

    out
}
