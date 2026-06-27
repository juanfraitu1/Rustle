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

/// Build a fuzzy-junction canonicalization map for read-chain collapse.
///
/// Alignment wobble places the same true splice site at slightly different
/// (donor, acceptor) coordinates across reads, which an exact intron-chain key
/// splits into near-duplicate isoforms. This clusters junctions whose donor AND
/// acceptor are both within `fuzz` bp of a higher-support junction and maps each
/// to that canonical (consensus) junction — the IsoSeq `--max-fuzzy-junction`
/// idea. `fuzz == 0` returns the identity map (exact behaviour preserved).
///
/// Greedy by descending support (deterministic, coord tie-break): the most-
/// supported junctions become canonicals first, lower-support neighbours snap to
/// them. Chain length is preserved (1:1 remap), so only same-length chains whose
/// every junction snaps together can merge.
fn build_fuzzy_junction_map(
    junctions_with_support: &[((u64, u64), f64)],
    fuzz: u64,
) -> std::collections::HashMap<(u64, u64), (u64, u64)> {
    use std::collections::HashMap;
    // Aggregate support per distinct junction.
    let mut support: HashMap<(u64, u64), f64> = HashMap::new();
    for (j, s) in junctions_with_support {
        *support.entry(*j).or_insert(0.0) += *s;
    }
    let mut map: HashMap<(u64, u64), (u64, u64)> = HashMap::with_capacity(support.len());
    if fuzz == 0 {
        for j in support.keys() {
            map.insert(*j, *j);
        }
        return map;
    }
    // Process distinct junctions by descending support (coord tie-break for
    // determinism): the most-supported become canonicals first; lower-support
    // neighbours within `fuzz` on BOTH coordinates snap to the nearest canonical.
    let mut distinct: Vec<((u64, u64), f64)> = support.into_iter().collect();
    distinct.sort_by(|a, b| {
        b.1.partial_cmp(&a.1)
            .unwrap_or(std::cmp::Ordering::Equal)
            .then(a.0.cmp(&b.0))
    });
    let mut canonicals: Vec<(u64, u64)> = Vec::new();
    for (j, _) in &distinct {
        let canon = canonicals.iter().find(|c| {
            j.0.abs_diff(c.0) <= fuzz && j.1.abs_diff(c.1) <= fuzz
        });
        match canon {
            Some(c) => {
                map.insert(*j, *c);
            }
            None => {
                canonicals.push(*j);
                map.insert(*j, *j);
            }
        }
    }
    map
}

/// Determine 5'-degradation containment folds for read-chain collapse.
///
/// A 5'-degraded read is the same molecule as a longer read but truncated at the
/// 5' cap, so its intron chain is a strand-appropriate sub-chain of the longer
/// "super" chain — a SUFFIX on the `+` strand (missing upstream/left introns), a
/// PREFIX on `-` (5' is the right end) — sharing the 3' terminus within `end3_tol`
/// bp. Such a chain is folded into the longest containing super (IsoSeq collapse
/// `--max-5-diff` / merge-5-shorter), rather than minting a truncated isoform.
///
/// `chains[i] = (intron_chain, abundance, min_start, max_end)`. Returns, per index,
/// `Some(parent_idx)` of the strictly-longer chain it folds into (longest-first;
/// callers resolve transitively to a root), or `None` if it is kept. Deterministic.
// Retained as the equivalence anchor for `legacy_wrapper_matches_5p_only` (proves the
// generalized `compute_degrade_folds(.., rc=false)` is byte-identical to the original
// 5'-only behavior). Production routes through `compute_degrade_folds` directly.
#[allow(dead_code)]
fn compute_5p_degrade_folds(
    chains: &[(Vec<(u64, u64)>, f64, u64, u64)],
    strand: char,
    end3_tol: u64,
) -> Vec<Option<usize>> {
    compute_degrade_folds(chains, strand, end3_tol, false)
}

/// `true` iff `sub` appears as a CONTIGUOUS sub-sequence of `sup` (i.e. there exists
/// an offset `k` with `sup[k..k+sub.len()] == sub`). Used by the read-coherence
/// degradation collapse to recognize internal fragments. `sub` empty => `false`.
fn is_contiguous_subrun<T: PartialEq>(sub: &[T], sup: &[T]) -> bool {
    if sub.is_empty() || sub.len() > sup.len() {
        return false;
    }
    sup.windows(sub.len()).any(|w| w == sub)
}

/// Determine degradation containment folds for read-chain collapse.
///
/// Default (`read_coherence = false`): the original 5'-degradation rule — a
/// 5'-truncated read is a strand-appropriate sub-chain of a longer "super" sharing
/// the 3' terminus within `end3_tol` (a SUFFIX on `+`, a PREFIX on `-`). This path is
/// BYTE-IDENTICAL to the historical `compute_5p_degrade_folds`.
///
/// Read-coherence (`read_coherence = true`): additionally fold
///   - 3'-truncated chains (the opposite-end truncation: a PREFIX on `+` / SUFFIX on
///     `-` sharing the 5' terminus within `end3_tol`), and
///   - internal-fragment chains (junctions are a CONTIGUOUS sub-run of a longer
///     chain; no terminus required — both ends were degraded).
/// A chain only folds when its junction list is a contiguous sub-sequence of the
/// super's (never a non-contiguous skip, which is a genuinely distinct isoform).
///
/// `chains[i] = (intron_chain, abundance, min_start, max_end)`. Returns, per index,
/// `Some(parent_idx)` of the strictly-longer chain it folds into (longest-first;
/// callers resolve transitively to a root), or `None` if it is kept. Deterministic.
fn compute_degrade_folds(
    chains: &[(Vec<(u64, u64)>, f64, u64, u64)],
    strand: char,
    end3_tol: u64,
    read_coherence: bool,
) -> Vec<Option<usize>> {
    let n = chains.len();
    let mut fold: Vec<Option<usize>> = vec![None; n];
    if n < 2 {
        return fold;
    }
    // Candidate-super order: longest chain first (so a shorter chain matches its
    // LONGEST container first), tie-break by abundance desc then chain coords.
    let mut order: Vec<usize> = (0..n).collect();
    order.sort_by(|&i, &j| {
        chains[j]
            .0
            .len()
            .cmp(&chains[i].0.len())
            .then(
                chains[j]
                    .1
                    .partial_cmp(&chains[i].1)
                    .unwrap_or(std::cmp::Ordering::Equal),
            )
            .then(chains[i].0.cmp(&chains[j].0))
    });
    for &b in &order {
        let (bchain, _, bmin, bmax) = &chains[b];
        for &a in &order {
            if a == b {
                continue;
            }
            let (achain, _, amin, amax) = &chains[a];
            if achain.len() <= bchain.len() {
                continue; // super must be STRICTLY longer
            }
            // 5'-degradation (default + rc): strand-appropriate truncation sharing 3'.
            let fold_5p = match strand {
                // 5' is the right end on '-': a 5'-truncation drops right/late
                // introns -> shorter chain is a PREFIX; 3' terminus is min_start.
                '-' => achain.starts_with(bchain.as_slice()) && bmin.abs_diff(*amin) <= end3_tol,
                // '+': 5'-truncation drops left/early introns -> SUFFIX; 3' is max_end.
                _ => achain.ends_with(bchain.as_slice()) && bmax.abs_diff(*amax) <= end3_tol,
            };
            let contained = if !read_coherence {
                fold_5p
            } else {
                // 3'-degradation: opposite-end truncation sharing the 5' terminus.
                let fold_3p = match strand {
                    // 5' is the left end on '+': a 3'-truncation drops late/right
                    // introns -> shorter chain is a PREFIX; 5' terminus is min_start.
                    '-' => achain.ends_with(bchain.as_slice()) && bmax.abs_diff(*amax) <= end3_tol,
                    _ => achain.starts_with(bchain.as_slice()) && bmin.abs_diff(*amin) <= end3_tol,
                };
                // Internal fragment: junctions are a contiguous sub-run of the super
                // (both ends degraded; no terminus required).
                let fold_internal = is_contiguous_subrun(bchain.as_slice(), achain.as_slice());
                fold_5p || fold_3p || fold_internal
            };
            if contained {
                fold[b] = Some(a); // longest-first order -> first match is the longest super
                break;
            }
        }
    }
    fold
}

/// Extract transcripts by grouping long-read transfrags by their exact intron chain.
///
/// Enable via `RUSTLE_READCHAIN=1`. Each unique (donor, acceptor) chain observed in
/// long-read transfrags becomes one transcript, provided total abundance >= readthr.
/// This avoids flow-composition artifacts: every emitted transcript is directly
/// supported by one or more actual reads with matching splice sites.
///
/// Cluster single-exon read spans into locus groups for `--read-chain-single`.
///
/// Single-exon reads carry no junctions to key on, so they are grouped by genomic
/// overlap: spans are merged when the gap between them is `<= min_gap` (overlap =
/// gap 0). Each cluster becomes one single-exon transcript spanning [min start,
/// max end] with summed abundance. Deterministic (sorts internally by coord).
fn cluster_single_exon_spans(
    spans: &[(u64, u64, f64)],
    min_gap: u64,
) -> Vec<(u64, u64, f64)> {
    if spans.is_empty() {
        return Vec::new();
    }
    let mut sorted: Vec<(u64, u64, f64)> = spans.to_vec();
    sorted.sort_by(|a, b| a.0.cmp(&b.0).then(a.1.cmp(&b.1)));
    let mut out: Vec<(u64, u64, f64)> = Vec::new();
    let mut cur = sorted[0];
    for &(s, e, ab) in &sorted[1..] {
        if s <= cur.1.saturating_add(min_gap) {
            // Overlapping (or within min_gap): extend the cluster, sum abundance.
            cur.1 = cur.1.max(e);
            cur.2 += ab;
        } else {
            out.push(cur);
            cur = (s, e, ab);
        }
    }
    out.push(cur);
    out
}

/// Single-exon transfrags (empty intron chain) are skipped unless `--read-chain-single`
/// (RUSTLE_READCHAIN_SINGLE), which clusters them by locus and emits one single-exon
/// transcript per cluster. NOT RECOMMENDED / strictly opt-in: measured on GGO chr19
/// HiFi every emitted single-exon call was a false positive (390/390 class u/i/p, 0
/// RefSeq matches) and it cost real multi-exon matches via downstream interference
/// (Sn 86.2->84.3, Pr 57.7->50.2) — single-exon loci on this data are spurious and no
/// abundance gate separates them. The switch is now live (previously dead) for inputs
/// where single-exon isoforms are real.
/// Effective terminal boundaries for a read-coherence chain group (opt-in trim). The collapse
/// otherwise sets the first-exon 5' and last-exon 3' boundaries to the EXTREME read extent
/// (`min_start`/`max_end`) over all reads sharing the intron chain, so a single run-on /
/// intra-primed molecule inflates the terminal exon. Here we cluster the per-read terminal
/// positions (`start_positions`/`end_positions`, each `(pos, weight)`) and snap each terminal to
/// the most-extreme cluster that clears a SUPPORT FLOOR, dropping minority runaway ends while
/// preserving genuinely-supported length.
///
/// The floor is RELATIVE to the group's OWN read weight (the summed position weights) with an
/// absolute minimum, so thin groups find no surviving cluster and fall back to the extreme
/// (conservative). `first_donor` / `last_acceptor` bound the terminal exons so the result always
/// keeps a positive-length terminal exon (`eff_start < first_donor`, `eff_end > last_acceptor`);
/// otherwise the extreme is kept. Returns `(eff_start, eff_end)`.
fn effective_terminal_boundaries(
    start_positions: &[(u64, f64)],
    end_positions: &[(u64, f64)],
    min_start: u64,
    max_end: u64,
    first_donor: u64,
    last_acceptor: u64,
    trim_frac: f64,
    trim_minabs: f64,
) -> (u64, u64) {
    let pos_w: f64 = start_positions.iter().map(|(_, w)| w).sum();
    let floor = ((trim_frac * pos_w).max(trim_minabs)).ceil().max(2.0) as u64;
    let sc =
        crate::tss_tts::cluster_positions_with_counts(start_positions, crate::tss_tts::CPAS_POS_BIN, floor);
    let ec =
        crate::tss_tts::cluster_positions_with_counts(end_positions, crate::tss_tts::CPAS_POS_BIN, floor);
    let mut eff_start = min_start;
    let mut eff_end = max_end;
    // 5' start: most-upstream surviving cluster (sc is ascending by position).
    if let Some(&(c, _)) = sc.first() {
        if c >= min_start && c < first_donor {
            eff_start = c;
        }
    }
    // 3' end: most-downstream surviving cluster.
    if let Some(&(c, _)) = ec.last() {
        if c <= max_end && c > last_acceptor {
            eff_end = c;
        }
    }
    (eff_start, eff_end)
}

pub fn extract_transcripts_readchain(
    graph: &Graph,
    transfrags: &mut [GraphTransfrag],
    bundle_chrom: &str,
    bundle_strand: char,
    config: &RunConfig,
) -> Vec<Transcript> {
    let debug = std::env::var_os("RUSTLE_READCHAIN_DEBUG").is_some();
    // Terminal effective-boundary trim (opt-in; default OFF => byte-identical). When on, the first/last
    // exon boundaries are snapped from the EXTREME read extent (min_start/max_end) to the most-extreme
    // per-read terminal CLUSTER that clears a support floor, dropping minority runaway ends. `frac` is
    // the relative support floor (fraction of the group's own read weight); `minabs` the absolute floor.
    let trim_terminal = std::env::var_os("RUSTLE_READCHAIN_TRIM_TERMINAL").is_some();
    let trim_frac: f64 = std::env::var("RUSTLE_READCHAIN_TRIM_FRAC")
        .ok()
        .and_then(|v| v.parse().ok())
        .unwrap_or(0.10);
    let trim_minabs: f64 = std::env::var("RUSTLE_READCHAIN_TRIM_MINABS")
        .ok()
        .and_then(|v| v.parse().ok())
        .unwrap_or(2.0);
    let source_id = graph.source_id;
    let sink_id = graph.sink_id;
    let min_flow = config.readthr.max(EPSILON);

    // Group long-read transfrags by intron chain.
    // Value: (total_abundance, min_start, max_end, representative_node_ids,
    //         per-read first-exon starts [(pos,weight)], per-read last-exon ends [(pos,weight)]).
    // The two position vectors are populated only when `trim_terminal` is set (otherwise they stay
    // empty and unused, keeping the default path byte-identical).
    let mut chain_groups: std::collections::HashMap<
        Vec<(u64, u64)>,
        (f64, u64, u64, Vec<usize>, Vec<(u64, f64)>, Vec<(u64, f64)>),
    > = std::collections::HashMap::new();

    // Pass 1: collect each long-read transfrag's intron chain + span + abundance.
    // --read-chain-single (RUSTLE_READCHAIN_SINGLE) also collects single-exon spans.
    let read_chain_single = std::env::var_os("RUSTLE_READCHAIN_SINGLE").is_some();
    let mut se_spans: Vec<(u64, u64, f64)> = Vec::new();
    let mut rc_entries: Vec<(Vec<(u64, u64)>, f64, u64, u64, Vec<usize>)> = Vec::new();
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

        // Single-exon reads (empty intron chain): the baseline handles these unless
        // --read-chain-single, which collects their spans for locus clustering below.
        if chain.is_empty() {
            if read_chain_single {
                let s = graph.nodes.get(inner[0]).map(|n| n.start).unwrap_or(0);
                let e = graph
                    .nodes
                    .get(*inner.last().unwrap())
                    .map(|n| n.end)
                    .unwrap_or(0);
                if e > s {
                    se_spans.push((s, e, tf.abundance));
                }
            }
            continue;
        }

        let first_node = inner[0];
        let last_node = *inner.last().unwrap();
        let tf_start = graph.nodes.get(first_node).map(|n| n.start).unwrap_or(0);
        let tf_end = graph.nodes.get(last_node).map(|n| n.end).unwrap_or(0);
        rc_entries.push((chain, tf.abundance, tf_start, tf_end, inner.to_vec()));
    }

    // Fuzzy-junction canonicalization (IsoSeq `--max-fuzzy-junction` analogue):
    // snap junctions wobbling within RUSTLE_READCHAIN_FUZZ bp to a per-cluster
    // consensus so alignment noise doesn't split one isoform into near-duplicate
    // chains. DEFAULT 0 (exact, no merging) — measured NET-NEGATIVE on GGO chr19
    // HiFi (deterministic): clean alignments have little junction wobble, and
    // near-junctions on paralog loci are often REAL distinct copies that merging
    // wrongly collapses — fuzz=5 lost RefSeq matches 1441->1424 (Sn 78.4->77.4)
    // with NO precision gain (Pr ~62.1 flat).
    // Kept as an OPT-IN knob (RUSTLE_READCHAIN_FUZZ=N) for noisier inputs (e.g. ONT)
    // where genuine wobble does split isoforms. Identity map when fuzz == 0.
    let fuzz: u64 = std::env::var("RUSTLE_READCHAIN_FUZZ")
        .ok()
        .and_then(|v| v.parse().ok())
        .unwrap_or(0);
    let jmap = if fuzz > 0 {
        let mut js: Vec<((u64, u64), f64)> = Vec::new();
        for (chain, abund, ..) in &rc_entries {
            for j in chain {
                js.push((*j, *abund));
            }
        }
        build_fuzzy_junction_map(&js, fuzz)
    } else {
        std::collections::HashMap::new()
    };

    // Pass 2: group by the (fuzzy-canonicalized) intron chain. The canonical chain
    // also supplies the emitted internal exon boundaries below (consensus coords);
    // first/last exon ends still float via min_start/max_end across the group.
    for (chain, abund, tf_start, tf_end, inner) in rc_entries {
        let key: Vec<(u64, u64)> = if fuzz > 0 {
            chain.iter().map(|j| jmap.get(j).copied().unwrap_or(*j)).collect()
        } else {
            chain
        };
        let entry = chain_groups
            .entry(key)
            .or_insert((0.0, tf_start, tf_end, inner, Vec::new(), Vec::new()));
        entry.0 += abund;
        if tf_start < entry.1 {
            entry.1 = tf_start;
        }
        if tf_end > entry.2 {
            entry.2 = tf_end;
        }
        if trim_terminal {
            // Keep this read's own terminal positions for effective-boundary clustering at emit time.
            entry.4.push((tf_start, abund));
            entry.5.push((tf_end, abund));
        }
    }

    // 5'-degradation containment collapse (DEFAULT ON, tol=100 bp — IsoSeq's own
    // default; RUSTLE_READCHAIN_DEGRADE_TOL=N to override, =0 to disable). Fold a
    // chain that is a strand-appropriate 5'-truncation of a longer chain (sharing
    // the 3' terminus within N bp) into that longer "super" chain, summing its
    // abundance in. Runs before the isofrac threshold so the super's boosted
    // abundance counts — this redistributes the high read-counts that 5'-degraded
    // fragments otherwise pile onto truncated keys (where they inflate max_abund and
    // suppress real minor isoforms). Measured on GGO chr19 HiFi: Sn 78.4->86.2,
    // matched 1441->1586 (+145 real RefSeq tx), Pr 62.1->57.7 (~F1-neutral),
    // deterministic. 4456 chains folded across 449 bundles.
    let degrade_tol: u64 = std::env::var("RUSTLE_READCHAIN_DEGRADE_TOL")
        .ok()
        .and_then(|v| v.parse::<u64>().ok())
        .unwrap_or(100);
    if degrade_tol > 0 {
        let tol = degrade_tol;
        let keys: Vec<Vec<(u64, u64)>> = chain_groups.keys().cloned().collect();
        let view: Vec<(Vec<(u64, u64)>, f64, u64, u64)> = keys
            .iter()
            .map(|k| {
                let v = &chain_groups[k];
                (k.clone(), v.0, v.1, v.2)
            })
            .collect();
        // Under --read-coherence (RUSTLE_READ_COHERENCE), generalize the collapse to
        // also fold 3'-truncated + internal-fragment chains. Default (env unset) keeps
        // the 5'-only behavior byte-for-byte.
        let rc_collapse = std::env::var_os("RUSTLE_READ_COHERENCE").is_some();
        let folds = compute_degrade_folds(&view, bundle_strand, tol, rc_collapse);
        // Resolve each folded chain to its transitive root super and accumulate.
        let mut add_to_root: Vec<f64> = vec![0.0; keys.len()];
        let mut removed = 0usize;
        for i in 0..keys.len() {
            if folds[i].is_some() {
                let mut r = i;
                while let Some(p) = folds[r] {
                    r = p;
                }
                add_to_root[r] += view[i].1;
                chain_groups.remove(&keys[i]);
                removed += 1;
            }
        }
        for i in 0..keys.len() {
            if add_to_root[i] > 0.0 {
                if let Some(entry) = chain_groups.get_mut(&keys[i]) {
                    entry.0 += add_to_root[i];
                }
            }
        }
        if debug && removed > 0 {
            eprintln!(
                "[READCHAIN] degradation collapse: folded {} truncated/fragment chain(s) into supers (tol={})",
                removed, tol
            );
        }
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
    // Tie-break on the intron-chain coordinates so equal-abundance chains have a
    // total order independent of HashMap iteration order — without this the
    // depletion order (and thus the emitted transcript set) varied run-to-run.
    let mut sorted_chains: Vec<(
        &Vec<(u64, u64)>,
        &(f64, u64, u64, Vec<usize>, Vec<(u64, f64)>, Vec<(u64, f64)>),
    )> = chain_groups.iter().collect();
    sorted_chains.sort_by(|a, b| {
        b.1 .0
            .partial_cmp(&a.1 .0)
            .unwrap_or(std::cmp::Ordering::Equal)
            .then_with(|| a.0.cmp(b.0))
    });

    for (chain, (total_abund, min_start, max_end, rep_nodes, start_positions, end_positions)) in
        sorted_chains
    {
        if *total_abund < isofrac_threshold {
            if debug {
                eprintln!(
                    "[READCHAIN] chain_len={} abund={:.2} < threshold={:.2} FILTERED",
                    chain.len(), total_abund, isofrac_threshold
                );
            }
            continue;
        }

        // Terminal effective-boundary trim (opt-in; default off => eff_* == min_start/max_end, so the
        // emitted exons are byte-identical). The collapse otherwise sets the first-exon 5' and last-exon
        // 3' boundaries to the EXTREME read extent over all reads sharing the intron chain, so one
        // run-on / intra-primed molecule inflates the terminal exon. Cluster the per-read terminal
        // positions and snap each terminal to the most-extreme cluster that clears a SUPPORT FLOOR,
        // dropping minority runaway ends while preserving genuinely-supported length. The floor is
        // RELATIVE to the group's own read weight (sum of the position-vector weights) -- NOT
        // total_abund, which the degrade-collapse inflates with folded-fragment abundance whose
        // positions are not in these vectors -- with an absolute minimum, so thin groups find no
        // surviving cluster and fall back to the extreme (conservative). Internal exons are untouched.
        let (mut eff_start, mut eff_end) = (*min_start, *max_end);
        if trim_terminal {
            let (s, e) = effective_terminal_boundaries(
                start_positions,
                end_positions,
                *min_start,
                *max_end,
                chain[0].0,
                chain.last().unwrap().1,
                trim_frac,
                trim_minabs,
            );
            eff_start = s;
            eff_end = e;
            if debug && (eff_start != *min_start || eff_end != *max_end) {
                let pos_w: f64 = start_positions.iter().map(|(_, w)| w).sum();
                eprintln!(
                    "[READCHAIN_TRIM] start {}->{} (-{}) end {}->{} (-{}) pos_w={:.1}",
                    min_start, eff_start, eff_start - *min_start,
                    max_end, eff_end, *max_end - eff_end, pos_w
                );
            }
        }

        // Build exon list from intron chain + span.
        // Internal exon boundaries are fully determined by (acceptor_{i-1}, donor_i);
        // only the first-exon start and last-exon end vary across reads.
        let mut exons: Vec<(u64, u64)> = Vec::new();
        exons.push((eff_start, chain[0].0)); // first exon
        for i in 1..chain.len() {
            exons.push((chain[i - 1].1, chain[i].0)); // internal exons
        }
        exons.push((chain.last().unwrap().1, eff_end)); // last exon

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

    // --read-chain-single: emit clustered single-exon transcripts. Single-exon
    // calls are FP-prone (the stranded-SE injector emitted 38/38 non-RefSeq), so
    // this is GATED behind --read-chain-single AND a strict abundance floor:
    // max(singlethr≈4.75, the bundle's isofrac threshold) — a single-exon cluster
    // must clear the same single-exon bar the flow assembler uses and not be a
    // minor fraction of a spliced locus. RUSTLE_READCHAIN_SINGLE_GAP merges spans
    // within N bp (default 0 = overlap only).
    if read_chain_single && !se_spans.is_empty() {
        let min_gap: u64 = std::env::var("RUSTLE_READCHAIN_SINGLE_GAP")
            .ok()
            .and_then(|v| v.parse().ok())
            .unwrap_or(0);
        let se_floor = config.singlethr.max(isofrac_threshold);
        let mut n_se = 0usize;
        for (s, e, abund) in cluster_single_exon_spans(&se_spans, min_gap) {
            if abund < se_floor {
                continue;
            }
            if len_half_open(s, e) < config.min_transcript_length {
                continue;
            }
            out.push(Transcript {
                chrom: bundle_chrom.to_string(),
                strand: bundle_strand,
                exon_cov: vec![abund],
                exons: vec![(s, e)],
                coverage: abund,
                source: Some("readchain_single".to_string()),
                is_longread: true,
                longcov: abund,
                ..Transcript::default()
            });
            n_se += 1;
        }
        if debug && n_se > 0 {
            eprintln!("[READCHAIN] single-exon: emitted {} clustered SE transcript(s)", n_se);
        }
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
    // Tie-break on the intron-chain coordinates so equal-abundance chains have a
    // total order independent of HashMap iteration order — without this the
    // depletion order (and thus the emitted transcript set) varied run-to-run.
    let mut sorted_chains: Vec<(&Vec<(u64, u64)>, &(f64, u64, u64, Vec<usize>))> =
        chain_groups.iter().collect();
    sorted_chains.sort_by(|a, b| {
        b.1 .0
            .partial_cmp(&a.1 .0)
            .unwrap_or(std::cmp::Ordering::Equal)
            .then_with(|| a.0.cmp(b.0))
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

#[cfg(test)]
mod fuzzy_junction_tests {
    use super::build_fuzzy_junction_map;

    // A lone junction maps to itself.
    #[test]
    fn identity_single() {
        let m = build_fuzzy_junction_map(&[((100, 200), 5.0)], 5);
        assert_eq!(m.get(&(100, 200)), Some(&(100, 200)));
    }

    // fuzz=0 preserves exact behaviour: every junction maps to itself, no merging.
    #[test]
    fn fuzz_zero_is_identity() {
        let js = [((100, 200), 5.0), ((102, 201), 1.0)];
        let m = build_fuzzy_junction_map(&js, 0);
        assert_eq!(m.get(&(100, 200)), Some(&(100, 200)));
        assert_eq!(m.get(&(102, 201)), Some(&(102, 201)));
    }

    // A low-support junction within `fuzz` of a higher-support one snaps to it;
    // a far junction stays its own canonical.
    #[test]
    fn snaps_wobble_to_higher_support() {
        let js = [
            ((100, 200), 10.0), // dominant
            ((102, 201), 2.0),  // +2/+1 wobble -> snaps to (100,200)
            ((500, 600), 3.0),  // far -> own canonical
        ];
        let m = build_fuzzy_junction_map(&js, 5);
        assert_eq!(m.get(&(100, 200)), Some(&(100, 200)));
        assert_eq!(m.get(&(102, 201)), Some(&(100, 200)));
        assert_eq!(m.get(&(500, 600)), Some(&(500, 600)));
    }

    // Outside the fuzz window on EITHER coordinate -> not merged.
    #[test]
    fn beyond_fuzz_not_merged() {
        let js = [((100, 200), 10.0), ((108, 200), 1.0)]; // donor +8 > fuzz 5
        let m = build_fuzzy_junction_map(&js, 5);
        assert_eq!(m.get(&(108, 200)), Some(&(108, 200)));
        // acceptor out of window
        let js2 = [((100, 200), 10.0), ((100, 207), 1.0)]; // acceptor +7 > 5
        let m2 = build_fuzzy_junction_map(&js2, 5);
        assert_eq!(m2.get(&(100, 207)), Some(&(100, 207)));
    }

    // Support is aggregated across duplicate observations before choosing canonical.
    #[test]
    fn aggregates_support_picks_dominant() {
        // (100,200) appears twice (total 4) vs (103,203) once (5). The single
        // (103,203) has higher support -> it becomes canonical, (100,200) snaps to it.
        let js = [((100, 200), 2.0), ((100, 200), 2.0), ((103, 203), 5.0)];
        let m = build_fuzzy_junction_map(&js, 5);
        assert_eq!(m.get(&(103, 203)), Some(&(103, 203)));
        assert_eq!(m.get(&(100, 200)), Some(&(103, 203)));
    }
}

#[cfg(test)]
mod degrade_collapse_tests {
    use super::{compute_5p_degrade_folds, compute_degrade_folds};

    // --- Generalized (read-coherence) collapse: 3' + internal + still-5' ---

    // 3'-truncated chain: shares the 5' terminus, junctions are a PREFIX of a longer
    // chain (on '+', 5' is the left end -> 3'-truncation drops late/right introns ->
    // PREFIX). Folds only under read-coherence.
    #[test]
    fn rc_plus_prefix_3p_trunc_folds() {
        let chains = vec![
            (vec![(100, 200), (300, 400), (500, 600)], 10.0, 50, 650), // super (idx0)
            (vec![(100, 200), (300, 400)], 2.0, 50, 450),              // 3'-trunc (prefix), same 5' (min_start)
        ];
        // Default 5'-only: a '+' PREFIX is NOT a 5'-truncation -> no fold.
        let f5 = compute_degrade_folds(&chains, '+', 100, false);
        assert_eq!(f5[1], None, "5'-only mode must not fold a 3'-truncation");
        // Read-coherence: prefix sharing the 5' terminus folds.
        let frc = compute_degrade_folds(&chains, '+', 100, true);
        assert_eq!(frc[0], None);
        assert_eq!(frc[1], Some(0), "rc mode folds the 3'-truncation prefix");
    }

    // Internal fragment: junctions are a contiguous sub-run of a longer chain but
    // share NEITHER terminus. Folds only under read-coherence (no terminus req).
    #[test]
    fn rc_internal_fragment_folds() {
        let chains = vec![
            (vec![(100, 200), (300, 400), (500, 600), (700, 800)], 10.0, 50, 850), // super idx0
            (vec![(300, 400), (500, 600)], 2.0, 250, 650), // internal sub-run, neither terminus shared
        ];
        let f5 = compute_degrade_folds(&chains, '+', 100, false);
        assert_eq!(f5[1], None, "5'-only mode must not fold an internal fragment");
        let frc = compute_degrade_folds(&chains, '+', 100, true);
        assert_eq!(frc[0], None);
        assert_eq!(frc[1], Some(0), "rc mode folds the internal fragment");
    }

    // The existing 5'-truncation (suffix, shares 3') must STILL fold under rc mode.
    #[test]
    fn rc_still_folds_5p_suffix() {
        let chains = vec![
            (vec![(100, 200), (300, 400), (500, 600)], 10.0, 50, 650),
            (vec![(300, 400), (500, 600)], 2.0, 250, 650), // 5'-trunc suffix, same max_end
        ];
        let frc = compute_degrade_folds(&chains, '+', 100, true);
        assert_eq!(frc[0], None);
        assert_eq!(frc[1], Some(0));
    }

    // A genuinely-distinct short isoform whose junctions are NOT a contiguous sub-run
    // of any longer chain must NOT fold, even under rc mode.
    #[test]
    fn rc_distinct_isoform_not_folded() {
        let chains = vec![
            (vec![(100, 200), (300, 400), (500, 600)], 10.0, 50, 650), // super
            (vec![(100, 200), (500, 600)], 5.0, 50, 650), // skips (300,400) -> NOT a contiguous sub-run
        ];
        let frc = compute_degrade_folds(&chains, '+', 100, true);
        assert_eq!(frc[0], None);
        assert_eq!(frc[1], None, "non-contiguous sub-chain is a distinct isoform, never folds");
    }

    // The old 5'-only wrapper must remain byte-equivalent to compute_degrade_folds
    // with read_coherence=false (default path unchanged).
    #[test]
    fn legacy_wrapper_matches_5p_only() {
        let chains = vec![
            (vec![(100, 200), (300, 400), (500, 600)], 10.0, 50, 650),
            (vec![(300, 400), (500, 600)], 2.0, 250, 650), // 5'-trunc suffix
            (vec![(100, 200), (300, 400)], 2.0, 50, 450),  // 3'-trunc prefix (must NOT fold here)
        ];
        let legacy = compute_5p_degrade_folds(&chains, '+', 100);
        let general = compute_degrade_folds(&chains, '+', 100, false);
        assert_eq!(legacy, general);
        assert_eq!(legacy[1], Some(0)); // suffix folds
        assert_eq!(legacy[2], None); // prefix does not (5'-only)
    }


    // + strand: a 2-junction SUFFIX of a 3-junction chain, same 3' (max_end),
    // folds into the longer chain.
    #[test]
    fn plus_suffix_folds() {
        let chains = vec![
            (vec![(100, 200), (300, 400), (500, 600)], 10.0, 50, 650), // super (idx0)
            (vec![(300, 400), (500, 600)], 2.0, 250, 650),              // 5'-trunc (idx1)
        ];
        let f = compute_5p_degrade_folds(&chains, '+', 100);
        assert_eq!(f[0], None);
        assert_eq!(f[1], Some(0));
    }

    // - strand: 5' is the right end, so a 5'-truncation is a PREFIX; shares the 3'
    // (min_start) terminus.
    #[test]
    fn minus_prefix_folds() {
        let chains = vec![
            (vec![(100, 200), (300, 400), (500, 600)], 10.0, 50, 650),
            (vec![(100, 200), (300, 400)], 2.0, 50, 450), // prefix, same min_start
        ];
        let f = compute_5p_degrade_folds(&chains, '-', 100);
        assert_eq!(f[0], None);
        assert_eq!(f[1], Some(0));
    }

    // A suffix whose 3' end is FAR from the super's 3' is NOT a degradation of it.
    #[test]
    fn rejects_3p_mismatch() {
        let chains = vec![
            (vec![(100, 200), (300, 400), (500, 600)], 10.0, 50, 650),
            (vec![(300, 400), (500, 600)], 2.0, 250, 900), // max_end 900 vs 650 -> >tol
        ];
        let f = compute_5p_degrade_folds(&chains, '+', 100);
        assert_eq!(f[1], None);
    }

    // Two structurally-unrelated chains do not fold.
    #[test]
    fn unrelated_not_folded() {
        let chains = vec![
            (vec![(100, 200), (300, 400)], 10.0, 50, 450),
            (vec![(700, 800), (900, 1000)], 5.0, 650, 1050),
        ];
        let f = compute_5p_degrade_folds(&chains, '+', 100);
        assert_eq!(f[0], None);
        assert_eq!(f[1], None);
    }

    // A nested chain folds into the LONGEST containing super (not an intermediate).
    #[test]
    fn folds_into_longest_super() {
        let chains = vec![
            (vec![(100, 200), (300, 400), (500, 600)], 10.0, 50, 650), // idx0 longest
            (vec![(300, 400), (500, 600)], 4.0, 250, 650),             // idx1 mid (suffix of 0)
            (vec![(500, 600)], 1.0, 450, 650),                         // idx2 short (suffix of 0 and 1)
        ];
        let f = compute_5p_degrade_folds(&chains, '+', 100);
        assert_eq!(f[0], None);
        assert_eq!(f[1], Some(0));
        assert_eq!(f[2], Some(0)); // longest container, not idx1
    }

    // Equal chains never fold (a strictly-longer super is required).
    #[test]
    fn equal_length_no_fold() {
        let chains = vec![
            (vec![(100, 200), (300, 400)], 10.0, 50, 450),
            (vec![(100, 200), (300, 400)], 3.0, 50, 450),
        ];
        let f = compute_5p_degrade_folds(&chains, '+', 100);
        assert_eq!(f[0], None);
        assert_eq!(f[1], None);
    }
}

#[cfg(test)]
mod single_exon_cluster_tests {
    use super::cluster_single_exon_spans;

    #[test]
    fn empty_input() {
        assert!(cluster_single_exon_spans(&[], 0).is_empty());
    }

    #[test]
    fn overlapping_merge_sum_abundance() {
        let c = cluster_single_exon_spans(&[(100, 200, 3.0), (150, 250, 2.0)], 0);
        assert_eq!(c, vec![(100, 250, 5.0)]);
    }

    #[test]
    fn disjoint_stay_separate() {
        let mut c = cluster_single_exon_spans(&[(500, 600, 2.0), (100, 200, 3.0)], 0);
        c.sort_by_key(|x| x.0);
        assert_eq!(c, vec![(100, 200, 3.0), (500, 600, 2.0)]);
    }

    #[test]
    fn gap_within_min_gap_merges() {
        // gap of 10 (200 -> 210) merges at min_gap=20, stays separate at min_gap=0.
        assert_eq!(
            cluster_single_exon_spans(&[(100, 200, 3.0), (210, 300, 2.0)], 20),
            vec![(100, 300, 5.0)]
        );
        assert_eq!(
            cluster_single_exon_spans(&[(100, 200, 3.0), (210, 300, 2.0)], 0).len(),
            2
        );
    }

    #[test]
    fn input_order_independent() {
        let a = cluster_single_exon_spans(&[(100, 200, 1.0), (150, 250, 1.0), (500, 600, 1.0)], 0);
        let b = cluster_single_exon_spans(&[(500, 600, 1.0), (150, 250, 1.0), (100, 200, 1.0)], 0);
        assert_eq!(a, b);
    }
}

#[cfg(test)]
mod terminal_trim_tests {
    use super::effective_terminal_boundaries;

    /// Build a per-read position vector from (position, read_count) groups (unit weights).
    fn positions(groups: &[(u64, usize)]) -> Vec<(u64, f64)> {
        let mut v = Vec::new();
        for &(p, n) in groups {
            for _ in 0..n {
                v.push((p, 1.0));
            }
        }
        v
    }

    #[test]
    fn no_trim_when_terminals_are_clean() {
        // 6 reads share one TSS and one TES -> a single supported cluster each -> no trim.
        let (s, e) = effective_terminal_boundaries(
            &positions(&[(1000, 6)]), &positions(&[(3000, 6)]),
            1000, 3000, 2000, 2500, 0.10, 2.0,
        );
        assert_eq!((s, e), (1000, 3000));
    }

    #[test]
    fn lone_5p_runaway_is_trimmed_to_supported_tss() {
        // 1 runaway read 400 bp upstream + 5 at the real TSS; floor=2 drops the singleton.
        let (s, e) = effective_terminal_boundaries(
            &positions(&[(600, 1), (1000, 5)]), &positions(&[(3000, 6)]),
            600, 3000, 2000, 2500, 0.10, 2.0,
        );
        assert_eq!(s, 1000, "5' boundary should snap to the supported cluster");
        assert_eq!(e, 3000, "3' boundary untouched");
    }

    #[test]
    fn thin_group_falls_back_to_extreme() {
        // 2 reads, both singleton clusters -> nothing clears floor=2 -> keep the extreme (conservative).
        let (s, e) = effective_terminal_boundaries(
            &positions(&[(500, 1), (1500, 1)]), &positions(&[(3000, 1), (4000, 1)]),
            500, 4000, 2000, 2500, 0.10, 2.0,
        );
        assert_eq!((s, e), (500, 4000));
    }

    #[test]
    fn supported_minority_alt_tss_is_preserved() {
        // upstream alt-TSS with 5 reads (>= floor=3) is real, not a runaway -> preserved, not trimmed.
        let (s, _e) = effective_terminal_boundaries(
            &positions(&[(700, 5), (1000, 25)]), &positions(&[(3000, 30)]),
            700, 3000, 2000, 2500, 0.10, 2.0,
        );
        assert_eq!(s, 700, "a terminal supported by >= floor reads must be preserved");
    }

    #[test]
    fn never_moves_boundary_past_the_junction() {
        // Degenerate: the only cluster center lands at/after the donor (and before the acceptor) ->
        // guard rejects the move and keeps the extreme, so the terminal exon stays well-defined.
        let (s, e) = effective_terminal_boundaries(
            &positions(&[(2500, 5)]), &positions(&[(2400, 5)]),
            2500, 2400, 2000, 2500, 0.10, 2.0,
        );
        assert_eq!((s, e), (2500, 2400));
    }

    #[test]
    fn input_order_independent() {
        let a = effective_terminal_boundaries(
            &positions(&[(600, 1), (1000, 5)]), &positions(&[(3000, 6)]),
            600, 3000, 2000, 2500, 0.10, 2.0,
        );
        let b = effective_terminal_boundaries(
            &positions(&[(1000, 5), (600, 1)]), &positions(&[(3000, 6)]),
            600, 3000, 2000, 2500, 0.10, 2.0,
        );
        assert_eq!(a, b);
    }
}
