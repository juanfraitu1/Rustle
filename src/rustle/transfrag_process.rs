//! Process transfrags: compatible_long merge, eliminate under threshold.
//!
//! ## Absorption
//! Absorption = a transfrag is merged into a _kept_ representative; its abundance is added to the
//! representative; the transfrag is marked absorbed (`weak = 1`), not kept as its own keeptrf entry.
//! The representative's `group` holds indices of transfrags absorbed into it; representative abundance
//! = sum of group (kept_cov).
//!
//! ## compatible_long return
//! - `ret == 0`: incompatible (do not merge).
//! - `ret == 1`: compatible, tf1 extends past tf2 (with introns).
//! - `ret == 2`: compatible, tf2 extends past tf1 (with introns).
//! - `ret == 3`: compatible, same or overlapping terminals (no intron-gap extension).
//! Incompatibility: different genes (if gene_id set), hard
//! boundary conflict, no overlap, containment, reachability, not continuous, middle gap filled, etc.
//!
//! ## Merge rules after compatibility
//! - **ret == 3**: Merge only if `lens[0] < edgedist` and `lens[2] < edgedist` (100). Opposing hard
//!   boundaries → keep both. Representative replacement: tf becomes rep if tf has hardstart+hardend or tf.guide.
//! - **ret == 1**: Absorb kept into tf (tf becomes rep) only if kept not guide, tf has longstart+longend,
//!   kept first/last not hard (or same node), `tf.abundance > DROP*kept.abundance`, `lens[1] < ssdist`, `lens[3] < ssdist`.
//! - **ret == 2**: Absorb tf into kept only if tf not guide, `lens[1] < ssdist`, `lens[3] < ssdist`, tf has no
//!   boundary evidence (or shares node with kept). Special: tf with longstart+longend and fewer nodes → separate keeptrf.
//!
//! ## Constants
//! - `ssdist` = 25: max intron-gap distance (lens[1], lens[3]) for ret=1/2 merge.
//! - `edgedist` = 100: max terminal distance (lens[0], lens[2]) for ret=3 merge.
//! - `DROP` = 0.5: for ret=1, tf must have abundance > DROP * kept.abundance to replace.
//! - `MAX_NODE`: unset lens[1]/lens[3] sentinel so they fail the ssdist check.

use crate::bitset::NodeSet;
use crate::graph::{FlowBranch, Graph, GraphNode, GraphTransfrag};
use crate::types::{DetHashMap as HashMap, DetHashSet as HashSet};
use std::fs::OpenOptions;
use std::io::{BufWriter, Write};
use std::sync::{Mutex, OnceLock};

/// Cross-bundle family-chain registry: intron-length-chain tuples from
/// transfrags that passed the strict boundary gate in earlier bundles.
/// Used by `RUSTLE_VG_FAMILY_RESCUE` to rescue single-boundary transfrags
/// whose intron chain matches a confirmed family member.
fn family_chain_registry() -> &'static Mutex<Vec<Vec<u32>>> {
    static REG: OnceLock<Mutex<Vec<Vec<u32>>>> = OnceLock::new();
    REG.get_or_init(|| Mutex::new(Vec::new()))
}

const FAMILY_CHAIN_MIN: usize = 5; // require >=5 introns to limit false matches

fn chain_register(chain: Vec<u32>) {
    if chain.len() < FAMILY_CHAIN_MIN {
        return;
    }
    let mut guard = family_chain_registry().lock().unwrap();
    if !guard.iter().any(|c| c == &chain) {
        guard.push(chain);
    }
}

/// Returns true if `candidate` is a contiguous subsequence of any registered
/// chain, with each intron length matching within `tol` bp.
fn chain_matches_registry(candidate: &[u32], tol: u32) -> bool {
    if candidate.len() < FAMILY_CHAIN_MIN {
        return false;
    }
    let guard = family_chain_registry().lock().unwrap();
    for reg in guard.iter() {
        if reg.len() < candidate.len() {
            continue;
        }
        for start in 0..=reg.len() - candidate.len() {
            let mut ok = true;
            for i in 0..candidate.len() {
                let a = candidate[i];
                let b = reg[start + i];
                let diff = if a > b { a - b } else { b - a };
                if diff > tol {
                    ok = false;
                    break;
                }
            }
            if ok {
                return true;
            }
        }
    }
    false
}

fn tf_intron_length_chain(tf: &GraphTransfrag, graph: &Graph) -> Vec<u32> {
    let mut lens: Vec<u32> = Vec::new();
    if tf.node_ids.len() < 2 {
        return lens;
    }
    for i in 0..tf.node_ids.len() - 1 {
        let Some(a) = graph.nodes.get(tf.node_ids[i]) else {
            continue;
        };
        let Some(b) = graph.nodes.get(tf.node_ids[i + 1]) else {
            continue;
        };
        // Only count true introns — gap between exons (b.start > a.end + 1).
        if b.start > a.end.saturating_add(1) {
            let ilen = (b.start - a.end - 1) as u32;
            lens.push(ilen);
        }
    }
    lens
}

/// Compute intron-length chain from a list of exon (start,end) tuples.
pub fn exons_intron_length_chain(exons: &[(u64, u64)]) -> Vec<u32> {
    let mut sorted = exons.to_vec();
    sorted.sort();
    let mut lens = Vec::new();
    for i in 0..sorted.len().saturating_sub(1) {
        let a_end = sorted[i].1;
        let b_start = sorted[i + 1].0;
        if b_start > a_end + 1 {
            lens.push((b_start - a_end - 1) as u32);
        }
    }
    lens
}

/// Find where `candidate` appears as a contiguous subsequence of `reference`
/// with each intron-length matching within `tol` bp. Returns the start offset
/// within `reference`, or None if not found.
pub fn chain_subsequence_offset(candidate: &[u32], reference: &[u32], tol: u32) -> Option<usize> {
    if candidate.is_empty() || candidate.len() > reference.len() {
        return None;
    }
    for start in 0..=reference.len() - candidate.len() {
        let mut ok = true;
        for i in 0..candidate.len() {
            let a = candidate[i];
            let b = reference[start + i];
            let diff = if a > b { a - b } else { b - a };
            if diff > tol {
                ok = false;
                break;
            }
        }
        if ok {
            return Some(start);
        }
    }
    None
}

/// Tolerant chain alignment — Needleman-Wunsch on intron-length sequences.
/// Returns the alignment score; higher (closer to 0) means more similar.
///
/// Scoring:
/// - Match: penalty `bp_penalty * |len_diff|` capped at `bp_cap`.
/// - Gap:   fixed penalty `gap_penalty` per skipped position.
///
/// Useful for detecting homology between distant paralogs whose splice graphs
/// have diverged past exact subsequence matching — the chr15 GOLGA6C cluster
/// (17-intron copies with per-position drift up to ~100 bp) scores around -5
/// to -640 pairwise, well above the -500 threshold a strict subseq misses.
pub fn chain_align_score(a: &[u32], b: &[u32]) -> f64 {
    let bp_penalty = 0.1f64;
    let gap_penalty = 500.0f64;
    let bp_cap = 2000.0f64;
    let n = a.len();
    let m = b.len();
    if n == 0 || m == 0 {
        return -(n.max(m) as f64 * gap_penalty);
    }
    let mut dp = vec![vec![0.0f64; m + 1]; n + 1];
    for i in 1..=n {
        dp[i][0] = -(i as f64) * gap_penalty;
    }
    for j in 1..=m {
        dp[0][j] = -(j as f64) * gap_penalty;
    }
    for i in 1..=n {
        for j in 1..=m {
            let diff = (a[i - 1] as i64 - b[j - 1] as i64).unsigned_abs() as f64;
            let pen = (diff * bp_penalty).min(bp_cap);
            let d = dp[i - 1][j - 1] - pen;
            let u = dp[i - 1][j] - gap_penalty;
            let l = dp[i][j - 1] - gap_penalty;
            dp[i][j] = d.max(u.max(l));
        }
    }
    dp[n][m]
}

/// Return true if candidate's chain aligns to any registered chain with
/// score ≥ threshold (more negative = looser). Used by the tolerant family
/// rescue path that catches distant sub-family members which strict
/// subsequence matching would miss.
pub fn chain_aligns_to_registry(candidate: &[u32], threshold: f64, min_len: usize) -> bool {
    if candidate.len() < min_len {
        return false;
    }
    let guard = family_chain_registry().lock().unwrap();
    for reg in guard.iter() {
        if reg.len() < min_len {
            continue;
        }
        if chain_align_score(candidate, reg) >= threshold {
            return true;
        }
    }
    false
}

const MAX_NODE: i64 = i64::MAX - 1; // unset lens sentinel
const SSDIST: i64 = 25;
const EDGEDIST: i64 = 100;
const DROP: f64 = 0.5; // abundance ratio for ret=1 replacement
const MAX_DIST: i64 = 200; // guide containment side check
const MAX_TRF_NUMBER: usize = 40000;

fn parse_trace_locus_env() -> Option<(u64, u64)> {
    let Ok(val) = std::env::var("RUSTLE_TRACE_LOCUS") else {
        return None;
    };
    let (a, b) = val.split_once('-')?;
    let start = a.trim().parse::<u64>().ok()?;
    let end = b.trim().parse::<u64>().ok()?;
    Some((start.min(end), start.max(end)))
}

fn parse_trace_intron_env() -> Option<(u64, u64)> {
    let Ok(val) = std::env::var("RUSTLE_TRACE_INTRON") else {
        return None;
    };
    let (a, b) = val.split_once('-')?;
    let donor = a.trim().parse::<u64>().ok()?;
    let acceptor = b.trim().parse::<u64>().ok()?;
    Some((donor, acceptor))
}

fn parse_trace_srcsink_nodes_env() -> Vec<usize> {
    let Ok(val) = std::env::var("RUSTLE_TRACE_SRCSINK_NODES") else {
        return Vec::new();
    };
    let mut nodes = Vec::new();
    for part in val.split(',') {
        let trimmed = part.trim();
        if trimmed.is_empty() {
            continue;
        }
        if let Ok(nid) = trimmed.parse::<usize>() {
            nodes.push(nid);
        }
    }
    nodes.sort_unstable();
    nodes.dedup();
    nodes
}

fn tf_overlaps_trace_graph(
    tf: &GraphTransfrag,
    graph: &Graph,
    trace_locus: Option<(u64, u64)>,
) -> bool {
    let Some((lo, hi)) = trace_locus else {
        return false;
    };
    tf.node_ids.iter().copied().any(|nid| {
        graph
            .nodes
            .get(nid)
            .map(|n| n.start <= hi && n.end >= lo)
            .unwrap_or(false)
    })
}

fn tf_span_graph(tf: &GraphTransfrag, graph: &Graph) -> (u64, u64) {
    let first = tf
        .node_ids
        .first()
        .and_then(|&nid| graph.nodes.get(nid))
        .map(|n| n.start)
        .unwrap_or(0);
    let last = tf
        .node_ids
        .last()
        .and_then(|&nid| graph.nodes.get(nid))
        .map(|n| n.end)
        .unwrap_or(0);
    (first, last)
}

fn tf_find_intron(
    tf: &GraphTransfrag,
    graph: &Graph,
    target: (u64, u64),
) -> Option<(usize, usize)> {
    for w in tf.node_ids.windows(2) {
        let a = w[0];
        let b = w[1];
        let Some(na) = graph.nodes.get(a) else {
            continue;
        };
        let Some(nb) = graph.nodes.get(b) else {
            continue;
        };
        if na.end < nb.start && (na.end, nb.start) == target {
            return Some((a, b));
        }
    }
    None
}

fn tf_junction_chain_coords(
    tf: &GraphTransfrag,
    graph: &Graph,
    source_id: usize,
    sink_id: usize,
) -> Vec<(u64, u64)> {
    let mut chain: Vec<(u64, u64)> = Vec::new();
    for w in tf.node_ids.windows(2) {
        let a = w[0];
        let b = w[1];
        if a == source_id || a == sink_id || b == source_id || b == sink_id {
            continue;
        }
        let Some(na) = graph.nodes.get(a) else {
            continue;
        };
        let Some(nb) = graph.nodes.get(b) else {
            continue;
        };
        if na.end < nb.start {
            chain.push((na.end, nb.start));
        }
    }
    chain
}

fn tf_has_novel_junction_vs_keeptrf(
    tf_chain: &[(u64, u64)],
    keeptrf: &[(usize, Vec<usize>, f64)],
    transfrags: &[GraphTransfrag],
    graph: &Graph,
    source_id: usize,
    sink_id: usize,
) -> bool {
    if tf_chain.is_empty() {
        return false;
    }
    let mut kept_junctions: HashSet<(u64, u64)> = Default::default();
    for (rep_idx, _, _) in keeptrf {
        let Some(rep) = transfrags.get(*rep_idx) else {
            continue;
        };
        for j in tf_junction_chain_coords(rep, graph, source_id, sink_id) {
            kept_junctions.insert(j);
        }
    }
    tf_chain.iter().any(|j| !kept_junctions.contains(j))
}

fn write_keeptrf_usepath_tsv(
    path: &str,
    keeptrf: &[(usize, Vec<usize>, f64)],
    trflong_insert: &[usize],
    transfrags: &[GraphTransfrag],
    graph: &Graph,
    mixed_mode: bool,
    keepsink: &[bool],
    hassink: &[Option<usize>],
) -> std::io::Result<()> {
    let existed = std::path::Path::new(path).exists();
    let file = OpenOptions::new().create(true).append(true).open(path)?;
    let write_header = if existed {
        file.metadata().map(|m| m.len() == 0).unwrap_or(true)
    } else {
        true
    };
    let mut writer = BufWriter::new(file);
    if write_header {
        writeln!(
            writer,
            "rep_idx\tgroup_cov\tgroup_size\tgroup_members\ttrflong_insert_rank\teffective_parse_rank\tseed_complete\tguide\tabundance\tread_count\tlongstart\tlongend\thardstart\thardend\thas_source_keep\thas_sink_keep\tweak\treal\ttrflong_seed\tusepath\tstart\tend\tnode_count\tjunction_chain\tnode_ids"
        )?;
    }

    let source_id = graph.source_id;
    let sink_id = graph.sink_id;
    let mut parse_order: Vec<usize> = transfrags
        .iter()
        .enumerate()
        .filter(|(_, tf)| {
            tf.trflong_seed
                && tf.weak == 0
                && !tf.node_ids.is_empty()
                && if mixed_mode { tf.usepath < -1 } else { tf.usepath >= 0 }
        })
        .map(|(idx, _)| idx)
        .collect();
    parse_order.sort_unstable_by_key(|&idx| transfrags[idx].usepath);
    let parse_rank: HashMap<usize, usize> = parse_order
        .into_iter()
        .enumerate()
        .map(|(rank, rep_idx)| (rep_idx, rank))
        .collect();
    let insert_rank: HashMap<usize, usize> = trflong_insert
        .iter()
        .copied()
        .enumerate()
        .map(|(rank, rep_idx)| (rep_idx, rank))
        .collect();

    for (rep_idx, group, group_cov) in keeptrf {
        let Some(tf) = transfrags.get(*rep_idx) else {
            continue;
        };
        if tf.node_ids.is_empty() {
            continue;
        }

        let first_real = tf
            .node_ids
            .iter()
            .copied()
            .find(|&nid| nid != source_id && nid != sink_id);
        let last_real = tf
            .node_ids
            .iter()
            .rev()
            .copied()
            .find(|&nid| nid != source_id && nid != sink_id);
        let start_node = first_real.and_then(|nid| graph.nodes.get(nid));
        let end_node = last_real.and_then(|nid| graph.nodes.get(nid));
        let start = start_node.map(|n| n.start).unwrap_or(0);
        let end = end_node.map(|n| n.end).unwrap_or(start);
        let hardstart = start_node.map(|n| n.hardstart).unwrap_or(false);
        let hardend = end_node.map(|n| n.hardend).unwrap_or(false);
        let has_source_keep = first_real
            .and_then(|nid| graph.nodes.get(source_id).map(|src| src.children.contains(nid)))
            .unwrap_or(false);
        let has_sink_keep = last_real
            .map(|nid| has_sink_completion(nid, graph, keepsink, hassink))
            .unwrap_or(false);
        let seed_complete =
            tf.guide || ((hardstart || has_source_keep) && (hardend || has_sink_keep));
        let junction_chain = tf
            .node_ids
            .iter()
            .copied()
            .filter(|&nid| nid != source_id && nid != sink_id)
            .collect::<Vec<_>>()
            .windows(2)
            .filter_map(|w| {
                let left = graph.nodes.get(w[0])?;
                let right = graph.nodes.get(w[1])?;
                if left.end < right.start {
                    Some(format!("{}-{}", left.end, right.start))
                } else {
                    None
                }
            })
            .collect::<Vec<_>>()
            .join(",");
        let node_ids = tf
            .node_ids
            .iter()
            .map(|nid| nid.to_string())
            .collect::<Vec<_>>()
            .join(",");
        let group_members = group
            .iter()
            .map(|idx| idx.to_string())
            .collect::<Vec<_>>()
            .join(",");
        writeln!(
            writer,
            "{}\t{:.4}\t{}\t{}\t{}\t{}\t{}\t{}\t{:.4}\t{:.4}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
            rep_idx,
            *group_cov,
            group.len(),
            group_members,
            insert_rank
                .get(rep_idx)
                .map(|v| v.to_string())
                .unwrap_or_default(),
            parse_rank
                .get(rep_idx)
                .map(|v| v.to_string())
                .unwrap_or_default(),
            seed_complete as u8,
            tf.guide as u8,
            tf.abundance,
            tf.read_count,
            tf.longstart,
            tf.longend,
            hardstart as u8,
            hardend as u8,
            has_source_keep as u8,
            has_sink_keep as u8,
            tf.weak,
            tf.real as u8,
            tf.trflong_seed as u8,
            tf.usepath,
            start,
            end,
            tf.node_ids.len(),
            junction_chain,
            node_ids,
        )?;
    }
    writer.flush()?;
    Ok(())
}

fn has_sink_completion(
    start_nid: usize,
    graph: &Graph,
    keepsink: &[bool],
    hassink: &[Option<usize>],
) -> bool {
    let sink_id = graph.sink_id;
    let mut current = start_nid;
    let mut steps = 0usize;
    let mut seen: HashSet<usize> = Default::default();
    loop {
        if current >= graph.nodes.len() || !seen.insert(current) {
            return false;
        }
        let node = &graph.nodes[current];
        if node.hardend {
            return true;
        }
        if current < hassink.len() && current < keepsink.len() && hassink[current].is_some() && keepsink[current] {
            return true;
        }
        if node.children.contains(sink_id) {
            return true;
        }
        let real_children: Vec<usize> = node
            .children
            .ones()
            .filter(|&cid| cid != sink_id)
            .collect();
        if real_children.len() != 1 {
            return false;
        }
        current = real_children[0];
        steps += 1;
        if steps > 64 {
            return false;
        }
    }
}

/// Promote terminal split-tail evidence into full-length long-read seed variants.
///
/// Rust already materializes 2-node synthetic tail links like `[82, 255]`, but some long-read
/// transfrags that end at the parent node keep only `... -> 82`. can still recover the
/// full terminal seed via later seed/path handling. To avoid losing that family before seed
/// extraction, synthesize a companion transfrag `[..., 82, 255]` when:
/// - the transfrag is long-read,
/// - it ends at a node with a unique non-sink contiguous 2-node long tail link,
/// - and that exact extended path does not already exist.
fn promote_terminal_tail_variants(
    transfrags: &mut Vec<GraphTransfrag>,
    graph: &mut Graph,
    trace_locus: Option<(u64, u64)>,
) {
    let source_id = graph.source_id;
    let sink_id = graph.sink_id;

    let mut tail_children: HashMap<usize, Vec<(usize, f64)>> = Default::default();
    for tf in transfrags.iter() {
        if !tf.longread || tf.node_ids.len() != 2 {
            continue;
        }
        let from = tf.node_ids[0];
        let to = tf.node_ids[1];
        if from == source_id || from == sink_id || to == source_id || to == sink_id {
            continue;
        }
        let Some(from_node) = graph.nodes.get(from) else {
            continue;
        };
        let Some(to_node) = graph.nodes.get(to) else {
            continue;
        };
        if from_node.end != to_node.start {
            continue;
        }
        tail_children
            .entry(from)
            .or_default()
            .push((to, tf.abundance.max(1.0)));
    }

    if tail_children.is_empty() {
        return;
    }

    // Only propagate a tail link when there is real long-read endpoint evidence that at least
    // one transfrag ending at this node extends beyond the current terminal node boundary.
    let mut family_has_tail_support: HashSet<usize> = Default::default();
    for tf in transfrags.iter() {
        if !tf.longread || tf.node_ids.len() < 2 {
            continue;
        }
        let Some(&last) = tf.node_ids.last() else {
            continue;
        };
        let Some(last_node) = graph.nodes.get(last) else {
            continue;
        };
        if tf.longend > last_node.end {
            family_has_tail_support.insert(last);
        }
    }

    let mut existing_paths: HashSet<Vec<usize>> =
        transfrags.iter().map(|tf| tf.node_ids.clone()).collect();
    let base_len = transfrags.len();
    let psize = graph.pattern_size();
    let mut picked: HashMap<(usize, usize, usize), usize> = Default::default();

    // Pick at most one base per (first node, last node, tail child): prefer the shortest path
    // (most boundary-specific) and then the highest abundance.
    for idx in 0..base_len {
        let tf = &transfrags[idx];
        if !tf.longread || tf.guide || tf.node_ids.len() < 2 || tf.node_ids.len() == 2 {
            continue;
        }
        let first = tf.node_ids[0];
        let Some(first_node) = graph.nodes.get(first) else {
            continue;
        };
        let first_has_source = first_node.parents.contains(source_id);
        let first_has_boundary = first_node.hardstart || first_has_source;
        if !first_has_boundary {
            continue;
        }
        let Some(&last) = tf.node_ids.last() else {
            continue;
        };
        if !family_has_tail_support.contains(&last) {
            continue;
        }
        let Some(children) = tail_children.get(&last) else {
            continue;
        };
        let Some(last_node) = graph.nodes.get(last) else {
            continue;
        };
        if !last_node.hardend && tf.longend == 0 {
            continue;
        }
        for &(child, _) in children {
            if tf.node_ids.contains(&child) {
                continue;
            }
            // Keep separate hardstart families, but collapse non-hardstart source-linked
            // internal starts to a single earliest-start representative per tail.
            let start_bucket = if first_node.hardstart {
                first
            } else {
                usize::MAX
            };
            let key = (last, start_bucket, child);
            match picked.get(&key).copied() {
                None => {
                    picked.insert(key, idx);
                }
                Some(prev_idx) => {
                    let prev = &transfrags[prev_idx];
                    let prev_first = prev.node_ids[0];
                    let prev_first_start = graph
                        .nodes
                        .get(prev_first)
                        .map(|n| n.start)
                        .unwrap_or(u64::MAX);
                    let this_first_start = first_node.start;
                    let better = this_first_start < prev_first_start
                        || (this_first_start == prev_first_start
                            && tf.node_ids.len() < prev.node_ids.len())
                        || (this_first_start == prev_first_start
                            && tf.node_ids.len() == prev.node_ids.len()
                            && tf.abundance > prev.abundance);
                    if better {
                        picked.insert(key, idx);
                    }
                }
            }
        }
    }

    if picked.is_empty() {
        if trace_locus.is_some() {
            eprintln!(
                "[TRACE_TAIL_PROMOTE] no_picks tail_bases={} tail_families={}",
                tail_children.len(),
                family_has_tail_support.len()
            );
        }
        return;
    }

    let mut added: Vec<GraphTransfrag> = Vec::new();
    let mut pick_list: Vec<((usize, usize, usize), usize)> = picked.into_iter().collect();
    pick_list.sort_unstable_by_key(|(_, idx)| *idx);
    for ((_, _, child), idx) in pick_list {
        let tf = &transfrags[idx];
        let Some(&last) = tf.node_ids.last() else {
            continue;
        };
        let Some(children) = tail_children.get(&last) else {
            continue;
        };
        let Some((_, tail_ab)) = children.iter().copied().find(|(c, _)| *c == child) else {
            if tf_overlaps_trace_graph(tf, graph, trace_locus) {
                eprintln!(
                    "[TRACE_TAIL_PROMOTE] skip_missing_tail_child base={} last={} want_child={}",
                    idx, last, child
                );
            }
            continue;
        };
        let mut new_nodes = tf.node_ids.clone();
        new_nodes.push(child);
        if !existing_paths.insert(new_nodes.clone()) {
            continue;
        }

        let mut ext = GraphTransfrag::new(new_nodes.clone(), psize);
        graph.set_pattern_edges_for_path(&mut ext.pattern, &new_nodes);
        ext.longread = true;
        ext.abundance = tail_ab.min(tf.abundance.max(1.0));
        ext.read_count = ext.abundance;
        ext.longstart = tf.longstart;
        ext.longend = graph.nodes.get(child).map(|n| n.end).unwrap_or(tf.longend);
        ext.poly_start_unaligned = tf.poly_start_unaligned;
        ext.poly_end_unaligned = tf.poly_end_unaligned;
        ext.killed_junction_orphan = tf.killed_junction_orphan;

        if tf_overlaps_trace_graph(tf, graph, trace_locus) {
            let (old_s, old_e) = tf_span_graph(tf, graph);
            let (new_s, new_e) = tf_span_graph(&ext, graph);
            eprintln!(
                "[TRACE_TAIL_PROMOTE] base={} span={}-{} -> extended span={}-{} add_path={:?} tail_ab={:.3}",
                idx,
                old_s,
                old_e,
                new_s,
                new_e,
                new_nodes,
                ext.abundance
            );
        }
        added.push(ext);
    }

    if added.is_empty() {
        return;
    }

    let start_idx = transfrags.len();
    transfrags.extend(added);
    for idx in start_idx..transfrags.len() {
        if transfrags[idx].node_ids.len() <= 1 {
            continue;
        }
        for &nid in &transfrags[idx].node_ids {
            if let Some(node) = graph.nodes.get_mut(nid) {
                if !node.trf_ids.contains(&idx) {
                    node.trf_ids.push(idx);
                }
            }
        }
    }
}

fn trace_tf_pattern_outgoing(label: &str, tf_idx: usize, tf: &GraphTransfrag, graph: &Graph) {
    let Some(&pivot) = tf.node_ids.last() else {
        return;
    };
    if !graph
        .nodes
        .get(pivot)
        .map(|n| trace_locus_active_span(n.start, n.end, parse_trace_locus_env()))
        .unwrap_or(false)
    {
        return;
    }
    let Some(node) = graph.nodes.get(pivot) else {
        return;
    };
    eprintln!(
        "[TRACE_PAT_TF] stage={} tf={} pivot={} coord={}-{} pop={} group={} longstart={} longend={}",
        label,
        tf_idx,
        pivot,
        node.start,
        node.end,
        tf.pattern.len_bits(),
        tf.group.len(),
        tf.longstart,
        tf.longend
    );
    for child in node.children.ones() {
        let ccoord = graph
            .nodes
            .get(child)
            .map(|n| (n.start, n.end))
            .unwrap_or((0, 0));
        let edge_on = graph
            .edge_bit_index(pivot, child)
            .map(|eid| tf.pattern.get_bit(eid))
            .unwrap_or(false);
        eprintln!(
            "[TRACE_PAT_TF]   child={} coord={}-{} edge_bit={} child_bit={} sink={}",
            child,
            ccoord.0,
            ccoord.1,
            edge_on,
            tf.pattern.get_bit(child),
            child == graph.sink_id
        );
    }
}

fn rebuild_node_trf_ids(graph: &mut Graph, transfrags: &[GraphTransfrag]) {
    let mut counts = vec![0usize; graph.nodes.len()];
    for tf in transfrags {
        if tf.node_ids.len() <= 1 {
            continue;
        }
        for &nid in &tf.node_ids {
            if nid < counts.len() {
                counts[nid] += 1;
            }
        }
    }

    for (nid, node) in graph.nodes.iter_mut().enumerate() {
        node.trf_ids.clear();
        let needed = counts[nid];
        if needed > node.trf_ids.capacity() {
            node.trf_ids.reserve(needed - node.trf_ids.capacity());
        }
    }

    for (tf_idx, tf) in transfrags.iter().enumerate() {
        if tf.node_ids.len() <= 1 {
            continue;
        }
        for &nid in &tf.node_ids {
            if let Some(node) = graph.nodes.get_mut(nid) {
                node.trf_ids.push(tf_idx);
            }
        }
    }
}

fn collect_active_transfrag_indices(
    transfrags: &[GraphTransfrag],
    source_id: usize,
    sink_id: usize,
) -> Vec<usize> {
    let mut active = Vec::with_capacity(transfrags.len());
    for (i, tf) in transfrags.iter().enumerate() {
        if tf.node_ids.is_empty() {
            continue;
        }
        if tf.node_ids[0] != source_id && tf.node_ids.last().copied() != Some(sink_id) {
            active.push(i);
        }
    }
    active
}

fn trace_locus_active_span(start: u64, end: u64, trace_locus: Option<(u64, u64)>) -> bool {
    let Some((lo, hi)) = trace_locus else {
        return false;
    };
    end >= lo && start <= hi
}

thread_local! {
    /// Last incompat-return site code inside compatible_long. Callers read this
    /// after a ret==0 to discover which gate fired. Paired with StringTie's
    /// PARITY_CL_TRACE to diff gate behavior per transfrag pair.
    pub static LAST_CL_INCOMPAT_SITE: std::cell::Cell<&'static str> =
        const { std::cell::Cell::new("") };
}

#[inline]
fn set_cl_site(site: &'static str) {
    LAST_CL_INCOMPAT_SITE.with(|c| c.set(site));
}

/// compatible_long return: 0 = incompatible, 1 = tf1 extends past tf2, 2 = tf2 extends past tf1, 3 = same/compatible.
/// lens: [len0, len1, len2, len3] for distance checks (len[]).
pub fn compatible_long(
    tf1: &GraphTransfrag,
    tf2: &GraphTransfrag,
    graph: &Graph,
) -> (u8, [i64; 4]) {
    let incompat = (0, [0i64, MAX_NODE, 0, MAX_NODE]);
    let n0 = &tf1.node_ids;
    let n1 = &tf2.node_ids;
    if n0.is_empty() || n1.is_empty() {
        set_cl_site("S0_EMPTY");
        return incompat;
    }

    let nodes = [n0, n1];
    let get_node = |idx: usize, i: usize| -> &GraphNode { &graph.nodes[nodes[idx][i]] };

    // Both starts can't be hard at different nodes
    if nodes[0][0] != nodes[1][0] {
        if get_node(0, 0).hardstart && get_node(1, 0).hardstart {
            set_cl_site("S1_HARDSTART_BOTH");
            return incompat;
        }
    }
    if nodes[0][nodes[0].len() - 1] != nodes[1][nodes[1].len() - 1] {
        if get_node(0, nodes[0].len() - 1).hardend && get_node(1, nodes[1].len() - 1).hardend {
            set_cl_site("S2_HARDEND_BOTH");
            return incompat;
        }
    }

    let tstart = [
        if tf1.longstart != 0 {
            tf1.longstart
        } else {
            get_node(0, 0).start
        },
        if tf2.longstart != 0 {
            tf2.longstart
        } else {
            get_node(1, 0).start
        },
    ];
    let tend = [
        if tf1.longend != 0 {
            tf1.longend
        } else {
            get_node(0, nodes[0].len() - 1).end
        },
        if tf2.longend != 0 {
            tf2.longend
        } else {
            get_node(1, nodes[1].len() - 1).end
        },
    ];

    let mut f = 0usize;
    let mut s = 1usize;
    if tstart[1] < tstart[0] {
        f = 1;
        s = 0;
    }
    if nodes[f][nodes[f].len() - 1] < nodes[s][0] {
        set_cl_site("S3_NO_OVERLAP");
        return incompat;
    }

    // compatible_long zeroes len[0..3] internally for compatible cases.
    // MAX_NODE is only a caller-side sentinel before the function is invoked.
    let mut lens = [0i64, 0, 0, 0];
    let mut rets = 3u8;
    let mut i = [0usize, 0usize];

    // START compatibility
    while i[f] < nodes[f].len() && nodes[f][i[f]] < nodes[s][0] {
        if i[f] + 1 < nodes[f].len() {
            let curr_end = get_node(f, i[f]).end;
            let next_start = get_node(f, i[f] + 1).start;
            if curr_end < next_start {
                rets = 1 + if f == 0 { 0 } else { 1 };
                if get_node(s, 0).hardstart {
                    set_cl_site("S4_HARDSTART_S");
                    return incompat;
                }
            }
        }
        if i[f] == 0 {
            lens[0] = (get_node(f, i[f]).end.saturating_sub(tstart[f])) as i64;
        } else {
            lens[0] += get_node(f, i[f]).length() as i64;
        }
        i[f] += 1;
    }
    if i[f] < nodes[f].len() && nodes[s][nodes[s].len() - 1] < nodes[f][i[f]] {
        set_cl_site("S5_CONTAINED");
        return incompat;
    }
    if i[f] == 0 {
        lens[0] = tstart[s] as i64 - tstart[f] as i64;
    }
    if i[f] < nodes[f].len() && nodes[s][0] < nodes[f][i[f]] {
        if !graph.can_reach(nodes[s][0], nodes[f][i[f]]) {
            set_cl_site("S6_CHILDPAT_START");
            return incompat;
        }
        lens[1] = (get_node(f, i[f]).start.saturating_sub(tstart[s])) as i64;
        while i[s] < nodes[s].len() && nodes[s][i[s]] < nodes[f][i[f]] {
            i[s] += 1;
            if i[s] < nodes[s].len() {
                let prev_end = get_node(s, i[s] - 1).end;
                let curr_start = get_node(s, i[s]).start;
                if prev_end < curr_start {
                    set_cl_site("S7_WALK_START_GAP");
                    return incompat;
                }
            }
        }
    } else if i[f] > 0 {
        lens[0] += (tstart[s].saturating_sub(get_node(s, 0).start)) as i64;
    }
    if f == 1 {
        lens[0] = -lens[0];
    }

    // END compatibility
    let mut rete = 3u8;
    f = 0;
    s = 1;
    if tend[1] > tend[0] {
        f = 1;
        s = 0;
    }
    let mut j = [nodes[0].len() - 1, nodes[1].len() - 1];

    while j[f] < nodes[f].len() && nodes[f][j[f]] > nodes[s][nodes[s].len() - 1] {
        if j[f] > 0 {
            let prev_end = get_node(f, j[f] - 1).end;
            let curr_start = get_node(f, j[f]).start;
            if prev_end < curr_start {
                rete = 1 + if f == 0 { 0 } else { 1 };
                if rets != 3 && rete != rets {
                    set_cl_site("S8_RETE_CONTRADICT");
                    return incompat;
                }
                if get_node(s, nodes[s].len() - 1).hardend {
                    set_cl_site("S9_HARDEND_S");
                    return incompat;
                }
            }
        }
        if j[f] < nodes[f].len() - 1 {
            lens[2] += get_node(f, j[f]).length() as i64;
        } else {
            lens[2] += (tend[f].saturating_sub(get_node(f, j[f]).start)) as i64;
        }
        if j[f] == 0 {
            break;
        }
        j[f] -= 1;
    }
    if j[f] < nodes[f].len() && nodes[s][0] > nodes[f][j[f]] {
        set_cl_site("S10_NO_OVERLAP_END");
        return incompat;
    }
    if j[f] == nodes[f].len() - 1 {
        lens[2] += (tend[f].saturating_sub(tend[s])) as i64;
    }
    if j[f] < nodes[f].len() && nodes[s][j[s]] > nodes[f][j[f]] {
        if !graph.can_reach(nodes[f][j[f]], nodes[s][j[s]]) {
            set_cl_site("S11_CHILDPAT_END");
            return incompat;
        }
        lens[3] = (tend[s] - get_node(f, j[f]).end) as i64;
        while j[s] < nodes[s].len() && nodes[s][j[s]] > nodes[f][j[f]] {
            if j[s] == 0 {
                break;
            }
            j[s] -= 1;
            if j[s] < nodes[s].len() {
                let curr_end = get_node(s, j[s]).end;
                let next_start = get_node(s, j[s] + 1).start;
                if curr_end < next_start {
                    set_cl_site("S12_WALK_END_GAP");
                    return incompat;
                }
            }
        }
    } else if j[f] < nodes[f].len() - 1 {
        let last_end = get_node(s, nodes[s].len() - 1).end;
        lens[2] += last_end.saturating_sub(tend[s]) as i64;
    }
    if i[0] > j[0] || i[1] > j[1] {
        set_cl_site("S13_I_GT_J");
        return incompat;
    }
    if f == 1 {
        lens[2] = -lens[2];
    }

    // Middle compatibility with hard-edge pattern check (strict 5738-5771 port)
    let mut fi = f;
    let mut si = s;
    let mut ii = [i[0], i[1]];
    let jj = [j[0], j[1]];
    let patterns = [&tf1.pattern, &tf2.pattern];
    while ii[fi] <= jj[fi] && ii[si] <= jj[si] {
        if nodes[fi][ii[fi]] == nodes[si][ii[si]] {
            ii[fi] += 1;
            ii[si] += 1;
        } else {
            if nodes[fi][ii[fi]] > nodes[si][ii[si]] {
                std::mem::swap(&mut fi, &mut si);
            }
            if ii[fi] == 0 || ii[si] == 0 {
                set_cl_site("S13b_MID_IDX_ZERO");
                return incompat;
            }
            while ii[fi] > 0
                && ii[fi] <= jj[fi]
                && ii[fi] < nodes[fi].len()
                && nodes[fi][ii[fi]] < nodes[si][ii[si]]
                && get_node(fi, ii[fi] - 1).end == get_node(fi, ii[fi]).start
            {
                ii[fi] += 1;
            }
            if ii[fi] <= jj[fi]
                && ii[fi] < nodes[fi].len()
                && ii[fi] > 0
                && get_node(fi, ii[fi] - 1).end == get_node(fi, ii[fi]).start
            {
                set_cl_site("S14_GAP_FILLED");
                return incompat; // gap filled
            }
            // Hard-edge pattern check (5756-5759):
            // if gap is hard in both transcripts, they are incompatible.
            if ii[fi] < nodes[fi].len() && ii[si] < nodes[si].len() {
                let f_prev = nodes[fi][ii[fi] - 1];
                let f_curr = nodes[fi][ii[fi]];
                let edge_key_f = (f_prev.min(f_curr), f_prev.max(f_curr));
                if let Some(&eid_f) = graph.gpos.get(&edge_key_f) {
                    if patterns[fi].get_bit(eid_f) {
                        let s_prev = nodes[si][ii[si] - 1];
                        let s_curr = nodes[si][ii[si]];
                        let edge_key_s = (s_prev.min(s_curr), s_prev.max(s_curr));
                        if let Some(&eid_s) = graph.gpos.get(&edge_key_s) {
                            if patterns[si].get_bit(eid_s) {
                                set_cl_site("S15_HARD_EDGE_BOTH");
                                return incompat;
                            }
                        }
                    }
                }
            }
            if ii[fi] > jj[fi] || ii[fi] >= nodes[fi].len() {
                break;
            }
            if nodes[fi][ii[fi]] > nodes[si][ii[si]] {
                std::mem::swap(&mut fi, &mut si);
            }
            while ii[fi] <= jj[fi]
                && ii[fi] < nodes[fi].len()
                && nodes[fi][ii[fi]] < nodes[si][ii[si]]
            {
                if ii[fi] + 1 < nodes[fi].len() {
                    if get_node(fi, ii[fi]).end < get_node(fi, ii[fi] + 1).start {
                        set_cl_site("S16_MID_WALK_GAP");
                        return incompat;
                    }
                }
                ii[fi] += 1;
            }
        }
    }

    if rets == 1 || rete == 1 {
        return (1, lens);
    }
    if rets == 2 || rete == 2 {
        return (2, lens);
    }
    (3, lens)
}

fn assign_incomplete_trf_to_nodes(tf_idx: usize, n1: usize, n2: usize, graph: &mut Graph) -> bool {
    if n2 <= n1 + 1 {
        return false;
    }
    let mut assigned = false;
    for i in (n1 + 1)..n2 {
        if graph.can_reach(n1, i) && graph.can_reach(i, n2) {
            if let Some(node) = graph.nodes.get_mut(i) {
                if !node.trf_ids.contains(&tf_idx) {
                    node.trf_ids.push(tf_idx);
                    node.trf_ids.sort_unstable();
                    assigned = true;
                }
            }
        }
    }
    assigned
}

fn trf_real(tf_idx: usize, graph: &Graph, transfrags: &mut [GraphTransfrag]) -> bool {
    let Some(tf_ro) = transfrags.get(tf_idx) else {
        return true;
    };
    if tf_ro.node_ids.len() < 2 {
        return true;
    }
    let tf_nodes = tf_ro.node_ids.clone();
    let tf_pattern = tf_ro.pattern.clone();
    let nt = tf_nodes.len() - 1;
    let mut i = 0usize;
    let mut n = tf_nodes[i];
    let mut totalcov = 0.0f64;
    let mut built_paths: Vec<FlowBranch> = Vec::new();
    while i < nt {
        let nextnode = tf_nodes[i + 1];
        if n + 1 == nextnode {
            n = nextnode;
            i += 1;
            continue;
        }
        if let Some(epos) = graph.edge_bit_index(n, nextnode) {
            if tf_pattern.get_bit(epos) {
                n = nextnode;
                i += 1;
                continue;
            }
        }
        let Some(node) = graph.nodes.get(n) else {
            return true;
        };
        let mut np = 0usize;
        totalcov = 0.0;
        for child in node.children.ones() {
            if child > nextnode {
                break;
            }
            if child == nextnode || graph.can_reach(child, nextnode) {
                let mut pab = 0.0;
                if let Some(epos) = graph.edge_bit_index(n, child) {
                    for &tid in &node.trf_ids {
                        if let Some(t2) = transfrags.get(tid) {
                            if t2.pattern.get_bit(epos) {
                                pab += t2.abundance;
                                totalcov += t2.abundance;
                            }
                        }
                    }
                }
                if pab > 0.0 {
                    built_paths.push(FlowBranch {
                        node: n,
                        contnode: child,
                        abundance: pab,
                    });
                    np += 1;
                }
            }
        }
        if totalcov <= 0.0 || np == 0 {
            return true;
        }
        if np == 1 {
            n = nextnode;
            i += 1;
            built_paths.clear();
        } else {
            break;
        }
    }
    if i == nt {
        if let Some(tf) = transfrags.get_mut(tf_idx) {
            tf.flow_paths.clear();
            tf.flow_path_idx = -1;
        }
        true
    } else {
        if totalcov > 0.0 {
            for p in &mut built_paths {
                p.abundance /= totalcov;
            }
        }
        if let Some(tf) = transfrags.get_mut(tf_idx) {
            tf.flow_paths = built_paths;
            tf.flow_path_idx = -1;
        }
        false
    }
}

/// Process one short-read fragment: redistribute srabund to compatible transfrags (process_srfrag; mixed mode).
/// u_idx indexes the transfrag with srabund to redistribute. Modifies transfrags in place.
pub fn process_srfrag(
    transfrags: &mut [GraphTransfrag],
    u_idx: usize,
    graph: &Graph,
    drop: f64,
    error_perc: f64,
    epsilon: f64,
) {
    let n_nodes = graph.nodes.len();
    let u = match transfrags.get(u_idx) {
        Some(t) => t,
        None => return,
    };
    let available = u.srabund * drop * error_perc - u.abundance;
    if available <= epsilon {
        return;
    }
    let n0 = *u.node_ids.first().unwrap_or(&0);
    if n0 >= n_nodes {
        return;
    }
    let u_last = *u.node_ids.last().unwrap_or(&0);
    let node0 = &graph.nodes[n0];
    let u_node_ids = u.node_ids.clone();
    let u_pattern = u.pattern.clone();

    let mut available = available;
    let mut equalinc: Vec<usize> = Vec::new();
    let mut equalext: Vec<usize> = Vec::new();

    for &t_idx in &node0.trf_ids {
        if t_idx >= transfrags.len() || t_idx == u_idx {
            continue;
        }
        let t = &transfrags[t_idx];
        let t_first = *t.node_ids.first().unwrap_or(&0);
        let t_last = *t.node_ids.last().unwrap_or(&0);
        if t_first > n0 || t_last < u_last {
            continue;
        }
        if t.pattern.contains_pattern(&u_pattern) {
            available -= t.abundance;
            if available <= epsilon {
                available = 0.0;
                break;
            }
        } else {
            if t_first == u.node_ids[0] && t_last == u_last {
                equalinc.push(t_idx);
            } else {
                equalext.push(t_idx);
            }
        }
    }
    if available <= epsilon {
        return;
    }

    let mut add_to_u = 0.0;
    for &t_idx in &equalinc {
        if available <= epsilon {
            break;
        }
        let solveabund = transfrags[t_idx].abundance.min(available);
        transfrags[t_idx].abundance -= solveabund;
        transfrags[t_idx].srabund -= solveabund;
        if transfrags[t_idx].abundance < epsilon {
            transfrags[t_idx].abundance = 0.0;
        }
        add_to_u += solveabund;
        available -= solveabund;
    }
    transfrags[u_idx].abundance += add_to_u;

    if available <= epsilon {
        return;
    }
    for &t_idx in &equalext {
        if available <= epsilon {
            break;
        }
        let t = &mut transfrags[t_idx];
        let mut merged: Vec<usize> = t.node_ids.clone();
        for &nid in &u_node_ids {
            if let Err(insert_at) = merged.binary_search(&nid) {
                merged.insert(insert_at, nid);
            }
        }
        t.node_ids = merged;
        t.node_id_set = NodeSet::from_iter(t.node_ids.iter().copied());
        t.pattern.or_assign(&u_pattern);
        available -= t.abundance;
    }

    if available <= epsilon {
        return;
    }

    // stitch_trf: assemble a chain of compatible included transfrags
    // to explain remaining sr abundance.
    let mut seltrfrag: Vec<usize> = Vec::new();
    for (t_idx, t) in transfrags.iter().enumerate() {
        if t_idx == u_idx || t.abundance <= epsilon {
            continue;
        }
        if t.pattern.contains_pattern(&u_pattern) {
            continue;
        }
        if u_pattern.contains_pattern(&t.pattern) {
            seltrfrag.push(t_idx);
        }
    }
    for &nid in u_node_ids.iter().skip(1) {
        if let Some(n) = graph.nodes.get(nid) {
            for &t_idx in &n.trf_ids {
                if t_idx >= transfrags.len() || t_idx == u_idx {
                    continue;
                }
                let t = &transfrags[t_idx];
                if t.abundance <= epsilon {
                    continue;
                }
                if t.node_ids.first().copied() == Some(nid)
                    && u_pattern.contains_pattern(&t.pattern)
                {
                    seltrfrag.push(t_idx);
                }
            }
        }
    }
    seltrfrag.sort_unstable();
    seltrfrag.dedup();
    seltrfrag.sort_unstable_by(|&a, &b| {
        transfrags[b]
            .abundance
            .partial_cmp(&transfrags[a].abundance)
            .unwrap_or(std::cmp::Ordering::Equal)
    });
    if seltrfrag.is_empty() {
        return;
    }

    let mut start_map: HashMap<usize, Vec<usize>> = Default::default();
    for (i, &tid) in seltrfrag.iter().enumerate() {
        let start = transfrags[tid].node_ids.first().copied().unwrap_or(0);
        start_map.entry(start).or_default().push(i);
    }

    fn stitch_trf_path(
        u_first: usize,
        u_last: usize,
        available: f64,
        i_start: usize,
        i_cur: usize,
        seltrfrag: &[usize],
        transfrags: &[GraphTransfrag],
        start_map: &HashMap<usize, Vec<usize>>,
        depth: usize,
    ) -> Option<(f64, Vec<usize>)> {
        if depth > 128 {
            return None;
        }
        let t = &transfrags[seltrfrag[i_cur]];
        let avail = available.min(t.abundance);
        if avail <= 0.0 {
            return None;
        }
        let mut nl = t.node_ids.last().copied().unwrap_or(u_first);
        if nl == u_last {
            nl = u_first;
        }
        let start0 = transfrags[seltrfrag[i_start]]
            .node_ids
            .first()
            .copied()
            .unwrap_or(u_first);
        if nl == start0 {
            return Some((avail, vec![i_cur]));
        }
        let nf = t.node_ids.first().copied().unwrap_or(u_first);
        if nf < start0 && nl > start0 {
            return None;
        }
        let Some(nexts) = start_map.get(&nl) else {
            return None;
        };
        for &j in nexts {
            if j <= i_start {
                continue;
            }
            if transfrags[seltrfrag[j]].abundance <= 0.0 {
                continue;
            }
            if let Some((maxab, mut path)) = stitch_trf_path(
                u_first,
                u_last,
                avail,
                i_start,
                j,
                seltrfrag,
                transfrags,
                start_map,
                depth + 1,
            ) {
                path.push(i_cur);
                return Some((maxab, path));
            }
        }
        None
    }

    while available > epsilon {
        let mut progressed = false;
        for i in 0..seltrfrag.len() {
            let tid = seltrfrag[i];
            if transfrags[tid].abundance <= epsilon {
                continue;
            }
            if let Some((maxabund, path_idx)) = stitch_trf_path(
                n0, u_last, available, i, i, &seltrfrag, transfrags, &start_map, 0,
            ) {
                if maxabund <= epsilon {
                    continue;
                }
                for pidx in path_idx {
                    let pt = seltrfrag[pidx];
                    transfrags[pt].abundance -= maxabund;
                    transfrags[pt].srabund -= maxabund;
                    if transfrags[pt].abundance < epsilon {
                        transfrags[pt].abundance = 0.0;
                    }
                }
                transfrags[u_idx].abundance += maxabund;
                available -= maxabund;
                progressed = true;
                if available <= epsilon {
                    break;
                }
            }
        }
        if !progressed {
            break;
        }
    }
}

/// Remove transfrags with abundance < threshold that don't touch source/sink; cap at max_count.
pub fn eliminate_transfrags_under_thr(
    transfrags: Vec<GraphTransfrag>,
    graph: &Graph,
    threshold: f64,
    max_count: usize,
) -> Vec<GraphTransfrag> {
    let source_id = graph.source_id;
    let sink_id = graph.sink_id;
    let mut current = transfrags;
    let mut thr = threshold;
    loop {
        let mut to_remove = Vec::new();
        for (i, tf) in current.iter().enumerate() {
            if tf.abundance >= thr {
                continue;
            }
            if tf.guide {
                // guide transfrags are preserved in threshold elimination.
                continue;
            }
            if tf.node_ids[0] == source_id || tf.node_ids.last().copied() == Some(sink_id) {
                continue;
            }
            to_remove.push(i);
        }
        if to_remove.is_empty() {
            if current.len() <= max_count {
                break;
            }
            thr += 1.0;
            continue;
        }
        let remove_set: HashSet<usize> = to_remove.into_iter().collect();
        current = current
            .into_iter()
            .enumerate()
            .filter(|(i, _)| !remove_set.contains(i))
            .map(|(_, tf)| tf)
            .collect();
        if current.len() <= max_count {
            break;
        }
        thr += 1.0;
    }
    current
}

/// Process transfrags: merge compatible (ret=1,2,3), zero absorbed, update representative abundance.
/// Returns updated transfrags (same vec, abundances modified). No guide mode.
pub fn process_transfrags(
    mut transfrags: Vec<GraphTransfrag>,
    graph: &mut Graph,
    min_abundance: f64,
    singlethr: f64,
    guided_mode: bool,
    long_mode: bool,
    mixed_mode: bool,
    verbose: bool,
    eonly: bool,
    keeptrf_export_path: Option<&str>,
) -> Vec<GraphTransfrag> {
    if transfrags.is_empty() {
        return transfrags;
    }
    let trace_locus = parse_trace_locus_env();
    let trace_intron = parse_trace_intron_env();
    let trace_srcsink_nodes = parse_trace_srcsink_nodes_env();
    if std::env::var_os("RUSTLE_PROMOTE_TERMINAL_TAIL_VARIANTS").is_some() {
        promote_terminal_tail_variants(&mut transfrags, graph, trace_locus);
    }
    let tf_weight = |tf: &GraphTransfrag| -> f64 {
        if mixed_mode {
            tf.abundance.max(tf.srabund)
        } else {
            tf.abundance
        }
    };
    let source_id = graph.source_id;
    let sink_id = graph.sink_id;
    if trace_locus.is_some() {
        for (idx, tf) in transfrags.iter().enumerate() {
            if !tf_overlaps_trace_graph(tf, graph, trace_locus) {
                continue;
            }
            let (span_s, span_e) = tf_span_graph(tf, graph);
            let first_real = tf
                .node_ids
                .iter()
                .copied()
                .find(|&n| n != source_id && n != sink_id)
                .unwrap_or(source_id);
            let last_real = tf
                .node_ids
                .iter()
                .rev()
                .copied()
                .find(|&n| n != source_id && n != sink_id)
                .unwrap_or(sink_id);
            let first_hard = graph
                .nodes
                .get(first_real)
                .map(|n| n.hardstart)
                .unwrap_or(false);
            let last_hard = graph
                .nodes
                .get(last_real)
                .map(|n| n.hardend)
                .unwrap_or(false);
            eprintln!(
                "[TRACE_RAWTRF] idx={} span={}-{} nodes={} abund={:.3} srabund={:.3} guide={} longstart={} longend={} hard={}/{} src_touch={} sink_touch={} real={}",
                idx,
                span_s,
                span_e,
                tf.node_ids.len(),
                tf.abundance,
                tf.srabund,
                tf.guide,
                tf.longstart,
                tf.longend,
                first_hard,
                last_hard,
                tf.node_ids.first().copied() == Some(source_id),
                tf.node_ids.last().copied() == Some(sink_id),
                tf.real
            );
            eprintln!("[TRACE_RAWTRF_NODES] idx={} {:?}", idx, tf.node_ids);
        }
    }
    if let Some(target) = trace_intron {
        for (idx, tf) in transfrags.iter().enumerate() {
            if let Some((a, b)) = tf_find_intron(tf, graph, target) {
                let (span_s, span_e) = tf_span_graph(tf, graph);
                eprintln!(
                    "[TRACE_INTRON_TF] stage=raw idx={} intron={}-{} edge={}=>{} span={}-{} nodes={} abund={:.3} longstart={} longend={}",
                    idx,
                    target.0,
                    target.1,
                    a,
                    b,
                    span_s,
                    span_e,
                    tf.node_ids.len(),
                    tf_weight(tf),
                    tf.longstart,
                    tf.longend
                );
            }
        }
    }
    transfrags = eliminate_transfrags_under_thr(transfrags, graph, min_abundance, MAX_TRF_NUMBER);
    if transfrags.is_empty() {
        return transfrags;
    }
    rebuild_node_trf_ids(graph, &transfrags);

    if guided_mode {
        let source_id = graph.source_id;
        let sink_id = graph.sink_id;
        let psize = graph.pattern_size();
        let mut to_add: Vec<GraphTransfrag> = Vec::new();
        for tf in &transfrags {
            if !tf.guide || tf.node_ids.is_empty() {
                continue;
            }
            let mut core: Vec<usize> = tf
                .node_ids
                .iter()
                .copied()
                .filter(|&n| n != source_id && n != sink_id)
                .collect();
            if core.is_empty() {
                continue;
            }
            let mut seen: HashSet<usize> = Default::default();
            core.retain(|n| seen.insert(*n));
            let first = core[0];
            let last = *core.last().unwrap_or(&first);
            if let Some(n) = graph.nodes.get_mut(first) {
                n.hardstart = true;
            }
            if let Some(n) = graph.nodes.get_mut(last) {
                n.hardend = true;
            }
            graph.add_edge(source_id, first);
            graph.add_edge(last, sink_id);
            let mut full = Vec::with_capacity(core.len() + 2);
            full.push(source_id);
            full.extend(core.iter().copied());
            full.push(sink_id);
            if transfrags.iter().any(|t| t.node_ids == full)
                || to_add.iter().any(|t| t.node_ids == full)
            {
                continue;
            }
            let mut gt = GraphTransfrag::new(full.clone(), psize);
            graph.set_pattern_edges_for_path(&mut gt.pattern, &full);
            gt.guide = true;
            gt.guide_tid = tf.guide_tid.clone();
            gt.longread = true;
            gt.longstart = graph.nodes[first].start;
            gt.longend = graph.nodes[last].end;
            to_add.push(gt);
        }
        if !to_add.is_empty() {
            transfrags.extend(to_add);
            rebuild_node_trf_ids(graph, &transfrags);
        }
    }

    graph.compute_reachability();

    let node_count = graph.nodes.len();
    let mut hassource: Vec<Option<usize>> = vec![None; node_count];
    let mut hassink: Vec<Option<usize>> = vec![None; node_count];
    let mut active_indices = collect_active_transfrag_indices(&transfrags, source_id, sink_id);

    // Incomplete transfrag handling (assign_incomplete_trf_to_nodes + trf_real).
    let mut incompletetrf: Vec<usize> = Vec::new();
    for &tf_idx in &active_indices {
        if tf_idx >= transfrags.len() || transfrags[tf_idx].node_ids.is_empty() {
            continue;
        }
        let mut incomplete = false;
        let nodes = transfrags[tf_idx].node_ids.clone();
        for w in nodes.windows(2) {
            let a = w[0];
            let b = w[1];
            if a == source_id || b == sink_id {
                continue;
            }
            let contiguous = graph
                .nodes
                .get(a)
                .zip(graph.nodes.get(b))
                .map(|(na, nb)| na.end == nb.start)
                .unwrap_or(false);
            if contiguous {
                continue;
            }
            let has_edge_in_pattern = graph
                .edge_bit_index(a, b)
                .map(|eid| transfrags[tf_idx].pattern.get_bit(eid))
                .unwrap_or(false);
            if !has_edge_in_pattern {
                let assigned = assign_incomplete_trf_to_nodes(tf_idx, a.min(b), a.max(b), graph);
                incomplete = incomplete || assigned;
            }
        }
        if incomplete {
            incompletetrf.push(tf_idx);
        } else if let Some(tf) = transfrags.get_mut(tf_idx) {
            tf.real = true;
        }
    }
    for &tf_idx in &incompletetrf {
        let real = trf_real(tf_idx, graph, &mut transfrags);
        if let Some(tf) = transfrags.get_mut(tf_idx) {
            tf.real = real;
        }
    }

    if long_mode && transfrags.len() > 1 {
        let mut order: Vec<usize> = (0..transfrags.len()).collect();
        order.sort_unstable_by(|&a, &b| {
            let ta = &transfrags[a];
            let tb = &transfrags[b];
            tb.guide
                .cmp(&ta.guide)
                .then_with(|| {
                    tf_weight(tb)
                        .partial_cmp(&tf_weight(ta))
                        .unwrap_or(std::cmp::Ordering::Equal)
                })
                .then_with(|| tb.node_ids.len().cmp(&ta.node_ids.len()))
                .then_with(|| tb.pattern.len_bits().cmp(&ta.pattern.len_bits()))
        });
        let already_sorted = order
            .iter()
            .enumerate()
            .all(|(new_idx, &old_idx)| new_idx == old_idx);
        if !already_sorted {
            let mut old = std::mem::take(&mut transfrags)
                .into_iter()
                .map(Some)
                .collect::<Vec<_>>();
            transfrags = order
                .into_iter()
                .map(|old_idx| old[old_idx].take().unwrap())
                .collect();
            rebuild_node_trf_ids(graph, &transfrags);
        }
        // iterates the longtrCmp-sorted transfrag array directly. Rebuild the
        // active non-source/non-sink index list after sorting so later keeptrf and
        // source/sink bookkeeping use the sorted identities, not stale pre-sort ids.
        active_indices = collect_active_transfrag_indices(&transfrags, source_id, sink_id);
    }

    // sorting:
    // - longreads: longtrCmp (guide first, then abundance, then structure)
    // - non-longreads: trCmp-like (larger node count first, then abundance)
    if long_mode {
        active_indices.sort_unstable_by(|&a, &b| {
            let ta = &transfrags[a];
            let tb = &transfrags[b];
            tb.guide
                .cmp(&ta.guide)
                .then_with(|| {
                    tf_weight(tb)
                        .partial_cmp(&tf_weight(ta))
                        .unwrap_or(std::cmp::Ordering::Equal)
                })
                .then_with(|| tb.node_ids.len().cmp(&ta.node_ids.len()))
                .then_with(|| tb.pattern.len_bits().cmp(&ta.pattern.len_bits()))
        });
    } else {
        // trCmp: reach first (last node - first node), then node count,
        // then pattern completeness, then abundance. 
        active_indices.sort_unstable_by(|&a, &b| {
            let ta = &transfrags[a];
            let tb = &transfrags[b];
            let reach = |t: &GraphTransfrag| -> usize {
                t.node_ids
                    .last()
                    .copied()
                    .unwrap_or(0)
                    .saturating_sub(*t.node_ids.first().unwrap_or(&0))
            };
            reach(tb)
                .cmp(&reach(ta))
                .then_with(|| tb.node_ids.len().cmp(&ta.node_ids.len()))
                .then_with(|| tb.pattern.len_bits().cmp(&ta.pattern.len_bits()))
                .then_with(|| {
                    tf_weight(tb)
                        .partial_cmp(&tf_weight(ta))
                        .unwrap_or(std::cmp::Ordering::Equal)
                })
        });
    }

    // in longread mode, hassource/hassink are assigned while
    // iterating transfrags in longtrCmp order (including source/sink transfrags).
    if long_mode {
        let mut all_order: Vec<usize> = (0..transfrags.len()).collect();
        all_order.sort_unstable_by(|&a, &b| {
            let ta = &transfrags[a];
            let tb = &transfrags[b];
            tb.guide
                .cmp(&ta.guide)
                .then_with(|| {
                    tf_weight(tb)
                        .partial_cmp(&tf_weight(ta))
                        .unwrap_or(std::cmp::Ordering::Equal)
                })
                .then_with(|| tb.node_ids.len().cmp(&ta.node_ids.len()))
                .then_with(|| tb.pattern.len_bits().cmp(&ta.pattern.len_bits()))
        });
        hassource.fill(None);
        hassink.fill(None);
        for i in all_order {
            if i >= transfrags.len() || transfrags[i].node_ids.len() < 2 {
                continue;
            }
            if transfrags[i].node_ids[0] == source_id {
                hassource[transfrags[i].node_ids[1]] = Some(i);
            } else if transfrags[i].node_ids.last().copied() == Some(sink_id) {
                hassink[transfrags[i].node_ids[0]] = Some(i);
            }
        }
    }
    if verbose {
        let mut top_abund: Vec<f64> = active_indices
            .iter()
            .take(5)
            .map(|&i| tf_weight(&transfrags[i]))
            .collect();
        top_abund.shrink_to_fit();
        eprintln!(
            "    process_transfrags: abundance_sort active (top={:?})",
            top_abund
        );
    }

    // keeptrf: (rep_idx, group_indices, group_cov)
    let mut keeptrf: Vec<(usize, Vec<usize>, f64)> = Vec::new();

    for &tf_idx in &active_indices {
        let tf = &transfrags[tf_idx];
        let first_nid = tf.node_ids[0];
        let last_nid = tf.node_ids[tf.node_ids.len() - 1];
        let first_node = &graph.nodes[first_nid];
        let last_node = &graph.nodes[last_nid];
        let tf_junction_chain = tf_junction_chain_coords(tf, graph, source_id, sink_id);

        let mut included = false;
        let mut included_via_group = false;
        let mut pending_guide_replace: Option<(usize, usize)> = None; // (kept_rep_idx, tf_idx)
        let mut mark_weak = false; // defer weak marking to avoid borrow conflict
        // Track whether tf REPLACED an existing keeptrf rep (ret=1 absorb-kept-into-tf
        // or ret=3 new_rep path). In those cases the "included=true" flag is
        // misleading — tf is not absorbed, it's the new representative — so we
        // must skip the mark_weak fallback that would zero out its trflong_seed
        // eligibility. This fixes the STRG.1.3/1.5 class of miss where the new
        // rep was erroneously dropped before parse_trflong ran.
        let mut tf_is_new_rep = false;
        for (kept_idx, group, kept_cov) in keeptrf.iter_mut() {
            let kept_tf = &transfrags[*kept_idx];
            let trace_pair = tf_overlaps_trace_graph(tf, graph, trace_locus)
                || tf_overlaps_trace_graph(kept_tf, graph, trace_locus);

            // guide containment branch (~6032-6081):
            // - if non-guide tf is fully contained in kept guide => included
            // - if kept guide is fully contained in non-guide tf and side distances are close => included;
            //   if kept guide has zero abundance, replace its abundance with tf abundance.
            if kept_tf.guide {
                if !tf.guide && kept_tf.pattern.contains_pattern(&tf.pattern) {
                    included = true;
                    break;
                } else if !tf.guide && tf.pattern.contains_pattern(&kept_tf.pattern) {
                    let mut contain = true;
                    let mut leftdist: i64 = 0;
                    let mut rightdist: i64 = 0;

                    let mut i = 0usize;
                    while i + 1 < tf.node_ids.len() && !kept_tf.pattern.get_bit(tf.node_ids[i]) {
                        let curr = &graph.nodes[tf.node_ids[i]];
                        let next = &graph.nodes[tf.node_ids[i + 1]];
                        if next.start > curr.end.saturating_add(1) {
                            contain = false;
                            break;
                        }
                        leftdist += next.length() as i64;
                        if i == 0 && tf.longstart > curr.start {
                            leftdist -= (tf.longstart - curr.start) as i64;
                        }
                        i += 1;
                    }
                    if contain {
                        let mut i = tf.node_ids.len() - 1;
                        let mut first = true;
                        while i > 0 && !kept_tf.pattern.get_bit(tf.node_ids[i]) {
                            let curr = &graph.nodes[tf.node_ids[i]];
                            let prev = &graph.nodes[tf.node_ids[i - 1]];
                            if curr.start > prev.end.saturating_add(1) {
                                contain = false;
                                break;
                            }
                            if first && tf.longend != 0 && tf.longend < curr.end {
                                rightdist -= (curr.end - tf.longend) as i64;
                                first = false;
                            }
                            i -= 1;
                        }
                    }
                    if contain {
                        if kept_tf.abundance <= 0.0 && leftdist < MAX_DIST && rightdist < MAX_DIST {
                            pending_guide_replace = Some((*kept_idx, tf_idx));
                            *kept_cov = tf_weight(tf);
                        }
                        included = true;
                        break;
                    }
                }
            }

            if guided_mode && !tf.guide && tf_weight(tf) < singlethr {
                included = true;
            }
            if included {
                continue;
            }

            let (ret, lens) = compatible_long(tf, kept_tf, graph);
            if trace_pair {
                let site = if ret == 0 {
                    LAST_CL_INCOMPAT_SITE.with(|c| c.get())
                } else {
                    ""
                };
                let (tf_s, tf_e) = tf_span_graph(tf, graph);
                let (kt_s, kt_e) = tf_span_graph(kept_tf, graph);
                eprintln!(
                    "[TRACE_SEEDMERGE] cand={} span={}-{} nodes={} abund={:.3} vs kept={} span={}-{} nodes={} abund={:.3} ret={} site={} lens=[{},{},{},{}] guide_pair={}/{} hard={}/{} kept_hard={}/{}",
                    tf_idx,
                    tf_s,
                    tf_e,
                    tf.node_ids.len(),
                    tf_weight(tf),
                    *kept_idx,
                    kt_s,
                    kt_e,
                    kept_tf.node_ids.len(),
                    tf_weight(kept_tf),
                    ret,
                    site,
                    lens[0],
                    lens[1],
                    lens[2],
                    lens[3],
                    tf.guide,
                    kept_tf.guide,
                    first_node.hardstart || tf.longstart != 0,
                    last_node.hardend || tf.longend != 0,
                    graph.nodes[kept_tf.node_ids[0]].hardstart || kept_tf.longstart != 0,
                    graph.nodes[*kept_tf.node_ids.last().unwrap()].hardend || kept_tf.longend != 0
                );
            }
            if ret == 0 {
                continue;
            }
            if ret == 3 {
                // Opposing hard boundaries → keep both ( reference 10759-10765)
                let tf_first_node = &graph.nodes[first_nid];
                let tf_last_node = &graph.nodes[last_nid];
                let kt_first_node = &graph.nodes[kept_tf.node_ids[0]];
                let kt_last_node = &graph.nodes[kept_tf.node_ids[kept_tf.node_ids.len() - 1]];
                if first_nid != kept_tf.node_ids[0]
                    && last_nid != kept_tf.node_ids[kept_tf.node_ids.len() - 1]
                {
                    let opposing = (tf_first_node.hardstart
                        && !tf_last_node.hardend
                        && !kt_first_node.hardstart
                        && kt_last_node.hardend)
                        || (!tf_first_node.hardstart
                            && tf_last_node.hardend
                            && kt_first_node.hardstart
                            && !kt_last_node.hardend);
                    if opposing {
                        if trace_pair {
                            eprintln!(
                                "[TRACE_SEEDMERGE]   decision=keep_both reason=opposing_hard_boundaries"
                            );
                        }
                        continue; // keep both — different transcript boundaries
                    }
                }
                // Merge only if both terminal distances < edgedist; signed lens[0], lens[2]
                let left_dist = lens[0];
                let right_dist = lens[2];
                if (!tf.guide || !kept_tf.guide) && left_dist < EDGEDIST && right_dist < EDGEDIST {
                    // 6144-6145:
                    // replace kept with current only if current is guide OR
                    // kept is non-guide and current has hardstart+hardend.
                    //
                    // Added guard: don't replace kept with a structurally poorer
                    // (fewer-node) cand. Otherwise a truncated hardstart+hardend
                    // cand can wipe out the richer kept rep — e.g. STRG.225.5
                    // retained-intron variant (21 nodes, abund=21) being replaced
                    // by a truncated 17-node cand ending before the retained
                    // intron region, dropping the retained-intron seed entirely.
                    let new_rep = tf.guide
                        || (!kept_tf.guide
                            && tf_first_node.hardstart
                            && tf_last_node.hardend
                            && tf.node_ids.len() >= kept_tf.node_ids.len());
                    if new_rep {
                        *kept_idx = tf_idx;
                        tf_is_new_rep = true;
                    }
                    group.push(tf_idx);
                    *kept_cov += tf_weight(tf);
                    included = true;
                    included_via_group = true;
                    if trace_pair {
                        eprintln!(
                            "[TRACE_SEEDMERGE]   decision=merge_ret3 new_rep={} group_size={} group_cov={:.3}",
                            new_rep,
                            group.len(),
                            *kept_cov
                        );
                    }
                    break;
                }
                if trace_pair {
                    eprintln!(
                        "[TRACE_SEEDMERGE]   decision=keep_both reason=ret3_edgedist left={} right={}",
                        left_dist, right_dist
                    );
                }
            }
            if ret == 1 {
                // tf extends past kept: absorb kept into tf (tf becomes new rep) only if ( reference 10828-10835)
                let kt_first = &graph.nodes[kept_tf.node_ids[0]];
                let kt_last = &graph.nodes[kept_tf.node_ids[kept_tf.node_ids.len() - 1]];
                let left_dist = lens[1];
                let right_dist = lens[3];
                if !kept_tf.guide
                    && tf.longstart != 0
                    && tf.longend != 0
                    && (!kt_first.hardstart || first_nid == kept_tf.node_ids[0])
                    && (!kt_last.hardend
                        || last_nid == kept_tf.node_ids[kept_tf.node_ids.len() - 1])
                    && tf_weight(tf) > DROP * tf_weight(kept_tf)
                    && left_dist < SSDIST
                    && right_dist < SSDIST
                {
                    *kept_idx = tf_idx;
                    tf_is_new_rep = true;
                    group.push(tf_idx);
                    *kept_cov += tf_weight(tf);
                    included = true;
                    included_via_group = true;
                    if trace_pair {
                        eprintln!(
                            "[TRACE_SEEDMERGE]   decision=merge_ret1 replace_kept=true group_size={} group_cov={:.3}",
                            group.len(),
                            *kept_cov
                        );
                    }
                    break;
                }
                if trace_pair {
                    eprintln!("[TRACE_SEEDMERGE]   decision=keep_both reason=ret1_gate");
                }
            }
            if ret == 2 {
                // kept extends past tf: absorb tf into kept only if ( reference 10869-10895)
                let left_dist = lens[1];
                let right_dist = lens[3];
                if !tf.guide && left_dist < SSDIST && right_dist < SSDIST {
                    let tf_first = &graph.nodes[first_nid];
                    let tf_last = &graph.nodes[last_nid];
                    if !(tf_first.hardstart || tf.longstart != 0)
                        && !(tf_last.hardend || tf.longend != 0)
                    {
                        let start_ok = !tf_first.hardstart || first_nid == kept_tf.node_ids[0];
                        let end_ok = !tf_last.hardend
                            || last_nid == kept_tf.node_ids[kept_tf.node_ids.len() - 1];
                        if start_ok && end_ok {
                            group.push(tf_idx);
                            *kept_cov += tf_weight(tf);
                            included = true;
                            included_via_group = true;
                            if trace_pair {
                                eprintln!(
                                    "[TRACE_SEEDMERGE]   decision=merge_ret2 group_size={} group_cov={:.3}",
                                    group.len(),
                                    *kept_cov
                                );
                            }
                            break;
                        }
                    }
                }
                if trace_pair {
                    eprintln!("[TRACE_SEEDMERGE]   decision=keep_both reason=ret2_gate");
                }
            }
        }
        // Long-read novel-splice rescue: if containment suppressed this transfrag but it carries
        // at least one splice junction not represented by current keeptrf representatives,
        // keep it eligible as its own seed even when it is only partially bounded.
        //
        // This matters for needy refs like STRG.29 where a left-anchored long-read fragment
        // carries the distinguishing early junctions, but keeptrf would otherwise absorb it into
        // a downstream representative that lacks those junctions.
        //
        // StringTie's rlink.cpp process_transfrags has no equivalent — contained tfs are
        // strictly marked weak=1. Disabled in RUSTLE_STRINGTIE_EXACT mode for parity.
        let novel_splice_rescue = !crate::stringtie_parity::stringtie_exact()
            && included
            && !included_via_group
            && long_mode
            && !tf.guide
            && (tf.longstart != 0 || tf.longend != 0)
            && tf_has_novel_junction_vs_keeptrf(
                &tf_junction_chain,
                &keeptrf,
                &transfrags,
                graph,
                source_id,
                sink_id,
            );
        if novel_splice_rescue {
            included = false;
        }
        if !included {
            // guides are added to keeptrf even without proper boundaries.
            let has_source_strict = first_node.hardstart || tf.longstart != 0;
            let has_sink = last_node.hardend || tf.longend != 0;
            // Short-tf stricter gate (opt-in RUSTLE_SHORT_TF_STRICT=1): for tfs with
            // ≤2 nodes, require BOTH boundary evidence. Prevents short mid-chain
            // tfs (like STRG.398's t=82) from becoming keeptrf seeds and depleting
            // dominant tfs via aggressive back/fwd extension. Opt-in because the
            // true STRG.398 divergence is deeper (graph alt-junction 10→14 lets
            // tf=42 avoid cov-drop weak marking and win over the preserved tf=8).
            let short_tf_stricter = !tf.guide
                && tf.node_ids.len() <= 2
                && std::env::var_os("RUSTLE_SHORT_TF_STRICT").is_some();
            let strict_pass = if short_tf_stricter {
                tf.longstart != 0 && tf.longend != 0
            } else {
                tf.guide || (has_source_strict && has_sink)
            };
            // Register chain for family-rescue lookup when a transfrag passes the
            // strict boundary gate. Later bundles' single-boundary transfrags that
            // match this chain (a confirmed family member) become eligible for
            // rescue under `RUSTLE_VG_FAMILY_RESCUE`.
            if strict_pass && long_mode && !tf.guide {
                let chain = tf_intron_length_chain(tf, graph);
                if chain.len() >= 3 {
                    chain_register(chain);
                }
            }
            // VG-family rescue: allow single-boundary long-read transfrags to become
            // keeptrf reps when their intron-length chain matches a registered
            // family member's chain. Two match modes:
            //   - default: exact contiguous-subsequence match within ±5 bp/intron
            //     (catches tight tandem paralogs like GOLGA6L7).
            //   - RUSTLE_VG_FAMILY_TOLERANT=1: Needleman-Wunsch alignment with
            //     score ≥ -500 (catches distant sub-family members whose chains
            //     have diverged past exact matching, e.g. GOLGA6C cluster).
            let vg_family_rescue = long_mode
                && !tf.guide
                && (has_source_strict || has_sink)
                && !strict_pass
                && std::env::var_os("RUSTLE_VG_FAMILY_RESCUE").is_some()
                && {
                    let chain = tf_intron_length_chain(tf, graph);
                    if std::env::var_os("RUSTLE_VG_FAMILY_TOLERANT").is_some() {
                        chain_aligns_to_registry(&chain, -500.0, FAMILY_CHAIN_MIN)
                    } else {
                        chain_matches_registry(&chain, 5)
                    }
                };
            if tf.guide || strict_pass || novel_splice_rescue || vg_family_rescue {
                if tf_overlaps_trace_graph(tf, graph, trace_locus) {
                    let (tf_s, tf_e) = tf_span_graph(tf, graph);
                    eprintln!(
                        "[TRACE_SEEDMERGE] keeptrf_add={} span={}-{} nodes={} abund={:.3} complete={} novel_splice_rescue={}",
                        tf_idx,
                        tf_s,
                        tf_e,
                        tf.node_ids.len(),
                        tf_weight(tf),
                        has_source_strict && has_sink,
                        novel_splice_rescue
                    );
                }
                if let Some(target) = trace_intron {
                    if let Some((a, b)) = tf_find_intron(tf, graph, target) {
                        let (tf_s, tf_e) = tf_span_graph(tf, graph);
                        eprintln!(
                            "[TRACE_INTRON_TF] stage=keeptrf_add idx={} intron={}-{} edge={}=>{} span={}-{} nodes={} abund={:.3} novel_splice_rescue={}",
                            tf_idx,
                            target.0,
                            target.1,
                            a,
                            b,
                            tf_s,
                            tf_e,
                            tf.node_ids.len(),
                            tf_weight(tf),
                            novel_splice_rescue
                        );
                    }
                }
                keeptrf.push((tf_idx, vec![tf_idx], tf_weight(tf)));
            } else {
                // transcript not included and lacks boundaries -> mark as weak
                transfrags[tf_idx].weak = 1;
            }
        } else {
            // ALL included transcripts are marked as weak.
            // Defer weak marking to avoid borrow conflict.
            mark_weak = true;
            if !included_via_group {
                // Included but not attached to any keeptrf group (e.g. guided low-abundance suppression).
                let trace_this = tf_overlaps_trace_graph(tf, graph, trace_locus);
                let trace_span = if trace_this {
                    Some(tf_span_graph(tf, graph))
                } else {
                    None
                };
                let tf_nodes_len = tf.node_ids.len();
                let tf_ab = tf_weight(tf);
                if let Some((tf_s, tf_e)) = trace_span {
                    eprintln!(
                        "[TRACE_SEEDMERGE] mark_weak={} span={}-{} nodes={} abund={:.3}",
                        tf_idx, tf_s, tf_e, tf_nodes_len, tf_ab
                    );
                }
            }
        }
        if let Some((rep_idx, src_idx)) = pending_guide_replace {
            let new_ab = transfrags[src_idx].abundance;
            transfrags[rep_idx].abundance = new_ab;
            transfrags[src_idx].abundance = 0.0;
        }
        // apply deferred weak marking to avoid borrow conflict. Skip when tf
        // replaced an existing keeptrf entry as the new rep — mark_weak is a
        // fallback for absorbed transfrags, not for promoted ones.
        if mark_weak && !tf_is_new_rep {
            transfrags[tf_idx].weak = 1;
        }
    }

    // compute keepsource/keepsink from the final keeptrf groups.
    //
    // Important nuance (matches the original algorithm behavior in practice): keep source/sink support only
    // for the representative start/end nodes that survived keeptrf selection. Computing these
    // flags over all active transfrags can incorrectly keep helpers at internal nodes (nodes
    // with non-source parents / non-sink children), which can siphon flux and cause misses.
    let mut keepsource = vec![false; node_count];
    let mut keepsink = vec![false; node_count];
    for (rep_idx, group, _group_cov) in &keeptrf {
        if *rep_idx >= transfrags.len() || transfrags[*rep_idx].node_ids.is_empty() {
            continue;
        }
        let rep_start = transfrags[*rep_idx].node_ids[0];
        let rep_end = *transfrags[*rep_idx].node_ids.last().unwrap();
        for &gid in group {
            if gid >= transfrags.len() || transfrags[gid].node_ids.is_empty() {
                continue;
            }
            let tf = &transfrags[gid];
            // eonly skips non-guide transfrags in keepsource/keepsink collection.
            if eonly && !tf.guide {
                continue;
            }
            // Keep flags are used to gate source/sink helpers whose abundance comes from
            // addsource/addsink, which skip guide transfrags; mirror that here to avoid
            // "kept but zero-abundance" helpers affecting completeness checks.
            if tf.guide {
                continue;
            }
            if tf.node_ids[0] == rep_start {
                if tf.longstart != 0
                    || (tf.longread
                        && graph.nodes.get(rep_start).map(|n| n.hardstart).unwrap_or(false))
                {
                    keepsource[rep_start] = true;
                }
            }
            if tf.node_ids.last().copied() == Some(rep_end) {
                if tf.longend != 0
                    || graph.nodes.get(rep_end).map(|n| n.hardend).unwrap_or(false)
                {
                    keepsink[rep_end] = true;
                }
            }
        }
    }

    // build trflong seed order from keeptrf before source/sink rewiring.
    // pass 1 (reverse keeptrf): incomplete non-guide
    // pass 2 (reverse keeptrf): complete or guide
    //
    // StringTie also runs pass 2 (rlink.cpp:6470-6477); Rustle matches. Pass 2
    // is DEFAULT-ON because disabling it regresses 1980→1136 predictions. Opt-out
    // via RUSTLE_TRFLONG_PASS2_OFF=1 for StringTie-parity experiments (does not
    // fix STRG.398 over-emission — divergence is actually in hassink[]
    // population, not pass-2 logic).
    let enable_pass2 = std::env::var_os("RUSTLE_TRFLONG_PASS2_OFF").is_none();
    let trace_trflong = std::env::var_os("RUSTLE_TRFLONG_TRACE").is_some();
    let rep_order_now: Vec<usize> = keeptrf.iter().map(|(r, _, _)| *r).collect();
    let mut trflong_insert: Vec<usize> = Vec::new();
    let mut n_complete = 0usize;
    let mut n_incomplete = 0usize;
    let mut n_guided = 0usize;
    for rid in rep_order_now.iter().rev().copied() {
        if transfrags[rid].node_ids.is_empty() {
            continue;
        }
        let n1 = transfrags[rid].node_ids[0];
        let n2 = *transfrags[rid].node_ids.last().unwrap();
        let hardstart = graph.nodes.get(n1).map(|n| n.hardstart).unwrap_or(false);
        let hardend = graph.nodes.get(n2).map(|n| n.hardend).unwrap_or(false);
        let has_source_keep = hassource[n1].is_some() && keepsource[n1];
        let has_sink_keep = hassink[n2].is_some() && keepsink[n2];
        let complete = transfrags[rid].guide
            || ((hardstart || has_source_keep)
                && (hardend || has_sink_keep || has_sink_completion(n2, graph, &keepsink, &hassink)));
        let incomplete = !transfrags[rid].guide
            && ((!hardstart && !has_source_keep) || (!hardend && !has_sink_keep));
        if transfrags[rid].guide {
            n_guided += 1;
        }
        if complete {
            n_complete += 1;
        }
        if trace_trflong && tf_overlaps_trace_graph(&transfrags[rid], graph, trace_locus) {
            let (sp_s, sp_e) = tf_span_graph(&transfrags[rid], graph);
            eprintln!(
                "[RUSTLE_TRFLONG] pass=1 tf={} first={} last={} coord={}-{} abund={:.3} longstart={} longend={} hardstart={} hardend={} hassource={:?} keepsource={} hassink={:?} keepsink={} decision={}",
                rid, n1, n2, sp_s, sp_e, tf_weight(&transfrags[rid]),
                transfrags[rid].longstart, transfrags[rid].longend,
                hardstart, hardend,
                hassource[n1], keepsource[n1],
                hassink[n2], keepsink[n2],
                if incomplete { "ADD" } else { "SKIP" }
            );
        }
        if incomplete {
            n_incomplete += 1;
            trflong_insert.push(rid);
        }
    }
    if enable_pass2 {
        for rid in rep_order_now.iter().rev().copied() {
            if transfrags[rid].node_ids.is_empty() {
                continue;
            }
            let n1 = transfrags[rid].node_ids[0];
            let n2 = *transfrags[rid].node_ids.last().unwrap();
            let hardstart = graph.nodes.get(n1).map(|n| n.hardstart).unwrap_or(false);
            let hardend = graph.nodes.get(n2).map(|n| n.hardend).unwrap_or(false);
            let has_source_keep = hassource[n1].is_some() && keepsource[n1];
            let has_sink_keep = hassink[n2].is_some() && keepsink[n2];
            let complete = transfrags[rid].guide
                || ((hardstart || has_source_keep)
                    && (hardend || has_sink_keep || has_sink_completion(n2, graph, &keepsink, &hassink)));
            if trace_trflong && tf_overlaps_trace_graph(&transfrags[rid], graph, trace_locus) {
                let (sp_s, sp_e) = tf_span_graph(&transfrags[rid], graph);
                eprintln!(
                    "[RUSTLE_TRFLONG] pass=2 tf={} first={} last={} coord={}-{} abund={:.3} decision={}",
                    rid, n1, n2, sp_s, sp_e, tf_weight(&transfrags[rid]),
                    if transfrags[rid].guide || complete { "ADD" } else { "SKIP" }
                );
            }
            if transfrags[rid].guide || complete {
                trflong_insert.push(rid);
            }
        }
    }

    // Diagnostic: when RUSTLE_NOSEED_DIAG is set, print transfrag details for bundles with 0 seeds.
    if std::env::var_os("RUSTLE_NOSEED_DIAG").is_some()
        && trflong_insert.is_empty()
        && transfrags
            .iter()
            .any(|tf| tf.longread && tf.abundance > 0.0)
    {
        let n_long = transfrags
            .iter()
            .filter(|tf| tf.longread && tf.abundance > 0.0)
            .count();
        eprintln!(
            "[NOSEED] active={} longread={} keeptrf={}",
            active_indices.len(),
            n_long,
            keeptrf.len()
        );
        for &tf_idx in &active_indices {
            let tf = &transfrags[tf_idx];
            if !tf.longread {
                continue;
            }
            let first_nid = tf.node_ids[0];
            let last_nid = *tf.node_ids.last().unwrap_or(&0);
            let hardstart = graph
                .nodes
                .get(first_nid)
                .map(|n| n.hardstart)
                .unwrap_or(false);
            let hardend = graph
                .nodes
                .get(last_nid)
                .map(|n| n.hardend)
                .unwrap_or(false);
            let (ns, ne) = graph
                .nodes
                .get(first_nid)
                .map(|n| (n.start, n.end))
                .unwrap_or((0, 0));
            let (ls, le) = graph
                .nodes
                .get(last_nid)
                .map(|n| (n.start, n.end))
                .unwrap_or((0, 0));
            let passes_gate = (tf.longstart != 0 || hardstart) && (tf.longend != 0 || hardend);
            eprintln!(
                "[NOSEED_TF] idx={} abund={:.2} longstart={} longend={} hardstart={} hardend={} first={}-{} last={}-{} weak={} guide={} included={} passes_keeptrf_gate={}",
                tf_idx, tf.abundance, tf.longstart, tf.longend, hardstart, hardend,
                ns, ne, ls, le, tf.weak, tf.guide, tf.weak != 0, passes_gate
            );
        }
    }

    // source/sink support rewiring after keeptrf selection.
    // Aggregate grouped abundance by representative start/end, then
    // create/update source->start and end->sink transfrags for hard boundaries.
    let mut addsource = vec![0.0; node_count];
    let mut addsink = vec![0.0; node_count];
    let mut source_contributors: Option<Vec<Vec<(usize, f64)>>> = if trace_srcsink_nodes.is_empty()
    {
        None
    } else {
        Some(vec![Vec::new(); node_count])
    };
    for (rep_idx, group, _group_cov) in &keeptrf {
        if transfrags[*rep_idx].node_ids.is_empty() {
            continue;
        }
        let n1 = transfrags[*rep_idx].node_ids[0];
        let n2 = *transfrags[*rep_idx].node_ids.last().unwrap();
        for &gid in group {
            if gid >= transfrags.len() || transfrags[gid].node_ids.is_empty() {
                continue;
            }
            // skip guide transfrags in addsource/addsink accumulation
            if transfrags[gid].guide {
                continue;
            }
            if transfrags[gid].node_ids[0] == n1 {
                let weight = tf_weight(&transfrags[gid]);
                addsource[n1] += weight;
                if let Some(src_groups) = source_contributors.as_mut() {
                    if trace_srcsink_nodes.binary_search(&n1).is_ok() {
                        src_groups[n1].push((gid, weight));
                    }
                }
            }
            if transfrags[gid].node_ids.last().copied() == Some(n2) {
                addsink[n2] += tf_weight(&transfrags[gid]);
            }
        }
    }

    if let Some(src_groups) = source_contributors.as_ref() {
        let mut contributes_to_source = vec![false; transfrags.len()];
        for &nid in &trace_srcsink_nodes {
            if nid >= src_groups.len() {
                continue;
            }
            for &(gid, _) in &src_groups[nid] {
                if gid < contributes_to_source.len() {
                    contributes_to_source[gid] = true;
                }
            }
        }

        for &nid in &trace_srcsink_nodes {
            if nid >= graph.nodes.len() {
                continue;
            }
            let node = &graph.nodes[nid];
            let mut rep_lines = Vec::new();
            for (rep_idx, group, group_cov) in &keeptrf {
                let Some(tf) = transfrags.get(*rep_idx) else {
                    continue;
                };
                if tf.node_ids.first().copied() != Some(nid) {
                    continue;
                }
                let mut members = Vec::new();
                let mut src_sum = 0.0;
                for &gid in group {
                    let Some(gtf) = transfrags.get(gid) else {
                        continue;
                    };
                    let first = gtf.node_ids.first().copied();
                    let contributes = first == Some(nid);
                    if contributes {
                        src_sum += tf_weight(gtf);
                    }
                    members.push(format!(
                        "{}:{:.3}:start={:?}:end={:?}:ls={}:le={}:weak={}:src={}",
                        gid,
                        tf_weight(gtf),
                        first,
                        gtf.node_ids.last().copied(),
                        gtf.longstart,
                        gtf.longend,
                        gtf.weak,
                        if contributes { 1 } else { 0 }
                    ));
                }
                rep_lines.push(format!(
                    "rep={} rep_abund={:.3} group_cov={:.3} src_sum={:.3} group=[{}]",
                    rep_idx,
                    tf_weight(tf),
                    *group_cov,
                    src_sum,
                    members.join(" ")
                ));
            }
            eprintln!(
                "[TRACE_SRCSINK_GROUP] node={}({}-{}) addsource={:.3} hassource={:?} keepsource={} reps={}",
                nid,
                node.start,
                node.end,
                addsource[nid],
                hassource[nid],
                keepsource[nid],
                rep_lines.len()
            );
            for line in rep_lines {
                eprintln!("[TRACE_SRCSINK_GROUP]   {}", line);
            }

            let mut missed = Vec::new();
            for (tf_idx, tf) in transfrags.iter().enumerate() {
                if tf.node_ids.first().copied() != Some(nid) {
                    continue;
                }
                if tf_weight(tf) <= 0.0 {
                    continue;
                }
                if contributes_to_source.get(tf_idx).copied().unwrap_or(false) {
                    continue;
                }
                missed.push(format!(
                    "{}:{:.3}:nodes={}:last={:?}:ls={}:le={}:weak={}:guide={}:usepath={}",
                    tf_idx,
                    tf_weight(tf),
                    tf.node_ids.len(),
                    tf.node_ids.last().copied(),
                    tf.longstart,
                    tf.longend,
                    tf.weak,
                    tf.guide,
                    tf.usepath
                ));
            }
            if !missed.is_empty() {
                eprintln!(
                    "[TRACE_SRCSINK_GROUP]   missed_start_tfs node={} -> [{}]",
                    nid,
                    missed.join(" ")
                );
            }
        }
    }

    if let Some((lo, hi)) = trace_locus {
        for nid in 1..sink_id {
            let Some(node) = graph.nodes.get(nid) else {
                continue;
            };
            if node.start > hi || node.end < lo {
                continue;
            }
            let mut start_reps: Vec<usize> = keeptrf
                .iter()
                .filter_map(|(rep_idx, _group, group_cov)| {
                    let tf = transfrags.get(*rep_idx)?;
                    if tf.node_ids.first().copied() == Some(nid) {
                        Some((*rep_idx, *group_cov))
                    } else {
                        None
                    }
                })
                .map(|(rep_idx, _group_cov)| {
                    eprintln!(
                        "[TRACE_SRCSINK_PRE] node={}({}-{}) rep={} addsource={:.3} addsink={:.3} hardstart={} hardend={} hassource={:?} hassink={:?} keepsource={} keepsink={} parents={:?} children={:?}",
                        nid,
                        node.start,
                        node.end,
                        rep_idx,
                        addsource[nid],
                        addsink[nid],
                        node.hardstart,
                        node.hardend,
                        hassource[nid],
                        hassink[nid],
                        keepsource[nid],
                        keepsink[nid],
                        node.parents.ones().collect::<Vec<_>>(),
                        node.children.ones().collect::<Vec<_>>()
                    );
                    rep_idx
                })
                .collect();
            start_reps.sort_unstable();
            if start_reps.is_empty()
                && ((addsource[nid] > 0.0)
                    || (addsink[nid] > 0.0)
                    || hassource[nid].is_some()
                    || hassink[nid].is_some()
                    || node.parents.contains(source_id)
                    || node.children.contains(sink_id))
            {
                eprintln!(
                    "[TRACE_SRCSINK_PRE] node={}({}-{}) rep=- addsource={:.3} addsink={:.3} hardstart={} hardend={} hassource={:?} hassink={:?} keepsource={} keepsink={} parents={:?} children={:?}",
                    nid,
                    node.start,
                    node.end,
                    addsource[nid],
                    addsink[nid],
                    node.hardstart,
                    node.hardend,
                    hassource[nid],
                    hassink[nid],
                    keepsource[nid],
                    keepsink[nid],
                    node.parents.ones().collect::<Vec<_>>(),
                    node.children.ones().collect::<Vec<_>>()
                );
            }
        }
    }

    for nid in 1..sink_id {
        let src_ab = addsource[nid];
        let snk_ab = addsink[nid];

        if let Some(src_tf_idx) = hassource[nid] {
            if !keepsource[nid] {
                transfrags[src_tf_idx].abundance = 0.0;
            } else {
                transfrags[src_tf_idx].abundance = src_ab;
            }
        } else if src_ab > 0.0 && graph.nodes.get(nid).map(|n| n.hardstart).unwrap_or(false) {
            graph.add_edge(source_id, nid);
            let mut tf = GraphTransfrag::new(vec![source_id, nid], graph.pattern_size());
            tf.abundance = src_ab;
            tf.longread = true;
            graph.set_pattern_edges_for_path(&mut tf.pattern, &tf.node_ids);
            let new_idx = transfrags.len();
            transfrags.push(tf);
            graph.nodes[source_id].trf_ids.push(new_idx);
            graph.nodes[nid].trf_ids.push(new_idx);
            hassource[nid] = Some(new_idx);
        }

        if let Some(snk_tf_idx) = hassink[nid] {
            if !keepsink[nid] {
                transfrags[snk_tf_idx].abundance = 0.0;
            } else {
                transfrags[snk_tf_idx].abundance = snk_ab;
            }
        } else if snk_ab > 0.0 && graph.nodes.get(nid).map(|n| n.hardend).unwrap_or(false) {
            graph.add_edge(nid, sink_id);
            let mut tf = GraphTransfrag::new(vec![nid, sink_id], graph.pattern_size());
            tf.abundance = snk_ab;
            tf.longread = true;
            graph.set_pattern_edges_for_path(&mut tf.pattern, &tf.node_ids);
            let new_idx = transfrags.len();
            transfrags.push(tf);
            graph.nodes[nid].trf_ids.push(new_idx);
            graph.nodes[sink_id].trf_ids.push(new_idx);
            hassink[nid] = Some(new_idx);
        }
    }

    if let Some((lo, hi)) = trace_locus {
        let source_trfs: Vec<String> = graph
            .nodes
            .get(source_id)
            .map(|n| n.trf_ids.clone())
            .unwrap_or_default()
            .into_iter()
            .filter_map(|t| {
                let tf = transfrags.get(t)?;
                let (span_s, span_e) = tf_span_graph(tf, graph);
                if span_e < lo || span_s > hi {
                    return None;
                }
                Some(format!(
                    "{}:{}-{} nodes={:?} abund={:.3} long={}",
                    t, span_s, span_e, tf.node_ids, tf.abundance, tf.longread
                ))
            })
            .collect();
        let sink_trfs: Vec<String> = graph
            .nodes
            .get(sink_id)
            .map(|n| n.trf_ids.clone())
            .unwrap_or_default()
            .into_iter()
            .filter_map(|t| {
                let tf = transfrags.get(t)?;
                let (span_s, span_e) = tf_span_graph(tf, graph);
                if span_e < lo || span_s > hi {
                    return None;
                }
                Some(format!(
                    "{}:{}-{} nodes={:?} abund={:.3} long={} longstart={} longend={}",
                    t,
                    span_s,
                    span_e,
                    tf.node_ids,
                    tf.abundance,
                    tf.longread,
                    tf.longstart,
                    tf.longend
                ))
            })
            .collect();
        eprintln!(
            "[TRACE_SRCSINK_POST] source_trfs={}",
            source_trfs.join(" | ")
        );
        eprintln!("[TRACE_SRCSINK_POST] sink_trfs={}", sink_trfs.join(" | "));
    }

    // Edges were potentially added/removed above (source/sink rewiring); refresh reachability
    // before long-read path extension stages consume parentpat/childpat.
    graph.compute_reachability();

    // keeps grouped support in keeptrf.cov and leaves the representative transfrag
    // abundance unchanged for parse_trflong/checktrf. Only mark absorbed members here.
    let rep_set: HashSet<usize> = keeptrf.iter().map(|(r, _, _)| *r).collect();
    for tf in transfrags.iter_mut() {
        tf.trflong_seed = false;
        tf.usepath = -1;
    }
    for (rep_idx, group, _group_cov) in &keeptrf {
        let absorbed: Vec<usize> = group
            .iter()
            .copied()
            .filter(|&idx| idx != *rep_idx)
            .collect();
        
        // update representative's longstart/longend to encompass
        // all absorbed transfrags. This ensures transcript boundaries reflect
        // the full extent of merged transfrags.
        let mut min_longstart = transfrags[*rep_idx].longstart;
        let mut max_longend = transfrags[*rep_idx].longend;
        for &idx in &absorbed {
            if transfrags[idx].longstart > 0 {
                if min_longstart == 0 || transfrags[idx].longstart < min_longstart {
                    min_longstart = transfrags[idx].longstart;
                }
            }
            if transfrags[idx].longend > max_longend {
                max_longend = transfrags[idx].longend;
            }
        }
        if min_longstart > 0 {
            transfrags[*rep_idx].longstart = min_longstart;
        }
        if max_longend > 0 {
            transfrags[*rep_idx].longend = max_longend;
        }
        
        transfrags[*rep_idx].group = absorbed.clone();
        for &idx in &absorbed {
            if !rep_set.contains(&idx) {
                transfrags[idx].weak = 1;
            }
        }
    }

    // parse_trflong consumes from tail, so effective processing order is reverse insertion.
    for &rid in &trflong_insert {
        if rid < transfrags.len() {
            transfrags[rid].trflong_seed = true;
        }
    }

    if mixed_mode {
        // Mixed-mode long pass uses negative usepath markers.
        // Later insertions are consumed first (tail semantics), so assign more-negative to later inserts.
        for (i, rid) in trflong_insert.iter().copied().enumerate() {
            transfrags[rid].usepath = -2 - i as i32;
        }
    } else {
        let total = trflong_insert.len();
        for (i, rid) in trflong_insert.iter().copied().enumerate() {
            transfrags[rid].usepath = (total.saturating_sub(1) - i) as i32;
        }
    }
    if trace_locus.is_some() {
        for &rid in &trflong_insert {
            if rid >= transfrags.len()
                || !tf_overlaps_trace_graph(&transfrags[rid], graph, trace_locus)
            {
                continue;
            }
            let tf = &transfrags[rid];
            let (tf_s, tf_e) = tf_span_graph(tf, graph);
            let mut origin_counts: HashMap<String, usize> = Default::default();
            let rep_origin = tf
                .origin_tag
                .clone()
                .unwrap_or_else(|| "unknown".to_string());
            *origin_counts.entry(rep_origin.clone()).or_insert(0) += 1;
            for &gid in &tf.group {
                if let Some(gtf) = transfrags.get(gid) {
                    let tag = gtf
                        .origin_tag
                        .clone()
                        .unwrap_or_else(|| "unknown".to_string());
                    *origin_counts.entry(tag).or_insert(0) += 1;
                }
            }
            let mut origin_summary: Vec<(String, usize)> = origin_counts.into_iter().collect();
            origin_summary.sort_by(|a, b| b.1.cmp(&a.1).then_with(|| a.0.cmp(&b.0)));
            let origin_summary = origin_summary
                .into_iter()
                .map(|(k, v)| format!("{k}:{v}"))
                .collect::<Vec<_>>()
                .join(",");
            eprintln!(
                "[TRACE_KEEPTRF] rep={} usepath={} seed={} span={}-{} nodes={} origin={} group={:?} origin_mix=[{}]",
                rid,
                tf.usepath,
                tf.trflong_seed,
                tf_s,
                tf_e,
                tf.node_ids.len(),
                rep_origin,
                tf.group
                ,
                origin_summary
            );
            trace_tf_pattern_outgoing("keeptrf.rep_pattern", rid, tf, graph);
            if let Some(target) = trace_intron {
                if let Some((a, b)) = tf_find_intron(tf, graph, target) {
                    eprintln!(
                        "[TRACE_INTRON_TF] stage=keeptrf_final idx={} intron={}-{} edge={}=>{} usepath={} seed={} group={:?}",
                        rid,
                        target.0,
                        target.1,
                        a,
                        b,
                        tf.usepath,
                        tf.trflong_seed,
                        tf.group
                    );
                }
            }
        }
    }
    if verbose {
        eprintln!(
            "    process_transfrags: trflong seeds={} complete={} incomplete={} guides={}",
            trflong_insert.len(),
            n_complete,
            n_incomplete,
            n_guided
        );
    }

    if verbose {
        let n_absorbed = transfrags.iter().filter(|tf| tf.weak != 0).count();
        eprintln!(
            "    process_transfrags: {} active, {} keeptrf, {} absorbed",
            active_indices.len(),
            keeptrf.len(),
            n_absorbed
        );
    }
    // fill disconnected parent-child edges with trthr transfrags.
    // For each parent-child edge in the graph, if no existing transfrag covers that edge,
    // create a minimal 2-node transfrag so path extraction can route through it.
    {
        let psize = graph.pattern_size();
        // Build allpat: union of all transfrag patterns
        let mut allpat = vec![false; psize];
        for tf in transfrags.iter() {
            for bi in 0..psize {
                if tf.pattern.get_bit(bi) {
                    allpat[bi] = true;
                }
            }
        }
        let source_id = graph.source_id;
        let sink_id = graph.sink_id;
        let mut new_tfs = Vec::new();
        for i in 1..graph.n_nodes.saturating_sub(1) {
            if i == source_id || i == sink_id {
                continue;
            }
            let children: Vec<usize> = graph.nodes[i].children.ones().collect::<Vec<_>>();
            for &c in &children {
                if let Some(&edge_pos) = graph.gpos.get(&(i, c)) {
                    if edge_pos < psize && !allpat[edge_pos] {
                        let mut tf = GraphTransfrag::new(vec![i, c], psize);
                        graph.set_pattern_edges_for_path(&mut tf.pattern, &tf.node_ids);
                        tf.abundance = min_abundance;
                        tf.read_count = min_abundance;
                        if long_mode {
                            tf.longread = true;
                        }
                        new_tfs.push(tf);
                    }
                }
            }
        }
        if verbose && !new_tfs.is_empty() {
            eprintln!(
                "    process_transfrags: filled {} disconnected parent-child edges",
                new_tfs.len()
            );
        }
        // Register new transfrags with their nodes' trf_ids
        let base_idx = transfrags.len();
        for (offset, tf) in new_tfs.iter().enumerate() {
            let tf_idx = base_idx + offset;
            for &nid in &tf.node_ids {
                if nid < graph.n_nodes {
                    graph.nodes[nid].trf_ids.push(tf_idx);
                }
            }
        }
        transfrags.extend(new_tfs);
        rebuild_node_trf_ids(graph, &transfrags);
    }

    // 
    // Normalize terminal helper abundances after compatibilities are built:
    // - for nodes whose only parent is source, source->node helper abundance becomes the
    //   total abundance of the other transfrags ending/continuing at that node
    // - for sink helpers with abundance <= 1, derive abundance from node bp coverage
    let source_children: Vec<usize> = graph
        .nodes
        .get(source_id)
        .map(|n| n.children.ones().collect())
        .unwrap_or_default();
    for child in source_children {
        let Some(cnode) = graph.nodes.get(child) else {
            continue;
        };
        if cnode.parents.count_ones() != 1 || !cnode.parents.contains(source_id) {
            continue;
        }
        let mut abundance = 0.0f64;
        let mut source_tf: Option<usize> = None;
        for &t in &cnode.trf_ids {
            let Some(tf) = transfrags.get(t) else {
                continue;
            };
            if tf.node_ids.last().copied() == Some(child) {
                source_tf = Some(t);
            } else {
                abundance += tf.abundance;
            }
        }
        if let Some(t0) = source_tf {
            if transfrags
                .get(t0)
                .map(|tf| tf.abundance > 0.0)
                .unwrap_or(false)
            {
                transfrags[t0].abundance = abundance;
            }
        }
    }

    let sink_parents: Vec<usize> = graph
        .nodes
        .get(sink_id)
        .map(|n| n.parents.ones().collect())
        .unwrap_or_default();
    for parent in sink_parents {
        let Some(pnode) = graph.nodes.get(parent) else {
            continue;
        };
        let mut abundance = 0.0f64;
        let mut sink_tf: Option<usize> = None;
        for &t in &pnode.trf_ids {
            let Some(tf) = transfrags.get(t) else {
                continue;
            };
            if tf.node_ids.last().copied() == Some(sink_id) {
                sink_tf = Some(t);
            } else {
                abundance += tf.abundance;
            }
        }
        let Some(t0) = sink_tf else {
            continue;
        };
        let Some(tf0) = transfrags.get_mut(t0) else {
            continue;
        };
        if tf0.abundance > 0.0 && tf0.abundance <= 1.0 {
            let plen = pnode.length() as f64;
            if plen > 0.0 {
                tf0.abundance = (pnode.coverage / plen) - abundance;
                if tf0.abundance < 1.0 {
                    tf0.abundance = 1.0;
                }
            }
        }
    }

    if let Some(path) = keeptrf_export_path {
        if let Err(err) = write_keeptrf_usepath_tsv(
            path,
            &keeptrf,
            &trflong_insert,
            &transfrags,
            graph,
            mixed_mode,
            &keepsink,
            &hassink,
        ) {
            eprintln!("warning: failed to write keeptrf/usepath TSV {}: {}", path, err);
        }
    }

    // StringTie-parity PARITY_TF dump: emit one line per transfrag in StringTie's
    // format. Gated via RUSTLE_PARITY_TF_DUMP=1 + RUSTLE_TRACE_LOCUS=start-end.
    // Enables side-by-side comparison of transfrag sets between the two pipelines.
    if std::env::var_os("RUSTLE_PARITY_TF_DUMP").is_some() {
        let locus: Option<(u64, u64)> = std::env::var("RUSTLE_TRACE_LOCUS")
            .ok()
            .and_then(|v| {
                let (s, e) = v.split_once('-')?;
                let s = s.trim().parse::<u64>().ok()?;
                let e = e.trim().parse::<u64>().ok()?;
                Some((s, e))
            });
        let in_locus = |nid: usize| -> bool {
            match locus {
                None => false,
                Some((lo, hi)) => graph
                    .nodes
                    .get(nid)
                    .map_or(false, |n| n.end >= lo && n.start <= hi),
            }
        };
        let any_in_locus =
            transfrags.iter().any(|tf| tf.node_ids.iter().any(|&nid| in_locus(nid)));
        if any_in_locus {
            for (t, tf) in transfrags.iter().enumerate() {
                let mut coord_parts = Vec::with_capacity(tf.node_ids.len());
                for &nid in &tf.node_ids {
                    let (s, e) = graph
                        .nodes
                        .get(nid)
                        .map(|n| (n.start, n.end))
                        .unwrap_or((0, 0));
                    coord_parts.push(format!("{}({}-{})", nid, s, e));
                }
                eprintln!(
                    "PARITY_TF t={} abund={:.4} longstart={} longend={} guide={} longread={} seed={} usepath={} nodes={}: {}",
                    t,
                    tf.abundance,
                    tf.longstart,
                    tf.longend,
                    tf.guide as u8,
                    tf.longread as u8,
                    tf.trflong_seed as u8,
                    tf.usepath,
                    tf.node_ids.len(),
                    coord_parts.join(" "),
                );
            }
        }
    }

    transfrags
}

/// Coalesce trflong_seed transfrags whose intron chains differ only by
/// alt-donor/acceptor shifts within a tolerance. Targets the STRG.309
/// class: 39 seeds producing 16 j-class variants from combinations of
/// alt-splice shifts at individually well-supported junctions.
///
/// Algorithm:
/// 1. Build global donor-cluster and acceptor-cluster maps: cluster all
///    donor (resp. acceptor) coords seen in seed transfrags using a
///    greedy sweep — coords within WINDOW bp of the previous cluster
///    center collapse into the same cluster ID.
/// 2. For each seed transfrag, compute a "blurred" intron chain
///    signature: (donor_cluster_id, acceptor_cluster_id) tuple per
///    intron. First/last node indices also included so distinct
///    TSS/TTS variants stay apart.
/// 3. Group seeds by signature. Within each group, keep the one with
///    max abundance as representative; merge others' abundance into
///    it; mark merged seeds as non-seeds (trflong_seed = false) and
///    zero their abundance (so max_flow doesn't see them).
///
/// Enabled by RUSTLE_TF_COALESCE=1. Tune with:
///   RUSTLE_TF_COALESCE_WINDOW=N  (default 50 bp)
///   RUSTLE_TF_COALESCE_DEBUG=1
///
/// Returns (#seeds_before, #seeds_after) for diagnostic tracking.
pub fn coalesce_alt_junc_seed_transfrags(
    transfrags: &mut [GraphTransfrag],
    graph: &Graph,
) -> (usize, usize) {
    // Default ON at lossless settings. Opt-out via RUSTLE_TF_COALESCE_OFF.
    if std::env::var_os("RUSTLE_TF_COALESCE_OFF").is_some() {
        return (0, 0);
    }
    // Defaults: W=15 R=0.1 is lossless on GGO_19 (preserves all 1674 matches,
    // reduces j-class 211→208, F1 87.14→87.20).
    let window: u64 = std::env::var("RUSTLE_TF_COALESCE_WINDOW")
        .ok().and_then(|v| v.parse().ok()).unwrap_or(15);
    let debug = std::env::var_os("RUSTLE_TF_COALESCE_DEBUG").is_some();

    // Collect all intron (donor_coord, acceptor_coord) pairs from seed
    // transfrags. Each intron corresponds to two consecutive NON-source,
    // NON-sink nodes where node[i].end != node[i+1].start (there's a gap).
    let source_id = graph.source_id;
    let sink_id = graph.sink_id;
    let mut all_donors: Vec<u64> = Vec::new();
    let mut all_acceptors: Vec<u64> = Vec::new();
    let mut seed_ids: Vec<usize> = Vec::new();
    for (idx, tf) in transfrags.iter().enumerate() {
        if !tf.trflong_seed || tf.abundance <= 0.0 { continue; }
        seed_ids.push(idx);
        let real_path: Vec<usize> = tf.node_ids.iter().copied()
            .filter(|&n| n != source_id && n != sink_id).collect();
        for w in real_path.windows(2) {
            let a = graph.nodes.get(w[0]);
            let b = graph.nodes.get(w[1]);
            if let (Some(a), Some(b)) = (a, b) {
                if a.end < b.start {
                    // intron a.end → b.start (half-open/exclusive ends)
                    all_donors.push(a.end);
                    all_acceptors.push(b.start);
                }
            }
        }
    }
    if seed_ids.is_empty() { return (0, 0); }

    // Cluster coord list: sort, sweep, assign cluster id each time a gap
    // >= window appears. Returns map coord → cluster_id.
    let cluster_coords = |coords: &[u64]| -> HashMap<u64, u32> {
        let mut unique: Vec<u64> = coords.to_vec();
        unique.sort_unstable();
        unique.dedup();
        let mut out: HashMap<u64, u32> = Default::default();
        let mut cur_id = 0u32;
        let mut last: Option<u64> = None;
        for c in unique {
            if let Some(l) = last {
                if c.saturating_sub(l) > window {
                    cur_id += 1;
                }
            }
            out.insert(c, cur_id);
            last = Some(c);
        }
        out
    };
    let donor_cluster = cluster_coords(&all_donors);
    let acceptor_cluster = cluster_coords(&all_acceptors);

    // Build signature for each seed. Signature = (first_real_node_coord_bucket,
    // [(donor_cluster, acceptor_cluster) per intron],
    //  last_real_node_coord_bucket).
    // Buckets ensure near-equal TSS/TTS coords collide but fundamentally
    // different TSS/TTS do not.
    let bucket_size = window.max(1);
    let bucket_of = |c: u64| -> u64 { c / bucket_size };
    let mut signatures: Vec<(u32, Vec<(u32, u32)>, u32)> = Vec::with_capacity(seed_ids.len());
    for &sid in &seed_ids {
        let tf = &transfrags[sid];
        let real_path: Vec<usize> = tf.node_ids.iter().copied()
            .filter(|&n| n != source_id && n != sink_id).collect();
        let first = real_path.first().copied().and_then(|n| graph.nodes.get(n))
            .map(|n| bucket_of(n.start) as u32).unwrap_or(0);
        let last = real_path.last().copied().and_then(|n| graph.nodes.get(n))
            .map(|n| bucket_of(n.end) as u32).unwrap_or(0);
        let mut introns: Vec<(u32, u32)> = Vec::new();
        for w in real_path.windows(2) {
            let a = graph.nodes.get(w[0]);
            let b = graph.nodes.get(w[1]);
            if let (Some(a), Some(b)) = (a, b) {
                if a.end < b.start {
                    let d = donor_cluster.get(&a.end).copied().unwrap_or(u32::MAX);
                    let ac = acceptor_cluster.get(&b.start).copied().unwrap_or(u32::MAX);
                    introns.push((d, ac));
                }
            }
        }
        signatures.push((first, introns, last));
    }

    // Group seeds by signature
    let mut groups: HashMap<(u32, Vec<(u32, u32)>, u32), Vec<usize>> = Default::default();
    for (i, &sid) in seed_ids.iter().enumerate() {
        let sig = signatures[i].clone();
        groups.entry(sig).or_default().push(sid);
    }

    // For each group with >1 seed, pick max-abundance representative;
    // merge only WEAK members (abundance < min_ratio * rep.abundance).
    // Preserves legit alt-splice isoforms with similar support.
    let max_ratio: f64 = std::env::var("RUSTLE_TF_COALESCE_MAX_RATIO")
        .ok().and_then(|v| v.parse().ok()).unwrap_or(0.1);
    let before_count = seed_ids.len();
    let mut merged_count = 0usize;
    for (_sig, members) in groups {
        if members.len() <= 1 { continue; }
        // Find max-abundance representative
        let rep = members.iter().copied().max_by(|&a, &b| {
            transfrags[a].abundance
                .partial_cmp(&transfrags[b].abundance)
                .unwrap_or(std::cmp::Ordering::Equal)
        }).unwrap();
        let rep_abund_orig = transfrags[rep].abundance;
        let mut rep_abund = rep_abund_orig;
        let mut rep_read_count = transfrags[rep].read_count;
        let mut rep_group: Vec<usize> = transfrags[rep].group.clone();
        let mut rep_longstart = transfrags[rep].longstart;
        let mut rep_longend = transfrags[rep].longend;
        for &m in &members {
            if m == rep { continue; }
            let (abund, rc, ls, le) = {
                let tm = &transfrags[m];
                (tm.abundance, tm.read_count, tm.longstart, tm.longend)
            };
            // Gate: only merge if weak relative to rep (< max_ratio)
            if rep_abund_orig > 0.0 && abund >= max_ratio * rep_abund_orig {
                // Similarly-abundant — leave alone (legit alt-splice)
                continue;
            }
            rep_abund += abund;
            rep_read_count += rc;
            rep_group.push(m);
            if ls > 0 && (rep_longstart == 0 || ls < rep_longstart) {
                rep_longstart = ls;
            }
            if le > rep_longend {
                rep_longend = le;
            }
            // Zero merged transfrag to remove from future consideration
            transfrags[m].abundance = 0.0;
            transfrags[m].trflong_seed = false;
            transfrags[m].weak = 1;
            merged_count += 1;
            if debug {
                eprintln!(
                    "TF_COALESCE merge seed idx={} abund={:.2} → rep={} rep_abund={:.2}",
                    m, abund, rep, rep_abund_orig
                );
            }
        }
        transfrags[rep].abundance = rep_abund;
        transfrags[rep].read_count = rep_read_count;
        transfrags[rep].group = rep_group;
        transfrags[rep].longstart = rep_longstart;
        transfrags[rep].longend = rep_longend;
    }
    if debug && merged_count > 0 {
        eprintln!(
            "TF_COALESCE seeds before={} after={} merged={} window={}bp",
            before_count, before_count - merged_count, merged_count, window
        );
    }
    (before_count, before_count - merged_count)
}
