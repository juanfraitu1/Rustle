//! Graph-based merge of transcripts (merge mode: CMTransfrag, process_merge_transfrags,
//! merge_transfrags, printMergeResults). Builds a splice graph from the union of input exons,
//! creates one transfrag per transcript, runs process_merge_transfrags and merge_transfrags,
//! then filters and writes GTF (printMergeResults-style).

use crate::graph::{Graph, GraphTransfrag};
use crate::graph_build::create_graph;
use crate::gtf::write_gtf;
use crate::path_extract::Transcript;
use crate::transfrag_process::compatible_long;
use crate::transcript_filter::filter_min_transcript_length;
use crate::types::{CTrInfo, DetHashSet as HashSet, Junction, RunConfig};
use std::collections::HashMap as StdHashMap;
use std::io::Write;

/// Merge transfrag: wraps a graph transfrag with optional source file/name for GTF attributes
/// CMTransfrag: transfrag + read indices; for GTF merge we store source_file_idx).
#[derive(Debug, Clone)]
pub struct MergeTransfrag {
    pub transfrag: GraphTransfrag,
    /// Index of the input GTF file this transfrag came from (for input_transcripts attribute).
    pub source_file_idx: Option<usize>,
    /// Optional source transcript name.
    pub source_name: Option<String>,
}

/// naming convention (`CMTransfrag` in merge mode).
pub type CMTransfrag = MergeTransfrag;

#[inline]
fn comptbl_pos(t1: usize, t2: usize, n: usize) -> usize {
    t2 + (t1 * (2 * n - t1 - 1)) / 2
}

#[inline]
fn is_compatible_pair(t1: usize, t2: usize, n: usize, compatible: &[bool]) -> bool {
    let (a, b) = if t1 < t2 { (t1, t2) } else { (t2, t1) };
    let pos = comptbl_pos(a, b, n);
    compatible.get(pos).copied().unwrap_or(false)
}

fn max_compon_size_with_penalty_rec(
    trnumber: usize,
    set: &[CTrInfo],
    compatible: &[bool],
    mark: &[bool],
    removable: &[bool],
    computed: &mut StdHashMap<String, (f64, Vec<usize>)>,
) -> Option<(f64, Vec<usize>)> {
    // MIN_VAL sentinel (effectively "invalid").
    const MIN_VAL: f64 = -1.0e12;

    if set.is_empty() {
        return None;
    }

    let mut result: Option<Vec<usize>> = None;
    let mut maxsize: f64 = MIN_VAL;
    let mut penalty = 0.0f64;

    for i in 0..set.len() {
        let mut size = set[i].abundance - penalty;
        let mut maxagreesize = 0.0f64;
        let mut agreeset: Vec<CTrInfo> = Vec::new();
        let mut key = String::new();

        for j in (i + 1)..set.len() {
            if is_compatible_pair(set[i].trno, set[j].trno, trnumber, compatible) {
                agreeset.push(set[j]);
                key.push_str(&set[j].trno.to_string());
                key.push(if set[j].abundance != 0.0 { ' ' } else { '.' });
                key.push(if set[j].penalty != 0.0 { ' ' } else { '.' });
                maxagreesize += set[j].abundance;
            } else if mark.get(set[j].trno).copied().unwrap_or(false) {
                if removable.get(set[j].trno).copied().unwrap_or(false) {
                    size -= set[j].penalty;
                } else {
                    size = MIN_VAL;
                    break;
                }
            }
        }

        let mut agreeresult: Option<(f64, Vec<usize>)> = None;
        if size > MIN_VAL && !agreeset.is_empty() && (size + maxagreesize > maxsize) {
            if let Some(cached) = computed.get(&key).cloned() {
                agreeresult = Some(cached);
            } else if let Some((agreesize, agreed_set)) = max_compon_size_with_penalty_rec(
                trnumber,
                &agreeset,
                compatible,
                mark,
                removable,
                computed,
            ) {
                computed.insert(key, (agreesize, agreed_set.clone()));
                agreeresult = Some((agreesize, agreed_set));
            }
            if let Some((agreesize, _)) = &agreeresult {
                size += *agreesize;
            }
        }

        if size > maxsize {
            let mut picked = vec![set[i].trno];
            if let Some((_, mut agreed)) = agreeresult {
                picked.append(&mut agreed);
            }
            maxsize = size;
            result = Some(picked);
        }

        if mark.get(set[i].trno).copied().unwrap_or(false) {
            if removable.get(set[i].trno).copied().unwrap_or(false) {
                penalty += set[i].penalty;
            } else {
                break;
            }
        }
    }

    result.map(|set| (maxsize, set))
}

/// Port of `max_compon_size_with_penalty`.
/// Returns (best_size, selected transcript indices).
fn max_compon_size_with_penalty(
    trnumber: usize,
    mut set: Vec<CTrInfo>,
    compatible: &[bool],
    mark: &[bool],
    removable: &[bool],
) -> Option<(f64, Vec<usize>)> {
    // Assumes `set` sorted by trno.
    set.sort_by_key(|s| s.trno);
    let mut memo: StdHashMap<String, (f64, Vec<usize>)> = StdHashMap::new();
    max_compon_size_with_penalty_rec(trnumber, &set, compatible, mark, removable, &mut memo)
}

fn select_max_compatible_transfrags(graph: &Graph, transfrags: &[GraphTransfrag]) -> Vec<usize> {
    let n = transfrags.len();
    if n <= 1 {
        return (0..n).collect();
    }

    let tbl_size = n * (n + 1) / 2;
    let mut compatible = vec![false; tbl_size];
    for i in 0..n {
        compatible[comptbl_pos(i, i, n)] = true;
    }
    if n > 1 {
        for i in 0..n {
            for j in (i + 1)..n {
                let c1 = compatible_long(&transfrags[i], &transfrags[j], graph).0 > 0;
                let c2 = compatible_long(&transfrags[j], &transfrags[i], graph).0 > 0;
                if c1 || c2 {
                    compatible[comptbl_pos(i, j, n)] = true;
                }
            }
        }
    }

    let set: Vec<CTrInfo> = transfrags
        .iter()
        .enumerate()
        .map(|(i, tf)| CTrInfo::new(i, tf.abundance.max(0.0), 0.0))
        .collect();
    let mark = vec![false; n];
    let removable = vec![true; n];
    max_compon_size_with_penalty(n, set, &compatible, &mark, &removable)
        .map(|(_, picked)| picked)
        .unwrap_or_else(|| (0..n).collect())
}

/// Build junction set from transcript intron chains.
fn junctions_from_transcripts(transcripts: &[Transcript]) -> Vec<Junction> {
    let mut jset: HashSet<Junction> = Default::default();
    for t in transcripts {
        if t.exons.len() < 2 {
            continue;
        }
        for i in 0..t.exons.len() - 1 {
            let donor = t.exons[i].1;
            let acceptor = t.exons[i + 1].0;
            jset.insert(Junction::new(donor, acceptor));
        }
    }
    let mut out: Vec<Junction> = jset.into_iter().collect();
    out.sort_by_key(|j| (j.donor, j.acceptor));
    out
}

/// Map a genomic (start, end) segment to graph node id(s) that overlap it.
/// Returns node ids in graph order (by start).
///
/// `create_graph` lays out real exon nodes between `source_id` and `sink_id` in increasing
/// genomic start order with disjoint intervals, so we binary-search the first overlapping
/// node and walk forward until `node.start >= end` — O(log N + K) instead of O(N) per segment.
fn segment_to_node_ids(graph: &Graph, start: u64, end: u64) -> Vec<usize> {
    if start >= end {
        return Vec::new();
    }
    let src = graph.source_id;
    let sink = graph.sink_id;
    let nodes = &graph.nodes;
    if nodes.len() <= 2 || sink <= src + 1 {
        return Vec::new();
    }
    let lo = src + 1;
    let hi = sink;
    if lo >= hi {
        return Vec::new();
    }
    let slice = &nodes[lo..hi];
    if slice.is_empty() {
        return Vec::new();
    }
    let i0 = slice.partition_point(|node| node.end <= start);
    let mut out = Vec::new();
    let mut i = i0;
    while i < slice.len() {
        let node = &slice[i];
        if node.start >= end {
            break;
        }
        if node.end > start && node.start < end {
            out.push(lo + i);
        }
        i += 1;
    }
    out
}

/// Build path (ordered node ids) for a transcript in the graph.
fn transcript_to_path(graph: &Graph, exons: &[(u64, u64)]) -> Vec<usize> {
    let mut path: Vec<usize> = Vec::new();
    let mut seen: HashSet<usize> = Default::default();
    for (s, e) in exons {
        for nid in segment_to_node_ids(graph, *s, *e) {
            if seen.insert(nid) {
                path.push(nid);
            }
        }
    }
    path.sort_by_key(|&nid| graph.nodes[nid].start);
    path
}

/// Build graph from transcripts (union of exons) for one (chrom, strand) group.
fn build_merge_graph(
    _transcripts: &[Transcript],
    bundle_start: u64,
    bundle_end: u64,
    junctions: &[Junction],
) -> Graph {
    use crate::types::CBundlenode;
    let bn = CBundlenode::new(bundle_start, bundle_end, 1.0, 0);
    create_graph(
        junctions,
        bundle_start,
        bundle_end,
        Some(&bn),
        0,
        None,
        '.',  // unstranded bundle in merge mode
        None, // no junction_stats in merge mode
        None, // no bpcov in merge mode
    )
}

/// Create one GraphTransfrag per transcript (path = nodes overlapping exons).
fn transcripts_to_transfrags(graph: &Graph, transcripts: &[Transcript]) -> Vec<GraphTransfrag> {
    let mut out = Vec::with_capacity(transcripts.len());
    for t in transcripts {
        let mut sorted_path = transcript_to_path(graph, &t.exons);
        if sorted_path.is_empty() {
            continue;
        }
        let psize = graph.pattern_size();
        // Tie-break on end / node id when starts coincide (rare).
        sorted_path.sort_unstable_by(|&a, &b| {
            let na = &graph.nodes[a];
            let nb = &graph.nodes[b];
            (na.start, na.end, a).cmp(&(nb.start, nb.end, b))
        });
        let mut tf = GraphTransfrag::new(sorted_path.clone(), psize);
        graph.set_pattern_edges_for_path(&mut tf.pattern, &sorted_path);
        tf.abundance = t.coverage.max(1.0);
        tf.real = true;
        out.push(tf);
    }
    out
}

/// Convert path (node ids) to exons; merge contiguous nodes into one exon (store_merge_prediction exons).
fn path_to_exons(graph: &Graph, path: &[usize]) -> Vec<(u64, u64)> {
    let mut nodes: Vec<(u64, u64)> = Vec::new();
    for &nid in path {
        if nid == graph.source_id || nid == graph.sink_id {
            continue;
        }
        if let Some(n) = graph.nodes.get(nid) {
            if n.end > n.start {
                nodes.push((n.start, n.end));
            }
        }
    }
    if nodes.is_empty() {
        return Vec::new();
    }
    nodes.sort_by_key(|(s, _)| *s);
    let mut exons = Vec::new();
    let (mut seg_s, mut seg_e) = nodes[0];
    for (s, e) in nodes.into_iter().skip(1) {
        if s == seg_e {
            seg_e = e;
        } else {
            exons.push((seg_s, seg_e));
            seg_s = s;
            seg_e = e;
        }
    }
    exons.push((seg_s, seg_e));
    exons
}

/// merge_transfrags-style: sort by abundance, extend paths with compatible (contained) transfrags,
/// produce one prediction per "path" with aggregated coverage (merge_transfrags + store_merge_prediction).
fn merge_transfrags_to_predictions(
    graph: &Graph,
    transfrags: &mut [GraphTransfrag],
    _isofrac: f64,
) -> Vec<Transcript> {
    if transfrags.is_empty() {
        return Vec::new();
    }
    // Sort by abundance desc, then node count desc (mgtrabundCmp).
    transfrags.sort_by(|a, b| {
        let ab = b
            .abundance
            .partial_cmp(&a.abundance)
            .unwrap_or(std::cmp::Ordering::Equal);
        if ab != std::cmp::Ordering::Equal {
            return ab;
        }
        b.node_ids.len().cmp(&a.node_ids.len())
    });

    let mut predictions: Vec<Transcript> = Vec::new();
    let mut used = vec![false; transfrags.len()];
    // Parity with component selection helper:
    // keep only the best compatible transfrag component before greedy path growth.
    let selected = select_max_compatible_transfrags(graph, transfrags);
    if !selected.is_empty() {
        let keep: HashSet<usize> = selected.into_iter().collect();
        for (i, u) in used.iter_mut().enumerate() {
            if !keep.contains(&i) {
                *u = true;
            }
        }
    }

    for t1 in 0..transfrags.len() {
        if used[t1] {
            continue;
        }
        let tf1 = &transfrags[t1];
        if tf1.node_ids.is_empty()
            || tf1.node_ids[0] == graph.source_id
            || tf1.node_ids.last().copied() == Some(graph.sink_id)
        {
            continue;
        }

        let mut path = tf1.node_ids.clone();
        let mut pathpat = tf1.pattern.clone();
        let mut cov = tf1.abundance * path_len_bp(graph, &path) as f64;
        used[t1] = true;

        for t2 in (t1 + 1)..transfrags.len() {
            if used[t2] {
                continue;
            }
            let tf2 = &transfrags[t2];
            if tf2.node_ids.is_empty() || tf2.abundance <= 0.0 {
                continue;
            }
            if pathpat.contains_pattern(&tf2.pattern) {
                cov += tf2.abundance * path_len_bp(graph, &tf2.node_ids) as f64;
                used[t2] = true;
                path = merge_paths(&path, &tf2.node_ids, graph);
                pathpat.reset();
                for &nid in &path {
                    pathpat.set_bit(nid);
                }
                for i in 0..path.len().saturating_sub(1) {
                    if let Some(eid) = graph.edge_bit_index(path[i], path[i + 1]) {
                        pathpat.set_bit(eid);
                    }
                }
            }
        }

        let path_len = path_len_bp(graph, &path);
        if path_len == 0 {
            continue;
        }
        let coverage = cov / path_len as f64;
        let exons = path_to_exons(graph, &path);
        let n_exons = exons.len();
        if n_exons == 0 {
            continue;
        }
        predictions.push(Transcript {
            chrom: String::new(),
            strand: '.',
            exons,
            coverage,
            exon_cov: vec![coverage; n_exons],
            tpm: 0.0,
            fpkm: 0.0,
            source: Some("merge".to_string()),
            is_longread: false,
            longcov: 0.0,
            bpcov_cov: 0.0,
            all_strand_cov: 0.0,
            transcript_id: None,
            gene_id: None,
            ref_transcript_id: None,
            ref_gene_id: None,
            hardstart: false,
            hardend: false,
        });
    }
    predictions
}

fn path_len_bp(graph: &Graph, path: &[usize]) -> u64 {
    path.iter()
        .filter(|&&nid| nid != graph.source_id && nid != graph.sink_id)
        .filter_map(|&nid| graph.nodes.get(nid))
        .map(|n| n.end.saturating_sub(n.start))
        .sum()
}

fn merge_paths(a: &[usize], b: &[usize], graph: &Graph) -> Vec<usize> {
    let mut set: HashSet<usize> = a.iter().copied().collect();
    for &n in b {
        set.insert(n);
    }
    let mut out: Vec<usize> = set.into_iter().collect();
    out.sort_by_key(|&nid| graph.nodes.get(nid).map(|n| n.start).unwrap_or(0));
    out
}

/// printMergeResults-style: filter overlapping predictions (isofrac, single-exon), assign gene/transcript IDs, write GTF.
pub fn print_merge_results<W: Write>(
    predictions: &mut [Transcript],
    chrom: &str,
    strand: char,
    label: &str,
    isofrac: f64,
    bundledist: u64,
    min_transcript_length: u64,
    writer: &mut W,
    gene_lines: bool,
) -> std::io::Result<usize> {
    for p in predictions.iter_mut() {
        p.chrom = chrom.to_string();
        p.strand = strand;
    }
    predictions.sort_by(|a, b| {
        (
            a.exons.first().map(|e| e.0).unwrap_or(0),
            a.exons.last().map(|e| e.1).unwrap_or(0),
        )
            .cmp(&(
                b.exons.first().map(|e| e.0).unwrap_or(0),
                b.exons.last().map(|e| e.1).unwrap_or(0),
            ))
    });

    let n = predictions.len();
    let mut keep = vec![true; n];
    for i in 0..n {
        if !keep[i] {
            continue;
        }
        let pi = &predictions[i];
        let pi_end = pi.exons.last().map(|e| e.1).unwrap_or(0);
        let mut j = i + 1;
        while j < n && predictions[j].exons.first().map(|e| e.0).unwrap_or(0) <= pi_end + bundledist {
            if keep[j] && predictions[j].strand == pi.strand {
                let pj = &predictions[j];
                if pj.exons.first().map(|e| e.0).unwrap_or(0) <= pi_end {
                    if pi.exons.len() > 1 && pj.coverage < isofrac * pi.coverage {
                        keep[i] = false;
                        break;
                    }
                    if pj.exons.len() > 1 && pi.coverage < isofrac * pj.coverage {
                        keep[j] = false;
                    }
                }
            }
            j += 1;
        }
    }

    let kept_indices: Vec<usize> = (0..n).filter(|&i| keep[i]).collect();
    let mut current_end = 0u64;
    let mut gene_no = 0usize;
    let mut tr_no = 0usize;
    for &i in &kept_indices {
        let t = &mut predictions[i];
        let start = t.exons.first().map(|e| e.0).unwrap_or(0);
        let end = t.exons.last().map(|e| e.1).unwrap_or(0);
        let tlen: u64 = t.exons.iter().map(|(s, e)| e.saturating_sub(*s)).sum();
        if tlen < min_transcript_length {
            continue;
        }
        if start > current_end {
            gene_no += 1;
            tr_no = 0;
        }
        tr_no += 1;
        t.transcript_id = Some(format!("{}.{}.{}", label, gene_no, tr_no));
        t.gene_id = Some(format!("{}.{}", label, gene_no));
        current_end = current_end.max(end);
    }

    let out_list: Vec<Transcript> = predictions
        .iter()
        .filter(|t| t.transcript_id.is_some())
        .cloned()
        .collect();
    write_gtf(&out_list, writer, label)?;
    Ok(out_list.len())
}

/// process_merge_transfrags: sort/merge/absorb contained transfrags for merge mode.
fn process_merge_transfrags(graph: &mut Graph, transfrags: &mut Vec<GraphTransfrag>) {
    if transfrags.is_empty() {
        return;
    }
    // (1) Sort by node count desc (most nodes first).
    transfrags.sort_by(|a, b| b.node_ids.len().cmp(&a.node_ids.len()));
    // (2) Absorb contained: if tf[i].pattern ⊆ tf[j].pattern, merge abundance and remove i.
    let mut i = 1;
    while i < transfrags.len() {
        let mut merged = false;
        for t2 in 0..i {
            if transfrags[t2]
                .pattern
                .contains_pattern(&transfrags[i].pattern)
            {
                transfrags[t2].abundance += transfrags[i].abundance;
                transfrags[t2].srabund += transfrags[i].srabund;
                transfrags.swap_remove(i);
                merged = true;
                break;
            }
        }
        if !merged {
            i += 1;
        }
    }
    // (3) Sort by abundance desc, then node count desc (mgtrabundCmp-style).
    transfrags.sort_by(|a, b| {
        let ab = b
            .abundance
            .partial_cmp(&a.abundance)
            .unwrap_or(std::cmp::Ordering::Equal);
        if ab != std::cmp::Ordering::Equal {
            return ab;
        }
        b.node_ids.len().cmp(&a.node_ids.len())
    });
    // (4) Rebuild graph.nodes[n].trf_ids.
    for node in &mut graph.nodes {
        node.trf_ids.clear();
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

/// Run graph-based merge for one (chrom, strand) group: build graph, create transfrags,
/// process_merge_transfrags, merge_transfrags, return predictions.
pub fn run_graph_merge_for_group(
    transcripts: &[Transcript],
    chrom: &str,
    strand: char,
    config: &RunConfig,
) -> Vec<Transcript> {
    if transcripts.is_empty() {
        return Vec::new();
    }
    let junctions = junctions_from_transcripts(transcripts);
    let bundle_start = transcripts
        .iter()
        .filter_map(|t| t.exons.first())
        .map(|e| e.0)
        .min()
        .unwrap_or(0);
    let bundle_end = transcripts
        .iter()
        .filter_map(|t| t.exons.last())
        .map(|e| e.1)
        .max()
        .unwrap_or(0);
    if bundle_start >= bundle_end {
        return Vec::new();
    }

    let mut graph = build_merge_graph(transcripts, bundle_start, bundle_end, &junctions);
    let mut tfs = transcripts_to_transfrags(&graph, transcripts);
    if tfs.is_empty() {
        return Vec::new();
    }

    process_merge_transfrags(&mut graph, &mut tfs);
    tfs.retain(|tf| !tf.node_ids.is_empty());

    let isofrac = config.transcript_isofrac;
    let mut predictions = merge_transfrags_to_predictions(&graph, &mut tfs, isofrac);
    for p in &mut predictions {
        p.chrom = chrom.to_string();
        p.strand = strand;
    }
    filter_min_transcript_length(predictions, config.min_transcript_length, config.verbose)
}
