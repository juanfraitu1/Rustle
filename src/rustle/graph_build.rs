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

/// Add sink edges at nodes whose `.end` coincides with a junction
/// donor position — these are candidate alt-TTS termini.
///
/// Motivation (STRG.294.3 class): StringTie emits alt-TTS isoforms
/// that terminate at an exon boundary which is ALSO a splice donor
/// for the dominant isoform. Rustle's graph treats such nodes as
/// pure "continue via junction" sites; without an explicit sink
/// edge the alt-TTS path is unreachable.
///
/// Two modes (env gated, independent, combinable):
///
/// 1. `RUSTLE_IMPLICIT_ALT_TTS_SINK=1` — heuristic: add sink edge at
///    every Primary node whose end matches any strand-filtered
///    junction donor.
///
/// 2. `RUSTLE_ALT_TTS_ORACLE=path/to/gtf` — oracle: parse a GTF and
///    extract per-strand transcript-end coords (last exon end for
///    `+`, first exon start for `-`). Add sink edges at Rustle
///    nodes whose end matches a coord on the same strand.
///    Diagnostic tool for isolating which alt-TTS sites Rustle
///    fails to detect natively; once Rustle can produce these
///    without oracle, the oracle can be removed.
///
/// Synthetic transfrag abundance: `RUSTLE_ALT_TTS_SINK_ABUND`
/// (default 1.0).
pub fn add_alt_tts_sink_edges(
    graph: &mut Graph,
    junctions: &[Junction],
    bundle_strand: char,
    junction_stats: Option<&JunctionStats>,
) -> Vec<GraphTransfrag> {
    let mut out: Vec<GraphTransfrag> = Vec::new();
    let implicit_on = std::env::var_os("RUSTLE_IMPLICIT_ALT_TTS_SINK").is_some();
    let oracle_path = std::env::var_os("RUSTLE_ALT_TTS_ORACLE");
    if !implicit_on && oracle_path.is_none() {
        return out;
    }

    // Build the set of candidate "end-of-tx" positions we'll honor.
    let mut candidate_ends: HashSet<u64> = HashSet::default();

    if implicit_on {
        let filtered = filter_junctions_for_bundle(junctions, bundle_strand, junction_stats);
        candidate_ends.extend(filtered.iter().map(|j| j.donor));
    }
    if let Some(path) = &oracle_path {
        if let Ok(content) = std::fs::read_to_string(path.to_str().unwrap_or("")) {
            // GTF parse: graph paths flow low-coord → high-coord
            // regardless of tx strand. A tx whose genomic range ends
            // at some coord X (a "terminus" in graph-path terms) has
            // its PATH last node ending at X — so we want both `+`
            // and `-` transcripts to register their highest-coord
            // exon end as a sink anchor. Likewise, each tx's lowest-
            // coord start is a source-anchor candidate (left out of
            // this pass to keep sink-only semantics).
            let mut cur_tid: Option<String> = None;
            let mut cur_strand: char = '.';
            let mut cur_start: u64 = u64::MAX;
            let mut cur_end: u64 = 0;
            let flush = |tid: &Option<String>, strand: char,
                         _cstart: u64, cend: u64,
                         candidates: &mut HashSet<u64>| {
                if tid.is_none() { return; }
                if strand == bundle_strand && cend > 0 {
                    candidates.insert(cend);
                }
            };
            for line in content.lines() {
                if line.starts_with('#') || line.is_empty() { continue; }
                let cols: Vec<&str> = line.split('\t').collect();
                if cols.len() < 9 || cols[2] != "exon" { continue; }
                let strand = cols[6].chars().next().unwrap_or('.');
                let es: u64 = cols[3].parse().unwrap_or(0);
                let ee: u64 = cols[4].parse().unwrap_or(0);
                // tx id
                let attrs = cols[8];
                let tid = attrs.find("transcript_id \"")
                    .and_then(|idx| {
                        let rest = &attrs[idx + "transcript_id \"".len()..];
                        rest.find('"').map(|end| rest[..end].to_string())
                    })
                    .unwrap_or_default();
                if cur_tid.as_ref() != Some(&tid) {
                    flush(&cur_tid, cur_strand, cur_start, cur_end, &mut candidate_ends);
                    cur_tid = Some(tid);
                    cur_strand = strand;
                    cur_start = u64::MAX;
                    cur_end = 0;
                }
                if es < cur_start { cur_start = es; }
                if ee > cur_end { cur_end = ee; }
            }
            flush(&cur_tid, cur_strand, cur_start, cur_end, &mut candidate_ends);
        }
    }

    if std::env::var_os("RUSTLE_ALT_TTS_DEBUG").is_some() {
        eprintln!(
            "ALT_TTS candidates for strand {} : {} entries (first 5: {:?})",
            bundle_strand,
            candidate_ends.len(),
            candidate_ends.iter().take(5).collect::<Vec<_>>()
        );
    }
    if candidate_ends.is_empty() {
        return out;
    }

    let source_id = graph.source_id;
    let sink_id = graph.sink_id;
    let abundance: f64 = std::env::var("RUSTLE_ALT_TTS_SINK_ABUND")
        .ok()
        .and_then(|v| v.parse().ok())
        .unwrap_or(1.0);
    let pattern_size = graph.pattern_size();
    let n_nodes = graph.n_nodes;
    let mut additions: Vec<usize> = Vec::new();
    // Nodes matching a candidate end but already having a sink edge —
    // we still want to mark them hardend so downstream fwd_to_sink
    // treats them as alt-TTS termini even without a new synthetic tf.
    let mut hardend_only: Vec<usize> = Vec::new();
    let diag_v = std::env::var_os("RUSTLE_ALT_TTS_DEBUG").is_some();
    for i in 0..n_nodes {
        if i == source_id || i == sink_id {
            continue;
        }
        let node = &graph.nodes[i];
        if node.role != NodeRole::Primary {
            continue;
        }
        if node.end <= node.start {
            continue;
        }
        // Graph paths flow low→high coord regardless of strand;
        // a path terminator node ENDS at the tx's highest-coord exon end.
        if !candidate_ends.contains(&node.end) {
            continue;
        }
        if diag_v {
            eprintln!(
                "ALT_TTS CANDIDATE node_id={} coord={}-{} has_sink_child={}",
                i, node.start, node.end, node.children.contains(sink_id)
            );
        }
        if node.children.contains(sink_id) {
            hardend_only.push(i);
            continue;
        }
        additions.push(i);
    }
    for nid in &additions {
        let nid = *nid;
        graph.add_edge(nid, sink_id);
        let mut tf = GraphTransfrag::new(vec![nid, sink_id], pattern_size);
        tf.abundance = abundance;
        tf.longread = true;
        // Mark the node with hardend so downstream filters/extensions
        // treat it as a legitimate terminus.
        crate::bump_hs!("graph_build.rs:392:hardend");
        graph.nodes[nid].hardend = true;
        if std::env::var_os("RUSTLE_ALT_TTS_DEBUG").is_some() {
            eprintln!(
                "ALT_TTS hardend SET (new sink): node_id={} coord={}-{}",
                nid, graph.nodes[nid].start, graph.nodes[nid].end
            );
        }
        out.push(tf);
    }
    for nid in &hardend_only {
        let nid = *nid;
        crate::bump_hs!("graph_build.rs:403:hardend");
        graph.nodes[nid].hardend = true;
        if std::env::var_os("RUSTLE_ALT_TTS_DEBUG").is_some() {
            eprintln!(
                "ALT_TTS hardend SET (existing sink): node_id={} coord={}-{}",
                nid, graph.nodes[nid].start, graph.nodes[nid].end
            );
        }
    }
    if std::env::var_os("RUSTLE_ALT_TTS_DEBUG").is_some() {
        eprintln!(
            "ALT_TTS added {} sink edges on strand {}",
            additions.len(),
            bundle_strand
        );
    }
    out
}

/// Discover hardend nodes from READ-SIGNAL: terminal-read cluster ending
/// within an exon whose high-coord boundary is a canonical junction donor.
///
/// Biology: alt-TTS isoforms terminate at the same canonical exon boundary
/// that longer siblings use as a splice donor. Reads supporting the shorter
/// variant align within the exon (with jitter, not exactly at the donor),
/// so no raw read-end clusters at the donor coord. But the CLUSTER of read
/// ends WITHIN the exon indicates biological termination at that exon.
///
/// For each primary node n whose end coord matches a junction donor on the
/// bundle strand, count reads whose last exon ends within [n.start, n.end].
/// If >= min_reads, mark n.hardend = true and add sink edge if missing.
///
/// Opt-in via RUSTLE_TERMINAL_DONOR_HARDEND=1. Threshold tunable via
/// RUSTLE_TERMINAL_DONOR_MIN (default 3).
pub fn discover_terminal_donor_hardends(
    graph: &mut Graph,
    reads: &[BundleRead],
    junctions: &[Junction],
    bundle_strand: char,
    junction_stats: Option<&JunctionStats>,
) -> Vec<GraphTransfrag> {
    let mut out: Vec<GraphTransfrag> = Vec::new();
    // Default ON. Opt-out via RUSTLE_TERMINAL_DONOR_OFF=1.
    if std::env::var_os("RUSTLE_TERMINAL_DONOR_OFF").is_some() {
        return out;
    }
    let min_reads: usize = std::env::var("RUSTLE_TERMINAL_DONOR_MIN")
        .ok().and_then(|v| v.parse().ok()).unwrap_or(3);
    let source_id = graph.source_id;
    let sink_id = graph.sink_id;
    let pattern_size = graph.pattern_size();

    // Build donor set from canonical bundle-strand junctions.
    let filtered = filter_junctions_for_bundle(junctions, bundle_strand, junction_stats);
    let donor_set: HashSet<u64> = filtered.iter().map(|j| j.donor).collect();
    if donor_set.is_empty() {
        return out;
    }

    // For each primary node whose end is a donor coord, count reads whose
    // LAST exon ends within [node.start, node.end]. (Reads with more exons
    // past this coord have their last exon elsewhere, not counted.)
    let n_nodes = graph.n_nodes;
    let mut additions: Vec<usize> = Vec::new();
    let mut hardend_only: Vec<usize> = Vec::new();
    let diag = std::env::var_os("RUSTLE_TERMINAL_DONOR_DEBUG").is_some();

    // Ratio gate: terminate-here / continue-past the donor. Alt-TTS
    // biology produces a visible split: some reads end in the exon,
    // some continue past via splice. Require terminate >= ratio * continue
    // to distinguish real alt-TTS from noise truncation.
    let min_ratio: f64 = std::env::var("RUSTLE_TERMINAL_DONOR_RATIO")
        .ok().and_then(|v| v.parse().ok()).unwrap_or(1.0);
    for i in 0..n_nodes {
        if i == source_id || i == sink_id {
            continue;
        }
        let node = &graph.nodes[i];
        if node.role != NodeRole::Primary {
            continue;
        }
        if node.end <= node.start {
            continue;
        }
        if !donor_set.contains(&node.end) {
            continue;
        }
        let (ns, ne) = (node.start, node.end);

        let mut count = 0usize;       // reads ending in this exon
        let mut continue_count = 0usize; // reads with splice past this donor
        for r in reads {
            let Some(last) = r.exons.last() else { continue };
            if last.1 >= ns && last.1 <= ne {
                count += 1;
                continue;
            }
            // Does any exon of this read end at or near this donor AND have a
            // subsequent exon? That indicates read continues past via splice.
            for (ei, ex) in r.exons.iter().enumerate() {
                if ei + 1 >= r.exons.len() { break; }
                if ex.1 == ne {
                    continue_count += 1;
                    break;
                }
            }
        }
        if count < min_reads {
            continue;
        }
        if min_ratio > 0.0 && continue_count > 0 {
            let r = count as f64 / continue_count as f64;
            if r < min_ratio {
                if diag {
                    eprintln!(
                        "TERMINAL_DONOR reject-ratio node_id={} coord={}-{} term={} cont={} ratio={:.2} min={}",
                        i, ns, ne, count, continue_count, r, min_ratio
                    );
                }
                continue;
            }
        }
        if diag {
            eprintln!(
                "TERMINAL_DONOR candidate node_id={} coord={}-{} term={} cont={} has_sink_child={}",
                i, ns, ne, count, continue_count, node.children.contains(sink_id)
            );
        }
        if node.children.contains(sink_id) {
            hardend_only.push(i);
        } else {
            additions.push(i);
        }
    }

    // Whether to also ADD sink edge (not just mark hardend). Off by default —
    // adding sink edges + synthetic transfrags disturbs flow decomposition
    // across many nodes that happen to be junction donors with some read
    // termination nearby.
    // By default, set ONLY alt_tts_end (path_extract HE-gate reads this).
    // Other filters (pairwise, RI, map_reads) continue reading `hardend`
    // unchanged — this decoupling is the whole point of the flag.
    // RUSTLE_TERMINAL_DONOR_SET_HARDEND=1 opts into setting hardend too
    // (the old behavior, kept for comparison / debugging).
    let set_hardend_too = std::env::var_os("RUSTLE_TERMINAL_DONOR_SET_HARDEND").is_some();
    let add_edges = std::env::var_os("RUSTLE_TERMINAL_DONOR_ADD_SINK").is_some();
    for &nid in &additions {
        graph.nodes[nid].alt_tts_end = true;
        if set_hardend_too {
            crate::bump_hs!("graph_build.rs:551:hardend");
            graph.nodes[nid].hardend = true;
        }
        if add_edges {
            graph.add_edge(nid, sink_id);
            let mut tf = GraphTransfrag::new(vec![nid, sink_id], pattern_size);
            tf.abundance = 1.0;
            tf.longread = true;
            out.push(tf);
        }
        if diag {
            eprintln!(
                "TERMINAL_DONOR alt_tts_end SET (hardend={} add_edge={}): node_id={} coord={}-{}",
                set_hardend_too, add_edges,
                nid, graph.nodes[nid].start, graph.nodes[nid].end
            );
        }
    }
    let _ = pattern_size;
    for &nid in &hardend_only {
        graph.nodes[nid].alt_tts_end = true;
        if set_hardend_too {
            crate::bump_hs!("graph_build.rs:572:hardend");
            graph.nodes[nid].hardend = true;
        }
        if diag {
            eprintln!(
                "TERMINAL_DONOR alt_tts_end SET (existing sink, hardend={}): node_id={} coord={}-{}",
                set_hardend_too, nid, graph.nodes[nid].start, graph.nodes[nid].end
            );
        }
    }
    out
}

/// Diagnostic MISSED-TX ORACLE: inject ref tx (from a StringTie GTF) as
/// high-abundance long-read seeds into the transfrag list. Classifies
/// per-tx why Rustle fails to emit: unmappable (exon has no matching
/// graph node = graph-build gap), partial (some exons missing), or
/// mappable but not-emitted (path-extract or filter gap).
///
/// Usage:
///   RUSTLE_MISSED_ORACLE=/path/to/stringtie.gtf \
///   [RUSTLE_MISSED_ORACLE_ONLY=tx1,tx2,...] \
///   [RUSTLE_MISSED_ORACLE_DEBUG=1] \
///   ./target/release/rustle -L in.bam -o out.gtf
///
/// If RUSTLE_MISSED_ORACLE_ONLY is set, restrict to that tx subset.
/// The injected transfrags have guide_tid = "oracle:<tid>" so downstream
/// emission can be traced to their oracle source.
pub fn inject_missed_tx_oracle(
    graph: &Graph,
    transfrags: &mut Vec<GraphTransfrag>,
    bundle_chrom: &str,
    bundle_start: u64,
    bundle_end: u64,
    bundle_strand: char,
) -> usize {
    let path = match std::env::var("RUSTLE_MISSED_ORACLE") {
        Ok(p) => p,
        Err(_) => return 0,
    };
    let content = match std::fs::read_to_string(&path) {
        Ok(c) => c,
        Err(_) => return 0,
    };
    let only: Option<std::collections::HashSet<String>> = std::env::var("RUSTLE_MISSED_ORACLE_ONLY")
        .ok()
        .map(|v| v.split(',').map(|s| s.trim().to_string()).collect());
    let debug = std::env::var_os("RUSTLE_MISSED_ORACLE_DEBUG").is_some();
    let abundance: f64 = std::env::var("RUSTLE_MISSED_ORACLE_ABUND")
        .ok().and_then(|v| v.parse().ok()).unwrap_or(100.0);

    let source_id = graph.source_id;
    let sink_id = graph.sink_id;
    let pattern_size = graph.pattern_size();

    // Parse GTF grouped by transcript
    let mut cur_tid: Option<String> = None;
    let mut cur_chrom: String = String::new();
    let mut cur_strand: char = '.';
    let mut cur_exons: Vec<(u64, u64)> = Vec::new();

    let mut injected = 0usize;

    let flush = |tid: &Option<String>,
                     chrom: &str,
                     strand: char,
                     exons: &[(u64, u64)],
                     transfrags: &mut Vec<GraphTransfrag>,
                     injected: &mut usize| {
        let Some(tid) = tid else { return };
        if let Some(only) = &only {
            if !only.contains(tid) { return; }
        }
        if chrom != bundle_chrom { return; }
        if exons.is_empty() { return; }
        let first = exons[0].0;
        let last = exons.last().unwrap().1;
        if first > bundle_end || last < bundle_start { return; }
        if strand != bundle_strand { return; }

        // Map each exon to graph nodes whose coord range is contained in or
        // overlaps the exon. Report per-exon hits to diagnose "unmappable".
        let mut nids: Vec<usize> = Vec::new();
        let mut exon_hit_counts: Vec<usize> = Vec::with_capacity(exons.len());
        for (es, ee) in exons {
            let mut this_exon_hits = 0;
            for (i, n) in graph.nodes.iter().enumerate() {
                if i == source_id || i == sink_id { continue; }
                if n.end <= n.start { continue; }
                // Inclusive GTF coords [es, ee]; node coords half-open [start, end).
                // Accept if node is fully inside OR node-end matches exon-end.
                let contained = n.start >= *es && n.end <= *ee + 1;
                let overlaps = n.start < *ee + 1 && n.end > *es;
                if contained || overlaps {
                    nids.push(i);
                    this_exon_hits += 1;
                }
            }
            exon_hit_counts.push(this_exon_hits);
        }
        let missing_exons = exon_hit_counts.iter().filter(|&&c| c == 0).count();
        let classification = if missing_exons == exons.len() {
            "UNMAPPABLE"
        } else if missing_exons > 0 {
            "PARTIAL"
        } else {
            "MAPPABLE"
        };

        if debug {
            eprintln!(
                "ORACLE_CLASSIFY tid={} bundle={}:{}-{} strand={} class={} exons_total={} exons_missing={} n_nodes_chained={}",
                tid, bundle_chrom, bundle_start, bundle_end, strand,
                classification, exons.len(), missing_exons, nids.len()
            );
        }
        if classification == "UNMAPPABLE" {
            return; // nothing to inject
        }

        nids.sort();
        nids.dedup();
        if nids.is_empty() { return; }

        let mut tf = GraphTransfrag::new(nids, pattern_size);
        tf.abundance = abundance;
        tf.trflong_seed = true;
        tf.longread = true;
        tf.real = true;
        tf.guide_tid = Some(format!("oracle:{tid}"));
        tf.longstart = first;
        tf.longend = last;
        transfrags.push(tf);
        *injected += 1;
    };

    for line in content.lines() {
        if line.starts_with('#') || line.is_empty() { continue; }
        let cols: Vec<&str> = line.split('\t').collect();
        if cols.len() < 9 { continue; }
        if cols[2] != "exon" { continue; }
        let chrom = cols[0];
        let strand = cols[6].chars().next().unwrap_or('.');
        let es: u64 = cols[3].parse().unwrap_or(0);
        let ee: u64 = cols[4].parse().unwrap_or(0);
        let attrs = cols[8];
        let tid = attrs.find("transcript_id \"")
            .and_then(|idx| {
                let rest = &attrs[idx + "transcript_id \"".len()..];
                rest.find('"').map(|end| rest[..end].to_string())
            })
            .unwrap_or_default();
        if cur_tid.as_ref() != Some(&tid) {
            flush(&cur_tid, &cur_chrom, cur_strand, &cur_exons,
                  transfrags, &mut injected);
            cur_tid = Some(tid);
            cur_chrom = chrom.to_string();
            cur_strand = strand;
            cur_exons.clear();
        }
        cur_exons.push((es, ee));
    }
    flush(&cur_tid, &cur_chrom, cur_strand, &cur_exons,
          transfrags, &mut injected);

    if debug && injected > 0 {
        eprintln!(
            "ORACLE_INJECTED total={} bundle={}:{}-{} strand={}",
            injected, bundle_chrom, bundle_start, bundle_end, bundle_strand
        );
    }
    injected
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
    let parity_on = crate::parity_decisions::is_enabled();
    let pre: Vec<&'a Junction> = junctions
        .iter()
        .filter(|j| {
            let mut reject_reason: Option<&'static str> = None;
            if let Some(js) = junction_stats {
                if let Some(stat) = js.get(j) {
                    if let Some(ts) = target_strand {
                        if stat.strand != Some(ts) {
                            reject_reason = Some("strand_mismatch");
                        }
                    }
                    if reject_reason.is_none() && stat.mm < 0.0 {
                        reject_reason = Some("mm_negative");
                    }
                    if reject_reason.is_none() && stat.strand == Some(0) {
                        reject_reason = Some("strand_zero");
                    }
                }
            }
            let accepted = reject_reason.is_none();
            if parity_on {
                let mm = junction_stats.and_then(|js| js.get(j)).map(|s| s.mm).unwrap_or(0.0);
                let strand_byte = junction_stats
                    .and_then(|js| js.get(j))
                    .and_then(|s| s.strand)
                    .map(|x| x as i64)
                    .unwrap_or(-1);
                let payload = format!(
                    r#""accepted":{},"bundle_strand":"{}","jstrand":{},"mm":{:.4},"reason":"{}""#,
                    accepted, bundle_strand, strand_byte, mm,
                    reject_reason.unwrap_or("ok"),
                );
                // Convention normalization: rustle stores 0-based half-open exon [s, e),
                // so j.donor = last-base+1 (0-based exclusive donor exon end) and j.acceptor =
                // first base of next exon (0-based). StringTie uses 1-based inclusive: jstart =
                // last-base of donor exon, jend = first base of acceptor exon. Donor coincidentally
                // matches numerically (0-based exclusive end == 1-based inclusive last base).
                // Acceptor needs +1 to align (rustle 0-based first base → StringTie 1-based first base).
                // Junction has no chrom field — emit None and rely on coord match for cross-tool diff.
                crate::parity_decisions::emit(
                    "junction_accept",
                    None,
                    j.donor,
                    j.acceptor + 1,
                    bundle_strand,
                    &payload,
                );
            }
            accepted
        })
        .collect();
    // Alt-junction coalescing: demote weak alt-donors/acceptors to reduce
    // graph over-segmentation. Targets the STRG.309 class of j-class
    // over-emission where Rustle explores 2^N combinations of alt-splice
    // shifts. Opt-in via RUSTLE_ALT_JUNC_COALESCE=1.
    if std::env::var_os("RUSTLE_ALT_JUNC_COALESCE").is_some() && junction_stats.is_some() {
        coalesce_weak_alt_junctions(pre, junction_stats.unwrap())
    } else {
        pre
    }
}

/// Drop alt-donor and alt-acceptor junctions that cluster within N bp of
/// a canonical (highest-read) junction AND have < min_ratio * canonical
/// read support. Reduces graph over-segmentation at alt-splice hot spots.
///
/// Algorithm:
/// 1. For each donor coord D, find the canonical junction (D, A_canon)
///    = max-read junction starting at D. Drop junctions (D, A') where
///    |A' - A_canon| <= window AND reads(D, A') < min_ratio * reads(D, A_canon).
///    This demotes weak alt-acceptors at the same donor.
/// 2. For each acceptor coord A, find the canonical junction (D_canon, A)
///    = max-read junction ending at A. Drop (D', A) where
///    |D' - D_canon| <= window AND reads(D', A) < min_ratio * reads(D_canon, A).
///    This demotes weak alt-donors at the same acceptor.
///
/// Does NOT drop junctions whose donor AND acceptor both differ from any
/// canonical — those are legitimate novel splice events to preserve.
///
/// Env vars:
///   RUSTLE_ALT_JUNC_WINDOW (default 100) — bp tolerance
///   RUSTLE_ALT_JUNC_MIN_RATIO (default 0.5) — alt must have at least this
///                                              fraction of canonical reads to be kept
fn coalesce_weak_alt_junctions<'a>(
    juncs: Vec<&'a Junction>,
    stats: &JunctionStats,
) -> Vec<&'a Junction> {
    let window: u64 = std::env::var("RUSTLE_ALT_JUNC_WINDOW")
        .ok().and_then(|v| v.parse().ok()).unwrap_or(100);
    let min_ratio: f64 = std::env::var("RUSTLE_ALT_JUNC_MIN_RATIO")
        .ok().and_then(|v| v.parse().ok()).unwrap_or(0.5);
    // Require minimum number of siblings (donors at same acceptor, or acceptors
    // at same donor) before demotion fires. STRG.309 class has 3-4+ siblings;
    // isolated alt-donor pairs are usually legitimate. Default 2 preserves prior
    // behavior; tuning to 3 targets the combinatorial-explosion pattern.
    let min_siblings: usize = std::env::var("RUSTLE_ALT_JUNC_MIN_SIBLINGS")
        .ok().and_then(|v| v.parse().ok()).unwrap_or(2);
    let debug = std::env::var_os("RUSTLE_ALT_JUNC_DEBUG").is_some();

    // Helper: get read support (mrcount) for a junction; 0 if no stats.
    let reads_of = |j: &Junction| -> f64 {
        stats.get(j).map(|s| s.mrcount.max(s.nreads_good)).unwrap_or(0.0)
    };

    // Group by acceptor: find canonical donor for each acceptor, drop weak alt-donors.
    use std::collections::HashMap as StdHashMap;
    let mut by_acceptor: StdHashMap<u64, Vec<&Junction>> = StdHashMap::new();
    for &j in &juncs {
        by_acceptor.entry(j.acceptor).or_default().push(j);
    }
    // Demote set: junctions to remove.
    let mut demoted: std::collections::HashSet<(u64, u64)> = Default::default();
    for (_acceptor, siblings) in &by_acceptor {
        if siblings.len() < min_siblings { continue; }
        // Find max-read canonical
        let canon = siblings.iter().copied().max_by(|a, b| {
            reads_of(a).partial_cmp(&reads_of(b))
                .unwrap_or(std::cmp::Ordering::Equal)
        });
        let Some(canon) = canon else { continue; };
        let canon_r = reads_of(canon);
        if canon_r <= 0.0 { continue; }
        for &alt in siblings {
            if alt.donor == canon.donor { continue; }
            let alt_r = reads_of(alt);
            let donor_dist = alt.donor.abs_diff(canon.donor);
            if donor_dist <= window && alt_r < min_ratio * canon_r {
                demoted.insert((alt.donor, alt.acceptor));
                if debug {
                    eprintln!(
                        "ALT_JUNC demote_at_acceptor alt={}-{} alt_reads={:.0} canon={}-{} canon_reads={:.0}",
                        alt.donor, alt.acceptor, alt_r,
                        canon.donor, canon.acceptor, canon_r
                    );
                }
            }
        }
    }

    // Group by donor: find canonical acceptor for each donor, drop weak alt-acceptors.
    let mut by_donor: StdHashMap<u64, Vec<&Junction>> = StdHashMap::new();
    for &j in &juncs {
        by_donor.entry(j.donor).or_default().push(j);
    }
    for (_donor, siblings) in &by_donor {
        if siblings.len() < min_siblings { continue; }
        let canon = siblings.iter().copied().max_by(|a, b| {
            reads_of(a).partial_cmp(&reads_of(b))
                .unwrap_or(std::cmp::Ordering::Equal)
        });
        let Some(canon) = canon else { continue; };
        let canon_r = reads_of(canon);
        if canon_r <= 0.0 { continue; }
        for &alt in siblings {
            if alt.acceptor == canon.acceptor { continue; }
            let alt_r = reads_of(alt);
            let acc_dist = alt.acceptor.abs_diff(canon.acceptor);
            if acc_dist <= window && alt_r < min_ratio * canon_r {
                demoted.insert((alt.donor, alt.acceptor));
                if debug {
                    eprintln!(
                        "ALT_JUNC demote_at_donor alt={}-{} alt_reads={:.0} canon={}-{} canon_reads={:.0}",
                        alt.donor, alt.acceptor, alt_r,
                        canon.donor, canon.acceptor, canon_r
                    );
                }
            }
        }
    }

    let before = juncs.len();
    let kept: Vec<&Junction> = juncs.into_iter()
        .filter(|j| !demoted.contains(&(j.donor, j.acceptor)))
        .collect();
    if debug {
        eprintln!(
            "ALT_JUNC coalesce: {} → {} ({} demoted, window={}bp ratio={:.2})",
            before, kept.len(), before - kept.len(), window, min_ratio
        );
    }
    kept
}

fn collect_bundlenode_events(
    filtered_juncs: &[&Junction],
    currentstart: u64,
    endbundle: u64,
    bundle_end: u64,
    demoted_donors: &HashSet<u64>,
    demoted_acceptors: &HashSet<u64>,
) -> Vec<(u64, JunctionEventType)> {
    let mut events: Vec<(u64, JunctionEventType)> = Vec::new();
    for j in filtered_juncs {
        if j.donor > currentstart && j.donor <= endbundle {
            if j.acceptor <= bundle_end {
                // Skip Start (donor) events at demoted alt-donor coords.
                // The junction still contributes as an edge from the node
                // containing this coord to the acceptor node, but we don't
                // create an extra graph-node boundary here. Reduces alt-
                // splice combinatorial path exploration at STRG.309 class.
                if !demoted_donors.contains(&j.donor) {
                    events.push((j.donor, JunctionEventType::Start));
                }
            }
        }
        if j.acceptor >= currentstart && j.acceptor <= endbundle {
            if !demoted_acceptors.contains(&j.acceptor) {
                events.push((j.acceptor, JunctionEventType::End));
            }
        }
    }
    events.sort_unstable();
    events.dedup_by_key(|e| (e.0, e.1));
    events
}

/// Compute alt-donor and alt-acceptor coords that should NOT create
/// graph-node boundaries. Targets the STRG.309 over-segmentation class:
/// at an acceptor shared by multiple donors, if an alt-donor has much
/// less read support than the canonical donor AND is within N bp of it,
/// demote the alt-donor coord so no extra graph segment is created.
///
/// Reads still map to junctions via exon boundaries; the demoted coord
/// just doesn't split the node. This keeps reads using alt-donors
/// coalesced into the same transfrag pattern as reads using canonical.
///
/// Opt-in via RUSTLE_GRAPH_ALT_COALESCE=1.
///
/// Skips junctions marked guide_match (preserves annotation-driven splits).
fn compute_demoted_alt_coords(
    filtered_juncs: &[&Junction],
    stats: &JunctionStats,
) -> (HashSet<u64>, HashSet<u64>) {
    let mut demoted_donors: HashSet<u64> = Default::default();
    let mut demoted_acceptors: HashSet<u64> = Default::default();
    if std::env::var_os("RUSTLE_GRAPH_ALT_COALESCE").is_none() {
        return (demoted_donors, demoted_acceptors);
    }
    let window: u64 = std::env::var("RUSTLE_GRAPH_ALT_WINDOW")
        .ok().and_then(|v| v.parse().ok()).unwrap_or(50);
    let min_ratio: f64 = std::env::var("RUSTLE_GRAPH_ALT_MIN_RATIO")
        .ok().and_then(|v| v.parse().ok()).unwrap_or(0.3);
    let debug = std::env::var_os("RUSTLE_GRAPH_ALT_DEBUG").is_some();

    let reads_of = |j: &Junction| -> f64 {
        stats.get(j).map(|s| s.mrcount.max(s.nreads_good)).unwrap_or(0.0)
    };
    let is_guide = |j: &Junction| -> bool {
        stats.get(j).map(|s| s.guide_match).unwrap_or(false)
    };

    use std::collections::HashMap as StdHashMap;

    // Demote alt-donors at the same acceptor.
    let mut by_acceptor: StdHashMap<u64, Vec<&Junction>> = StdHashMap::new();
    for &j in filtered_juncs {
        by_acceptor.entry(j.acceptor).or_default().push(j);
    }
    for (_acc, sibs) in &by_acceptor {
        if sibs.len() < 2 { continue; }
        let canon = sibs.iter().copied().max_by(|a, b| {
            reads_of(a).partial_cmp(&reads_of(b))
                .unwrap_or(std::cmp::Ordering::Equal)
        });
        let Some(canon) = canon else { continue; };
        if is_guide(canon) { continue; }
        let canon_r = reads_of(canon);
        if canon_r <= 0.0 { continue; }
        for &alt in sibs {
            if alt.donor == canon.donor { continue; }
            if is_guide(alt) { continue; }
            let alt_r = reads_of(alt);
            let dist = alt.donor.abs_diff(canon.donor);
            if dist <= window && alt_r < min_ratio * canon_r {
                demoted_donors.insert(alt.donor);
                if debug {
                    eprintln!(
                        "GRAPH_ALT demote_donor={} alt_reads={:.0} canon_donor={} canon_reads={:.0} dist={}bp",
                        alt.donor, alt_r, canon.donor, canon_r, dist
                    );
                }
            }
        }
    }

    // Demote alt-acceptors at the same donor.
    let mut by_donor: StdHashMap<u64, Vec<&Junction>> = StdHashMap::new();
    for &j in filtered_juncs {
        by_donor.entry(j.donor).or_default().push(j);
    }
    for (_d, sibs) in &by_donor {
        if sibs.len() < 2 { continue; }
        let canon = sibs.iter().copied().max_by(|a, b| {
            reads_of(a).partial_cmp(&reads_of(b))
                .unwrap_or(std::cmp::Ordering::Equal)
        });
        let Some(canon) = canon else { continue; };
        if is_guide(canon) { continue; }
        let canon_r = reads_of(canon);
        if canon_r <= 0.0 { continue; }
        for &alt in sibs {
            if alt.acceptor == canon.acceptor { continue; }
            if is_guide(alt) { continue; }
            let alt_r = reads_of(alt);
            let dist = alt.acceptor.abs_diff(canon.acceptor);
            if dist <= window && alt_r < min_ratio * canon_r {
                demoted_acceptors.insert(alt.acceptor);
                if debug {
                    eprintln!(
                        "GRAPH_ALT demote_acceptor={} alt_reads={:.0} canon_acceptor={} canon_reads={:.0} dist={}bp",
                        alt.acceptor, alt_r, canon.acceptor, canon_r, dist
                    );
                }
            }
        }
    }

    if debug && !(demoted_donors.is_empty() && demoted_acceptors.is_empty()) {
        eprintln!(
            "GRAPH_ALT demoted donors={} acceptors={} (W={}bp R={:.2})",
            demoted_donors.len(), demoted_acceptors.len(), window, min_ratio
        );
    }
    (demoted_donors, demoted_acceptors)
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
    let (demoted_donors, demoted_acceptors) = junction_stats
        .map(|s| compute_demoted_alt_coords(&filtered_juncs, s))
        .unwrap_or_else(|| (Default::default(), Default::default()));
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

        let events = collect_bundlenode_events(&filtered_juncs, bn.start, bn.end, bundle_end,
            &demoted_donors, &demoted_acceptors);
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
    // In pure-overlap mode, Primary disjoint splits don't exist and
    // read routing must include overlap nodes (OverlapAnchor,
    // JunctionEntry). Gated behind RUSTLE_PURE_OVERLAP=1.
    let pure_overlap = std::env::var_os("RUSTLE_PURE_OVERLAP").is_some();
    for (nid, node) in graph.nodes.iter().enumerate() {
        if nid == graph.source_id || nid == graph.sink_id {
            continue;
        }
        let accept = if pure_overlap {
            node.role.accepts_reads_pure_overlap()
        } else {
            node.role.accepts_reads()
        };
        if !accept {
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
    let mut _sink_futuretr: Vec<(usize, f64)> = Vec::new();
    create_graph_inner(junctions, _bundle_start, bundle_end, bundlenodes,
        junction_support, _reads, bundle_strand, junction_stats, bpcov, None, &[], &[],
        &mut _sink_futuretr, "")
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
    sink_futuretr_out: &mut Vec<(usize, f64)>,
    bundle_chrom: &str,
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
    // Graph-segment alt-coalesce: compute demoted alt-donor/acceptor coords.
    // When enabled via RUSTLE_GRAPH_ALT_COALESCE=1, weak alt positions within
    // N bp of a canonical donor/acceptor (at the same other-endpoint) do NOT
    // create graph-node boundaries. Reduces combinatorial path exploration.
    let (demoted_donors, demoted_acceptors) = junction_stats
        .map(|s| compute_demoted_alt_coords(&filtered_juncs, s))
        .unwrap_or_else(|| (Default::default(), Default::default()));

    // ends[] hashmap: tracks junction endpoints for cross-junction linking.
    // Persists across bundlenodes ( 4080-4088).
    let mut ends: HashMap<u64, Vec<usize>> = Default::default();

    let mut sink_parents: Vec<usize> = Vec::new();
    let sink_futuretr: &mut Vec<(usize, f64)> = sink_futuretr_out;
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
            collect_bundlenode_events(&filtered_juncs, currentstart, endbundle, bundle_end,
                &demoted_donors, &demoted_acceptors);


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
        let mut bnode_lstart_owned: Vec<ReadBoundary> = Vec::new();
        let mut bnode_lend_owned: Vec<ReadBoundary> = Vec::new();
        let mut use_bnode_boundaries = false;
        // compute boundaries per-bundlenode when longtrim is active
        // and no external boundaries were provided. External boundaries come from
        // CPAS/poly-A evidence; when absent, use coverage-derivative detection.
        let bnode_span = endbundle.saturating_sub(currentstart);
        let min_lt_span = 2 * (100 + 50) + 50; // 2*(CHI_WIN+CHI_THR) + margin = 350
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
                    use_bnode_boundaries = true;
            }
        }
        let (effective_lstart, effective_lend): (&[ReadBoundary], &[ReadBoundary]) =
            if use_bnode_boundaries {
                (bnode_lstart_owned.as_slice(), bnode_lend_owned.as_slice())
            } else {
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
                        _bundle_start, bundle_end, bundle_chrom, nodeend,
                        has_junction_start, has_junction_end, source_bid,
                        &mut sink_parents,
                        sink_futuretr,
                    );
                }
            }

            match ev_type {
                JunctionEventType::Start => {
                    // Junction start at donor position
                    // Set current graphnode end to donor.
                    let pure_overlap_active = std::env::var_os("RUSTLE_PURE_OVERLAP").is_some();
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
                    // In pure-overlap mode, the donor is a shared exon
                    // boundary across ALL parallels in this bundlenode
                    // (OverlapAnchor + each JunctionEntry). Each parallel's
                    // end must be shrunk to pos too.
                    if pure_overlap_active {
                        for &par in &bundlenode_aliases {
                            if par == graphnode_id { continue; }
                            if par < graph.nodes.len() && graph.nodes[par].end > pos {
                                graph.nodes[par].end = pos;
                            }
                        }
                    }

                    // Record all junctions starting at this position in ends[].
                    // In pure-overlap mode, register ALL parallels so downstream
                    // bundlenodes link to each; otherwise only graphnode_id.
                    for j in &filtered_juncs {
                        if j.donor == pos {
                            if pure_overlap_active && !bundlenode_aliases.is_empty() {
                                for &par in &bundlenode_aliases {
                                    ends.entry(j.acceptor).or_default().push(par);
                                }
                                // Also register graphnode_id if not already in aliases
                                if !bundlenode_aliases.contains(&graphnode_id) {
                                    ends.entry(j.acceptor).or_default().push(graphnode_id);
                                }
                            } else {
                                ends.entry(j.acceptor).or_default().push(graphnode_id);
                            }
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

                    // PURE OVERLAP mode (RUSTLE_PURE_OVERLAP=1):
                    // Don't shrink the current node; create a JunctionEntry
                    // parallel spanning [acceptor, endbundle) and link
                    // ends[acceptor] to it. Initial bundlenode node becomes
                    // an OverlapAnchor (set on first acceptor only).
                    if std::env::var_os("RUSTLE_PURE_OVERLAP").is_some()
                        && graph.nodes[graphnode_id].start < acceptor
                        && acceptor < endbundle
                    {
                        // Promote the initial bundlenode node (still Primary)
                        // to OverlapAnchor the first time we fork.
                        if graph.nodes[graphnode_id].role == NodeRole::Primary
                            && !bundlenode_aliases.contains(&graphnode_id)
                        {
                            graph.nodes[graphnode_id].role = NodeRole::OverlapAnchor;
                            bundlenode_aliases.push(graphnode_id);
                        }
                        // Create JunctionEntry [acceptor, endbundle)
                        let entry = graph.add_node(acceptor, endbundle);
                        entry.source_bnode = Some(source_bid);
                        entry.coverage = 0.0;
                        entry.role = NodeRole::JunctionEntry;
                        let entry_id = entry.node_id;
                        create_graph_new_node(
                            trace_s,
                            trace_g,
                            entry_id,
                            acceptor,
                            endbundle,
                            "junction_entry(pure_overlap)",
                        );
                        // Link junction parents → entry
                        if let Some(parent_ids) = ends.get(&acceptor).cloned() {
                            for pid in parent_ids {
                                graph.add_edge(pid, entry_id);
                            }
                        }
                        // Track entry for end-of-bundlenode children patch
                        bundlenode_aliases.push(entry_id);
                        // Advance graphnode_id to new entry so subsequent donor
                        // events attach here (acceptors fork further from entry).
                        graphnode_id = entry_id;
                        continue;
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
                    _bundle_start, bundle_end, bundle_chrom, endbundle,
                    true, true, source_bid,
                    &mut sink_parents,
                    sink_futuretr,
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
            // In pure-overlap mode, all parallels share this bundle-end
            // terminus; each must also be a sink parent so downstream
            // bundlenodes can link to any of them via ends[].
            if std::env::var_os("RUSTLE_PURE_OVERLAP").is_some() {
                for &par in &bundlenode_aliases {
                    if par != graphnode_id && par < graph.nodes.len() {
                        sink_parents.push(par);
                    }
                }
            }
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
                        crate::bump_hs!("graph_build.rs:2013:hardstart");
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
                        crate::bump_hs!("graph_build.rs:2055:hardend");
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
    bundle_end: u64,
    bundle_chrom: &str,
    nodeend: u64,
    startcov: bool, // true if junction exists at current node start
    endcov: bool,   // true if junction exists at nodeend
    source_bid: usize,
    sink_parents: &mut Vec<usize>,
    // Sink-connector futuretr: (prev_node, tmpcov) entries for synthetic 2-node
    // [prev, sink] transfrags. StringTie creates these inside its longtrim()
    // via futuretr with n2=-1; they populate hassink[prev] and make the
    // terminal node "complete" at trflong.Add time so chimeric tfs ending at
    // this node SKIP pass 1 (matching StringTie behavior).
    sink_futuretr: &mut Vec<(usize, f64)>,
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

    // RUSTLE_PARITY_LONGTRIM_SPLIT_TSV: emit per-event longtrim decisions
    // for cross-tool diff against ST's PARITY_LONGTRIM_SPLIT_TSV.
    // Schema matches ST: source, chrom, bundle_start, bundle_end, sno,
    // bnode_start, bnode_end, pos, direction, abundance, tmpcov, accepted, reason.
    let lt_emit = |graph: &Graph, graphnode_id: usize, direction: &str, pos: u64,
                   abundance: f64, tmpcov: f64, accepted: u8, reason: &str| {
        if let Ok(p) = std::env::var("RUSTLE_PARITY_LONGTRIM_SPLIT_TSV") {
            if !p.is_empty() {
                use std::io::Write;
                use std::sync::{Mutex, OnceLock};
                static WR: OnceLock<Mutex<Option<std::fs::File>>> = OnceLock::new();
                static HDR: OnceLock<Mutex<bool>> = OnceLock::new();
                let wr = WR.get_or_init(|| Mutex::new(
                    std::fs::OpenOptions::new()
                        .create(true).append(true).open(&p).ok()));
                let hdr = HDR.get_or_init(|| Mutex::new(false));
                if let (Ok(mut f_opt), Ok(mut hdr_w)) = (wr.lock(), hdr.lock()) {
                    if let Some(f) = f_opt.as_mut() {
                        if !*hdr_w {
                            let _ = writeln!(f, "source\tchrom\tbundle_start\tbundle_end\tsno\tbnode_start\tbnode_end\tpos\tdirection\tabundance\ttmpcov\taccepted\treason");
                            *hdr_w = true;
                        }
                        let sno: i32 = match bundle_strand {
                            '-' => 0,
                            '+' => 2,
                            _ => 1,
                        };
                        let bn_start = graph.nodes[graphnode_id].start;
                        let bn_end = graph.nodes[graphnode_id].end.saturating_sub(1);
                        let _ = writeln!(f,
                            "rustle\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{:.4}\t{:.4}\t{}\t{}",
                            bundle_chrom, bundle_start, bundle_end, sno,
                            bn_start, bn_end, pos, direction,
                            abundance, tmpcov, accepted, reason);
                    }
                }
            }
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
        // Rustle default: positions with <3 reads are alignment jitter, not
        // real TSS/TES. StringTie has NO such gate — any lend/lstart within
        // the node is considered, decided solely on tmpcov > 0.
        //
        // RUSTLE_LONGTRIM_STRICT_LEND=1 drops both gates to match ST semantics
        // (targets the 5096 last_node_end drift cases measured in session 5,
        // 2026-04-24, where Rustle's last graph node extends 20-100bp past ST's
        // because single-read ends never triggered a split).
        let strict_lend = std::env::var_os("RUSTLE_LONGTRIM_STRICT_LEND").is_some();
        // RUSTLE_LONGTRIM_STRICT_LEND: plumb real read-end boundaries through
        // (rather than coverage-derivative only). Default tuning chosen to
        // minimize regression vs baseline — sweep on GGO_19 showed no setting
        // that BOTH reduces last_node_end drift AND preserves the 1655 match
        // baseline. Best at MR=5 MT=20 (-34 matches). See memory for table.
        let min_boundary_reads: f64 = if strict_lend {
            std::env::var("RUSTLE_LONGTRIM_STRICT_MIN_READS")
                .ok()
                .and_then(|v| v.parse::<f64>().ok())
                .unwrap_or(5.0)
        } else {
            3.0
        };
        // Minimum tmpcov (window coverage contrast) to accept a split.
        // StringTie has no threshold (splits on any positive tmpcov); keeping 25
        // here avoids false splits at ambiguous boundaries. Override via
        // RUSTLE_LONGTRIM_MIN_TMPCOV; strict_lend forces 0.0.
        let min_tmpcov: f64 = if strict_lend {
            std::env::var("RUSTLE_LONGTRIM_STRICT_MIN_TMPCOV")
                .ok()
                .and_then(|v| v.parse::<f64>().ok())
                .unwrap_or(20.0)
        } else {
            std::env::var("RUSTLE_LONGTRIM_MIN_TMPCOV")
                .ok()
                .and_then(|v| v.parse::<f64>().ok())
                .unwrap_or(25.0)
        };

        if use_start {
            let pos = lstart[*nls].pos;
            let boundary_cov = lstart[*nls].cov;
            let cur_start = graph.nodes[*graphnode_id].start;

            // Skip boundaries with insufficient read support.
            if boundary_cov.abs() < min_boundary_reads {
                lt_emit(graph, *graphnode_id, "source", pos, boundary_cov, 0.0, 0, "low_reads");
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
            // CPAS-style force: when a boundary cluster has strong support
            // (>= cpas_force_cov reads, default 15) it represents a real
            // alt-TSS even if the coverage transition isn't sharp. Mirrors
            // ST's CPAS path in find_all_trims (rlink.cpp:2712) where the
            // tstartend match sets lastdrop=0 to force trimpoint emission
            // regardless of coverage drop. Off by default — opt in via
            // RUSTLE_CPAS_FORCE_COV=N (set 0 to disable; recommended 15-20).
            let cpas_force_cov: f64 = std::env::var("RUSTLE_CPAS_FORCE_COV")
                .ok().and_then(|v| v.parse().ok()).unwrap_or(0.0);
            if cpas_force_cov > 0.0 && boundary_cov >= cpas_force_cov && tmpcov < min_tmpcov + 1.0 {
                tmpcov = min_tmpcov + 1.0;
            }
            // ST find_all_trims-style localdrop ratio gate (rlink.cpp:2718).
            // ST's source-trim accepts lstart splits only when
            //   covleft_avg / covright_avg < localdrop (default 0.2 = 5× drop).
            // Rustle's longtrim_inline previously used only the absolute
            // tmpcov gate, which fires on flat-coverage internal positions
            // where avg coverage is high. Result: 6× hardstart over-marking
            // vs ST. Default ON (opt-out via RUSTLE_LONGTRIM_RATIO_GATE_OFF=1).
            let ratio_gate_on =
                std::env::var_os("RUSTLE_LONGTRIM_RATIO_GATE").is_some();
            let ratio_threshold: f64 = std::env::var("RUSTLE_LONGTRIM_RATIO")
                .ok().and_then(|v| v.parse().ok()).unwrap_or(0.2);
            let mut ratio_pass = true;
            if ratio_gate_on && tmpcov > min_tmpcov && start_ok && end_ok && pos > cur_start {
                let startpos = (pos - bundle_start) as i64;
                let winstart = (startpos - CHI_THR).max(0);
                let winend = (startpos + CHI_THR - 1).min(bpcov_len as i64 - 1);
                let left_window = (startpos - winstart).max(1) as f64;
                let right_window = (winend - startpos + 1).max(1) as f64;
                let avg_left = get_cov(winstart, startpos - 1) / left_window;
                let avg_right = get_cov(startpos, winend) / right_window;
                if avg_right > 0.0 && avg_left / avg_right >= ratio_threshold {
                    ratio_pass = false;
                }
            }
            // CPAS-force boundaries bypass the ratio gate (matches ST's
            // lastdrop=0 force at find_all_trims:2712).
            if cpas_force_cov > 0.0 && boundary_cov >= cpas_force_cov {
                ratio_pass = true;
            }
            if tmpcov > min_tmpcov && ratio_pass {
                let cur_end = graph.nodes[*graphnode_id].end;
                if pos > graph.nodes[*graphnode_id].start && pos < cur_end {
                    lt_emit(graph, *graphnode_id, "source", pos, boundary_cov, tmpcov, 1, "tmpcov_pos");
                    // Split: [cur_start, pos) and [pos, cur_end)
                    let prev_id = *graphnode_id;
                    graph.nodes[prev_id].end = pos;
                    let new_node = graph.add_node(pos, cur_end);
                    new_node.source_bnode = Some(source_bid);
                    crate::bump_hs!("graph_build.rs:2275:hardstart");
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
                        crate::bump_hs!("graph_build.rs:2288:hardend");
                        graph.nodes[prev_id].hardend = true;
                        sink_parents.push(prev_id);
                        // Record synthetic [prev_id, sink] connector for
                        // post-graph transfrag materialization. Mirrors
                        // StringTie longtrim() futuretr with n2=-1. Gated by
                        // RUSTLE_SINK_CONNECTOR_TF (opt-in) — StringTie parity
                        // candidate; validation pending.
                        if std::env::var_os("RUSTLE_SINK_CONNECTOR_TF").is_some() {
                            sink_futuretr.push((prev_id, tmpcov));
                        }
                    }
                    *graphnode_id = new_id;
                } else {
                    lt_emit(graph, *graphnode_id, "source", pos, boundary_cov, tmpcov, 0, "out_of_range");
                }
            } else {
                let reason = if !ratio_pass { "ratio_fail" }
                             else if !(start_ok && end_ok) { "out_of_anchor" }
                             else { "tmpcov_le_min" };
                lt_emit(graph, *graphnode_id, "source", pos, boundary_cov, tmpcov, 0, reason);
            }
            *nls += 1;
        } else {
            let pos = lend[*nle].pos;
            let boundary_cov = lend[*nle].cov;
            let cur_start = graph.nodes[*graphnode_id].start;

            if boundary_cov.abs() < min_boundary_reads {
                lt_emit(graph, *graphnode_id, "sink", pos, boundary_cov, 0.0, 0, "low_reads");
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
            // CPAS-style force (mirror of lstart): force a split at
            // strongly-supported alt-TTS boundaries even when surrounding
            // coverage is flat. See lstart branch for rationale.
            let cpas_force_cov: f64 = std::env::var("RUSTLE_CPAS_FORCE_COV")
                .ok().and_then(|v| v.parse().ok()).unwrap_or(0.0);
            let boundary_end_cov = lend[*nle].cov;
            if cpas_force_cov > 0.0 && boundary_end_cov >= cpas_force_cov && tmpcov < min_tmpcov + 1.0 {
                tmpcov = min_tmpcov + 1.0;
            }
            // ST localdrop ratio gate (mirror of lstart): for sink-trim,
            // covright_avg / covleft_avg < localdrop (right is 5× lower).
            // Bypassed at strong CPAS clusters.
            let ratio_gate_on =
                std::env::var_os("RUSTLE_LONGTRIM_RATIO_GATE").is_some();
            let ratio_threshold: f64 = std::env::var("RUSTLE_LONGTRIM_RATIO")
                .ok().and_then(|v| v.parse().ok()).unwrap_or(0.2);
            let mut ratio_pass = true;
            if ratio_gate_on && tmpcov > min_tmpcov && start_ok && end_ok && pos > cur_start {
                let endpos = (pos - bundle_start) as i64;
                let winstart = (endpos - CHI_THR + 1).max(0);
                let winend = (endpos + CHI_THR).min(bpcov_len as i64 - 1);
                let left_window = (endpos - winstart + 1).max(1) as f64;
                let right_window = (winend - endpos).max(1) as f64;
                let avg_left = get_cov(winstart, endpos) / left_window;
                let avg_right = get_cov(endpos + 1, winend) / right_window;
                if avg_left > 0.0 && avg_right / avg_left >= ratio_threshold {
                    ratio_pass = false;
                }
            }
            if cpas_force_cov > 0.0 && boundary_end_cov >= cpas_force_cov {
                ratio_pass = true;
            }
            if tmpcov > min_tmpcov && ratio_pass {
                let split = pos + 1; // half-open boundary
                let cur_end = graph.nodes[*graphnode_id].end;
                if split > graph.nodes[*graphnode_id].start && split < cur_end {
                    lt_emit(graph, *graphnode_id, "sink", pos, boundary_cov, tmpcov, 1, "tmpcov_pos");
                    // Split: [cur_start, split) and [split, cur_end)
                    graph.nodes[*graphnode_id].end = split;
                    crate::bump_hs!("graph_build.rs:2365:hardend");
                    graph.nodes[*graphnode_id].hardend = true;
                    let new_node = graph.add_node(split, cur_end);
                    new_node.source_bnode = Some(source_bid);
                    let new_id = new_node.node_id;
                    // Prev → new (contiguous)
                    graph.add_edge(*graphnode_id, new_id);
                    // Prev → sink (hardend)
                    let prev_for_sink = *graphnode_id;
                    sink_parents.push(prev_for_sink);
                    // Record synthetic [prev, sink] connector for lend event.
                    // Same rationale as lstart branch above.
                    if std::env::var_os("RUSTLE_SINK_CONNECTOR_TF").is_some() {
                        sink_futuretr.push((prev_for_sink, tmpcov));
                    }
                    *graphnode_id = new_id;
                } else {
                    lt_emit(graph, *graphnode_id, "sink", pos, boundary_cov, tmpcov, 0, "out_of_range");
                }
            } else {
                let reason = if !ratio_pass { "ratio_fail" }
                             else if !(start_ok && end_ok) { "out_of_anchor" }
                             else { "tmpcov_le_min" };
                lt_emit(graph, *graphnode_id, "sink", pos, boundary_cov, tmpcov, 0, reason);
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
    bundle_chrom: &str,
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
        // Default: pass empty externals so per-bundlenode coverage-derivative
        // detection runs inside create_graph_inner.
        //
        // RUSTLE_LONGTRIM_STRICT_LEND=1 overrides: feed real read-end boundaries
        // directly so `longtrim_inline` can split nodes at actual read 3' ends.
        // This is the ST-faithful behavior — StringTie's longtrim takes raw
        // lstart/lend per-read and splits on coverage contrast alone. Rustle's
        // default (coverage-derivative-only) misses 5096 cases where the last
        // graph node extends 20-100bp past the read's actual end.
        if std::env::var_os("RUSTLE_LONGTRIM_STRICT_LEND").is_some() {
            (lstart, lend)
        } else {
            (&[] as &[ReadBoundary], &[] as &[ReadBoundary])
        }
    } else {
        (&[] as &[ReadBoundary], &[] as &[ReadBoundary])
    };
    let mut sink_futuretr: Vec<(usize, f64)> = Vec::new();
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
        &mut sink_futuretr,
        bundle_chrom,
    );

    // Zero-width node compaction (DEFAULT-ON): junction events at shared
    // donor/acceptor coords leave `[X, X)` empty nodes. Treat them as transparent
    // via transitive closure — for every path u→[zw*]→v, add direct edge u→v.
    // Reduces graph node count ~1300 / edges ~4500 globally on GGO_19. Net +12
    // fewer novel chains, -3 valid chains (20% FP rate, far below the 40-60% FP
    // rate of generic filters).
    // Opt-out: RUSTLE_COMPACT_ZW_NODES_OFF=1.
    //
    // Extension (opt-in): RUSTLE_COMPACT_PASS_THROUGH=1 also merges pass-through
    // nodes (exactly one contiguous parent + one contiguous child, non-source/sink)
    // via the same transitive-closure approach. Addresses remaining non-zw
    // over-segmentation at ~6,500 nodes globally.
    if std::env::var_os("RUSTLE_COMPACT_ZW_NODES_OFF").is_none() {
        let include_pt = std::env::var_os("RUSTLE_COMPACT_PASS_THROUGH").is_some();
        let removed = graph.compact_transparent_nodes(include_pt);
        if std::env::var_os("RUSTLE_COMPACT_ZW_NODES_TRACE").is_some() && removed > 0 {
            eprintln!(
                "COMPACT_ZW bundle={}-{} strand={} removed={} remaining_nodes={} pass_through={}",
                bundle_start, bundle_end, bundle_strand, removed, graph.n_nodes, include_pt
            );
        }
    }

    // Proper chain-compression pass (opt-in via RUSTLE_COMPACT_CHAIN=1).
    // Extends head-node coord to absorb downstream pass-through members.
    //
    // EXPERIMENTAL — DO NOT ENABLE. Measured on GGO_19: removes 1219 nodes
    // across 351 bundles but regresses catastrophically:
    //   matching tx 1637 → 729 (-908)
    //   transcript F1  87 → 36
    //   base sens 94.3% → 82.5%
    // Filter excludes source/sink/hardstart/hardend/longtrim_cov>0 nodes,
    // but the surviving pass-throughs (exon-interior segments created by
    // junction events elsewhere in the locus) still carry signal that the
    // downstream pipeline (pathpat compatibility, seed selection, flow
    // path enumeration) depends on. Merging them yields a more permissive
    // graph with many more emitted paths, most novel/spurious.
    //
    // Kept for future investigation: the right approach would require
    // simultaneous recalibration of downstream pathpat heuristics.
    if std::env::var_os("RUSTLE_COMPACT_CHAIN").is_some() {
        let removed = graph.compact_chain_pass_through();
        if std::env::var_os("RUSTLE_COMPACT_CHAIN_TRACE").is_some() && removed > 0 {
            eprintln!(
                "COMPACT_CHAIN bundle={}-{} strand={} removed={} remaining_nodes={}",
                bundle_start, bundle_end, bundle_strand, removed, graph.n_nodes
            );
        }
    }

    let mut stats = CreateGraphLongtrimStats {
        lstart_events: lstart.len(),
        lend_events: lend.len(),
        ..Default::default()
    };
    let mut synthetic: Vec<GraphTransfrag> = Vec::new();

    // NOTE: sink-connector synthesis moved to pipeline.rs, just before
    // process_transfrags runs. Graph node IDs can be remapped by post-
    // create_graph steps (prune/split/redirect), which would invalidate
    // [n, sink] connectors created here. See `synthesize_sink_connectors` in
    // pipeline.rs for the active implementation. sink_futuretr still collected
    // here for future reference (e.g. to seed abundance) but not materialized.
    let _ = (&sink_futuretr, bundle_strand, bundle_start, bundle_end);

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
        // RUSTLE_APPLY_LONGTRIM_DIRECT_OFF=1: skip the bpcov-derivative
        // longtrim split pass that runs on top of longtrim_inline. This
        // pass is the dominant hardstart/hardend over-marking source —
        // 1001 + 158 + 19 hardstart writes per default GGO_19 run come
        // from here regardless of whether external lstart/lend are
        // supplied (the boundary_map's own bpcov-diff scan generates
        // candidates internally).
        let skip_apply =
            std::env::var_os("RUSTLE_APPLY_LONGTRIM_DIRECT_OFF").is_some();
        let (lt_synth, lt_stats) = if skip_apply {
            (
                Vec::<GraphTransfrag>::new(),
                LongtrimStats::default(),
            )
        } else {
            apply_longtrim_direct(
                &mut graph,
                &boundary_map,
                &schedule,
                bpcov,
                longtrim_min_boundary_cov,
                trace_strand_index(bundle_strand),
                bundle_start,
                bundle_end,
            )
        };
        stats.applied = !skip_apply;
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
