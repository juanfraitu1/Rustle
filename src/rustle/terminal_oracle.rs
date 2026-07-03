//! DEAD CODE -- StringTie-era assembler; slated for removal.
//! NOT part of the multi-copy-family thesis (O1 family-def / O2 copy-assign / O3 ASJ / O4 absent-copy).
//! See docs/RETIREMENT_AND_MIGRATION.md. Do not extend.
//! Injection oracle (analysis-only, default OFF via RUSTLE_TERMINAL_ORACLE=<path>).
//! On identical-node bundles, override Rustle's synthetic source/sink edge set with
//! StringTie's captured `graph_edge` set, reconciling the synthetic transfrag abundances,
//! to bound the FP/TP impact of the terminal-wiring divergence before any production change.
//!
//! CONFOUND: flow's build_lr_edge_capacities creates virtual source/sink capacity edges
//! for partial real transfrags, so this measures the synthetic-coverage-edge contribution
//! — a LOWER BOUND on full ST-terminal suppression.

use std::collections::{HashMap, HashSet};
use std::sync::OnceLock;

use crate::bpcov::BpcovStranded;
use crate::graph::{Graph, GraphTransfrag};

/// Coverage-drop divisor used when reconciling synthetic terminal-edge abundance.
/// Mirrors `graph_build::DROP` (private there); kept in sync by value.
const DROP: f64 = 0.5;

/// Per-bundle ST terminal edges: bundle key (start1, end, strand) -> set of (from_token, to_token).
/// Tokens: "A-B" for a real node (1-based inclusive), "SRC"/"SNK" for source/sink.
pub type OracleMap = HashMap<(u64, u64, char), HashSet<(String, String)>>;

static ORACLE: OnceLock<OracleMap> = OnceLock::new();

/// Load ST graph_edge events once from the oracle file path.
pub fn load_oracle(path: &str) -> &'static OracleMap {
    ORACLE.get_or_init(|| {
        let mut map: OracleMap = HashMap::new();
        let Ok(content) = std::fs::read_to_string(path) else {
            return map;
        };
        for line in content.lines() {
            let Ok(v) = serde_json::from_str::<serde_json::Value>(line) else {
                continue;
            };
            if v.get("step").and_then(|s| s.as_str()) != Some("graph_edge") {
                continue;
            }
            let start = v.get("start").and_then(|x| x.as_u64()).unwrap_or(0);
            let end = v.get("end").and_then(|x| x.as_u64()).unwrap_or(0);
            let strand = v
                .get("strand")
                .and_then(|x| x.as_str())
                .and_then(|s| s.chars().next())
                .unwrap_or('.');
            let edges_str = v
                .get("payload")
                .and_then(|p| p.get("edges"))
                .and_then(|e| e.as_str())
                .unwrap_or("");
            let set = map.entry((start, end, strand)).or_default();
            for e in edges_str.split(',') {
                if let Some((f, t)) = e.split_once('>') {
                    set.insert((f.to_string(), t.to_string()));
                }
            }
        }
        map
    })
}

/// Node-coordinate token, matching the graph_edge emit convention (1-based inclusive).
fn node_token(g: &Graph, nid: usize) -> String {
    if nid == g.source_id {
        return "SRC".to_string();
    }
    if nid == g.sink_id {
        return "SNK".to_string();
    }
    let n = &g.nodes[nid];
    format!("{}-{}", n.start + 1, n.end)
}

/// Build a token -> node_id map for real nodes (for mapping ST tokens back to RU ids).
fn token_to_id(g: &Graph) -> HashMap<String, usize> {
    let mut m = HashMap::new();
    for i in 0..g.n_nodes {
        if i == g.source_id || i == g.sink_id {
            continue;
        }
        m.insert(node_token(g, i), i);
    }
    m
}

/// Average coverage of a node (icov), mirroring add_coverage_source_sink_edges.
fn node_avg_cov(g: &Graph, nid: usize, bpcov: &BpcovStranded, sno: usize) -> f64 {
    let n = &g.nodes[nid];
    if n.end <= n.start {
        return 0.0;
    }
    let len = (n.end - n.start) as f64;
    let si = bpcov.plus.idx(n.start);
    let ei = bpcov.plus.idx(n.end);
    bpcov.get_cov_range(sno, si, ei) / len
}

/// Override RU's source/sink edges on an identical-node bundle to match ST's captured set.
/// Returns (edges_removed, edges_installed); mutates `graph` and filters/extends `coverage_synth`.
/// Returns None (no-op) if ST's referenced nodes are not a subset of RU's (cannot inject cleanly).
pub fn override_terminal_edges(
    graph: &mut Graph,
    bundle_key: (u64, u64, char),
    coverage_synth: &mut Vec<GraphTransfrag>,
    bpcov: &BpcovStranded,
    sno: usize,
    oracle: &OracleMap,
) -> Option<(usize, usize)> {
    let st_edges = oracle.get(&bundle_key)?;
    let tok2id = token_to_id(graph);

    // Identical-node check: ST's referenced real nodes (tokens appearing in ST's edges,
    // minus SRC/SNK) must be a subset of RU's real-node token set. Bundles failing this
    // are skipped (can't inject cleanly).
    let ru_nodes: HashSet<String> = tok2id.keys().cloned().collect();
    let st_nodes: HashSet<String> = st_edges
        .iter()
        .flat_map(|(f, t)| [f.clone(), t.clone()])
        .filter(|x| x != "SRC" && x != "SNK")
        .collect();
    if !st_nodes.is_subset(&ru_nodes) {
        return None;
    }

    let src = graph.source_id;
    let snk = graph.sink_id;

    // Current RU source/sink edges as node ids.
    let mut ru_src_edges: Vec<usize> = graph.nodes[src].children.ones().collect();
    let mut ru_snk_edges: Vec<usize> = graph.nodes[snk].parents.ones().collect();

    // ST's desired source/sink edges, mapped to RU ids.
    let st_src_targets: HashSet<usize> = st_edges
        .iter()
        .filter(|(f, _)| f == "SRC")
        .filter_map(|(_, t)| tok2id.get(t).copied())
        .collect();
    let st_snk_sources: HashSet<usize> = st_edges
        .iter()
        .filter(|(_, t)| t == "SNK")
        .filter_map(|(f, _)| tok2id.get(f).copied())
        .collect();

    let mut removed = 0usize;
    let mut installed = 0usize;

    // 1) REMOVE RU source edges ST does not have.
    for &nid in &ru_src_edges {
        if !st_src_targets.contains(&nid) {
            graph.remove_edge(src, nid);
            removed += 1;
        }
    }
    // 2) REMOVE RU sink edges ST does not have.
    for &nid in &ru_snk_edges {
        if !st_snk_sources.contains(&nid) {
            graph.remove_edge(nid, snk);
            removed += 1;
        }
    }
    // Recompute after removals for the "install missing" pass.
    ru_src_edges = graph.nodes[src].children.ones().collect();
    ru_snk_edges = graph.nodes[snk].parents.ones().collect();
    let ru_src_set: HashSet<usize> = ru_src_edges.iter().copied().collect();
    let ru_snk_set: HashSet<usize> = ru_snk_edges.iter().copied().collect();

    let pattern_size = graph.pattern_size();
    // 3) INSTALL ST source edges RU lacks (+ synth transfrag with reconciled abundance).
    for &nid in &st_src_targets {
        if !ru_src_set.contains(&nid) {
            let icov = node_avg_cov(graph, nid, bpcov, sno);
            let parcov: f64 = graph.nodes[nid]
                .parents
                .ones()
                .filter(|&p| p != src && p != snk)
                .map(|p| node_avg_cov(graph, p, bpcov, sno))
                .sum();
            graph.add_edge(src, nid);
            let mut tf = GraphTransfrag::new(vec![src, nid], pattern_size);
            tf.abundance = ((icov - parcov) / DROP).max(0.0);
            coverage_synth.push(tf);
            installed += 1;
        }
    }
    // 4) INSTALL ST sink edges RU lacks (+ synth transfrag with reconciled abundance).
    for &nid in &st_snk_sources {
        if !ru_snk_set.contains(&nid) {
            let icov = node_avg_cov(graph, nid, bpcov, sno);
            let chcov: f64 = graph.nodes[nid]
                .children
                .ones()
                .filter(|&c| c != src && c != snk)
                .map(|c| node_avg_cov(graph, c, bpcov, sno))
                .sum();
            graph.add_edge(nid, snk);
            let mut tf = GraphTransfrag::new(vec![nid, snk], pattern_size);
            tf.abundance = ((icov - chcov) / DROP).max(0.0);
            coverage_synth.push(tf);
            installed += 1;
        }
    }

    // 5) DROP synth transfrags whose 2-node source/sink edge was removed.
    coverage_synth.retain(|tf| {
        if tf.node_ids.len() != 2 {
            return true;
        }
        let (a, b) = (tf.node_ids[0], tf.node_ids[1]);
        if a == src {
            return graph.nodes[src].children.contains(b);
        }
        if b == snk {
            return graph.nodes[snk].parents.contains(a);
        }
        true
    });

    Some((removed, installed))
}
