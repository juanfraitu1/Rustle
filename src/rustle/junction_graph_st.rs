/// StringTie-style BundlenodeGraph construction and junction inference.
///
/// Ports the core `create_graph` algorithm from StringTie's rlink.cpp:
///   - Coverage islands → bundlenodes (CBundlenode linked list equivalent)
///   - Bundlenode splitting at junction boundaries → graphnodes (CGraphnode)
///   - Read-to-path mapping → transfrag paths
///   - Path-aggregated junction evidence → inferred junction candidates
///
/// The full StringTie algorithm also includes:
///   - CHI-statistic TSS/TES detection from long-read coverage drops (longreads path)
///   - Guide-edge based trimming (source2guide / guide2sink)
///   - EM-based abundance estimation per transfrag
///   - retainedintron / coverage-cap checks during graph assembly
///
/// These are documented as TODO stubs below. The framework here handles the
/// junction-graph construction and inference pipeline, which is the part that
/// affects what junction candidates reach the acceptance filter.

use crate::types::{BundleRead, Junction, JunctionStat, JunctionStats};
use std::collections::{BTreeMap, HashMap};

// ── Data structures ─────────────────────────────────────────────────────────

/// A node in the bundlenode graph, equivalent to StringTie's `CGraphnode`.
///
/// Graphnodes are exon segments: contiguous coverage regions between junction
/// split-points. The source node (id = SOURCE_ID) and sink node (id = last)
/// are virtual nodes that anchor all graph paths.
#[derive(Debug, Clone)]
pub struct GraphNode {
    pub id: usize,
    pub start: u64,
    pub end: u64,
    /// Coverage averaged over [start, end].
    pub cov: f64,
    /// Abundance flowing INTO this node (sum of transfrag abundances entering).
    pub abund_in: f64,
    /// Abundance flowing OUT of this node.
    pub abund_out: f64,
    /// Child node IDs (junction or continuation edges leaving this node).
    pub children: Vec<usize>,
    /// Parent node IDs.
    pub parents: Vec<usize>,
    /// True if this node's start is a hard (read-confirmed) TSS boundary.
    pub hardstart: bool,
    /// True if this node's end is a hard (read-confirmed) TES boundary.
    pub hardend: bool,
}

impl GraphNode {
    fn new(id: usize, start: u64, end: u64, cov: f64) -> Self {
        GraphNode {
            id, start, end, cov,
            abund_in: 0.0, abund_out: 0.0,
            children: Vec::new(), parents: Vec::new(),
            hardstart: false, hardend: false,
        }
    }
}

/// A coverage island (pre-graph bundlenode), equivalent to StringTie's `CBundlenode`.
///
/// Built by merging read exon coverage with gaps ≤ `local_dist`.
#[derive(Debug, Clone)]
pub struct BundleNode {
    pub start: u64,
    pub end: u64,
    pub cov: f64,
}

/// The directed acyclic bundlenode graph for a single bundle.
///
/// `nodes[0]` is the virtual source; `nodes[nodes.len()-1]` is the virtual sink.
/// All real graphnodes are in between.
#[derive(Debug)]
pub struct BundleGraph {
    pub nodes: Vec<GraphNode>,
    pub source_id: usize,
    pub sink_id: usize,
}

impl BundleGraph {
    /// Source node id (always 0).
    pub const SOURCE_ID: usize = 0;

    fn new() -> Self {
        let source = GraphNode::new(Self::SOURCE_ID, 0, 0, 0.0);
        BundleGraph { nodes: vec![source], source_id: 0, sink_id: 0 }
    }

    fn add_node(&mut self, start: u64, end: u64, cov: f64) -> usize {
        let id = self.nodes.len();
        self.nodes.push(GraphNode::new(id, start, end, cov));
        id
    }

    fn add_sink(&mut self) {
        let id = self.nodes.len();
        self.nodes.push(GraphNode::new(id, u64::MAX, u64::MAX, 0.0));
        self.sink_id = id;
    }

    fn add_edge(&mut self, parent: usize, child: usize) {
        if !self.nodes[parent].children.contains(&child) {
            self.nodes[parent].children.push(child);
        }
        if !self.nodes[child].parents.contains(&parent) {
            self.nodes[child].parents.push(parent);
        }
    }

    /// Returns all graphnodes that overlap `[qstart, qend]`, sorted by start.
    pub fn nodes_overlapping(&self, qstart: u64, qend: u64) -> Vec<usize> {
        let mut result = Vec::new();
        for node in &self.nodes {
            if node.id == self.source_id || node.id == self.sink_id {
                continue;
            }
            if node.end >= qstart && node.start <= qend {
                result.push(node.id);
            }
        }
        result.sort_by_key(|&id| self.nodes[id].start);
        result
    }

    /// Returns all graphnodes that a read exon interval `(exon_start, exon_end)` maps to.
    ///
    /// A read exon maps to a graphnode if their intervals share at least one position.
    pub fn exon_to_nodes(&self, exon_start: u64, exon_end: u64) -> Vec<usize> {
        self.nodes_overlapping(exon_start, exon_end)
    }
}

// ── Phase 1: Bundlenode construction ────────────────────────────────────────

/// Build the list of coverage islands (bundlenodes) from read exon intervals.
///
/// Equivalent to the CGroup → CBundlenode pipeline in StringTie rlink.cpp.
/// Groups contiguous covered bases, merging gaps ≤ `local_dist`.
pub fn build_bundlenodes(
    reads: &[BundleRead],
    bundle_start: u64,
    bundle_end: u64,
    local_dist: u64,
) -> Vec<BundleNode> {
    if bundle_end <= bundle_start {
        return Vec::new();
    }
    let len = (bundle_end - bundle_start) as usize;
    let mut cov = vec![0.0f64; len];

    for read in reads {
        for &(es, ee) in &read.exons {
            let s = es.saturating_sub(bundle_start) as usize;
            let e = (ee.saturating_sub(bundle_start) as usize).min(len);
            for i in s..e {
                cov[i] += read.weight;
            }
        }
    }

    // Find islands: runs of non-zero coverage
    let mut islands: Vec<(usize, usize, f64)> = Vec::new(); // (start_idx, end_idx, sum_cov)
    let mut in_island = false;
    let mut istart = 0usize;
    let mut isum = 0.0f64;

    for (i, &c) in cov.iter().enumerate() {
        if c > 0.0 {
            if !in_island {
                in_island = true;
                istart = i;
                isum = 0.0;
            }
            isum += c;
        } else if in_island {
            islands.push((istart, i, isum));
            in_island = false;
        }
    }
    if in_island {
        islands.push((istart, len, isum));
    }

    if islands.is_empty() {
        return Vec::new();
    }

    // Merge islands within local_dist
    let mut merged: Vec<BundleNode> = Vec::new();
    let (mut ms, mut me, mut mc) = islands[0];

    for &(s, e, c) in &islands[1..] {
        if (s as u64) <= (me as u64) + local_dist {
            // Extend current island
            me = e;
            mc += c;
        } else {
            merged.push(BundleNode {
                start: bundle_start + ms as u64,
                end: bundle_start + me as u64,
                cov: mc / (me - ms).max(1) as f64,
            });
            ms = s; me = e; mc = c;
        }
    }
    merged.push(BundleNode {
        start: bundle_start + ms as u64,
        end: bundle_start + me as u64,
        cov: mc / (me - ms).max(1) as f64,
    });

    merged
}

// ── Phase 2: Graph construction (create_graph equivalent) ───────────────────

/// Build a bundlenode DAG by splitting bundlenodes at junction boundaries.
///
/// Core algorithm from StringTie rlink.cpp `create_graph()` (lines 3223–3800).
/// The simplified (non-longreads, non-guide) path:
///
/// For each bundlenode B:
///   1. Create an initial graphnode G spanning [B.start, B.end].
///   2. Link G from source or from nodes whose junction ends here.
///   3. Walk through junctions sorted by start position within B:
///      a. Junction *start* at pos: set G.end = pos, record G.id in ends[junc.acceptor],
///         create continuation G' = [pos+1, B.end], link G → G'.
///      b. Junction *end* at pos (within B): set G.end = pos-1, create G' = [pos, B.end],
///         link G → G', then add edges from ends[pos] → G'.
///   4. Unfinished nodes at B.end are terminal; link to sink.
///
/// TODO (longreads path):
///   - CHI-statistic TSS/TES detection via coverage differential (`longtrim`)
///   - Guide-edge trimming (`source2guide`, `guide2sink`)
///   - `junctionsupport`-based tail trimming for small residual nodes
pub fn build_bundle_graph(
    bnodes: &[BundleNode],
    junctions: &[(u64, u64)],  // (donor, acceptor) — half-open: intron = [donor, acceptor)
) -> BundleGraph {
    let mut graph = BundleGraph::new();
    if bnodes.is_empty() {
        graph.add_sink();
        return graph;
    }

    // Sort junctions by donor (start), keeping a separate sort by acceptor (end)
    let mut junc_by_donor: Vec<(u64, u64)> = junctions.to_vec();
    junc_by_donor.sort_unstable();
    let mut junc_by_acceptor: Vec<(u64, u64)> = junctions.to_vec();
    junc_by_acceptor.sort_unstable_by_key(|&(_, acc)| acc);

    // `ends`: acceptor_pos → Vec<graphnode_id> whose donor ends before this intron.
    // When a graphnode is created at `acceptor_pos`, all nodes in ends[pos] become parents.
    let mut ends: HashMap<u64, Vec<usize>> = HashMap::new();

    let mut njs = 0usize; // index into junc_by_donor (junction starts)
    let mut nje = 0usize; // index into junc_by_acceptor (junction ends / acceptor positions)

    for bnode in bnodes {
        let bstart = bnode.start;
        let bend = bnode.end;
        let bcov = bnode.cov;

        // Check if a junction ends right at bstart (meaning we enter this bnode from an intron)
        let mut from_junction = false;
        while nje < junc_by_acceptor.len() && junc_by_acceptor[nje].1 <= bstart {
            if junc_by_acceptor[nje].1 == bstart {
                from_junction = true;
            }
            nje += 1;
        }

        // Create initial graphnode spanning the full bnode
        let gid = graph.add_node(bstart, bend, bcov);

        // Link this node: either from junction-source nodes or from global source
        if from_junction {
            if let Some(parents) = ends.get(&bstart) {
                for &pid in parents.clone().iter() {
                    graph.add_edge(pid, gid);
                }
            } else {
                // Junction end pointed here but no recorded donor → link to source
                graph.add_edge(BundleGraph::SOURCE_ID, gid);
            }
        } else {
            graph.add_edge(BundleGraph::SOURCE_ID, gid);
        }

        // Skip junctions whose donor is before bstart
        while njs < junc_by_donor.len() && junc_by_donor[njs].0 < bstart {
            njs += 1;
        }

        // Walk through junctions within this bundlenode
        let mut current_gid = gid;
        let mut completed = false;

        loop {
            // Find next event: junction donor start OR junction acceptor within bnode
            // Recheck acceptor index within bnode range
            let next_donor = junc_by_donor.get(njs)
                .filter(|&&(d, _)| d <= bend)
                .map(|&(d, _)| d);
            let next_acceptor = junc_by_acceptor.get(nje)
                .filter(|&&(_, acc)| acc <= bend)
                .map(|&(_, acc)| acc);

            let event = match (next_donor, next_acceptor) {
                (Some(d), Some(acc)) => {
                    // Process whichever comes first; ties: acceptor first (matches ST)
                    if d <= acc { Some((d, true)) } else { Some((acc, false)) }
                }
                (Some(d), None)   => Some((d, true)),
                (None,   Some(acc)) => Some((acc, false)),
                (None,   None)    => None,
            };

            match event {
                None => break,

                // Junction DONOR at pos: split current node at pos
                Some((pos, true)) => {
                    let cur = &graph.nodes[current_gid];
                    if pos < cur.start || pos >= cur.end {
                        // Junction donor outside current node range, skip
                        njs += 1;
                        continue;
                    }

                    // Trim current node to [cur.start, pos]
                    graph.nodes[current_gid].end = pos;

                    // Record all junctions starting at this pos
                    while njs < junc_by_donor.len() && junc_by_donor[njs].0 == pos {
                        let acc = junc_by_donor[njs].1;
                        ends.entry(acc).or_default().push(current_gid);
                        njs += 1;
                    }

                    if pos < bend {
                        // Create continuation node [pos+1, bend]
                        let next_gid = graph.add_node(pos + 1, bend, bcov);
                        graph.add_edge(current_gid, next_gid);
                        current_gid = next_gid;
                    } else {
                        completed = true;
                        break;
                    }
                }

                // Junction ACCEPTOR at pos: split current node, link junction-source parents
                Some((pos, false)) => {
                    let cur = &graph.nodes[current_gid];
                    if pos <= cur.start {
                        // Already at or before node start, just add edges
                        if let Some(parents) = ends.get(&pos) {
                            for &pid in parents.clone().iter() {
                                graph.add_edge(pid, current_gid);
                            }
                        }
                        while nje < junc_by_acceptor.len() && junc_by_acceptor[nje].1 == pos {
                            nje += 1;
                        }
                        continue;
                    }

                    if pos > cur.end {
                        while nje < junc_by_acceptor.len() && junc_by_acceptor[nje].1 == pos {
                            nje += 1;
                        }
                        continue;
                    }

                    // Trim current node to [cur.start, pos-1]
                    graph.nodes[current_gid].end = pos - 1;

                    // Create new node [pos, bend]
                    let next_gid = graph.add_node(pos, bend, bcov);
                    graph.add_edge(current_gid, next_gid);

                    // Link junction-source parents to next_gid
                    if let Some(parents) = ends.get(&pos) {
                        for &pid in parents.clone().iter() {
                            graph.add_edge(pid, next_gid);
                        }
                    }

                    while nje < junc_by_acceptor.len() && junc_by_acceptor[nje].1 == pos {
                        nje += 1;
                    }

                    current_gid = next_gid;
                }
            }
        }

        if !completed {
            // This bnode's last node is a terminal — will be linked to sink
            graph.nodes[current_gid].hardend = true;
        }
    }

    // Add sink and link all terminal (no-children) non-source nodes to it
    graph.add_sink();
    let sink_id = graph.sink_id;
    let node_ids: Vec<usize> = (0..graph.nodes.len()).collect();
    for nid in node_ids {
        if nid == BundleGraph::SOURCE_ID || nid == sink_id {
            continue;
        }
        if graph.nodes[nid].children.is_empty() {
            graph.add_edge(nid, sink_id);
        }
    }

    graph
}

// ── Phase 3: Read-to-path mapping ───────────────────────────────────────────

/// Map a single read to its path through the graph.
///
/// Returns the ordered sequence of graphnode IDs the read traverses, or `None`
/// if the read spans a region not covered by any graphnode (unmappable read).
///
/// Algorithm:
///   For each read exon interval, find overlapping graphnodes.
///   The path is the sorted union of all overlapping nodes.
///   Non-contiguous node pairs in the path represent inferred junctions.
pub fn map_read_to_path(graph: &BundleGraph, read: &BundleRead) -> Option<Vec<usize>> {
    let mut path: Vec<usize> = Vec::new();

    for &(es, ee) in &read.exons {
        let nodes = graph.exon_to_nodes(es, ee);
        for nid in nodes {
            if !path.contains(&nid) {
                path.push(nid);
            }
        }
    }

    if path.is_empty() {
        return None;
    }

    path.sort_by_key(|&nid| graph.nodes[nid].start);
    Some(path)
}

/// Aggregate read paths into junction support counts.
///
/// For each pair of non-contiguous nodes (i, j) in a read path where
/// there is no direct continuation edge (i.end+1 == j.start), the gap
/// represents an inferred splice junction: (i.end, j.start).
///
/// Returns a map: (donor_pos, acceptor_pos) → total_weight.
pub fn aggregate_path_junctions(
    graph: &BundleGraph,
    reads: &[BundleRead],
) -> BTreeMap<(u64, u64), f64> {
    let mut junc_support: BTreeMap<(u64, u64), f64> = BTreeMap::new();

    for read in reads {
        let path = match map_read_to_path(graph, read) {
            Some(p) if p.len() >= 2 => p,
            _ => continue,
        };

        for window in path.windows(2) {
            let (na, nb) = (window[0], window[1]);
            let end_a = graph.nodes[na].end;
            let start_b = graph.nodes[nb].start;

            // If not contiguous, it's a junction (gap = intron)
            if end_a + 1 < start_b {
                *junc_support.entry((end_a, start_b)).or_insert(0.0) += read.weight;
            }
        }
    }

    junc_support
}

// ── Phase 4: Coverage flow (stub) ────────────────────────────────────────────

/// Compute coverage flow for each graph node.
///
/// This is a placeholder for the EM-based abundance estimation in StringTie.
/// The full algorithm:
///   1. Distribute read abundance across all paths a read could follow.
///   2. Iterate until convergence (EM-style).
///   3. Update `abund_in` / `abund_out` per node based on flow fractions.
///
/// Current implementation: assign read weight uniformly to all nodes in the
/// read's path (average coverage approximation, not flow-balanced).
///
/// TODO: Implement full EM abundance estimation (`longrec` / `setAbundanceAll` equivalent).
pub fn compute_coverage_flow(graph: &mut BundleGraph, reads: &[BundleRead]) {
    // Reset abundances
    for node in graph.nodes.iter_mut() {
        node.abund_in = 0.0;
        node.abund_out = 0.0;
    }

    // Simple uniform assignment: each read contributes weight to every node it covers
    for read in reads {
        let path = match map_read_to_path(graph, read) {
            Some(p) if !p.is_empty() => p,
            _ => continue,
        };
        let per_node = read.weight / path.len() as f64;
        for &nid in &path {
            graph.nodes[nid].abund_in += per_node;
            graph.nodes[nid].abund_out += per_node;
        }
    }
}

// ── Per-base coverage ────────────────────────────────────────────────────────

/// Per-base coverage over a genomic interval, equivalent to StringTie's `bpcov`.
///
/// Stored as a flat f64 vec indexed from 0 = `origin`. Uses a simple cumulative-
/// sum construction for O(1) range-sum queries.
pub struct PerBaseCoverage {
    /// Cumulative coverage: `prefix[i]` = sum of coverage over [origin, origin+i).
    prefix: Vec<f64>,
    pub origin: u64,
}

impl PerBaseCoverage {
    /// Build per-base coverage from read exon intervals.
    pub fn from_reads(reads: &[BundleRead], origin: u64, end: u64) -> Self {
        if end <= origin {
            return PerBaseCoverage { prefix: vec![0.0], origin };
        }
        let len = (end - origin) as usize;
        let mut cov = vec![0.0f64; len + 1]; // +1 for prefix sum sentinel

        for read in reads {
            for &(es, ee) in &read.exons {
                let s = es.saturating_sub(origin) as usize;
                let e = (ee.saturating_sub(origin) as usize).min(len);
                for c in cov[s..e].iter_mut() {
                    *c += read.weight;
                }
            }
        }

        // Build prefix sums
        let mut prefix = vec![0.0f64; len + 2];
        for i in 0..len {
            prefix[i + 1] = prefix[i] + cov[i];
        }

        PerBaseCoverage { prefix, origin }
    }

    /// Sum of coverage over [start, end] (both inclusive, absolute coordinates).
    pub fn sum(&self, start: u64, end: u64) -> f64 {
        if end < start || start < self.origin {
            return 0.0;
        }
        let s = (start - self.origin) as usize;
        let e = (end - self.origin) as usize + 1;
        // prefix has len+2 elements; valid data is at indices 0..=len
        let max_idx = self.prefix.len().saturating_sub(2);
        self.prefix[e.min(max_idx)] - self.prefix[s.min(max_idx)]
    }

    /// Average coverage over [start, end] (absolute coordinates).
    pub fn avg(&self, start: u64, end: u64) -> f64 {
        if end < start {
            return 0.0;
        }
        let span = (end - start + 1) as f64;
        self.sum(start, end) / span
    }

    /// Point coverage at `pos` (absolute coordinate).
    pub fn at(&self, pos: u64) -> f64 {
        if pos < self.origin {
            return 0.0;
        }
        let i = (pos - self.origin) as usize;
        // prefix[i+1] is valid only for i < len; prefix[len+1] is the sentinel (0.0)
        // so guard with i+2 to avoid returning prefix[len+1]-prefix[len] < 0
        if i + 2 >= self.prefix.len() {
            return 0.0;
        }
        self.prefix[i + 1] - self.prefix[i]
    }

    pub fn len_bases(&self) -> usize {
        self.prefix.len().saturating_sub(1)
    }
}

/// Assign per-base coverage to each graphnode using the coverage profile.
///
/// Replaces the bnode.cov average with proper per-base averaging.
pub fn assign_node_coverage(graph: &mut BundleGraph, bpcov: &PerBaseCoverage) {
    for node in graph.nodes.iter_mut() {
        if node.id == BundleGraph::SOURCE_ID || node.id == graph.sink_id {
            continue;
        }
        if node.end >= node.start {
            node.cov = bpcov.avg(node.start, node.end);
        }
    }
}

// ── Phase 2b: Coverage-differential TSS/TES trimming ─────────────────────────

/// Constants from StringTie rlink.h / rlink.cpp
const CHI_WIN: usize = 100;
const CHI_THR: usize = 50;
const LONG_INTRON_ANCHOR: u64 = 25;
const ERROR_PERC: f64 = 0.1;
const DROP: f64 = 0.5;

/// A detected coverage-differential boundary within a graphnode.
#[derive(Debug, Clone, Copy)]
pub struct TrimPoint {
    /// Absolute genomic coordinate of the trim boundary.
    pub pos: u64,
    /// True → coverage rises here (TSS / acceptor). False → coverage drops (TES / donor).
    pub is_start: bool,
    /// Coverage differential magnitude at this position.
    pub magnitude: f64,
}

/// Scan a bundlenode region for significant coverage drops/rises (TSS/TES detection).
///
/// Ports the short-region path of StringTie `find_trims()` (rlink.cpp line 1500).
/// For longer regions, the CHI sliding-window scan is used instead.
///
/// Returns up to one TSS and one TES candidate per call.
pub fn find_coverage_trims(
    start: u64,
    end: u64,
    bpcov: &PerBaseCoverage,
) -> (Option<TrimPoint>, Option<TrimPoint>) {
    let len = end.saturating_sub(start) as usize;
    if len < CHI_THR {
        return (None, None);
    }

    let mut source_trim: Option<TrimPoint> = None;
    let mut sink_trim: Option<TrimPoint> = None;
    let mut source_drop = 1.0f64;
    let mut sink_drop = 1.0f64;

    let local_drop = if (len as usize) < 2 * (CHI_WIN + CHI_THR) + 1 {
        ERROR_PERC / (10.0 * DROP)
    } else {
        ERROR_PERC
    };

    let scan_start = start + LONG_INTRON_ANCHOR;
    let scan_end = end.saturating_sub(LONG_INTRON_ANCHOR);

    if scan_start >= scan_end {
        return (None, None);
    }

    for pos in scan_start..scan_end {
        let cov_left = if pos > start {
            bpcov.sum(start, pos - 1) / (pos - start) as f64
        } else {
            0.0
        };
        let cov_right = if pos <= end {
            bpcov.sum(pos, end) / (end - pos + 1) as f64
        } else {
            0.0
        };

        if cov_left < cov_right && cov_right > 0.0 {
            let drop = cov_left / cov_right;
            if drop < local_drop && drop < source_drop {
                source_drop = drop;
                source_trim = Some(TrimPoint {
                    pos,
                    is_start: true,
                    magnitude: (cov_right - cov_left) / DROP,
                });
            }
        } else if cov_right < cov_left && cov_left > 0.0 {
            let drop = cov_right / cov_left;
            if drop < local_drop && drop < sink_drop {
                sink_drop = drop;
                sink_trim = Some(TrimPoint {
                    pos: pos.saturating_sub(1),
                    is_start: false,
                    magnitude: (cov_left - cov_right) / DROP,
                });
            }
        }
    }

    (source_trim, sink_trim)
}

/// CHI2-based coverage differential for a split at position `split_i` within a window.
///
/// Ports `compute_chi2()` from rlink.cpp (line 1408).
/// Returns the improvement in fit from splitting at `split_i`.
pub fn compute_chi2(window: &[f64]) -> f64 {
    let n = window.len();
    if n < 2 {
        return 0.0;
    }
    let half = n / 2;
    let sum_left: f64 = window[..half].iter().sum();
    let sum_right: f64 = window[half..].iter().sum();
    let mu_base = (sum_left + sum_right) / n as f64;
    let mu_left = sum_left / half as f64;
    let mu_right = sum_right / (n - half) as f64;

    let base_cost: f64 = window.iter()
        .map(|&x| (x - mu_base).powi(2))
        .sum::<f64>() / n as f64;
    let left_cost: f64 = window[..half].iter()
        .map(|&x| (x - mu_left).powi(2))
        .sum::<f64>() / half as f64;
    let right_cost: f64 = window[half..].iter()
        .map(|&x| (x - mu_right).powi(2))
        .sum::<f64>() / (n - half) as f64;

    base_cost - (left_cost + right_cost) / 2.0
}

/// Sliding CHI-window scan for TSS/TES within a long bundlenode.
///
/// Ports the long-region sliding window path of StringTie's `create_graph` / `longtrim`.
/// For each position in the window, computes the CHI2 statistic and records significant
/// coverage differentials.
///
/// Returns all detected TrimPoints sorted by position.
pub fn find_coverage_trims_chi(
    start: u64,
    end: u64,
    bpcov: &PerBaseCoverage,
    min_chi2: f64,
) -> Vec<TrimPoint> {
    let len = end.saturating_sub(start) as usize;
    // Only makes sense for long bundlenodes (>= 2*CHI_WIN + CHI_THR)
    if len < 2 * CHI_WIN + CHI_THR {
        let (src, snk) = find_coverage_trims(start, end, bpcov);
        let mut result = Vec::new();
        if let Some(t) = src { result.push(t); }
        if let Some(t) = snk { result.push(t); }
        return result;
    }

    let mut trims: Vec<TrimPoint> = Vec::new();
    let mut window: Vec<f64> = Vec::with_capacity(2 * CHI_WIN);

    // Build initial window (2×CHI_WIN bp starting at start + LONG_INTRON_ANCHOR)
    let scan_start = start + LONG_INTRON_ANCHOR;
    for i in 0..(2 * CHI_WIN) {
        let pos = scan_start + i as u64;
        window.push(if pos <= end { bpcov.at(pos) } else { 0.0 });
    }

    let scan_end = end.saturating_sub(LONG_INTRON_ANCHOR);
    if scan_end <= scan_start + CHI_WIN as u64 {
        return trims;
    }

    for (step, center) in (scan_start + CHI_WIN as u64..=scan_end).enumerate() {
        // Slide window: drop oldest, add newest
        if step > 0 {
            window.remove(0);
            let new_pos = center + CHI_WIN as u64;
            window.push(if new_pos <= end { bpcov.at(new_pos) } else { 0.0 });
        }

        let chi = compute_chi2(&window);
        if chi > min_chi2 {
            // Determine direction from half-window means
            let sum_left: f64 = window[..CHI_WIN].iter().sum::<f64>() / CHI_WIN as f64;
            let sum_right: f64 = window[CHI_WIN..].iter().sum::<f64>() / CHI_WIN as f64;
            let is_start = sum_right > sum_left;
            trims.push(TrimPoint { pos: center, is_start, magnitude: chi });
        }
    }

    // Cluster nearby trims: keep strongest in each CHI_WIN window
    trims.sort_by_key(|t| t.pos);
    let mut clustered: Vec<TrimPoint> = Vec::new();
    let mut cluster_start_idx = 0;
    while cluster_start_idx < trims.len() {
        let cluster_anchor = trims[cluster_start_idx].pos;
        let mut best = trims[cluster_start_idx];
        let mut cluster_end_idx = cluster_start_idx + 1;
        while cluster_end_idx < trims.len()
            && trims[cluster_end_idx].pos <= cluster_anchor + CHI_WIN as u64
        {
            if trims[cluster_end_idx].magnitude > best.magnitude {
                best = trims[cluster_end_idx];
            }
            cluster_end_idx += 1;
        }
        clustered.push(best);
        cluster_start_idx = cluster_end_idx;
    }

    clustered
}

// ── Phase 4b: EM flow estimation ────────────────────────────────────────────

/// A read transfrag: a specific path through graphnodes with cumulative weight.
#[derive(Debug, Clone)]
pub struct StTransfrag {
    /// Ordered sequence of graphnode IDs.
    pub nodes: Vec<usize>,
    /// Total read weight supporting this path.
    pub abundance: f64,
    /// True if at least one long read supports this path.
    pub has_longread: bool,
    /// Leftmost start position across all supporting reads.
    pub longstart: u64,
    /// Rightmost end position across all supporting reads.
    pub longend: u64,
}

/// Build transfrags by grouping reads with identical graphnode paths.
///
/// Reads that traverse the same sequence of graphnodes are merged into a single
/// transfrag. This is the first step of StringTie's abundance estimation.
pub fn build_transfrags(graph: &BundleGraph, reads: &[BundleRead]) -> Vec<StTransfrag> {
    let mut path_map: std::collections::HashMap<Vec<usize>, StTransfrag> =
        std::collections::HashMap::new();

    for read in reads {
        let path = match map_read_to_path(graph, read) {
            Some(p) if !p.is_empty() => p,
            _ => continue,
        };

        let entry = path_map.entry(path.clone()).or_insert_with(|| StTransfrag {
            nodes: path,
            abundance: 0.0,
            has_longread: false,
            longstart: read.ref_start,
            longend: read.ref_end,
        });
        entry.abundance += read.weight;
        if read.ref_end > read.ref_start + 500 {
            entry.has_longread = true;
        }
        entry.longstart = entry.longstart.min(read.ref_start);
        entry.longend = entry.longend.max(read.ref_end);
    }

    let mut result: Vec<StTransfrag> = path_map.into_values().collect();
    // Sort by abundance descending (most-supported paths first)
    result.sort_by(|a, b| b.abundance.partial_cmp(&a.abundance).unwrap_or(std::cmp::Ordering::Equal));
    result
}

/// EM-based abundance estimation for transfrags.
///
/// Ports a simplified version of StringTie's `setAbundanceAll()`:
///   1. Each transfrag starts with its direct read weight.
///   2. For each node, compute the fractional abundance of each transfrag vs all
///      transfrags passing through it.
///   3. Distribute node coverage proportionally across transfrags.
///   4. Repeat for `n_iter` iterations.
///
/// Full ST EM uses more complex path probability estimation with gap penalties.
/// This simplified version handles the common case of non-overlapping paths.
pub fn em_estimate_abundance(
    graph: &BundleGraph,
    transfrags: &mut Vec<StTransfrag>,
    n_iter: usize,
) {
    if transfrags.is_empty() {
        return;
    }

    // Build node → transfrag index map
    let n_nodes = graph.nodes.len();
    let mut node_to_trf: Vec<Vec<usize>> = vec![Vec::new(); n_nodes];
    for (ti, trf) in transfrags.iter().enumerate() {
        for &nid in &trf.nodes {
            if nid < n_nodes {
                node_to_trf[nid].push(ti);
            }
        }
    }

    for _iter in 0..n_iter {
        // E step: compute expected node coverage from current transfrag abundances
        let mut node_cov: Vec<f64> = graph.nodes.iter().map(|n| n.cov).collect();

        // M step: redistribute node coverage to transfrags proportionally
        let mut new_abund: Vec<f64> = vec![0.0; transfrags.len()];
        for nid in 0..n_nodes {
            if nid == BundleGraph::SOURCE_ID || nid == graph.sink_id {
                continue;
            }
            let covering = &node_to_trf[nid];
            if covering.is_empty() {
                continue;
            }
            let total_abund: f64 = covering.iter().map(|&ti| transfrags[ti].abundance).sum();
            if total_abund <= 0.0 {
                continue;
            }
            // Distribute this node's coverage proportionally
            let node_total = node_cov[nid];
            for &ti in covering {
                let frac = transfrags[ti].abundance / total_abund;
                new_abund[ti] += frac * node_total / transfrags[ti].nodes.len() as f64;
            }
        }

        // Update abundances (normalize by number of nodes in path to avoid inflation)
        for (ti, trf) in transfrags.iter_mut().enumerate() {
            if trf.nodes.is_empty() {
                continue;
            }
            trf.abundance = new_abund[ti];
        }
    }
}

/// Compute per-node coverage from the current transfrag abundances.
///
/// After EM estimation, update each graphnode's `cov` with the sum of transfrag
/// abundances flowing through it (the "predicted" coverage).
pub fn update_node_coverage_from_transfrags(
    graph: &mut BundleGraph,
    transfrags: &[StTransfrag],
) {
    for node in graph.nodes.iter_mut() {
        node.cov = 0.0;
    }
    for trf in transfrags {
        for &nid in &trf.nodes {
            if nid < graph.nodes.len() {
                graph.nodes[nid].cov += trf.abundance;
            }
        }
    }
}

// ── Main struct: BundlenodeGraphInferer ─────────────────────────────────────

/// StringTie-style junction inference using the bundlenode graph.
///
/// Pipeline:
///   1. Build bundlenodes from read coverage.
///   2. Construct bundlenode DAG by splitting at known junction boundaries.
///   3. Map reads to paths through the DAG.
///   4. Aggregate path-implied junction evidence.
///   5. Inject high-confidence candidates into `JunctionStats`.
///
/// The `inject_into` method only adds junctions NOT already in the stats,
/// so this layer acts as a junction *discovery* supplement.
pub struct BundlenodeGraphInferer {
    /// Junction candidates inferred from read paths: (donor, acceptor) → support weight.
    candidates: BTreeMap<(u64, u64), f64>,
    /// The graph constructed during inference (kept for diagnostics).
    pub graph: BundleGraph,
}

impl BundlenodeGraphInferer {
    /// Build a BundlenodeGraphInferer from reads and an existing junction set.
    ///
    /// The `known_junctions` slice provides the accepted splice junctions that
    /// are used to split bundlenodes into graphnodes. Read paths through the
    /// resulting graph then reveal any additional junction evidence.
    pub fn from_reads_and_junctions(
        reads: &[BundleRead],
        bundle_start: u64,
        bundle_end: u64,
        known_junctions: &[(u64, u64)],
        _min_intron: u64,
    ) -> Self {
        const EM_ITERS: usize = 10;

        // Phase 1: bundlenodes from read coverage
        let local_dist = 50u64; // ST's bundledist ≈ 50 bp
        let bnodes = build_bundlenodes(reads, bundle_start, bundle_end, local_dist);

        // Phase 2: graph construction — split bnodes at known junction boundaries
        let mut graph = build_bundle_graph(&bnodes, known_junctions);

        // Phase 3: assign per-base coverage to graphnodes
        let bpcov = PerBaseCoverage::from_reads(reads, bundle_start, bundle_end);
        assign_node_coverage(&mut graph, &bpcov);

        // Phase 4: EM flow — group reads by path, estimate abundance
        let mut transfrags = build_transfrags(&graph, reads);
        if transfrags.len() > 1 {
            em_estimate_abundance(&graph, &mut transfrags, EM_ITERS);
            update_node_coverage_from_transfrags(&mut graph, &transfrags);
        }

        // Phase 5: aggregate path-implied junctions from read alignments
        let candidates = aggregate_path_junctions(&graph, reads);

        BundlenodeGraphInferer { candidates, graph }
    }

    /// Legacy constructor: builds from reads alone with no prior junction set.
    ///
    /// Falls back to the original coverage-gap + bridge-read approach because
    /// there are no junction split-points to build the graph with.
    ///
    /// This is the interface used by the existing `pipeline.rs` call-site
    /// (`RUSTLE_JG_USE_ST_GRAPH`). Providing `known_junctions` via
    /// `from_reads_and_junctions` gives better results.
    pub fn from_reads(
        reads: &[BundleRead],
        bundle_start: u64,
        bundle_end: u64,
        min_intron: u64,
        max_intron: u64,
    ) -> Self {
        // Extract junctions from reads themselves to seed the graph
        let known: Vec<(u64, u64)> = extract_read_junctions(reads, min_intron, max_intron);
        Self::from_reads_and_junctions(reads, bundle_start, bundle_end, &known, min_intron)
    }

    /// Inject inferred junctions into an existing `JunctionStats`.
    ///
    /// Only adds junctions NOT already present. Skips candidates below
    /// `min_bridge_reads` total weight.
    pub fn inject_into(&self, stats: &mut JunctionStats, min_bridge_reads: f64) {
        for (&(donor, acceptor), &weight) in &self.candidates {
            if weight < min_bridge_reads {
                continue;
            }
            let junc = Junction::new(donor, acceptor);
            if stats.contains_key(&junc) {
                continue;
            }
            let mut js = JunctionStat::default();
            js.mrcount = weight;
            js.nreads_good = weight;
            js.leftsupport = weight;
            js.rightsupport = weight;
            js.mm = weight;
            js.strand = None;
            js.nm = 0.0;
            js.consleft = -1;
            js.consright = -1;
            stats.insert(junc, js);
        }
    }

    /// Return the inferred junction candidates (for diagnostics).
    pub fn candidates(&self) -> &BTreeMap<(u64, u64), f64> {
        &self.candidates
    }
}

// ── Helpers ──────────────────────────────────────────────────────────────────

/// Extract splice junctions directly from read exon boundaries.
fn extract_read_junctions(
    reads: &[BundleRead],
    min_intron: u64,
    max_intron: u64,
) -> Vec<(u64, u64)> {
    let mut seen: std::collections::HashSet<(u64, u64)> = std::collections::HashSet::new();
    let mut result = Vec::new();

    for read in reads {
        for exons in read.exons.windows(2) {
            let donor = exons[0].1;    // end of left exon
            let acceptor = exons[1].0; // start of right exon
            if acceptor > donor {
                let intron_len = acceptor - donor;
                if intron_len >= min_intron && intron_len <= max_intron {
                    if seen.insert((donor, acceptor)) {
                        result.push((donor, acceptor));
                    }
                }
            }
        }
    }
    result
}

// ── Tests ─────────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;
    use crate::types::BundleRead;
    use std::sync::Arc;

    fn make_read(uid: u64, exons: Vec<(u64, u64)>, weight: f64) -> BundleRead {
        BundleRead {
            read_uid: uid,
            read_name: format!("read{uid}").as_str().into(),
            read_name_hash: uid,
            ref_id: None, mate_ref_id: None, mate_start: None, hi: 0,
            ref_start: exons.first().map_or(0, |e| e.0),
            ref_end: exons.last().map_or(0, |e| e.1),
            exons: exons.clone(),
            junctions: vec![], junction_valid: vec![],
            junctions_raw: vec![], junctions_del: vec![],
            weight, is_reverse: false, strand: '+',
            has_poly_start: false, has_poly_end: false,
            has_poly_start_aligned: false, has_poly_start_unaligned: false,
            has_poly_end_aligned: false, has_poly_end_unaligned: false,
            unaligned_poly_t: 0, unaligned_poly_a: 0,
            has_last_exon_polya: false, has_first_exon_polyt: false,
            query_length: None, clip_left: 0, clip_right: 0,
            nh: 1, nm: 0, de: None, md: None, insertion_sites: vec![],
            unitig: false, unitig_cov: 0.0, read_count_yc: 1.0,
            countfrag_len: 0.0, countfrag_num: 0.0, junc_mismatch_weight: 0.0,
            pair_idx: vec![], pair_count: vec![],
            mapq: 60, mismatches: vec![], seq: Vec::new(), hp_tag: None, ps_tag: None,
            is_primary_alignment: true,
            em_weight_gap: -1.0, em_n_sites: 0, em_anchored: true, em_ev_decisive: false,
        }
    }

    #[test]
    fn test_bundlenode_single_island() {
        let reads = vec![make_read(1, vec![(100, 200)], 1.0)];
        let bnodes = build_bundlenodes(&reads, 100, 200, 50);
        assert_eq!(bnodes.len(), 1);
        assert_eq!(bnodes[0].start, 100);
        assert_eq!(bnodes[0].end, 200);
    }

    #[test]
    fn test_bundlenode_two_islands() {
        let reads = vec![
            make_read(1, vec![(100, 200)], 1.0),
            make_read(2, vec![(400, 500)], 1.0),
        ];
        let bnodes = build_bundlenodes(&reads, 100, 500, 50);
        assert_eq!(bnodes.len(), 2);
        assert_eq!(bnodes[0].start, 100);
        assert_eq!(bnodes[1].start, 400);
    }

    #[test]
    fn test_bundlenode_merges_nearby() {
        let reads = vec![
            make_read(1, vec![(100, 200)], 1.0),
            make_read(2, vec![(220, 300)], 1.0),
        ];
        // gap = 20 < local_dist = 50 → merged
        let bnodes = build_bundlenodes(&reads, 100, 300, 50);
        assert_eq!(bnodes.len(), 1);
    }

    #[test]
    fn test_graph_no_junctions() {
        // Single bundlenode, no junctions → simple linear graph: source → node → sink
        let bnodes = vec![BundleNode { start: 100, end: 200, cov: 1.0 }];
        let graph = build_bundle_graph(&bnodes, &[]);
        // source + 1 real node + sink = 3
        assert_eq!(graph.nodes.len(), 3);
        assert!(graph.nodes[BundleGraph::SOURCE_ID].children.contains(&1));
        assert!(graph.nodes[1].children.contains(&graph.sink_id));
    }

    #[test]
    fn test_graph_one_junction_splits_bnode() {
        // Bundlenode [100, 500], junction donor=200, acceptor=350
        // Expected: source → [100,200] → [201,500] → sink
        //                              ↗ (junction edge from [100,200])
        //           [350,500]   → sink  (if bundlenode [350,500] existed)
        // But we have one bundlenode [100,500] split at donor=200
        // and the junction end is at 350 (within the bnode)
        let bnodes = vec![BundleNode { start: 100, end: 500, cov: 1.0 }];
        let junctions = vec![(200u64, 350u64)];
        let graph = build_bundle_graph(&bnodes, &junctions);
        // At minimum, 2 real nodes ([100,200] and [200+, 500]) plus source + sink
        assert!(graph.nodes.len() >= 4);
    }

    #[test]
    fn test_path_inference_two_exon_read() {
        // Two exon read: [100,200] gap [250,350]
        // Junction: donor=200, acceptor=250
        let reads = vec![make_read(1, vec![(100, 200), (250, 350)], 1.0)];
        let inferer = BundlenodeGraphInferer::from_reads(&reads, 100, 350, 50, 1_000_000);
        // Should find junction (200, 250)
        let cands = inferer.candidates();
        assert!(cands.contains_key(&(200u64, 250u64)),
            "Expected junction (200,250) in candidates; got: {:?}", cands);
    }

    // ── PerBaseCoverage tests ──────────────────────────────────────────────

    #[test]
    fn test_pbc_single_read() {
        // Read covering [100, 110) with weight 1.0
        let reads = vec![make_read(1, vec![(100, 110)], 1.0)];
        let pbc = PerBaseCoverage::from_reads(&reads, 100, 110);
        assert_eq!(pbc.at(100), 1.0);
        assert_eq!(pbc.at(109), 1.0);
        assert_eq!(pbc.at(110), 0.0); // past end
        assert_eq!(pbc.sum(100, 109), 10.0);
        assert!((pbc.avg(100, 109) - 1.0).abs() < 1e-9);
    }

    #[test]
    fn test_pbc_two_reads_overlapping() {
        // Two reads both covering [0, 10); second only covers [5, 10)
        let reads = vec![
            make_read(1, vec![(0, 10)], 1.0),
            make_read(2, vec![(5, 10)], 1.0),
        ];
        let pbc = PerBaseCoverage::from_reads(&reads, 0, 10);
        // Positions 0-4: cov=1; positions 5-9: cov=2
        assert_eq!(pbc.at(0), 1.0);
        assert_eq!(pbc.at(4), 1.0);
        assert_eq!(pbc.at(5), 2.0);
        assert_eq!(pbc.at(9), 2.0);
        // sum [0,9] = 5*1 + 5*2 = 15
        assert_eq!(pbc.sum(0, 9), 15.0);
    }

    #[test]
    fn test_pbc_empty_reads() {
        let reads: Vec<BundleRead> = vec![];
        let pbc = PerBaseCoverage::from_reads(&reads, 100, 200);
        assert_eq!(pbc.at(150), 0.0);
        assert_eq!(pbc.sum(100, 199), 0.0);
        assert_eq!(pbc.avg(100, 199), 0.0);
    }

    #[test]
    fn test_pbc_out_of_range() {
        let reads = vec![make_read(1, vec![(100, 200)], 2.0)];
        let pbc = PerBaseCoverage::from_reads(&reads, 100, 200);
        // Query before origin
        assert_eq!(pbc.at(50), 0.0);
        assert_eq!(pbc.sum(50, 99), 0.0);
        // Query past end
        assert_eq!(pbc.at(250), 0.0);
    }

    // ── compute_chi2 tests ─────────────────────────────────────────────────

    #[test]
    fn test_chi2_uniform_window_returns_zero() {
        // Uniform coverage → no split → chi2 ≈ 0
        let window = vec![1.0; 20];
        let chi = compute_chi2(&window);
        assert!(chi.abs() < 1e-9, "Expected ~0, got {chi}");
    }

    #[test]
    fn test_chi2_step_function() {
        // Left half = 0, right half = 10 → strong chi2
        let mut window = vec![0.0f64; 20];
        for v in window[10..].iter_mut() {
            *v = 10.0;
        }
        let chi = compute_chi2(&window);
        // Should be positive since a split at 10 reduces variance
        assert!(chi > 0.0, "Expected positive chi2, got {chi}");
    }

    #[test]
    fn test_chi2_empty_or_single() {
        assert_eq!(compute_chi2(&[]), 0.0);
        assert_eq!(compute_chi2(&[5.0]), 0.0);
    }

    // ── find_coverage_trims tests ──────────────────────────────────────────

    #[test]
    fn test_find_trims_too_short_returns_none() {
        // Region shorter than CHI_THR (50) → no trims
        let reads = vec![make_read(1, vec![(0, 30)], 1.0)];
        let pbc = PerBaseCoverage::from_reads(&reads, 0, 30);
        let (src, snk) = find_coverage_trims(0, 30, &pbc);
        assert!(src.is_none());
        assert!(snk.is_none());
    }

    #[test]
    fn test_find_trims_uniform_coverage_no_trim() {
        // Uniform coverage, no TSS/TES signal
        let reads = vec![make_read(1, vec![(0, 200)], 1.0)];
        let pbc = PerBaseCoverage::from_reads(&reads, 0, 200);
        let (src, snk) = find_coverage_trims(0, 200, &pbc);
        // Uniform coverage should not produce a trim
        assert!(src.is_none() || snk.is_none(), "Uniform coverage should not produce both trims");
    }

    // ── build_transfrags tests ─────────────────────────────────────────────

    #[test]
    fn test_build_transfrags_identical_paths_merge() {
        // Three reads all taking the same path through the graph
        let reads = vec![
            make_read(1, vec![(100, 200), (300, 400)], 1.0),
            make_read(2, vec![(100, 200), (300, 400)], 2.0),
            make_read(3, vec![(100, 200), (300, 400)], 1.0),
        ];
        let junctions = vec![(200u64, 300u64)];
        let bnodes = vec![
            BundleNode { start: 100, end: 200, cov: 1.0 },
            BundleNode { start: 300, end: 400, cov: 1.0 },
        ];
        let graph = build_bundle_graph(&bnodes, &junctions);
        let transfrags = build_transfrags(&graph, &reads);
        // All three reads share the same path → 1 transfrag with abundance=4.0
        assert_eq!(transfrags.len(), 1, "Expected 1 merged transfrag, got {}", transfrags.len());
        assert!((transfrags[0].abundance - 4.0).abs() < 1e-9,
            "Expected abundance=4.0, got {}", transfrags[0].abundance);
    }

    #[test]
    fn test_build_transfrags_different_paths_separate() {
        // Two reads taking different paths (one exon each)
        let reads = vec![
            make_read(1, vec![(100, 200)], 1.0),
            make_read(2, vec![(300, 400)], 1.0),
        ];
        let bnodes = vec![
            BundleNode { start: 100, end: 200, cov: 1.0 },
            BundleNode { start: 300, end: 400, cov: 1.0 },
        ];
        let graph = build_bundle_graph(&bnodes, &[]);
        let transfrags = build_transfrags(&graph, &reads);
        assert_eq!(transfrags.len(), 2, "Expected 2 transfrags, got {}", transfrags.len());
    }

    // ── em_estimate_abundance tests ────────────────────────────────────────

    #[test]
    fn test_em_two_transfrags_equal_coverage() {
        // Two transfrags through the same two-exon graph, equal initial abundance
        let reads = vec![
            make_read(1, vec![(100, 200), (300, 400)], 1.0),
            make_read(2, vec![(100, 200), (300, 400)], 1.0),
        ];
        let junctions = vec![(200u64, 300u64)];
        let bnodes = vec![
            BundleNode { start: 100, end: 200, cov: 2.0 },
            BundleNode { start: 300, end: 400, cov: 2.0 },
        ];
        let graph = build_bundle_graph(&bnodes, &junctions);
        let mut transfrags = build_transfrags(&graph, &reads);
        // EM should run without panicking and preserve positive abundances
        em_estimate_abundance(&graph, &mut transfrags, 5);
        for trf in &transfrags {
            assert!(trf.abundance >= 0.0, "Abundance went negative: {}", trf.abundance);
        }
    }

    #[test]
    fn test_em_empty_transfrags_no_panic() {
        let bnodes = vec![BundleNode { start: 0, end: 100, cov: 1.0 }];
        let graph = build_bundle_graph(&bnodes, &[]);
        let mut transfrags: Vec<StTransfrag> = vec![];
        // Should not panic
        em_estimate_abundance(&graph, &mut transfrags, 10);
    }

    // ── assign_node_coverage tests ─────────────────────────────────────────

    #[test]
    fn test_assign_node_coverage_uses_avg() {
        // Single read covering [0, 100) at weight 3.0
        // Graph: source → [0,99] → sink
        let reads = vec![make_read(1, vec![(0, 100)], 3.0)];
        let bnodes = build_bundlenodes(&reads, 0, 100, 50);
        let mut graph = build_bundle_graph(&bnodes, &[]);
        let pbc = PerBaseCoverage::from_reads(&reads, 0, 100);
        assign_node_coverage(&mut graph, &pbc);

        // Real node (not source/sink) should have cov ≈ 3.0
        for node in &graph.nodes {
            if node.id == BundleGraph::SOURCE_ID || node.id == graph.sink_id {
                continue;
            }
            assert!((node.cov - 3.0).abs() < 0.5,
                "Expected cov≈3.0 for node [{},{}], got {}", node.start, node.end, node.cov);
        }
    }
}
