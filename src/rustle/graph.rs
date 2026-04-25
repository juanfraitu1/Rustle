//! Splice graph: nodes = exon segments, edges = collinear + junction edges.

use crate::bitset::NodeSet;
use crate::bitvec::GBitVec;
use crate::coord::len_half_open;
use crate::types::{AssemblyMode, DetHashMap as HashMap, DetHashSet as HashSet};

/// Role classification for a splice-graph node.
///
/// In the default (disjoint-split) build, every node is `Primary` — it
/// represents an exclusive genomic interval with normal participation in
/// all downstream systems.
///
/// Overlap-model builds (StringTie-parity; gated behind
/// `RUSTLE_OVERLAP_NODE_SCAFFOLD` / `RUSTLE_PURE_OVERLAP` variants) add
/// nodes whose genomic range intentionally overlaps with `Primary`
/// siblings. Each overlap role selects different participation flags:
///
/// - `accepts_reads` — whether read-to-node mapping (bundle2graph,
///   `substitute_alias_for_spanning_run`) may route reads here.
/// - `accrues_coverage` — whether per-base coverage aggregation sums bp
///   mass onto this node. Set to false for overlap nodes so the same bp
///   isn't counted twice.
/// - `prune_autoattach` — whether `prune_graph_nodes*` source/sink
///   auto-attach loops consider this node. Overlap nodes receive edges
///   programmatically; auto-attach would create spurious short paths.
/// - `longtrim_schedule` — whether per-bundlenode longtrim scheduling
///   iterates this node. Overlap nodes share a bundlenode with their
///   Primary siblings; longtrim runs per Primary.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum NodeRole {
    /// Normal graph node. Participates in everything.
    Primary,
    /// StringTie-node-33 equivalent: full-span anchor OVERLAPPING one or
    /// more Primary disjoint splits within the same bundlenode.
    OverlapAnchor,
    /// StringTie-node-34/35 equivalent: junction-acceptor entry node,
    /// spans [acceptor, endbundle) and is entered only via the junction
    /// landing at its start. Overlaps downstream Primary splits.
    JunctionEntry,
}

impl Default for NodeRole {
    fn default() -> Self {
        NodeRole::Primary
    }
}

impl NodeRole {
    /// Whether STANDARD read-to-node mapping (via bundle2graph in
    /// coord-sorted order) considers this node. Returns true only for
    /// `Primary` because overlap nodes would coord-intersect the same
    /// bp as Primary siblings; including both would route each read to
    /// multiple nodes and break exon-to-node exclusivity.
    ///
    /// In pure-overlap mode (no Primary siblings), use
    /// [`NodeRole::accepts_reads_pure_overlap`] instead.
    #[inline]
    pub fn accepts_reads(self) -> bool {
        matches!(self, NodeRole::Primary)
    }

    /// Pure-overlap-mode read routing: overlap nodes DO accept reads
    /// because they're the only nodes covering their bp range in this
    /// mode. Used by a separate pure-overlap bundle2graph variant.
    #[inline]
    pub fn accepts_reads_pure_overlap(self) -> bool {
        true
    }

    /// Whether per-base coverage aggregation sums onto this node.
    ///
    /// - `Primary`: yes
    /// - `OverlapAnchor`: no — its genomic range is covered by Primary
    ///   splits (overlap-alongside) or by JunctionEntry children
    ///   (pure-overlap). Aggregating would double-count.
    /// - `JunctionEntry`: no — same reasoning: its range overlaps the
    ///   OverlapAnchor's, which aggregates (or in overlap-alongside
    ///   mode the Primary splits do).
    #[inline]
    pub fn accrues_coverage(self) -> bool {
        matches!(self, NodeRole::Primary)
    }

    /// Whether prune source/sink auto-attach considers this node.
    /// Overlap nodes receive edges programmatically (at bundlenode
    /// creation time); auto-attach would create spurious single-node
    /// paths.
    #[inline]
    pub fn prune_autoattach(self) -> bool {
        matches!(self, NodeRole::Primary)
    }

    /// Whether longtrim bundle scheduling iterates this node.
    /// Longtrim runs per Primary split; overlap nodes share their
    /// bundlenode with a Primary and don't need their own schedule.
    #[inline]
    pub fn longtrim_schedule(self) -> bool {
        matches!(self, NodeRole::Primary)
    }

    /// Whether this node is any flavor of overlap (not Primary).
    #[inline]
    pub fn is_overlap(self) -> bool {
        !matches!(self, NodeRole::Primary)
    }
}

/// Single node in the splice graph (exon segment).
#[derive(Debug, Clone)]
pub struct GraphNode {
    pub node_id: usize,
    pub start: u64,
    pub end: u64,
    /// The bundlenode that created this graph node (create_graphnode -> bundle2graph).
    pub source_bnode: Option<usize>,
    /// Total base-coverage mass on this node.
    ///
    /// This is the Rust/`CGraphnode::cov` quantity used by `update_abundance`,
    /// `process_transfrags`, and `noderate = coverage / nodecov`. It lives in
    /// bp-weight units: sum(read_weight * overlapped_bp) across mapped reads, or the
    /// equivalent sum rematerialized from `bpcov`.
    ///
    /// It is intentionally not the bundlenode's aggregate coverage estimate.
    pub coverage: f64,
    pub children: crate::bitset::SmallBitset,
    pub parents: crate::bitset::SmallBitset,
    pub hardstart: bool,
    pub hardend: bool,
    /// Read-signal alt-TTS signal (terminal-read cluster at junction-donor
    /// coord). Consumed ONLY by path_extract's HARDEND_TERMINATE gate to
    /// avoid disturbing hardend-reading filters (pairwise, RI, map_reads
    /// split decisions). Set by discover_terminal_donor_hardends.
    pub alt_tts_end: bool,
    /// Sum of poly_start_unaligned over transfrags that start at this node (annotate).
    pub poly_start_unaligned_total: u32,
    /// Sum of poly_start_aligned over transfrags that start at this node .
    pub poly_start_aligned_total: u32,
    /// Sum of poly_end_unaligned over transfrags that end at this node (annotate).
    pub poly_end_unaligned_total: u32,
    /// Sum of poly_end_aligned over transfrags that end at this node .
    pub poly_end_aligned_total: u32,
    /// compute_nodecov: incoming abundance through this node.
    pub abundin: f64,
    /// compute_nodecov: outgoing abundance through this node.
    pub abundout: f64,
    /// `nodecov = max(abundin, abundout)` in abundance units.
    pub nodecov: f64,
    /// Sum of long-read transfrag abundance through this node.
    pub longcov: f64,
    /// `noderate = coverage / nodecov` (bp-weight / abundance).
    pub noderate: f64,
    pub trf_ids: Vec<usize>,
    pub childpat: Option<crate::bitset::SmallBitset>,
    pub parentpat: Option<crate::bitset::SmallBitset>,
    /// Coverage-drop contrast at a longtrim split point.
    /// Set on both sides of a longtrim split (the hardend upstream node carries
    /// this value on its .longtrim_cov, and the hardstart downstream node does too).
    /// Used by map_reads read-split gate to distinguish strong gene-end boundaries
    /// from weak alt-TTS/alt-TSS boundaries.
    /// 0.0 when no longtrim split applies to this node.
    pub longtrim_cov: f64,
    /// Role classification for participation in read routing, coverage
    /// aggregation, prune auto-attach, and longtrim scheduling. Default
    /// `Primary` matches the historical behavior of all nodes.
    /// See [`NodeRole`] for semantics.
    pub role: NodeRole,
}

impl GraphNode {
    pub fn new(node_id: usize, start: u64, end: u64) -> Self {
        Self {
            node_id,
            start,
            end,
            source_bnode: None,
            coverage: 0.0,
            children: crate::bitset::SmallBitset::empty(),
            parents: crate::bitset::SmallBitset::empty(),
            hardstart: false,
            hardend: false,
            alt_tts_end: false,
            poly_start_unaligned_total: 0,
            poly_start_aligned_total: 0,
            poly_end_unaligned_total: 0,
            poly_end_aligned_total: 0,
            abundin: 0.0,
            abundout: 0.0,
            nodecov: 0.0,
            longcov: 0.0,
            noderate: 1.0,
            trf_ids: Vec::new(),
            childpat: None,
            parentpat: None,
            longtrim_cov: 0.0,
            role: NodeRole::Primary,
        }
    }

    /// Convenience: alias for the common check.
    #[inline]
    pub fn is_overlap(&self) -> bool {
        self.role.is_overlap()
    }

    pub fn length(&self) -> u64 {
        len_half_open(self.start, self.end)
    }
}

/// Transfrag: path through graph (list of node IDs) with abundance.
#[derive(Debug, Clone)]
pub struct FlowBranch {
    pub node: usize,
    pub contnode: usize,
    pub abundance: f64,
}

/// the original algorithm naming convention (`CPath` in header).
pub type CPath = FlowBranch;

#[derive(Debug, Clone)]
pub struct GraphTransfrag {
    pub node_ids: Vec<usize>,
    /// Bitset representation of node_ids for fast membership testing (
    pub node_id_set: NodeSet,
    pub pattern: GBitVec,
    /// Long-read / mixed-mode flow mass for this transfrag.
    ///
    /// This is the quantity carried through `process_transfrags`, `long_max_flow`,
    /// `nodecov`, and `longcov` output. It is not per-base coverage.
    pub abundance: f64,
    /// Total supporting read mass independent of long/short flow buckets.
    pub read_count: f64,
    pub longstart: u64,
    pub longend: u64,
    pub longread: bool,
    /// `shortread` flag: at least one short-read contributes to this transfrag.
    pub shortread: bool,
    /// Short-read abundance (`srabund`); used for short-read / mixed flow.
    pub srabund: f64,
    /// Count of reads with aligned polyA/T at 5' end.
    pub poly_start_aligned: u16,
    /// Count of reads with aligned polyA/T at 3' end.
    pub poly_end_aligned: u16,
    /// Count of reads with polyA/T at 5' (poly_start_unaligned)
    pub poly_start_unaligned: u16,
    /// Count of reads with polyA/T at 3' (poly_end_unaligned)
    pub poly_end_unaligned: u16,
    /// True if this transfrag came from a read split at a killed junction (V99 killed_junction_orphan).
    pub killed_junction_orphan: bool,
    /// True if transfrag has a weak link (coverage drop between consecutive contiguous nodes; compute_weak).
    pub coverage_weak: bool,
    /// 1 if absorbed into another (merged into representative); 0 otherwise (tf.weak, reference line 11014).
    pub weak: u8,
    /// `real` flag: true when transfrag is considered complete/solid after incomplete handling.
    pub real: bool,
    /// Indices of transfrags absorbed into this one (representative only; group, reference line 2978).
    pub group: Vec<usize>,
    /// True if from reference/guide (guide); blocks absorption when kept is guide (ret=1), allows rep replacement (ret=3).
    pub guide: bool,
    /// Guide transcript id when this transfrag is guide-derived (t_eq equivalent anchor).
    pub guide_tid: Option<String>,
    /// True when selected as a long-read extraction seed (trflong behavior).
    pub trflong_seed: bool,
    /// Mixed-mode parse order marker (usepath=-2-ntrf convention; -1 means unset).
    pub usepath: i32,
    /// transfrag.path equivalent for unresolved/incomplete transfrags.
    pub flow_paths: Vec<FlowBranch>,
    /// Selected branch index for current path-based flow step (usepath for path branches).
    pub flow_path_idx: i32,
    /// Nascent marker (nascent handling scaffold).
    pub nascent: bool,
    /// Lightweight provenance tag for seed-inventory tracing.
    pub origin_tag: Option<String>,
}

impl GraphTransfrag {
    pub fn new(node_ids: Vec<usize>, pattern_size: usize) -> Self {
        let mut node_id_set = NodeSet::with_capacity(pattern_size);
        for &nid in &node_ids {
            node_id_set.insert(nid);
        }
        let mut pattern = GBitVec::new(pattern_size);
        for &nid in &node_ids {
            pattern.set_bit(nid);
        }
        Self {
            node_ids,
            node_id_set,
            pattern,
            abundance: 0.0,
            read_count: 0.0,
            longstart: 0,
            longend: 0,
            longread: false,
            shortread: false,
            srabund: 0.0,
            poly_start_aligned: 0,
            poly_end_aligned: 0,
            poly_start_unaligned: 0,
            poly_end_unaligned: 0,
            killed_junction_orphan: false,
            coverage_weak: false,
            weak: 0,
            real: false,
            group: Vec::new(),
            guide: false,
            guide_tid: None,
            trflong_seed: false,
            usepath: -1,
            flow_paths: Vec::new(),
            flow_path_idx: -1,
            nascent: false,
            origin_tag: None,
        }
    }

    /// Rebuild pattern after graph modifications (pruning, trimming) to fix stale edge bits.
    /// transfrag patterns must reflect current graph topology for capacity network.
    pub fn rebuild_pattern(&mut self, graph: &Graph) {
        // Grow pattern if graph has expanded
        let psize = graph.pattern_size();
        if psize > self.pattern.len_bits() + (psize - self.pattern.len_bits()) {
            self.pattern.grow(psize);
        }
        // Clear and rebuild from current node_ids
        self.pattern.reset();
        for &nid in &self.node_ids {
            self.pattern.set_bit(nid);
        }
        graph.set_pattern_edges_for_path(&mut self.pattern, &self.node_ids);
    }
}

/// the original algorithm naming convention (`CGraphnode` in header).
pub type CGraphnode = GraphNode;
/// the original algorithm naming convention (`CTransfrag` in header).
pub type CTransfrag = GraphTransfrag;
/// the original algorithm naming convention (`CMTransfrag` in header merge mode scaffold).
pub type CMTransfrag = GraphTransfrag;

/// Splice graph: source(0), real nodes 1..n-1, sink(n-1).
#[derive(Debug)]
pub struct Graph {
    pub nodes: Vec<GraphNode>,
    pub source_id: usize,
    pub sink_id: usize,
    pub n_nodes: usize,
    /// (min(from,to), max(from,to)) -> edge_id for pattern bits
    pub gpos: HashMap<(usize, usize), usize>,
    pub edgeno: usize,
    /// Track the next available edge ID to avoid O(|E|) scans.
    next_edge_id: usize,
}

impl Graph {
    pub fn new() -> Self {
        Self {
            nodes: Vec::new(),
            source_id: 0,
            sink_id: 0,
            n_nodes: 0,
            gpos: Default::default(),
            edgeno: 0,
            next_edge_id: 0,
        }
    }

    pub fn add_node(&mut self, start: u64, end: u64) -> &mut GraphNode {
        let node_id = self.nodes.len();
        self.nodes.push(GraphNode::new(node_id, start, end));
        self.n_nodes = self.nodes.len();
        // Keep edge-bit namespace strictly above node-id namespace.
        // Use tracked next_edge_id instead of O(|E|) scan.
        if self.next_edge_id < self.n_nodes && !self.gpos.is_empty() {
            self.reindex_edge_bits_dense();
        }
        self.nodes.last_mut().unwrap()
    }

    pub fn add_edge(&mut self, from_id: usize, to_id: usize) {
        if from_id >= self.nodes.len() || to_id >= self.nodes.len() {
            return;
        }
        if !self.nodes[from_id].children.contains(to_id) {
            self.nodes[from_id].children.insert_grow(to_id);
        }
        if !self.nodes[to_id].parents.contains(from_id) {
            self.nodes[to_id].parents.insert_grow(from_id);
        }
        let key = (from_id.min(to_id), from_id.max(to_id));
        if !self.gpos.contains_key(&key) {
            let eid = self.next_edge_id.max(self.n_nodes);
            self.gpos.insert(key, eid);
            self.next_edge_id = eid.saturating_add(1);
            self.edgeno += 1;
        }
    }

    /// Connected components of the splice graph excluding source/sink.
    /// Uses undirected adjacency over existing edges.
    pub fn junction_connected_components(&self) -> Vec<Vec<usize>> {
        let mut comps: Vec<Vec<usize>> = Vec::new();
        if self.nodes.is_empty() {
            return comps;
        }
        let mut seen = vec![false; self.nodes.len()];
        for start in 0..self.nodes.len() {
            if start == self.source_id || start == self.sink_id {
                continue;
            }
            if seen[start] {
                continue;
            }
            let mut stack = vec![start];
            seen[start] = true;
            let mut comp = Vec::new();
            while let Some(u) = stack.pop() {
                comp.push(u);
                let node = &self.nodes[u];
                for v in node.children.ones().chain(node.parents.ones()) {
                    if v == self.source_id || v == self.sink_id {
                        continue;
                    }
                    if !seen[v] {
                        seen[v] = true;
                        stack.push(v);
                    }
                }
            }
            if !comp.is_empty() {
                comps.push(comp);
            }
        }
        comps
    }

    /// Build a subgraph containing source/sink and the specified internal nodes.
    /// Returns (subgraph, new_to_old, old_to_new).
    pub fn subgraph_from_nodes(&self, keep_nodes: &[usize]) -> (Graph, Vec<usize>, Vec<usize>) {
        let mut keep: Vec<usize> = keep_nodes
            .iter()
            .copied()
            .filter(|&n| n != self.source_id && n != self.sink_id)
            .collect();
        keep.sort_unstable();
        keep.dedup();

        let mut sub = Graph::new();
        let mut new_to_old: Vec<usize> = Vec::new();
        let mut old_to_new: Vec<usize> = vec![usize::MAX; self.nodes.len()];

        let push_clone = |old_id: usize,
                          sub: &mut Graph,
                          new_to_old: &mut Vec<usize>,
                          old_to_new: &mut [usize]| {
            if let Some(node) = self.nodes.get(old_id) {
                let mut cloned = node.clone();
                cloned.node_id = sub.nodes.len();
                cloned.parents.clear();
                cloned.children.clear();
                sub.nodes.push(cloned);
                sub.n_nodes = sub.nodes.len();
                new_to_old.push(old_id);
                if old_id < old_to_new.len() {
                    old_to_new[old_id] = sub.nodes.len() - 1;
                }
            }
        };

        push_clone(self.source_id, &mut sub, &mut new_to_old, &mut old_to_new);
        for &oid in &keep {
            push_clone(oid, &mut sub, &mut new_to_old, &mut old_to_new);
        }
        push_clone(self.sink_id, &mut sub, &mut new_to_old, &mut old_to_new);

        sub.source_id = 0;
        sub.sink_id = sub.nodes.len().saturating_sub(1);

        for &old_id in &new_to_old {
            let new_u = match old_to_new.get(old_id) {
                Some(&v) if v != usize::MAX => v,
                _ => continue,
            };
            if let Some(old_node) = self.nodes.get(old_id) {
                for child in old_node.children.ones().collect::<Vec<_>>() {
                    let new_v = match old_to_new.get(child) {
                        Some(&v) if v != usize::MAX => v,
                        _ => continue,
                    };
                    sub.add_edge(new_u, new_v);
                }
            }
        }
        sub.reindex_edge_bits_dense();
        (sub, new_to_old, old_to_new)
    }

    /// Remove edge from_id -> to_id (delete_connection).
    pub fn remove_edge(&mut self, from_id: usize, to_id: usize) {
        if from_id >= self.nodes.len() || to_id >= self.nodes.len() {
            return;
        }
        self.nodes[from_id].children.remove(to_id);
        self.nodes[to_id].parents.remove(from_id);
        let key = (from_id.min(to_id), from_id.max(to_id));
        if self.gpos.remove(&key).is_some() {
            self.edgeno = self.edgeno.saturating_sub(1);
        }
    }

    /// Pattern-bit index for an edge if it exists in graph (min/max node order).
    pub fn edge_bit_index(&self, from: usize, to: usize) -> Option<usize> {
        let key = if from < to { (from, to) } else { (to, from) };
        self.gpos.get(&key).copied()
    }

    pub fn pattern_size(&self) -> usize {
        let max_bit = self
            .gpos
            .values()
            .copied()
            .max()
            .unwrap_or(self.n_nodes.saturating_sub(1));
        (max_bit.saturating_add(65)).max(self.n_nodes + 64)
    }

    pub fn node(&self, id: usize) -> Option<&GraphNode> {
        self.nodes.get(id)
    }

    pub fn node_mut(&mut self, id: usize) -> Option<&mut GraphNode> {
        self.nodes.get_mut(id)
    }

    /// Compute childpat/parentpat (reachability) for compatibility checks.
    /// Is a non-source/sink node a "pass-through": exactly one non-source/sink
    /// parent p AND exactly one non-source/sink child c, AND contiguous with
    /// both (p.end == node.start AND node.end == c.start).
    ///
    /// Pass-through nodes add nothing path-structurally — they can be merged
    /// into their neighbors without changing the set of emittable paths.
    /// Conservative definition: skips any node with branching or junction gap.
    fn is_pass_through(&self, i: usize) -> bool {
        if i == self.source_id || i == self.sink_id {
            return false;
        }
        let node = &self.nodes[i];
        let parents_nonspecial: Vec<usize> = node
            .parents
            .ones()
            .filter(|&p| p != self.source_id && p != self.sink_id)
            .collect();
        let children_nonspecial: Vec<usize> = node
            .children
            .ones()
            .filter(|&c| c != self.source_id && c != self.sink_id)
            .collect();
        if parents_nonspecial.len() != 1 || children_nonspecial.len() != 1 {
            return false;
        }
        let p = parents_nonspecial[0];
        let c = children_nonspecial[0];
        // Must be contiguous with both neighbors (no intron gap).
        self.nodes[p].end == node.start && node.end == self.nodes[c].start
    }

    /// Remove zero-width graph nodes (nodes where `end <= start`, excluding
    /// source/sink). These are artifacts of junction events firing at shared
    /// donor/acceptor coords. Zero-width nodes are treated as TRANSPARENT:
    /// every path that went u→[zw…]→v in the original graph becomes u→v in
    /// the compacted graph (transitive closure through zw chains).
    ///
    /// When `include_pass_through` is true, also compacts non-zero-width
    /// pass-through nodes (exactly one contiguous parent and one contiguous
    /// child, non-source/sink). These carry no path-branching and can be
    /// merged without changing the set of emittable chains.
    ///
    /// Safe to call after `create_graph_inner` but BEFORE transfrag construction
    /// (pathpat bit layout depends on node IDs).
    ///
    /// Returns the number of nodes removed.
    pub fn compact_zero_width_nodes(&mut self) -> usize {
        self.compact_transparent_nodes(false)
    }

    /// Chain-compression: merge maximal linear chains of contiguous pass-through
    /// nodes into a single "head" node. Each chain has the shape:
    ///
    ///   N0 → N1 → N2 → ... → Nk  (all nodes contiguous, Ni.end == Ni+1.start)
    ///
    /// where N1..Nk-1 are pass-through (each has 1 non-source/sink parent and
    /// 1 non-source/sink child). N0 is the "anchor" head of the chain (no
    /// pass-through constraint). Nk-1 is the "tail" pass-through. Nk is where
    /// the chain ends (either branches out or reaches sink).
    ///
    /// Merge semantics (unlike edge-rewriting): the head node N0 has its `end`
    /// extended to Nk-1.end, and coverage + other scalar fields aggregated from
    /// the absorbed nodes. N1..Nk-1 are removed. Reads that previously mapped
    /// to N1..Nk-1 will map to the extended N0 because its coord range now
    /// covers theirs.
    ///
    /// CRITICAL: must be called BEFORE `map_reads` so read-to-node coords match.
    ///
    /// Returns the number of nodes merged away.
    pub fn compact_chain_pass_through(&mut self) -> usize {
        let n_old = self.nodes.len();
        if n_old <= 2 {
            return 0;
        }

        // Helper: check if node is pass-through.
        // STRICT definition:
        //   - exactly 1 parent, not source, not sink (no alt entry point)
        //   - exactly 1 child, not source, not sink (no alt exit point)
        //   - contiguous with both (no intron gap)
        //   - NOT a hardstart (not a TSS)
        //   - NOT a hardend (not a TTS)
        // Excluding nodes with source/sink neighbors prevents compacting
        // hardstart/hardend positions which carry biological significance.
        let is_pt = |i: usize, nodes: &[GraphNode], src: usize, snk: usize| -> bool {
            if i == src || i == snk {
                return false;
            }
            let node = &nodes[i];
            // Exclude nodes that are hardstart/hardend (TSS/TTS).
            if node.hardstart || node.hardend {
                return false;
            }
            // Exclude nodes carrying soft-boundary signal. longtrim_cov is set
            // on coverage-drop split boundaries (below the hardend threshold).
            // alt_tts_end is a terminal-read cluster at a junction donor. Both
            // are consumed by downstream filters that rely on the node boundary
            // being visible. Merging destroys that signal.
            if node.longtrim_cov > 0.0 || node.alt_tts_end {
                return false;
            }
            // ALL parents must be non-special, exactly 1 such.
            let all_parents: Vec<usize> = node.parents.ones().collect();
            if all_parents.len() != 1 {
                return false;
            }
            let p = all_parents[0];
            if p == src || p == snk {
                return false;
            }
            // ALL children must be non-special, exactly 1 such.
            let all_children: Vec<usize> = node.children.ones().collect();
            if all_children.len() != 1 {
                return false;
            }
            let c = all_children[0];
            if c == src || c == snk {
                return false;
            }
            // Do not merge across longtrim boundaries on the parent side:
            // parent.longtrim_cov > 0 means parent ends AT a coverage-drop
            // split. Extending over that drop would erase it.
            if nodes[p].longtrim_cov > 0.0 {
                return false;
            }
            // Contiguous: no intron gap either side.
            nodes[p].end == node.start && node.end == nodes[c].start
        };

        // 1. Identify all pass-through nodes.
        let is_passthrough: Vec<bool> = (0..n_old)
            .map(|i| is_pt(i, &self.nodes, self.source_id, self.sink_id))
            .collect();

        // 2. For each pass-through, walk backward to find the chain HEAD
        //    (the non-passthrough node that starts the chain). Heads must
        //    themselves be non-passthrough, but their child is pass-through.
        //
        //    We compute: head_of[i] = the chain head for each pass-through i.
        let mut head_of: Vec<Option<usize>> = vec![None; n_old];
        for i in 0..n_old {
            if !is_passthrough[i] {
                continue;
            }
            // Walk backward: parent of i (must exist since is_pt verified 1 parent)
            let parents_ns: Vec<usize> = self.nodes[i]
                .parents
                .ones()
                .filter(|&p| p != self.source_id && p != self.sink_id)
                .collect();
            if parents_ns.len() != 1 {
                continue;
            }
            let mut cur = parents_ns[0];
            let mut guard = 0usize;
            while is_passthrough[cur] && guard < n_old {
                let ps: Vec<usize> = self.nodes[cur]
                    .parents
                    .ones()
                    .filter(|&p| p != self.source_id && p != self.sink_id)
                    .collect();
                if ps.len() != 1 {
                    break;
                }
                cur = ps[0];
                guard += 1;
            }
            head_of[i] = Some(cur);
        }

        // 3. Build chain assignments: for each chain, head + tail + [pass-through members].
        //    Walk forward from each head through pass-through children.
        use std::collections::HashMap as StdMap;
        let mut chain_members: StdMap<usize, Vec<usize>> = StdMap::new(); // head → [pt members in order]
        for i in 0..n_old {
            if !is_passthrough[i] {
                continue;
            }
            if let Some(head) = head_of[i] {
                chain_members.entry(head).or_default().push(i);
            }
        }

        if chain_members.is_empty() {
            return 0;
        }

        // 4. For each chain, extend head's coords and transfer coverage.
        let mut removed_ids: std::collections::HashSet<usize> = Default::default();
        for (&head, members) in chain_members.iter() {
            // Sort members by genomic start (they should form a linear chain).
            let mut pts: Vec<usize> = members.clone();
            pts.sort_by_key(|&i| self.nodes[i].start);

            // Extend head.end to last pass-through's end.
            let last_end = pts.last().map(|&i| self.nodes[i].end).unwrap_or(0);
            if last_end > self.nodes[head].end {
                // Transfer coverage/abundance from absorbed nodes to head.
                let mut added_cov = 0.0f64;
                let mut added_longcov = 0.0f64;
                let mut added_nodecov = 0.0f64;
                for &pt_id in &pts {
                    added_cov += self.nodes[pt_id].coverage;
                    added_longcov += self.nodes[pt_id].longcov;
                    added_nodecov += self.nodes[pt_id].nodecov;
                }
                self.nodes[head].end = last_end;
                self.nodes[head].coverage += added_cov;
                self.nodes[head].longcov += added_longcov;
                self.nodes[head].nodecov += added_nodecov;

                // Propagate hardend if the last pass-through had it.
                if let Some(&last_pt) = pts.last() {
                    if self.nodes[last_pt].hardend {
                        crate::bump_hs!("graph.rs:752:hardend");
                        self.nodes[head].hardend = true;
                    }
                }

                // Mark pass-through members for removal.
                for &pt_id in &pts {
                    removed_ids.insert(pt_id);
                }
            }
        }

        if removed_ids.is_empty() {
            return 0;
        }

        // 5. Rebuild: for each chain, head → [head's original children minus first pt] + [tail's children].
        //    Simpler: collect all edges from the OLD graph. For each edge (u, v):
        //    - if u and v are both kept (not in removed_ids): keep as-is
        //    - if u is kept, v is removed: trace v's chain-tail's children, add (u, tail_child)
        //    - if u is removed: skip (will be added from head's iteration)
        //
        //    But we need to remap: the chain's head gets the edges FROM pts (head → C in the graph).

        // Easier: enumerate all non-removed nodes; for each, rebuild children.
        //   - Filter out children that are in removed_ids (those are chain members).
        //   - Instead: if child is in removed_ids, replace with child's chain's final child (what comes after the chain).
        let mut final_child_of = std::collections::HashMap::<usize, Option<usize>>::new();
        for (head, members) in chain_members.iter() {
            let mut pts: Vec<usize> = members.clone();
            pts.sort_by_key(|&i| self.nodes[i].start);
            if let Some(&last_pt) = pts.last() {
                let tail_children_ns: Vec<usize> = self.nodes[last_pt]
                    .children
                    .ones()
                    .filter(|&c| c != self.source_id && c != self.sink_id)
                    .collect();
                // The chain's successor (non-special) is the final child
                let succ = tail_children_ns.first().copied();
                // Also map all intermediate pts to this successor so any edge
                // into them gets redirected (though that shouldn't happen in linear chain).
                for &pt in &pts {
                    final_child_of.insert(pt, succ);
                }
                let _ = head;
            }
        }

        // 6. Rewrite edges.
        // Collect new edges based on surviving nodes.
        let mut new_edges: Vec<(usize, usize)> = Vec::new();
        for u in 0..n_old {
            if removed_ids.contains(&u) {
                continue;
            }
            for v in self.nodes[u].children.ones() {
                if removed_ids.contains(&v) {
                    // v is chain member; redirect to chain's final successor
                    if let Some(Some(succ)) = final_child_of.get(&v) {
                        if !removed_ids.contains(succ) {
                            new_edges.push((u, *succ));
                        }
                    }
                    // If no successor (chain ends at sink), drop edge
                } else {
                    new_edges.push((u, v));
                }
            }
        }

        // Also: chain HEADS need edges to their final successor (bypassing absorbed pts).
        for (head, members) in chain_members.iter() {
            if removed_ids.contains(head) {
                continue;
            }
            let mut pts: Vec<usize> = members.clone();
            pts.sort_by_key(|&i| self.nodes[i].start);
            if let Some(&last_pt) = pts.last() {
                let tail_children: Vec<usize> = self.nodes[last_pt].children.ones().collect();
                for c in tail_children {
                    if !removed_ids.contains(&c) {
                        new_edges.push((*head, c));
                    }
                }
            }
        }
        new_edges.sort_unstable();
        new_edges.dedup();

        // 7. Compact node array: build new IDs for surviving nodes.
        let mut new_id: Vec<Option<usize>> = vec![None; n_old];
        let mut new_nodes: Vec<GraphNode> = Vec::with_capacity(n_old);
        for (old_id, node) in self.nodes.iter().enumerate() {
            if removed_ids.contains(&old_id) {
                continue;
            }
            new_id[old_id] = Some(new_nodes.len());
            let mut cloned = node.clone();
            cloned.node_id = new_nodes.len();
            cloned.parents.clear();
            cloned.children.clear();
            cloned.childpat = None;
            cloned.parentpat = None;
            new_nodes.push(cloned);
        }
        let removed = removed_ids.len();

        // 8. Apply.
        self.nodes = new_nodes;
        self.n_nodes = self.nodes.len();
        self.source_id = new_id
            .get(self.source_id)
            .copied()
            .flatten()
            .unwrap_or(0);
        self.sink_id = new_id
            .get(self.sink_id)
            .copied()
            .flatten()
            .unwrap_or(self.n_nodes.saturating_sub(1));
        self.gpos.clear();
        self.edgeno = 0;
        self.next_edge_id = self.n_nodes;

        // 9. Rebuild edges via new IDs.
        for (u, v) in new_edges {
            if let (Some(nu), Some(nv)) = (new_id[u], new_id[v]) {
                if nu != nv {
                    self.add_edge(nu, nv);
                }
            }
        }

        // 10. Recompute reachability.
        self.compute_reachability();

        removed
    }

    /// Extended version: optionally include pass-through non-zw nodes.
    pub fn compact_transparent_nodes(&mut self, include_pass_through: bool) -> usize {
        let n_old = self.nodes.len();
        // 1. Identify transparent nodes:
        //    - zero-width (start == end), excluding source/sink
        //    - pass-through (1 contig parent, 1 contig child, non-source/sink)
        //      — only if include_pass_through is true.
        let is_zw: Vec<bool> = (0..n_old)
            .map(|i| {
                i != self.source_id
                    && i != self.sink_id
                    && self.nodes[i].end <= self.nodes[i].start
            })
            .collect();
        // Recompute including pass-through if requested.
        let is_zw: Vec<bool> = if include_pass_through {
            (0..n_old)
                .map(|i| is_zw[i] || self.is_pass_through(i))
                .collect()
        } else {
            is_zw
        };
        if !is_zw.iter().any(|&b| b) {
            return 0;
        }

        // 2. Build new IDs (compact, zero-widths dropped).
        let mut new_id: Vec<Option<usize>> = vec![None; n_old];
        let mut new_nodes: Vec<GraphNode> = Vec::with_capacity(n_old);
        for (old_id, node) in self.nodes.iter().enumerate() {
            if is_zw[old_id] {
                continue;
            }
            new_id[old_id] = Some(new_nodes.len());
            let mut cloned = node.clone();
            cloned.node_id = new_nodes.len();
            cloned.parents.clear();
            cloned.children.clear();
            cloned.childpat = None;
            cloned.parentpat = None;
            new_nodes.push(cloned);
        }
        let removed = n_old - new_nodes.len();

        // 3. Transitive closure: for each non-zw node u, find all non-zw v such
        // that there's a path u → [zw*] → v. This captures edges through zero-
        // width chains. DFS through zw-only intermediate nodes.
        //
        // We walk forward from each non-zw u through any zw children, collecting
        // all non-zw descendants reachable via zw-only paths.
        let mut edges: Vec<(usize, usize)> = Vec::new();
        for u in 0..n_old {
            if is_zw[u] {
                continue;
            }
            let Some(new_u) = new_id[u] else { continue };
            // BFS from u's children. Skip zw nodes but continue through their children.
            let mut stack: Vec<usize> = self.nodes[u].children.ones().collect();
            let mut visited: HashSet<usize> = Default::default();
            while let Some(v) = stack.pop() {
                if !visited.insert(v) {
                    continue;
                }
                if v == u {
                    continue;
                }
                if is_zw[v] {
                    // Transit through zw node: enqueue its children.
                    for c in self.nodes[v].children.ones() {
                        if !visited.contains(&c) {
                            stack.push(c);
                        }
                    }
                } else if let Some(new_v) = new_id[v] {
                    if new_u != new_v {
                        edges.push((new_u, new_v));
                    }
                    // Don't walk past non-zw — the edge u→v is enough; downstream
                    // edges originate from v's own iteration.
                }
            }
        }
        edges.sort_unstable();
        edges.dedup();

        // Resolve helper (used only for source/sink remap).
        let resolve_sink_source = |old: usize| -> Option<usize> {
            if old < new_id.len() {
                new_id[old]
            } else {
                None
            }
        };

        // 4. Apply new node list, remap source/sink.
        self.nodes = new_nodes;
        self.n_nodes = self.nodes.len();
        self.source_id = resolve_sink_source(self.source_id).unwrap_or(0);
        self.sink_id = resolve_sink_source(self.sink_id).unwrap_or(self.n_nodes.saturating_sub(1));
        self.gpos.clear();
        self.edgeno = 0;
        self.next_edge_id = self.n_nodes;

        // 5. Rebuild edges.
        for (u, v) in edges {
            self.add_edge(u, v);
        }

        // 6. Recompute reachability.
        self.compute_reachability();

        removed
    }

    pub fn compute_reachability(&mut self) {
        let n = self.n_nodes;
        let reachable = self.source_reachable_mask();
        let children: Vec<Vec<usize>> = self
            .nodes
            .iter()
            .map(|node| node.children.ones().collect())
            .collect();
        let parents: Vec<Vec<usize>> = self
            .nodes
            .iter()
            .map(|node| node.parents.ones().collect())
            .collect();

        fn fill_childpat(
            i: usize,
            graph: &Graph,
            reachable: &crate::bitset::SmallBitset,
            children: &[Vec<usize>],
            state: &mut [u8],
            memo: &mut [HashSet<usize>],
        ) {
            if state[i] == 2 {
                return;
            }
            if state[i] == 1 {
                return;
            }
            state[i] = 1;
            let mut to_add: HashSet<usize> = Default::default();
            for &c in &children[i] {
                if !reachable.contains(c) {
                    continue;
                }
                fill_childpat(c, graph, reachable, children, state, memo);
                to_add.insert(c);
                if let Some(eid) = graph.edge_bit_index(i, c) {
                    to_add.insert(eid);
                }
                to_add.extend(memo[c].iter().copied());
            }
            memo[i] = to_add;
            state[i] = 2;
        }

        fn fill_parentpat(
            i: usize,
            graph: &Graph,
            reachable: &crate::bitset::SmallBitset,
            parents: &[Vec<usize>],
            state: &mut [u8],
            memo: &mut [HashSet<usize>],
        ) {
            if state[i] == 2 {
                return;
            }
            if state[i] == 1 {
                return;
            }
            state[i] = 1;
            let mut to_add: HashSet<usize> = Default::default();
            for &p in &parents[i] {
                if !reachable.contains(p) {
                    continue;
                }
                fill_parentpat(p, graph, reachable, parents, state, memo);
                to_add.insert(p);
                if let Some(eid) = graph.edge_bit_index(p, i) {
                    to_add.insert(eid);
                }
                to_add.extend(memo[p].iter().copied());
            }
            memo[i] = to_add;
            state[i] = 2;
        }

        let mut child_state = vec![0u8; n];
        let mut parent_state = vec![0u8; n];
        let mut child_memo: Vec<HashSet<usize>> = (0..n).map(|_| Default::default()).collect();
        let mut parent_memo: Vec<HashSet<usize>> = (0..n).map(|_| Default::default()).collect();

        for i in 0..n {
            if !reachable.contains(i) {
                self.nodes[i].childpat = None;
                self.nodes[i].parentpat = None;
                continue;
            }
            fill_childpat(
                i,
                self,
                &reachable,
                &children,
                &mut child_state,
                &mut child_memo,
            );
            fill_parentpat(
                i,
                self,
                &reachable,
                &parents,
                &mut parent_state,
                &mut parent_memo,
            );
        }

        for i in 0..n {
            if reachable.contains(i) {
                let cm = std::mem::take(&mut child_memo[i]);
                let pm = std::mem::take(&mut parent_memo[i]);
                let mut cbs = crate::bitset::SmallBitset::empty();
                for id in cm {
                    cbs.insert_grow(id);
                }
                let mut pbs = crate::bitset::SmallBitset::empty();
                for id in pm {
                    pbs.insert_grow(id);
                }
                self.nodes[i].childpat = Some(cbs);
                self.nodes[i].parentpat = Some(pbs);
            }
        }
    }

    pub fn can_reach(&self, from_id: usize, to_id: usize) -> bool {
        self.nodes
            .get(from_id)
            .and_then(|n| n.childpat.as_ref())
            .map(|pat| pat.contains(to_id))
            .unwrap_or(false)
    }

    /// Set pattern edge bits for consecutive nodes in path (pattern edges).
    pub fn set_pattern_edges_for_path(&self, pattern: &mut GBitVec, path: &[usize]) {
        for i in 0..path.len().saturating_sub(1) {
            let (a, b) = (path[i], path[i + 1]);
            let key = (a.min(b), a.max(b));
            if let Some(&eid) = self.gpos.get(&key) {
                pattern.set_bit(eid);
            }
        }
    }

    /// Rebuild edge-bit IDs densely above node-id range.
    /// This mirrors `lastgpos` growth semantics where edge bits are assigned
    /// after graph shape stabilizes and remain disjoint from node bits.
    pub fn reindex_edge_bits_dense(&mut self) {
        let mut keys: Vec<(usize, usize)> = self.gpos.keys().copied().collect();
        keys.sort_unstable();
        let mut next = self.n_nodes;
        let mut new_map: HashMap<(usize, usize), usize> =
            HashMap::with_capacity_and_hasher(keys.len(), Default::default());
        for key in keys {
            new_map.insert(key, next);
            next = next.saturating_add(1);
        }
        self.gpos = new_map;
        self.edgeno = self.gpos.len();
        self.next_edge_id = next;
    }

    /// futuretr sink-proximity suppression analogue:
    /// return true if a downstream path from `from` reaches sink within `anchor` genomic distance.
    pub fn has_sink_within_anchor_downstream(&self, from: usize, anchor: u64) -> bool {
        let sink = self.sink_id;
        let Some(start_node) = self.nodes.get(from) else {
            return false;
        };
        let mut q: std::collections::VecDeque<(usize, u64, u64)> =
            std::collections::VecDeque::new();
        let mut seen: HashSet<usize> = Default::default();
        q.push_back((from, 0, start_node.end));
        seen.insert(from);
        while let Some((nid, dist, prev_end)) = q.pop_front() {
            let Some(node) = self.nodes.get(nid) else {
                continue;
            };
            for c in node.children.ones().collect::<Vec<_>>() {
                if c == sink {
                    return true;
                }
                let Some(cn) = self.nodes.get(c) else {
                    continue;
                };
                let gap = cn.start.saturating_sub(prev_end).saturating_sub(1);
                let add = gap + cn.length();
                let nd = dist.saturating_add(add);
                if nd > anchor {
                    continue;
                }
                if seen.insert(c) {
                    q.push_back((c, nd, cn.end));
                }
            }
        }
        false
    }

    /// traverse_dfs helper:
    /// collect deferred source->node and node->sink transfrags without mutating the graph.
    ///
    /// the original algorithm adds these terminal transfrags during later traversal, after read->transfrag
    /// mapping has already happened. If Rust materializes the source/sink edges here, long-read
    /// mapping sees synthetic source parents too early and can trim starts against them.
    pub fn synthesize_terminal_transfrags(
        &mut self,
        _mode: AssemblyMode,
        trthr: f64,
    ) -> Vec<GraphTransfrag> {
        let mut synth = Vec::new();
        let source_id = self.source_id;
        let sink_id = self.sink_id;
        // traverse_dfs: only source-reachable nodes are traversed and get
        // source/sink terminal links.
        let reachable = self.source_reachable_mask();

        for nid in 1..sink_id {
            let only_source_parent = {
                let p = &self.nodes[nid].parents;
                p.count_ones() == 1 && p.contains(source_id)
            };
            if !reachable.contains(nid) && !only_source_parent {
                continue;
            }
            if only_source_parent {
                let mut tf = GraphTransfrag::new(vec![source_id, nid], self.pattern_size());
                tf.abundance = trthr;
                tf.read_count = trthr;
                tf.longread = true;
                synth.push(tf);
            }

            let tail_needs_sink_tf = {
                let children = &self.nodes[nid].children;
                // futuretr:
                // a node->sink synthetic transfrag can exist even when the node is already
                // terminally connected to sink. The later futuretr materialization stage
                // decides whether to suppress that link based on downstream sink proximity.
                // We therefore need the deferred sink transfrag for:
                // - orphaned nodes
                // - nodes whose only child is sink
                children.is_empty() || (children.count_ones() == 1 && children.contains(sink_id))
            };
            if tail_needs_sink_tf {
                let mut tf = GraphTransfrag::new(vec![nid, sink_id], self.pattern_size());
                tf.abundance = trthr;
                tf.read_count = trthr;
                tf.longread = true;
                synth.push(tf);
            }
        }

        synth
    }

    #[inline]
    fn source_reachable_mask(&self) -> crate::bitset::SmallBitset {
        let mut visited = crate::bitset::SmallBitset::with_capacity(self.nodes.len().min(64));
        if self.source_id >= self.n_nodes {
            return visited;
        }
        let mut stack = vec![self.source_id];
        while let Some(nid) = stack.pop() {
            if visited.contains(nid) {
                continue;
            }
            visited.insert_grow(nid);
            for child in self.nodes[nid].children.ones().collect::<Vec<_>>() {
                if !visited.contains(child) {
                    stack.push(child);
                }
            }
        }
        visited
    }
}
