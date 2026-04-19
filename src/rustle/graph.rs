//! Splice graph: nodes = exon segments, edges = collinear + junction edges.

use crate::bitset::NodeSet;
use crate::bitvec::GBitVec;
use crate::coord::len_half_open;
use crate::types::{AssemblyMode, DetHashMap as HashMap, DetHashSet as HashSet};

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
        }
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
