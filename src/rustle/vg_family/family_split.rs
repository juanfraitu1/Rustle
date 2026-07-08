//! De-novo family DECOMPOSITION — the `bench/denovo_family_split.py` core.
//!
//! Connected components of the POA-homology edge graph transitively close a SUPERFAMILY (e.g. KRAB-ZNF)
//! into one blob, because genes that merely share a tandem DOMAIN pairwise-confirm. A real recent-duplicate
//! family is a DENSE, mutually-homologous subgraph; an over-merge is a SPARSE chain bridged by lower
//! `core_recip` domain edges. We keep the same POA edges but split large components by **weighted
//! modularity** (uses the `core_recip` edge weight; no arbitrary similarity cutoff), then FLAG sparse
//! homology webs (`size >= web_min_size & density < web_max_density`) rather than claim them as families.
//!
//! Community detection is a **hand-rolled, deterministic** weighted-modularity Louvain (the Python uses
//! `networkx.community.louvain_communities`, seeded). Byte-matching networkx is infeasible (its greedy
//! order + RNG); the Python pipeline established the partition is robust (resolution 1→8 holds 96–96.7%
//! concordance), so the faithful target is the same modularity DEFINITION with deterministic tie-breaks
//! (sorted node/community order, strict-improvement moves), not a byte-identical partition.

use std::collections::{BTreeMap, BTreeSet};

/// Only decompose connected components this size or larger (small families left intact).
pub const MIN_DECOMP: usize = 6;
/// Modularity resolution (`gamma`); 1.0 = standard.
pub const RESOLUTION: f64 = 1.0;
/// A final community of `>= WEB_MIN_SIZE` nodes ...
pub const WEB_MIN_SIZE: usize = 10;
/// ... and `< WEB_MAX_DENSITY` density is flagged a homology web (a sparse domain-sharing over-merge),
/// not a family — and Web families are EXCLUDED from copy-assignment (denovo_pipeline.rs).
///
/// Set to 0.30 to ALIGN with the validated DNA cDNA-homology manifest bar (`make_dna_family_manifest.py`
/// drops `overmerge_sparse` at `density < 0.30`). Measured on the real de-novo split (698 families): the
/// family-density distribution is BIMODAL — median 1.000 (cliques) with NOTHING legitimate in the
/// [0.15, 0.30) band; the only 4 large (n>=10) families there are multi-chromosome domain-sharing
/// over-merges (DSFAM0 = a 164-member ZNF spanning 19 chromosomes, etc.) that the old 0.15 bar let
/// through as "families." `WEB_MIN_SIZE` stays 10 (not the manifest's 4) ON PURPOSE: a SMALL sparse
/// group can be a real divergent family (e.g. a 7-copy single-chrom MAGEB at density 0.24), so only
/// LARGE-and-sparse is called a web. See loose end L2 in bench/LOOSE_ENDS_AUDIT.md.
pub const WEB_MAX_DENSITY: f64 = 0.30;

/// Tunable decomposition parameters (defaults mirror `denovo_family_split.py`).
#[derive(Clone, Copy, Debug)]
pub struct SplitParams {
    pub min_decomp: usize,
    pub resolution: f64,
    pub web_min_size: usize,
    pub web_max_density: f64,
}

impl Default for SplitParams {
    fn default() -> Self {
        SplitParams {
            min_decomp: MIN_DECOMP,
            resolution: RESOLUTION,
            web_min_size: WEB_MIN_SIZE,
            web_max_density: WEB_MAX_DENSITY,
        }
    }
}

/// A discrete recent-duplicate family vs a non-discretizable homology web.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum FamilyClass {
    Family,
    Web,
}

/// Structural diagnostics of a community in the edge graph.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct CommunityStats {
    pub n: usize,
    pub n_edges: usize,
    pub density: f64,
    pub avg_core_recip: f64,
    pub n_articulation: usize,
}

/// A final decomposed family with its structural diagnostics and class.
#[derive(Clone, Debug)]
pub struct SplitFamily {
    pub members: Vec<usize>,
    pub stats: CommunityStats,
    pub class: FamilyClass,
}

/// Connected components (each a sorted node list) of the edge graph, keeping components of `>= min_size`.
/// Nodes are the ids appearing in `edges`; isolated nodes (no edge) are not represented.
fn uf_find(parent: &mut [usize], mut x: usize) -> usize {
    while parent[x] != x {
        parent[x] = parent[parent[x]];
        x = parent[x];
    }
    x
}
fn uf_union(parent: &mut [usize], a: usize, b: usize) {
    let ra = uf_find(parent, a);
    let rb = uf_find(parent, b);
    if ra != rb {
        parent[ra] = rb;
    }
}

pub fn connected_components(edges: &[(usize, usize, f64)], min_size: usize) -> Vec<Vec<usize>> {
    let max_id = match edges.iter().map(|&(a, b, _)| a.max(b)).max() {
        Some(m) => m,
        None => return Vec::new(),
    };
    let mut parent: Vec<usize> = (0..=max_id).collect();
    let mut present = vec![false; max_id + 1];
    for &(a, b, _) in edges {
        present[a] = true;
        present[b] = true;
        uf_union(&mut parent, a, b);
    }
    let mut groups: BTreeMap<usize, Vec<usize>> = BTreeMap::new();
    for node in 0..=max_id {
        if present[node] {
            let r = uf_find(&mut parent, node);
            groups.entry(r).or_default().push(node);
        }
    }
    let mut comps: Vec<Vec<usize>> = groups
        .into_values()
        .filter(|c| c.len() >= min_size)
        .collect();
    for c in &mut comps {
        c.sort_unstable();
    }
    // deterministic: size desc, then smallest member.
    comps.sort_by(|a, b| b.len().cmp(&a.len()).then_with(|| a[0].cmp(&b[0])));
    comps
}

/// Hand-rolled deterministic weighted-modularity Louvain over a graph with local node ids `0..n`.
/// Returns the communities (each a sorted Vec of local node ids). `resolution` is the modularity gamma.
pub fn louvain_communities(n: usize, edges: &[(usize, usize, f64)], resolution: f64) -> Vec<Vec<usize>> {
    if n == 0 {
        return Vec::new();
    }
    let labels = louvain_labels(n, edges, resolution);
    let mut groups: BTreeMap<usize, Vec<usize>> = BTreeMap::new();
    for (node, &lab) in labels.iter().enumerate() {
        groups.entry(lab).or_default().push(node);
    }
    let mut comms: Vec<Vec<usize>> = groups.into_values().collect();
    for c in &mut comms {
        c.sort_unstable();
    }
    comms.sort_by(|a, b| a[0].cmp(&b[0])); // deterministic: by smallest member
    comms
}

/// Multi-level Louvain: run local-moving to convergence, aggregate communities into super-nodes, repeat
/// until a level makes no moves. Returns a community label (renumbered `0..k`) for each original node.
fn louvain_labels(n: usize, edges: &[(usize, usize, f64)], resolution: f64) -> Vec<usize> {
    let mut node2super: Vec<usize> = (0..n).collect();
    let mut cur_n = n;
    let mut cur_edges: Vec<(usize, usize, f64)> = edges.to_vec();
    loop {
        let (labels, moved) = louvain_one_level(cur_n, &cur_edges, resolution);
        for s in node2super.iter_mut() {
            *s = labels[*s];
        }
        let k = labels.iter().copied().max().map(|m| m + 1).unwrap_or(0);
        if !moved || k >= cur_n {
            break;
        }
        // aggregate: each community becomes a super-node; internal edges become self-loops.
        let mut agg: BTreeMap<(usize, usize), f64> = BTreeMap::new();
        for &(i, j, w) in &cur_edges {
            let (ci, cj) = (labels[i], labels[j]);
            *agg.entry((ci.min(cj), ci.max(cj))).or_insert(0.0) += w;
        }
        cur_edges = agg.into_iter().map(|((a, b), w)| (a, b, w)).collect();
        cur_n = k;
    }
    // renumber the final super-node ids to a compact 0..K.
    let mut remap: BTreeMap<usize, usize> = BTreeMap::new();
    let mut next = 0usize;
    node2super
        .iter()
        .map(|&s| {
            *remap.entry(s).or_insert_with(|| {
                let v = next;
                next += 1;
                v
            })
        })
        .collect()
}

/// One Louvain level: greedy local moving to a local modularity optimum. Returns `(labels 0..k, moved)`.
/// Weighted-modularity gain (python-louvain formulation): moving node `i` into community `C` gains
/// `remove_cost + w(i,C) - γ·Σ_tot(C)·k_i/(2m)`, where `remove_cost = -w(i,C_old) + γ·(Σ_tot(C_old)-k_i)·k_i/(2m)`.
/// Determinism: nodes scanned in id order; neighbour communities in sorted (BTreeMap) order; only STRICT
/// improvements move, so the lowest-id community wins ties and the node stays put on a non-positive best.
fn louvain_one_level(n: usize, edges: &[(usize, usize, f64)], resolution: f64) -> (Vec<usize>, bool) {
    let mut deg = vec![0.0f64; n];
    let mut neigh: Vec<Vec<(usize, f64)>> = vec![Vec::new(); n];
    let mut m = 0.0f64;
    for &(i, j, w) in edges {
        m += w;
        if i == j {
            deg[i] += 2.0 * w; // a self-loop counts twice toward degree
        } else {
            deg[i] += w;
            deg[j] += w;
            neigh[i].push((j, w));
            neigh[j].push((i, w));
        }
    }
    if m <= 0.0 {
        return ((0..n).collect(), false);
    }
    let two_m = 2.0 * m;
    let mut com: Vec<usize> = (0..n).collect();
    let mut com_deg: Vec<f64> = deg.clone(); // Σ_tot of each community (community id ∈ 0..n)
    let mut any_moved = false;
    let mut improved = true;
    while improved {
        improved = false;
        for node in 0..n {
            let c_old = com[node];
            let dctw = deg[node] / two_m;
            let mut ncw: BTreeMap<usize, f64> = BTreeMap::new();
            for &(nb, w) in &neigh[node] {
                *ncw.entry(com[nb]).or_insert(0.0) += w;
            }
            let dnc_old = *ncw.get(&c_old).unwrap_or(&0.0);
            let remove_cost = -dnc_old + resolution * (com_deg[c_old] - deg[node]) * dctw;
            com_deg[c_old] -= deg[node]; // isolate the node
            let mut best_com = c_old;
            let mut best_inc = 0.0f64;
            for (&c, &dnc) in &ncw {
                let inc = remove_cost + dnc - resolution * com_deg[c] * dctw;
                if inc > best_inc {
                    best_inc = inc;
                    best_com = c;
                }
            }
            com_deg[best_com] += deg[node];
            com[node] = best_com;
            if best_com != c_old {
                improved = true;
                any_moved = true;
            }
        }
    }
    let mut remap: BTreeMap<usize, usize> = BTreeMap::new();
    let mut next = 0usize;
    let labels: Vec<usize> = com
        .iter()
        .map(|&c| {
            *remap.entry(c).or_insert_with(|| {
                let v = next;
                next += 1;
                v
            })
        })
        .collect();
    (labels, any_moved)
}

/// Articulation (cut) points of an undirected graph with local node ids `0..n`. Returns sorted ids.
pub fn articulation_points(n: usize, edges: &[(usize, usize)]) -> Vec<usize> {
    let mut adj: Vec<Vec<usize>> = vec![Vec::new(); n];
    for &(a, b) in edges {
        if a != b {
            adj[a].push(b);
            adj[b].push(a);
        }
    }
    let mut visited = vec![false; n];
    let mut disc = vec![0usize; n];
    let mut low = vec![0usize; n];
    let mut is_ap = vec![false; n];
    let mut timer = 0usize;
    // iterative Tarjan (avoids recursion-depth blowups on large families).
    for start in 0..n {
        if visited[start] {
            continue;
        }
        let mut stack: Vec<(usize, isize, usize)> = vec![(start, -1, 0)]; // (node, parent, child cursor)
        let mut root_children = 0usize;
        while let Some(&(u, parent, ci)) = stack.last() {
            if ci == 0 {
                visited[u] = true;
                timer += 1;
                disc[u] = timer;
                low[u] = timer;
            }
            if ci < adj[u].len() {
                stack.last_mut().unwrap().2 += 1;
                let v = adj[u][ci];
                if v as isize == parent {
                    continue;
                }
                if !visited[v] {
                    if parent == -1 {
                        root_children += 1;
                    }
                    stack.push((v, u as isize, 0));
                } else {
                    low[u] = low[u].min(disc[v]);
                }
            } else {
                stack.pop();
                if let Some(&(p, pparent, _)) = stack.last() {
                    low[p] = low[p].min(low[u]);
                    if pparent != -1 && low[u] >= disc[p] {
                        is_ap[p] = true; // non-root cut vertex
                    }
                }
            }
        }
        if root_children > 1 {
            is_ap[start] = true;
        }
    }
    (0..n).filter(|&i| is_ap[i]).collect()
}

/// Structural diagnostics of the subgraph induced by `members` (global node ids) over `edges`.
pub(crate) fn community_stats(members: &[usize], edges: &[(usize, usize, f64)]) -> CommunityStats {
    let set: BTreeSet<usize> = members.iter().copied().collect();
    let n = set.len();
    let mut n_edges = 0usize;
    let mut wsum = 0.0f64;
    let mut internal: Vec<(usize, usize)> = Vec::new();
    for &(a, b, w) in edges {
        if a != b && set.contains(&a) && set.contains(&b) {
            n_edges += 1;
            wsum += w;
            internal.push((a, b));
        }
    }
    let density = if n > 1 {
        2.0 * n_edges as f64 / (n as f64 * (n as f64 - 1.0))
    } else {
        1.0
    };
    let avg_core_recip = if n_edges > 0 { wsum / n_edges as f64 } else { 0.0 };
    // articulation only for n > 2 (python `arts = ... if n > 2 else 0`); relabel members to 0..n.
    let n_articulation = if n > 2 {
        let idx: BTreeMap<usize, usize> = set.iter().enumerate().map(|(i, &node)| (node, i)).collect();
        let local: Vec<(usize, usize)> = internal.iter().map(|&(a, b)| (idx[&a], idx[&b])).collect();
        articulation_points(n, &local).len()
    } else {
        0
    };
    CommunityStats { n, n_edges, density, avg_core_recip, n_articulation }
}

/// Classify a community by size + density (web iff `n >= web_min_size && density < web_max_density`).
pub fn classify(n: usize, density: f64, p: &SplitParams) -> FamilyClass {
    if n >= p.web_min_size && density < p.web_max_density {
        FamilyClass::Web
    } else {
        FamilyClass::Family
    }
}

/// Decompose the POA-homology edge graph into final families: connected components, with components
/// `>= min_decomp` split by weighted modularity (communities of `>= 2` kept), each annotated with
/// structural diagnostics and a family/web class. Sorted by size desc, then smallest member.
pub fn decompose_families(edges: &[(usize, usize, f64)], p: &SplitParams) -> Vec<SplitFamily> {
    let comps = connected_components(edges, 2);
    let mut final_members: Vec<Vec<usize>> = Vec::new();
    for comp in comps {
        if comp.len() < p.min_decomp {
            final_members.push(comp);
            continue;
        }
        // decompose: relabel the component to local ids 0..k, Louvain, map communities (>= 2) back.
        let set: BTreeSet<usize> = comp.iter().copied().collect();
        let idx: BTreeMap<usize, usize> = comp.iter().enumerate().map(|(i, &nd)| (nd, i)).collect();
        let local_edges: Vec<(usize, usize, f64)> = edges
            .iter()
            .filter(|&&(a, b, _)| a != b && set.contains(&a) && set.contains(&b))
            .map(|&(a, b, w)| (idx[&a], idx[&b], w))
            .collect();
        for community in louvain_communities(comp.len(), &local_edges, p.resolution) {
            if community.len() >= 2 {
                final_members.push(community.iter().map(|&li| comp[li]).collect());
            }
        }
    }
    let mut out: Vec<SplitFamily> = final_members
        .into_iter()
        .map(|mut members| {
            members.sort_unstable();
            let stats = community_stats(&members, edges);
            let class = classify(stats.n, stats.density, p);
            SplitFamily { members, stats, class }
        })
        .collect();
    // deterministic output order: size desc, then smallest member.
    out.sort_by(|a, b| {
        b.members
            .len()
            .cmp(&a.members.len())
            .then_with(|| a.members[0].cmp(&b.members[0]))
    });
    out
}

/// Internal edge density of a node block over the induced subgraph: 2|E|/(|C|(|C|-1)); <=1 node = 1.0.
fn induced_density(members: &[usize], edges: &[(usize, usize, f64)]) -> f64 {
    let set: std::collections::HashSet<usize> = members.iter().copied().collect();
    let n = members.len();
    if n <= 1 {
        return 1.0;
    }
    let m = edges.iter().filter(|(a, b, _)| set.contains(a) && set.contains(b)).count();
    2.0 * m as f64 / (n as f64 * (n as f64 - 1.0))
}

/// Restrict edges to those with both endpoints in `members`, remapped to local indices 0..members.len().
fn induced_edges(members: &[usize], edges: &[(usize, usize, f64)]) -> (usize, Vec<(usize, usize, f64)>) {
    let idx: std::collections::HashMap<usize, usize> =
        members.iter().enumerate().map(|(i, &g)| (g, i)).collect();
    let local = edges
        .iter()
        .filter_map(|&(a, b, w)| match (idx.get(&a), idx.get(&b)) {
            (Some(&la), Some(&lb)) => Some((la, lb, w)),
            _ => None,
        })
        .collect();
    (members.len(), local)
}

/// Guaranteed-progress splitter: adaptive-resolution Louvain, else component split, else deterministic halving.
/// Returns global-index blocks; never a single block equal to the input when |members| > 2.
fn split_once(members: &[usize], edges: &[(usize, usize, f64)]) -> Vec<Vec<usize>> {
    let (n, local) = induced_edges(members, edges);
    for res in [1.0, 2.0, 4.0, 8.0] {
        let parts = louvain_communities(n, &local, res);
        if parts.len() >= 2 {
            return parts
                .into_iter()
                .map(|p| p.into_iter().map(|l| members[l]).collect())
                .collect();
        }
    }
    // connected-component fallback (also catches a disconnected block).
    let comps = crate::vg_family::read_conflict::conflict_families(
        n,
        &local.iter().map(|&(a, b, _)| (a, b, 1usize)).collect::<Vec<_>>(),
    );
    if comps.len() >= 2 {
        return comps
            .into_iter()
            .map(|c| c.into_iter().map(|l| members[l]).collect())
            .collect();
    }
    // deterministic halving.
    let h = members.len() / 2;
    vec![members[..h].to_vec(), members[h..].to_vec()]
}

/// gamma-quasi-clique partition: keep a block whole iff it is already a gamma-quasi-clique (or <=2 nodes),
/// else split (guaranteed-progress) and recurse. Blocks partition 0..n. Deterministic (Louvain is
/// deterministic here).
pub fn gamma_quasi_clique_partition(n: usize, edges: &[(usize, usize, f64)], gamma: f64) -> Vec<Vec<usize>> {
    // start from raw connected components, then refine each.
    let comps = crate::vg_family::read_conflict::conflict_families(
        n,
        &edges.iter().map(|&(a, b, _)| (a, b, 1usize)).collect::<Vec<_>>(),
    );
    let mut out = Vec::new();
    let mut stack: Vec<Vec<usize>> = comps;
    while let Some(block) = stack.pop() {
        if block.len() <= 2 || induced_density(&block, edges) >= gamma {
            out.push(block);
            continue;
        }
        let parts = split_once(&block, edges);
        if parts.len() == 1 && parts[0].len() == block.len() {
            out.push(block); // no-progress guard (shouldn't happen)
        } else {
            stack.extend(parts);
        }
    }
    out
}

#[cfg(test)]
mod tests {
    use super::*;

    fn we(a: usize, b: usize) -> (usize, usize, f64) {
        (a, b, 1.0)
    }

    // ---- connected_components ----

    #[test]
    fn cc_two_separate_edges_two_components() {
        let e = [we(0, 1), we(2, 3)];
        let mut c = connected_components(&e, 2);
        c.sort();
        assert_eq!(c, vec![vec![0, 1], vec![2, 3]]);
    }

    #[test]
    fn cc_transitive_one_component() {
        let e = [we(0, 1), we(1, 2)];
        assert_eq!(connected_components(&e, 2), vec![vec![0, 1, 2]]);
    }

    // ---- louvain_communities ----

    #[test]
    fn louvain_single_edge_one_community() {
        assert_eq!(louvain_communities(2, &[we(0, 1)], 1.0), vec![vec![0, 1]]);
    }

    #[test]
    fn louvain_single_clique_one_community() {
        let e = [we(0, 1), we(1, 2), we(0, 2)];
        assert_eq!(louvain_communities(3, &e, 1.0), vec![vec![0, 1, 2]]);
    }

    #[test]
    fn louvain_two_cliques_bridge_splits() {
        // two triangles {0,1,2} and {3,4,5}, strong internal weight, bridged by ONE weak edge.
        let e = [
            we(0, 1), we(1, 2), we(0, 2),
            we(3, 4), we(4, 5), we(3, 5),
            (0, 3, 0.05),
        ];
        let mut comms = louvain_communities(6, &e, 1.0);
        comms.sort();
        assert_eq!(comms, vec![vec![0, 1, 2], vec![3, 4, 5]]);
    }

    // ---- articulation_points ----

    #[test]
    fn artic_path_middle_is_cut() {
        assert_eq!(articulation_points(3, &[(0, 1), (1, 2)]), vec![1]);
    }

    #[test]
    fn artic_triangle_none() {
        assert_eq!(articulation_points(3, &[(0, 1), (1, 2), (0, 2)]), Vec::<usize>::new());
    }

    #[test]
    fn artic_path_of_four_two_cuts() {
        assert_eq!(articulation_points(4, &[(0, 1), (1, 2), (2, 3)]), vec![1, 2]);
    }

    // ---- community_stats ----

    #[test]
    fn stats_clique_density_one() {
        let e = [(0, 1, 0.5), (1, 2, 0.7), (0, 2, 0.9)];
        let s = community_stats(&[0, 1, 2], &e);
        assert_eq!(s.n, 3);
        assert_eq!(s.n_edges, 3);
        assert!((s.density - 1.0).abs() < 1e-9);
        assert!((s.avg_core_recip - (0.5 + 0.7 + 0.9) / 3.0).abs() < 1e-9);
        assert_eq!(s.n_articulation, 0);
    }

    #[test]
    fn stats_path_density_and_articulation() {
        // path 0-1-2-3: 3 edges, density = 2*3/(4*3) = 0.5, articulation points {1,2}.
        let e = [(0, 1, 0.2), (1, 2, 0.2), (2, 3, 0.2)];
        let s = community_stats(&[0, 1, 2, 3], &e);
        assert_eq!(s.n_edges, 3);
        assert!((s.density - 0.5).abs() < 1e-9);
        assert_eq!(s.n_articulation, 2);
    }

    #[test]
    fn stats_only_counts_internal_edges() {
        // members {0,1} but an edge 1-2 leaves the set -> not counted.
        let e = [(0, 1, 0.4), (1, 2, 0.9)];
        let s = community_stats(&[0, 1], &e);
        assert_eq!(s.n_edges, 1);
        assert!((s.avg_core_recip - 0.4).abs() < 1e-9);
    }

    // ---- classify ----

    #[test]
    fn classify_web_vs_family() {
        let p = SplitParams::default(); // web_max_density = 0.30, web_min_size = 10
        assert_eq!(classify(10, 0.10, &p), FamilyClass::Web); // size>=10 & sparse
        assert_eq!(classify(12, 0.167, &p), FamilyClass::Web); // the DSFAM0-class over-merge (164 ZNF/19 chr): was Family at the old 0.15 bar
        assert_eq!(classify(10, 0.29, &p), FamilyClass::Web); // just under the aligned 0.30 bar
        assert_eq!(classify(10, 0.40, &p), FamilyClass::Family); // dense -> family
        assert_eq!(classify(10, 1.00, &p), FamilyClass::Family); // clique
        assert_eq!(classify(5, 0.10, &p), FamilyClass::Family); // too SMALL -> family even if sparse (small divergent fam, e.g. MAGEB)
        assert_eq!(classify(7, 0.24, &p), FamilyClass::Family); // n<10 stays family (the MAGEB-class protection)
    }

    // ---- decompose_families ----

    #[test]
    fn decompose_small_component_kept_intact() {
        let e = [we(0, 1), we(1, 2), we(0, 2)]; // 3 < MIN_DECOMP
        let fams = decompose_families(&e, &SplitParams::default());
        assert_eq!(fams.len(), 1);
        assert_eq!(fams[0].members, vec![0, 1, 2]);
        assert_eq!(fams[0].class, FamilyClass::Family);
    }

    #[test]
    fn decompose_large_component_splits() {
        // two K4 cliques bridged by one weak edge; 8 nodes >= MIN_DECOMP -> split into two size-4 families.
        let mut e = vec![];
        for &(a, b) in &[(0, 1), (0, 2), (0, 3), (1, 2), (1, 3), (2, 3)] {
            e.push(we(a, b));
        }
        for &(a, b) in &[(4, 5), (4, 6), (4, 7), (5, 6), (5, 7), (6, 7)] {
            e.push(we(a, b));
        }
        e.push((0, 4, 0.05));
        let fams = decompose_families(&e, &SplitParams::default());
        assert_eq!(fams.len(), 2);
        let mut sizes: Vec<usize> = fams.iter().map(|f| f.members.len()).collect();
        sizes.sort();
        assert_eq!(sizes, vec![4, 4]);
        assert!(fams.iter().all(|f| f.class == FamilyClass::Family));
    }

    #[test]
    fn decompose_flags_sparse_web() {
        // a star: center 0 with 14 leaves. Louvain keeps it one community; n=15, density ~0.133 -> web.
        let e: Vec<(usize, usize, f64)> = (1..15).map(|leaf| we(0, leaf)).collect();
        let fams = decompose_families(&e, &SplitParams::default());
        assert_eq!(fams.len(), 1);
        assert_eq!(fams[0].members.len(), 15);
        assert!(fams[0].stats.density < 0.15);
        assert_eq!(fams[0].class, FamilyClass::Web);
    }

    #[test]
    fn louvain_resolution_controls_granularity() {
        // K6 clique. Low resolution keeps it whole; high resolution shatters it into singletons -- this
        // directly pins the `gamma` term in the modularity gain.
        let mut e = vec![];
        for a in 0..6 {
            for b in (a + 1)..6 {
                e.push(we(a, b));
            }
        }
        assert_eq!(louvain_communities(6, &e, 0.5).len(), 1, "low resolution keeps the clique whole");
        assert_eq!(louvain_communities(6, &e, 2.0).len(), 6, "high resolution shatters K6 into singletons");
    }

    #[test]
    fn louvain_multilevel_aggregation_merges_supernodes() {
        // Hierarchical: two DISCONNECTED groups, each = 4 triangles whose anchor nodes form a strong K4.
        // Level-1 local-moving forms the per-triangle communities; only the AGGREGATION level coalesces the
        // 4 triangle super-nodes of a group into one community. The correct partition is 2 groups of 12.
        let mut e: Vec<(usize, usize, f64)> = vec![];
        for g in 0..2 {
            let base = g * 12;
            let anchors = [base, base + 3, base + 6, base + 9];
            for t in 0..4 {
                let (x, y, z) = (base + 3 * t, base + 3 * t + 1, base + 3 * t + 2);
                e.push(we(x, y));
                e.push(we(y, z));
                e.push(we(x, z)); // triangle (weight 1.0)
            }
            for i in 0..anchors.len() {
                for j in (i + 1)..anchors.len() {
                    e.push((anchors[i], anchors[j], 5.0)); // strong inter-triangle K4 among anchors
                }
            }
        }
        let mut comms = louvain_communities(24, &e, 1.0);
        comms.sort();
        assert_eq!(comms.len(), 2, "two groups");
        assert!(comms.iter().all(|c| c.len() == 12), "each group is its 12 nodes: {comms:?}");
        assert_eq!(comms[0], (0..12).collect::<Vec<_>>());
        assert_eq!(comms[1], (12..24).collect::<Vec<_>>());
    }

    #[test]
    fn louvain_beats_weakest_edge_cut() {
        // The single WEAKEST edge (0-2, weight 0.1) is INSIDE community A; the correct split is across the
        // bridge (2-3, weight 0.3). A "cut the weakest edge" heuristic would wrongly split A; only a real
        // modularity optimiser recovers {0,1,2},{3,4,5}.
        let e = [
            (0, 1, 1.0), (1, 2, 1.0), (0, 2, 0.1), // community A: one weak intra edge
            (3, 4, 1.0), (4, 5, 1.0), (3, 5, 1.0), // community B: strong triangle
            (2, 3, 0.3),                           // bridge (stronger than the weakest intra edge)
        ];
        let mut comms = louvain_communities(6, &e, 1.0);
        comms.sort();
        assert_eq!(comms, vec![vec![0, 1, 2], vec![3, 4, 5]]);
    }

    #[test]
    fn decompose_drops_singleton_communities() {
        // K7 (strong, weight 10) plus a pendant node 7 weakly attached (0-7, weight 0.05). At resolution
        // 1.1 the clique survives but the pendant detaches into its OWN singleton community, which
        // decompose drops (`len >= 2`), so node 7 vanishes from all families.
        let mut e = vec![];
        for a in 0..7 {
            for b in (a + 1)..7 {
                e.push((a, b, 10.0));
            }
        }
        e.push((0, 7, 0.05));
        let p = SplitParams { resolution: 1.1, ..SplitParams::default() };
        let fams = decompose_families(&e, &p);
        assert_eq!(fams.len(), 1, "only the clique family survives");
        assert_eq!(fams[0].members, (0..7).collect::<Vec<_>>());
        assert!(!fams[0].members.contains(&7), "the dropped singleton node is gone");
    }

    // ---- gamma_quasi_clique_partition ----

    #[test]
    fn gamma_quasi_clique_keeps_array_whole_splits_repeat_chain() {
        // A 5-node dense clique (a tandem array) stays ONE block.
        let clique: Vec<(usize, usize, f64)> = {
            let mut e = Vec::new();
            for i in 0..5 { for j in (i + 1)..5 { e.push((i, j, 1.0)); } }
            e
        };
        let blocks = gamma_quasi_clique_partition(5, &clique, 0.20);
        assert_eq!(blocks.len(), 1, "a dense array is one gamma-quasi-clique");
        assert_eq!(blocks[0].len(), 5);

        // A LONG sparse bridge chain (12-node path: density 2*11/(12*11)=0.167 < gamma=0.20) must split.
        let chain: Vec<(usize, usize, f64)> = (0..11).map(|i| (i, i + 1, 1.0)).collect();
        let blocks = gamma_quasi_clique_partition(12, &chain, 0.20);
        assert!(blocks.len() >= 2, "a sparse repeat-bridge chain is split, got {:?}", blocks);
    }

    #[test]
    fn community_stats_relabels_noncontiguous_ids() {
        // members are NON-contiguous global ids forming a path 10-3-27-5; exercises the global->local
        // relabel the integration layer relies on. Articulation points are the two middle nodes.
        let e = [(10, 3, 0.2), (3, 27, 0.2), (27, 5, 0.2)];
        let s = community_stats(&[10, 3, 27, 5], &e);
        assert_eq!(s.n, 4);
        assert_eq!(s.n_edges, 3);
        assert!((s.density - 0.5).abs() < 1e-9);
        assert_eq!(s.n_articulation, 2);
    }
}
