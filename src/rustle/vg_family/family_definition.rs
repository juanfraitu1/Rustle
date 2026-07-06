//! O1 multi-copy family-definition predicate `distinct_loci`.
//!
//! **Terminology.** This module uses "locus" in the **physical-span** sense:
//! a `(chrom, start, end)` genomic interval. It is NOT the gene-locus / splice-
//! junction community used to collapse isoforms (that is
//! `family_detect::collapse_loci`). See `docs/VG_FAMILY_TERMS.md` for the
//! canonical vocabulary.
//!
//! Faithful Rust port of the Python reference
//! `bench/genome_family_def.py::distinct_loci` (LOCUS_OVERLAP = 0.50). This is
//! the first Python->Rust migration of the family-definition layer: it counts
//! the DISTINCT PHYSICAL LOCI among a family's member gene spans, collapsing
//! redundant/nested annotations at the same physical locus (e.g. MAGEA9 ==
//! a nested LOC copy) before the ">= 2 distinct loci" multi-copy test.
//!
//! The DNA multi-copy family predicate `|distinct_loci(C)| >= 2` is defined in
//! `bench/DNA_FAMILY_DEFINITION_FORMAL.md` (the >=2-loci CERTIFICATE); the
//! collapse step is the `LOCUS_OVERLAP` reciprocal-span rule documented there.
//!
//! Semantics preserved EXACTLY from the Python reference:
//!   (a) every member counts as its own locus unless collapsed (the Python
//!       `uf.find(i)` singleton guarantee) -- a member that overlaps no other
//!       is one locus;
//!   (b) two members collapse to one locus iff SAME contig AND overlap > 0 AND
//!       overlap >= 0.50 * len_a AND overlap >= 0.50 * len_b (RECIPROCAL);
//!   (c) collapse is transitive (union-find), so a chain A~B~C where A and C do
//!       not directly overlap enough still collapses to one locus;
//!   (d) members on different contigs never collapse.
//!
//! The count is invariant to member order and to union order, so the per-contig
//! `sort_by_key(start)` and the naive-union union-find below reproduce the
//! Python `len(uf.groups())` byte-for-byte. Verified against a golden fixture
//! generated from the shipped Python build (see `#[test] distinct_loci_parity`).

use std::collections::{BTreeMap, BTreeSet, HashMap, HashSet, VecDeque};

/// Two member spans on the SAME contig are the SAME physical locus iff their
/// spans reciprocally overlap by >= this fraction of BOTH lengths. Mirrors
/// `genome_family_def.LOCUS_OVERLAP`.
pub const LOCUS_OVERLAP: f64 = 0.50;

/// Number of distinct physical loci among `members`, each `(chrom, start, end)`
/// with 0-based half-open coords (matching the Python `genes` dict).
///
/// Members on the same contig whose spans reciprocally overlap by at least
/// `LOCUS_OVERLAP` of both lengths are collapsed (union-find) into one locus;
/// the return value is the count of distinct groups. See the module docs for
/// the exact semantics preserved from `genome_family_def.py::distinct_loci`.
pub fn distinct_loci(members: &[(String, i64, i64)]) -> usize {
    let mut uf = UnionFind::new(members.len());

    // Group member indices by contig (cross-contig pairs never collapse).
    let mut by_contig: HashMap<&str, Vec<usize>> = HashMap::new();
    for (i, (chrom, _, _)) in members.iter().enumerate() {
        by_contig.entry(chrom.as_str()).or_default().push(i);
    }

    for idxs in by_contig.values_mut() {
        // Start-sorted so the `sb >= ea` break is a valid overlap-prune, exactly
        // as in the Python `idxs.sort(key=lambda i: genes[i]["start"])`.
        idxs.sort_by_key(|&i| members[i].1);
        for a in 0..idxs.len() {
            let ia = idxs[a];
            let (sa, ea) = (members[ia].1, members[ia].2);
            for &ib in &idxs[a + 1..] {
                let (sb, eb) = (members[ib].1, members[ib].2);
                if sb >= ea {
                    break; // start-sorted: nothing further on this contig can overlap A
                }
                let ov = ea.min(eb) - sa.max(sb);
                // RECIPROCAL >= LOCUS_OVERLAP of both lengths (f64 to match Python).
                if ov > 0
                    && ov as f64 >= LOCUS_OVERLAP * (ea - sa) as f64
                    && ov as f64 >= LOCUS_OVERLAP * (eb - sb) as f64
                {
                    uf.union(ia, ib);
                }
            }
        }
    }
    uf.num_groups()
}

/// Minimal union-find (path-compression + naive union), a faithful port of the
/// Python `genome_family_def.UF`. Naive union (`parent[root_a] = root_b`) is
/// kept deliberately: the returned GROUP COUNT is union-order invariant, so the
/// exact merge direction is irrelevant to `distinct_loci`.
struct UnionFind {
    parent: Vec<usize>,
}

impl UnionFind {
    fn new(n: usize) -> Self {
        UnionFind {
            parent: (0..n).collect(),
        }
    }

    fn find(&mut self, x: usize) -> usize {
        let mut root = x;
        while self.parent[root] != root {
            root = self.parent[root];
        }
        // Path compression (as in the Python `while self.p[x] != r` walk).
        let mut cur = x;
        while self.parent[cur] != root {
            let next = self.parent[cur];
            self.parent[cur] = root;
            cur = next;
        }
        root
    }

    fn union(&mut self, a: usize, b: usize) {
        let ra = self.find(a);
        let rb = self.find(b);
        if ra != rb {
            self.parent[ra] = rb;
        }
    }

    /// Number of distinct groups among all `n` members (each member is a node,
    /// matching the Python `uf.find(i)` singleton guarantee).
    fn num_groups(&mut self) -> usize {
        let n = self.parent.len();
        let mut roots = std::collections::HashSet::with_capacity(n);
        for i in 0..n {
            let r = self.find(i);
            roots.insert(r);
        }
        roots.len()
    }
}

// ===========================================================================
// REFINEMENT operator R  (Python port of genome_family_def.refine_families /
// _refine_component / _induced_density / _split_once).
//
// R splits each raw single-linkage component (connected component of the kept
// homology graph) into gamma-quasi-clique COMMUNITIES via a density-gated
// recursion, then keeps the blocks that pass the >= 2-distinct-loci multi-copy
// predicate (reusing `distinct_loci` above). This is the operator the shipped
// RNA family definition (`bench/family_rna_refine.py`) applies at gamma=0.20.
//
// OPTION B -- DETERMINISTIC SPLITTER (no networkx, no RNG). The Python reference
// uses a seeded networkx Louvain splitter; here `split_once` is a fully
// deterministic, RNG-free splitter:
//   1. connected components of the induced subgraph  (deterministic; matches the
//      Python `nx.connected_components` branch);
//   2. otherwise deterministic GREEDY-MODULARITY agglomeration (Clauset-Newman-
//      Moore, with a sorted community-pair tie-break) -- NATURAL communities,
//      preferred over crude halving where it applies;
//   3. otherwise a sanctioned deterministic HALVING by sorted node order (the
//      guaranteed-progress last resort; see DNA_FAMILY_DEFINITION_FORMAL.md).
// Every branch yields >= 2 non-empty parts for n > 2, so the recursion provably
// terminates with every leaf a gamma-quasi-clique or |C| <= 2.
//
// PORT CONTRACT (verified against a golden fixture in
// `#[test] refine_families_parity_and_certificate`):
//   * on the 1932 raw components that never invoke the splitter (n <= 2 OR
//     induced density >= gamma) the refined family SET is BYTE-IDENTICAL to
//     Python's (the deterministic gate path is a faithful port);
//   * on the 1 component that DOES invoke the splitter (a 213-node blob, density
//     < gamma) every refined family is a VALID witness -- a gamma-quasi-clique
//     (or size <= 2) with >= 2 distinct loci -- but the specific family SET may
//     differ from Python's (the exact-max-gamma-quasi-clique partition is
//     NP-hard, so both are certificate-passing witnesses; the deterministic
//     Rust community split is one such witness, the seeded Louvain another).
// ===========================================================================

/// gamma cohesion parameter: a refined family must be a gamma-quasi-clique
/// (internal edge density >= GAMMA). Mirrors `genome_family_def.GAMMA`.
pub const GAMMA: f64 = 0.20;

/// A certificate-passing family whose dominant-biotype fraction is < PURITY is
/// DESCRIPTIVELY flagged repeat_mosaic and sorted to the END of the catalog (it
/// is reported, not filtered). Mirrors `genome_family_def.PURITY`.
pub const PURITY: f64 = 0.50;

/// A family member: a node's gene metadata, mirroring one entry of the Python
/// `genes` dict. `chrom/start/end` drive `distinct_loci`; `name/biotype` drive
/// the `refine_families` catalog sort (`biotype_purity` + name tie-break).
#[derive(Clone, Debug)]
pub struct Gene {
    pub chrom: String,
    pub start: i64,
    pub end: i64,
    pub name: String,
    pub biotype: String,
}

/// Internal edge density `rho_in(C) = 2|E(C)| / (|C|(|C|-1))` of the subgraph
/// induced by `nodes` under the (symmetric, within-family) adjacency `adj`. A
/// singleton/empty subset is 1.0 (trivially cohesive). Faithful port of
/// `genome_family_def._induced_density`.
pub fn induced_density(nodes: &BTreeSet<usize>, adj: &HashMap<usize, BTreeSet<usize>>) -> f64 {
    let n = nodes.len();
    if n <= 1 {
        return 1.0;
    }
    let mut e = 0usize;
    for &u in nodes {
        if let Some(nb) = adj.get(&u) {
            for &w in nb {
                // count each induced edge once (u < w), matching the Python `sub` build
                if w > u && nodes.contains(&w) {
                    e += 1;
                }
            }
        }
    }
    2.0 * e as f64 / (n as f64 * (n as f64 - 1.0))
}

/// Connected components of the subgraph induced by `nodes` (deterministic BFS in
/// ascending node order). Matches the Python `nx.connected_components(sub)` split
/// branch; the ascending `BTreeSet` iteration makes the component order stable.
fn connected_components(
    nodes: &BTreeSet<usize>,
    adj: &HashMap<usize, BTreeSet<usize>>,
) -> Vec<BTreeSet<usize>> {
    let mut visited: HashSet<usize> = HashSet::new();
    let mut comps: Vec<BTreeSet<usize>> = Vec::new();
    for &start in nodes {
        if visited.contains(&start) {
            continue;
        }
        let mut comp: BTreeSet<usize> = BTreeSet::new();
        let mut queue: VecDeque<usize> = VecDeque::new();
        queue.push_back(start);
        visited.insert(start);
        while let Some(u) = queue.pop_front() {
            comp.insert(u);
            if let Some(nb) = adj.get(&u) {
                for &w in nb {
                    if nodes.contains(&w) && !visited.contains(&w) {
                        visited.insert(w);
                        queue.push_back(w);
                    }
                }
            }
        }
        comps.push(comp);
    }
    comps
}

/// Deterministic greedy-modularity (Clauset-Newman-Moore) community partition of
/// the subgraph induced by `nodes`. Starts with every node its own community and
/// repeatedly merges the community pair with the greatest positive modularity
/// gain `dQ = k_xy/m - D_x*D_y/(2 m^2)` (D = sum of induced degrees), stopping
/// when no merge is strictly positive. DETERMINISTIC: community id = min node in
/// the community; ties in `dQ` break to the lexicographically smallest community
/// pair (sorted keys). No RNG.
fn greedy_modularity(
    nodes: &BTreeSet<usize>,
    adj: &HashMap<usize, BTreeSet<usize>>,
) -> Vec<BTreeSet<usize>> {
    let verts: Vec<usize> = nodes.iter().copied().collect(); // ascending -> pos = index
    let k = verts.len();
    let pos_of: HashMap<usize, usize> = verts.iter().enumerate().map(|(p, &v)| (v, p)).collect();

    // induced edge list (pi < pj) + induced degree per position
    let mut edges: Vec<(usize, usize)> = Vec::new();
    let mut degree = vec![0i64; k];
    for (pi, &v) in verts.iter().enumerate() {
        if let Some(nb) = adj.get(&v) {
            for &w in nb {
                if let Some(&pj) = pos_of.get(&w) {
                    degree[pi] += 1; // each incident induced edge counted once per endpoint
                    if pi < pj {
                        edges.push((pi, pj));
                    }
                }
            }
        }
    }
    let m = edges.len();
    if m == 0 {
        // no induced edges -> each node its own community (unreachable for a
        // connected n>2 subgraph, but keeps the function total)
        return verts
            .iter()
            .map(|&v| {
                let mut s = BTreeSet::new();
                s.insert(v);
                s
            })
            .collect();
    }
    let mf = m as f64;
    let mut comm: Vec<usize> = (0..k).collect(); // community id = min pos, invariant-maintained

    loop {
        // tot[c] = sum of induced degrees in community c
        let mut tot: HashMap<usize, i64> = HashMap::new();
        for p in 0..k {
            *tot.entry(comm[p]).or_default() += degree[p];
        }
        // cross[(ca,cb)] = # induced edges between distinct communities ca < cb
        let mut cross: HashMap<(usize, usize), usize> = HashMap::new();
        for &(pi, pj) in &edges {
            let (ca, cb) = (comm[pi], comm[pj]);
            if ca != cb {
                let key = if ca < cb { (ca, cb) } else { (cb, ca) };
                *cross.entry(key).or_default() += 1;
            }
        }
        if cross.is_empty() {
            break;
        }
        // deterministic best: max dQ, ties -> smallest (ca,cb)
        let mut keys: Vec<(usize, usize)> = cross.keys().copied().collect();
        keys.sort_unstable();
        let mut best: Option<((usize, usize), f64)> = None;
        for key in keys {
            let (ca, cb) = key;
            let ke = cross[&key] as f64;
            let dq = ke / mf - (tot[&ca] as f64 * tot[&cb] as f64) / (2.0 * mf * mf);
            match best {
                None => best = Some((key, dq)),
                Some((_, bq)) => {
                    if dq > bq {
                        best = Some((key, dq));
                    }
                }
            }
        }
        let (key, bq) = best.unwrap();
        if bq <= 0.0 {
            break;
        }
        let (keep, drop) = key; // key.0 < key.1 -> merge drop into the smaller id keep
        for c in comm.iter_mut() {
            if *c == drop {
                *c = keep;
            }
        }
    }

    let mut groups: BTreeMap<usize, BTreeSet<usize>> = BTreeMap::new();
    for p in 0..k {
        groups.entry(comm[p]).or_default().insert(verts[p]);
    }
    groups.into_values().collect()
}

/// Deterministic halving of `nodes` by sorted node order into two non-empty
/// blocks (guaranteed-progress last resort; the Python `_split_once` deterministic
/// halving). For n > 2 both halves are strictly smaller than n.
fn halve(nodes: &BTreeSet<usize>) -> Vec<BTreeSet<usize>> {
    let verts: Vec<usize> = nodes.iter().copied().collect();
    let mid = verts.len() / 2;
    let a: BTreeSet<usize> = verts[..mid].iter().copied().collect();
    let b: BTreeSet<usize> = verts[mid..].iter().copied().collect();
    vec![a, b]
}

/// GUARANTEED-PROGRESS deterministic splitter (OPTION B; no networkx / no RNG).
/// Never returns a single block equal to the input for n > 2, so the recursion
/// in `refine_component` provably terminates. See the module-level R note.
fn split_once(
    nodes: &BTreeSet<usize>,
    adj: &HashMap<usize, BTreeSet<usize>>,
) -> Vec<BTreeSet<usize>> {
    // 1. disconnected induced subgraph -> split by connected components
    let comps = connected_components(nodes, adj);
    if comps.len() >= 2 {
        return comps;
    }
    // 2. connected -> deterministic greedy-modularity natural communities
    let communities = greedy_modularity(nodes, adj);
    if communities.len() >= 2 {
        return communities;
    }
    // 3. modularity-stuck -> sanctioned deterministic halving (guaranteed >= 2 parts)
    halve(nodes)
}

/// R applied to ONE raw component's node set -> a PARTITION into gamma-quasi-clique
/// (or size <= 2) blocks. Density-gated recursion: keep a block whole iff it is
/// already a gamma-quasi-clique (or |C| <= 2), else split (deterministically) and
/// recurse. Faithful port of `genome_family_def._refine_component`.
pub fn refine_component(
    nodes: &BTreeSet<usize>,
    adj: &HashMap<usize, BTreeSet<usize>>,
    gamma: f64,
) -> Vec<BTreeSet<usize>> {
    let n = nodes.len();
    if n <= 2 || induced_density(nodes, adj) >= gamma {
        return vec![nodes.clone()];
    }
    let mut out: Vec<BTreeSet<usize>> = Vec::new();
    for part in split_once(nodes, adj) {
        if part.len() == n {
            // no-progress guard (cannot fire with the guaranteed splitter; matches Python)
            return vec![nodes.clone()];
        }
        out.extend(refine_component(&part, adj, gamma));
    }
    out
}

/// Distinct physical loci among member indices (reuses `distinct_loci` over the
/// members' `(chrom, start, end)` coords). Port of `distinct_loci(members, genes)`.
fn distinct_loci_members(members: &[usize], genes: &[Gene]) -> usize {
    let coords: Vec<(String, i64, i64)> = members
        .iter()
        .map(|&i| (genes[i].chrom.clone(), genes[i].start, genes[i].end))
        .collect();
    distinct_loci(&coords)
}

/// Dominant-biotype fraction of a family: `max_b |{i in C : biotype(i)=b}| / |C|`.
/// Faithful port of `genome_family_def.biotype_purity`.
pub fn biotype_purity(members: &[usize], genes: &[Gene]) -> f64 {
    if members.is_empty() {
        return 1.0;
    }
    let mut counts: HashMap<&str, usize> = HashMap::new();
    for &i in members {
        *counts.entry(genes[i].biotype.as_str()).or_default() += 1;
    }
    let max = counts.values().copied().max().unwrap_or(0);
    max as f64 / members.len() as f64
}

/// R over a whole raw catalog: split each raw family (connected component) into
/// gamma-quasi-clique communities, keep the blocks passing the >= 2-distinct-loci
/// multi-copy predicate, then sort `(biotype_purity < PURITY, -len, name)` so
/// repeat_mosaic hubs sort to the end and genuine (largest-first) families lead.
/// Faithful port of `genome_family_def.refine_families`.
///
/// `fams` are raw components as node-index lists; `all_edges` are the homology
/// edges (index pairs); `genes` is indexed by node id (assign indices in node-id
/// sort order to reproduce the exact Python catalog order). No RNG; deterministic.
pub fn refine_families(
    fams: &[Vec<usize>],
    all_edges: &[(usize, usize)],
    genes: &[Gene],
    gamma: f64,
) -> Vec<Vec<usize>> {
    // node -> raw-family index
    let mut fam_of: HashMap<usize, usize> = HashMap::new();
    for (k, m) in fams.iter().enumerate() {
        for &i in m {
            fam_of.insert(i, k);
        }
    }
    // adjacency restricted to same-raw-family pairs (matches the Python `adj` build)
    let mut adj: HashMap<usize, BTreeSet<usize>> = HashMap::new();
    for &(u, v) in all_edges {
        if let (Some(&fu), Some(&fv)) = (fam_of.get(&u), fam_of.get(&v)) {
            if fu == fv {
                adj.entry(u).or_default().insert(v);
                adj.entry(v).or_default().insert(u);
            }
        }
    }

    let mut refined: Vec<Vec<usize>> = Vec::new();
    for m in fams {
        let nodes: BTreeSet<usize> = m.iter().copied().collect();
        for block in refine_component(&nodes, &adj, gamma) {
            let members: Vec<usize> = block.into_iter().collect(); // BTreeSet -> sorted (== Python sorted(block))
            if distinct_loci_members(&members, genes) >= 2 {
                refined.push(members);
            }
        }
    }

    // repeat_mosaic hubs (purity < PURITY) sort to the END; genuine families
    // first, largest-first, name tie-break. Stable sort (deterministic).
    refined.sort_by(|a, b| {
        let pa = (biotype_purity(a, genes) < PURITY) as u8;
        let pb = (biotype_purity(b, genes) < PURITY) as u8;
        pa.cmp(&pb)
            .then_with(|| b.len().cmp(&a.len())) // -len => longest first
            .then_with(|| genes[a[0]].name.cmp(&genes[b[0]].name))
    });
    refined
}

#[cfg(test)]
mod tests {
    use super::*;

    fn m(chrom: &str, start: i64, end: i64) -> (String, i64, i64) {
        (chrom.to_string(), start, end)
    }

    #[test]
    fn empty_and_singleton() {
        assert_eq!(distinct_loci(&[]), 0);
        assert_eq!(distinct_loci(&[m("c1", 100, 200)]), 1);
    }

    #[test]
    fn reciprocal_and_boundary() {
        // exact 50% overlap of both -> collapse
        assert_eq!(distinct_loci(&[m("c1", 100, 200), m("c1", 150, 250)]), 1);
        // just under 50% -> distinct
        assert_eq!(distinct_loci(&[m("c1", 100, 200), m("c1", 151, 251)]), 2);
        // nested but reciprocal FAILS on the big span (10% of it) -> distinct
        assert_eq!(distinct_loci(&[m("c1", 100, 1100), m("c1", 100, 200)]), 2);
        // different contig, identical coords -> never collapse
        assert_eq!(distinct_loci(&[m("c1", 100, 200), m("c2", 100, 200)]), 2);
        // zero-length member never overlaps (ov > 0 fails) -> distinct
        assert_eq!(distinct_loci(&[m("c1", 100, 100), m("c1", 90, 110)]), 2);
        // transitive chain A~B~C (A,C don't directly overlap enough) -> one locus
        assert_eq!(
            distinct_loci(&[m("c1", 0, 100), m("c1", 30, 130), m("c1", 60, 160)]),
            1
        );
    }

    /// Byte-parity against the golden fixture emitted by
    /// `bench/gen_distinct_loci_fixture.py` from the shipped Python build. Every
    /// row's Rust `distinct_loci` MUST equal the Python golden count.
    #[test]
    fn distinct_loci_parity() {
        let fixture = include_str!("testdata/distinct_loci_fixture.tsv");
        let mut n_rows = 0usize;
        let mut n_collapse = 0usize;
        let mut n_xchrom = 0usize;
        for (i, line) in fixture.lines().enumerate() {
            if line.is_empty() || line.starts_with('#') {
                continue;
            }
            let cols: Vec<&str> = line.split('\t').collect();
            assert_eq!(
                cols.len(),
                4,
                "fixture line {}: expected 4 tab-columns, got {}",
                i + 1,
                cols.len()
            );
            let label = cols[0];
            let n_members: usize = cols[1].parse().expect("n_members int");
            let golden: usize = cols[2].parse().expect("golden int");
            let members: Vec<(String, i64, i64)> = cols[3]
                .split(';')
                .map(|tok| {
                    let mut p = tok.split('|');
                    let chrom = p.next().expect("chrom").to_string();
                    let start: i64 = p.next().expect("start").parse().expect("start int");
                    let end: i64 = p.next().expect("end").parse().expect("end int");
                    assert!(p.next().is_none(), "member token has extra fields: {tok}");
                    (chrom, start, end)
                })
                .collect();
            assert_eq!(
                members.len(),
                n_members,
                "fixture line {} ({}): declared n_members {} != parsed {}",
                i + 1,
                label,
                n_members,
                members.len()
            );

            let got = distinct_loci(&members);
            assert_eq!(
                got, golden,
                "PARITY MISMATCH family '{}' (fixture line {}): rust distinct_loci={} \
                 != python golden={} (n_members={})",
                label,
                i + 1,
                got,
                golden,
                n_members
            );

            n_rows += 1;
            if golden < n_members {
                n_collapse += 1;
            }
            let chroms: std::collections::HashSet<&str> =
                members.iter().map(|(c, _, _)| c.as_str()).collect();
            if chroms.len() >= 2 {
                n_xchrom += 1;
            }
        }
        // Guard the fixture stays representative (diverse, exercises collapse).
        assert!(n_rows >= 200, "fixture too small: {n_rows} rows (< 200)");
        assert!(
            n_collapse >= 10,
            "fixture under-exercises same-locus collapse: only {n_collapse} collapse rows"
        );
        assert!(
            n_xchrom >= 10,
            "fixture under-exercises cross-chrom families: only {n_xchrom} rows"
        );
        eprintln!(
            "distinct_loci parity: {n_rows} rows OK ({n_collapse} collapse, {n_xchrom} cross-chrom)"
        );
    }

    // ======================================================================
    // Refinement operator R (refine_families / refine_component / splitter).
    // ======================================================================

    fn g(chrom: &str, start: i64, end: i64) -> Gene {
        Gene {
            chrom: chrom.to_string(),
            start,
            end,
            name: format!("{chrom}:{start}"),
            biotype: "protein_coding".to_string(),
        }
    }

    fn build_adj(edges: &[(usize, usize)]) -> HashMap<usize, BTreeSet<usize>> {
        let mut adj: HashMap<usize, BTreeSet<usize>> = HashMap::new();
        for &(u, v) in edges {
            adj.entry(u).or_default().insert(v);
            adj.entry(v).or_default().insert(u);
        }
        adj
    }

    /// Run R on one component (nodes 0..n) + adjacency, returning the SET of
    /// refined member-sets that pass the >= 2-distinct-loci filter.
    fn refine_component_set(
        n: usize,
        adj: &HashMap<usize, BTreeSet<usize>>,
        genes: &[Gene],
        gamma: f64,
    ) -> HashSet<Vec<usize>> {
        let nodes: BTreeSet<usize> = (0..n).collect();
        let mut out: HashSet<Vec<usize>> = HashSet::new();
        for block in refine_component(&nodes, adj, gamma) {
            let mem: Vec<usize> = block.into_iter().collect();
            if distinct_loci_members(&mem, genes) >= 2 {
                out.insert(mem);
            }
        }
        out
    }

    #[test]
    fn induced_density_basics() {
        let nodes: BTreeSet<usize> = (0..3).collect();
        // triangle -> density 1.0
        assert!((induced_density(&nodes, &build_adj(&[(0, 1), (1, 2), (0, 2)])) - 1.0).abs() < 1e-12);
        // path 0-1-2 -> 2E/(n(n-1)) = 2*2/6
        assert!(
            (induced_density(&nodes, &build_adj(&[(0, 1), (1, 2)])) - (2.0 * 2.0 / 6.0)).abs()
                < 1e-12
        );
        // singleton -> 1.0
        let one: BTreeSet<usize> = [0].into_iter().collect();
        assert!((induced_density(&one, &build_adj(&[])) - 1.0).abs() < 1e-12);
    }

    #[test]
    fn split_disconnected_by_components() {
        // two triangles, no cross edge -> connected-components branch, 2 parts
        let adj = build_adj(&[(0, 1), (1, 2), (0, 2), (3, 4), (4, 5), (3, 5)]);
        let nodes: BTreeSet<usize> = (0..6).collect();
        let parts = split_once(&nodes, &adj);
        assert_eq!(parts.len(), 2);
        assert_eq!(parts.iter().map(|p| p.len()).collect::<Vec<_>>(), vec![3, 3]);
        // deterministic
        assert_eq!(split_once(&nodes, &adj), parts);
    }

    #[test]
    fn greedy_modularity_barbell() {
        // two K4 cliques joined by a single bridge -> >= 2 natural communities
        let mut edges = Vec::new();
        for a in 0..4 {
            for b in (a + 1)..4 {
                edges.push((a, b));
            }
        }
        for a in 4..8 {
            for b in (a + 1)..8 {
                edges.push((a, b));
            }
        }
        edges.push((3, 4)); // bridge
        let adj = build_adj(&edges);
        let nodes: BTreeSet<usize> = (0..8).collect();
        let parts = greedy_modularity(&nodes, &adj);
        assert!(
            parts.len() >= 2,
            "greedy modularity should split the barbell, got {} communities",
            parts.len()
        );
        // the two K4 are recovered as the natural communities
        assert!(parts.contains(&(0..4).collect::<BTreeSet<usize>>()));
        assert!(parts.contains(&(4..8).collect::<BTreeSet<usize>>()));
        // deterministic across reruns
        assert_eq!(greedy_modularity(&nodes, &adj), parts);
    }

    #[test]
    fn halve_guarantees_progress() {
        // guaranteed-progress last resort: two non-empty halves, each < n for n>2
        let nodes: BTreeSet<usize> = (0..5).collect();
        let parts = halve(&nodes);
        assert_eq!(parts.len(), 2);
        assert!(parts.iter().all(|p| !p.is_empty() && p.len() < 5));
        assert_eq!(parts[0].len() + parts[1].len(), 5);
    }

    #[test]
    fn refine_component_leaves_are_cohesive() {
        // 5 triangles in a chain (overall density 0.181 < gamma) -> must recurse;
        // every leaf a gamma-quasi-clique (or size<=2) and the blocks PARTITION.
        let mut edges = Vec::new();
        for t in 0..5 {
            let b = t * 3;
            edges.push((b, b + 1));
            edges.push((b + 1, b + 2));
            edges.push((b, b + 2));
        }
        for t in 0..4 {
            edges.push((t * 3 + 2, t * 3 + 3)); // bridges
        }
        let adj = build_adj(&edges);
        let nodes: BTreeSet<usize> = (0..15).collect();
        let gamma = 0.20;
        assert!(
            induced_density(&nodes, &adj) < gamma,
            "test graph must be below gamma to force a split"
        );
        let blocks = refine_component(&nodes, &adj, gamma);
        // partition check
        let mut union: BTreeSet<usize> = BTreeSet::new();
        for blk in &blocks {
            for &x in blk {
                assert!(union.insert(x), "blocks overlap on {x}");
            }
        }
        assert_eq!(union, nodes, "blocks do not partition the component");
        // certificate: every leaf cohesive
        for blk in &blocks {
            assert!(
                blk.len() <= 2 || induced_density(blk, &adj) >= gamma,
                "leaf not gamma-cohesive: size {} density {}",
                blk.len(),
                induced_density(blk, &adj)
            );
        }
        // determinism
        assert_eq!(refine_component(&nodes, &adj, gamma), blocks);
    }

    #[test]
    fn refine_families_small_and_sort() {
        // two disjoint cohesive families on distinct loci; largest-first order
        let genes = vec![
            g("c1", 0, 100),
            g("c1", 200, 300),
            g("c1", 400, 500), // fam A (3 loci)
            g("c2", 0, 100),
            g("c2", 200, 300), // fam B (2 loci)
        ];
        let fams = vec![vec![0, 1, 2], vec![3, 4]];
        let edges = vec![(0, 1), (1, 2), (0, 2), (3, 4)];
        let refined = refine_families(&fams, &edges, &genes, 0.20);
        assert_eq!(refined.len(), 2);
        assert_eq!(refined[0], vec![0, 1, 2]); // 3-locus family first
        assert_eq!(refined[1], vec![3, 4]);
        // determinism
        assert_eq!(refine_families(&fams, &edges, &genes, 0.20), refined);
    }

    #[test]
    fn refine_families_repeat_mosaic_sorts_last() {
        let gb = |chrom: &str, s: i64, e: i64, bt: &str| Gene {
            chrom: chrom.to_string(),
            start: s,
            end: e,
            name: format!("{chrom}:{s}"),
            biotype: bt.to_string(),
        };
        // A: pure 2-locus family; B: larger 4-locus mosaic (purity 0.25 < PURITY)
        let genes = vec![
            gb("c1", 0, 100, "protein_coding"),
            gb("c1", 200, 300, "protein_coding"), // A: 0,1
            gb("c2", 0, 100, "protein_coding"),
            gb("c2", 200, 300, "lncRNA"),
            gb("c2", 400, 500, "snRNA"),
            gb("c2", 600, 700, "pseudogene"), // B: 2,3,4,5
        ];
        let fams = vec![vec![0, 1], vec![2, 3, 4, 5]];
        // B is a K4 (density 1.0 >= gamma)
        let edges = vec![
            (0, 1),
            (2, 3),
            (3, 4),
            (4, 5),
            (2, 4),
            (2, 5),
            (3, 5),
        ];
        let refined = refine_families(&fams, &edges, &genes, 0.20);
        assert_eq!(refined.len(), 2);
        assert!((biotype_purity(&[2, 3, 4, 5], &genes) - 0.25).abs() < 1e-12);
        // mosaic hub sorts LAST even though it is the larger family
        assert_eq!(refined[0], vec![0, 1], "pure family should lead");
        assert_eq!(refined[1], vec![2, 3, 4, 5], "repeat_mosaic hub sorts last");
    }

    fn parse_gene_row(nd: &serde_json::Value) -> Gene {
        let a = nd.as_array().expect("node row");
        Gene {
            chrom: a[0].as_str().expect("chrom").to_string(),
            start: a[1].as_i64().expect("start"),
            end: a[2].as_i64().expect("end"),
            name: a[3].as_str().expect("name").to_string(),
            biotype: a[4].as_str().expect("biotype").to_string(),
        }
    }

    /// GOLDEN-FIXTURE parity + certificate against the shipped RNA build
    /// (bench/gen_refine_families_fixture.py): SET byte-parity on the 1932 clean
    /// (non-splitter) components, and gamma-quasi-clique + >= 2-loci CERTIFICATE
    /// validity on the 1 splitter blob (213 nodes). Deterministic (rerun-stable).
    #[test]
    fn refine_families_parity_and_certificate() {
        let raw = include_str!("testdata/refine_families_fixture.json");
        let fx: serde_json::Value = serde_json::from_str(raw).expect("parse fixture json");
        let gamma = fx["gamma"].as_f64().expect("gamma");
        assert!(
            (gamma - GAMMA).abs() < 1e-12,
            "fixture gamma {gamma} != GAMMA {GAMMA}"
        );
        let comps = fx["components"].as_array().expect("components array");
        assert_eq!(comps.len(), fx["n_components"].as_u64().expect("n_components") as usize);

        let mut n_split = 0usize;
        let mut n_clean = 0usize;
        let mut clean_fam_total = 0usize;
        let mut blob_nodes = 0usize;
        let mut rust_blob_fams = 0usize;
        let mut py_blob_fams = 0usize;

        for comp in comps {
            let was_split = comp["was_split"].as_bool().expect("was_split");
            let genes: Vec<Gene> = comp["nodes"]
                .as_array()
                .expect("nodes")
                .iter()
                .map(parse_gene_row)
                .collect();
            let n = genes.len();
            let mut adj: HashMap<usize, BTreeSet<usize>> = HashMap::new();
            for e in comp["edges"].as_array().expect("edges") {
                let ep = e.as_array().expect("edge");
                let i = ep[0].as_u64().unwrap() as usize;
                let j = ep[1].as_u64().unwrap() as usize;
                adj.entry(i).or_default().insert(j);
                adj.entry(j).or_default().insert(i);
            }
            let py_set: HashSet<Vec<usize>> = comp["python_families"]
                .as_array()
                .expect("python_families")
                .iter()
                .map(|f| {
                    f.as_array()
                        .unwrap()
                        .iter()
                        .map(|x| x.as_u64().unwrap() as usize)
                        .collect::<Vec<usize>>()
                })
                .collect();

            let rust_set = refine_component_set(n, &adj, &genes, gamma);

            if was_split {
                n_split += 1;
                blob_nodes = n;
                rust_blob_fams = rust_set.len();
                py_blob_fams = py_set.len();
                // CERTIFICATE: gamma-quasi-clique (or size<=2) AND >= 2 distinct loci
                for fam in &rust_set {
                    assert!(
                        distinct_loci_members(fam, &genes) >= 2,
                        "blob family fails the >= 2-loci predicate: {fam:?}"
                    );
                    let fset: BTreeSet<usize> = fam.iter().copied().collect();
                    let dens = induced_density(&fset, &adj);
                    assert!(
                        fam.len() <= 2 || dens >= gamma,
                        "blob family not gamma-cohesive: size {} density {}",
                        fam.len(),
                        dens
                    );
                }
                // determinism on the splitter blob
                assert_eq!(
                    rust_set,
                    refine_component_set(n, &adj, &genes, gamma),
                    "blob refinement not deterministic across reruns"
                );
            } else {
                // BYTE-PARITY (set-based) on the deterministic clean components
                assert_eq!(
                    rust_set, py_set,
                    "CLEAN-COMPONENT PARITY MISMATCH (n={n}): rust {} vs python {} families",
                    rust_set.len(),
                    py_set.len()
                );
                n_clean += 1;
                clean_fam_total += rust_set.len();
            }
        }

        assert_eq!(n_split, 1, "expected exactly 1 splitter component, got {n_split}");
        assert_eq!(blob_nodes, 213, "expected the 213-node blob, got {blob_nodes}");
        assert!(n_clean >= 1900, "too few clean components: {n_clean}");
        assert!(clean_fam_total >= 500, "too few clean families: {clean_fam_total}");
        eprintln!(
            "refine_families parity+certificate: {n_clean} clean comps SET-parity EXACT \
             ({clean_fam_total} families); blob(n={blob_nodes}) certificate OK -> rust \
             {rust_blob_fams} families vs python {py_blob_fams} (valid non-unique witnesses)"
        );
    }

    /// DETERMINISM (criterion 3): the FULL refine_families (adjacency build +
    /// filter + sort + splitter) over the whole reconstructed graph reruns
    /// byte-identical, and every family passes the >= 2-loci predicate.
    #[test]
    fn refine_families_full_deterministic() {
        let raw = include_str!("testdata/refine_families_fixture.json");
        let fx: serde_json::Value = serde_json::from_str(raw).expect("parse fixture json");
        let gamma = fx["gamma"].as_f64().unwrap();

        let mut genes: Vec<Gene> = Vec::new();
        let mut fams: Vec<Vec<usize>> = Vec::new();
        let mut edges: Vec<(usize, usize)> = Vec::new();
        for comp in fx["components"].as_array().unwrap() {
            let base = genes.len();
            let nodes_json = comp["nodes"].as_array().unwrap();
            for nd in nodes_json {
                genes.push(parse_gene_row(nd));
            }
            fams.push((base..base + nodes_json.len()).collect());
            for e in comp["edges"].as_array().unwrap() {
                let ep = e.as_array().unwrap();
                edges.push((
                    base + ep[0].as_u64().unwrap() as usize,
                    base + ep[1].as_u64().unwrap() as usize,
                ));
            }
        }

        let r1 = refine_families(&fams, &edges, &genes, gamma);
        let r2 = refine_families(&fams, &edges, &genes, gamma);
        assert_eq!(r1, r2, "refine_families is not deterministic across reruns");
        for fam in &r1 {
            assert!(
                distinct_loci_members(fam, &genes) >= 2,
                "refined family fails the >= 2-loci predicate"
            );
        }
        assert!(
            r1.len() >= 584 && r1.len() <= 650,
            "refined family count {} out of the expected band [584,650]",
            r1.len()
        );
        eprintln!(
            "refine_families full-graph deterministic: {} families (byte-identical across reruns)",
            r1.len()
        );
    }
}
