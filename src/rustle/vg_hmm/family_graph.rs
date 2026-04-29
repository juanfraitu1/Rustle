//! Family graph: union of per-copy splice graphs across a gene family.
//!
//! Nodes are exon-equivalence classes (groups of exons across copies that
//! we treat as the same exon). Edges are junctions observed in any copy.
//! Copy-specific exons become singleton nodes — that is what makes a
//! family graph different from the per-copy splice graphs it is built
//! from.

use std::ops::Index;

#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, PartialOrd, Ord)]
pub struct NodeIdx(pub usize);

pub type CopyId = usize;

#[derive(Debug, Clone)]
pub struct ExonClass {
    pub idx: NodeIdx,
    pub chrom: String,
    pub span: (u64, u64),
    pub strand: char,
    /// (CopyId, raw genomic exon sequence on the transcript strand).
    pub per_copy_sequences: Vec<(CopyId, Vec<u8>)>,
    /// True iff exactly one copy contributes — i.e. a bubble branch.
    pub copy_specific: bool,
}

#[derive(Debug, Clone, Copy)]
pub struct JunctionEdge {
    pub from: NodeIdx,
    pub to: NodeIdx,
    pub family_support: u32,
    pub strand: char,
}

#[derive(Debug, Clone)]
pub struct FamilyGraph {
    pub family_id: usize,
    pub nodes: Vec<ExonClass>,
    pub edges: Vec<JunctionEdge>,
}

impl FamilyGraph {
    pub fn n_nodes(&self) -> usize { self.nodes.len() }
    pub fn n_edges(&self) -> usize { self.edges.len() }
}

impl Index<NodeIdx> for FamilyGraph {
    type Output = ExonClass;
    fn index(&self, idx: NodeIdx) -> &ExonClass { &self.nodes[idx.0] }
}

use crate::types::Bundle;

/// Collect, dedup, and sort exonic intervals from a bundle's reads.
/// Returns half-open (start, end) pairs on the bundle's chromosome.
pub fn extract_copy_exons(bundle: &Bundle) -> Vec<(u64, u64)> {
    let mut all: Vec<(u64, u64)> = bundle.reads.iter()
        .flat_map(|r| r.exons.iter().copied())
        .collect();
    all.sort_unstable();
    all.dedup();
    all
}

/// One member of an exon cluster: (copy id, exon index within copy).
pub type ExonRef = (CopyId, usize);

/// Cluster exons across copies by reciprocal-overlap fraction on the same
/// (chrom, strand). `min_recip` is the minimum reciprocal-overlap (e.g. 0.30
/// = 30%) for two exons to join the same cluster. Single-copy exons emerge
/// as singleton clusters.
pub fn cluster_by_position(
    copies: &[(&str, char, Vec<(u64, u64)>)],
    min_recip: f64,
) -> Vec<Vec<ExonRef>> {
    // Flat list of (copy_id, exon_idx, chrom, strand, start, end).
    let mut all: Vec<(CopyId, usize, &str, char, u64, u64)> = Vec::new();
    for (cid, (chrom, strand, exons)) in copies.iter().enumerate() {
        for (ei, &(s, e)) in exons.iter().enumerate() {
            all.push((cid, ei, chrom, *strand, s, e));
        }
    }

    // Union-Find over the flat list.
    let n = all.len();
    let mut parent: Vec<usize> = (0..n).collect();
    fn find(p: &mut [usize], x: usize) -> usize {
        if p[x] == x { x } else { let r = find(p, p[x]); p[x] = r; r }
    }
    fn union(p: &mut [usize], a: usize, b: usize) {
        let ra = find(p, a); let rb = find(p, b);
        if ra != rb { p[ra] = rb; }
    }

    for i in 0..n {
        for j in (i + 1)..n {
            let (_, _, ci, si, ai, bi) = all[i];
            let (_, _, cj, sj, aj, bj) = all[j];
            if ci != cj || si != sj { continue; }
            // Reciprocal overlap fraction.
            let inter_s = ai.max(aj); let inter_e = bi.min(bj);
            if inter_e <= inter_s { continue; }
            let inter = (inter_e - inter_s) as f64;
            let li = (bi - ai) as f64; let lj = (bj - aj) as f64;
            let r = (inter / li).min(inter / lj);
            if r >= min_recip { union(&mut parent, i, j); }
        }
    }

    // Group flat indices by root.
    use std::collections::BTreeMap;
    let mut groups: BTreeMap<usize, Vec<ExonRef>> = BTreeMap::new();
    for i in 0..n {
        let r = find(&mut parent, i);
        groups.entry(r).or_default().push((all[i].0, all[i].1));
    }
    groups.into_values().collect()
}

use std::collections::HashSet;

fn fnv1a(bytes: &[u8]) -> u64 {
    let mut h: u64 = 0xcbf29ce484222325;
    for &b in bytes { h ^= b as u64; h = h.wrapping_mul(0x100000001b3); }
    h
}

fn minimizers(seq: &[u8], k: usize, w: usize) -> HashSet<u64> {
    let mut out = HashSet::new();
    if seq.len() < k { return out; }
    let n = seq.len() - k + 1;
    for win_start in 0..n.saturating_sub(w).max(1) {
        let win_end = (win_start + w).min(n);
        let mut best: Option<u64> = None;
        for i in win_start..win_end {
            let h = fnv1a(&seq[i..i + k]);
            best = Some(best.map_or(h, |b| b.min(h)));
        }
        if let Some(h) = best { out.insert(h); }
    }
    out
}

/// Within a position-overlap cluster, split by minimizer-Jaccard similarity.
pub fn refine_by_minimizer_jaccard(
    cluster: &[(CopyId, Vec<u8>)],
    min_jaccard: f64,
    k: usize,
    w: usize,
) -> Vec<Vec<CopyId>> {
    let n = cluster.len();
    let mins: Vec<HashSet<u64>> = cluster.iter().map(|(_, s)| minimizers(s, k, w)).collect();

    let mut parent: Vec<usize> = (0..n).collect();
    fn find(p: &mut [usize], x: usize) -> usize {
        if p[x] == x { x } else { let r = find(p, p[x]); p[x] = r; r }
    }
    fn union(p: &mut [usize], a: usize, b: usize) {
        let ra = find(p, a); let rb = find(p, b);
        if ra != rb { p[ra] = rb; }
    }

    for i in 0..n {
        for j in (i + 1)..n {
            let inter = mins[i].intersection(&mins[j]).count() as f64;
            let union_sz = mins[i].union(&mins[j]).count() as f64;
            if union_sz == 0.0 { continue; }
            if inter / union_sz >= min_jaccard { union(&mut parent, i, j); }
        }
    }
    use std::collections::BTreeMap;
    let mut g: BTreeMap<usize, Vec<CopyId>> = BTreeMap::new();
    for i in 0..n {
        let r = find(&mut parent, i);
        g.entry(r).or_default().push(cluster[i].0);
    }
    g.into_values().collect()
}
