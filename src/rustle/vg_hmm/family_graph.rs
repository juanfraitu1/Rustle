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
