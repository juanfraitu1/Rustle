//! Copy-graph objects (v1): every copy of a family is a tagged, corroborable PATH in one GFA 1.1
//! variation graph. A REFERENCE walk makes a reference-absent copy visibly an arm the reference does
//! not take. Pure builder — no I/O; the caller fills the parallel vectors and writes the strings.

use std::collections::BTreeSet;

/// Per-copy status across the (in-genome / annotated) axes and the absent subtypes.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum CopyStatus {
    Reference,
    InGenomeAnnotated,
    InGenomeUnannotated,
    AnnotationUnknown,
    AbsentCollapsed,
    AbsentDivergent,
}

impl CopyStatus {
    /// `ST:Z:` tag value.
    pub fn tag(&self) -> &'static str {
        match self {
            CopyStatus::Reference => "reference",
            CopyStatus::InGenomeAnnotated => "in-genome-annotated",
            CopyStatus::InGenomeUnannotated => "in-genome-unannotated",
            CopyStatus::AnnotationUnknown => "annotation-unknown",
            CopyStatus::AbsentCollapsed => "absent-collapsed",
            CopyStatus::AbsentDivergent => "absent-divergent",
        }
    }
    pub fn is_absent(&self) -> bool {
        matches!(self, CopyStatus::AbsentCollapsed | CopyStatus::AbsentDivergent)
    }
    /// Bandage colour for arms unique to this status.
    pub fn colour(&self) -> &'static str {
        match self {
            CopyStatus::Reference => "#9aa0a6",
            CopyStatus::AbsentCollapsed | CopyStatus::AbsentDivergent => "#d93025",
            CopyStatus::InGenomeUnannotated => "#188038",
            CopyStatus::InGenomeAnnotated => "#1a73e8",
            CopyStatus::AnnotationUnknown => "#a142f4",
        }
    }
}

/// Corroboration evidence carried as GFA tags. `None` => tag omitted (never faked).
#[derive(Clone, Debug, Default)]
pub struct Corrob {
    pub reads: Option<u32>,        // RC:i:
    pub suns: Option<u32>,         // SU:i: (filled by the builder if left None)
    pub map_identity: Option<f64>, // MI:f:
}

/// One PSV column, already known to be usable (genome_pos + ref_allele both Some).
#[derive(Clone, Debug)]
pub struct PsvColumn {
    pub col: usize,             // original column index (for provenance only)
    pub genome_pos: Option<u64>,
    pub ref_allele: Option<u8>,
}

/// One copy as a path: its allele per column (None = gap => routes through the reference allele node).
#[derive(Clone, Debug)]
pub struct CopyPath {
    pub id: String,
    pub alleles: Vec<Option<u8>>,
    pub status: CopyStatus,
    pub corrob: Corrob,
}

/// One read as a walk over the columns it observed (None = unobserved).
#[derive(Clone, Debug)]
pub struct ReadWalk {
    pub name: String,
    pub obs: Vec<Option<u8>>,
    pub assigned_copy: Option<usize>, // index into CopyGraph.copies; None = tied/K=0 (grey)
}

/// A whole family's variation graph. columns, every copy.alleles, every read.obs are length M and
/// share column order; backbone is length M+1.
#[derive(Clone, Debug)]
pub struct CopyGraph {
    pub family: String,
    pub columns: Vec<PsvColumn>,
    pub backbone: Vec<Vec<u8>>,
    pub copies: Vec<CopyPath>,
    pub reads: Vec<ReadWalk>,
}

/// Assembled GFA line groups (dedup + ordering handled by the caller/writer).
#[derive(Default, Debug)]
pub struct GfaLines {
    pub header: String,
    pub segs: Vec<String>,
    pub links: Vec<String>,
    pub paths: Vec<String>,
    pub walks: Vec<String>,
}

impl CopyGraph {
    fn m(&self) -> usize { self.columns.len() }
    fn bb(&self, i: usize) -> String { format!("{}_bb{}", self.family, i) }
    fn allele_node(&self, ci: usize, b: u8) -> String {
        format!("{}_c{}_{}", self.family, ci, b as char)
    }

    /// The set of alleles present at column `ci` (reference ∪ all copies ∪ all reads), sorted.
    fn alleles_at(&self, ci: usize) -> BTreeSet<u8> {
        let mut set = BTreeSet::new();
        if let Some(b) = self.columns[ci].ref_allele { set.insert(b); }
        for c in &self.copies {
            if let Some(Some(b)) = c.alleles.get(ci) { set.insert(*b); }
        }
        for r in &self.reads {
            if let Some(Some(b)) = r.obs.get(ci) { set.insert(*b); }
        }
        set
    }

    /// Walk string "bb0+,c0_x+,bb1+,...,bbM+" given the allele taken at each column (`taken[ci]`).
    /// A `None` in `taken` routes through the reference allele node at that column.
    fn walk_tokens(&self, taken: &[Option<u8>]) -> Vec<String> {
        let m = self.m();
        let mut toks = Vec::with_capacity(2 * m + 1);
        for ci in 0..m {
            toks.push(format!("{}+", self.bb(ci)));
            let b = taken.get(ci).and_then(|o| *o).or(self.columns[ci].ref_allele);
            if let Some(b) = b {
                toks.push(format!("{}+", self.allele_node(ci, b)));
            }
        }
        toks.push(format!("{}+", self.bb(m)));
        toks
    }

    pub fn gfa_lines(&self) -> GfaLines {
        let mut out = GfaLines { header: "H\tVN:Z:1.1".into(), ..Default::default() };
        let m = self.m();
        // backbone spacer S-nodes bb0..=bbM
        for i in 0..=m {
            let seq = String::from_utf8_lossy(&self.backbone[i]).to_string();
            out.segs.push(format!("S\t{}\t{}\tSN:Z:spacer", self.bb(i), seq));
        }
        // allele S-nodes + L-lines bb{ci} -> allele -> bb{ci+1}
        for ci in 0..m {
            let pos = self.columns[ci].genome_pos.unwrap_or(0);
            for b in self.alleles_at(ci) {
                let nid = self.allele_node(ci, b);
                out.segs.push(format!("S\t{}\t{}\tPO:i:{}", nid, b as char, pos));
                out.links.push(format!("L\t{}\t+\t{}\t+\t0M", self.bb(ci), nid));
                out.links.push(format!("L\t{}\t+\t{}\t+\t0M", nid, self.bb(ci + 1)));
            }
        }
        // REFERENCE walk: the genome's own allele at each column.
        let ref_taken: Vec<Option<u8>> = self.columns.iter().map(|c| c.ref_allele).collect();
        let ref_walk = self.walk_tokens(&ref_taken).join(",");
        out.paths.push(format!("P\t{}_REFERENCE\t{}\t*\tST:Z:reference", self.family, ref_walk));
        out
    }

    /// One self-contained GFA string (header + this family's lines). Convenience for tests / single-family use.
    pub fn to_gfa(&self) -> String {
        let g = self.gfa_lines();
        let mut s = String::new();
        s.push_str(&g.header); s.push('\n');
        for l in g.segs.iter().chain(g.links.iter()).chain(g.paths.iter()).chain(g.walks.iter()) {
            s.push_str(l); s.push('\n');
        }
        s
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn tiny_graph() -> CopyGraph {
        // 2 columns, backbone spacers of len 3, one reference + one copy
        CopyGraph {
            family: "FAM1".into(),
            columns: vec![
                PsvColumn { col: 0, genome_pos: Some(100), ref_allele: Some(b'A') },
                PsvColumn { col: 1, genome_pos: Some(200), ref_allele: Some(b'C') },
            ],
            backbone: vec![b"NNN".to_vec(), b"NNN".to_vec(), b"NNN".to_vec()],
            copies: vec![CopyPath {
                id: "FAM1_copy0".into(),
                alleles: vec![Some(b'A'), Some(b'G')],
                status: CopyStatus::InGenomeAnnotated,
                corrob: Corrob { reads: Some(5), suns: None, map_identity: Some(0.99) },
            }],
            reads: vec![],
        }
    }

    #[test]
    fn constructs_and_reports_shape() {
        let g = tiny_graph();
        assert_eq!(g.columns.len(), 2);
        assert_eq!(g.backbone.len(), 3);
        assert_eq!(g.copies[0].status.tag(), "in-genome-annotated");
        assert!(!g.copies[0].status.is_absent());
        assert!(CopyStatus::AbsentCollapsed.is_absent());
    }

    fn parse_line_prefixes(gfa: &str) -> (usize, usize, usize) {
        let (mut s, mut l, mut h) = (0, 0, 0);
        for line in gfa.lines() {
            match line.chars().next() {
                Some('S') => s += 1,
                Some('L') => l += 1,
                Some('H') => h += 1,
                _ => {}
            }
        }
        (h, s, l)
    }

    #[test]
    fn skeleton_has_backbone_alleles_and_links() {
        let g = tiny_graph(); // col0 alleles {A(ref), A(copy)} => {A}; col1 alleles {C(ref), G(copy)} => {C,G}
        let gfa = g.to_gfa();
        let (h, s, _l) = parse_line_prefixes(&gfa);
        assert_eq!(h, 1, "one header");
        // backbone: bb0,bb1,bb2 (3) + allele nodes: col0 {A}=1, col1 {C,G}=2 => 3 => total S = 6
        assert_eq!(s, 6, "3 backbone + 3 allele S-nodes");
        assert!(gfa.contains("S\tFAM1_bb0\tNNN"));
        assert!(gfa.contains("S\tFAM1_c0_A\tA\tPO:i:100"));
        assert!(gfa.contains("S\tFAM1_c1_G\tG\tPO:i:200"));
        // every allele node linked to its flanking backbone (no dangling by construction)
        assert!(gfa.contains("L\tFAM1_bb0\t+\tFAM1_c0_A\t+\t0M"));
        assert!(gfa.contains("L\tFAM1_c0_A\t+\tFAM1_bb1\t+\t0M"));
        assert!(gfa.contains("L\tFAM1_bb1\t+\tFAM1_c1_G\t+\t0M"));
    }

    #[test]
    fn reference_walk_threads_reference_alleles() {
        let g = tiny_graph();
        let gfa = g.to_gfa();
        // reference alleles are A (col0) and C (col1)
        assert!(gfa.contains(
            "P\tFAM1_REFERENCE\tFAM1_bb0+,FAM1_c0_A+,FAM1_bb1+,FAM1_c1_C+,FAM1_bb2+\t*\tST:Z:reference"
        ), "reference P-line missing or wrong:\n{}", gfa);
    }
}
