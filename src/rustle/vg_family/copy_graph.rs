//! Copy-graph objects (v1): every copy of a family is a tagged, corroborable PATH in one GFA 1.1
//! variation graph. A REFERENCE walk makes a reference-absent copy visibly an arm the reference does
//! not take. Pure builder — no I/O; the caller fills the parallel vectors and writes the strings.

use std::collections::BTreeSet;

/// Neutral faint colour for allele nodes observed ONLY in reads (carried by no CopyPath and unequal to
/// the reference allele). Distinct from the backbone light-grey so read-only arms stay legible in Bandage.
const READ_ONLY_COLOUR: &str = "#e8eaed";

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

/// A shared exon node in the exon presence/absence graph — one genomic exon interval.
#[derive(Clone, Debug)]
pub struct ExonClass { pub chrom: String, pub start: u64, pub end: u64 }

/// A copy as an ordered walk over exon-class indices.
#[derive(Clone, Debug)]
pub struct CopyExonPath { pub id: String, pub exon_nodes: Vec<usize>, pub status: CopyStatus, pub corrob: Corrob }

/// Family exon presence/absence graph. `nodes` sorted by genomic start; each copy walks a subset.
#[derive(Clone, Debug)]
pub struct ExonGraph { pub family: String, pub nodes: Vec<ExonClass>, pub copies: Vec<CopyExonPath> }

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

    /// Number of columns where copy `c`'s allele is BOTH unique among the family's copies AND differs
    /// from the reference allele — a private *divergent* marker (SUN). A copy identical to the
    /// reference scores 0. (This reference-exclusion filter deviates from the plan's illustrative
    /// snippet but matches the acceptance test's SUN semantics.)
    fn private_columns(&self, c: usize) -> u32 {
        let mut n = 0u32;
        for ci in 0..self.m() {
            let Some(Some(b)) = self.copies[c].alleles.get(ci) else { continue };
            // Only count columns where this allele differs from the reference
            if self.columns[ci].ref_allele == Some(*b) { continue; }
            // Check if no other copy has this same allele
            let unique = self.copies.iter().enumerate()
                .all(|(k, other)| k == c || other.alleles.get(ci).and_then(|o| *o) != Some(*b));
            if unique { n += 1; }
        }
        n
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

        // copy P-lines with corroboration tags
        for (copy_idx, cp) in self.copies.iter().enumerate() {
            let walk = self.walk_tokens(&cp.alleles).join(",");
            let name = if cp.status.is_absent() {
                format!("{}_copy{}_ABSENT", self.family, copy_idx)
            } else {
                format!("{}_copy{}", self.family, copy_idx)
            };
            let mut tags = String::new();
            if let Some(rc) = cp.corrob.reads { tags.push_str(&format!("\tRC:i:{}", rc)); }
            let su = cp.corrob.suns.unwrap_or_else(|| self.private_columns(copy_idx));
            tags.push_str(&format!("\tSU:i:{}", su));
            if let Some(mi) = cp.corrob.map_identity { tags.push_str(&format!("\tMI:f:{:.3}", mi)); }
            tags.push_str(&format!("\tST:Z:{}", cp.status.tag()));
            out.paths.push(format!("P\t{}\t{}\t*{}", name, walk, tags));
        }

        // read W-lines over each read's observed span; gaps within span route through the reference node.
        for r in &self.reads {
            let first = r.obs.iter().position(|o| o.is_some());
            let last = r.obs.iter().rposition(|o| o.is_some());
            let (Some(first), Some(last)) = (first, last) else { continue };
            let mut toks: Vec<String> = Vec::new();
            for ci in first..=last {
                toks.push(format!(">{}", self.bb(ci)));
                let b = r.obs[ci].or(self.columns[ci].ref_allele);
                if let Some(b) = b {
                    toks.push(format!(">{}", self.allele_node(ci, b)));
                }
            }
            toks.push(format!(">{}", self.bb(last + 1)));
            let w = toks.join("");
            let hap = r.assigned_copy.map(|c| c as i64).unwrap_or(-1).max(0);
            out.walks.push(format!("W\t{}\t{}\t{}\t0\t{}\t{}", r.name, hap, self.family, toks.len(), w));
        }

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

    /// Bandage node colours (keyed on SEGMENT names): reference-walk nodes grey, absent-only divergent
    /// nodes red, other copy-divergent nodes their copy's status colour, backbone light grey.
    pub fn colours_csv(&self) -> String {
        use std::collections::BTreeMap;
        let mut colour: BTreeMap<String, &'static str> = BTreeMap::new();
        // backbone
        for i in 0..=self.m() {
            colour.insert(self.bb(i), "#dadce0");
        }
        for ci in 0..self.m() {
            let refb = self.columns[ci].ref_allele;
            for b in self.alleles_at(ci) {
                let nid = self.allele_node(ci, b);
                if Some(b) == refb {
                    colour.insert(nid, CopyStatus::Reference.colour());
                    continue;
                }
                // walked by any absent copy? (and not the reference allele) — absent wins over non-absent.
                if let Some(c) = self.copies.iter()
                    .find(|c| c.status.is_absent() && c.alleles.get(ci).and_then(|o| *o) == Some(b)) {
                    colour.insert(nid, c.status.colour());
                } else if let Some(c) = self.copies.iter()
                    .find(|c| c.alleles.get(ci).and_then(|o| *o) == Some(b)) {
                    colour.insert(nid, c.status.colour());
                } else {
                    // observed only in reads (no copy carries it, not the reference) — neutral read-only.
                    colour.insert(nid, READ_ONLY_COLOUR);
                }
            }
        }
        let mut s = String::new();
        for (k, v) in colour { s.push_str(&format!("{},{}\n", k, v)); }
        s
    }

    /// Legend: each status actually present (plus reference) → its colour.
    pub fn legend_tsv(&self) -> String {
        use std::collections::BTreeSet;
        let mut statuses: BTreeSet<&'static str> = BTreeSet::new();
        statuses.insert("reference");
        let mut rows: Vec<(&'static str, &'static str)> = vec![("reference", CopyStatus::Reference.colour())];
        for c in &self.copies {
            if statuses.insert(c.status.tag()) {
                rows.push((c.status.tag(), c.status.colour()));
            }
        }
        let mut s = String::new();
        for (st, col) in rows { s.push_str(&format!("{}\t{}\n", st, col)); }
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

    #[test]
    fn copy_paths_carry_tags_and_absent_diverges() {
        // 3 columns; reference = A,A,A. copy0 in-genome matches ref. copy1 ABSENT diverges at col1 & col2.
        let g = CopyGraph {
            family: "FAM2".into(),
            columns: (0..3).map(|i| PsvColumn {
                col: i, genome_pos: Some(100 + i as u64), ref_allele: Some(b'A'),
            }).collect(),
            backbone: vec![b"NN".to_vec(); 4],
            copies: vec![
                CopyPath { id: "FAM2_copy0".into(), alleles: vec![Some(b'A'), Some(b'A'), Some(b'A')],
                    status: CopyStatus::InGenomeAnnotated,
                    corrob: Corrob { reads: Some(8), suns: None, map_identity: Some(0.998) } },
                CopyPath { id: "FAM2_copy1".into(), alleles: vec![Some(b'A'), Some(b'G'), Some(b'T')],
                    status: CopyStatus::AbsentDivergent,
                    corrob: Corrob { reads: Some(12), suns: None, map_identity: Some(0.952) } },
            ],
            reads: vec![],
        };
        let gfa = g.to_gfa();
        // absent copy P-line named with _ABSENT and tagged
        let absent = gfa.lines().find(|l| l.starts_with("P\tFAM2_copy1_ABSENT")).expect("absent P-line");
        assert!(absent.contains("RC:i:12"));
        assert!(absent.contains("MI:f:0.952"));
        assert!(absent.contains("ST:Z:absent-divergent"));
        // SU (private columns): copy1's allele is unique (vs copy0) at col1(G) and col2(T) => SU:i:2
        assert!(absent.contains("SU:i:2"), "expected SU:i:2 in: {}", absent);
        // it walks the divergent nodes the reference walk does NOT (c1_G, c2_T)
        assert!(absent.contains("FAM2_c1_G+"));
        assert!(absent.contains("FAM2_c2_T+"));
        // in-genome copy0 has SU:i:0 (never unique) and is not _ABSENT
        let c0 = gfa.lines().find(|l| l.starts_with("P\tFAM2_copy0\t")).expect("copy0 P-line");
        assert!(c0.contains("SU:i:0"));
        assert!(c0.contains("ST:Z:in-genome-annotated"));
    }

    #[test]
    fn omits_unknown_corrob_tags() {
        // Honesty rule NEGATIVE path: when reads/map_identity are None, RC:i: and MI:f: are OMITTED,
        // while SU:i: (always computed) and ST:Z: (always emitted) remain.
        let g = CopyGraph {
            family: "FAM3".into(),
            columns: vec![
                PsvColumn { col: 0, genome_pos: Some(100), ref_allele: Some(b'A') },
                PsvColumn { col: 1, genome_pos: Some(200), ref_allele: Some(b'C') },
            ],
            backbone: vec![b"NN".to_vec(); 3],
            copies: vec![CopyPath {
                id: "FAM3_copy0".into(),
                alleles: vec![Some(b'A'), Some(b'G')],
                status: CopyStatus::AnnotationUnknown,
                corrob: Corrob { reads: None, suns: None, map_identity: None },
            }],
            reads: vec![],
        };
        let gfa = g.to_gfa();
        let cp = gfa.lines().find(|l| l.starts_with("P\tFAM3_copy0\t")).expect("copy0 P-line");
        // unknown values => tags omitted (never faked)
        assert!(!cp.contains("RC:i:"), "RC:i: must be omitted when reads is None: {}", cp);
        assert!(!cp.contains("MI:f:"), "MI:f: must be omitted when map_identity is None: {}", cp);
        // always-present tags
        assert!(cp.contains("SU:i:"), "SU:i: must always be present: {}", cp);
        assert!(cp.contains("ST:Z:annotation-unknown"), "ST:Z: must always be present: {}", cp);
    }

    // Assert every P-line and W-line step is backed by an L-line (parses walks, checks adjacency set).
    fn assert_no_dangling(gfa: &str) {
        use std::collections::HashSet;
        let mut links: HashSet<(String, String)> = HashSet::new();
        for l in gfa.lines().filter(|l| l.starts_with("L\t")) {
            let f: Vec<&str> = l.split('\t').collect(); // L from + to + 0M
            links.insert((f[1].to_string(), f[3].to_string()));
        }
        let node = |tok: &str| tok.trim_start_matches(['>', '<']).trim_end_matches(['>', '<', '+', '-']).to_string();
        for l in gfa.lines() {
            let seq: Vec<String> = if l.starts_with("P\t") {
                l.split('\t').nth(2).unwrap().split(',').map(node).collect()
            } else if l.starts_with("W\t") {
                let w = l.split('\t').nth(6).unwrap();
                w.split_inclusive(['>', '<']).filter(|s| s.len() > 1).map(node).collect()
            } else { continue };
            for pair in seq.windows(2) {
                assert!(links.contains(&(pair[0].clone(), pair[1].clone())),
                    "dangling walk edge {}->{} in line: {}", pair[0], pair[1], l);
            }
        }
    }

    #[test]
    fn reads_walk_with_backing_links() {
        let mut g = tiny_graph(); // 2 cols, ref A,C
        g.reads = vec![
            ReadWalk { name: "readX".into(), obs: vec![Some(b'A'), Some(b'C')], assigned_copy: Some(0) },
            ReadWalk { name: "readY".into(), obs: vec![None, Some(b'C')], assigned_copy: None },
        ];
        let gfa = g.to_gfa();
        assert!(gfa.lines().any(|l| l.starts_with("W\treadX")), "readX walk missing");
        assert!(gfa.lines().any(|l| l.starts_with("W\treadY")), "readY walk missing");
        assert_no_dangling(&gfa);
    }

    #[test]
    fn read_walk_internal_gap_routes_through_reference() {
        // 3 columns, reference A,A,A. A read observes col0 (A) and col2 (T) but NOT col1 (internal gap):
        // the gap at col1 must route through the reference allele node c1_A. The walk spans
        // bb(first_obs=0) .. bb(last_obs+1=3).
        let g = CopyGraph {
            family: "FAM4".into(),
            columns: (0..3).map(|i| PsvColumn {
                col: i, genome_pos: Some(100 + i as u64), ref_allele: Some(b'A'),
            }).collect(),
            backbone: vec![b"NN".to_vec(); 4],
            copies: vec![],
            reads: vec![ReadWalk {
                name: "readG".into(),
                obs: vec![Some(b'A'), None, Some(b'T')],
                assigned_copy: None,
            }],
        };
        let gfa = g.to_gfa();
        // exact W-line: gap col1 routes through reference node c1_A; starts bb0, ends bb3; 7 tokens.
        assert!(gfa.lines().any(|l| l ==
            "W\treadG\t0\tFAM4\t0\t7\t>FAM4_bb0>FAM4_c0_A>FAM4_bb1>FAM4_c1_A>FAM4_bb2>FAM4_c2_T>FAM4_bb3"
        ), "internal-gap W-line missing or wrong:\n{}", gfa);
        assert_no_dangling(&gfa);
    }

    #[test]
    fn read_observing_zero_columns_emits_no_walk() {
        // A read with no observations (all None) must emit NO W-line (early continue).
        let mut g = tiny_graph(); // 2 cols
        g.reads = vec![ReadWalk {
            name: "readEmpty".into(),
            obs: vec![None, None],
            assigned_copy: None,
        }];
        let gfa = g.to_gfa();
        assert!(!gfa.lines().any(|l| l.starts_with("W\t")), "no W-line expected for a read with zero observations:\n{}", gfa);
    }

    #[test]
    fn colours_mark_absent_red_reference_grey() {
        // reuse the 3-column absent-copy graph
        let g = CopyGraph {
            family: "FAM3".into(),
            columns: (0..3).map(|i| PsvColumn { col: i, genome_pos: Some(10 + i as u64), ref_allele: Some(b'A') }).collect(),
            backbone: vec![b"NN".to_vec(); 4],
            copies: vec![
                CopyPath { id: "FAM3_copy0".into(), alleles: vec![Some(b'A'), Some(b'A'), Some(b'A')],
                    status: CopyStatus::InGenomeAnnotated, corrob: Corrob::default() },
                CopyPath { id: "FAM3_copy1".into(), alleles: vec![Some(b'A'), Some(b'G'), Some(b'T')],
                    status: CopyStatus::AbsentDivergent, corrob: Corrob::default() },
            ],
            reads: vec![],
        };
        let csv = g.colours_csv();
        // reference allele node grey
        assert!(csv.lines().any(|l| l == "FAM3_c0_A,#9aa0a6"), "ref node not grey:\n{}", csv);
        // absent-only divergent nodes red
        assert!(csv.lines().any(|l| l == "FAM3_c1_G,#d93025"), "absent node not red:\n{}", csv);
        assert!(csv.lines().any(|l| l == "FAM3_c2_T,#d93025"));
        // legend lists the two statuses in use
        let legend = g.legend_tsv();
        assert!(legend.contains("reference\t#9aa0a6"));
        assert!(legend.contains("absent-divergent\t#d93025"));
    }

    #[test]
    fn colours_read_only_allele_gets_neutral() {
        // A base observed ONLY in a read (differs from ref_allele AND carried by no CopyPath) still
        // gets a GFA allele segment via alleles_at — it must receive the neutral read-only colour, not
        // be silently dropped from colours.csv (which would render uncoloured in Bandage).
        let g = CopyGraph {
            family: "FAM5".into(),
            columns: vec![
                PsvColumn { col: 0, genome_pos: Some(100), ref_allele: Some(b'A') },
            ],
            backbone: vec![b"NN".to_vec(); 2],
            copies: vec![CopyPath {
                id: "FAM5_copy0".into(), alleles: vec![Some(b'A')],
                status: CopyStatus::InGenomeAnnotated, corrob: Corrob::default(),
            }],
            // read observes 'G' at col0 — neither the reference (A) nor any copy (A) carries it.
            reads: vec![ReadWalk {
                name: "readR".into(), obs: vec![Some(b'G')], assigned_copy: None,
            }],
        };
        let csv = g.colours_csv();
        // the read-only node exists in the GFA (alleles_at folds it in)…
        assert!(g.to_gfa().contains("S\tFAM5_c0_G\tG"), "read-only allele node missing from GFA:\n{}", g.to_gfa());
        // …and it must be deliberately coloured neutral, distinct from the backbone light-grey.
        assert!(csv.lines().any(|l| l == "FAM5_c0_G,#e8eaed"), "read-only node not neutral:\n{}", csv);
        // reference allele still grey
        assert!(csv.lines().any(|l| l == "FAM5_c0_A,#9aa0a6"), "ref node not grey:\n{}", csv);
    }

    #[test]
    fn colours_absent_wins_over_non_absent_at_shared_node() {
        // A single non-reference allele node walked by BOTH an absent copy and a non-absent copy at the
        // same column/base must come out RED (absent precedence), not the non-absent status colour.
        let g = CopyGraph {
            family: "FAM6".into(),
            columns: vec![
                PsvColumn { col: 0, genome_pos: Some(100), ref_allele: Some(b'A') },
            ],
            backbone: vec![b"NN".to_vec(); 2],
            copies: vec![
                // non-absent copy walks G at col0…
                CopyPath { id: "FAM6_copy0".into(), alleles: vec![Some(b'G')],
                    status: CopyStatus::InGenomeAnnotated, corrob: Corrob::default() },
                // …and an absent copy walks the SAME G at col0.
                CopyPath { id: "FAM6_copy1".into(), alleles: vec![Some(b'G')],
                    status: CopyStatus::AbsentDivergent, corrob: Corrob::default() },
            ],
            reads: vec![],
        };
        let csv = g.colours_csv();
        // absent wins: the shared node is red, NOT the in-genome blue (#1a73e8).
        assert!(csv.lines().any(|l| l == "FAM6_c0_G,#d93025"),
            "shared absent/non-absent node must be red (absent wins):\n{}", csv);
        assert!(!csv.lines().any(|l| l == "FAM6_c0_G,#1a73e8"),
            "shared node must NOT take the non-absent colour:\n{}", csv);
    }

    #[test]
    fn exon_graph_constructs() {
        let g = ExonGraph {
            family: "F".into(),
            nodes: vec![ExonClass { chrom: "c".into(), start: 0, end: 100 }],
            copies: vec![CopyExonPath { id: "F_copy0".into(), exon_nodes: vec![0],
                status: CopyStatus::InGenomeAnnotated, corrob: Corrob::default() }],
        };
        assert_eq!(g.nodes.len(), 1);
        assert_eq!(g.copies[0].exon_nodes, vec![0]);
    }
}
