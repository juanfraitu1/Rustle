//! Copy-graph objects (v1): every copy of a family is a tagged, corroborable PATH in one GFA 1.1
//! variation graph. A REFERENCE walk makes a reference-absent copy visibly an arm the reference does
//! not take. Pure builder — no I/O; the caller fills the parallel vectors and writes the strings.
//!
//! **STATUS:** OPT-IN — --phase (src/bin/copy_assign.rs:235-236, default_value_t = false)  (docs/MODULE_STATUS.md; assigned by reachability, not by this header)

use std::collections::{BTreeMap, BTreeSet};

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

/// Per-read significance certificate carried onto the audit W-line.
#[derive(Clone, Debug)]
pub struct ReadCert {
    pub p_value: f64,
    pub min_p_value: f64,
    pub status: super::copy_assign::AssignStatus,
}

/// One read as a walk over the columns it observed (None = unobserved).
#[derive(Clone, Debug)]
pub struct ReadWalk {
    pub name: String,
    pub obs: Vec<Option<u8>>,
    pub assigned_copy: Option<usize>, // index into CopyGraph.copies; None = tied/K=0 (grey)
    /// Significance certificate from copy_assign (p_value/min_p_value/status). `None` => the old,
    /// untagged W-line (backward-compatible / opt-in); `Some` appends CP/PV/MP/ST tags.
    pub cert: Option<ReadCert>,
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
            let mut line = format!("W\t{}\t{}\t{}\t0\t{}\t{}", r.name, hap, self.family, toks.len(), w);
            if let Some(c) = &r.cert {
                use super::copy_assign::AssignStatus::*;
                let st = match c.status { Assigned => "Assigned", Ambiguous => "Ambiguous", Tied => "Tied" };
                let cp = r.assigned_copy.map(|c| format!("copy{c}")).unwrap_or_else(|| "none".into());
                line.push_str(&format!("\tCP:Z:{cp}\tPV:f:{}\tMP:f:{}\tST:Z:{st}", c.p_value, c.min_p_value));
            }
            out.walks.push(line);
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

impl ExonGraph {
    fn node(&self, k: usize) -> String {
        format!("{}_E{}", self.family, k)
    }

    /// Bandage node colours (keyed on exon node names): a class walked by ≥1 non-absent copy → grey;
    /// else if only absent copies walk it → the owner's colour (red). Skip classes with no walkers.
    pub fn colours_csv(&self) -> String {
        let n = self.nodes.len();
        let mut colour: BTreeMap<String, &'static str> = BTreeMap::new();
        for k in 0..n {
            let walkers: Vec<&CopyExonPath> = self.copies.iter().filter(|c| c.exon_nodes.contains(&k)).collect();
            let on_ref = walkers.iter().any(|c| !c.status.is_absent());
            let col = if on_ref {
                CopyStatus::Reference.colour()               // grey shared/reference exon
            } else if let Some(c) = walkers.first() {
                c.status.colour()                            // copy-specific arm -> owner colour (red for absent)
            } else { continue };
            colour.insert(self.node(k), col);
        }
        let mut s = String::new();
        for (kk, v) in colour { s.push_str(&format!("{},{}\n", kk, v)); }
        s
    }

    /// Legend: each status actually present (plus reference) → its colour.
    pub fn legend_tsv(&self) -> String {
        let mut seen: BTreeSet<&'static str> = BTreeSet::new();
        let mut s = format!("reference\t{}\n", CopyStatus::Reference.colour());
        seen.insert("reference");
        for c in &self.copies {
            if seen.insert(c.status.tag()) {
                s.push_str(&format!("{}\t{}\n", c.status.tag(), c.status.colour()));
            }
        }
        s
    }

    /// Reciprocal overlap = min(inter/len_a, inter/len_b); 0 if disjoint or different chrom.
    fn recip_overlap(a: (&str, u64, u64), b: (&str, u64, u64)) -> f64 {
        if a.0 != b.0 { return 0.0; }
        let lo = a.1.max(b.1); let hi = a.2.min(b.2);
        if hi <= lo { return 0.0; }
        let inter = (hi - lo) as f64;
        (inter / (a.2 - a.1) as f64).min(inter / (b.2 - b.1) as f64)
    }

    pub fn to_gfa(&self, exon_seq: impl Fn(&ExonClass) -> Vec<u8>) -> String {
        // reference = classes present in >=1 non-absent copy; fallback = present in all copies
        let n = self.nodes.len();
        let non_absent: Vec<&CopyExonPath> = self.copies.iter().filter(|c| !c.status.is_absent()).collect();
        let mut ref_nodes: Vec<usize> = (0..n).filter(|&k|
            if non_absent.is_empty() { self.copies.iter().all(|c| c.exon_nodes.contains(&k)) }
            else { non_absent.iter().any(|c| c.exon_nodes.contains(&k)) }
        ).collect();
        ref_nodes.sort();

        // per-class RC = sum of reads over copies walking it
        let rc = |k: usize| -> u32 {
            self.copies.iter().filter(|c| c.exon_nodes.contains(&k)).filter_map(|c| c.corrob.reads).sum()
        };

        let mut s = String::from("H\tVN:Z:1.1\n");
        for (k, ec) in self.nodes.iter().enumerate() {
            let seq = String::from_utf8_lossy(&exon_seq(ec)).to_string();
            s.push_str(&format!("S\t{}\t{}\tPO:i:{}\tRC:i:{}\n", self.node(k), seq, ec.start, rc(k)));
        }
        // L-lines = union of consecutive adjacencies across reference + all copy walks
        use std::collections::BTreeSet;
        let mut links: BTreeSet<(usize, usize)> = BTreeSet::new();
        let mut add_walk = |walk: &[usize], set: &mut BTreeSet<(usize,usize)>| {
            for w in walk.windows(2) { set.insert((w[0], w[1])); }
        };
        add_walk(&ref_nodes, &mut links);
        for c in &self.copies { add_walk(&c.exon_nodes, &mut links); }
        for (a, b) in &links { s.push_str(&format!("L\t{}\t+\t{}\t+\t0M\n", self.node(*a), self.node(*b))); }
        // REFERENCE P-line
        let rwalk: String = ref_nodes.iter().map(|k| format!("{}+", self.node(*k))).collect::<Vec<_>>().join(",");
        s.push_str(&format!("P\t{}_REFERENCE\t{}\t*\tST:Z:reference\n", self.family, rwalk));
        // copy P-lines
        for (ci, c) in self.copies.iter().enumerate() {
            let name = if c.status.is_absent() { format!("{}_copy{}_ABSENT", self.family, ci) } else { format!("{}_copy{}", self.family, ci) };
            let walk: String = c.exon_nodes.iter().map(|k| format!("{}+", self.node(*k))).collect::<Vec<_>>().join(",");
            let mut tags = String::new();
            if let Some(r) = c.corrob.reads { tags.push_str(&format!("\tRC:i:{}", r)); }
            if let Some(mi) = c.corrob.map_identity { tags.push_str(&format!("\tMI:f:{:.3}", mi)); }
            tags.push_str(&format!("\tST:Z:{}", c.status.tag()));
            s.push_str(&format!("P\t{}\t{}\t*{}\n", name, walk, tags));
        }
        s
    }

    pub fn from_copies(family: &str, copies: &[(String, CopyStatus, Corrob, String, Vec<(u64,u64)>)]) -> ExonGraph {
        // flatten every exon as (copy_idx, chrom, start, end); skip malformed/zero-length (end <= start)
        // so the `(a.2 - a.1)` length divisor never underflows.
        let mut flat: Vec<(usize, String, u64, u64)> = Vec::new();
        for (ci, (_, _, _, chrom, exons)) in copies.iter().enumerate() {
            for &(s, e) in exons { if e > s { flat.push((ci, chrom.clone(), s, e)); } }
        }

        // union-find over flat exon items
        let n = flat.len();
        let mut parent: Vec<usize> = (0..n).collect();

        fn find(p: &mut Vec<usize>, x: usize) -> usize {
            if p[x] != x { let r = find(p, p[x]); p[x] = r; }
            p[x]
        }

        for i in 0..n {
            for j in (i+1)..n {
                if Self::recip_overlap((&flat[i].1, flat[i].2, flat[i].3), (&flat[j].1, flat[j].2, flat[j].3)) >= 0.30 {
                    let (a, b) = (find(&mut parent, i), find(&mut parent, j));
                    if a != b { parent[a] = b; }
                }
            }
        }

        // group items by root, build one ExonClass per group (min start, max end)
        let mut groups: BTreeMap<usize, Vec<usize>> = BTreeMap::new();
        for i in 0..n {
            let r = find(&mut parent, i);
            groups.entry(r).or_default().push(i);
        }

        let mut classes: Vec<(String, u64, u64, Vec<usize>)> = groups.values().map(|idxs| {
            let start = idxs.iter().map(|&i| flat[i].2).min().unwrap();
            let end = idxs.iter().map(|&i| flat[i].3).max().unwrap();
            let members: Vec<usize> = idxs.iter().map(|&i| flat[i].0).collect();
            (flat[idxs[0]].1.clone(), start, end, members)
        }).collect();

        classes.sort_by_key(|c| (c.1, c.2));   // by genomic (start, end) -> canonical E0..En order

        let nodes: Vec<ExonClass> = classes.iter().map(|(chrom, s, e, _)|
            ExonClass { chrom: chrom.clone(), start: *s, end: *e }
        ).collect();

        let copy_paths: Vec<CopyExonPath> = copies.iter().enumerate().map(|(ci, (id, status, corrob, _chrom, _exons))| {
            let mut exon_nodes: Vec<usize> = (0..classes.len()).filter(|&k| classes[k].3.contains(&ci)).collect();
            exon_nodes.sort();
            CopyExonPath { id: id.clone(), exon_nodes, status: *status, corrob: corrob.clone() }
        }).collect();

        ExonGraph { family: family.to_string(), nodes, copies: copy_paths }
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
            ReadWalk { name: "readX".into(), obs: vec![Some(b'A'), Some(b'C')], assigned_copy: Some(0), cert: None },
            ReadWalk { name: "readY".into(), obs: vec![None, Some(b'C')], assigned_copy: None, cert: None },
        ];
        let gfa = g.to_gfa();
        assert!(gfa.lines().any(|l| l.starts_with("W\treadX")), "readX walk missing");
        assert!(gfa.lines().any(|l| l.starts_with("W\treadY")), "readY walk missing");
        assert_no_dangling(&gfa);
    }

    #[test]
    fn read_walk_cert_tags_emitted_when_present() {
        use super::super::copy_assign::AssignStatus;
        let mut g = tiny_graph(); // 2 cols, ref A,C; copy0 alleles A,G
        g.reads = vec![
            ReadWalk { name: "r1".into(), obs: vec![Some(b'A'), Some(b'C')], assigned_copy: Some(0),
                       cert: Some(ReadCert { p_value: 0.001, min_p_value: 0.0005, status: AssignStatus::Assigned }) },
            ReadWalk { name: "r2".into(), obs: vec![Some(b'A'), Some(b'C')], assigned_copy: None, cert: None },
        ];
        let gfa = g.to_gfa();
        let r1 = gfa.lines().find(|l| l.starts_with("W\tr1")).expect("r1 walk");
        assert!(r1.contains("CP:Z:copy0") && r1.contains("PV:f:0.001") && r1.contains("MP:f:0.0005") && r1.contains("ST:Z:Assigned"));
        let r2 = gfa.lines().find(|l| l.starts_with("W\tr2")).expect("r2 walk");
        assert!(!r2.contains("CP:Z") && !r2.contains("PV:f"), "no cert -> no tags (backward-compatible)");
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
                cert: None,
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
            cert: None,
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
                name: "readR".into(), obs: vec![Some(b'G')], assigned_copy: None, cert: None,
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

    #[test]
    fn from_copies_clusters_and_flags_copy_specific_exon() {
        // copy0 exons E1,E3 ; copy1 exons E1,E2(extra),E3 — E2 is copy1-specific.
        let copies = vec![
            ("F_copy0".to_string(), CopyStatus::InGenomeAnnotated, Corrob::default(), "chr1".to_string(),
                vec![(0u64,100u64), (300,400)]),
            ("F_copy1".to_string(), CopyStatus::AbsentDivergent, Corrob::default(), "chr1".to_string(),
                vec![(0,100), (150,250), (300,400)]),
        ];
        let g = ExonGraph::from_copies("F", &copies);
        assert_eq!(g.nodes.len(), 3, "E1,E2,E3");
        // find the class only copy1 walks (the extra exon ~150-250)
        let owners: Vec<Vec<usize>> = (0..g.nodes.len())
            .map(|k| g.copies.iter().enumerate().filter(|(_, c)| c.exon_nodes.contains(&k)).map(|(i,_)| i).collect())
            .collect();
        let copy_specific: Vec<usize> = (0..g.nodes.len()).filter(|&k| owners[k] == vec![1]).collect();
        assert_eq!(copy_specific.len(), 1, "exactly one copy1-specific exon");
        // copy0 walks 2 classes, copy1 walks 3
        assert_eq!(g.copies[0].exon_nodes.len(), 2);
        assert_eq!(g.copies[1].exon_nodes.len(), 3);
    }

    #[test]
    fn from_copies_respects_overlap_threshold() {
        // A=(0,100), B=(70,170): inter=30 over len 100 => recip=0.30 => AT threshold => MERGE (1 class).
        let merge = vec![
            ("F_copy0".to_string(), CopyStatus::InGenomeAnnotated, Corrob::default(), "chr1".to_string(),
                vec![(0u64,100u64)]),
            ("F_copy1".to_string(), CopyStatus::InGenomeAnnotated, Corrob::default(), "chr1".to_string(),
                vec![(70u64,170u64)]),
        ];
        let g = ExonGraph::from_copies("F", &merge);
        assert_eq!(g.nodes.len(), 1, "recip overlap exactly 0.30 merges into one class");
        assert_eq!(g.copies[0].exon_nodes, vec![0]);
        assert_eq!(g.copies[1].exon_nodes, vec![0]);

        // A=(0,100), B'=(71,171): inter=29 over len 100 => recip=0.29 => below threshold => SEPARATE (2 classes).
        let split = vec![
            ("F_copy0".to_string(), CopyStatus::InGenomeAnnotated, Corrob::default(), "chr1".to_string(),
                vec![(0u64,100u64)]),
            ("F_copy1".to_string(), CopyStatus::InGenomeAnnotated, Corrob::default(), "chr1".to_string(),
                vec![(71u64,171u64)]),
        ];
        let g = ExonGraph::from_copies("F", &split);
        assert_eq!(g.nodes.len(), 2, "recip overlap 0.29 (< 0.30) stays two separate classes");
        // classes sorted by start: E0=(0,100) is copy0's, E1=(71,171) is copy1's.
        assert_eq!(g.copies[0].exon_nodes, vec![0]);
        assert_eq!(g.copies[1].exon_nodes, vec![1]);
    }

    #[test]
    fn exon_gfa_has_reference_skip_and_arm_no_dangling() {
        let copies = vec![
            ("F_copy0".to_string(), CopyStatus::InGenomeAnnotated, Corrob { reads: Some(10), suns: None, map_identity: None },
                "chr1".to_string(), vec![(0u64,100u64),(300,400)]),
            ("F_copy1".to_string(), CopyStatus::AbsentDivergent, Corrob { reads: Some(5), suns: None, map_identity: Some(0.95) },
                "chr1".to_string(), vec![(0,100),(150,250),(300,400)]),
        ];
        let g = ExonGraph::from_copies("F", &copies);
        let gfa = g.to_gfa(|ec| vec![b'A'; (ec.end - ec.start) as usize]);
        // reference exists and is the shared backbone (2 classes), copy1 absent walks 3
        assert!(gfa.contains("P\tF_REFERENCE"));
        let c1 = gfa.lines().find(|l| l.starts_with("P\tF_copy1_ABSENT")).unwrap();
        assert!(c1.contains("MI:f:0.950"));
        assert!(c1.contains("ST:Z:absent-divergent"));
        // the copy1-specific exon node exists with RC:i:5 (only copy1, 5 reads)
        let arm = g.copies[1].exon_nodes.iter().find(|&&k| !g.copies[0].exon_nodes.contains(&k)).copied().unwrap();
        assert!(gfa.contains(&format!("RC:i:5")));
        assert!(gfa.lines().any(|l| l.starts_with(&format!("S\tF_E{}", arm))));
        // no dangling: every P-line step is backed by an L-line
        assert_no_dangling(&gfa);
    }

    #[test]
    fn exon_colours_arm_red_shared_grey() {
        let copies = vec![
            ("F_copy0".to_string(), CopyStatus::InGenomeAnnotated, Corrob::default(), "chr1".to_string(), vec![(0u64,100u64),(300,400)]),
            ("F_copy1".to_string(), CopyStatus::AbsentDivergent, Corrob::default(), "chr1".to_string(), vec![(0,100),(150,250),(300,400)]),
        ];
        let g = ExonGraph::from_copies("F", &copies);
        let csv = g.colours_csv();
        let arm = g.copies[1].exon_nodes.iter().find(|&&k| !g.copies[0].exon_nodes.contains(&k)).copied().unwrap();
        assert!(csv.lines().any(|l| l == format!("F_E{},#d93025", arm)), "arm not red:\n{}", csv);
        // a shared class (walked by the in-genome copy0) is grey
        let shared = g.copies[0].exon_nodes[0];
        assert!(csv.lines().any(|l| l == format!("F_E{},#9aa0a6", shared)));
        assert!(g.legend_tsv().contains("absent-divergent\t#d93025"));
    }
}

// ---- merged 2026-10-05: was `vg_family/copy_discovery.rs`, now the inline module below (one component) ----
#[allow(clippy::all)]
pub mod copy_discovery {
//! Discovery of candidate gene-family copies from read alignment ties.
//!
//! **STATUS:** OPT-IN
//!
//! Reached from `copy_assign` behind `--discover-copies` (`src/bin/copy_assign.rs`, `default_value_t =
//! false`): `discover_copies_for_family` calls [`cluster_tie_partners`] once per family inside that
//! flag's own `if args.discover_copies {` block. REPORT ONLY -- the result is written to
//! `<out>.discovered_copies.tsv` and never mutates the input catalog or this run's own assignments.
//!
//! ⚠ A [`DiscoveredCopy`] row is a candidate ALIGNED-BLOCK CLUSTER, not necessarily a whole candidate
//! copy: clustering is per-block (see [`aligned_blocks`]), so a genuine multi-exon copy supported by only
//! a few reads can fragment into several rows, one per exon, each independently only needing its own
//! `TIE_PARTNER_MIN_SUPPORT` reads -- more chances for a coincidental cluster than a per-placement scheme
//! would give. Not detected or merged here (docs/o1_ledger.md §6l6 addendum, "Named limitation"); a
//! follow-up would regroup blocks sharing a supporting placement into one candidate with its own
//! `exon_blocks` column, mirroring `catalog_input::exon_blocks_str`.

use std::collections::{HashMap, HashSet};
use crate::vg_family::copy_split::AlignedRead;
use crate::vg_family::denovo_assemble::BamRead;

pub const TIE_PARTNER_MERGE_DISTANCE_BP: u64 = 500;
pub const TIE_PARTNER_MIN_SUPPORT: usize = 2;

#[derive(Clone, Debug, PartialEq)]
pub struct DiscoveredCopy {
    pub family_id: String,
    pub chrom: String,
    pub start: u64,
    pub end: u64,
    /// MAJORITY strand across the placements supporting this cluster (`false` -> `+`, `true` -> `-` from
    /// each record's own SAM FLAG 0x10), an exact tie defaulting to `+` -- the rule the design doc's own
    /// "Open Question" section resolved, matching `majority_read_strand`'s existing convention
    /// (`denovo_pipeline.rs`) and `build_footprint_seq`'s documented `+`-on-tie placeholder caveat.
    pub strand: char,
    pub n_supporting_reads: usize,
    pub read_names: Vec<String>,
    pub nearest_copy_tid: String,
    /// `None` when this family has NO copy on the candidate's own chromosome (reachable for a
    /// cross-chromosome family) -- written as `NA`, never as a `u64::MAX` sentinel.
    pub nearest_copy_distance: Option<u64>,
}

/// One max-AS placement of one AS-tied read, carried as ALIGNED BLOCKS rather than a bounding span.
///
/// ⚠ The blocks, not a `(ref_start, ref_end)` pair, are the load-bearing part. A bounding span counts a
/// spliced-out intron (`N`) as if the read covered it, which made `inside_any_copy` call a placement
/// "inside" a catalog copy that only an intron spans (no aligned base ever landing there) and inflated
/// reported cluster widths by whole intron lengths. This is the same bug class `block_overlap`
/// (`src/bin/copy_assign.rs`) was introduced to fix at the O3 truth gate.
#[derive(Clone, Debug, PartialEq)]
pub struct TiePlacement {
    pub chrom: String,
    /// One `(start, end)` per `M`/`=`/`X` run, in reference order. `D`/`N` advance the reference position
    /// without producing a block (exactly `block_overlap`'s own accounting).
    pub blocks: Vec<(u64, u64)>,
    /// SAM FLAG 0x10 of the record this placement came from -- the strand majority vote's input.
    pub reverse: bool,
}

/// Every `M`/`=`/`X` run of an alignment as its own reference interval `[pos, pos + n)`. The copy_assign
/// binary's `block_overlap` (pysam `get_blocks()` semantics) is computed from these blocks.
/// `D` and `N` advance the reference cursor without emitting a block, so a deletion splits a run here too
/// (matching `block_overlap`) -- usually immaterial since alignment `D` runs are short, but a `D` longer
/// than `TIE_PARTNER_MERGE_DISTANCE_BP` would split one placement into two reported clusters, same as a
/// real gap between two placements. Not observed on the real substrate this module was built against.
pub fn aligned_blocks(read: &AlignedRead) -> Vec<(u64, u64)> {
    let mut pos = read.ref_start;
    let mut out = Vec::new();
    for &(op, n) in &read.cigar {
        match op {
            'M' | '=' | 'X' => {
                out.push((pos, pos + n));
                pos += n;
            }
            'D' | 'N' => pos += n,
            _ => {}
        }
    }
    out
}

/// Does ANY aligned block of this placement overlap a copy already in the catalog?
///
/// Per-block, not a bounding-box test: a read whose intron merely SPANS a catalog copy has no aligned
/// base inside it and is NOT "inside" that copy.
fn inside_any_copy(existing_copies: &[(String, u64, u64, String)], p: &TiePlacement) -> bool {
    existing_copies.iter().any(|(c_chrom, c_start, c_end, _)| {
        *c_chrom == p.chrom && p.blocks.iter().any(|(b_start, b_end)| b_start < c_end && b_end > c_start)
    })
}

/// Nearest catalog copy of this family ON THE SAME CHROMOSOME, and its distance in bp (`0` when the
/// candidate overlaps it). `("NA", None)` when the family has no copy on that chromosome at all.
fn nearest_copy(existing_copies: &[(String, u64, u64, String)], chrom: &str, start: u64, end: u64) -> (String, Option<u64>) {
    existing_copies
        .iter()
        .filter(|(c_chrom, ..)| c_chrom == chrom)
        .map(|(_, c_start, c_end, tid)| {
            let d = if end <= *c_start { c_start - end } else if start >= *c_end { start - c_end } else { 0 };
            (tid.clone(), Some(d))
        })
        .min_by_key(|(_, d)| *d)
        .unwrap_or(("NA".to_string(), None))
}

/// Cluster out-of-catalog AS-tie-partner positions into candidate copies for one family.
///
/// `tied_reads`: `(read_name, that read's max-AS placements)` -- ALREADY RESTRICTED BY THE CALLER to the
/// reads this family actually considered (`discover_copies_for_family`, `src/bin/copy_assign.rs`). Passing
/// a whole region's tied reads to every family is the cross-family pooling bug the final whole-branch
/// review caught: it attributed one identical site, with an identical read list, to 3-8 different
/// `family_id`s.
/// `existing_copies`: `(chrom, start, end, tid)` for every copy already catalogued in this family (the
/// caller zips `FamilyAssignment::copy_spans` with `copy_tids`).
/// Positions already inside an existing copy span are defensively re-excluded here (never surfaced), even
/// though the caller is expected to have filtered them out already.
///
/// Pure: no BAM, no catalog types, no I/O.
pub fn cluster_tie_partners(
    tied_reads: &[(String, Vec<TiePlacement>)],
    family_id: &str,
    existing_copies: &[(String, u64, u64, String)],
    merge_distance: u64,
    min_support: usize,
) -> Vec<DiscoveredCopy> {
    // Flatten to one site per ALIGNED BLOCK of every placement that lands outside every existing copy.
    // `pid` identifies the contributing placement so the strand vote counts each placement once, not once
    // per exon block. Exclusion is per-PLACEMENT (a placement with any block inside a catalog copy is
    // dropped whole); clustering is per-BLOCK, so a spliced read's intron never chains two real loci into
    // one candidate.
    struct Site {
        pid: usize,
        name: String,
        chrom: String,
        start: u64,
        end: u64,
        reverse: bool,
    }
    let mut sites: Vec<Site> = Vec::new();
    let mut pid = 0usize;
    for (name, placements) in tied_reads {
        for p in placements {
            if !inside_any_copy(existing_copies, p) {
                for &(start, end) in &p.blocks {
                    sites.push(Site { pid, name: name.clone(), chrom: p.chrom.clone(), start, end, reverse: p.reverse });
                }
            }
            pid += 1;
        }
    }
    // Sort by (chrom, start, end) so overlap/proximity clustering is a single linear pass. `sort_by` is
    // stable, so equal keys keep `tied_reads`' own (already deterministic) order.
    sites.sort_by(|a, b| (a.chrom.as_str(), a.start, a.end).cmp(&(b.chrom.as_str(), b.start, b.end)));

    struct Cluster {
        chrom: String,
        start: u64,
        end: u64,
        read_names: Vec<String>,
        seen_names: HashSet<String>,
        voted: HashSet<usize>,
        n_fwd: usize,
        n_rev: usize,
    }
    let mut clusters: Vec<Cluster> = Vec::new();
    for s in sites {
        let merge = match clusters.last() {
            Some(last) => last.chrom == s.chrom && s.start <= last.end + merge_distance,
            None => false,
        };
        if !merge {
            clusters.push(Cluster {
                chrom: s.chrom.clone(),
                start: s.start,
                end: s.end,
                read_names: Vec::new(),
                seen_names: HashSet::new(),
                voted: HashSet::new(),
                n_fwd: 0,
                n_rev: 0,
            });
        }
        let c = clusters.last_mut().expect("just pushed or matched");
        c.end = c.end.max(s.end);
        // `seen_names` keeps membership O(1) while `read_names` keeps insertion order for the report.
        if c.seen_names.insert(s.name.clone()) {
            c.read_names.push(s.name);
        }
        if c.voted.insert(s.pid) {
            if s.reverse {
                c.n_rev += 1;
            } else {
                c.n_fwd += 1;
            }
        }
    }

    clusters
        .into_iter()
        .filter(|c| c.read_names.len() >= min_support)
        .map(|c| {
            let (nearest_copy_tid, nearest_copy_distance) = nearest_copy(existing_copies, &c.chrom, c.start, c.end);
            DiscoveredCopy {
                family_id: family_id.to_string(),
                chrom: c.chrom,
                start: c.start,
                end: c.end,
                // Majority vote; an exact tie defaults to `+` (design doc's Open Question, resolved).
                strand: if c.n_rev > c.n_fwd { '-' } else { '+' },
                n_supporting_reads: c.read_names.len(),
                read_names: c.read_names,
                nearest_copy_tid,
                nearest_copy_distance,
            }
        })
        .collect()
}

/// Every AS-tied read in the region, with all of its own max-AS placements as aligned blocks.
///
/// A read is AS-tied here iff at least two of its non-supplementary records share its maximum `AS`.
pub fn tie_partner_placements(bam_reads: &[BamRead]) -> Vec<(String, Vec<TiePlacement>)> {
    let mut by_name: HashMap<&str, Vec<&BamRead>> = HashMap::new();
    for br in bam_reads.iter().filter(|b| !b.is_supplementary) {
        by_name.entry(br.name.as_str()).or_default().push(br);
    }
    let mut out = Vec::new();
    for (name, placements) in by_name {
        if placements.len() < 2 {
            continue;
        }
        let max_as = placements.iter().map(|b| b.as_score).max().unwrap();
        let tied: Vec<&BamRead> = placements.iter().copied().filter(|b| b.as_score == max_as).collect();
        if tied.len() >= 2 {
            out.push((
                name.to_string(),
                tied.into_iter()
                    .map(|b| TiePlacement {
                        chrom: b.chrom.clone(),
                        blocks: aligned_blocks(&b.read),
                        reverse: b.reverse,
                    })
                    .collect(),
            ));
        }
    }
    // Sort by read name for deterministic output order (HashMap iteration is randomized).
    out.sort_by(|a, b| a.0.cmp(&b.0));
    out
}

#[cfg(test)]
mod tests {
    use super::*;

    fn one_copy(tid: &str, chrom: &str, start: u64, end: u64) -> Vec<(String, u64, u64, String)> {
        vec![(chrom.to_string(), start, end, tid.to_string())]
    }

    /// A single-block (unspliced) placement on the forward strand.
    fn pl(chrom: &str, start: u64, end: u64) -> TiePlacement {
        TiePlacement { chrom: chrom.to_string(), blocks: vec![(start, end)], reverse: false }
    }

    fn pl_rev(chrom: &str, start: u64, end: u64) -> TiePlacement {
        TiePlacement { chrom: chrom.to_string(), blocks: vec![(start, end)], reverse: true }
    }

    #[test]
    fn defensively_excludes_positions_inside_a_catalog_copy() {
        let existing = one_copy("c0", "chr1", 1000, 2000);
        // this "tied" position sits INSIDE c0's span -- must never surface as a discovery
        let tied = vec![("read1".to_string(), vec![pl("chr1", 1200, 1300)])];
        let out = cluster_tie_partners(&tied, "FAM0", &existing, 500, 2);
        assert!(out.is_empty(), "a position already inside a catalog copy must never be reported");
    }

    #[test]
    fn merges_positions_within_merge_distance_and_respects_min_support() {
        let existing = one_copy("c0", "chr1", 1000, 2000);
        let tied = vec![
            ("read1".to_string(), vec![pl("chr1", 5000, 5100)]),
            ("read2".to_string(), vec![pl("chr1", 5050, 5150)]), // within 500bp of read1's site
        ];
        let out = cluster_tie_partners(&tied, "FAM0", &existing, 500, 2);
        assert_eq!(out.len(), 1, "two nearby out-of-catalog positions with 2 supporting reads = 1 cluster");
        assert_eq!(out[0].n_supporting_reads, 2);
        assert_eq!(out[0].nearest_copy_tid, "c0");
        assert_eq!(out[0].nearest_copy_distance, Some(3000)); // 5000 - 2000
    }

    #[test]
    fn keeps_clusters_separate_beyond_merge_distance() {
        let existing = one_copy("c0", "chr1", 1000, 2000);
        let tied = vec![
            ("read1".to_string(), vec![pl("chr1", 5000, 5100)]),
            ("read2".to_string(), vec![pl("chr1", 5100, 5100)]),
            ("read3".to_string(), vec![pl("chr1", 9000, 9100)]),
            ("read4".to_string(), vec![pl("chr1", 9050, 9150)]),
        ];
        let out = cluster_tie_partners(&tied, "FAM0", &existing, 500, 2);
        assert_eq!(out.len(), 2, "two far-apart pairs must stay two separate clusters");
    }

    #[test]
    fn drops_clusters_below_min_support() {
        let existing = one_copy("c0", "chr1", 1000, 2000);
        let tied = vec![("read1".to_string(), vec![pl("chr1", 5000, 5100)])];
        let out = cluster_tie_partners(&tied, "FAM0", &existing, 500, 2);
        assert!(out.is_empty(), "a single supporting read must not clear min_support=2");
    }

    #[test]
    fn a_copy_inside_a_spliced_out_intron_does_not_exclude_the_placement() {
        // Mirrors `block_overlap_ignores_an_intron_that_merely_spans_the_window`
        // (src/bin/copy_assign.rs) at the containment layer. ref_start=1000, "10M5000N10M": aligned blocks
        // are [1000,1010) and [6010,6020); the intron covers [1010,6010) with no aligned base in it. A
        // catalog copy sitting entirely inside that intron must NOT swallow the placement -- with the old
        // bounding-span test (`ref_start..ref_end` = 1000..6020) it did, silently killing the candidate.
        let spliced = TiePlacement {
            chrom: "chr1".to_string(),
            blocks: aligned_blocks(&AlignedRead {
                ref_start: 1000,
                cigar: vec![('M', 10), ('N', 5000), ('M', 10)],
                seq: vec![],
                qual: vec![],
            }),
            reverse: false,
        };
        assert_eq!(spliced.blocks, vec![(1000, 1010), (6010, 6020)]);
        let copy_in_the_intron = one_copy("c0", "chr1", 3000, 3100);
        let tied = vec![
            ("read1".to_string(), vec![spliced.clone()]),
            ("read2".to_string(), vec![spliced.clone()]),
        ];
        let out = cluster_tie_partners(&tied, "FAM0", &copy_in_the_intron, 500, 2);
        assert_eq!(out.len(), 2, "both aligned blocks survive as their own candidates; neither is 'inside' c0");
        assert_eq!((out[0].start, out[0].end), (1000, 1010));
        assert_eq!((out[1].start, out[1].end), (6010, 6020), "the 5kb intron never chains the two blocks");

        // Contrast: a copy that genuinely overlaps one of the ALIGNED blocks does exclude the placement.
        let copy_on_the_block = one_copy("c0", "chr1", 1005, 1100);
        let out2 = cluster_tie_partners(&tied, "FAM0", &copy_on_the_block, 500, 2);
        assert!(out2.is_empty(), "a real aligned-base overlap must still exclude the whole placement");
    }

    #[test]
    fn strand_is_the_majority_vote_and_a_tie_defaults_to_plus() {
        let existing = one_copy("c0", "chr1", 1000, 2000);
        // 2 forward + 1 reverse -> '+'
        let fwd_majority = vec![
            ("read1".to_string(), vec![pl("chr1", 5000, 5100)]),
            ("read2".to_string(), vec![pl("chr1", 5050, 5150)]),
            ("read3".to_string(), vec![pl_rev("chr1", 5060, 5160)]),
        ];
        let out = cluster_tie_partners(&fwd_majority, "FAM0", &existing, 500, 2);
        assert_eq!(out.len(), 1);
        assert_eq!(out[0].strand, '+', "2 forward vs 1 reverse");

        // 2 reverse + 1 forward -> '-'
        let rev_majority = vec![
            ("read1".to_string(), vec![pl_rev("chr1", 5000, 5100)]),
            ("read2".to_string(), vec![pl_rev("chr1", 5050, 5150)]),
            ("read3".to_string(), vec![pl("chr1", 5060, 5160)]),
        ];
        let out = cluster_tie_partners(&rev_majority, "FAM0", &existing, 500, 2);
        assert_eq!(out.len(), 1);
        assert_eq!(out[0].strand, '-', "2 reverse vs 1 forward");

        // exact tie -> '+' by the documented default
        let tie = vec![
            ("read1".to_string(), vec![pl("chr1", 5000, 5100)]),
            ("read2".to_string(), vec![pl_rev("chr1", 5050, 5150)]),
        ];
        let out = cluster_tie_partners(&tie, "FAM0", &existing, 500, 2);
        assert_eq!(out.len(), 1);
        assert_eq!(out[0].strand, '+', "an exact strand tie defaults to '+'");
    }

    #[test]
    fn a_multi_exon_placement_votes_once_for_strand() {
        // One reverse placement with 3 exon blocks close enough to land in ONE cluster must not outvote
        // two forward single-block reads: the vote is per PLACEMENT, not per block.
        let existing = one_copy("c0", "chr1", 1000, 2000);
        let three_exons = TiePlacement {
            chrom: "chr1".to_string(),
            blocks: vec![(5000, 5050), (5100, 5150), (5200, 5250)],
            reverse: true,
        };
        let tied = vec![
            ("read1".to_string(), vec![pl("chr1", 5000, 5100)]),
            ("read2".to_string(), vec![pl("chr1", 5050, 5150)]),
            ("read3".to_string(), vec![three_exons]),
        ];
        let out = cluster_tie_partners(&tied, "FAM0", &existing, 500, 2);
        assert_eq!(out.len(), 1);
        assert_eq!(out[0].n_supporting_reads, 3, "distinct read NAMES, not blocks");
        assert_eq!(out[0].strand, '+', "2 forward placements outvote 1 reverse placement's 3 blocks");
    }

    #[test]
    fn no_copy_on_this_chromosome_yields_na_not_a_sentinel() {
        // A cross-chromosome family: the candidate lands on chr2, every catalog copy is on chr1.
        let existing = one_copy("c0", "chr1", 1000, 2000);
        let tied = vec![
            ("read1".to_string(), vec![pl("chr2", 5000, 5100)]),
            ("read2".to_string(), vec![pl("chr2", 5050, 5150)]),
        ];
        let out = cluster_tie_partners(&tied, "FAM0", &existing, 500, 2);
        assert_eq!(out.len(), 1);
        assert_eq!(out[0].nearest_copy_tid, "NA");
        assert_eq!(out[0].nearest_copy_distance, None, "no copy on chr2 -> NA, never u64::MAX");
    }

    #[test]
    fn tie_partner_placements_finds_reads_tied_at_their_own_max_as() {
        use crate::vg_family::denovo_assemble::BamRead;
        use crate::vg_family::copy_split::AlignedRead;
        let mk = |name: &str, chrom: &str, start: u64, as_score: i32| BamRead {
            chrom: chrom.into(),
            read: AlignedRead { ref_start: start, cigar: vec![('M', 100)], seq: vec![], qual: vec![] },
            mapq: 0, name: name.into(), as_score, de: 0.0,
            is_supplementary: false, is_secondary: as_score != 200, reverse: false, ts: None,
        };
        let reads = vec![
            mk("tied_read", "chr1", 1000, 200),   // best
            mk("tied_read", "chr1", 5000, 200),   // tied with the above
            mk("tied_read", "chr1", 9000, 150),   // worse, not part of the tie
            mk("unique_read", "chr1", 2000, 300), // only one placement, never tied
        ];
        let out = tie_partner_placements(&reads);
        assert_eq!(out.len(), 1);
        assert_eq!(out[0].0, "tied_read");
        assert_eq!(out[0].1.len(), 2, "only the 2 max-scoring placements, not the 150-scoring one");
        assert_eq!(out[0].1[0].blocks, vec![(1000, 1100)], "placements carry aligned blocks, not a span");
    }

    #[test]
    fn tie_partner_placements_carries_the_strand_flag_and_splits_on_introns() {
        use crate::vg_family::denovo_assemble::BamRead;
        use crate::vg_family::copy_split::AlignedRead;
        let mk = |start: u64, reverse: bool| BamRead {
            chrom: "chr1".into(),
            read: AlignedRead { ref_start: start, cigar: vec![('M', 10), ('N', 500), ('M', 10)], seq: vec![], qual: vec![] },
            mapq: 0, name: "r".into(), as_score: 200, de: 0.0,
            is_supplementary: false, is_secondary: false, reverse, ts: None,
        };
        let out = tie_partner_placements(&[mk(100, false), mk(9000, true)]);
        assert_eq!(out.len(), 1);
        let p = &out[0].1;
        assert_eq!(p.len(), 2);
        assert_eq!(p[0].blocks, vec![(100, 110), (610, 620)], "the intron is not part of any block");
        assert!(!p[0].reverse);
        assert!(p[1].reverse, "FLAG 0x10 is carried through per placement");
    }
}
}
