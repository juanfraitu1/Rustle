# Copy-Graph Objects (v1) Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Make every copy of a family a tagged, corroborable path in one Bandage-loadable family variation graph emitted by `copy_assign --phase`, with a REFERENCE walk so a reference-absent copy is visibly a path the reference does not take.

**Architecture:** A new pure builder module `src/rustle/vg_family/copy_graph.rs` owns graph construction and GFA/colour serialization from plain parallel vectors (columns, per-copy alleles, per-read observations). `copy_assign --phase` fills those vectors from the existing `FamilyAssignment` + `GenomeIndex` and calls the builder, replacing the inline emitter. The builder has zero I/O and zero domain dependencies, so it is fully unit-testable in memory.

**Tech Stack:** Rust (crate `rustle`), existing `rustle::genome::GenomeIndex` (`fetch_sequence`), GFA 1.1, Bandage.

## Global Constraints

- Default `copy_assign` behaviour (no `--phase`) MUST be byte-identical to today (the builder runs only inside the existing `if args.phase` block).
- No new k-mer computation anywhere in this feature.
- Output MUST be valid GFA 1.1 and open in Bandage; EVERY P-line and W-line walk step MUST be backed by an `L` line (no dangling edges).
- Honesty rule: a corroboration tag is OMITTED when its value is unknown; never write a fabricated value.
- v2 is OUT OF SCOPE: exon presence/absence bubbles, ORF/protein (`OR`) tag, real per-transcript backbone sequence.
- Existing `--phase` TSV outputs (`phase_block_lines`, `phased_hap_lines`, `phased_read_lines`) MUST stay unchanged; only the GFA/colour emission is replaced.

## File Structure

- **Create** `src/rustle/vg_family/copy_graph.rs` — the pure builder: types (`CopyStatus`, `Corrob`, `PsvColumn`, `CopyPath`, `ReadWalk`, `CopyGraph`) + methods (`to_gfa`, `gfa_lines`, `colours_csv`, `legend_tsv`) + unit tests.
- **Modify** `src/rustle/vg_family/mod.rs` — add `pub mod copy_graph;`.
- **Modify** `src/bin/copy_assign.rs` — inside the `if args.phase` block, replace the inline GFA building (approx lines 946–994) and the GFA writer (approx lines 1176–1210) with a call into `copy_graph`; add the optional `--gff` arg + annotation tagging.

Node id scheme (all ids sanitized whitespace-free; family prefix keeps ids unique across families in one file):
- backbone spacer node `i`: `{fam}_bb{i}` (i in `0..=M`, where `M = columns.len()`)
- allele node at graph-column `ci` with base `b`: `{fam}_c{ci}_{b}`

**Parallelism contract (critical):** `columns`, every `CopyPath.alleles`, and every `ReadWalk.obs` are all length `M` and share the same column order. `backbone` is length `M+1`. The caller (Task 7) is responsible for filtering out unusable columns (missing genome position or reference allele) and building these parallel vectors; the builder assumes they are already aligned.

---

### Task 1: Module scaffold + types

**Files:**
- Create: `src/rustle/vg_family/copy_graph.rs`
- Modify: `src/rustle/vg_family/mod.rs` (add `pub mod copy_graph;` alongside the other `pub mod` lines)
- Test: in `src/rustle/vg_family/copy_graph.rs` (`#[cfg(test)] mod tests`)

**Interfaces:**
- Produces: the public types below, used by all later tasks and by Task 7's wiring.

- [ ] **Step 1: Write the failing test**

Add to the bottom of the new file:

```rust
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
}
```

- [ ] **Step 2: Run test to verify it fails**

Run: `cargo test -p rustle copy_graph::tests::constructs_and_reports_shape 2>&1 | tail -20`
Expected: FAIL — `copy_graph` module / types not found (compile error).

- [ ] **Step 3: Write minimal implementation**

At the top of `src/rustle/vg_family/copy_graph.rs`:

```rust
//! Copy-graph objects (v1): every copy of a family is a tagged, corroborable PATH in one GFA 1.1
//! variation graph. A REFERENCE walk makes a reference-absent copy visibly an arm the reference does
//! not take. Pure builder — no I/O; the caller fills the parallel vectors and writes the strings.

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
```

And add `pub mod copy_graph;` to `src/rustle/vg_family/mod.rs`.

- [ ] **Step 4: Run test to verify it passes**

Run: `cargo test -p rustle copy_graph::tests::constructs_and_reports_shape 2>&1 | tail -20`
Expected: PASS (1 passed).

- [ ] **Step 5: Commit**

```bash
git add src/rustle/vg_family/copy_graph.rs src/rustle/vg_family/mod.rs
git commit -m "feat(copy_graph): scaffold copy-graph object types"
```

---

### Task 2: GFA skeleton — backbone spacers + allele bubbles + L-lines

**Files:**
- Modify: `src/rustle/vg_family/copy_graph.rs`
- Test: same file's `tests` module

**Interfaces:**
- Produces: `CopyGraph::gfa_lines(&self) -> GfaLines` and `CopyGraph::to_gfa(&self) -> String`. `GfaLines { header, segs, links, paths, walks }` (all `Vec<String>` except header `String`). This task fills `header`, `segs`, `links`; `paths`/`walks` come in Tasks 3–5.

- [ ] **Step 1: Write the failing test**

Add to `tests`:

```rust
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
```

- [ ] **Step 2: Run test to verify it fails**

Run: `cargo test -p rustle copy_graph::tests::skeleton_has_backbone_alleles_and_links 2>&1 | tail -20`
Expected: FAIL — `to_gfa` / `gfa_lines` not found.

- [ ] **Step 3: Write minimal implementation**

Add to `copy_graph.rs` (an `impl CopyGraph` block):

```rust
use std::collections::BTreeSet;

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
```

- [ ] **Step 4: Run test to verify it passes**

Run: `cargo test -p rustle copy_graph::tests::skeleton_has_backbone_alleles_and_links 2>&1 | tail -20`
Expected: PASS.

- [ ] **Step 5: Commit**

```bash
git add src/rustle/vg_family/copy_graph.rs
git commit -m "feat(copy_graph): GFA skeleton — backbone spacers, allele bubbles, backing L-lines"
```

---

### Task 3: REFERENCE walk P-line

**Files:**
- Modify: `src/rustle/vg_family/copy_graph.rs`
- Test: same file's `tests` module

**Interfaces:**
- Consumes: `gfa_lines` from Task 2.
- Produces: a REFERENCE P-line `{fam}_REFERENCE` in `GfaLines.paths`, walking `bb0, c{ci}_{ref}, ..., bbM`, tagged `ST:Z:reference`.

- [ ] **Step 1: Write the failing test**

```rust
    #[test]
    fn reference_walk_threads_reference_alleles() {
        let g = tiny_graph();
        let gfa = g.to_gfa();
        // reference alleles are A (col0) and C (col1)
        assert!(gfa.contains(
            "P\tFAM1_REFERENCE\tFAM1_bb0+,FAM1_c0_A+,FAM1_bb1+,FAM1_c1_C+,FAM1_bb2+\t*\tST:Z:reference"
        ), "reference P-line missing or wrong:\n{}", gfa);
    }
```

- [ ] **Step 2: Run test to verify it fails**

Run: `cargo test -p rustle copy_graph::tests::reference_walk_threads_reference_alleles 2>&1 | tail -20`
Expected: FAIL — no `FAM1_REFERENCE` path.

- [ ] **Step 3: Write minimal implementation**

Add a helper and call it inside `gfa_lines` (after the allele-node loop, before returning `out`):

```rust
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
```

Then in `gfa_lines`, before `out` is returned, push the reference path:

```rust
        // REFERENCE walk: the genome's own allele at each column.
        let ref_taken: Vec<Option<u8>> = self.columns.iter().map(|c| c.ref_allele).collect();
        let ref_walk = self.walk_tokens(&ref_taken).join(",");
        out.paths.push(format!("P\t{}_REFERENCE\t{}\t*\tST:Z:reference", self.family, ref_walk));
```

- [ ] **Step 4: Run test to verify it passes**

Run: `cargo test -p rustle copy_graph::tests::reference_walk_threads_reference_alleles 2>&1 | tail -20`
Expected: PASS. Also re-run Task 2's test to confirm no regression: `cargo test -p rustle copy_graph:: 2>&1 | tail -20`.

- [ ] **Step 5: Commit**

```bash
git add src/rustle/vg_family/copy_graph.rs
git commit -m "feat(copy_graph): emit the REFERENCE walk P-line"
```

---

### Task 4: Copy P-lines with RC/SU/MI/ST tags + _ABSENT naming

**Files:**
- Modify: `src/rustle/vg_family/copy_graph.rs`
- Test: same file's `tests` module

**Interfaces:**
- Consumes: `walk_tokens` (Task 3), `Corrob`, `CopyStatus`.
- Produces: one P-line per copy in `GfaLines.paths`, name `{fam}_copy{ci}` (suffix `_ABSENT` when `status.is_absent()`), tags `RC:i: SU:i: MI:f: ST:Z:` (each omitted when its value is unknown). `SU` is computed here (private-column count) when `corrob.suns` is `None`.

- [ ] **Step 1: Write the failing test**

```rust
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
```

- [ ] **Step 2: Run test to verify it fails**

Run: `cargo test -p rustle copy_graph::tests::copy_paths_carry_tags_and_absent_diverges 2>&1 | tail -20`
Expected: FAIL — no copy P-lines emitted.

- [ ] **Step 3: Write minimal implementation**

Add a private-column helper:

```rust
    /// Number of columns where copy `c`'s (Some) allele is unique among the family's copies (reference excluded).
    fn private_columns(&self, c: usize) -> u32 {
        let mut n = 0u32;
        for ci in 0..self.m() {
            let Some(Some(b)) = self.copies[c].alleles.get(ci) else { continue };
            let unique = self.copies.iter().enumerate()
                .all(|(k, other)| k == c || other.alleles.get(ci).and_then(|o| *o) != Some(*b));
            if unique { n += 1; }
        }
        n
    }
```

Then in `gfa_lines`, after the reference path push, emit copy paths:

```rust
        // copy P-lines with corroboration tags
        for (ci, cp) in self.copies.iter().enumerate() {
            let walk = self.walk_tokens(&cp.alleles).join(",");
            let name = if cp.status.is_absent() {
                format!("{}_copy{}_ABSENT", self.family, ci)
            } else {
                format!("{}_copy{}", self.family, ci)
            };
            let mut tags = String::new();
            if let Some(rc) = cp.corrob.reads { tags.push_str(&format!("\tRC:i:{}", rc)); }
            let su = cp.corrob.suns.unwrap_or_else(|| self.private_columns(ci));
            tags.push_str(&format!("\tSU:i:{}", su));
            if let Some(mi) = cp.corrob.map_identity { tags.push_str(&format!("\tMI:f:{:.3}", mi)); }
            tags.push_str(&format!("\tST:Z:{}", cp.status.tag()));
            out.paths.push(format!("P\t{}\t{}\t*{}", name, walk, tags));
        }
```

- [ ] **Step 4: Run test to verify it passes**

Run: `cargo test -p rustle copy_graph::tests::copy_paths_carry_tags_and_absent_diverges 2>&1 | tail -20`
Expected: PASS.

- [ ] **Step 5: Commit**

```bash
git add src/rustle/vg_family/copy_graph.rs
git commit -m "feat(copy_graph): copy P-lines with RC/SU/MI/ST tags and _ABSENT naming"
```

---

### Task 5: Read W-lines with backing L-lines (dangling-free)

**Files:**
- Modify: `src/rustle/vg_family/copy_graph.rs`
- Test: same file's `tests` module

**Interfaces:**
- Consumes: `walk_tokens`, `alleles_at`.
- Produces: one W-line per read that observed ≥1 column, spanning `[first_obs ..= last_obs]`, routing gaps within the span through the reference node. Because every walk step is a backbone↔allele adjacency for which Task 2 already emitted an L-line, no walk step is dangling.

- [ ] **Step 1: Write the failing test**

```rust
    // Assert every P-line and W-line step is backed by an L-line (parses walks, checks adjacency set).
    fn assert_no_dangling(gfa: &str) {
        use std::collections::HashSet;
        let mut links: HashSet<(String, String)> = HashSet::new();
        for l in gfa.lines().filter(|l| l.starts_with("L\t")) {
            let f: Vec<&str> = l.split('\t').collect(); // L from + to + 0M
            links.insert((f[1].to_string(), f[3].to_string()));
        }
        let node = |tok: &str| tok.trim_start_matches(['>', '<']).trim_end_matches('+').trim_end_matches('-').to_string();
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
```

- [ ] **Step 2: Run test to verify it fails**

Run: `cargo test -p rustle copy_graph::tests::reads_walk_with_backing_links 2>&1 | tail -20`
Expected: FAIL — no W-lines emitted.

- [ ] **Step 3: Write minimal implementation**

In `gfa_lines`, after the copy paths, emit read walks:

```rust
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
```

Note: the read's observed allele node is always in `alleles_at(ci)` (reads contribute to that set in Task 2), and the reference fallback node is too — so both flanking L-lines already exist. No dangling by construction.

- [ ] **Step 4: Run test to verify it passes**

Run: `cargo test -p rustle copy_graph::tests::reads_walk_with_backing_links 2>&1 | tail -20`
Expected: PASS. Re-run the whole module: `cargo test -p rustle copy_graph:: 2>&1 | tail -20` (all pass).

- [ ] **Step 5: Commit**

```bash
git add src/rustle/vg_family/copy_graph.rs
git commit -m "feat(copy_graph): read W-lines with backing L-lines (no dangling edges)"
```

---

### Task 6: colours.csv (segment-keyed) + legend.tsv

**Files:**
- Modify: `src/rustle/vg_family/copy_graph.rs`
- Test: same file's `tests` module

**Interfaces:**
- Produces: `CopyGraph::colours_csv(&self) -> String` (rows `segment,#RRGGBB`) and `CopyGraph::legend_tsv(&self) -> String` (rows `status\t#RRGGBB`). Colour rule: a non-reference allele node walked by any absent copy => red; a node on the reference walk => grey; any other non-reference allele node => the walking copy's status colour; backbone nodes => light grey.

- [ ] **Step 1: Write the failing test**

```rust
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
```

- [ ] **Step 2: Run test to verify it fails**

Run: `cargo test -p rustle copy_graph::tests::colours_mark_absent_red_reference_grey 2>&1 | tail -20`
Expected: FAIL — `colours_csv` not found.

- [ ] **Step 3: Write minimal implementation**

```rust
impl CopyGraph {
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
                // walked by any absent copy? (and not the reference allele)
                let by_absent = self.copies.iter()
                    .any(|c| c.status.is_absent() && c.alleles.get(ci).and_then(|o| *o) == Some(b));
                if by_absent {
                    colour.insert(nid, CopyStatus::AbsentDivergent.colour());
                } else if let Some(c) = self.copies.iter()
                    .find(|c| c.alleles.get(ci).and_then(|o| *o) == Some(b)) {
                    colour.insert(nid, c.status.colour());
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
```

- [ ] **Step 4: Run test to verify it passes**

Run: `cargo test -p rustle copy_graph::tests::colours_mark_absent_red_reference_grey 2>&1 | tail -20`
Expected: PASS. Full module green: `cargo test -p rustle copy_graph:: 2>&1 | tail -20`.

- [ ] **Step 5: Commit**

```bash
git add src/rustle/vg_family/copy_graph.rs
git commit -m "feat(copy_graph): segment-keyed colours.csv + legend.tsv"
```

---

### Task 7: Wire into `copy_assign --phase`

**Files:**
- Modify: `src/bin/copy_assign.rs` (inside the `if args.phase` block, approx lines 946–994 and the GFA writer approx 1176–1210)
- Test: `src/bin/copy_assign.rs` `#[cfg(test)] mod tests` (a wiring helper unit test) + a manual byte-identical check

**Interfaces:**
- Consumes: `rustle::vg_family::copy_graph::{CopyGraph, CopyPath, ReadWalk, PsvColumn, Corrob, CopyStatus, GfaLines}`, `FamilyAssignment`, `rustle::genome::GenomeIndex::fetch_sequence`.
- Produces: `<out>.phase.gfa`, `<out>.phase.gfa.colours.csv`, `<out>.phase.gfa.legend.tsv` built from the per-family `CopyGraph`s.

**Design of the wiring:** factor a pure helper `fn build_copy_graph(fid: &str, fa: &FamilyAssignment, ref_base: impl Fn(&str, u64) -> Option<u8>, bam_reads: &[BamRead]) -> CopyGraph` so it is unit-testable with an injected `ref_base` closure (no genome I/O in the test). In `main`, call it with `ref_base = |chrom, pos| genome_for(chrom).ok().and_then(|g| g.fetch_sequence(chrom, pos, pos+1)).and_then(|v| v.first().copied())`. Accumulate each family's `gfa_lines()`/`colours_csv()`/`legend_tsv()` into the file writers (dedup `segs`/`links` across families with the existing `HashSet` approach). Keep `phase_block_lines`/`phased_hap_lines`/`phased_read_lines` exactly as they are.

**Status derivation (in `build_copy_graph`):** a copy `ci` is absent when any of its assigned reads is `discovery_coupled`, OR `ci` is among the collapsed/rescued tail (`ci >= fa.copy_tids.len() - (fa.collapsed_copies + fa.rescued_copies)`); subtype = `AbsentCollapsed` if the copy came from within-locus collapse, else `AbsentDivergent` (v1: default absent copies to `AbsentDivergent`, upgrade to `AbsentCollapsed` when `fa.collapsed_copies > 0` and the copy index is in the collapsed tail). In-genome copies are `AnnotationUnknown` (Task 8 refines to annotated/unannotated when `--gff` is given). `MI` comes from the copy's remap identity if the family carries one, else `None`.

- [ ] **Step 1: Write the failing test**

Add to `copy_assign.rs` tests (build a minimal `FamilyAssignment` via `FamilyAssignment::empty()` then set fields):

```rust
    #[test]
    fn build_copy_graph_maps_family_to_graph() {
        use rustle::vg_family::denovo_pipeline::FamilyAssignment;
        use rustle::vg_family::copy_graph::CopyStatus;
        let mut fa = FamilyAssignment::empty();
        fa.chrom = "chr1".into();
        fa.n_copies = 2;
        fa.copy_tids = vec!["c0".into(), "c1".into()];
        fa.psv_col_pos = vec![Some(100), Some(200)];
        fa.copy_psv_alleles = vec![vec![Some(b'A'), Some(b'A')], vec![Some(b'A'), Some(b'G')]];
        fa.read_psv_obs = vec![];
        fa.assignments = vec![];
        // injected reference: A at both positions
        let ref_base = |_c: &str, _p: u64| Some(b'A');
        let g = build_copy_graph("CAFAM0", &fa, ref_base, &[]);
        assert_eq!(g.columns.len(), 2);
        assert_eq!(g.copies.len(), 2);
        assert_eq!(g.columns[0].ref_allele, Some(b'A'));
        // graph renders the reference walk + a divergent copy
        let gfa = g.to_gfa();
        assert!(gfa.contains("P\tCAFAM0_REFERENCE"));
        assert!(gfa.contains("CAFAM0_c1_G+"));
    }
```

- [ ] **Step 2: Run test to verify it fails**

Run: `cargo test -p rustle --bin copy_assign build_copy_graph_maps_family_to_graph 2>&1 | tail -20`
Expected: FAIL — `build_copy_graph` not defined.

- [ ] **Step 3: Write minimal implementation**

Add `build_copy_graph` to `copy_assign.rs` (free function). Filter to usable columns (genome_pos + ref_base both Some), build parallel vectors, derive status:

```rust
fn build_copy_graph(
    fid: &str,
    fa: &rustle::vg_family::denovo_pipeline::FamilyAssignment,
    ref_base: impl Fn(&str, u64) -> Option<u8>,
    bam_reads: &[BamRead],
) -> rustle::vg_family::copy_graph::CopyGraph {
    use rustle::vg_family::copy_graph::*;
    // usable columns (both a genome position and a reference base), remembering the original index
    let mut cols: Vec<PsvColumn> = Vec::new();
    let mut keep: Vec<usize> = Vec::new();
    for (j, p) in fa.psv_col_pos.iter().enumerate() {
        if let Some(pos) = p {
            if let Some(rb) = ref_base(&fa.chrom, *pos) {
                cols.push(PsvColumn { col: j, genome_pos: Some(*pos), ref_allele: Some(rb) });
                keep.push(j);
            }
        }
    }
    let sel = |row: &Vec<Option<u8>>| keep.iter().map(|&j| row.get(j).copied().flatten()).collect::<Vec<_>>();
    let backbone = vec![b"NNNNNNNNNN".to_vec(); cols.len() + 1];

    // which copies are absent (discovery_coupled reads) — collapsed tail => AbsentCollapsed else AbsentDivergent
    let n = fa.copy_tids.len();
    let absent_tail_start = n.saturating_sub(fa.collapsed_copies + fa.rescued_copies);
    let mut coupled = vec![false; n];
    for (_ri, a) in &fa.assignments {
        if a.discovery_coupled && a.best_copy < n { coupled[a.best_copy] = true; }
    }
    let copies = (0..n).map(|ci| {
        let is_absent = coupled[ci] || ci >= absent_tail_start && (fa.collapsed_copies + fa.rescued_copies) > 0;
        let status = if is_absent {
            if fa.collapsed_copies > 0 && ci >= absent_tail_start { CopyStatus::AbsentCollapsed }
            else { CopyStatus::AbsentDivergent }
        } else {
            CopyStatus::AnnotationUnknown
        };
        let reads = fa.assignments.iter()
            .filter(|(_, a)| a.best_copy == ci && a.status_is_assigned()).count() as u32;
        CopyPath {
            id: format!("{}_copy{}", fid, ci),
            alleles: fa.copy_psv_alleles.get(ci).map(|r| sel(r)).unwrap_or_default(),
            status,
            corrob: Corrob { reads: Some(reads), suns: None, map_identity: None },
        }
    }).collect();

    let reads = fa.assignments.iter().zip(fa.read_psv_obs.iter()).map(|((ri, a), obs)| {
        ReadWalk {
            name: sanitize_gfa_id(&format!("{}_{}", fid, bam_reads[*ri])),
            obs: sel(obs),
            assigned_copy: if a.status_is_assigned() { Some(a.best_copy) } else { None },
        }
    }).collect();

    CopyGraph { family: fid.to_string(), columns: cols, backbone, copies, reads }
}
```

If `Assignment` has no `status_is_assigned()` helper, inline `matches!(a.status, AssignStatus::Assigned)` (the type already imported in this file). Then **replace** the inline GFA-building blocks (the two `for` loops that push to `gfa_segs/gfa_links/gfa_paths/gfa_walks/gfa_colors`, approx lines 949–994) with:

```rust
                    let ref_base = |chrom: &str, pos: u64| {
                        genome_for(chrom).ok()
                            .and_then(|g| g.fetch_sequence(chrom, pos, pos + 1))
                            .and_then(|v| v.first().copied())
                    };
                    let cg = build_copy_graph(&fid, fa, ref_base, bam_reads);
                    let gl = cg.gfa_lines();
                    for s in gl.segs { gfa_segs.insert(s); }
                    for l in gl.links { gfa_links.insert(l); }   // gfa_links becomes HashSet<String>
                    gfa_paths.extend(gl.paths);
                    gfa_walks.extend(gl.walks);
                    for row in cg.colours_csv().lines() { gfa_colors.push(row.to_string()); }
                    legend_rows.extend(cg.legend_tsv().lines().map(|s| s.to_string()));
```

Adjust the declarations near line 576: change `gfa_links: HashSet<(String,String)>` to `HashSet<String>`, and add `let mut legend_rows: Vec<String> = Vec::new();`. Update the writer (approx 1176–1210): the L-line write becomes `for l in &links { writeln!(gf, "{}", l)?; }` (already full `L\t...` strings), and after writing `.colours.csv`, also write `.phase.gfa.legend.tsv` from de-duplicated `legend_rows`.

- [ ] **Step 4: Run test to verify it passes**

Run: `cargo test -p rustle --bin copy_assign build_copy_graph_maps_family_to_graph 2>&1 | tail -20`
Expected: PASS. Then build the binary: `cargo build -p rustle --bin copy_assign 2>&1 | tail -5` (Expected: Finished).

- [ ] **Step 5: Byte-identical guard (no `--phase`) + commit**

Confirm the non-phase path is untouched by diffing an assignments run before/after (small region, foreground, per the WSL2 rule):

Run:
```bash
B=/home/juanfra/winloci_scratch/GGO_mm.bam; F=/home/juanfra/winloci_scratch/GGO.fasta
R=NC_073229.2:49040000-49075000
./target/release/copy_assign --bam $B --fasta $F --region $R --out /home/juanfra/winloci_scratch/cg_nophase 2>/dev/null
md5sum /home/juanfra/winloci_scratch/cg_nophase.assignments.tsv
```
Expected: the md5 matches a run of the pre-change binary on the same region (record it before starting; no `--phase` means the new code never executes).

```bash
git add src/bin/copy_assign.rs
git commit -m "feat(copy_assign): --phase emits the copy_graph (REFERENCE walk + tagged copy paths + legend)"
```

---

### Task 8: Optional `--gff` — annotation axis

**Files:**
- Modify: `src/bin/copy_assign.rs` (add `--gff` arg; pass annotation status into `build_copy_graph`)
- Test: `src/bin/copy_assign.rs` tests

**Interfaces:**
- Consumes: `build_copy_graph` (Task 7), `fa.copy_spans` (genomic span per copy).
- Produces: in-genome copies tagged `InGenomeAnnotated` / `InGenomeUnannotated` by overlap with GFF/BED features; `AnnotationUnknown` when no `--gff`.

- [ ] **Step 1: Write the failing test**

```rust
    #[test]
    fn annotation_axis_from_intervals() {
        use rustle::vg_family::denovo_pipeline::FamilyAssignment;
        use rustle::vg_family::copy_graph::CopyStatus;
        let mut fa = FamilyAssignment::empty();
        fa.chrom = "chr1".into();
        fa.copy_tids = vec!["c0".into(), "c1".into()];
        fa.copy_spans = vec![("chr1".into(), 1000, 2000), ("chr1".into(), 5000, 6000)];
        fa.psv_col_pos = vec![Some(1500)];
        fa.copy_psv_alleles = vec![vec![Some(b'A')], vec![Some(b'A')]];
        // annotation covers only copy0's span
        let ann = vec![("chr1".to_string(), 900u64, 2100u64)];
        let s0 = annotation_status(&fa, 0, Some(&ann));
        let s1 = annotation_status(&fa, 1, Some(&ann));
        assert_eq!(s0, CopyStatus::InGenomeAnnotated);
        assert_eq!(s1, CopyStatus::InGenomeUnannotated);
        assert_eq!(annotation_status(&fa, 1, None), CopyStatus::AnnotationUnknown);
    }
```

- [ ] **Step 2: Run test to verify it fails**

Run: `cargo test -p rustle --bin copy_assign annotation_axis_from_intervals 2>&1 | tail -20`
Expected: FAIL — `annotation_status` not defined.

- [ ] **Step 3: Write minimal implementation**

Add a free function and wire it into the in-genome branch of `build_copy_graph` (replace the `CopyStatus::AnnotationUnknown` for non-absent copies with `annotation_status(fa, ci, ann)`, threading an `ann: Option<&[(String,u64,u64)]>` param through `build_copy_graph`):

```rust
/// In-genome annotation axis: overlap of copy `ci`'s span with any annotated interval.
/// `None` intervals => AnnotationUnknown (we never claim "unannotated" unchecked).
fn annotation_status(
    fa: &rustle::vg_family::denovo_pipeline::FamilyAssignment,
    ci: usize,
    ann: Option<&[(String, u64, u64)]>,
) -> rustle::vg_family::copy_graph::CopyStatus {
    use rustle::vg_family::copy_graph::CopyStatus;
    let Some(ann) = ann else { return CopyStatus::AnnotationUnknown };
    let Some((c, s, e)) = fa.copy_spans.get(ci) else { return CopyStatus::AnnotationUnknown };
    let hit = ann.iter().any(|(ac, as_, ae)| ac == c && *as_ < *e && *s < *ae);
    if hit { CopyStatus::InGenomeAnnotated } else { CopyStatus::InGenomeUnannotated }
}
```

Add the CLI arg to `struct Args` (near the other options):

```rust
    /// Optional gene annotation (GFF3/GTF/BED) to tag in-genome copies annotated vs unannotated in the
    /// --phase copy graph. Without it, in-genome copies are tagged `annotation-unknown`.
    #[arg(long)]
    gff: Option<String>,
```

Parse it once before the sweep into `Vec<(String,u64,u64)>` (chrom, start0, end) — accept BED (0-based, 3 cols) and GFF/GTF (1-based, cols 1/4/5; convert start to 0-based). Pass `args.gff`-derived intervals into the `build_copy_graph` call. Minimal parser:

```rust
fn parse_annotation(path: &str) -> anyhow::Result<Vec<(String, u64, u64)>> {
    let mut out = Vec::new();
    for line in std::fs::read_to_string(path)?.lines() {
        if line.starts_with('#') || line.trim().is_empty() { continue; }
        let f: Vec<&str> = line.split('\t').collect();
        if f.len() >= 5 && f[3].parse::<u64>().is_ok() && f[4].parse::<u64>().is_ok() {
            out.push((f[0].to_string(), f[3].parse::<u64>()?.saturating_sub(1), f[4].parse()?)); // GFF/GTF 1-based
        } else if f.len() >= 3 && f[1].parse::<u64>().is_ok() && f[2].parse::<u64>().is_ok() {
            out.push((f[0].to_string(), f[1].parse()?, f[2].parse()?)); // BED 0-based
        }
    }
    Ok(out)
}
```

- [ ] **Step 4: Run test to verify it passes**

Run: `cargo test -p rustle --bin copy_assign annotation_axis_from_intervals 2>&1 | tail -20`
Expected: PASS. Build: `cargo build -p rustle --bin copy_assign 2>&1 | tail -5`.

- [ ] **Step 5: Commit**

```bash
git add src/bin/copy_assign.rs
git commit -m "feat(copy_assign): optional --gff tags in-genome copies annotated/unannotated in the copy graph"
```

---

### Task 9: DSFAM26 render smoke (data-gated)

**Files:**
- Test: manual command (data-dependent; not a CI unit test)

**Interfaces:**
- Consumes: the full feature.

- [ ] **Step 1: Render DSFAM26 and eyeball in Bandage**

Run (FOREGROUND, output under winloci_scratch — WSL2 crash rule):
```bash
./target/release/copy_assign \
  --bam /home/juanfra/winloci_scratch/GGO_mm.bam \
  --fasta /home/juanfra/winloci_scratch/GGO.fasta \
  --region NC_073229.2:49040000-49075000 \
  --out /home/juanfra/winloci_scratch/dsfam26_o4 \
  --phase --dump-psv --collapse
grep -c "^S" /home/juanfra/winloci_scratch/dsfam26_o4.phase.gfa
grep -c "REFERENCE" /home/juanfra/winloci_scratch/dsfam26_o4.phase.gfa
grep -c "_ABSENT" /home/juanfra/winloci_scratch/dsfam26_o4.phase.gfa
```
Expected: a non-zero S count, exactly one `_REFERENCE` P-line, and ≥1 `_ABSENT` copy P-line (the collapsed MHC copy). If `_ABSENT` is 0, the collapsed copy was not surfaced — confirm the collapse/absent-copy flag that populates the second copy is enabled for the region (check `--collapse` vs the current flag name in `copy_assign --help`; the emitter renders whatever copies the assignment produced).

- [ ] **Step 2: Load in Bandage**

Load `dsfam26_o4.phase.gfa`, apply `dsfam26_o4.phase.gfa.colours.csv`. Confirm the `_ABSENT` copy path threads red nodes the grey REFERENCE walk skips. No commit (verification only).

---

## Self-Review

**Spec coverage:**
- REFERENCE walk → Task 3. ✓
- Backbone spacers (not intronic seq) → Task 2. ✓
- PSV allele bubbles + L-lines (no dangling) → Tasks 2 & 5 (dangling assertion). ✓
- Copy P-lines + RC/SU/MI/ST + `_ABSENT` naming + honesty (omit unknown) → Task 4. ✓
- SU computed as private-column count → Task 4. ✓
- colours.csv segment-keyed (absent red / reference grey) + legend → Task 6. ✓
- Wiring, default byte-identical, keep existing TSVs → Task 7. ✓
- Optional `--gff` annotation axis (unknown when absent) → Task 8. ✓
- DSFAM26 smoke → Task 9. ✓
- Status taxonomy enum → Task 1. ✓
- v2 (exon presence/absence, ORF) excluded → not in any task. ✓

**Placeholder scan:** no TBD/TODO; every code step shows complete code; test bodies are concrete.

**Type consistency:** `CopyGraph`/`CopyPath`/`PsvColumn`/`ReadWalk`/`Corrob`/`CopyStatus`/`GfaLines` used identically across Tasks 1–8; `gfa_lines`/`to_gfa`/`colours_csv`/`legend_tsv`/`walk_tokens`/`private_columns`/`alleles_at`/`build_copy_graph`/`annotation_status` names are stable. One integration caveat flagged in Task 7 (`gfa_links` type change `HashSet<(String,String)>` → `HashSet<String>`, `Assignment::status_is_assigned` may need inlining) — resolve against the real code during implementation.
