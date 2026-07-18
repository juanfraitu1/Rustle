# Copy-Graph v2 Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Add the `MI` genome-map-identity tag and a sibling `<out>.exon.gfa` exon presence/absence graph to `copy_assign --phase`, so a copy differing by exon content (esp. a reference-absent/unannotated copy) is a visible arm nothing else walks.

**Architecture:** One shared upstream fix threads genomic per-copy data (`copy_introns`, `copy_map_identity`) into `FamilyAssignment`; then `MI` reads the identity, and a new `ExonGraph` type (in `copy_graph.rs`) builds genomic exon classes from `copy_introns` and emits a separate GFA. The v1 `CopyGraph`/`.phase.gfa` is untouched.

**Tech Stack:** Rust (crate `rustle`), `rustle::genome::GenomeIndex::fetch_sequence`, GFA 1.1, Bandage.

## Global Constraints

- Default `copy_assign` (no `--phase`) MUST be byte-identical; the new `FamilyAssignment` fields are consumed only by v2 code.
- New `FamilyAssignment` fields (`copy_map_identity`, `copy_introns`) MUST be length-aligned with `copy_tids` at EVERY construction/push site.
- Valid GFA 1.1; every P-line walk step backed by an L-line (skip + detour edges both emitted via union-of-adjacencies) — no dangling.
- Honesty rule: `MI` omitted when unknown, never faked (in-genome copies stay omitted; no fake `1.0`).
- No new k-mer computation. v1 `CopyGraph` and `.phase.gfa`/`.colours.csv`/`.legend.tsv` output unchanged.
- OUT OF SCOPE: individual read W-lines in the exon graph, in-genome self-identity, reviving the dead `FamilyGraph`/`consensus.rs` layer.

## File Structure

- **Modify** `src/rustle/vg_family/absent_copy.rs` — `Admission::Copy` carries the remap identity (Task 1).
- **Modify** `src/rustle/vg_family/denovo_pipeline.rs` — `FamilyAssignment` gains `copy_map_identity` + `copy_introns`, populated at all sites (Task 2).
- **Modify** `src/bin/copy_assign.rs` — `build_copy_graph` fills `Corrob.map_identity` (Task 3); new `build_exon_graph` + `.exon.gfa` writer (Task 8).
- **Modify** `src/rustle/vg_family/copy_graph.rs` — new `ExonClass`/`CopyExonPath`/`ExonGraph` types + clustering + emitters (Tasks 4–7).

---

### Task 1: `Admission::Copy` carries the remap identity

**Files:**
- Modify: `src/rustle/vg_family/absent_copy.rs` (enum `:72`, gate-5 `:152-155`)
- Modify: `src/rustle/vg_family/denovo_pipeline.rs` (match arms `:1233`, `:1514`; test literal `:3417`)
- Test: `absent_copy.rs` tests

**Interfaces:**
- Produces: `Admission::Copy(DenovoTranscript, Option<f64>)` — the second field is the genome remap identity (`Some(id)`), `None` only where an admitted copy has no computed identity.

- [ ] **Step 1: Write the failing test**

Add to `absent_copy.rs` tests (adapt to the existing hermetic admit-test harness if one exists; otherwise this enum-shape test):
```rust
#[test]
fn admission_copy_carries_identity() {
    let t = DenovoTranscript { tid: "t".into(), chrom: "c".into(), start: 0, end: 10, n_reads: 5,
        strand: '+', introns: vec![], seq: b"ACGTACGTAC".to_vec() };
    let a = Admission::Copy(t.clone(), Some(0.95));
    match a { Admission::Copy(_, id) => assert_eq!(id, Some(0.95)), _ => panic!("wrong variant") }
}
```

- [ ] **Step 2: Run test to verify it fails**

Run: `cargo test -p rustle --lib admission_copy_carries_identity 2>&1 | tail -20`
Expected: FAIL — `Admission::Copy` takes one field.

- [ ] **Step 3: Write minimal implementation**

In `absent_copy.rs:72`, change the variant:
```rust
    Copy(DenovoTranscript, Option<f64>),
```
At gate 5 (`:152-155`), pass the bound `id`:
```rust
    match remap_identity(&t.seq) {
        Some(id) if id < p.remap_max_identity => Admission::Copy(t, Some(id)),
        Some(_) => dna_needs(cand, ">=98% remap identity (paralog-leak or het)"),
        None => dna_needs(cand, "no homology on remap"),
    }
```
Search `absent_copy.rs` for EVERY other `Admission::Copy(` construction and add the second field (`Some(id)` if an identity is in scope, else `None`). Update the match arms:
- `denovo_pipeline.rs:1233`: `Admission::Copy(t) => t,` → `Admission::Copy(t, _) => t,` (Task 2 will use the id here).
- `denovo_pipeline.rs:1514`: `Admission::Copy(t) => admitted.push(t),` → `Admission::Copy(t, _) => admitted.push(t),` (Task 2 will capture the id).
- `denovo_pipeline.rs:3417` (test literal): add the second field.

- [ ] **Step 4: Run test to verify it passes**

Run: `cargo test -p rustle --lib admission_copy_carries_identity 2>&1 | tail -20` (PASS) and `cargo build -p rustle 2>&1 | tail -5` (Finished — all arms updated).

- [ ] **Step 5: Commit**
```bash
git add src/rustle/vg_family/absent_copy.rs src/rustle/vg_family/denovo_pipeline.rs
git commit -m "feat(absent_copy): Admission::Copy carries the genome remap identity"
```

---

### Task 2: `FamilyAssignment` gains `copy_map_identity` + `copy_introns`

**Files:**
- Modify: `src/rustle/vg_family/denovo_pipeline.rs` (struct `:377`, `empty()` `:449`, main build `:1623-1640`, stage-2 admission `:1508-1537`, admit_novel_pools `:1240-1244`)
- Test: `denovo_pipeline.rs` tests

**Interfaces:**
- Produces: `FamilyAssignment.copy_map_identity: Vec<Option<f64>>` and `FamilyAssignment.copy_introns: Vec<Vec<(u64,u64)>>`, both length-aligned with `copy_tids`. Absent copies → `Some(remap_id)`; in-genome copies → `None`. `copy_introns[k]` = the copy's genomic `(donor,acceptor)` chain.

- [ ] **Step 1: Write the failing test**

Add to `denovo_pipeline.rs` tests:
```rust
#[test]
fn family_assignment_has_aligned_intron_and_identity_fields() {
    let fa = FamilyAssignment::empty();
    assert_eq!(fa.copy_introns.len(), 0);
    assert_eq!(fa.copy_map_identity.len(), 0);
    // a built fa (via the existing detect_and_assign test harness) must keep these aligned with copy_tids:
    // reuse whatever test already builds a FamilyAssignment and assert:
    //   assert_eq!(fa.copy_introns.len(), fa.copy_tids.len());
    //   assert_eq!(fa.copy_map_identity.len(), fa.copy_tids.len());
}
```
Also strengthen an EXISTING `detect_and_assign`/admission test (e.g. near `:2440-2540`, the freeze/admit tests) to assert `fa.copy_introns.len() == fa.copy_tids.len()` and `fa.copy_map_identity.len() == fa.copy_tids.len()`, and that an admitted absent copy's `copy_map_identity` entry `.is_some()`.

- [ ] **Step 2: Run test to verify it fails**

Run: `cargo test -p rustle --lib family_assignment_has_aligned 2>&1 | tail -20`
Expected: FAIL — fields don't exist.

- [ ] **Step 3: Write minimal implementation**

Add the fields to the struct (`denovo_pipeline.rs:377`, after `copy_junctions`):
```rust
    /// genomic (donor,acceptor) intron chain per copy, parallel to copy_tids. For v2 exon graph.
    pub copy_introns: Vec<Vec<(u64,u64)>>,
    /// per-copy genome remap identity, parallel to copy_tids: Some for reference-ABSENT copies
    /// (their discovery remap identity), None for in-genome copies (never remapped). For MI:f:.
    pub copy_map_identity: Vec<Option<f64>>,
```
In `empty()` (`:449`): add `copy_introns: Vec::new(),` and `copy_map_identity: Vec::new(),`.

Thread a tid→identity map so `copy_map_identity` survives admission/pruning (keyed on tid, aligned with `copy_tids` which is built from `all_copies`). Near the top of the family loop where `all_copies` is first bound, add:
```rust
    let mut remap_id_by_tid: std::collections::HashMap<String, f64> = std::collections::HashMap::new();
```
In the stage-2 admission block (`:1508-1537`), capture the id from the updated arm:
```rust
    match absent_copy::admit_candidate(cand, host, genome, fasta_path, &AbsentCopyParams::default()) {
        Admission::Copy(t, id) => { if let Some(v) = id { remap_id_by_tid.insert(t.tid.clone(), v); } admitted.push(t); }
        Admission::DnaNeeds(r) => dna_needs.push(r),
    }
```
In `admit_novel_pools_with_admitter` (`:1233` + the push block `:1240-1244`), capture the id and push ALL parallel fields (this site currently pushes `copy_tids` without `copy_spans` — fix that gap too so every per-copy vector stays aligned):
```rust
    let (t, adm_id) = match admitted { Admission::Copy(t, id) => (t, id), Admission::DnaNeeds(_) => continue };
    // ... after fa.copy_tids.push(t.tid.clone()); add:
    fa.copy_spans.push((t.chrom.clone(), t.start, t.end));
    fa.copy_introns.push(t.introns.clone());
    fa.copy_map_identity.push(adm_id);
```
In the main struct build (`:1640`, alongside `copy_tids`/`copy_spans` which map over `all_copies`):
```rust
    copy_introns: all_copies.iter().map(|c| c.introns.clone()).collect(),
    copy_map_identity: all_copies.iter().map(|c| remap_id_by_tid.get(&c.tid).copied()).collect(),
```
(Any other place that constructs a `FamilyAssignment` literal — e.g. `gated_family` — must add the two fields; `cargo build` will flag them. Use `Vec::new()` for empty families, or `vec![None; n]`/`vec![Vec::new(); n]` where a copy roster exists.)

- [ ] **Step 4: Run test to verify it passes**

Run: `cargo test -p rustle --lib family_assignment_has_aligned 2>&1 | tail -20` (PASS) + `cargo build -p rustle 2>&1 | tail -5` (Finished) + `cargo test -p rustle --lib denovo_pipeline:: 2>&1 | grep "test result"` (no regressions).

- [ ] **Step 5: Commit**
```bash
git add src/rustle/vg_family/denovo_pipeline.rs
git commit -m "feat(denovo): FamilyAssignment carries per-copy genomic introns + map identity"
```

---

### Task 3: MI tag — `build_copy_graph` fills `Corrob.map_identity`

**Files:**
- Modify: `src/bin/copy_assign.rs` (`build_copy_graph`, `Corrob` construction `:420`)
- Test: `copy_assign.rs` tests

**Interfaces:**
- Consumes: `fa.copy_map_identity` (Task 2).

- [ ] **Step 1: Write the failing test**
```rust
#[test]
fn build_copy_graph_fills_mi_from_copy_map_identity() {
    use rustle::vg_family::denovo_pipeline::FamilyAssignment;
    let mut fa = FamilyAssignment::empty();
    fa.chrom = "chr1".into();
    fa.copy_tids = vec!["c0".into(), "c1".into()];
    fa.psv_col_pos = vec![Some(100)];
    fa.copy_psv_alleles = vec![vec![Some(b'A')], vec![Some(b'G')]];
    fa.copy_map_identity = vec![None, Some(0.952)];
    let g = build_copy_graph("CAFAM0", &fa, |_c, _p| Some(b'A'), &[]);
    let gfa = g.to_gfa();
    let c1 = gfa.lines().find(|l| l.starts_with("P\tCAFAM0_copy1")).unwrap();
    assert!(c1.contains("MI:f:0.952"), "copy1 MI missing: {}", c1);
    let c0 = gfa.lines().find(|l| l.starts_with("P\tCAFAM0_copy0")).unwrap();
    assert!(!c0.contains("MI:f:"), "copy0 must omit MI: {}", c0);
}
```
(If `build_copy_graph`'s signature already has the `ann` param from v1's Task 8, pass `None` for it.)

- [ ] **Step 2: Run test to verify it fails**

Run: `cargo test -p rustle --bin copy_assign build_copy_graph_fills_mi 2>&1 | tail -20`
Expected: FAIL — `MI:f:` absent (still hardcoded `None`).

- [ ] **Step 3: Write minimal implementation**

In `build_copy_graph`, at the `Corrob` construction (`:420`), replace `map_identity: None` with:
```rust
    map_identity: fa.copy_map_identity.get(copy_idx).copied().flatten(),
```
(Use the copy-loop index variable — `copy_idx` per v1's Task 4 rename.)

- [ ] **Step 4: Run test to verify it passes**

Run: `cargo test -p rustle --bin copy_assign build_copy_graph_fills_mi 2>&1 | tail -20` (PASS) + `cargo test -p rustle --bin copy_assign 2>&1 | grep "test result"`.

- [ ] **Step 5: Commit**
```bash
git add src/bin/copy_assign.rs
git commit -m "feat(copy_assign): MI:f: tag from FamilyAssignment.copy_map_identity"
```

---

### Task 4: Exon-graph types + skeleton

**Files:**
- Modify: `src/rustle/vg_family/copy_graph.rs` (new types; reuse v1 `CopyStatus`/`Corrob`)
- Test: same file's `tests`

**Interfaces:**
- Produces: `ExonClass { chrom: String, start: u64, end: u64 }`, `CopyExonPath { id: String, exon_nodes: Vec<usize>, status: CopyStatus, corrob: Corrob }`, `ExonGraph { family: String, nodes: Vec<ExonClass>, copies: Vec<CopyExonPath> }`.

- [ ] **Step 1: Write the failing test**
```rust
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
```

- [ ] **Step 2: Run test to verify it fails**

Run: `cargo test -p rustle --lib exon_graph_constructs 2>&1 | tail -20`
Expected: FAIL — types absent.

- [ ] **Step 3: Write minimal implementation**

Add to `copy_graph.rs`:
```rust
/// A shared exon node in the exon presence/absence graph — one genomic exon interval.
#[derive(Clone, Debug)]
pub struct ExonClass { pub chrom: String, pub start: u64, pub end: u64 }

/// A copy as an ordered walk over exon-class indices.
#[derive(Clone, Debug)]
pub struct CopyExonPath { pub id: String, pub exon_nodes: Vec<usize>, pub status: CopyStatus, pub corrob: Corrob }

/// Family exon presence/absence graph. `nodes` sorted by genomic start; each copy walks a subset.
#[derive(Clone, Debug)]
pub struct ExonGraph { pub family: String, pub nodes: Vec<ExonClass>, pub copies: Vec<CopyExonPath> }
```

- [ ] **Step 4: Run test to verify it passes**

Run: `cargo test -p rustle --lib exon_graph_constructs 2>&1 | tail -20` (PASS).

- [ ] **Step 5: Commit**
```bash
git add src/rustle/vg_family/copy_graph.rs
git commit -m "feat(copy_graph): exon-graph object types"
```

---

### Task 5: Exon-class clustering (`ExonGraph::from_copies`)

**Files:**
- Modify: `src/rustle/vg_family/copy_graph.rs`
- Test: same file's `tests`

**Interfaces:**
- Consumes: per-copy genomic exons.
- Produces: `ExonGraph::from_copies(family: &str, copies: &[(String, CopyStatus, Corrob, String, Vec<(u64,u64)>)]) -> ExonGraph` — each tuple is `(copy_id, status, corrob, chrom, exons)`. Clusters exons by reciprocal-overlap (≥0.30, same chrom) into position-sorted classes; each `CopyExonPath.exon_nodes` = the class indices that copy contributes an exon to.

- [ ] **Step 1: Write the failing test**
```rust
    #[test]
    fn from_copies_clusters_and_flags_copy_specific_exon() {
        // copy0 exons E1,E3 ; copy1 exons E1,E2(extra),E3 — E2 is copy1-specific.
        let copies = vec![
            ("F_copy0".to_string(), CopyStatus::InGenomeAnnotated, Corrob::default(),
                vec![(0u64,100u64), (300,400)]),
            ("F_copy1".to_string(), CopyStatus::AbsentDivergent, Corrob::default(),
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
```

- [ ] **Step 2: Run test to verify it fails**

Run: `cargo test -p rustle --lib from_copies_clusters 2>&1 | tail -20`
Expected: FAIL — `from_copies` not found.

- [ ] **Step 3: Write minimal implementation**
```rust
impl ExonGraph {
    /// reciprocal overlap = min(inter/len_a, inter/len_b); 0 if disjoint or different chrom.
    fn recip_overlap(a: (&str,u64,u64), b: (&str,u64,u64)) -> f64 {
        if a.0 != b.0 { return 0.0; }
        let lo = a.1.max(b.1); let hi = a.2.min(b.2);
        if hi <= lo { return 0.0; }
        let inter = (hi - lo) as f64;
        (inter / (a.2 - a.1) as f64).min(inter / (b.2 - b.1) as f64)
    }

    pub fn from_copies(family: &str, copies: &[(String, CopyStatus, Corrob, String, Vec<(u64,u64)>)]) -> ExonGraph {
        // flatten every exon as (copy_idx, chrom, start, end)
        let mut flat: Vec<(usize, String, u64, u64)> = Vec::new();
        for (ci, (_, _, _, chrom, exons)) in copies.iter().enumerate() {
            for &(s, e) in exons { flat.push((ci, chrom.clone(), s, e)); }
        }
        // union-find over flat exon items
        let n = flat.len();
        let mut parent: Vec<usize> = (0..n).collect();
        fn find(p: &mut Vec<usize>, x: usize) -> usize { if p[x]!=x { let r=find(p,p[x]); p[x]=r; } p[x] }
        for i in 0..n { for j in (i+1)..n {
            if Self::recip_overlap((&flat[i].1, flat[i].2, flat[i].3), (&flat[j].1, flat[j].2, flat[j].3)) >= 0.30 {
                let (a,b)=(find(&mut parent,i),find(&mut parent,j)); if a!=b { parent[a]=b; }
            }
        }}
        // group items by root, build one ExonClass per group (min start, max end)
        use std::collections::BTreeMap;
        let mut groups: BTreeMap<usize, Vec<usize>> = BTreeMap::new();
        for i in 0..n { let r = find(&mut parent, i); groups.entry(r).or_default().push(i); }
        let mut classes: Vec<(String, u64, u64, Vec<usize>)> = groups.values().map(|idxs| {
            let start = idxs.iter().map(|&i| flat[i].2).min().unwrap();
            let end = idxs.iter().map(|&i| flat[i].3).max().unwrap();
            let members: Vec<usize> = idxs.iter().map(|&i| flat[i].0).collect();
            (flat[idxs[0]].1.clone(), start, end, members)
        }).collect();
        classes.sort_by_key(|c| c.1);   // by genomic start -> E0..En
        let nodes: Vec<ExonClass> = classes.iter().map(|(chrom,s,e,_)| ExonClass { chrom: chrom.clone(), start:*s, end:*e }).collect();
        let copy_paths: Vec<CopyExonPath> = copies.iter().enumerate().map(|(ci, (id, status, corrob, _chrom, _exons))| {
            let mut exon_nodes: Vec<usize> = (0..classes.len()).filter(|&k| classes[k].3.contains(&ci)).collect();
            exon_nodes.sort();
            CopyExonPath { id: id.clone(), exon_nodes, status: *status, corrob: corrob.clone() }
        }).collect();
        ExonGraph { family: family.to_string(), nodes, copies: copy_paths }
    }
}
```

- [ ] **Step 4: Run test to verify it passes**

Run: `cargo test -p rustle --lib from_copies_clusters 2>&1 | tail -20` (PASS).

- [ ] **Step 5: Commit**
```bash
git add src/rustle/vg_family/copy_graph.rs
git commit -m "feat(copy_graph): cluster per-copy exons into shared exon classes"
```

---

### Task 6: `ExonGraph::to_gfa` — REFERENCE walk + copy walks + skip edges + RC

**Files:**
- Modify: `src/rustle/vg_family/copy_graph.rs`
- Test: same file's `tests`

**Interfaces:**
- Produces: `ExonGraph::to_gfa(&self, exon_seq: impl Fn(&ExonClass) -> Vec<u8>) -> String` — GFA 1.1: `S` per class (`{fam}_E{k}` seq `PO:i:start` `RC:i:{sum of walking copies' corrob.reads}`), a `{fam}_REFERENCE` P-line over classes present in ≥1 non-absent copy (fallback: all copies) `ST:Z:reference`, one `{fam}_copy{ci}[_ABSENT]` P-line per copy with `RC`/`MI`/`ST` tags, and L-lines = union of consecutive-class adjacencies over the REFERENCE + all copy walks.

- [ ] **Step 1: Write the failing test**
```rust
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
        assert_no_dangling(&gfa);   // reuse the v1 helper (parses P/W walks vs L set)
    }
```

- [ ] **Step 2: Run test to verify it fails**

Run: `cargo test -p rustle --lib exon_gfa_has_reference 2>&1 | tail -20`
Expected: FAIL — `to_gfa` not found (and `from_copies` signature must be the chrom-tagged 5-tuple from Task 5's note).

- [ ] **Step 3: Write minimal implementation**
```rust
impl ExonGraph {
    fn node(&self, k: usize) -> String { format!("{}_E{}", self.family, k) }

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
}
```
(`from_copies` from Task 5 must be the chrom-tagged 5-tuple version so `ExonClass.chrom` is real.)

- [ ] **Step 4: Run test to verify it passes**

Run: `cargo test -p rustle --lib exon_gfa_has_reference 2>&1 | tail -20` (PASS) + `cargo test -p rustle --lib copy_graph:: 2>&1 | grep "test result"` (no regressions).

- [ ] **Step 5: Commit**
```bash
git add src/rustle/vg_family/copy_graph.rs
git commit -m "feat(copy_graph): ExonGraph::to_gfa — REFERENCE walk, copy walks, skip edges, per-class RC"
```

---

### Task 7: `ExonGraph` colours + legend

**Files:**
- Modify: `src/rustle/vg_family/copy_graph.rs`
- Test: same file's `tests`

**Interfaces:**
- Produces: `ExonGraph::colours_csv(&self) -> String` (rows `segment,#RRGGBB`) and `ExonGraph::legend_tsv(&self) -> String`. Rule: a class walked by exactly one copy that is absent → that copy's status colour (red); a class on the reference (≥1 non-absent copy) → grey `#9aa0a6`; else the walking copy's status colour.

- [ ] **Step 1: Write the failing test**
```rust
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
```

- [ ] **Step 2: Run test to verify it fails**

Run: `cargo test -p rustle --lib exon_colours_arm_red 2>&1 | tail -20`
Expected: FAIL — `colours_csv` not found.

- [ ] **Step 3: Write minimal implementation**
```rust
impl ExonGraph {
    pub fn colours_csv(&self) -> String {
        use std::collections::BTreeMap;
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
    pub fn legend_tsv(&self) -> String {
        use std::collections::BTreeSet;
        let mut seen: BTreeSet<&'static str> = BTreeSet::new();
        let mut s = format!("reference\t{}\n", CopyStatus::Reference.colour()); seen.insert("reference");
        for c in &self.copies { if seen.insert(c.status.tag()) { s.push_str(&format!("{}\t{}\n", c.status.tag(), c.status.colour())); } }
        s
    }
}
```

- [ ] **Step 4: Run test to verify it passes**

Run: `cargo test -p rustle --lib exon_colours_arm_red 2>&1 | tail -20` (PASS).

- [ ] **Step 5: Commit**
```bash
git add src/rustle/vg_family/copy_graph.rs
git commit -m "feat(copy_graph): ExonGraph colours (arm by owner, shared grey) + legend"
```

---

### Task 8: Wire `build_exon_graph` into `copy_assign --phase`

**Files:**
- Modify: `src/bin/copy_assign.rs` (inside the `if args.phase` block; new `build_exon_graph` free fn; `.exon.gfa` writer)
- Test: `copy_assign.rs` tests

**Interfaces:**
- Consumes: `fa.copy_introns`, `fa.copy_spans`, `fa.copy_map_identity`, `fa.assignments`, `fa.copy_tids`; the per-copy status derivation from v1's `build_copy_graph`; `genome_for(chrom).fetch_sequence`.
- Produces: `<out>.exon.gfa` + `.exon.gfa.colours.csv` + `.exon.gfa.legend.tsv`.

- [ ] **Step 1: Write the failing test**
```rust
#[test]
fn build_exon_graph_makes_copy_specific_arm() {
    use rustle::vg_family::denovo_pipeline::FamilyAssignment;
    let mut fa = FamilyAssignment::empty();
    fa.chrom = "chr1".into();
    fa.copy_tids = vec!["c0".into(), "c1".into()];
    fa.copy_spans = vec![("chr1".into(), 0, 400), ("chr1".into(), 0, 400)];
    // copy0 introns skip 100-300 (exons 0-100, 300-400); copy1 has an extra exon (exons 0-100,150-250,300-400)
    fa.copy_introns = vec![ vec![(100,300)], vec![(100,150),(250,300)] ];
    fa.copy_map_identity = vec![None, Some(0.95)];
    fa.assignments = vec![]; // no reads needed for structure
    let g = build_exon_graph("CAFAM0", &fa, |_c,_s,_e| vec![b'A'; 10]);
    // copy1 walks one more class than copy0
    assert!(g.copies[1].exon_nodes.len() > g.copies[0].exon_nodes.len());
    let gfa = g.to_gfa(|ec| vec![b'A'; (ec.end-ec.start) as usize]);
    assert!(gfa.contains("P\tCAFAM0_REFERENCE"));
    assert!(gfa.lines().any(|l| l.starts_with("P\tCAFAM0_copy1_ABSENT")));
}
```
(Absent status for copy1 comes from `fa.copy_map_identity[1].is_some()` OR the same discovery-coupled logic v1 uses — see Step 3; here `Some(0.95)` marks copy1 absent.)

- [ ] **Step 2: Run test to verify it fails**

Run: `cargo test -p rustle --bin copy_assign build_exon_graph_makes 2>&1 | tail -20`
Expected: FAIL — `build_exon_graph` not defined.

- [ ] **Step 3: Write minimal implementation**

Add `build_exon_graph` (free fn) to `copy_assign.rs`. Reconstruct each copy's genomic exons from `copy_introns[ci]` + `copy_spans[ci]`, derive status the SAME way `build_copy_graph` does (reuse a shared helper or replicate: absent iff a discovery_coupled read assigned to it OR `fa.copy_map_identity[ci].is_some()`; subtype by span co-location; else in-genome via `annotation_status`), build `Corrob { reads: <assigned count>, suns: None, map_identity: fa.copy_map_identity[ci] }`, and call `ExonGraph::from_copies`:
```rust
fn build_exon_graph(
    fid: &str,
    fa: &rustle::vg_family::denovo_pipeline::FamilyAssignment,
    _exon_seq_probe: impl Fn(&str, u64, u64) -> Vec<u8>,   // reserved; to_gfa takes the real fetch
) -> rustle::vg_family::copy_graph::ExonGraph {
    use rustle::vg_family::copy_graph::*;
    let n = fa.copy_tids.len();
    // reuse build_copy_graph's status derivation (extract it into a shared `copy_status(fa, ci)` helper).
    let copies: Vec<(String, CopyStatus, Corrob, String, Vec<(u64,u64)>)> = (0..n).map(|ci| {
        let (chrom, start, end) = fa.copy_spans.get(ci).cloned().unwrap_or((fa.chrom.clone(), 0, 0));
        let introns = fa.copy_introns.get(ci).cloned().unwrap_or_default();
        // genomic exons from intron chain + outer bounds:
        let mut exons = Vec::with_capacity(introns.len()+1); let mut prev = start;
        for (d,a) in &introns { exons.push((prev,*d)); prev = *a; } exons.push((prev,end));
        let reads = fa.assignments.iter().filter(|(_,a)| a.best_copy==ci && matches!(a.status, AssignStatus::Assigned)).count() as u32;
        let status = copy_status(fa, ci);   // shared helper extracted from build_copy_graph
        let corrob = Corrob { reads: Some(reads), suns: None, map_identity: fa.copy_map_identity.get(ci).copied().flatten() };
        (format!("{}_copy{}", fid, ci), status, corrob, chrom, exons)
    }).collect();
    ExonGraph::from_copies(fid, &copies)
}
```
(Extract the v1 status derivation into `fn copy_status(fa, ci) -> CopyStatus` and call it from both `build_copy_graph` and here — DRY.) Then, inside the `if args.phase` block after the v1 `.phase.gfa` is written, build+write the exon graph per family and fold into `<out>.exon.gfa` files, using `genome_for(&ec.chrom).fetch_sequence(&ec.chrom, ec.start, ec.end)` as the `to_gfa` `exon_seq` closure (fall back to `vec![b'N'; (end-start) as usize]` on `None`, logging a count). Write `<out>.exon.gfa` (header once + folded S/L/P), `<out>.exon.gfa.colours.csv`, `<out>.exon.gfa.legend.tsv`.

- [ ] **Step 4: Run test to verify it passes**

Run: `cargo test -p rustle --bin copy_assign build_exon_graph_makes 2>&1 | tail -20` (PASS) + `cargo build --release -p rustle --bin copy_assign 2>&1 | tail -5` (Finished — binary rebuilt for Task 9).

- [ ] **Step 5: Byte-identical guard + commit**

Confirm no-`--phase` output unchanged (structural: all edits inside `if args.phase` + new free fns). Then:
```bash
git add src/bin/copy_assign.rs
git commit -m "feat(copy_assign): --phase emits <out>.exon.gfa (exon presence/absence arms)"
```

---

### Task 9: GWFAM61 render smoke + final review

**Files:** manual (data-gated)

- [ ] **Step 1: Render GWFAM61 and inspect the exon arm**

Run FOREGROUND, output under winloci_scratch (WSL2 crash rule). GWFAM61 = human CNTNAP3/CNTNAP3B, chr9; the family forms with `--homology-primary`. Data on `/mnt/linuxdisk` (mount first if needed):
```bash
RUSTLE_SKIP_POA_DIAGNOSTIC=1 ./target/release/copy_assign \
  --bam /mnt/linuxdisk/home/juanfraitu/winloci_data/soto_reads.bam \
  --fasta /mnt/linuxdisk/home/juanfraitu/winloci_data/chm13v2.0.fa \
  --region chr9:39086152-43680114 --out /home/juanfra/winloci_scratch/gwfam61 \
  --phase --dump-psv --homology-primary
grep -c "^S" /home/juanfra/winloci_scratch/gwfam61.exon.gfa
grep -c "REFERENCE" /home/juanfra/winloci_scratch/gwfam61.exon.gfa
grep "^P" /home/juanfra/winloci_scratch/gwfam61.exon.gfa | awk -F'\t' '{print $2}'
```
Expected: a non-empty exon graph with a `_REFERENCE` walk and two copy walks differing by ≥1 exon node. If 0 families, the region/flags need adjusting (the family region may need widening or the two loci run separately then linked — note this in the report). ⚠The chr9 loci are ~4.5 Mb apart; if a single region is too heavy, run each locus (`chr9:39086152-39301963` and `chr9:43536410-43680114`) — but then they are separate families, so for the demo prefer whichever surfaces a single family with both copies, else document that the exon-arm is validated by the unit tests and a co-located family (e.g. a compact PCDHB/GSTM-like family with an exon difference).

- [ ] **Step 2: Load in Bandage; verify the arm**

Load `gwfam61.exon.gfa`, apply `gwfam61.exon.gfa.colours.csv`. Confirm the copy-specific exon is an arm one copy walks that the grey REFERENCE walk skips. No commit (verification only). Then run the final whole-branch review.

---

## Self-Review

**Spec coverage:**
- Admission carries identity → Task 1. ✓
- `copy_map_identity` + `copy_introns` aligned at all sites → Task 2 (incl. the admit_novel_pools copy_spans gap fix). ✓
- MI tag → Task 3. ✓
- ExonClass/CopyExonPath/ExonGraph types → Task 4. ✓
- Exon-class clustering (reciprocal-overlap) → Task 5. ✓
- REFERENCE walk (shared classes, fallback all) + copy walks + union-of-adjacencies L-lines (no dangling) + per-class RC → Task 6. ✓
- Colours (arm by owner / shared grey) + legend → Task 7. ✓
- `build_exon_graph` wiring + `.exon.gfa` under `--phase`, default byte-identical → Task 8. ✓
- Demo render → Task 9. ✓
- Out-of-scope (read W-lines, in-genome self-identity, FamilyGraph revival) → excluded. ✓

**Placeholder scan:** Task 5's `from_copies` has a documented two-pass sketch; the implementer note directs to the chrom-tagged 5-tuple (used consistently in Tasks 6–8). No TBD/TODO; every code step has complete code.

**Type consistency:** `ExonGraph`/`ExonClass`/`CopyExonPath`, `from_copies` (chrom-tagged 5-tuple), `to_gfa(exon_seq)`, `colours_csv`, `legend_tsv`, `build_exon_graph`, `copy_status` (shared helper), `copy_introns`/`copy_map_identity` are stable across Tasks 4–8. `Admission::Copy(DenovoTranscript, Option<f64>)` consistent Tasks 1–2. One flagged item: Task 8 relies on extracting v1's status derivation into a shared `copy_status(fa, ci)` — the implementer must factor it out of `build_copy_graph` (a small refactor, keep `build_copy_graph` behavior identical).
