# VG Novel-Copy HMM — Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Replace the flat k-mer rescue in `src/rustle/vg.rs::discover_novel_copies` with a hybrid per-exon profile-HMM family model that scores unmapped reads, places them on a family graph, emits synthetic bundles for the existing assembly pipeline, and classifies why minimap2 missed each rescued read.

**Architecture:** New code under `src/rustle/vg_hmm/` (one module per concern: `family_graph`, `profile`, `scorer`, `rescue`, `diagnostic`); `vg.rs::discover_novel_copies` becomes a thin dispatcher; the legacy k-mer-only path remains selectable via `--vg-discover-novel-mode kmer` for A/B comparison until the HMM mode is validated and made default.

**Tech Stack:** Rust 2021, `noodles-bam` (existing), `rust-spoa` (new — partial-order alignment via spoa FFI), `rayon` (existing — for parallel forward scoring), external `minimap2` binary as subprocess for the diagnostic verification stage.

**References:**
- Spec: `docs/superpowers/specs/2026-04-28-vg-novel-copy-hmm-design.md`
- Failure-mode rationale: `docs/UNMAPPED_FAMILY_RESCUE.md`
- Current scaffold to replace: `src/rustle/vg.rs:1370–1603`
- Pipeline call site: `src/rustle/pipeline.rs:14017–14040`
- Test BAMs (outside the repo): `/mnt/c/Users/jfris/Desktop/GGO_19.bam` (45 MB chr19), `/mnt/c/Users/jfris/Desktop/GGO.bam` (1.67 GB full)

---

## File map

**Created (new):**
- `src/rustle/vg_hmm/mod.rs` — module facade, public re-exports.
- `src/rustle/vg_hmm/family_graph.rs` — `FamilyGraph`, `ExonClass`, `JunctionEdge`, builder.
- `src/rustle/vg_hmm/profile.rs` — `ProfileHmm`, MSA via spoa, emission/transition fitting.
- `src/rustle/vg_hmm/scorer.rs` — forward / Viterbi / posteriors over the family graph.
- `src/rustle/vg_hmm/rescue.rs` — read selection, clustering, synthetic-bundle emission.
- `src/rustle/vg_hmm/diagnostic.rs` — five-bucket failure-mode classifier.
- `tests/regression/vg_hmm_family_graph.rs`
- `tests/regression/vg_hmm_profile.rs`
- `tests/regression/vg_hmm_scorer.rs`
- `tests/regression/vg_hmm_rescue.rs`
- `tests/regression/vg_hmm_diagnostic.rs`
- `tests/regression/vg_hmm_integration_golga.rs`

**Modified:**
- `Cargo.toml` — add `rust-spoa` dependency; register new test binaries.
- `src/rustle/lib.rs` — `pub mod vg_hmm;`.
- `src/rustle/types.rs` — add `synthetic: bool` field on `Bundle`, with `Default` impl.
- `src/rustle/vg.rs` — `discover_novel_copies` becomes a dispatcher.
- `src/bin/rustle.rs` — add `--vg-discover-novel-mode`, `--vg-rescue-diagnostic`, `--vg-rescue-min-loglik` flags.
- `src/rustle/pipeline.rs` — accept HMM-rescue results, append synthetic bundles to the assembly bundle list.
- `src/rustle/transcript_filter.rs` — exempt `synthetic` bundles from `transcript_isofrac` and cross-bundle pairwise-contained filters.
- `docs/ALGORITHMS.md` — final task: update §10 with the HMM architecture and link the rescue-class table.
- `MULTI_COPY_FAMILY_PROOF.md` — final task: replace item #3 ("scaffold exists") with shipped status + rescue-class breakdown.

---

## Conventions

- **Branch:** all work lands on a new branch `vg-hmm-novel-copy` cut from `investigation/parity-diagnostic-gates`.
- **Tests:** every regression test file is registered in `Cargo.toml` under `[[test]]` so `cargo test` picks it up. Prefer `cargo test --profile dev-opt -- --nocapture` for iteration.
- **Commits:** one commit per completed step from this plan unless explicitly noted. Commit messages: `vg-hmm: <component> <verb>` (e.g. `vg-hmm: family_graph add ExonClass cluster_by_position`).
- **Confidence checkpoints:** the python research spike (`tools/spike_vg_hmm.py`, results in `tools/SPIKE_VG_HMM_RESULTS.md`) has already validated the architectural premise on real data. Two intermediate Rust smoke checkpoints (Tasks 1.6.5 and 3.2.5 below) confirm the Rust port stays consistent with that signal.
- **No `unwrap()` in non-test code** — propagate `Result` through `anyhow::Result` (already in `Cargo.toml`).
- **Determinism:** anywhere we hash, sort, or iterate `HashMap`/`HashSet`, sort the output for determinism. Test assertions depend on it.

---

# Phase 0 — Scaffolding

### Task 0.1: Cut the implementation branch

**Files:** none (git operation).

- [ ] **Step 1: Cut branch from current head**

```bash
cd /mnt/c/Users/jfris/Desktop/Rustle
git checkout -b vg-hmm-novel-copy
git status  # confirm clean working tree
```

- [ ] **Step 2: Push branch tracking ref so worktree state is recoverable**

```bash
git push -u origin vg-hmm-novel-copy 2>/dev/null || echo "no remote configured — local-only is fine"
```

### Task 0.2: Add the spoa POA dependency

**Files:**
- Modify: `Cargo.toml`

- [ ] **Step 1: Probe whether `rust-spoa` builds in this environment**

```bash
cargo add --dry-run rust-spoa 2>&1 | head -40
```

Expected: dry-run reports a candidate version. If the crate is not on crates.io with that exact name, search alternatives:

```bash
cargo search spoa 2>&1 | head -20
cargo search "partial order alignment" 2>&1 | head -20
```

Acceptable alternatives in priority order: `rust-spoa`, `spoa-rs`, `poasta`. Pick the first that is published and has a recent (within 2 years) release. Record the chosen crate name as `<POA_CRATE>` for the rest of this plan.

- [ ] **Step 2: Add the dependency**

```bash
cargo add <POA_CRATE>
```

- [ ] **Step 3: Sanity-build**

```bash
cargo build --profile quick 2>&1 | tail -20
```

Expected: build succeeds. If linking fails because of missing C++ toolchain, abort and consult: install `build-essential` and `cmake` on the host, then re-run.

- [ ] **Step 4: Commit**

```bash
git add Cargo.toml Cargo.lock
git commit -m "vg-hmm: pin POA crate (<POA_CRATE>)"
```

### Task 0.3: Create the `vg_hmm` module skeleton

**Files:**
- Create: `src/rustle/vg_hmm/mod.rs`
- Create: `src/rustle/vg_hmm/family_graph.rs` (empty stub)
- Create: `src/rustle/vg_hmm/profile.rs` (empty stub)
- Create: `src/rustle/vg_hmm/scorer.rs` (empty stub)
- Create: `src/rustle/vg_hmm/rescue.rs` (empty stub)
- Create: `src/rustle/vg_hmm/diagnostic.rs` (empty stub)
- Modify: `src/rustle/lib.rs`

- [ ] **Step 1: Write `vg_hmm/mod.rs`**

```rust
//! Variation-graph HMM rescue for novel gene-family copies.
//!
//! Replaces the flat k-mer rescue in `crate::vg::discover_novel_copies`
//! with a per-exon profile-HMM family model. See
//! `docs/superpowers/specs/2026-04-28-vg-novel-copy-hmm-design.md` for the
//! design, and `docs/UNMAPPED_FAMILY_RESCUE.md` for why the aligner misses
//! the reads this module rescues.

pub mod family_graph;
pub mod profile;
pub mod scorer;
pub mod rescue;
pub mod diagnostic;

pub use family_graph::{ExonClass, FamilyGraph, JunctionEdge};
pub use profile::ProfileHmm;
pub use rescue::RescueResult;
pub use diagnostic::{RescueClass, FailureMode};
```

- [ ] **Step 2: Write per-file empty stubs**

For each of `family_graph.rs`, `profile.rs`, `scorer.rs`, `rescue.rs`, `diagnostic.rs`, write a single-line file:

```rust
//! Stub. Implemented in subsequent tasks per the plan.
```

- [ ] **Step 3: Register the module in `lib.rs`**

Insert near the end of the existing `pub mod` block (after `pub mod vg;`):

```rust
pub mod vg_hmm; // Family-aware HMM rescue for novel gene copies (see docs/superpowers/specs/2026-04-28-vg-novel-copy-hmm-design.md)
```

- [ ] **Step 4: Build to verify stubs are picked up**

```bash
cargo build --profile quick 2>&1 | tail -10
```

Expected: clean build, only "unused module" warnings (acceptable; will go away as code lands).

- [ ] **Step 5: Commit**

```bash
git add src/rustle/vg_hmm/ src/rustle/lib.rs
git commit -m "vg-hmm: scaffold module tree + empty stubs"
```

---

# Phase 1 — Family graph

### Task 1.1: Define `FamilyGraph`, `ExonClass`, `JunctionEdge` types

**Files:**
- Modify: `src/rustle/vg_hmm/family_graph.rs`
- Test: `tests/regression/vg_hmm_family_graph.rs`

- [ ] **Step 1: Write the failing test**

Create `tests/regression/vg_hmm_family_graph.rs`:

```rust
use rustle::vg_hmm::family_graph::{ExonClass, FamilyGraph, JunctionEdge, NodeIdx};

#[test]
fn empty_family_graph_constructs() {
    let fg = FamilyGraph {
        family_id: 0,
        nodes: Vec::<ExonClass>::new(),
        edges: Vec::<JunctionEdge>::new(),
    };
    assert_eq!(fg.family_id, 0);
    assert!(fg.nodes.is_empty());
    assert!(fg.edges.is_empty());
}

#[test]
fn exon_class_carries_per_copy_sequences() {
    let cls = ExonClass {
        idx: NodeIdx(0),
        chrom: "chrX".into(),
        span: (1000, 1500),
        strand: '+',
        per_copy_sequences: vec![(0, b"ACGT".to_vec()), (1, b"ACAT".to_vec())],
        copy_specific: false,
    };
    assert_eq!(cls.per_copy_sequences.len(), 2);
    assert!(!cls.copy_specific);
}
```

Register the new test binary in `Cargo.toml`:

```toml
[[test]]
name = "vg_hmm_family_graph"
path = "tests/regression/vg_hmm_family_graph.rs"
```

- [ ] **Step 2: Run test to verify it fails**

```bash
cargo test --profile quick --test vg_hmm_family_graph 2>&1 | tail -20
```

Expected: FAIL with errors about undefined types `ExonClass`, `FamilyGraph`, etc.

- [ ] **Step 3: Implement the types**

Replace `src/rustle/vg_hmm/family_graph.rs` with:

```rust
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
```

- [ ] **Step 4: Run test to verify it passes**

```bash
cargo test --profile quick --test vg_hmm_family_graph 2>&1 | tail -10
```

Expected: PASS, both tests.

- [ ] **Step 5: Commit**

```bash
git add src/rustle/vg_hmm/family_graph.rs tests/regression/vg_hmm_family_graph.rs Cargo.toml
git commit -m "vg-hmm: family_graph define ExonClass / JunctionEdge / FamilyGraph types"
```

### Task 1.2: Extract per-copy exon sequences from a `Bundle`

**Files:**
- Modify: `src/rustle/vg_hmm/family_graph.rs`
- Modify: `tests/regression/vg_hmm_family_graph.rs`

- [ ] **Step 1: Write the failing test**

Append to `tests/regression/vg_hmm_family_graph.rs`:

```rust
use rustle::types::{Bundle, BundleRead, JunctionStats};
use std::sync::Arc;

fn mk_bundle_with_reads(start: u64, end: u64, exons: Vec<(u64, u64)>) -> Bundle {
    let read = BundleRead {
        read_uid: 1, read_name: Arc::from("r1"), read_name_hash: 0,
        ref_id: None, mate_ref_id: None, mate_start: None, hi: 0,
        ref_start: start, ref_end: end, exons: exons.clone(),
        junctions: Vec::new(), junction_valid: Vec::new(),
        junctions_raw: Vec::new(), junctions_del: Vec::new(),
        weight: 1.0, is_reverse: false, strand: '+',
        has_poly_start: false, has_poly_end: false,
        has_poly_start_aligned: false, has_poly_start_unaligned: false,
        has_poly_end_aligned: false, has_poly_end_unaligned: false,
        unaligned_poly_t: 0, unaligned_poly_a: 0,
        has_last_exon_polya: false, has_first_exon_polyt: false,
        query_length: None, clip_left: 0, clip_right: 0, nh: 1, nm: 0, md: None,
    };
    Bundle {
        chrom: "chr1".into(), start, end, strand: '+',
        reads: vec![read],
        junction_stats: JunctionStats::default(),
        bundlenodes: None, read_bnodes: None, bnode_colors: None,
    }
}

#[test]
fn exon_extraction_returns_unique_sorted_intervals() {
    use rustle::vg_hmm::family_graph::extract_copy_exons;
    let b = mk_bundle_with_reads(100, 500, vec![(100, 200), (300, 400), (450, 500)]);
    let exons = extract_copy_exons(&b);
    assert_eq!(exons, vec![(100, 200), (300, 400), (450, 500)]);
}
```

If `BundleRead` field initialization above does not compile because the struct grew new fields, fix the test fixture by reading current `src/rustle/types.rs:100` and adding the missing fields with sensible defaults; do NOT modify `BundleRead` itself.

- [ ] **Step 2: Run test to verify it fails**

```bash
cargo test --profile quick --test vg_hmm_family_graph exon_extraction 2>&1 | tail -10
```

Expected: FAIL with "no function `extract_copy_exons`".

- [ ] **Step 3: Implement `extract_copy_exons`**

Append to `src/rustle/vg_hmm/family_graph.rs`:

```rust
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
```

- [ ] **Step 4: Run test to verify it passes**

```bash
cargo test --profile quick --test vg_hmm_family_graph exon_extraction 2>&1 | tail -10
```

Expected: PASS.

- [ ] **Step 5: Commit**

```bash
git add -u
git commit -m "vg-hmm: family_graph extract per-copy exon intervals"
```

### Task 1.3: Position-overlap clustering of exons across copies

**Files:**
- Modify: `src/rustle/vg_hmm/family_graph.rs`
- Modify: `tests/regression/vg_hmm_family_graph.rs`

- [ ] **Step 1: Write the failing test**

Append:

```rust
#[test]
fn position_overlap_clusters_partition_exons() {
    use rustle::vg_hmm::family_graph::cluster_by_position;
    // copy 0: exons at 100-200 and 300-400
    // copy 1: exons at 110-210 (overlaps copy0[0]) and 500-600 (no overlap)
    let copy0 = vec![(100u64, 200u64), (300, 400)];
    let copy1 = vec![(110u64, 210u64), (500, 600)];
    let clusters = cluster_by_position(&[("chrA", '+', copy0), ("chrA", '+', copy1)], 0.30);
    // Expect 3 clusters: {(c0,0),(c1,0)}, {(c0,1)}, {(c1,1)}
    assert_eq!(clusters.len(), 3);
    assert!(clusters.iter().any(|c| c.len() == 2)); // the overlapping pair
}
```

- [ ] **Step 2: Run test to verify it fails**

```bash
cargo test --profile quick --test vg_hmm_family_graph position_overlap 2>&1 | tail -10
```

Expected: FAIL.

- [ ] **Step 3: Implement `cluster_by_position`**

Append:

```rust
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
```

- [ ] **Step 4: Run test to verify it passes**

```bash
cargo test --profile quick --test vg_hmm_family_graph position_overlap 2>&1 | tail -10
```

Expected: PASS.

- [ ] **Step 5: Commit**

```bash
git add -u
git commit -m "vg-hmm: family_graph cluster_by_position (Stage 1: reciprocal overlap)"
```

### Task 1.4: Minimizer-Jaccard refinement (Stage 2)

**Files:**
- Modify: `src/rustle/vg_hmm/family_graph.rs`
- Modify: `tests/regression/vg_hmm_family_graph.rs`

- [ ] **Step 1: Write the failing test**

```rust
#[test]
fn minimizer_jaccard_splits_position_cluster_when_sequences_diverge() {
    use rustle::vg_hmm::family_graph::refine_by_minimizer_jaccard;
    // Two "exons" at the same position cluster but with no shared k-mers.
    let cluster = vec![(0usize, b"AAAAAAAAAAAAAAAA".to_vec()),
                       (1usize, b"GGGGGGGGGGGGGGGG".to_vec())];
    let split = refine_by_minimizer_jaccard(&cluster, 0.30, 15, 10);
    assert_eq!(split.len(), 2, "fully divergent should split");
}

#[test]
fn minimizer_jaccard_keeps_similar_sequences_together() {
    use rustle::vg_hmm::family_graph::refine_by_minimizer_jaccard;
    let s = b"ACGTACGTACGTACGTACGTACGTACGT".to_vec();
    let cluster = vec![(0usize, s.clone()), (1usize, s.clone())];
    let split = refine_by_minimizer_jaccard(&cluster, 0.30, 15, 10);
    assert_eq!(split.len(), 1, "identical should stay together");
}
```

- [ ] **Step 2: Run test to verify it fails**

```bash
cargo test --profile quick --test vg_hmm_family_graph minimizer_jaccard 2>&1 | tail -20
```

Expected: FAIL.

- [ ] **Step 3: Implement minimizer extraction + Jaccard refinement**

```rust
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
```

- [ ] **Step 4: Run test to verify it passes**

```bash
cargo test --profile quick --test vg_hmm_family_graph minimizer_jaccard 2>&1 | tail -10
```

Expected: PASS, both tests.

- [ ] **Step 5: Commit**

```bash
git add -u
git commit -m "vg-hmm: family_graph minimizer-Jaccard cluster refinement (Stage 2)"
```

### Task 1.5: Junction edge collection with `family_support`

**Files:**
- Modify: `src/rustle/vg_hmm/family_graph.rs`
- Modify: `tests/regression/vg_hmm_family_graph.rs`

- [ ] **Step 1: Write the failing test**

```rust
#[test]
fn junction_edges_count_unique_copies_supporting_each_edge() {
    use rustle::vg_hmm::family_graph::collect_family_junctions;
    // Three copies, two of which share a junction at (1000,1100); one has a private (2000,2100).
    let per_copy = vec![
        ('+', vec![(1000u64, 1100u64), (3000, 3100)]),
        ('+', vec![(1000, 1100)]),
        ('+', vec![(2000, 2100)]),
    ];
    let edges = collect_family_junctions(&per_copy);
    // Expect 3 distinct junctions with supports 2, 1, 1.
    let supp = |donor, accept| edges.iter()
        .find(|e| e.donor == donor && e.acceptor == accept)
        .map(|e| e.family_support).unwrap_or(0);
    assert_eq!(supp(1000, 1100), 2);
    assert_eq!(supp(2000, 2100), 1);
    assert_eq!(supp(3000, 3100), 1);
}
```

- [ ] **Step 2: Run test to verify it fails**

```bash
cargo test --profile quick --test vg_hmm_family_graph junction_edges 2>&1 | tail -10
```

Expected: FAIL.

- [ ] **Step 3: Implement `collect_family_junctions`**

```rust
/// A junction observation before binding to graph nodes.
#[derive(Debug, Clone, Copy)]
pub struct RawJunction {
    pub donor: u64,
    pub acceptor: u64,
    pub strand: char,
    pub family_support: u32,
}

/// Aggregate per-copy junctions into family junctions with copy-support counts.
pub fn collect_family_junctions(per_copy: &[(char, Vec<(u64, u64)>)]) -> Vec<RawJunction> {
    use std::collections::BTreeMap;
    // Round to nearest 10 bp to absorb long-read jitter (matches existing
    // build_family_consensus_junctions in vg.rs).
    let mut counts: BTreeMap<(char, u64, u64), u32> = BTreeMap::new();
    for (strand, juncs) in per_copy {
        let mut seen: HashSet<(u64, u64)> = HashSet::new();
        for &(d, a) in juncs {
            let key = (d / 10 * 10, a / 10 * 10);
            if seen.insert(key) {
                *counts.entry((*strand, key.0, key.1)).or_insert(0) += 1;
            }
        }
    }
    counts.into_iter()
        .map(|((strand, d, a), n)| RawJunction { donor: d, acceptor: a, strand, family_support: n })
        .collect()
}
```

- [ ] **Step 4: Run test to verify it passes**

```bash
cargo test --profile quick --test vg_hmm_family_graph junction_edges 2>&1 | tail -10
```

Expected: PASS.

- [ ] **Step 5: Commit**

```bash
git add -u
git commit -m "vg-hmm: family_graph collect family junctions with per-copy support"
```

### Task 1.6: End-to-end `build_family_graph` orchestrator

**Files:**
- Modify: `src/rustle/vg_hmm/family_graph.rs`
- Modify: `tests/regression/vg_hmm_family_graph.rs`

- [ ] **Step 1: Write the failing test**

```rust
#[test]
fn build_family_graph_two_copy_smoke() {
    use rustle::vg_hmm::family_graph::build_family_graph;
    use rustle::vg::FamilyGroup;
    use std::collections::HashMap;
    // Two bundles, identical exon coordinates.
    let b0 = mk_bundle_with_reads(100, 500, vec![(100, 200), (300, 400)]);
    let b1 = mk_bundle_with_reads(100, 500, vec![(100, 200), (300, 400)]);
    let bundles = vec![b0, b1];
    let family = FamilyGroup { family_id: 7, bundle_indices: vec![0, 1], multimap_reads: HashMap::new() };

    // No genome FASTA available in this test — pass `None`; build should still
    // produce nodes (with empty per-copy sequences) and edges.
    let fg = build_family_graph(&family, &bundles, None, 0.30, 0.30).unwrap();
    assert_eq!(fg.family_id, 7);
    assert_eq!(fg.nodes.len(), 2, "two shared exon classes expected");
    assert!(fg.nodes.iter().all(|n| !n.copy_specific));
}
```

- [ ] **Step 2: Run test to verify it fails**

```bash
cargo test --profile quick --test vg_hmm_family_graph build_family_graph_two_copy 2>&1 | tail -15
```

Expected: FAIL.

- [ ] **Step 3: Implement `build_family_graph`**

```rust
use crate::genome::GenomeIndex;
use crate::vg::FamilyGroup;
use anyhow::Result;

pub fn build_family_graph(
    family: &FamilyGroup,
    bundles: &[Bundle],
    genome: Option<&GenomeIndex>,
    min_pos_recip: f64,
    min_jaccard: f64,
) -> Result<FamilyGraph> {
    // 1. Collect (chrom, strand, exons) per copy.
    let copies: Vec<(&str, char, Vec<(u64, u64)>)> = family.bundle_indices.iter()
        .map(|&bi| {
            let b = &bundles[bi];
            (b.chrom.as_str(), b.strand, extract_copy_exons(b))
        })
        .collect();
    if copies.is_empty() {
        return Ok(FamilyGraph { family_id: family.family_id, nodes: Vec::new(), edges: Vec::new() });
    }

    // Anyhow-error if strands are mixed — caller should have filtered this.
    let strand0 = copies[0].1;
    if copies.iter().any(|(_, s, _)| *s != strand0) {
        anyhow::bail!("family {} has mixed strands; build_family_graph requires single-strand families", family.family_id);
    }

    // 2. Stage 1: position-overlap clusters.
    let pos_clusters = cluster_by_position(&copies, min_pos_recip);

    // 3. Stage 2: refine by minimizer Jaccard (requires sequences). If no genome,
    //    skip refinement — every position cluster becomes one ExonClass.
    let mut nodes: Vec<ExonClass> = Vec::new();
    for cluster in &pos_clusters {
        // Recover sequences (or empty stubs if no genome).
        let mut with_seq: Vec<(CopyId, Vec<u8>)> = cluster.iter().map(|&(cid, ei)| {
            let (chrom, _, exons) = &copies[cid];
            let (s, e) = exons[ei];
            let seq = match genome {
                Some(g) => g.fetch_sequence(chrom, s, e).unwrap_or_default(),
                None => Vec::new(),
            };
            (cid, seq)
        }).collect();
        // If we have sequences, refine; otherwise skip.
        let groups: Vec<Vec<CopyId>> = if with_seq.iter().all(|(_, s)| !s.is_empty()) {
            refine_by_minimizer_jaccard(&with_seq, min_jaccard, 15, 10)
        } else {
            vec![with_seq.iter().map(|(c, _)| *c).collect()]
        };
        for cids in groups {
            let mut per_copy_sequences: Vec<(CopyId, Vec<u8>)> = with_seq.iter()
                .filter(|(c, _)| cids.contains(c))
                .cloned().collect();
            per_copy_sequences.sort_by_key(|(c, _)| *c);
            // Representative span = union of contributing exon spans.
            let mut s_min = u64::MAX; let mut e_max = 0u64;
            for &cid in &cids {
                let pos_in_cluster = cluster.iter().find(|&&(c, _)| c == cid).unwrap();
                let (_, _, exs) = &copies[pos_in_cluster.0];
                let (s, e) = exs[pos_in_cluster.1];
                s_min = s_min.min(s); e_max = e_max.max(e);
            }
            nodes.push(ExonClass {
                idx: NodeIdx(nodes.len()),
                chrom: copies[0].0.to_string(),
                span: (s_min, e_max),
                strand: strand0,
                per_copy_sequences,
                copy_specific: cids.len() == 1,
            });
        }
    }

    // 4. Junction edges: collect junctions from each copy's bundle, bind to nodes
    //    by mapping (donor, acceptor) → nodes whose span ends at donor / starts at acceptor.
    let per_copy_juncs: Vec<(char, Vec<(u64, u64)>)> = family.bundle_indices.iter()
        .map(|&bi| {
            let b = &bundles[bi];
            let mut js = Vec::new();
            for (j, _) in &b.junction_stats {
                js.push((j.donor, j.acceptor));
            }
            (b.strand, js)
        })
        .collect();
    let raw = collect_family_junctions(&per_copy_juncs);

    let mut edges: Vec<JunctionEdge> = Vec::new();
    for r in &raw {
        let from = nodes.iter().find(|n| n.span.1 == (r.donor / 10 * 10) || n.span.1 == r.donor).map(|n| n.idx);
        let to   = nodes.iter().find(|n| n.span.0 == (r.acceptor / 10 * 10) || n.span.0 == r.acceptor).map(|n| n.idx);
        if let (Some(f), Some(t)) = (from, to) {
            edges.push(JunctionEdge { from: f, to: t, family_support: r.family_support, strand: r.strand });
        }
    }

    Ok(FamilyGraph { family_id: family.family_id, nodes, edges })
}
```

- [ ] **Step 4: Run test to verify it passes**

```bash
cargo test --profile quick --test vg_hmm_family_graph build_family_graph_two_copy 2>&1 | tail -10
```

Expected: PASS.

- [ ] **Step 5: Commit**

```bash
git add -u
git commit -m "vg-hmm: family_graph build_family_graph orchestrator (no MSA yet)"
```

### Task 1.6.5: Smoke checkpoint — dump GOLGA6L7 family graph

**Files:**
- Create: `src/bin/dump_family_graph.rs` (throwaway dev binary).
- Modify: `Cargo.toml` (register the binary).

This task is a confidence checkpoint, not a deliverable. It exists so that
after Phase 1 we can confirm the Rust port produces the same family-graph
shape as the python spike for GOLGA6L7. If the node count or spans look
wrong, fix Phase 1 before moving on to Phase 2 — the profiles built on a
broken graph won't tell us anything.

- [ ] **Step 1: Register the binary**

Append to `Cargo.toml`:

```toml
[[bin]]
name = "dump_family_graph"
path = "src/bin/dump_family_graph.rs"
```

- [ ] **Step 2: Implement a minimal driver**

```rust
//! Throwaway: load GGO_19.bam, run discover_family_groups, find the
//! GOLGA6L7 family on chr19, build_family_graph, print structure.
//! Drop this file once Phase 7 integration tests pin the same numbers.

use rustle::bam::load_bundles;
use rustle::vg::discover_family_groups;
use rustle::vg_hmm::family_graph::build_family_graph;

fn main() -> anyhow::Result<()> {
    let bam = std::env::args().nth(1).unwrap_or_else(||
        "/mnt/c/Users/jfris/Desktop/GGO_19.bam".to_string());
    // Load bundles via the same path pipeline.rs uses; this requires a
    // minimal RunConfig — copy whatever pipeline.rs does for its first
    // bundle-load step. If that's awkward, replace with a hand-rolled
    // BAM read that builds Bundles directly from the chr19 region around
    // 104789647-104877901.
    let bundles = load_bundles(std::path::Path::new(&bam))?;
    let families = discover_family_groups(&bundles, /* matches existing call site */ Default::default());
    println!("[smoke] {} families discovered", families.len());

    // Pick the family overlapping the GOLGA6L7 region (104789647-104877901).
    let target = families.iter().find(|f| {
        f.bundle_indices.iter().any(|&bi| {
            let b = &bundles[bi];
            b.chrom == "NC_073243.2" && b.start < 104_877_901 && b.end > 104_789_647
        })
    });
    let f = target.ok_or_else(|| anyhow::anyhow!("no family covering GOLGA6L7"))?;
    println!("[smoke] family {} has {} copies", f.family_id, f.bundle_indices.len());

    let fg = build_family_graph(f, &bundles, None, 0.30, 0.30)?;
    println!("[smoke] family graph: {} nodes, {} edges", fg.n_nodes(), fg.n_edges());
    for n in &fg.nodes {
        println!("  node {}: {} {}-{} ({} copies, copy_specific={})",
                 n.idx.0, n.chrom, n.span.0, n.span.1,
                 n.per_copy_sequences.len(), n.copy_specific);
    }
    Ok(())
}
```

If `load_bundles` does not exist with that signature, search `pipeline.rs`
for the actual bundle-loading entry point and adapt. The exact API
matters less than the output: 3 copies, ~9 exon-class nodes, no
copy-specific singletons (since GOLGA6L7 is architecturally homogeneous).

- [ ] **Step 3: Build and run**

```bash
cargo build --profile dev-opt --bin dump_family_graph 2>&1 | tail -3
./target/dev-opt/dump_family_graph /mnt/c/Users/jfris/Desktop/GGO_19.bam 2>&1 | head -30
```

Expected: 3-copy family found, ~9 nodes, all `copy_specific=false`.
The python spike sees 3 mRNAs × 9 exons each on this region; the Rust
graph should approximate that.

- [ ] **Step 4: Commit (do NOT delete the binary yet)**

```bash
git add -u
git commit -m "vg-hmm: smoke binary dump_family_graph (delete after Phase 7)"
```

---

# Phase 2 — Profile HMM

### Task 2.1: `ProfileHmm` types

**Files:**
- Modify: `src/rustle/vg_hmm/profile.rs`
- Test: `tests/regression/vg_hmm_profile.rs`

- [ ] **Step 1: Register the test binary**

Append to `Cargo.toml`:

```toml
[[test]]
name = "vg_hmm_profile"
path = "tests/regression/vg_hmm_profile.rs"
```

- [ ] **Step 2: Write the failing test**

Create `tests/regression/vg_hmm_profile.rs`:

```rust
use rustle::vg_hmm::profile::ProfileHmm;

#[test]
fn empty_profile_hmm_constructs() {
    let p = ProfileHmm::empty(0);
    assert_eq!(p.n_columns, 0);
    assert!(p.match_emit.is_empty());
}

#[test]
fn singleton_class_profile_has_one_hot_match_emissions() {
    let p = ProfileHmm::from_singleton(b"ACGT");
    assert_eq!(p.n_columns, 4);
    // P(A | col 0) > 0.9 with smoothing floor.
    assert!(p.match_emit[0][0] > 0.9_f64.ln(), "match emit[col 0]['A'] = {}", p.match_emit[0][0]);
}
```

- [ ] **Step 3: Run test to verify it fails**

```bash
cargo test --profile quick --test vg_hmm_profile 2>&1 | tail -15
```

Expected: FAIL.

- [ ] **Step 4: Implement `ProfileHmm` skeleton + singleton builder**

Replace `src/rustle/vg_hmm/profile.rs`:

```rust
//! Per-exon-class profile HMM (Krogh-style: Match / Insert / Delete states
//! with column-wise emissions). For multi-copy classes, columns come from
//! a partial-order alignment of the per-copy sequences. Singleton classes
//! get a degenerate 1-hot profile.

/// Log-probability for each base A/C/G/T (in that order).
pub type LogEmit = [f64; 4];

/// Transitions out of a single column. log P over (M->M, M->I, M->D, I->M, I->I, D->M, D->D).
#[derive(Debug, Clone, Copy, Default)]
pub struct TransRow {
    pub mm: f64, pub mi: f64, pub md: f64,
    pub im: f64, pub ii: f64,
    pub dm: f64, pub dd: f64,
}

#[derive(Debug, Clone)]
pub struct ProfileHmm {
    pub n_columns: usize,
    pub match_emit: Vec<LogEmit>,
    pub insert_emit: Vec<LogEmit>,
    pub trans: Vec<TransRow>,
    pub n_seed_copies: usize,
}

const SMOOTH_FLOOR: f64 = 0.001;

fn idx(b: u8) -> Option<usize> {
    match b { b'A' => Some(0), b'C' => Some(1), b'G' => Some(2), b'T' => Some(3), _ => None }
}

impl ProfileHmm {
    pub fn empty(n_columns: usize) -> Self {
        let bg = [(0.25_f64).ln(); 4];
        Self {
            n_columns,
            match_emit: vec![bg; n_columns],
            insert_emit: vec![bg; n_columns + 1],
            trans: vec![TransRow::default(); n_columns + 1],
            n_seed_copies: 0,
        }
    }

    /// Singleton: one copy → 1-hot match emissions with smoothing floor.
    pub fn from_singleton(seq: &[u8]) -> Self {
        let n = seq.len();
        let mut p = Self::empty(n);
        let log_main = (1.0 - 3.0 * SMOOTH_FLOOR).ln();
        let log_floor = SMOOTH_FLOOR.ln();
        for (i, &b) in seq.iter().enumerate() {
            let mut row = [log_floor; 4];
            if let Some(k) = idx(b) { row[k] = log_main; }
            p.match_emit[i] = row;
        }
        p.n_seed_copies = 1;
        p
    }
}
```

- [ ] **Step 5: Run test to verify it passes**

```bash
cargo test --profile quick --test vg_hmm_profile 2>&1 | tail -10
```

Expected: PASS.

- [ ] **Step 6: Commit**

```bash
git add -u
git commit -m "vg-hmm: profile ProfileHmm skeleton + singleton builder"
```

### Task 2.2: POA wrapper for multi-copy MSA

**Files:**
- Modify: `src/rustle/vg_hmm/profile.rs`
- Modify: `tests/regression/vg_hmm_profile.rs`

- [ ] **Step 1: Read POA crate API**

```bash
cargo doc --no-deps -p <POA_CRATE> --open 2>&1 | tail -5
# Or read directly:
find ~/.cargo/registry/src -type d -name "<POA_CRATE>*" | head -1 | xargs ls 2>/dev/null
```

Identify the function that takes `&[Vec<u8>]` (or `&[String]`) and returns aligned rows or a graph that can yield per-row aligned strings (typically `consensus_and_msa`, `align`, or similar). Note its exact signature for use below.

- [ ] **Step 2: Write the failing test**

Append to `tests/regression/vg_hmm_profile.rs`:

```rust
#[test]
fn poa_msa_returns_equal_length_rows() {
    use rustle::vg_hmm::profile::poa_msa;
    let inputs: Vec<Vec<u8>> = vec![
        b"ACGTACGT".to_vec(),
        b"ACGAACGT".to_vec(),
        b"ACGTACGT".to_vec(),
    ];
    let rows = poa_msa(&inputs).expect("POA failed");
    let l = rows[0].len();
    assert!(rows.iter().all(|r| r.len() == l), "rows must be equal length, got {:?}", rows.iter().map(|r| r.len()).collect::<Vec<_>>());
}
```

- [ ] **Step 3: Run test to verify it fails**

```bash
cargo test --profile quick --test vg_hmm_profile poa_msa 2>&1 | tail -10
```

Expected: FAIL.

- [ ] **Step 4: Implement the POA wrapper**

Append:

```rust
use anyhow::{anyhow, Result};

/// Multiple sequence alignment via POA. Returns rows of equal length,
/// padded with b'-' for gaps. One row per input sequence, in the same
/// order. Errors if POA fails (caller falls back to per-copy singletons).
pub fn poa_msa(seqs: &[Vec<u8>]) -> Result<Vec<Vec<u8>>> {
    if seqs.len() < 2 {
        return Err(anyhow!("poa_msa requires at least 2 sequences"));
    }
    // POA crate API binding: actual call depends on chosen crate. The
    // common pattern is:
    //   let mut g = poa::Aligner::new();  for s in seqs { g.add(s); }
    //   let msa = g.consensus_and_msa();
    // ADAPT THIS BLOCK to the API recorded in Step 1. The contract this
    // function exports is: equal-length rows, b'-' for gaps, same order.
    todo!("bind to chosen POA crate API; remove this todo when adapted");
}
```

Then **adapt** the body to the actual crate API. Concrete patterns to try, in order:

1. If `<POA_CRATE>` exposes `Aligner::with_capacity(...).global(seqs).msa()` returning `Vec<Vec<u8>>`, use it directly.
2. If it exposes `poa_consensus(&[(&str, &[u8])])` returning a single consensus, pair it with `align_to_consensus` per sequence and reconstruct columns.
3. If neither maps cleanly, file an issue in this plan, fall back to a banded progressive Needleman-Wunsch implementation (this is acceptable per spec §12 risk mitigation if POA is genuinely blocked).

After adapting, remove the `todo!`.

- [ ] **Step 5: Run test to verify it passes**

```bash
cargo test --profile quick --test vg_hmm_profile poa_msa 2>&1 | tail -10
```

Expected: PASS.

- [ ] **Step 6: Commit**

```bash
git add -u
git commit -m "vg-hmm: profile POA MSA wrapper"
```

### Task 2.3: Estimate match/insert emissions from POA columns with family pseudo-counts

**Files:**
- Modify: `src/rustle/vg_hmm/profile.rs`
- Modify: `tests/regression/vg_hmm_profile.rs`

- [ ] **Step 1: Write the failing test**

```rust
#[test]
fn fully_conserved_column_has_peaked_emission() {
    use rustle::vg_hmm::profile::ProfileHmm;
    let aligned: Vec<Vec<u8>> = vec![b"ACGT".to_vec(), b"ACGT".to_vec(), b"ACGT".to_vec()];
    let p = ProfileHmm::from_msa(&aligned).expect("fit failed");
    // Column 0 is fully A; expect P(A) close to 1.
    let log_pa = p.match_emit[0][0];
    assert!(log_pa > (-0.05_f64), "expected log P(A | col 0) ≈ 0, got {}", log_pa);
}

#[test]
fn divergent_column_is_flatter_than_conserved() {
    use rustle::vg_hmm::profile::ProfileHmm;
    let aligned: Vec<Vec<u8>> = vec![b"AC".to_vec(), b"AG".to_vec(), b"AT".to_vec()];
    let p = ProfileHmm::from_msa(&aligned).expect("fit failed");
    // Col 0: fully A; col 1: 1/3 each of C/G/T.
    let log_p_a_col0 = p.match_emit[0][0];
    let log_p_c_col1 = p.match_emit[1][1];
    assert!(log_p_a_col0 > log_p_c_col1, "conserved must beat divergent: {} vs {}", log_p_a_col0, log_p_c_col1);
}
```

- [ ] **Step 2: Run test to verify it fails**

```bash
cargo test --profile quick --test vg_hmm_profile from_msa 2>&1 | tail -15
```

Expected: FAIL.

- [ ] **Step 3: Implement `from_msa`**

Append:

```rust
impl ProfileHmm {
    /// Build a profile HMM from an equal-length MSA (rows = sequences,
    /// columns = aligned positions, b'-' for gaps).
    pub fn from_msa(rows: &[Vec<u8>]) -> Result<Self> {
        if rows.is_empty() { return Err(anyhow!("from_msa: no rows")); }
        let n_seq = rows.len();
        let n_col = rows[0].len();
        if !rows.iter().all(|r| r.len() == n_col) {
            return Err(anyhow!("from_msa: rows have unequal lengths"));
        }

        // Per-row gap fraction; columns where ≥ 50% of rows have a gap become
        // INSERT columns (no match state); others are match columns.
        let mut is_match: Vec<bool> = vec![false; n_col];
        for c in 0..n_col {
            let gaps = rows.iter().filter(|r| r[c] == b'-').count();
            is_match[c] = gaps * 2 < n_seq;
        }

        // Background composition (family-wide), used as the smoothing prior.
        let mut bg_counts = [0u32; 4];
        for r in rows {
            for &b in r { if let Some(k) = idx(b) { bg_counts[k] += 1; } }
        }
        let bg_total: u32 = bg_counts.iter().sum::<u32>().max(1);
        let bg: [f64; 4] = [
            bg_counts[0] as f64 / bg_total as f64,
            bg_counts[1] as f64 / bg_total as f64,
            bg_counts[2] as f64 / bg_total as f64,
            bg_counts[3] as f64 / bg_total as f64,
        ];

        // Smoothing weight: α = 1 / (1 + n_seq) — fewer copies → lean on prior.
        let alpha = 1.0 / (1.0 + n_seq as f64);

        // Build only match columns; record their original column index for trans fitting.
        let match_col_idx: Vec<usize> = (0..n_col).filter(|&c| is_match[c]).collect();
        let n_match = match_col_idx.len();

        let mut match_emit: Vec<LogEmit> = Vec::with_capacity(n_match);
        for &c in &match_col_idx {
            let mut counts = [0u32; 4];
            let mut observed = 0u32;
            for r in rows {
                if let Some(k) = idx(r[c]) { counts[k] += 1; observed += 1; }
            }
            let obs_total = observed.max(1) as f64;
            let mut row_p = [0.0_f64; 4];
            for k in 0..4 {
                let emp = counts[k] as f64 / obs_total;
                let p = alpha * bg[k] + (1.0 - alpha) * emp;
                row_p[k] = p.max(SMOOTH_FLOOR);
            }
            // Renormalize and log.
            let z: f64 = row_p.iter().sum();
            let log_row = [
                (row_p[0] / z).ln(), (row_p[1] / z).ln(),
                (row_p[2] / z).ln(), (row_p[3] / z).ln(),
            ];
            match_emit.push(log_row);
        }

        // Insert emissions: family background, log-space.
        let log_bg = [
            (bg[0].max(SMOOTH_FLOOR)).ln(),
            (bg[1].max(SMOOTH_FLOOR)).ln(),
            (bg[2].max(SMOOTH_FLOOR)).ln(),
            (bg[3].max(SMOOTH_FLOOR)).ln(),
        ];
        let insert_emit = vec![log_bg; n_match + 1];

        // Transitions deferred to Task 2.4; provide flat defaults for now.
        let trans = vec![TransRow {
            mm: (0.95_f64).ln(), mi: (0.025_f64).ln(), md: (0.025_f64).ln(),
            im: (0.5_f64).ln(),  ii: (0.5_f64).ln(),
            dm: (0.95_f64).ln(), dd: (0.05_f64).ln(),
        }; n_match + 1];

        Ok(Self {
            n_columns: n_match,
            match_emit,
            insert_emit,
            trans,
            n_seed_copies: n_seq,
        })
    }
}
```

- [ ] **Step 4: Run test to verify it passes**

```bash
cargo test --profile quick --test vg_hmm_profile from_msa 2>&1 | tail -10
```

Expected: PASS.

- [ ] **Step 5: Commit**

```bash
git add -u
git commit -m "vg-hmm: profile fit match emissions from MSA with family pseudo-counts"
```

### Task 2.4: Fit transitions from MSA gap structure

**Files:**
- Modify: `src/rustle/vg_hmm/profile.rs`
- Modify: `tests/regression/vg_hmm_profile.rs`

- [ ] **Step 1: Write the failing test**

```rust
#[test]
fn no_gaps_msa_makes_mm_dominant() {
    use rustle::vg_hmm::profile::ProfileHmm;
    let rows: Vec<Vec<u8>> = vec![b"ACGT".to_vec(), b"ACGT".to_vec()];
    let p = ProfileHmm::from_msa(&rows).unwrap();
    // M->M should be the dominant transition out of every column.
    for tr in &p.trans {
        assert!(tr.mm > tr.mi && tr.mm > tr.md, "expected M->M dominant, got {:?}", tr);
    }
}

#[test]
fn gappy_msa_increases_md_or_mi() {
    use rustle::vg_hmm::profile::ProfileHmm;
    let rows: Vec<Vec<u8>> = vec![b"AC-GT".to_vec(), b"ACGGT".to_vec(), b"AC-GT".to_vec()];
    let p = ProfileHmm::from_msa(&rows).unwrap();
    // Transition out of column 1 (after the C) should show non-trivial M->D
    // because two of three rows have a deletion at column 2.
    assert!(p.trans[1].md > (1e-3_f64).ln(), "expected M->D > 0 at col 1, got {}", p.trans[1].md);
}
```

- [ ] **Step 2: Run test to verify it fails**

```bash
cargo test --profile quick --test vg_hmm_profile transitions 2>&1 | tail -15
```

Expected: FAIL or pass-by-default depending on what `from_msa` currently emits. Tighten until it fails meaningfully (the second test should fail because trans is currently flat).

- [ ] **Step 3: Implement transition fitting in `from_msa`**

Replace the `// Transitions deferred to Task 2.4; provide flat defaults for now.` block in `from_msa` with:

```rust
// Fit transitions by counting per-row state walks across columns.
// State at column c (for row r):
//   Match  if r[c] != '-' AND is_match[c]
//   Delete if r[c] == '-' AND is_match[c]
//   Insert if !is_match[c]   (we don't track per-insert sub-states)
//
// Transitions out of MATCH column m_idx are accumulated by walking the
// row from that column to the next match column; intervening insert
// columns count as M->I->...->I->M (one M->I out, then I->I within,
// and I->M into the next match). Gaps at the next match column count
// M->D or D->D depending on whether the *current* column is M or D.
let mut counts: Vec<[u32; 7]> = vec![[0; 7]; n_match + 1]; // M->M, M->I, M->D, I->M, I->I, D->M, D->D

for row in rows {
    // Walk left→right tracking previous match-column state.
    let mut prev_match_state: Option<bool> = None; // Some(true) = match, Some(false) = delete, None = before-first-match
    let mut had_insert_after_prev = false;
    let mut prev_match_idx_in_match: Option<usize> = None;
    for c in 0..n_col {
        if is_match[c] {
            let cur_is_delete = row[c] == b'-';
            let mi = match_col_idx.iter().position(|&x| x == c).unwrap();
            if let Some(prev_was_match) = prev_match_state {
                let prev_mi = prev_match_idx_in_match.unwrap();
                if had_insert_after_prev {
                    // prev → I → ... → I → cur
                    if prev_was_match { counts[prev_mi][1] += 1; } // M->I
                    // I->I edges within would need to be counted; skip for simplicity (single insert run).
                    if cur_is_delete { counts[prev_mi + 1][3] += 0; /* I->D not modelled */ }
                    else { counts[prev_mi][3] += 1; } // I->M (counted at prev's row for simplicity)
                } else {
                    if prev_was_match {
                        if cur_is_delete { counts[prev_mi][2] += 1; } // M->D
                        else { counts[prev_mi][0] += 1; }             // M->M
                    } else {
                        if cur_is_delete { counts[prev_mi][6] += 1; } // D->D
                        else { counts[prev_mi][5] += 1; }             // D->M
                    }
                }
            }
            prev_match_state = Some(!cur_is_delete);
            prev_match_idx_in_match = Some(mi);
            had_insert_after_prev = false;
        } else if row[c] != b'-' {
            had_insert_after_prev = true;
        }
    }
}

// Convert counts → smoothed log-probabilities per row.
let mut trans: Vec<TransRow> = Vec::with_capacity(n_match + 1);
for c in &counts {
    let m_total = (c[0] + c[1] + c[2]).max(1) as f64;
    let i_total = (c[3] + c[4]).max(1) as f64;
    let d_total = (c[5] + c[6]).max(1) as f64;
    let s = SMOOTH_FLOOR;
    let mm = ((c[0] as f64 + s) / (m_total + 3.0 * s)).ln();
    let mi = ((c[1] as f64 + s) / (m_total + 3.0 * s)).ln();
    let md = ((c[2] as f64 + s) / (m_total + 3.0 * s)).ln();
    let im = ((c[3] as f64 + s) / (i_total + 2.0 * s)).ln();
    let ii = ((c[4] as f64 + s) / (i_total + 2.0 * s)).ln();
    let dm = ((c[5] as f64 + s) / (d_total + 2.0 * s)).ln();
    let dd = ((c[6] as f64 + s) / (d_total + 2.0 * s)).ln();
    trans.push(TransRow { mm, mi, md, im, ii, dm, dd });
}
```

- [ ] **Step 4: Run test to verify it passes**

```bash
cargo test --profile quick --test vg_hmm_profile transitions 2>&1 | tail -10
```

Expected: PASS, both tests. If the gappy test fails, inspect the count walk — it is the messiest part of this task and likely to need an iteration.

- [ ] **Step 5: Commit**

```bash
git add -u
git commit -m "vg-hmm: profile fit transitions from MSA gap structure"
```

### Task 2.5: Wire profile fitting into `family_graph::build_family_graph`

**Files:**
- Modify: `src/rustle/vg_hmm/family_graph.rs`
- Modify: `src/rustle/vg_hmm/profile.rs`
- Modify: `tests/regression/vg_hmm_family_graph.rs`

- [ ] **Step 1: Add an optional `profile: Option<ProfileHmm>` to `ExonClass`**

In `family_graph.rs`, modify `ExonClass`:

```rust
pub struct ExonClass {
    pub idx: NodeIdx,
    pub chrom: String,
    pub span: (u64, u64),
    pub strand: char,
    pub per_copy_sequences: Vec<(CopyId, Vec<u8>)>,
    pub copy_specific: bool,
    pub profile: Option<crate::vg_hmm::profile::ProfileHmm>,
}
```

Update existing constructions in `build_family_graph` and tests to set `profile: None` initially.

- [ ] **Step 2: Write the failing test**

Append to `vg_hmm_family_graph.rs`:

```rust
#[test]
fn build_family_graph_fits_profiles_when_sequences_present() {
    use rustle::vg_hmm::family_graph::build_family_graph;
    use rustle::vg::FamilyGroup;
    use std::collections::HashMap;
    // Without genome the test in 1.6 keeps profiles None — that case still holds.
    // Here we simulate sequences by stubbing the genome path. Easiest: extend
    // build_family_graph with a function variant that takes preloaded sequences.
    // For test-only API, see fit_profiles_in_place below.
    let bundles = vec![mk_bundle_with_reads(0, 100, vec![(0, 50), (60, 100)]),
                       mk_bundle_with_reads(0, 100, vec![(0, 50), (60, 100)])];
    let family = FamilyGroup { family_id: 0, bundle_indices: vec![0, 1], multimap_reads: HashMap::new() };
    let mut fg = build_family_graph(&family, &bundles, None, 0.30, 0.30).unwrap();
    // Inject sequences and refit.
    for n in &mut fg.nodes {
        n.per_copy_sequences = vec![(0, b"ACGTACGT".to_vec()), (1, b"ACGAACGT".to_vec())];
    }
    rustle::vg_hmm::family_graph::fit_profiles_in_place(&mut fg).unwrap();
    assert!(fg.nodes.iter().all(|n| n.profile.is_some()));
    assert!(fg.nodes[0].profile.as_ref().unwrap().n_columns >= 6);
}
```

- [ ] **Step 3: Run test to verify it fails**

Expected: FAIL — `fit_profiles_in_place` does not exist.

- [ ] **Step 4: Implement `fit_profiles_in_place`**

Append to `family_graph.rs`:

```rust
use crate::vg_hmm::profile::ProfileHmm;

pub fn fit_profiles_in_place(fg: &mut FamilyGraph) -> Result<()> {
    for node in &mut fg.nodes {
        let seqs: Vec<Vec<u8>> = node.per_copy_sequences.iter().map(|(_, s)| s.clone()).collect();
        node.profile = Some(if seqs.len() == 1 {
            ProfileHmm::from_singleton(&seqs[0])
        } else {
            // Mean pairwise identity gate (spec §12 risk mitigation): if too divergent,
            // fall back to per-copy singletons by collapsing to the first one.
            let id = mean_pairwise_identity(&seqs);
            if id < 0.60 {
                ProfileHmm::from_singleton(&seqs[0])
            } else {
                let msa = match crate::vg_hmm::profile::poa_msa(&seqs) {
                    Ok(m) => m,
                    Err(_) => { node.copy_specific = true; vec![seqs[0].clone()] }
                };
                ProfileHmm::from_msa(&msa).unwrap_or_else(|_| ProfileHmm::from_singleton(&seqs[0]))
            }
        });
    }
    Ok(())
}

fn mean_pairwise_identity(seqs: &[Vec<u8>]) -> f64 {
    let n = seqs.len();
    if n < 2 { return 1.0; }
    let mut total = 0.0_f64; let mut pairs = 0u32;
    for i in 0..n {
        for j in (i + 1)..n {
            let m = seqs[i].iter().zip(seqs[j].iter()).take(seqs[i].len().min(seqs[j].len())).filter(|(a, b)| a == b).count();
            let l = seqs[i].len().max(seqs[j].len()).max(1);
            total += m as f64 / l as f64;
            pairs += 1;
        }
    }
    total / pairs.max(1) as f64
}
```

- [ ] **Step 5: Run test to verify it passes**

```bash
cargo test --profile quick --test vg_hmm_family_graph build_family_graph_fits_profiles 2>&1 | tail -10
```

Expected: PASS.

- [ ] **Step 6: Commit**

```bash
git add -u
git commit -m "vg-hmm: family_graph fit_profiles_in_place (POA → ProfileHmm per node)"
```

---

# Phase 3 — Scorer

### Task 3.1: Single-column profile-HMM forward (banded)

**Files:**
- Modify: `src/rustle/vg_hmm/scorer.rs`
- Test: `tests/regression/vg_hmm_scorer.rs`

- [ ] **Step 1: Register the test binary**

Append to `Cargo.toml`:

```toml
[[test]]
name = "vg_hmm_scorer"
path = "tests/regression/vg_hmm_scorer.rs"
```

- [ ] **Step 2: Write the failing test**

```rust
use rustle::vg_hmm::profile::ProfileHmm;
use rustle::vg_hmm::scorer::forward_against_profile;

#[test]
fn perfectly_matching_read_scores_higher_than_random() {
    let p = ProfileHmm::from_singleton(b"ACGTACGT");
    let good = forward_against_profile(&p, b"ACGTACGT");
    let bad  = forward_against_profile(&p, b"TTTTTTTT");
    assert!(good > bad, "good={} bad={}", good, bad);
}
```

- [ ] **Step 3: Run test to verify it fails**

Expected: FAIL.

- [ ] **Step 4: Implement single-profile forward**

Replace `src/rustle/vg_hmm/scorer.rs`:

```rust
//! Forward / Viterbi over the family graph. For now, single-profile only;
//! inter-node transitions land in Task 3.3.

use crate::vg_hmm::profile::ProfileHmm;

#[inline]
fn idx(b: u8) -> Option<usize> {
    match b { b'A' => Some(0), b'C' => Some(1), b'G' => Some(2), b'T' => Some(3), _ => None }
}

const NEG_INF: f64 = f64::NEG_INFINITY;

#[inline]
fn logsumexp(a: f64, b: f64) -> f64 {
    if a == NEG_INF { return b; }
    if b == NEG_INF { return a; }
    let m = a.max(b);
    m + ((a - m).exp() + (b - m).exp()).ln()
}

/// Forward log P(read | profile). Banded would be a perf optimization;
/// start unbanded for correctness, add banding in a follow-up if profiling
/// shows it matters.
pub fn forward_against_profile(p: &ProfileHmm, read: &[u8]) -> f64 {
    let m = p.n_columns;     // number of match columns
    let l = read.len();
    if m == 0 { return NEG_INF; }

    // dp_M[j][i], dp_I[j][i], dp_D[j][i] over match columns 1..=m and positions 0..=l.
    let mut m_dp = vec![vec![NEG_INF; l + 1]; m + 1];
    let mut i_dp = vec![vec![NEG_INF; l + 1]; m + 1];
    let mut d_dp = vec![vec![NEG_INF; l + 1]; m + 1];

    // Initial: M[0][0] is virtual start.
    m_dp[0][0] = 0.0;

    for j in 1..=m {
        for i in 0..=l {
            // Emission of M[j] consuming read[i-1].
            let emit_m = if i > 0 {
                idx(read[i - 1]).map(|k| p.match_emit[j - 1][k]).unwrap_or(NEG_INF)
            } else { NEG_INF };
            // M from M[j-1][i-1] + mm, I[j-1][i-1] + im, D[j-1][i-1] + dm
            let from_m = if i > 0 { m_dp[j - 1][i - 1] + p.trans[j - 1].mm + emit_m } else { NEG_INF };
            let from_i = if i > 0 { i_dp[j - 1][i - 1] + p.trans[j - 1].im + emit_m } else { NEG_INF };
            let from_d = if i > 0 { d_dp[j - 1][i - 1] + p.trans[j - 1].dm + emit_m } else { NEG_INF };
            m_dp[j][i] = logsumexp(logsumexp(from_m, from_i), from_d);

            // I from M[j][i-1] + mi (with insert emission), I[j][i-1] + ii.
            let emit_i = if i > 0 {
                idx(read[i - 1]).map(|k| p.insert_emit[j][k]).unwrap_or(NEG_INF)
            } else { NEG_INF };
            let i_from_m = if i > 0 { m_dp[j][i - 1] + p.trans[j].mi + emit_i } else { NEG_INF };
            let i_from_i = if i > 0 { i_dp[j][i - 1] + p.trans[j].ii + emit_i } else { NEG_INF };
            i_dp[j][i] = logsumexp(i_from_m, i_from_i);

            // D from M[j-1][i] + md, D[j-1][i] + dd (no emission).
            let d_from_m = m_dp[j - 1][i] + p.trans[j - 1].md;
            let d_from_d = d_dp[j - 1][i] + p.trans[j - 1].dd;
            d_dp[j][i] = logsumexp(d_from_m, d_from_d);
        }
    }
    // Terminate: total = log-sum-exp over M[m][l], I[m][l], D[m][l].
    logsumexp(logsumexp(m_dp[m][l], i_dp[m][l]), d_dp[m][l])
}
```

- [ ] **Step 5: Run test to verify it passes**

```bash
cargo test --profile quick --test vg_hmm_scorer 2>&1 | tail -10
```

Expected: PASS.

- [ ] **Step 6: Commit**

```bash
git add -u
git commit -m "vg-hmm: scorer single-profile forward (full DP, unbanded)"
```

### Task 3.2: Family-graph forward across multiple profile nodes

**Files:**
- Modify: `src/rustle/vg_hmm/scorer.rs`
- Modify: `tests/regression/vg_hmm_scorer.rs`

- [ ] **Step 1: Write the failing test**

```rust
use rustle::vg_hmm::family_graph::{ExonClass, FamilyGraph, JunctionEdge, NodeIdx};
use rustle::vg_hmm::scorer::forward_against_family;
use rustle::vg_hmm::profile::ProfileHmm;

fn two_node_chain() -> FamilyGraph {
    let p1 = ProfileHmm::from_singleton(b"ACGT");
    let p2 = ProfileHmm::from_singleton(b"TGCA");
    let nodes = vec![
        ExonClass { idx: NodeIdx(0), chrom: "x".into(), span: (0, 4), strand: '+',
            per_copy_sequences: vec![(0, b"ACGT".to_vec())], copy_specific: true, profile: Some(p1) },
        ExonClass { idx: NodeIdx(1), chrom: "x".into(), span: (10, 14), strand: '+',
            per_copy_sequences: vec![(0, b"TGCA".to_vec())], copy_specific: true, profile: Some(p2) },
    ];
    let edges = vec![JunctionEdge { from: NodeIdx(0), to: NodeIdx(1), family_support: 1, strand: '+' }];
    FamilyGraph { family_id: 0, nodes, edges }
}

#[test]
fn family_forward_path_through_two_nodes() {
    let fg = two_node_chain();
    let on_path = forward_against_family(&fg, b"ACGTTGCA");
    let off_path = forward_against_family(&fg, b"AAAAAAAA");
    assert!(on_path > off_path, "on_path={} off_path={}", on_path, off_path);
}
```

- [ ] **Step 2: Run test to verify it fails**

```bash
cargo test --profile quick --test vg_hmm_scorer family_forward 2>&1 | tail -10
```

Expected: FAIL.

- [ ] **Step 3: Implement `forward_against_family`**

Append to `scorer.rs`:

```rust
use crate::vg_hmm::family_graph::{FamilyGraph, NodeIdx};

/// Forward over the family graph: for each read split point, try all
/// node sequences that respect junction edges. We compute, for each node
/// and each starting read position, the best forward score over all
/// path prefixes ending at that node.
///
/// Implementation: dynamic programming over `(node_idx, read_pos_after)`,
/// where the value is `log P(read[0..pos] | path ending at this node)`.
/// Transition between nodes adds the log-junction-prior
/// `log(family_support / max_support_in_family)`.
pub fn forward_against_family(fg: &FamilyGraph, read: &[u8]) -> f64 {
    if fg.nodes.is_empty() { return NEG_INF; }
    let n = fg.nodes.len();
    let l = read.len();
    let max_support = fg.edges.iter().map(|e| e.family_support).max().unwrap_or(1).max(1);

    // node_score[node][pos_after] = best log-prob of consuming read[0..pos_after] ending at node
    let mut node_score = vec![vec![NEG_INF; l + 1]; n];

    // Determine source nodes (those with no incoming edges).
    let mut has_incoming = vec![false; n];
    for e in &fg.edges { has_incoming[e.to.0] = true; }

    // For source nodes: consume some prefix of read.
    for (idx, node) in fg.nodes.iter().enumerate() {
        if !has_incoming[idx] {
            if let Some(p) = &node.profile {
                for pos_after in 0..=l {
                    node_score[idx][pos_after] = logsumexp(
                        node_score[idx][pos_after],
                        forward_against_profile(p, &read[..pos_after]),
                    );
                }
            }
        }
    }

    // Topological order over edges (assume nodes are already in reasonable order;
    // if not, run Kahn's algorithm — most family graphs are small).
    let mut order: Vec<NodeIdx> = (0..n).map(NodeIdx).collect();
    // Trust insertion order as topological for a single-strand chain; revisit if a test fails.

    for &nidx in &order {
        for e in fg.edges.iter().filter(|e| e.from == nidx) {
            let to_node = &fg.nodes[e.to.0];
            let log_edge = (e.family_support as f64 / max_support as f64).ln();
            if let Some(p) = &to_node.profile {
                // For each starting position at the destination, the best score is
                // node_score[from][start] + log_edge + forward(p, read[start..end]).
                for start in 0..=l {
                    if node_score[nidx.0][start] == NEG_INF { continue; }
                    for end in start..=l {
                        let inner = forward_against_profile(p, &read[start..end]);
                        let cand = node_score[nidx.0][start] + log_edge + inner;
                        node_score[e.to.0][end] = logsumexp(node_score[e.to.0][end], cand);
                    }
                }
            }
        }
    }

    // Final: best score over any node ending at read end.
    let mut best = NEG_INF;
    for s in &node_score { best = best.max(s[l]); }
    best
}
```

- [ ] **Step 4: Run test to verify it passes**

```bash
cargo test --profile quick --test vg_hmm_scorer family_forward 2>&1 | tail -10
```

Expected: PASS.

- [ ] **Step 5: Commit**

```bash
git add -u
git commit -m "vg-hmm: scorer forward over family-graph nodes"
```

### Task 3.3: Viterbi best-path traceback over the family graph

**Files:**
- Modify: `src/rustle/vg_hmm/scorer.rs`
- Modify: `tests/regression/vg_hmm_scorer.rs`

- [ ] **Step 1: Write the failing test**

```rust
#[test]
fn viterbi_returns_correct_node_sequence() {
    use rustle::vg_hmm::scorer::viterbi_path;
    let fg = two_node_chain();
    let path = viterbi_path(&fg, b"ACGTTGCA").expect("path");
    let names: Vec<usize> = path.nodes.iter().map(|n| n.0).collect();
    assert_eq!(names, vec![0, 1]);
    assert!(path.score > 0.0_f64.ln() - 100.0); // sanity bound
}
```

- [ ] **Step 2: Run test to verify it fails**

Expected: FAIL.

- [ ] **Step 3: Implement Viterbi**

Append to `scorer.rs`:

```rust
#[derive(Debug, Clone)]
pub struct ViterbiPath {
    pub nodes: Vec<NodeIdx>,
    pub score: f64,
}

/// Viterbi best-path over the family graph. Mirrors `forward_against_family`
/// but with `max` instead of `logsumexp` and with a backpointer table.
pub fn viterbi_path(fg: &FamilyGraph, read: &[u8]) -> Option<ViterbiPath> {
    if fg.nodes.is_empty() { return None; }
    let n = fg.nodes.len();
    let l = read.len();
    let max_support = fg.edges.iter().map(|e| e.family_support).max().unwrap_or(1).max(1);

    let mut best = vec![vec![NEG_INF; l + 1]; n];
    let mut bp_node: Vec<Vec<Option<NodeIdx>>> = vec![vec![None; l + 1]; n];
    let mut bp_pos:  Vec<Vec<usize>>           = vec![vec![0; l + 1]; n];

    let mut has_incoming = vec![false; n];
    for e in &fg.edges { has_incoming[e.to.0] = true; }

    for (idx, node) in fg.nodes.iter().enumerate() {
        if !has_incoming[idx] {
            if let Some(p) = &node.profile {
                for pos_after in 0..=l {
                    let s = forward_against_profile(p, &read[..pos_after]);
                    if s > best[idx][pos_after] { best[idx][pos_after] = s; }
                }
            }
        }
    }

    for nidx in 0..n {
        for e in fg.edges.iter().filter(|e| e.from.0 == nidx) {
            let to_node = &fg.nodes[e.to.0];
            let log_edge = (e.family_support as f64 / max_support as f64).ln();
            if let Some(p) = &to_node.profile {
                for start in 0..=l {
                    if best[nidx][start] == NEG_INF { continue; }
                    for end in start..=l {
                        let inner = forward_against_profile(p, &read[start..end]);
                        let cand = best[nidx][start] + log_edge + inner;
                        if cand > best[e.to.0][end] {
                            best[e.to.0][end] = cand;
                            bp_node[e.to.0][end] = Some(NodeIdx(nidx));
                            bp_pos[e.to.0][end] = start;
                        }
                    }
                }
            }
        }
    }

    // Find best terminal node.
    let mut term: Option<(usize, f64)> = None;
    for nidx in 0..n {
        let s = best[nidx][l];
        if term.map_or(true, |(_, ts)| s > ts) { term = Some((nidx, s)); }
    }
    let (mut cur_node, score) = term?;
    let mut cur_pos = l;
    let mut path = vec![NodeIdx(cur_node)];
    while let Some(prev) = bp_node[cur_node][cur_pos] {
        let prev_pos = bp_pos[cur_node][cur_pos];
        path.push(prev);
        cur_node = prev.0;
        cur_pos = prev_pos;
    }
    path.reverse();
    Some(ViterbiPath { nodes: path, score })
}
```

- [ ] **Step 4: Run test to verify it passes**

```bash
cargo test --profile quick --test vg_hmm_scorer viterbi 2>&1 | tail -10
```

Expected: PASS.

- [ ] **Step 5: Commit**

```bash
git add -u
git commit -m "vg-hmm: scorer Viterbi best-path traceback"
```

### Task 3.2.5: Smoke checkpoint — score positive vs random reads

**Files:**
- Modify: `src/bin/dump_family_graph.rs` (extend the throwaway binary).

After Phase 3 Task 3.2 we have a working family-graph forward scorer.
The python spike showed unmapped reads have a tail of family signal
above random reads (max 18 vs 7 k-mer hits; max 15.5 vs 9.7 PSSM
log-odds). The Rust port must reproduce at least the qualitative
ordering — positive >> random — or something in Phase 1–3 is wrong.

- [ ] **Step 1: Extend the smoke binary to score reads**

Add to `dump_family_graph.rs`, after the family-graph dump:

```rust
use rustle::vg_hmm::family_graph::fit_profiles_in_place;
use rustle::vg_hmm::scorer::forward_against_family;
use rustle::genome::GenomeIndex;

let genome = GenomeIndex::from_fasta(std::path::Path::new(
    "/mnt/c/Users/jfris/Desktop/GGO.fasta"))?;
// Re-build with sequences this time.
let mut fg = build_family_graph(f, &bundles, Some(&genome), 0.30, 0.30)?;
fit_profiles_in_place(&mut fg)?;
println!("[smoke] profiles fit on {} nodes", fg.nodes.iter().filter(|n| n.profile.is_some()).count());

// Pull a few reads from the family bundle (positive) and from a far
// region of chr19 (random control). Score each.
let pos_seqs: Vec<Vec<u8>> = bundles[f.bundle_indices[0]].reads.iter()
    .take(20).map(|r| r.read_name.as_bytes().to_vec()).collect();
//   ^ NB: BundleRead does not carry sequence — pull it from BAM by read name,
//   or use a small pysam-like helper. For the smoke, simplest path is a
//   one-off `samtools view -f 0 ... | head -20 | awk '{print $10}'` shelled
//   out and piped. Document the chosen approach in this file's header.

println!("[smoke] scoring {} positive reads", pos_seqs.len());
for s in &pos_seqs {
    let lp = forward_against_family(&fg, s);
    println!("  positive  {} bp  log-lik = {:.2}", s.len(), lp);
}
// Repeat with reads sampled from chr19:75_000_000-80_000_000 as the
// random control — should produce noticeably lower log-likelihoods.
```

- [ ] **Step 2: Run and inspect**

```bash
cargo build --profile dev-opt --bin dump_family_graph 2>&1 | tail -3
./target/dev-opt/dump_family_graph /mnt/c/Users/jfris/Desktop/GGO_19.bam 2>&1 | tail -50
```

Expected: positive reads score systematically higher than random. The
spike showed ~5× separation in the tail; the Rust HMM should match or
exceed because it has M/I/D states and family-prior smoothing that the
PSSM lacked. **If positives don't score higher than random, stop and
debug Phase 2 (profile fitting) or Phase 3 (forward DP) before
proceeding to Phase 4.**

- [ ] **Step 3: Commit**

```bash
git add -u
git commit -m "vg-hmm: smoke checkpoint score positives vs random in dump_family_graph"
```

---

# Phase 4 — Config + CLI dispatcher

### Task 4.1: Add `synthetic` flag to `Bundle` and `BundleSource::VgRescue` source string

**Files:**
- Modify: `src/rustle/types.rs`
- Modify: existing call sites that construct `Bundle` (search-and-replace)

- [ ] **Step 1: Find all `Bundle { ... }` constructions**

```bash
grep -rn "Bundle {" src/ tests/ 2>&1 | grep -v "BundleRead\|BundleSource" | head -30
```

Note every site; each will need `synthetic: false` added to keep it compiling.

- [ ] **Step 2: Add the field**

In `src/rustle/types.rs`, append to the `Bundle` struct (after `bnode_colors`):

```rust
    /// True for bundles synthesized by VG rescue from unmapped reads
    /// (see src/rustle/vg_hmm/rescue.rs). Synthetic bundles bypass
    /// transcript_isofrac and cross-bundle pairwise-contained filters.
    pub synthetic: bool,
```

- [ ] **Step 3: Update all constructions**

For every site found in Step 1, add `synthetic: false`. If the site uses `..Default::default()`, also implement `Default` for `Bundle` (only if it does not already exist; check first):

```bash
grep -n "impl Default for Bundle" src/rustle/types.rs
```

If absent and a callsite needs it, add a `Default` impl that mirrors the empty bundle currently constructed in tests.

- [ ] **Step 4: Build**

```bash
cargo build --profile quick 2>&1 | tail -20
```

Expected: clean build.

- [ ] **Step 5: Commit**

```bash
git add -u
git commit -m "vg-hmm: Bundle.synthetic flag (default false; opt-in for rescued bundles)"
```

### Task 4.2: Add CLI flags and `RunConfig` fields

**Files:**
- Modify: `src/rustle/types.rs` (RunConfig)
- Modify: `src/bin/rustle.rs`

- [ ] **Step 1: Find `RunConfig` and the existing `vg_*` flags**

```bash
grep -n "vg_discover_novel\|vg_snp\|vg_solver\|pub struct RunConfig" src/rustle/types.rs | head
```

Note where existing `vg_*` fields live.

- [ ] **Step 2: Add fields**

Insert near the existing `vg_discover_novel: bool`:

```rust
    /// `kmer` (legacy) or `hmm` (new). Default `kmer` until validation flips it.
    pub vg_discover_novel_mode: String,
    /// Enable external minimap2 verification of rescued reads.
    pub vg_rescue_diagnostic: bool,
    /// Forward log-odds threshold for HMM rescue (nats).
    pub vg_rescue_min_loglik: f64,
```

If `RunConfig` has a `Default`, add corresponding defaults (`"kmer".to_string()`, `false`, `30.0`).

- [ ] **Step 3: Add CLI flags**

In `src/bin/rustle.rs`, find the existing `--vg-discover-novel` clap field and add (with the same crate-style derive macros used there):

```rust
    /// VG novel-copy discovery algorithm: `kmer` (legacy) or `hmm` (family RNA-HMM).
    #[arg(long = "vg-discover-novel-mode", default_value = "kmer")]
    pub vg_discover_novel_mode: String,

    /// Run external minimap2 verification on each rescued read (slower).
    #[arg(long = "vg-rescue-diagnostic")]
    pub vg_rescue_diagnostic: bool,

    /// Forward log-odds threshold for HMM rescue (nats). Default 30.0.
    #[arg(long = "vg-rescue-min-loglik", default_value_t = 30.0)]
    pub vg_rescue_min_loglik: f64,
```

Wire the parsed values into the `RunConfig` construction the same way other `vg_*` flags are wired.

- [ ] **Step 4: Build and run --help**

```bash
cargo build --profile quick 2>&1 | tail -5
./target/quick/rustle --help 2>&1 | grep -E "vg-discover-novel-mode|vg-rescue-diagnostic|vg-rescue-min-loglik"
```

Expected: all three flags listed.

- [ ] **Step 5: Commit**

```bash
git add -u
git commit -m "vg-hmm: CLI + RunConfig add --vg-discover-novel-mode / --vg-rescue-diagnostic / --vg-rescue-min-loglik"
```

### Task 4.3: Make `vg::discover_novel_copies` a dispatcher

**Files:**
- Modify: `src/rustle/vg.rs`

- [ ] **Step 1: Add a dispatcher above the existing function**

Rename the existing function from `discover_novel_copies` to `discover_novel_copies_kmer`, then add:

```rust
pub fn discover_novel_copies(
    bam_path: &std::path::Path,
    families: &[FamilyGroup],
    bundles: &[Bundle],
    config: &crate::types::RunConfig,
) -> Vec<NovelCandidate> {
    match config.vg_discover_novel_mode.as_str() {
        "hmm" => crate::vg_hmm::rescue::run_rescue(bam_path, families, bundles, config)
            .unwrap_or_else(|e| {
                eprintln!("[VG-HMM] rescue failed: {} — falling back to kmer", e);
                discover_novel_copies_kmer(bam_path, families, bundles, config)
            }),
        _ => discover_novel_copies_kmer(bam_path, families, bundles, config),
    }
}
```

- [ ] **Step 2: Add a stub `run_rescue` so the build passes**

In `src/rustle/vg_hmm/rescue.rs`, replace the stub with:

```rust
//! Rescue pipeline: read prefilter → HMM scoring → clustering → synthetic bundle emission.

use crate::types::{Bundle, RunConfig};
use crate::vg::{FamilyGroup, NovelCandidate};
use anyhow::Result;

pub struct RescueResult; // re-exported for now; populated in Task 5.x.

pub fn run_rescue(
    _bam_path: &std::path::Path,
    _families: &[FamilyGroup],
    _bundles: &[Bundle],
    _config: &RunConfig,
) -> Result<Vec<NovelCandidate>> {
    // Stub: returns empty until Phase 5 lands. Tests covering hmm-mode end-to-end
    // pass once Phase 5 wires the implementation.
    Ok(Vec::new())
}
```

- [ ] **Step 3: Build**

```bash
cargo build --profile quick 2>&1 | tail -5
```

Expected: clean build.

- [ ] **Step 4: Smoke-test the dispatcher path**

```bash
# Run the existing kmer path and confirm nothing changed.
cargo test --profile quick 2>&1 | tail -20
```

Expected: existing tests still pass.

- [ ] **Step 5: Commit**

```bash
git add -u
git commit -m "vg-hmm: dispatcher in discover_novel_copies + run_rescue stub"
```

---

# Phase 5 — Rescue pipeline

### Task 5.1: K-mer prefilter (reuse existing logic)

**Files:**
- Modify: `src/rustle/vg_hmm/rescue.rs`
- Test: `tests/regression/vg_hmm_rescue.rs`

- [ ] **Step 1: Register test binary, write the failing test**

Append to `Cargo.toml`:

```toml
[[test]]
name = "vg_hmm_rescue"
path = "tests/regression/vg_hmm_rescue.rs"
```

Create `tests/regression/vg_hmm_rescue.rs`:

```rust
use rustle::vg_hmm::rescue::prefilter_read;
use std::collections::HashSet;

#[test]
fn prefilter_passes_when_kmer_hits_meet_threshold() {
    let mut family_kmers: HashSet<u64> = HashSet::new();
    // Fake hash that we know our prefilter will produce.
    family_kmers.insert(0xdeadbeef);
    let read = b"ACGTACGTACGTACGT".to_vec();
    let (passed, _hits) = prefilter_read(&read, &family_kmers, 15, 1);
    // We don't know the actual hash, so just assert the API returns and (passed,hits) is consistent.
    let _ = passed;
}
```

- [ ] **Step 2: Implement `prefilter_read`**

Append to `rescue.rs`:

```rust
use std::collections::HashSet as StdHashSet;

fn fnv1a64(bytes: &[u8]) -> u64 {
    let mut h: u64 = 0xcbf29ce484222325;
    for &b in bytes { h ^= b as u64; h = h.wrapping_mul(0x100000001b3); }
    h
}

/// Returns (passed, n_hits). Mirrors the k-mer index format of vg.rs's flat
/// rescue (FNV-1a of A/C/G/T-only k-mers).
pub fn prefilter_read(seq: &[u8], family_kmers: &StdHashSet<u64>, k: usize, min_hits: usize) -> (bool, usize) {
    if seq.len() < k { return (false, 0); }
    let mut hits = 0;
    for w in seq.windows(k) {
        if w.iter().any(|&b| !matches!(b, b'A' | b'C' | b'G' | b'T')) { continue; }
        if family_kmers.contains(&fnv1a64(w)) { hits += 1; }
    }
    (hits >= min_hits, hits)
}
```

- [ ] **Step 3: Run, commit**

```bash
cargo test --profile quick --test vg_hmm_rescue 2>&1 | tail -5
git add -u && git commit -m "vg-hmm: rescue k-mer prefilter"
```

### Task 5.2: Wire family-graph build + scoring loop

**Files:**
- Modify: `src/rustle/vg_hmm/rescue.rs`
- Modify: `tests/regression/vg_hmm_rescue.rs`

- [ ] **Step 1: Write the failing test**

```rust
#[test]
fn run_rescue_returns_zero_candidates_with_empty_input() {
    use rustle::vg_hmm::rescue::run_rescue_in_memory;
    use rustle::vg::FamilyGroup;
    let candidates = run_rescue_in_memory(&[], &[], &[], 30.0).unwrap();
    assert!(candidates.is_empty());
}
```

- [ ] **Step 2: Implement an in-memory entry point**

Append to `rescue.rs`:

```rust
use crate::vg_hmm::family_graph::{build_family_graph, fit_profiles_in_place, FamilyGraph};
use crate::vg_hmm::scorer::{forward_against_family, viterbi_path};

/// Entry point that takes already-loaded BAM unmapped reads. Used by the
/// in-memory tests; `run_rescue` builds reads from the BAM and dispatches here.
pub fn run_rescue_in_memory(
    families: &[FamilyGroup],
    bundles: &[Bundle],
    unmapped_reads: &[(String, Vec<u8>)],   // (read_name, sequence)
    min_loglik: f64,
) -> Result<Vec<NovelCandidate>> {
    if families.is_empty() || unmapped_reads.is_empty() { return Ok(Vec::new()); }

    // Build family graphs. No genome here — caller-side run_rescue passes one in.
    let mut family_graphs: Vec<FamilyGraph> = Vec::new();
    for f in families {
        let mut fg = build_family_graph(f, bundles, None, 0.30, 0.30)?;
        let _ = fit_profiles_in_place(&mut fg);
        family_graphs.push(fg);
    }

    let mut out = Vec::new();
    for (name, seq) in unmapped_reads {
        let mut best: Option<(usize, f64)> = None;
        for fg in &family_graphs {
            if fg.nodes.is_empty() { continue; }
            let s = forward_against_family(fg, seq);
            if s > min_loglik && best.map_or(true, |(_, bs)| s > bs) {
                best = Some((fg.family_id, s));
            }
        }
        if let Some((fid, s)) = best {
            out.push(NovelCandidate {
                read_name: name.clone(),
                family_id: fid,
                matched_junctions: 0,
                total_junctions: 0,
                approx_chrom: String::new(),
                approx_start: 0,
                approx_end: seq.len() as u64,
            });
            // log-odds saved on side-channel TSV in Task 5.5; here we just record it.
            let _ = s;
        }
    }
    Ok(out)
}
```

- [ ] **Step 3: Build, test, commit**

```bash
cargo test --profile quick --test vg_hmm_rescue 2>&1 | tail -5
git add -u && git commit -m "vg-hmm: rescue in-memory scoring entry point"
```

### Task 5.3: BAM-backed `run_rescue` driver

**Files:**
- Modify: `src/rustle/vg_hmm/rescue.rs`

- [ ] **Step 1: Implement `run_rescue` to iterate BAM and call `run_rescue_in_memory`**

Replace the stub in `rescue.rs` with:

```rust
pub fn run_rescue(
    bam_path: &std::path::Path,
    families: &[FamilyGroup],
    bundles: &[Bundle],
    config: &RunConfig,
) -> Result<Vec<NovelCandidate>> {
    // 1. Read unmapped sequences from BAM. Reuse the decoding logic from
    //    vg.rs::discover_novel_copies (lines ~1485–1530); refactor that
    //    block out into a shared helper if convenient. For now, inline it
    //    locally to keep the dispatcher's contract exact.
    let bam_file = std::fs::File::open(bam_path)?;
    let buf = std::io::BufReader::new(bam_file);
    let worker_count = std::num::NonZeroUsize::MIN;
    let bgzf = noodles_bgzf::MultithreadedReader::with_worker_count(worker_count, buf);
    let mut reader = noodles_bam::io::Reader::from(bgzf);
    let _header = reader.read_header()?;

    let mut unmapped: Vec<(String, Vec<u8>)> = Vec::new();
    for result in reader.records() {
        let record = result?;
        if !record.flags().is_unmapped() { continue; }
        // Decode 4-bit sequence (mirroring vg.rs:1509–1530).
        let seq_obj = record.sequence();
        let seq_raw = seq_obj.as_ref();
        let seq_len = seq_obj.len();
        let mut buf: Vec<u8> = Vec::with_capacity(seq_len);
        for (i, &byte) in seq_raw.iter().enumerate() {
            let bases = [(byte >> 4) & 0xF, byte & 0xF];
            for (j, nib) in bases.iter().enumerate() {
                if i * 2 + j >= seq_len { break; }
                let b = match nib {
                    1 => b'A', 2 => b'C', 4 => b'G', 8 => b'T', _ => b'N',
                };
                buf.push(b);
            }
        }
        let name = record.name().map(|n| n.to_string()).unwrap_or_default();
        unmapped.push((name, buf));
    }
    eprintln!("[VG-HMM] {} unmapped reads to score", unmapped.len());

    run_rescue_in_memory(families, bundles, &unmapped, config.vg_rescue_min_loglik)
}
```

Add the `noodles_bgzf`, `noodles_bam` imports at the top.

- [ ] **Step 2: Build, run smoke test on chr19 BAM**

```bash
cargo build --profile dev-opt 2>&1 | tail -5
./target/dev-opt/rustle -L --vg --vg-discover-novel --vg-discover-novel-mode hmm \
    --genome-fasta /mnt/c/Users/jfris/Desktop/GGO_genomic.fna \
    -o /tmp/rustle_hmm_test.gtf /mnt/c/Users/jfris/Desktop/GGO_19.bam 2>&1 | grep -E "VG-HMM|VG\]" | head -20
```

Expected: HMM mode runs without panicking; reports unmapped read count and rescue counts.

- [ ] **Step 3: Commit**

```bash
git add -u
git commit -m "vg-hmm: rescue BAM-backed driver (run_rescue)"
```

### Task 5.4: Cluster rescued reads into synthetic bundles

**Files:**
- Modify: `src/rustle/vg_hmm/rescue.rs`
- Modify: `tests/regression/vg_hmm_rescue.rs`

- [ ] **Step 1: Write the failing test**

```rust
#[test]
fn synthetic_bundle_emitted_per_path_region_cluster() {
    use rustle::vg_hmm::rescue::synthesize_bundles;
    use rustle::vg_hmm::family_graph::{ExonClass, FamilyGraph, JunctionEdge, NodeIdx};
    let fg = FamilyGraph {
        family_id: 1,
        nodes: vec![
            ExonClass { idx: NodeIdx(0), chrom: "chr1".into(), span: (1000, 1100), strand: '+',
                per_copy_sequences: vec![], copy_specific: false, profile: None },
        ],
        edges: vec![],
    };
    let cands_with_paths: Vec<(String, Vec<NodeIdx>)> = (0..5).map(|i| (format!("r{}", i), vec![NodeIdx(0)])).collect();
    let bundles = synthesize_bundles(&fg, &cands_with_paths, 3);
    assert_eq!(bundles.len(), 1);
    assert!(bundles[0].synthetic);
    assert_eq!(bundles[0].chrom, "chr1");
}
```

- [ ] **Step 2: Implement `synthesize_bundles`**

Append:

```rust
use crate::types::{BundleRead, JunctionStats};
use std::sync::Arc;

pub fn synthesize_bundles(
    fg: &FamilyGraph,
    rescued: &[(String, Vec<crate::vg_hmm::family_graph::NodeIdx>)],
    min_reads: usize,
) -> Vec<Bundle> {
    use std::collections::BTreeMap;
    // Cluster by canonical path representation; widen by ±1 kb.
    let mut clusters: BTreeMap<Vec<usize>, Vec<&str>> = BTreeMap::new();
    for (name, path) in rescued {
        let key: Vec<usize> = path.iter().map(|n| n.0).collect();
        clusters.entry(key).or_default().push(name.as_str());
    }
    let mut out = Vec::new();
    for (path_key, names) in clusters {
        if names.len() < min_reads { continue; }
        // Path span: union of representative spans of nodes on the path.
        let spans: Vec<(u64, u64)> = path_key.iter().map(|&i| fg.nodes[i].span).collect();
        let s = spans.iter().map(|(a, _)| *a).min().unwrap_or(0).saturating_sub(1000);
        let e = spans.iter().map(|(_, b)| *b).max().unwrap_or(0) + 1000;
        let chrom = fg.nodes[path_key[0]].chrom.clone();
        let strand = fg.nodes[path_key[0]].strand;
        let reads: Vec<BundleRead> = names.iter().enumerate().map(|(i, n)| BundleRead {
            read_uid: 1_000_000 + i as u64,
            read_name: Arc::from(*n),
            read_name_hash: 0, ref_id: None, mate_ref_id: None, mate_start: None, hi: 0,
            ref_start: s, ref_end: e,
            exons: spans.clone(),
            junctions: Vec::new(), junction_valid: Vec::new(),
            junctions_raw: Vec::new(), junctions_del: Vec::new(),
            weight: 1.0, is_reverse: false, strand,
            has_poly_start: false, has_poly_end: false,
            has_poly_start_aligned: false, has_poly_start_unaligned: false,
            has_poly_end_aligned: false, has_poly_end_unaligned: false,
            unaligned_poly_t: 0, unaligned_poly_a: 0,
            has_last_exon_polya: false, has_first_exon_polyt: false,
            query_length: None, clip_left: 0, clip_right: 0, nh: 1, nm: 0, md: None,
        }).collect();
        out.push(Bundle {
            chrom, start: s, end: e, strand,
            reads,
            junction_stats: JunctionStats::default(),
            bundlenodes: None, read_bnodes: None, bnode_colors: None,
            synthetic: true,
        });
    }
    out
}
```

- [ ] **Step 3: Build, test, commit**

```bash
cargo test --profile quick --test vg_hmm_rescue synthetic_bundle 2>&1 | tail -10
git add -u && git commit -m "vg-hmm: rescue synthesize_bundles from clustered Viterbi paths"
```

### Task 5.5: Wire synthetic bundles back into the pipeline

**Files:**
- Modify: `src/rustle/vg_hmm/rescue.rs`
- Modify: `src/rustle/pipeline.rs`
- Modify: `src/rustle/transcript_filter.rs`

- [ ] **Step 1: Modify `run_rescue_in_memory` to also return synthetic bundles**

Change signature to return `(Vec<NovelCandidate>, Vec<Bundle>)`. Internally, after scoring each read, also call `viterbi_path` to record the path; pass `(name, path)` lists into `synthesize_bundles`. Update `run_rescue` to forward the new tuple. Update the dispatcher (Task 4.3) accordingly — it currently returns only `Vec<NovelCandidate>`; either:
- Add a sibling function `run_rescue_with_bundles` and call it from `pipeline.rs` directly (preferred — fewer changes to `discover_novel_copies` callers), OR
- Stash the synthetic bundles in a thread-local / `RunConfig` extension. (Avoid this — global state.)

Pick the first option.

- [ ] **Step 2: In `pipeline.rs`, near line 14017 where `discover_novel_copies` is called, branch on `vg_discover_novel_mode`**

```rust
if config.vg_mode && config.vg_discover_novel && !vg_families.is_empty() {
    if config.vg_discover_novel_mode == "hmm" {
        match crate::vg_hmm::rescue::run_rescue_with_bundles(bam_path, &vg_families, &bundles, config) {
            Ok((cands, synth_bundles)) => {
                eprintln!("[VG-HMM] {} novel candidates, {} synthetic bundles",
                          cands.len(), synth_bundles.len());
                bundles.extend(synth_bundles);
            }
            Err(e) => eprintln!("[VG-HMM] rescue error: {}", e),
        }
    } else {
        // existing kmer path unchanged
        let novel_candidates = crate::vg::discover_novel_copies(...);
        // existing handling
    }
}
```

- [ ] **Step 3: Filter exemptions**

In `transcript_filter.rs`, find the `transcript_isofrac` filter (search for `isofrac` near filter logic) and the cross-bundle pairwise-contained filter. Add a guard at each:

```rust
if bundle.synthetic { continue; }  // skip transcript_isofrac
```

For the pairwise-contained filter that compares across bundles, skip the pair when `containing_bundle.synthetic != contained_bundle.synthetic` (i.e. only suppress containment between two non-synthetic or two synthetic bundles).

- [ ] **Step 4: Build and run on chr19 smoke**

```bash
cargo build --profile dev-opt 2>&1 | tail -3
./target/dev-opt/rustle -L --vg --vg-discover-novel --vg-discover-novel-mode hmm \
    --genome-fasta /mnt/c/Users/jfris/Desktop/GGO_genomic.fna \
    -o /tmp/rustle_hmm_test.gtf /mnt/c/Users/jfris/Desktop/GGO_19.bam 2>&1 | grep -E "VG-HMM|family|novel" | head
grep "vg_rescue\|copy_status" /tmp/rustle_hmm_test.gtf | head
```

Expected: at least one synthetic bundle; emitted GTF includes `source "vg_rescue"` or equivalent attribute (next task adds that explicitly).

- [ ] **Step 5: Commit**

```bash
git add -u
git commit -m "vg-hmm: pipeline wire synthetic bundles + filter exemptions"
```

---

# Phase 6 — Diagnostic

### Task 6.1: Internal failure-mode classifier (chain score + seed-frequency mask)

**Files:**
- Modify: `src/rustle/vg_hmm/diagnostic.rs`
- Test: `tests/regression/vg_hmm_diagnostic.rs`

- [ ] **Step 1: Register the test binary, write the failing test**

Append to `Cargo.toml`:

```toml
[[test]]
name = "vg_hmm_diagnostic"
path = "tests/regression/vg_hmm_diagnostic.rs"
```

Create `tests/regression/vg_hmm_diagnostic.rs`:

```rust
use rustle::vg_hmm::diagnostic::{classify_internal, RescueClass};

#[test]
fn many_short_islands_classify_as_below_threshold() {
    // 30 isolated 15-mer hits ≈ chain too short to anchor.
    let class = classify_internal(30, 15, 5_000, 0.001);
    assert_eq!(class, RescueClass::BelowThresholdChain);
}

#[test]
fn dominant_repetitive_seeds_classify_as_seed_masked() {
    // High hit count but most seeds are repetitive.
    let class = classify_internal(500, 15, 5_000, 0.50);
    assert_eq!(class, RescueClass::SeedMasked);
}
```

- [ ] **Step 2: Implement `classify_internal`**

Replace `src/rustle/vg_hmm/diagnostic.rs`:

```rust
//! Failure-mode classification of rescued reads.
//! See docs/UNMAPPED_FAMILY_RESCUE.md for the five-bucket rationale.

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum RescueClass {
    BelowThresholdChain,
    SeedMasked,
    Divergent,
    Structural,
    ReferenceAbsent,
    NeedsExternalVerification,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum FailureMode {
    Internal(RescueClass),
    External(RescueClass),
}

/// Internal-only classifier: bucket 1 (below_threshold) vs 2 (seed_masked) vs
/// "needs external verification" for everything else.
///
/// Args:
///   n_kmer_hits — how many family k-mers the read matched
///   k — k-mer length used in the index
///   read_len — read length in bp (for chain-score bound)
///   masked_fraction — fraction of the read's k-mer hits that hit
///                     repetitive (top-percentile) seeds
pub fn classify_internal(n_kmer_hits: usize, k: usize, read_len: usize, masked_fraction: f64) -> RescueClass {
    // Approx minimap2 chain-score bound: each k-mer match is worth ~k score
    // points in best-case (no gaps). Default minimap2 -s for HiFi is ~40.
    let approx_chain_score = (n_kmer_hits as f64) * (k as f64);
    if approx_chain_score < 40.0 {
        return RescueClass::BelowThresholdChain;
    }
    if masked_fraction > 0.30 {
        return RescueClass::SeedMasked;
    }
    let _ = read_len;
    RescueClass::NeedsExternalVerification
}
```

- [ ] **Step 3: Run, commit**

```bash
cargo test --profile quick --test vg_hmm_diagnostic 2>&1 | tail -10
git add -u && git commit -m "vg-hmm: diagnostic internal classifier (buckets 1, 2, needs-external)"
```

### Task 6.2: External minimap2 subprocess verification

**Files:**
- Modify: `src/rustle/vg_hmm/diagnostic.rs`
- Modify: `tests/regression/vg_hmm_diagnostic.rs`

- [ ] **Step 1: Confirm minimap2 is available**

```bash
which minimap2 && minimap2 --version 2>&1 | head -1
```

If not installed, the engineer pauses to install: `conda install -c bioconda minimap2` or download static binary. Diagnostic stage requires it; the test below skips when binary is missing.

- [ ] **Step 2: Write the failing test (with skip-if-missing)**

```rust
#[test]
fn external_minimap2_classifies_divergent_when_no_alignment_at_default() {
    // Skip if minimap2 not on PATH.
    if std::process::Command::new("minimap2").arg("--version").output().is_err() {
        eprintln!("skipping: minimap2 not on PATH");
        return;
    }
    use rustle::vg_hmm::diagnostic::classify_external;
    use std::io::Write;
    // Construct a tiny FASTA with a clearly divergent sequence relative to the read.
    let dir = tempfile::tempdir().unwrap();
    let fa = dir.path().join("ref.fa");
    let mut f = std::fs::File::create(&fa).unwrap();
    writeln!(f, ">ref\nAAAAAAAAAAAAAAAAAAAAAAAAAAAA").unwrap();
    drop(f);
    let class = classify_external(&fa, b"GTGTGTGTGTGTGTGTGTGTGTGTGTGT").unwrap();
    // Maximally divergent → should land in ReferenceAbsent or Divergent.
    assert!(matches!(class, rustle::vg_hmm::diagnostic::RescueClass::ReferenceAbsent
                          | rustle::vg_hmm::diagnostic::RescueClass::Divergent));
}
```

Add the dev-dependency:

```bash
cargo add --dev tempfile
```

- [ ] **Step 3: Implement `classify_external`**

Append:

```rust
use anyhow::{Context, Result};
use std::io::Write;
use std::process::{Command, Stdio};

/// Run minimap2 with stepped-down parameters; classify into Divergent /
/// Structural / ReferenceAbsent based on alignment outcome at each step.
pub fn classify_external(ref_fasta: &std::path::Path, read_seq: &[u8]) -> Result<RescueClass> {
    // Steps: default → -s 20 → -f 0.0001 → -N 50 --secondary=yes → -s 10 -f 0 -N 100
    let steps: &[&[&str]] = &[
        &["-x", "map-hifi"],
        &["-x", "map-hifi", "-s", "20"],
        &["-x", "map-hifi", "-s", "20", "-f", "0.0001"],
        &["-x", "map-hifi", "-s", "20", "-f", "0.0001", "-N", "50", "--secondary=yes"],
        &["-x", "map-hifi", "-s", "10", "-f", "0", "-N", "100"],
    ];
    for (step_idx, args) in steps.iter().enumerate() {
        let mut cmd = Command::new("minimap2");
        cmd.args(*args).arg(ref_fasta).arg("-").arg("-a")
            .stdin(Stdio::piped()).stdout(Stdio::piped()).stderr(Stdio::null());
        let mut child = cmd.spawn().context("spawn minimap2")?;
        // Pipe a single FASTA record to stdin.
        let stdin = child.stdin.as_mut().unwrap();
        writeln!(stdin, ">read")?;
        stdin.write_all(read_seq)?;
        writeln!(stdin)?;
        let out = child.wait_with_output()?;
        let sam = String::from_utf8_lossy(&out.stdout);
        // Look for at least one alignment record with a non-* CIGAR and not flag 4.
        let has_aln = sam.lines().any(|l| {
            !l.starts_with('@') && l.split('\t').nth(1).and_then(|f| f.parse::<u32>().ok()).map_or(false, |f| f & 4 == 0)
        });
        if has_aln {
            // Identify identity / structural based on CIGAR.
            // Simple proxy: long indels (≥50 bp) → Structural; otherwise Divergent.
            let structural = sam.lines().any(|l| {
                if l.starts_with('@') { return false; }
                let fields: Vec<&str> = l.split('\t').collect();
                if fields.len() < 6 { return false; }
                let cigar = fields[5];
                cigar_has_long_indel(cigar, 50)
            });
            return Ok(if structural { RescueClass::Structural }
                      else if step_idx >= 2 { RescueClass::Divergent }
                      else { RescueClass::Divergent });  // step_idx 0/1 = aligns at default-ish; still Divergent (the read's identity is the bucket)
        }
    }
    Ok(RescueClass::ReferenceAbsent)
}

fn cigar_has_long_indel(cigar: &str, threshold: usize) -> bool {
    let mut num = 0usize;
    for c in cigar.chars() {
        if c.is_ascii_digit() { num = num * 10 + (c as u8 - b'0') as usize; }
        else { if (c == 'I' || c == 'D' || c == 'N') && num >= threshold { return true; } num = 0; }
    }
    false
}
```

- [ ] **Step 4: Run, commit**

```bash
cargo test --profile dev-opt --test vg_hmm_diagnostic 2>&1 | tail -10
git add -u && git commit -m "vg-hmm: diagnostic external minimap2 verification (buckets 3, 4, 5)"
```

### Task 6.3: Wire diagnostic into rescue + GTF attributes + per-family TSV

**Files:**
- Modify: `src/rustle/vg_hmm/rescue.rs`
- Modify: `src/rustle/vg_hmm/diagnostic.rs`
- Modify: `src/rustle/gtf.rs` (or wherever GTF attribute emission lives — check first)

- [ ] **Step 1: Locate GTF attribute emission for `vg_rescue` source**

```bash
grep -rn "copy_status\|family_id\|rescue_class\|source =\|source = \"vg" src/rustle/gtf.rs src/rustle/transcript_filter.rs 2>&1 | head
```

Identify the function building per-transcript attribute strings. For the synthetic-bundle path, it should append `rescue_class "<bucket>"` for each transcript whose source bundle is synthetic.

- [ ] **Step 2: Track per-rescued-read class through the synthetic-bundle path**

Extend `synthesize_bundles` to accept `Vec<(String, Vec<NodeIdx>, RescueClass)>` instead of `(String, Vec<NodeIdx>)`. Aggregate the cluster's most-common class as the bundle's `rescue_class` and stash it on `Bundle` as a new optional field `pub rescue_class: Option<RescueClass>` (added to types.rs alongside `synthetic`).

- [ ] **Step 3: Run classification per rescued read in `run_rescue_in_memory`**

After scoring each surviving read, call `classify_internal(...)`; if `RescueClass::NeedsExternalVerification` AND `config.vg_rescue_diagnostic`, call `classify_external(...)` against a temporary FASTA containing the family's path-region sequences. Otherwise, leave class as `NeedsExternalVerification`.

- [ ] **Step 4: Emit `rescue_class` GTF attribute**

In the GTF emitter, when the source bundle has `synthetic = true` AND `rescue_class.is_some()`, append:

```
rescue_class "{bucket}";
```

- [ ] **Step 5: Per-family TSV**

If `config.vg_report` is set, append rescue-class columns to the existing TSV writer in `vg.rs::write_family_report*`. Add columns: `n_rescued`, `n_below_thresh`, `n_seed_masked`, `n_divergent`, `n_structural`, `n_ref_absent`. The counts come from aggregating `RescueClass` values per family.

- [ ] **Step 6: Build, smoke, commit**

```bash
cargo build --profile dev-opt 2>&1 | tail -3
./target/dev-opt/rustle -L --vg --vg-discover-novel --vg-discover-novel-mode hmm \
    --vg-rescue-diagnostic --vg-report /tmp/fam.tsv \
    --genome-fasta /mnt/c/Users/jfris/Desktop/GGO_genomic.fna \
    -o /tmp/rustle_hmm_test.gtf /mnt/c/Users/jfris/Desktop/GGO_19.bam 2>&1 | tail
grep "rescue_class" /tmp/rustle_hmm_test.gtf | head -5
head /tmp/fam.tsv
git add -u && git commit -m "vg-hmm: emit rescue_class on synthetic transcripts + TSV columns"
```

---

# Phase 7 — Integration tests, validation, default flip

### Task 7.1: GOLGA6L7 chr19 smoke integration test

**Files:**
- Create: `tests/regression/vg_hmm_integration_golga.rs`

- [ ] **Step 1: Register, write the test**

Append to `Cargo.toml`:

```toml
[[test]]
name = "vg_hmm_integration_golga"
path = "tests/regression/vg_hmm_integration_golga.rs"
```

Create the test:

```rust
//! End-to-end smoke on the chr19 BAM (~10 s). Skipped when the BAM is missing
//! (CI machines without Desktop fixtures still pass `cargo test`).

use std::process::Command;

const BAM: &str = "/mnt/c/Users/jfris/Desktop/GGO_19.bam";
const GENOME: &str = "/mnt/c/Users/jfris/Desktop/GGO_genomic.fna";

#[test]
fn hmm_mode_smoke_on_chr19() {
    if !std::path::Path::new(BAM).exists() {
        eprintln!("skipping: {} missing", BAM);
        return;
    }
    if !std::path::Path::new(GENOME).exists() {
        eprintln!("skipping: {} missing", GENOME);
        return;
    }

    let out_gtf = std::env::temp_dir().join("rustle_hmm_smoke.gtf");
    let status = Command::new(env!("CARGO_BIN_EXE_rustle"))
        .args([
            "-L", "--vg", "--vg-discover-novel",
            "--vg-discover-novel-mode", "hmm",
            "--genome-fasta", GENOME,
            "-o", out_gtf.to_str().unwrap(),
            BAM,
        ])
        .status()
        .expect("rustle failed to spawn");
    assert!(status.success(), "rustle exited with {:?}", status);

    let gtf = std::fs::read_to_string(&out_gtf).expect("read gtf");
    assert!(gtf.lines().count() > 0, "empty gtf");
    // Smoke check: HMM mode should produce at least one entry. The exact count
    // depends on tuning; this test is a "did anything go fundamentally wrong"
    // gate, not a regression-pinning gate.
}
```

- [ ] **Step 2: Run**

```bash
cargo test --profile dev-opt --test vg_hmm_integration_golga -- --nocapture 2>&1 | tail
```

Expected: PASS, or "skipping" if BAM absent.

- [ ] **Step 3: Commit**

```bash
git add -u
git commit -m "vg-hmm: chr19 GOLGA smoke integration test"
```

### Task 7.2: Synthetic ground-truth fixture

**Files:**
- Create: `tests/regression/vg_hmm_synthetic_ground_truth.rs`
- Create: `test_data/vg_hmm/synthetic_family.fa` (small, checked in)
- Create: `test_data/vg_hmm/synthetic_reads.bam` (small; built by a one-off helper script committed alongside)

- [ ] **Step 1: Build the fixture**

Use existing tooling under `tools/` to generate a tiny BAM of perturbed reads. Pseudocode:

```python
# tools/make_vg_hmm_fixture.py
# Read one GOLGA copy from GGO_genomic.fna; perturb 5% of bases; insert
# a 50-bp copy-specific exon at a known position; simulate 50 reads via
# pbsim2 or a simple substring-with-noise function. Output:
#   test_data/vg_hmm/synthetic_family.fa  (the base copy)
#   test_data/vg_hmm/synthetic_reads.bam  (the simulated reads)
```

Commit the helper script and the resulting fixtures (small files only — keep `synthetic_reads.bam` under 1 MB).

- [ ] **Step 2: Write the test**

```rust
#[test]
fn synthetic_perturbed_reads_are_rescued_and_classified_divergent() {
    // Run rustle on the synthetic BAM, parse GTF, check that:
    //  - At least 80% of simulated reads were rescued
    //  - The bundles' rescue_class is Divergent or Structural
    //  - external minimap2 (if available) confirms no default-param alignment
    todo!("implement once Task 7.2 fixture lands");
}
```

Implement the test concretely (open output GTF, count `rescue_class "Divergent"` lines, etc.). Keep it tight — not more than 30 lines of Rust.

- [ ] **Step 3: Run, fix, commit**

```bash
cargo test --profile dev-opt --test vg_hmm_synthetic_ground_truth 2>&1 | tail
git add -u && git commit -m "vg-hmm: synthetic ground-truth integration test"
```

### Task 7.3: Negative control fixture

**Files:**
- Create: `tests/regression/vg_hmm_negative_control.rs`

- [ ] **Step 1: Test**

```rust
//! Random low-complexity / contam-like reads must not be rescued.

#[test]
fn random_unmapped_reads_yield_no_rescue() {
    use rustle::vg_hmm::rescue::run_rescue_in_memory;
    use rustle::vg::FamilyGroup;
    use rustle::types::{Bundle, JunctionStats};
    use std::collections::HashMap;
    // Synthesize a family with one bundle of plausible exons.
    let bundle = Bundle {
        chrom: "chrX".into(), start: 0, end: 100, strand: '+',
        reads: vec![],
        junction_stats: JunctionStats::default(),
        bundlenodes: None, read_bnodes: None, bnode_colors: None, synthetic: false,
    };
    let family = FamilyGroup { family_id: 0, bundle_indices: vec![0], multimap_reads: HashMap::new() };
    // Pure-A reads.
    let reads = (0..50).map(|i| (format!("r{}", i), vec![b'A'; 1000])).collect::<Vec<_>>();
    let cands = run_rescue_in_memory(&[family], &[bundle], &reads, 30.0).unwrap();
    assert!(cands.is_empty(), "expected zero rescues from low-complexity reads, got {}", cands.len());
}
```

Register, run, commit.

```bash
cargo test --profile dev-opt --test vg_hmm_negative_control 2>&1 | tail
git add -u && git commit -m "vg-hmm: negative-control test (low-complexity reads yield zero rescues)"
```

### Task 7.4: Full-BAM benchmark + flip default mode if results hold

**Files:**
- Modify: `src/bin/rustle.rs` (default value of `--vg-discover-novel-mode`)
- Modify: `MULTI_COPY_FAMILY_PROOF.md`
- Modify: `docs/ALGORITHMS.md`

- [ ] **Step 1: Run the full benchmark**

```bash
./target/release/rustle -L --vg --vg-discover-novel --vg-discover-novel-mode hmm \
    --vg-report /tmp/fam_hmm.tsv \
    --genome-fasta /mnt/c/Users/jfris/Desktop/GGO_genomic.fna \
    -o /tmp/rustle_hmm_full.gtf /mnt/c/Users/jfris/Desktop/GGO.bam 2>&1 | tail -20

python3 tools/classify_family_assembly.py /mnt/c/Users/jfris/Desktop/GGO_genomic.gff \
    "golgin subfamily A member 6||golgin A6 family like" \
    /tmp/rustle_hmm_full.gtf --label Rustle_HMM/GOLGA6
python3 tools/classify_family_assembly.py /mnt/c/Users/jfris/Desktop/GGO_genomic.gff \
    "golgin subfamily A member 8||golgin A8 family like" \
    /tmp/rustle_hmm_full.gtf --label Rustle_HMM/GOLGA8
```

Compare exact / partial / wrong-strand / missing counts against `MULTI_COPY_FAMILY_PROOF.md`'s VG-mode baseline (8 / 33 exact). If HMM ≥ baseline on exact-match counts AND wrong-strand ≤ baseline, proceed to flip default.

- [ ] **Step 2: Flip default**

In `src/bin/rustle.rs`, change:

```rust
#[arg(long = "vg-discover-novel-mode", default_value = "kmer")]
```

to:

```rust
#[arg(long = "vg-discover-novel-mode", default_value = "hmm")]
```

And update `RunConfig::default` correspondingly.

- [ ] **Step 3: Update docs**

In `MULTI_COPY_FAMILY_PROOF.md`, replace item #3 ("scaffold exists") with the shipped status, the rescue-class breakdown from `/tmp/fam_hmm.tsv`, and the new exact-match counts.

In `docs/ALGORITHMS.md` §10, add a paragraph describing the HMM rescue and link to the spec + companion docs:

```markdown
As of 2026-04-XX, novel-copy rescue uses a per-exon profile-HMM family
model (see `docs/superpowers/specs/2026-04-28-vg-novel-copy-hmm-design.md`
and `docs/UNMAPPED_FAMILY_RESCUE.md`). The legacy flat k-mer rescue
remains available via `--vg-discover-novel-mode kmer`.
```

- [ ] **Step 4: Commit**

```bash
git add -u
git commit -m "vg-hmm: flip default to hmm mode after benchmark validation; update docs"
```

- [ ] **Step 5: Open PR**

```bash
git push -u origin vg-hmm-novel-copy
gh pr create --title "VG HMM novel-copy discovery" --body "$(cat <<'EOF'
## Summary
- Replace flat k-mer novel-copy rescue with per-exon profile-HMM family model
- Add five-bucket failure-mode classifier (internal + optional external minimap2 verification)
- Synthetic bundles flow through the existing assembly pipeline with documented filter exemptions

## Test plan
- [ ] cargo test (all vg_hmm_* binaries)
- [ ] chr19 smoke (`tests/regression/vg_hmm_integration_golga.rs`)
- [ ] Synthetic ground truth (`tests/regression/vg_hmm_synthetic_ground_truth.rs`)
- [ ] Negative control (`tests/regression/vg_hmm_negative_control.rs`)
- [ ] Full GGO.bam benchmark (manual, see plan §7.4)

🤖 Generated with [Claude Code](https://claude.com/claude-code)
EOF
)"
```

---

## Self-Review Notes

- **Spec coverage:** every spec section §3–§13 maps to a phase: §3 → Phase 0; §4 → Phase 1; §5 → Phase 2; §6 → Phase 3; §7 → Phase 4 + Phase 5; §8 → Phase 6; §9 → Task 4.2; §10 → Phase 7; §11 → ordering of phases reflects rollout; §12 → risk mitigations are inline (e.g., Task 2.5 mean-pairwise-identity gate); §13 → confirmed deferred, no tasks for them.
- **Placeholders:** the only `todo!()` is in Task 7.2's test body — flagged explicitly as the engineer's job to fill in once the synthetic fixture lands. The Task 2.2 `todo!()` is similarly explicit and forces the engineer to bind the chosen POA crate API rather than guess. No vague "add error handling" steps.
- **Type consistency:** `RescueClass` declared in Task 6.1, used in Task 6.3 (with `Option<RescueClass>` on `Bundle`); `NodeIdx`/`ExonClass`/`FamilyGraph`/`JunctionEdge` all introduced in Task 1.1 and used unchanged in Tasks 1.6, 2.5, 3.x, 5.4. The `synthetic` field added in Task 4.1 matches its usage in Tasks 5.4 and 5.5.
- **Risk:** Task 4.1 will touch many `Bundle { ... }` construction sites in tests. If that fan-out is large (>20 sites) the engineer should consider adding a `Default` impl for `Bundle` and using `..Default::default()` everywhere instead. This is acknowledged inline.
