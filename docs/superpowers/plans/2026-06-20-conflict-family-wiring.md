# Conflict-Graph Family Wiring Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Replace POA-sequence-similarity as the authoritative family-membership criterion in `detect_and_assign` with the read-conflict graph, so domain-sharers (reads map to one locus only → no conflict edge) are never grouped with true multi-copy paralog families (reads genuinely cross-map with tied AS scores → conflict edges).

**Architecture:** After `collapse_loci` produces `reps`, run `build_read_placements` + `conflict_edges` + `conflict_families` on `bam_reads` (which include primary + secondary alignments). Convert the resulting `Vec<Vec<usize>>` families to `Vec<SplitFamily>` using existing `community_stats` + `classify` so `colocated_families` can consume them unchanged. POA detection (`detect_edges_reporting`) is kept running purely for diagnostics but its output is no longer used to drive `colocated_families`.

**Tech Stack:** Rust, existing `family_split::{community_stats, classify, SplitFamily, CommunityStats, FamilyClass, SplitParams}`, existing `read_conflict::{conflict_edges, conflict_families, ConflictParams}`, existing `denovo_pipeline::{build_read_placements, colocated_families}`.

## Global Constraints

- All new behaviour is gated through `DenovoConfig.conflict: ConflictParams` (default: `as_tie=0.9, min_reads=2`) — no hidden global state.
- `colocated_families` signature is **unchanged** — it still takes `&[SplitFamily]`.
- `detect_families` (no-BAM variant) is **unchanged** — it stays POA-only, that function does not have access to `bam_reads`.
- Every structural change is TDD: write the failing test, run to see red, implement, run to see green, commit.
- `cargo test --lib` must be 0 failures at every commit.
- Default `min_reads` for conflict is **2** (the real-data smoke revealed 1-read edges are noise; the AS-tie check already handles the principal filtering, but spurious single-read secondaries still slip through).

---

## File map

| File | Change |
|------|--------|
| `src/rustle/vg_family/denovo_pipeline.rs` | Add `conflict: ConflictParams` to `DenovoConfig`; add `conflict_to_split_families()`; rewire `detect_and_assign` to use conflict split; update test fixture |
| `src/rustle/vg_family/family_split.rs` | `pub(crate)` the `community_stats` function (currently `pub(super)` or private — verify and expose) |
| No new files | Everything fits in the existing modules |

---

## Task 1 — Expose `community_stats` and confirm `classify` is `pub`

`community_stats` computes `CommunityStats` (density, avg_core_recip, n_articulation) for an arbitrary member-set over a weighted edge list. We need it in `denovo_pipeline.rs`.

**Files:**
- Modify: `src/rustle/vg_family/family_split.rs`

**Interfaces:**
- Produces: `pub fn community_stats(members: &[usize], edges: &[(usize, usize, f64)]) -> CommunityStats` (crate-visible)
- Produces: `pub fn classify(n: usize, density: f64, p: &SplitParams) -> FamilyClass` (confirm already `pub`)

- [ ] **Step 1: Check current visibility of `community_stats` and `classify`**

```bash
grep -n "^fn community_stats\|^pub fn community_stats\|^pub(crate) fn community_stats\|^pub fn classify" \
  src/rustle/vg_family/family_split.rs
```

Expected: `community_stats` will be `fn` (private) or `pub(super)`. `classify` should already be `pub`.

- [ ] **Step 2: Make `community_stats` pub (no test needed — it is already tested indirectly via `decompose_families`)**

Find the line with `fn community_stats` and change it to `pub(crate) fn community_stats`. Example diff:

```rust
// before
fn community_stats(members: &[usize], edges: &[(usize, usize, f64)]) -> CommunityStats {

// after
pub(crate) fn community_stats(members: &[usize], edges: &[(usize, usize, f64)]) -> CommunityStats {
```

- [ ] **Step 3: Verify build is clean**

```bash
cargo build --lib 2>&1 | grep "^error"
```

Expected: no output (no errors).

- [ ] **Step 4: Commit**

```bash
git add src/rustle/vg_family/family_split.rs
git commit -m "refactor(family_split): pub(crate) community_stats for conflict-wiring reuse"
```

---

## Task 2 — Add `conflict: ConflictParams` to `DenovoConfig`

`DenovoConfig` is the single config struct threaded through the pipeline. Adding `conflict` here keeps all knobs in one place.

**Files:**
- Modify: `src/rustle/vg_family/denovo_pipeline.rs`

**Interfaces:**
- Consumes: `ConflictParams { as_tie: f64, min_reads: usize }` (already in `read_conflict.rs`, already imported via `use super::read_conflict::{..., ConflictParams, ...}`)
- Produces: `DenovoConfig.conflict: ConflictParams`

- [ ] **Step 1: Write the failing test**

Add this test to the `tests` module in `denovo_pipeline.rs`:

```rust
#[test]
fn denovoconfig_default_conflict_params_are_sane() {
    let cfg = DenovoConfig::default();
    // as_tie=0.9 means a 10% AS margin still counts as tied; min_reads=2 guards single spurious secondaries.
    assert!((cfg.conflict.as_tie - 0.9).abs() < 1e-9);
    assert_eq!(cfg.conflict.min_reads, 2);
}
```

- [ ] **Step 2: Run test to verify it fails**

```bash
cargo test --lib -- vg_family::denovo_pipeline::tests::denovoconfig_default_conflict_params_are_sane 2>&1 | tail -5
```

Expected: compile error or `FAILED` (field does not exist yet).

- [ ] **Step 3: Add `conflict` field to `DenovoConfig`**

In `denovo_pipeline.rs`, find `DenovoConfig` and its `Default` impl:

```rust
// before
pub struct DenovoConfig {
    pub pass1_min_reads: u32,
    pub gate: GateParams,
    pub detect: DetectParams,
    pub split: SplitParams,
}

impl Default for DenovoConfig {
    fn default() -> Self {
        DenovoConfig {
            pass1_min_reads: PASS1_MIN_READS,
            gate: GateParams::default(),
            detect: DetectParams::default(),
            split: SplitParams::default(),
        }
    }
}

// after
pub struct DenovoConfig {
    pub pass1_min_reads: u32,
    pub gate: GateParams,
    pub detect: DetectParams,
    pub split: SplitParams,
    pub conflict: ConflictParams,
}

impl Default for DenovoConfig {
    fn default() -> Self {
        DenovoConfig {
            pass1_min_reads: PASS1_MIN_READS,
            gate: GateParams::default(),
            detect: DetectParams::default(),
            split: SplitParams::default(),
            conflict: ConflictParams { as_tie: 0.9, min_reads: 2 },
        }
    }
}
```

- [ ] **Step 4: Run test to verify it passes**

```bash
cargo test --lib -- vg_family::denovo_pipeline::tests::denovoconfig_default_conflict_params_are_sane 2>&1 | tail -5
```

Expected: `test ... ok`.

- [ ] **Step 5: Full suite green**

```bash
cargo test --lib 2>&1 | tail -3
```

Expected: `N passed; 0 failed`.

- [ ] **Step 6: Commit**

```bash
git add src/rustle/vg_family/denovo_pipeline.rs
git commit -m "feat(denovo_pipeline): add conflict: ConflictParams to DenovoConfig (as_tie=0.9, min_reads=2)"
```

---

## Task 3 — Add `conflict_to_split_families` converter

This converts `conflict_families()`'s `Vec<Vec<usize>>` output to the `Vec<SplitFamily>` that `colocated_families` already expects, by running `community_stats` (structural diagnostics) and `classify` (Family vs Web) over each conflict component.

Edge weights in the conflict kernel are `usize` (read counts); `community_stats` takes `f64`. The conversion is `w as f64`.

**Files:**
- Modify: `src/rustle/vg_family/denovo_pipeline.rs`

**Interfaces:**
- Consumes: `community_stats` (Task 1), `classify`, `SplitParams` from `family_split`
- Consumes: `conflict_families: &[Vec<usize>]`, `c_edges: &[(usize, usize, usize)]` (integer-weight conflict edges), `n_loci: usize`, `p: &SplitParams`
- Produces: `fn conflict_to_split_families(families: &[Vec<usize>], c_edges: &[(usize, usize, usize)], p: &SplitParams) -> Vec<SplitFamily>`

The function lives beside `build_read_placements` in `denovo_pipeline.rs`.

- [ ] **Step 1: Add the import for `community_stats` and `classify` to `denovo_pipeline.rs`**

```rust
// Find this line:
use super::family_split::{decompose_families, FamilyClass, SplitFamily, SplitParams};

// Change to:
use super::family_split::{classify, community_stats, decompose_families, FamilyClass, SplitFamily, SplitParams};
```

- [ ] **Step 2: Write the failing test**

Add to the `tests` module:

```rust
#[test]
fn conflict_to_split_families_produces_one_family_per_component() {
    // Two disjoint conflict components: {0,1} linked by 5 reads, {2,3} by 3 reads.
    let families = vec![vec![0usize, 1], vec![2usize, 3]];
    let c_edges = vec![(0usize, 1usize, 5usize), (2, 3, 3)];
    let p = SplitParams::default();
    let split = conflict_to_split_families(&families, &c_edges, &p);
    assert_eq!(split.len(), 2);
    assert!(split.iter().all(|sf| sf.class == FamilyClass::Family));
    // members are sorted ascending within each family.
    assert_eq!(split[0].members, vec![0, 1]);
    assert_eq!(split[1].members, vec![2, 3]);
    // density of a 2-node clique (1 edge / 1 possible) = 1.0.
    for sf in &split {
        assert!((sf.stats.density - 1.0).abs() < 1e-9, "2-node clique density must be 1.0");
    }
}
```

- [ ] **Step 3: Run to verify it fails**

```bash
cargo test --lib -- vg_family::denovo_pipeline::tests::conflict_to_split_families_produces_one_family_per_component 2>&1 | tail -5
```

Expected: compile error (function not defined yet).

- [ ] **Step 4: Implement `conflict_to_split_families`**

Add this function to `denovo_pipeline.rs` (immediately after `build_read_placements`):

```rust
/// Convert the read-conflict kernel's component output to `SplitFamily` objects that `colocated_families`
/// consumes. Edge weights (read counts) are cast to `f64` for `community_stats`; class follows the same
/// size+density rule as the POA path so large, sparse conflict components can still be flagged as Webs.
fn conflict_to_split_families(
    families: &[Vec<usize>],
    c_edges: &[(usize, usize, usize)],
    p: &SplitParams,
) -> Vec<SplitFamily> {
    let float_edges: Vec<(usize, usize, f64)> =
        c_edges.iter().map(|&(a, b, w)| (a, b, w as f64)).collect();
    let mut out: Vec<SplitFamily> = families
        .iter()
        .map(|members| {
            let mut m = members.clone();
            m.sort_unstable();
            let stats = community_stats(&m, &float_edges);
            let class = classify(stats.n, stats.density, p);
            SplitFamily { members: m, stats, class }
        })
        .collect();
    // deterministic: size desc, then smallest member.
    out.sort_by(|a, b| {
        b.members.len().cmp(&a.members.len()).then_with(|| a.members[0].cmp(&b.members[0]))
    });
    out
}
```

- [ ] **Step 5: Run test to verify it passes**

```bash
cargo test --lib -- vg_family::denovo_pipeline::tests::conflict_to_split_families_produces_one_family_per_component 2>&1 | tail -5
```

Expected: `test ... ok`.

- [ ] **Step 6: Full suite green**

```bash
cargo test --lib 2>&1 | tail -3
```

- [ ] **Step 7: Commit**

```bash
git add src/rustle/vg_family/denovo_pipeline.rs src/rustle/vg_family/family_split.rs
git commit -m "feat(denovo_pipeline): add conflict_to_split_families — convert read-conflict components to SplitFamily"
```

---

## Task 4 — Fix the end-to-end test fixture: add a genuine cross-mapping BAM record

The existing test `detect_and_assign_resolves_multimapper_end_to_end` uses `two_paralogs_with_psvs()`. That fixture produces one `BamRead` at locus 0 (copyA's region, ref_start=0), with no secondary at locus 1. After we switch to conflict-based families, the conflict graph sees no cross-locus read → no family → the assignment step is skipped and the test assertion on `fas.len() == 1` fails.

Fix: add a SECONDARY `BamRead` record for the same read at locus 1's coordinates. Both records carry the same read name (mimicking a real multimapper: the aligner outputs primary at copyA and secondary at copyB with tied AS).

**Files:**
- Modify: `src/rustle/vg_family/denovo_pipeline.rs` (the `two_paralogs_with_psvs` helper and possibly the existing test assertions)

**Context:** In `two_paralogs_with_psvs`:
- `copyA` gene: exon1 [0,200) intron [200,220) exon2 [220,420) → rep at locus 0
- `copyB` gene: exon1 [1000,1200) intron [1200,1220) exon2 [1220,1420) → rep at locus 1
- The existing bam_read has `ref_start=0`, `cigar=[('M',200),('N',20),('M',200)]`, `seq=copyB_spliced`
- It crosses loci: the read's ref-span is [0,420), which is locus 0 only (locus 1 starts at 1000)

For the conflict graph to fire, we need the SAME read name to also appear at locus 1. A minimal secondary record:
- `name`: same as the primary read
- `ref_start = 1000`, `cigar = [('M', 200), ('N', 20), ('M', 200)]`
- `seq`: copyA_spliced (the secondary placement is the read at the paralog's locus)
- `as_score` ≈ the primary's score (tied), `mapq = 0` (secondary convention)

- [ ] **Step 1: Update `two_paralogs_with_psvs` to return the secondary as well**

In `two_paralogs_with_psvs`, the return `bam` vec currently has one entry. Add a second entry with `ref_start=1000`, same name, same (or 1 lower) `as_score`:

```rust
// find this block in two_paralogs_with_psvs:
let copyb_spliced = cat(&[&fb1, &core_b, &fb2]);
let ar = AlignedRead { ref_start: 0, cigar: vec![('M', 200), ('N', 20), ('M', 200)], seq: copyb_spliced };
let bam = vec![BamRead { chrom: "c1".into(), read: ar, mapq: 60, name: "readB".into(), as_score: 380 }];

// replace with:
let copyb_spliced = cat(&[&fb1, &core_b, &fb2]);
let copya_spliced = cat(&[&fa1, &core_a, &fa2]);
// primary: readB maps to copyA's genomic region (aligner primary placement).
let ar_primary = AlignedRead { ref_start: 0, cigar: vec![('M', 200), ('N', 20), ('M', 200)], seq: copyb_spliced };
// secondary: same read also maps (with tied AS) to copyB's genomic region.
let ar_secondary = AlignedRead { ref_start: 1000, cigar: vec![('M', 200), ('N', 20), ('M', 200)], seq: copya_spliced };
let bam = vec![
    BamRead { chrom: "c1".into(), read: ar_primary,   mapq: 60, name: "readB".into(), as_score: 380 },
    BamRead { chrom: "c1".into(), read: ar_secondary, mapq:  0, name: "readB".into(), as_score: 379 },
];
```

`as_score: 379` vs `380` — ratio 0.997 ≥ 0.9 (the default `as_tie`), so these are tied.

- [ ] **Step 2: Run the existing test to see it now fails (without the wiring change — conflict graph doesn't fire yet because we haven't switched the path)**

```bash
cargo test --lib -- vg_family::denovo_pipeline::tests::detect_and_assign_resolves_multimapper_end_to_end 2>&1 | tail -10
```

Expected at this point: the test still PASSES (we haven't wired conflict in yet — POA still drives families). If it passes, good; we are just preparing the fixture.

- [ ] **Step 3: Full suite green with updated fixture**

```bash
cargo test --lib 2>&1 | tail -3
```

Expected: same pass count as before (the secondary record doesn't break the POA-based path).

- [ ] **Step 4: Commit fixture change separately (so it's reviewable before the wiring)**

```bash
git add src/rustle/vg_family/denovo_pipeline.rs
git commit -m "test(denovo_pipeline): add secondary BAM record in two_paralogs_with_psvs fixture for conflict-wiring readiness"
```

---

## Task 5 — Wire conflict families into `detect_and_assign`

Replace `decompose_families(&edges, &cfg.split)` with `conflict_to_split_families(...)` as the source of `split` that feeds into `colocated_families`. Keep POA detection running for the fallback-edge audit log; stop using it for family membership.

**Files:**
- Modify: `src/rustle/vg_family/denovo_pipeline.rs`

**The current code block to replace** (lines ~290–344 after earlier edits):

```rust
// [diagnostic block]  { let placements = ...; let c_edges = ...; let c_fams = ...; eprintln!(...); }

let (edges, fallback_pairs) = detect_edges_reporting(&reps, &cfg.detect);
eprintln!("[detect_and_assign] {} homology edges ...", ...);
let fallback: Vec<FallbackEdge> = ...;
let split = decompose_families(&edges, &cfg.split);
```

**Target code:**

```rust
// Conflict-graph: the AUTHORITATIVE family criterion (reads genuinely ambiguous → family;
// domain-sharer exon reads map to one locus → no edge).
let placements = build_read_placements(bam_reads, &reps);
let c_edges = conflict_edges(reps.len(), &placements, &cfg.conflict);
let c_fams = conflict_families(reps.len(), &c_edges);
eprintln!(
    "[detect_and_assign] conflict-graph: {} edges -> {} families",
    c_edges.len(), c_fams.len(),
);
for (fi, fam) in c_fams.iter().enumerate() {
    let members: Vec<&str> = fam.iter().map(|&i| reps[i].tid.as_str()).collect();
    let coords: Vec<String> = fam.iter()
        .map(|&i| format!("{}:{}-{}", reps[i].chrom, reps[i].start, reps[i].end))
        .collect();
    let edge_weights: Vec<usize> = c_edges.iter()
        .filter(|&&(a, b, _)| fam.contains(&a) && fam.contains(&b))
        .map(|&(_, _, w)| w)
        .collect();
    eprintln!(
        "  conflict-fam{fi} n={} reads_linking={:?}: {} @ {}",
        fam.len(), edge_weights, members.join(","), coords.join(" | "),
    );
}
let split = conflict_to_split_families(&c_fams, &c_edges, &cfg.split);

// POA homology edges — kept for the large-seq fallback audit log; no longer drives family membership.
let (poa_edges, fallback_pairs) = detect_edges_reporting(&reps, &cfg.detect);
eprintln!(
    "[detect_and_assign] POA homology (diagnostic): {} edges ({} via large-seq fallback)",
    poa_edges.len(), fallback_pairs.len(),
);
let fallback: Vec<FallbackEdge> = fallback_pairs
    .iter()
    .map(|&(a, b)| FallbackEdge {
        chrom: reps[a].chrom.clone(),
        tid_a: reps[a].tid.clone(), start_a: reps[a].start, end_a: reps[a].end, len_a: reps[a].seq.len(),
        tid_b: reps[b].tid.clone(), start_b: reps[b].start, end_b: reps[b].end, len_b: reps[b].seq.len(),
    })
    .collect();
```

Note: `decompose_families` is no longer called in this path. The `detect_edges` import may become unused if `detect_edges` is not used elsewhere — verify and remove if so.

- [ ] **Step 1: Apply the code replacement described above**

Make the edit to `denovo_pipeline.rs`. Do NOT change `colocated_families` call below — it still receives `&split` and `win` unchanged.

- [ ] **Step 2: Run the existing end-to-end test**

```bash
cargo test --lib -- vg_family::denovo_pipeline::tests::detect_and_assign_resolves_multimapper_end_to_end 2>&1 | tail -10
```

Expected: PASS — the fixture now has a proper secondary record so the conflict graph fires, and assignment proceeds as before.

- [ ] **Step 3: Full suite**

```bash
cargo test --lib 2>&1 | tail -5
```

Expected: same pass count (591 + new tests from Tasks 2–4 = 594+).

- [ ] **Step 4: Clean up unused imports if any**

```bash
cargo build --lib 2>&1 | grep "unused import.*detect_edges\b"
```

If `detect_edges` (the non-reporting variant) is now unused in `denovo_pipeline.rs`, remove it from the import line:
```rust
// remove `detect_edges` from the use statement, keep `detect_edges_reporting`
use super::family_detect::{collapse_loci, detect_edges_reporting, DenovoTranscript, DetectParams};
```

Also remove `decompose_families` from the import if it is no longer called here:
```rust
// remove `decompose_families` from:
use super::family_split::{classify, community_stats, decompose_families, FamilyClass, SplitFamily, SplitParams};
// becomes:
use super::family_split::{classify, community_stats, FamilyClass, SplitFamily, SplitParams};
```

- [ ] **Step 5: Full suite green after cleanup**

```bash
cargo test --lib 2>&1 | tail -3
```

- [ ] **Step 6: Commit**

```bash
git add src/rustle/vg_family/denovo_pipeline.rs
git commit -m "feat(denovo_pipeline): wire conflict-graph families as authoritative family criterion in detect_and_assign"
```

---

## Task 6 — Smoke test on real data to confirm domain-sharer discrimination

Run the smoke test on three regions with known ground truth:

| Region | Expected outcome | Why |
|--------|-----------------|-----|
| RFPL: `NC_073243.2:104789647-104877901` | 0 conflict families | Copies on different chromosomes; within this region reads are not genuinely ambiguous |
| MAGEA: `NC_073247.2:161251228-164865959` | ≥2 conflict families with ≥38 reads each | Real tandem array; reads cross-map with tied AS |
| A known domain-sharer region (pick a gene with a shared exon — e.g. use the POA-only family from the earlier diagnostic that had 0 conflict edges) | 0 conflict families | Only the shared exon maps; no alternative placement → no edge |

- [ ] **Step 1: Run RFPL smoke**

```bash
RUSTLE_DENOVO_SMOKE_BAM=/home/juanfra/winloci_scratch/GGO.bam \
RUSTLE_DENOVO_SMOKE_FASTA=/home/juanfra/winloci_scratch/GGO.fasta \
RUSTLE_DENOVO_SMOKE_REGION=NC_073243.2:104789647-104877901 \
cargo test --release --lib -- --ignored smoke_detect_and_assign_real --nocapture 2>&1 \
  | grep -E "conflict-graph|co-located|AGGREGATE"
```

Expected:
```
[detect_and_assign] conflict-graph: 0 edges -> 0 families
co-located families with assignments: 0
```

- [ ] **Step 2: Run MAGEA smoke**

```bash
RUSTLE_DENOVO_SMOKE_BAM=/home/juanfra/winloci_scratch/GGO.bam \
RUSTLE_DENOVO_SMOKE_FASTA=/home/juanfra/winloci_scratch/GGO.fasta \
RUSTLE_DENOVO_SMOKE_REGION=NC_073247.2:161251228-164865959 \
cargo test --release --lib -- --ignored smoke_detect_and_assign_real --nocapture 2>&1 \
  | grep -E "conflict-fam|co-located|AGGREGATE"
```

Expected: 2–4 conflict families, co-located families ≥ 2, AGGREGATE unique-mapper agreement ≥ 95%.

- [ ] **Step 3: Confirm unique-mapper agreement holds**

The AGGREGATE unique-mapper agreement was 100% (52/52) before the wiring change. Verify it is still ≥ 95% after (conflict families may be slightly different from POA families, but assignment within families should be equally accurate).

- [ ] **Step 4: Commit smoke results as a note in the plan or as a commit message annotation**

No code change — just document in the commit message:
```bash
git commit --allow-empty -m "smoke(denovo_pipeline): conflict-graph families validated on MAGEA (≥95% unique-mapper agreement) and RFPL (0 edges) regions"
```

---

## Self-Review

**Spec coverage:**
- ✓ Domain-sharers excluded: single-placement reads make no edge (existing test + new test Task 2)
- ✓ Real copies found: multimapper with secondary record fires edge (Task 4 fixture + Task 5 wire)
- ✓ Config exposed: `DenovoConfig.conflict` (Task 2)
- ✓ No interface breakage: `colocated_families` signature unchanged (Task 3 produces `SplitFamily`)
- ✓ POA kept for diagnostics and fallback-edge audit (Task 5)
- ✓ Smoke validated on real data (Task 6)

**Placeholder scan:** None found — all steps have concrete code.

**Type consistency:**
- `conflict_to_split_families` takes `&[(usize, usize, usize)]` (integer weights from `conflict_edges`) — consistent with `read_conflict::conflict_edges` return type `Vec<(usize, usize, usize)>`.
- `community_stats` takes `&[(usize, usize, f64)]` — conversion `w as f64` inside `conflict_to_split_families` ✓
- `SplitFamily.members: Vec<usize>` — filled from `Vec<Vec<usize>>` returned by `conflict_families` ✓
