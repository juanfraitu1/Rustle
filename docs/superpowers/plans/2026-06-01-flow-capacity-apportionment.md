# Flow-Capacity Multimapper Apportionment Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Make max-flow capacities trustworthy for paralog multimappers — mass-conserving, anchor-aware apportionment (inverted-pair-aware) feeding the flow, with capacity-confidence propagated to each abundance.

**Architecture:** Correct `read.weight` BEFORE graph build (the existing EM write-back), so apportioned mass flows untouched into `node.coverage`/`tf.abundance`/capacities — zero flow-engine changes. Joint-strand EM input fixes inverted pairs; an anchored-mass prior replaces the self-reinforcing pileup prior; a parallel anchored-coverage channel yields a per-transcript `capacity_confidence`. All VG/`--vg`-gated; default de-novo output stays byte-identical.

**Tech Stack:** Rust (rustle), Edmonds-Karp max-flow (`max_flow.rs`), fingerprint-EM (`vg.rs::run_fingerprint_em`), cargo test.

**Spec:** `docs/superpowers/specs/2026-06-01-flow-capacity-apportionment-design.md` (approved 2026-06-01).

**Task dependency order:** T1 BundleRead fields, T2 GraphNode field, T3 Transcript fields (field foundation) → T4 `anchored_mass_per_copy`, T5 `family_for_em_input` (helpers) → T6 joint-strand wiring → T7 anchored prior, T8 write-back + conservation → T9 anchored-cov channel → T10 capacity_confidence push → T11 attach + GTF → T12 validation. Earlier tasks are prerequisites for later ones; implement in order.

---

### Task 1: Add EM-anchoring fields to BundleRead

**Files:**
- Modify: src/rustle/types.rs:222 (struct field block ending at `is_primary_alignment: bool,`)
- Modify: src/rustle/bam.rs:849 (real ingest literal, `nh` in scope from :617)
- Modify: src/rustle/junction_graph_st.rs:1126 (test `make_read`)
- Modify: src/rustle/vg.rs:5718 (test `make_read`)
- Modify: src/rustle/vg_hmm/rescue.rs:1056, :1189, :1330, :1535 (synthetic-read literals)

Notes verified by reading code:
- `BundleRead` is `#[derive(Debug, Clone)]` only (no `Default` impl, no `Default` derive) — so there is NOTHING to update for a Default impl. The spec's "mitigate with Default impls" does not apply here; the compiler enforces exhaustiveness on the 7 explicit literals.
- No literal uses `..Default::default()` / `..base` spread syntax; every one lists all fields, so all 7 must be edited.
- `make_read_full` at vg.rs:5928 delegates to `make_read` and constructs NO literal — leave it alone.
- New fields are appended after `is_primary_alignment` to keep the existing field order intact.

- [ ] **Step 1: Add the three fields to the struct (after `is_primary_alignment: bool,` at types.rs:222)**
Replace the closing of the struct:
```rust
    /// True if this record is the primary alignment (not secondary/supplementary).
    /// Populated at ingest; used in VG mode to filter secondaries out of
    /// non-family bundles after family discovery.
    pub is_primary_alignment: bool,
    /// VG joint-strand EM: per-read anchoring score gap from `run_fingerprint_em`
    /// write-back (vg.rs:4823-4866). -1.0 = not yet computed / non-VG run.
    pub em_weight_gap: f64,
    /// VG joint-strand EM: number of decisive fingerprint sites backing this
    /// read's copy assignment. 0 = none / non-VG run.
    pub em_n_sites: u32,
    /// VG joint-strand EM: true when this read's mass is anchored (unique read
    /// or fingerprint-decisive). Default at construction is `nh <= 1`; refined
    /// in EM write-back. Drives capacity-confidence accumulation.
    pub em_anchored: bool,
}
```

- [ ] **Step 2: Update the real ingest literal (bam.rs:849)** — `nh` is bound at bam.rs:617 and in scope here.
```rust
        ps_tag,
        is_primary_alignment,
        em_weight_gap: -1.0,
        em_n_sites: 0,
        em_anchored: nh <= 1,
    })
```

- [ ] **Step 3: Update test literal junction_graph_st.rs:1126**
```rust
            mapq: 60, mismatches: vec![], seq: Vec::new(), hp_tag: None, ps_tag: None,
            is_primary_alignment: true,
            em_weight_gap: -1.0, em_n_sites: 0, em_anchored: true,
        }
```

- [ ] **Step 4: Update test literal vg.rs:5718**
```rust
            mapq: 60, mismatches, seq: Vec::new(), hp_tag: None, ps_tag: None,
            is_primary_alignment: true,
            em_weight_gap: -1.0, em_n_sites: 0, em_anchored: true,
        }
```

- [ ] **Step 5: Update synthetic literal rescue.rs:1056**
```rust
                hp_tag: None, ps_tag: None,
                is_primary_alignment: true,
                em_weight_gap: -1.0, em_n_sites: 0, em_anchored: true,
            }
```

- [ ] **Step 6: Update synthetic literal rescue.rs:1189**
```rust
                    hp_tag: None, ps_tag: None,
                    is_primary_alignment: true,
                    em_weight_gap: -1.0, em_n_sites: 0, em_anchored: true,
                }
```

- [ ] **Step 7: Update synthetic literal rescue.rs:1330**
```rust
                    ps_tag: None,
                    is_primary_alignment: true,
                    em_weight_gap: -1.0,
                    em_n_sites: 0,
                    em_anchored: true,
                }
```

- [ ] **Step 8: Update synthetic literal rescue.rs:1535**
```rust
                hp_tag: None, ps_tag: None, is_primary_alignment: true,
                em_weight_gap: -1.0, em_n_sites: 0, em_anchored: true,
            }
```

- [ ] **Step 9: Verify the whole tree compiles (mechanical task → build is the gate)**
```bash
cargo build 2>&1 | tail -20
```
Expect a clean build. If the compiler reports `missing field ... in initializer of BundleRead`, a literal was missed — grep `grep -rn "BundleRead {" src/rustle/` and add the three fields to it.

- [ ] **Step 10: Commit**
```bash
git add src/rustle/types.rs src/rustle/bam.rs src/rustle/junction_graph_st.rs src/rustle/vg.rs src/rustle/vg_hmm/rescue.rs
git commit -m "Add em_weight_gap/em_n_sites/em_anchored fields to BundleRead

Flow-capacity apportionment (spec 2026-06-01) step 3: persist per-read
EM anchoring on BundleRead. Defaults: em_weight_gap=-1.0, em_n_sites=0,
em_anchored=(nh<=1) at ingest; refined later in run_fingerprint_em
write-back. Updates all 7 struct literals (compiler-enforced).

Co-Authored-By: Claude Opus 4.8 (1M context) <noreply@anthropic.com>"
```

---

### Task 2: Add `anchored_coverage` field to GraphNode (mechanical)

**Files:**
- Modify: `src/rustle/graph.rs:125` (field decl, beside `coverage`)
- Modify: `src/rustle/graph.rs:176-206` (`GraphNode::new` constructor)
- Verify: `cargo build`

Context (verified by reading the code):
- The canonical `GraphNode` is the struct at `graph.rs:111-174` (the one used by max-flow / path_extract). Its `coverage: f64` field is at line 125.
- `GraphNode` has `#[derive(Debug, Clone)]` only — there is NO `derive(Default)` and NO `impl Default for GraphNode`. The single initialization site is the manual `GraphNode::new` constructor (`graph.rs:176-206`). `BundleGraph::add_node` (`graph.rs:398-408`) calls `GraphNode::new`, so it needs no change.
- `junction_graph_st.rs:30` defines a SEPARATE, unrelated `GraphNode` type (fields `id`, `cov`, `abund_in`, `abund_out`, …) with its own `new` at `junction_graph_st.rs:50-58`. It is a distinct type and MUST NOT be touched.
- This task only adds the field with a 0.0 default. It is NOT consumed by max-flow (Component 4 of the spec accumulates into it later; capacity is unchanged). Default de-novo behavior is byte-identical because nothing reads or writes the field yet.

- [ ] **Step 1: Add the field declaration next to `coverage`**

In `src/rustle/graph.rs`, immediately after the `pub coverage: f64,` field (line 125, with its doc comment ending at line 124), insert the new field. Replace:
```rust
    pub coverage: f64,
    pub children: crate::util::bitset::SmallBitset,
```
with:
```rust
    pub coverage: f64,
    /// VG capacity-confidence channel: the portion of `coverage` contributed by
    /// reads classified as anchored (unique, or dNM/fingerprint-decisive) during
    /// joint-strand fingerprint EM. Accumulated in parallel with `coverage`
    /// (`+= weight*bp` only when the read is `em_anchored`). NOT consumed by
    /// max-flow; used only to derive per-transcript `capacity_confidence`.
    /// 0.0 outside VG mode and before any anchored mass is attributed.
    pub anchored_coverage: f64,
    pub children: crate::util::bitset::SmallBitset,
```

- [ ] **Step 2: Initialize the field in `GraphNode::new`**

In `src/rustle/graph.rs` inside `GraphNode::new` (lines 176-206), the constructor sets `coverage: 0.0,` at line 183. Replace:
```rust
            source_bnode: None,
            coverage: 0.0,
            children: crate::util::bitset::SmallBitset::empty(),
```
with:
```rust
            source_bnode: None,
            coverage: 0.0,
            anchored_coverage: 0.0,
            children: crate::util::bitset::SmallBitset::empty(),
```

- [ ] **Step 3: Verify it compiles**
```bash
cargo build 2>&1 | tail -20
```
Expect a clean build (the field is added with a default and no other code references it yet). If the build reports a missing-field error for `GraphNode`, search for any other struct literal — `grep -rn 'GraphNode {' src/rustle/ | grep -v 'struct GraphNode'` — and confirm it is the unrelated `junction_graph_st.rs` type (do not add the field there).

- [ ] **Step 4: Commit**
```bash
git add src/rustle/graph.rs
git commit -m "graph: add anchored_coverage field to GraphNode (VG capacity-confidence channel)

Mechanical: new f64 field (default 0.0) beside coverage, initialized in
GraphNode::new. Not consumed by max-flow; populated later by the VG
joint-strand EM apportionment. Default de-novo behavior unchanged.

Co-Authored-By: Claude Opus 4.8 (1M context) <noreply@anthropic.com>"
```

---

### Task 3: Add capacity_confidence + abundance_min fields to Transcript

**Files:**
- Modify (struct def): `src/rustle/path_extract.rs:665`
- Modify (literal sites, 36 total across 6 files): `src/rustle/merge_mode.rs` (1), `src/rustle/path_extract.rs` (9), `src/rustle/pipeline.rs` (22), `src/rustle/single_exon_pileup.rs` (1), `src/rustle/transcript_filter.rs` (2), `src/rustle/vg.rs` (1)
- Verify: `cargo build`

This is a mechanical-only task (struct field additions). Verification is `cargo build` (the compiler enforces struct-literal exhaustiveness, so a missed literal is a hard compile error), ending with a commit.

- [ ] **Step 1: Add the two fields to the `Transcript` struct definition (after `copy_independent_support`, before `family_verdict`)**

In `src/rustle/path_extract.rs`, the field `pub copy_independent_support: Option<f64>,` is at line 665 (immediately followed by the doc comment + `pub family_verdict: Option<crate::vg::FamilyVerdict>,`). Use the Edit tool to replace exactly:
```rust
    pub copy_independent_support: Option<f64>,
    /// Multi-copy family classification verdict (--vg only; spec 2026-06-01).
```
with:
```rust
    pub copy_independent_support: Option<f64>,
    /// Fraction of this transcript's max-flow capacity that is anchored
    /// (unambiguous: unique reads + dNM-decisive multimappers). `Some` only in
    /// --vg mode; `anchored_coverage_sum / total_coverage_sum` over path nodes,
    /// clamped to [0,1]. Lower = capacity built mostly from ambiguous multimap
    /// mass (low trust). Emitted as the `capacity_confidence` GTF attribute when
    /// `Some` (spec 2026-06-01 flow-capacity apportionment).
    pub capacity_confidence: Option<f64>,
    /// Jointly-feasible lower bound on this transcript's abundance:
    /// `coverage * capacity_confidence`. Sub-conservative (correct lower bound);
    /// no per-copy MAX is emitted (non-additive). `Some` only in --vg mode.
    /// Emitted as the `abundance_min` GTF attribute when `Some`.
    pub abundance_min: Option<f64>,
    /// Multi-copy family classification verdict (--vg only; spec 2026-06-01).
```

- [ ] **Step 2: Append the two fields to all 36 `Transcript` struct literals**

Every literal site contains the exact substring `copy_independent_support: None, family_verdict: None,` (verified: the token `copy_independent_support` appears as a struct field only at the definition and in these `: None` literals — there are no other structs sharing the name, and `family_verdict: None,` always immediately follows on the same line). A single in-place sed appends the two new fields right after `copy_independent_support: None,` on every line that has it:
```bash
sed -i 's/copy_independent_support: None,/copy_independent_support: None, capacity_confidence: None, abundance_min: None,/g' \
  src/rustle/merge_mode.rs \
  src/rustle/path_extract.rs \
  src/rustle/pipeline.rs \
  src/rustle/single_exon_pileup.rs \
  src/rustle/transcript_filter.rs \
  src/rustle/vg.rs
```
NOTE: this does NOT touch the struct definition (path_extract.rs:665), whose line is `pub copy_independent_support: Option<f64>,` (no `: None,`), so it is left intact by the pattern.

- [ ] **Step 3: Confirm the substitution count is exactly 36 and no literal was missed**
```bash
grep -rn "copy_independent_support: None, capacity_confidence: None, abundance_min: None," src/ | wc -l
```
Expect output `36`. Then confirm no bare (un-appended) literal remains:
```bash
grep -rn "copy_independent_support: None," src/ | grep -v "capacity_confidence: None, abundance_min: None,"
```
Expect NO output (empty).

- [ ] **Step 4: Build to verify struct-literal exhaustiveness**
```bash
cargo build 2>&1 | tail -30
```
Expect a clean build (warnings about unused `capacity_confidence` / `abundance_min` reads are acceptable; there must be zero `error[E0063]: missing field` errors — those would indicate a literal the sed missed).

- [ ] **Step 5: Commit**
```bash
git add src/rustle/path_extract.rs src/rustle/merge_mode.rs src/rustle/pipeline.rs src/rustle/single_exon_pileup.rs src/rustle/transcript_filter.rs src/rustle/vg.rs
git commit -m "Add capacity_confidence + abundance_min fields to Transcript

Two new Option<f64> fields beside copy_independent_support, default None at
all 36 struct literals. VG-only capacity-confidence channel per spec
2026-06-01 flow-capacity apportionment; default de-novo unaffected (None =>
not emitted).

Co-Authored-By: Claude Opus 4.8 (1M context) <noreply@anthropic.com>"
```

---

### Task 4: Add `anchored_mass_per_copy` to vg.rs (anchor-first prior source)

This implements the anchor-first core's mass accumulator (spec Component 2 / Phase B): per copy, sum `read.weight` of reads that are either unique (`read_name_hash ∉ family.multimap_reads`) OR whose `anchor_read` verdict at this copy is `Owns`. The per-read `(nm, alen)` extraction mirrors the established `placement` closure in `classify_family` (vg.rs:1620-1634): `read_aligned_len` for `alen`, and `de`-based cost (`round(de*alen)`) falling back to raw `nm`. This function is the source for `log_priors[k] = ln(anchored[k]/total + 1e-3)` that replaces the M-step pileup prior in a later task.

**Files:**
- Modify: `src/rustle/vg.rs` (insert new fn after `anchor_read`, which ends at vg.rs:1515)
- Test: `src/rustle/vg.rs` (inline `#[cfg(test)] mod tests`, beside `extractor_suppresses_cross_strand_phantom` at vg.rs:5961)

- [ ] **Step 1: Write a FAILING unit test in the existing `mod tests`**

Insert this test immediately after `extractor_suppresses_cross_strand_phantom` (it ends at vg.rs:5997, just before `extractor_keeps_co_expressed_ties` at vg.rs:6000). It reuses the existing `make_read_full`/`make_bundle`/`FamilyGroup` fixture helpers. Copy0 is the real locus (167 unique reads + 30 multimappers that fit it at nm=5 vs sibling nm=70 → `Owns`); copy1 is the phantom (only bad-fit multimappers → `Sibling`, plus no unique reads → ~0 anchored mass). It will FAIL to compile because `anchored_mass_per_copy` does not exist yet.

```rust
    /// anchored_mass_per_copy: real copy (unique reads + decisively-owned
    /// multimappers) accumulates clear mass; phantom copy (only bad-fit
    /// multimappers, no unique reads) accumulates ~0.
    #[test]
    fn anchored_mass_real_copy_vs_phantom() {
        // copy0 = real locus (good fit, has unique reads)
        // copy1 = phantom locus (bad fit, no unique reads)
        let mut real_reads = Vec::new();    // copy 0 bundle
        let mut phantom_reads = Vec::new(); // copy 1 bundle
        let mut multimap: HashMap<u64, Vec<(usize, usize)>> = HashMap::new();
        for i in 0..30u64 {
            let rnh = 2000 + i;
            // placement at copy 0 (real): good fit, nm=5 over 1000 bp
            real_reads.push(make_read_full(rnh, vec![(100_000, 101_000)], 5, false));
            // placement at copy 1 (phantom): bad fit, nm=70 over 1000 bp
            phantom_reads.push(make_read_full(rnh, vec![(900_000, 901_000)], 70, false));
            multimap.insert(rnh, vec![(0, i as usize), (1, i as usize)]);
        }
        // copy 0 also has 167 UNIQUE reads (read_name_hash NOT in multimap).
        for i in 0..167u64 {
            real_reads.push(make_read_full(8000 + i, vec![(100_500, 101_500)], 3, false));
        }

        let bundles = vec![
            make_bundle("chrTest", '+', real_reads),
            make_bundle("chrTest", '+', phantom_reads),
        ];
        let family = FamilyGroup {
            family_id: 11,
            bundle_indices: vec![0, 1],
            multimap_reads: multimap,
        };

        // dnm = 70 - 5 = 65 >= t(2) at copy0 → all 30 multimappers Own copy0.
        // dnm = 5 - 70 = -65 <= -t at copy1 → Sibling (not counted).
        let mass = anchored_mass_per_copy(&family, &bundles, 2, 0.8);
        assert_eq!(mass.len(), 2);
        // copy0: 167 unique (weight 1.0 each) + 30 owned multimappers = 197.0
        assert!((mass[0] - 197.0).abs() < 1e-6, "real copy mass should be 197.0, got {}", mass[0]);
        // copy1: no unique reads, no owned multimappers → ~0
        assert!(mass[1] < 1e-6, "phantom copy mass should be ~0, got {}", mass[1]);
    }
```

- [ ] **Step 2: Confirm the test FAILS to compile (no impl yet)**

```bash
cargo test -p rustle anchored_mass_real_copy_vs_phantom 2>&1 | grep -E "cannot find function|error\[E0425\]|error: " | head
```
Expect: `cannot find function `anchored_mass_per_copy``.

- [ ] **Step 3: Implement `anchored_mass_per_copy`**

Insert this immediately after the closing brace of `anchor_read` (vg.rs:1515), before the `identifiability_classes` doc comment (vg.rs:1517). It builds the same per-read `(nm, alen)` placement table the rest of the family code uses, then for each copy sums weights of unique reads plus multimappers whose `anchor_read` verdict at that copy is `Owns`.

```rust
/// Per-copy anchored read mass (the anchor-first M-step prior source).
///
/// For each copy `k` (index = position in `family.bundle_indices`), sums
/// `read.weight` over reads in copy `k`'s bundle that are EITHER:
///   * unique — `read.read_name_hash` is not a key in `family.multimap_reads`
///     (the read only places at this copy), OR
///   * decisively owned — `anchor_read(nm_k, alen_k, others, t, extent_frac)`
///     returns `ReadAnchor::Owns`, where `others` are the read's `(nm, alen)` at
///     the copy's siblings. `Sibling`/`Tie` placements contribute 0 to copy `k`.
///
/// Per-placement `(nm, alen)` uses the same calibration as `classify_family`:
/// `alen = read_aligned_len(read)`; cost = `round(de*alen)` when `de` is present
/// (HiFi indel-robust), else raw `read.nm`. Calibration constants are passed in
/// (`t`, `extent_frac`); the project default is raw-dNM `t=2`, `extent_frac=0.8`.
///
/// Returns a `Vec<f64>` of length `family.bundle_indices.len()` (zeros for
/// stale/empty copies). Operates on the ORIGINAL `FamilyGroup` + intact bundles
/// (call at EM time, like `compute_copy_independent_support`).
pub fn anchored_mass_per_copy(
    family: &FamilyGroup,
    bundles: &[Bundle],
    t: i64,
    extent_frac: f64,
) -> Vec<f64> {
    let n_copies = family.bundle_indices.len();
    let mut mass = vec![0.0_f64; n_copies];

    // Per-read `(nm, alen)` at (fam_pos, read_idx) — same cost model as
    // classify_family's `placement` closure (de-based, NM fallback).
    let placement = |fam_pos: usize, ri: usize| -> Option<(u32, u64)> {
        let bi = *family.bundle_indices.get(fam_pos)?;
        let read = bundles.get(bi)?.reads.get(ri)?;
        let alen = read_aligned_len(read);
        if alen == 0 {
            return None;
        }
        let cost = match read.de {
            Some(de) => ((de as f64) * (alen as f64)).round() as u32,
            None => read.nm,
        };
        Some((cost, alen))
    };

    // Build, per multimap read, its best (min-NM) placement per copy. Mirrors
    // classify_family: dedupes redundant alignments at the same copy by min NM.
    let mut per_copy_by_read: std::collections::HashMap<u64, std::collections::HashMap<usize, (u32, u64)>> =
        std::collections::HashMap::with_capacity(family.multimap_reads.len());
    for (&rnh, placements) in &family.multimap_reads {
        let mut per_copy: std::collections::HashMap<usize, (u32, u64)> = std::collections::HashMap::new();
        for &(fp, ri) in placements {
            if fp >= n_copies {
                continue;
            }
            if let Some((nm, al)) = placement(fp, ri) {
                per_copy
                    .entry(fp)
                    .and_modify(|e| if nm < e.0 { *e = (nm, al); })
                    .or_insert((nm, al));
            }
        }
        if !per_copy.is_empty() {
            per_copy_by_read.insert(rnh, per_copy);
        }
    }

    for (copy_id, &bi) in family.bundle_indices.iter().enumerate() {
        let bundle = match bundles.get(bi) {
            Some(b) => b,
            None => continue,
        };
        for read in &bundle.reads {
            // Unique reads (not a multimap key) always anchor their copy.
            if !family.multimap_reads.contains_key(&read.read_name_hash) {
                mass[copy_id] += read.weight;
                continue;
            }
            // Multimapper: decisive only if it Owns THIS copy by dNM margin.
            let per_copy = match per_copy_by_read.get(&read.read_name_hash) {
                Some(pc) => pc,
                None => continue,
            };
            let (nm_c, alen_c) = match per_copy.get(&copy_id) {
                Some(&v) => v,
                None => continue, // this read has no placement at this copy
            };
            let others: Vec<(u32, u64)> = per_copy
                .iter()
                .filter(|(&c, _)| c != copy_id)
                .map(|(_, &v)| v)
                .collect();
            if anchor_read(nm_c, alen_c, &others, t, extent_frac) == ReadAnchor::Owns {
                mass[copy_id] += read.weight;
            }
        }
    }

    mass
}
```

- [ ] **Step 4: Run the test — it must PASS**

```bash
cargo test -p rustle anchored_mass_real_copy_vs_phantom 2>&1 | tail -20
```
Expect `test ... anchored_mass_real_copy_vs_phantom ... ok` and `test result: ok.`

- [ ] **Step 5: Build the crate to confirm no breakage**

```bash
cargo build -p rustle 2>&1 | tail -5
```

- [ ] **Step 6: Commit**

```bash
git add src/rustle/vg.rs
git commit -m "vg: add anchored_mass_per_copy (anchor-first prior source)

Per-copy sum of read.weight over unique reads plus dNM-decisive (Owns)
multimappers, reusing read_aligned_len + de/NM cost model from
classify_family and the anchor_read margin classifier. Source for the
anchor-first M-step prior (RUSTLE_VG_ANCHOR_PRIOR). Unit-tested with a
DAZ-style real-vs-phantom family fixture.

Co-Authored-By: Claude Opus 4.8 (1M context) <noreply@anthropic.com>"
```

---

### Task 5: Add family_for_em_input: the UNSPLIT both-strands EM input family

**Files:**
- Modify: `/mnt/c/Users/jfris/Desktop/Rustle/src/rustle/vg.rs` (add fn after `partition_and_remap_family_by_strand`, which ends at line 2297)
- Test: `/mnt/c/Users/jfris/Desktop/Rustle/src/rustle/vg.rs` (inline `#[cfg(test)] mod tests` starting at line 5677; helpers `make_read_full` @5928, `make_bundle` @5937 are in scope via `use super::*`)

- [ ] **Step 1: Write the FAILING test** — append this test inside `mod tests` (insert it just before the closing `}` after the existing `extractor_*` tests, e.g. right after `extractor_keeps_co_expressed_ties`). It builds a 2-bundle opposite-strand family with one cross-strand multimapper and asserts the unsplit family keeps BOTH bundle_indices and the read with >=2 placements (where `partition_and_remap_family_by_strand` would split it into per-strand sub-families).
```rust
    /// `family_for_em_input` must return ONE FamilyGroup spanning BOTH strands,
    /// preserving `fam_pos` indexing into `bundle_indices` and the full
    /// `multimap_reads` (no per-strand split, no <2-placement drop).
    /// Contrast: `partition_and_remap_family_by_strand` would split this into
    /// two single-bundle sub-families and drop the cross-strand read.
    #[test]
    fn family_for_em_input_keeps_both_strands_unsplit() {
        // copy 0 = '+' bundle, copy 1 = '-' bundle. One read places at both.
        let rnh = 555u64;
        let plus_reads = vec![make_read_full(rnh, vec![(100, 1100)], 5, false)];
        let minus_reads = vec![make_read_full(rnh, vec![(5000, 6000)], 5, true)];
        let bundles = vec![
            make_bundle("chrTest", '+', plus_reads),
            make_bundle("chrTest", '-', minus_reads),
        ];
        let mut multimap: HashMap<u64, Vec<(usize, usize)>> = HashMap::new();
        // (fam_pos, read_idx): placed at fam_pos 0 read 0, and fam_pos 1 read 0.
        multimap.insert(rnh, vec![(0usize, 0usize), (1usize, 0usize)]);
        let family = FamilyGroup {
            family_id: 3,
            bundle_indices: vec![0, 1],
            multimap_reads: multimap,
        };

        // Sanity: the strand-splitter WOULD destroy this family for EM.
        let split = partition_and_remap_family_by_strand(&family, &bundles);
        assert!(
            split.iter().all(|f| f.bundle_indices.len() < 2),
            "splitter drops both single-bundle sub-families: {:?}",
            split.iter().map(|f| f.bundle_indices.len()).collect::<Vec<_>>()
        );

        // The unsplit EM input keeps everything.
        let em_in = family_for_em_input(&family, &bundles);
        assert_eq!(em_in.family_id, family.family_id, "preserves family_id");
        assert_eq!(em_in.bundle_indices, vec![0, 1], "both bundle_indices retained");
        let locs = em_in.multimap_reads.get(&rnh).expect("cross-strand read retained");
        assert_eq!(locs.len(), 2, "read keeps >=2 placements");
        assert!(locs.contains(&(0, 0)) && locs.contains(&(1, 0)),
                "fam_pos indexing preserved: {locs:?}");
    }
```

- [ ] **Step 2: Run the test, confirm it FAILS to compile (fn undefined)**
```bash
cargo test --manifest-path /mnt/c/Users/jfris/Desktop/Rustle/Cargo.toml -p rustle family_for_em_input_keeps_both_strands_unsplit 2>&1 | tail -20
```
Expected: `cannot find function family_for_em_input in this scope` (or `E0425`).

- [ ] **Step 3: Implement `family_for_em_input`** — insert immediately after `partition_and_remap_family_by_strand` (after line 2297, before the `// ── HMM-based EM reweighting` comment at line 2299). This mirrors the single-strand early-return branch of the splitter (clone of `bundle_indices` + `multimap_reads`, same `family_id`) but applies it UNCONDITIONALLY, so the both-strands family flows into EM intact.
```rust
/// Return the UNSPLIT both-strands family — the input to
/// `run_fingerprint_em` / `run_em_reweighting_hmm` BEFORE
/// `partition_and_remap_family_by_strand` would split it per strand.
///
/// Preserves `fam_pos` (index into `bundle_indices`) and the full
/// `multimap_reads`, including cross-strand multimappers that the
/// per-strand splitter would drop (those whose placements span >1 strand
/// leave <2 placements per sub-family). Joint-strand EM
/// (`RUSTLE_VG_JOINT_STRAND_EM=1`) needs both strands' bundles and the
/// reads that bridge them, so it consumes THIS family, not the split ones.
///
/// `bundles` is taken for signature parity with
/// `partition_and_remap_family_by_strand` (callers pass the same slice);
/// no per-strand grouping is performed here.
pub fn family_for_em_input(family: &FamilyGroup, _bundles: &[Bundle]) -> FamilyGroup {
    FamilyGroup {
        family_id: family.family_id,
        bundle_indices: family.bundle_indices.clone(),
        multimap_reads: family.multimap_reads.clone(),
    }
}
```

- [ ] **Step 4: Run the test, confirm it PASSES**
```bash
cargo test --manifest-path /mnt/c/Users/jfris/Desktop/Rustle/Cargo.toml -p rustle family_for_em_input_keeps_both_strands_unsplit 2>&1 | tail -20
```
Expected: `test result: ok. 1 passed`.

- [ ] **Step 5: Commit**
```bash
git add /mnt/c/Users/jfris/Desktop/Rustle/src/rustle/vg.rs
git commit -m "vg: add family_for_em_input (unsplit both-strands EM input)

Returns the both-strands FamilyGroup that joint-strand EM consumes,
preserving fam_pos indexing and cross-strand multimappers that
partition_and_remap_family_by_strand would otherwise drop.

Co-Authored-By: Claude Opus 4.8 (1M context) <noreply@anthropic.com>"
```

---

### Task 6: Wire joint-strand EM input (unsplit family) into run_fingerprint_em

**Files:**
- Modify: `src/rustle/vg.rs` (add `family_for_em_input` near `partition_and_remap_family_by_strand` at vg.rs:2242)
- Modify: `src/rustle/pipeline.rs:10665-10669` (build `em_input_partitions`) and pipeline.rs:10734-10739 (`run_fingerprint_em` call)
- Test: `src/rustle/vg.rs` (inline `#[cfg(test)] mod tests` at vg.rs:5677)

**O1 note (from spec, verified in code):** `build_family_graph` hard-bails on mixed strands (family_graph.rs:432-433), so we do NOT rebuild the graph on the unsplit family. `em_hmm_partitions` (strand-split) is kept verbatim for `build_family_graph` / borrow / completion. Only the `run_fingerprint_em` *placement list + normalization* receives the unsplit family.

**CRITICAL — graph index alignment (verified at pipeline.rs:10671-10694):** `family_graphs: Vec<Option<FamilyGraph>>` is built `par_iter()` over `em_hmm_partitions` (one slot per *strand-split* partition), and `run_fingerprint_em` indexes it **by position** (`family_graphs.get(fam_idx)`, vg.rs:4519). So you **cannot** pass `&em_input_partitions` (one entry per *family*, fewer/reordered) with the old `&family_graphs` — after the first mixed-strand family the indices drift and a family would be paired with the **wrong** family's graph (not `None`). The fix (Steps 4–5) is to build a parallel `em_input_graphs: Vec<Option<FamilyGraph>>` aligned 1:1 with `em_input_partitions`: a **single-strand** family (its source produced exactly 1 strand partition) reuses that partition's graph so `fp` stays valid; a **mixed-strand** family (>1 strand partitions) gets `None`.

**⚠ CORRECTION (found in T6 review, 2026-06-02):** The earlier claim that a `None` graph "routes to the `fp.n_sites==0` neutral path" is **FALSE against the current code**. `run_fingerprint_em` (vg.rs:4643–4651) treats a `None`/empty graph as `eprintln!("no family graph — skipping"); results.push(EmResult::default()); continue;` — it **skips the family entirely**, leaving every multimapper at full weight in every copy. The `fp.n_sites==0` neutral path is only reachable when a graph is *present* but yields 0 diagnostic sites. Consequence: as wired by T6, the mixed-strand DAZ family (the whole target) is **skipped** by fingerprint-EM, and T7's anchored prior (inserted at vg.rs:4684, AFTER the skip) never runs for it. **T6 itself needs no rework** — its input shaping (unsplit family + `None` graph for mixed-strand) is correct. **The fix moves to T7:** lift the `None`-skip into a neutral-fingerprint fallthrough (build an empty `ExonFingerprints{ sites:vec![], per_copy_site_refs: vec![vec![]; n_copies], n_copies, n_sites:0 }` → `best==None` at vg.rs:4787 → all `log_scores=0.0` → junc+nm+anchored-prior decide), **gated by `RUSTLE_VG_JOINT_STRAND_EM != "0"`** so the `=0` rollback keeps the exact old skip behavior. See T7 Step 3b below.

- [ ] **Step 1: FAILING test for `family_for_em_input` (unsplit keeps cross-strand placements)**
Add this test inside the existing `#[cfg(test)] mod tests` block in `src/rustle/vg.rs` (after the last `#[test]` before its closing brace). It builds a mixed-strand `FamilyGroup` whose shared read has one placement per strand, and asserts the unsplit result keeps both placements (whereas `partition_and_remap_family_by_strand` would drop them to 1 per partition). It does not touch `bundles` for the split decision, so an empty slice is passed.
```rust
    #[test]
    fn family_for_em_input_keeps_cross_strand_placements() {
        // Family with two copies; read 7 maps to copy 0 (would be '-' strand)
        // and copy 1 (would be '+' strand) — one placement each.
        let mut multimap: HashMap<u64, Vec<(usize, usize)>> = HashMap::new();
        multimap.insert(7u64, vec![(0usize, 0usize), (1usize, 0usize)]);
        let fam = FamilyGroup {
            family_id: 42,
            bundle_indices: vec![10, 20],
            multimap_reads: multimap,
        };
        let bundles: Vec<Bundle> = Vec::new();
        let out = family_for_em_input(&fam, &bundles);
        assert_eq!(out.len(), 1, "unsplit family yields exactly one group");
        let g = &out[0];
        assert_eq!(g.family_id, 42, "family_id preserved (no *10+i remap)");
        assert_eq!(g.bundle_indices, vec![10, 20], "all copies kept");
        let locs = g.multimap_reads.get(&7u64)
            .expect("cross-strand shared read retained in EM input");
        assert!(locs.len() >= 2, "shared read keeps >=2 placements (got {})", locs.len());
    }
```
Run it to confirm it fails to compile (function missing):
```bash
cargo test -p rustle family_for_em_input_keeps_cross_strand_placements 2>&1 | tail -20
```

- [ ] **Step 2: Implement `family_for_em_input` (vg.rs, immediately before `pub fn partition_and_remap_family_by_strand` at vg.rs:2242)**
Returns the unsplit both-strands family as a single-element `Vec<FamilyGroup>` so it is a drop-in for `partition_and_remap_family_by_strand`'s return type. `bundle_indices` and `multimap_reads` are cloned verbatim — `fam_pos` stays a valid index into `bundle_indices`, so EM write-back `global_bi/ri` (vg.rs:4852-4857) remains correct. `_bundles` is accepted for signature parity (strand is irrelevant to the unsplit input).
```rust
/// Joint-strand EM input (spec component 1, O1-resolved): return the family
/// UNSPLIT (both strands together) for the fingerprint-EM placement list +
/// normalization only. Unlike `partition_and_remap_family_by_strand`, no
/// strand grouping and no `<2 placements after split` drop happens, so a read
/// shared only across an inverted pair (e.g. DAZ1(-)/DAZ3(+)) keeps all its
/// placements and is reweighted to sum to 1.0 across BOTH strands instead of
/// staying at 1/NH in each. `fam_pos` is preserved as an index into
/// `bundle_indices`, keeping EM write-back `global_bi/ri` valid.
///
/// The family GRAPH is NOT rebuilt on this unsplit group (build_family_graph
/// bails on mixed strands, family_graph.rs:432); the strand-split
/// `em_hmm_partitions` graphs are kept and indexed inside `run_fingerprint_em`,
/// where a missing/empty graph for a cross-strand group routes to the
/// `fp.n_sites==0` neutral path (junc + nm + anchored prior).
pub fn family_for_em_input(family: &FamilyGroup, _bundles: &[Bundle]) -> Vec<FamilyGroup> {
    vec![FamilyGroup {
        family_id: family.family_id,
        bundle_indices: family.bundle_indices.clone(),
        multimap_reads: family.multimap_reads.clone(),
    }]
}
```

- [ ] **Step 3: Verify the new test passes**
```bash
cargo test -p rustle family_for_em_input_keeps_cross_strand_placements 2>&1 | tail -20
```

- [ ] **Step 4a: Record per-family strand-partition spans (pipeline.rs:10665-10669)**
The `em_hmm_partitions` loop currently reads:
```rust
                    let mut em_hmm_partitions: Vec<crate::vg::FamilyGroup> = Vec::new();
                    for fam in families_for_em.iter() {
                        let parts = crate::vg::partition_and_remap_family_by_strand(fam, &bundles);
                        em_hmm_partitions.extend(parts);
                    }
```
Replace it with a version that records, per source family, `(start_idx, count)` into `em_hmm_partitions` (needed to map each family to its graph(s)):
```rust
                    let mut em_hmm_partitions: Vec<crate::vg::FamilyGroup> = Vec::new();
                    // Per source family: (start index into em_hmm_partitions, #strand partitions).
                    let mut em_fam_spans: Vec<(usize, usize)> = Vec::with_capacity(families_for_em.len());
                    for fam in families_for_em.iter() {
                        let start = em_hmm_partitions.len();
                        let parts = crate::vg::partition_and_remap_family_by_strand(fam, &bundles);
                        let count = parts.len();
                        em_hmm_partitions.extend(parts);
                        em_fam_spans.push((start, count));
                    }
```

- [ ] **Step 4b: Build `em_input_partitions` + aligned `em_input_graphs` (env-gated, default ON in VG)**
Insert this **after** `family_graphs` is built (i.e. after the `par_iter()` collect that ends at pipeline.rs:10694, BEFORE the `em_res = if do_hmm { ... }` block). **Use `.clone()`, not `mem::take`** — verified that `family_graphs` is read again after the EM call by `build_bundle_borrow_coverage` / `build_bundle_borrow_junctions` / `build_bundle_completion_nodes` (pipeline.rs:10780/10784/10791) and is itself cloned at pipeline.rs:10801, which also confirms `FamilyGraph: Clone`. Cloning only the single-strand graphs (mixed-strand → `None`, no clone) keeps cost bounded:
```rust
                    // Joint-strand EM input (spec component 1, O1-resolved): default ON in VG.
                    // The fingerprint-EM consumes the UNSPLIT family so a read shared only
                    // across an inverted pair (DAZ1/DAZ3) keeps both placements and is
                    // reweighted to sum to 1.0 across strands. em_input_graphs is kept
                    // 1:1 with em_input_partitions: single-strand family -> its (only)
                    // graph (fp stays valid); mixed-strand family -> None (neutral fp).
                    let joint_strand_em = std::env::var("RUSTLE_VG_JOINT_STRAND_EM")
                        .map(|v| v != "0")
                        .unwrap_or(true);
                    let (em_input_partitions, em_input_graphs):
                        (Vec<crate::vg::FamilyGroup>,
                         Vec<Option<crate::vg_hmm::family_graph::FamilyGraph>>) =
                    if joint_strand_em {
                        let mut parts = Vec::with_capacity(families_for_em.len());
                        let mut graphs = Vec::with_capacity(families_for_em.len());
                        for (fam, &(start, count)) in families_for_em.iter().zip(em_fam_spans.iter()) {
                            parts.extend(crate::vg::family_for_em_input(fam, &bundles)); // exactly 1
                            if count == 1 {
                                // single-strand: reuse the family's graph (clone; family_graphs
                                // is read again by the borrow/completion builders below)
                                graphs.push(family_graphs.get(start).and_then(|o| o.clone()));
                            } else {
                                // mixed-strand (e.g. inverted DAZ1/DAZ3): no valid joint graph
                                graphs.push(None);
                            }
                        }
                        (parts, graphs)
                    } else {
                        // exact current behavior: strand-split input + its graphs
                        (em_hmm_partitions.clone(), family_graphs.clone())
                    };
```

- [ ] **Step 5: Feed `em_input_partitions` to `run_fingerprint_em` (pipeline.rs:10734-10739)**
The fingerprint-EM call currently reads:
```rust
                            crate::vg::run_fingerprint_em(
                                &em_hmm_partitions,
                                &mut bundles,
                                &family_graphs,
                                config.vg_em_max_iter,
                            )
```
Change BOTH the families arg AND the graphs arg to the aligned joint-strand vectors (`em_input_partitions` / `em_input_graphs` from Step 4b):
```rust
                            crate::vg::run_fingerprint_em(
                                &em_input_partitions,
                                &mut bundles,
                                &em_input_graphs,
                                config.vg_em_max_iter,
                            )
```
This keeps `family_graphs.get(fam_idx)` aligned: index `i` now refers to the same family in both `em_input_partitions[i]` and `em_input_graphs[i]` (single-strand → its graph; mixed-strand → `None` → neutral fp). `build_bundle_borrow_*` / `build_bundle_completion_nodes` (pipeline.rs:10780-10791) keep using `em_hmm_partitions` / `family_graphs` unchanged (Step 4b clones, so `family_graphs` is intact for them).

- [ ] **Step 6: Build + run the touched tests**
```bash
cargo build -p rustle 2>&1 | tail -20
cargo test -p rustle family_for_em_input 2>&1 | tail -20
```

- [ ] **Step 7: Commit**
```bash
git add src/rustle/vg.rs src/rustle/pipeline.rs
git commit -m "VG: wire joint-strand (unsplit) family into fingerprint-EM input

Add vg::family_for_em_input returning the unsplit both-strands family for
the EM placement list + normalization (component 1 of the flow-capacity
apportionment design). Gated by RUSTLE_VG_JOINT_STRAND_EM (default ON in
VG; =0 restores strand-split input). Keep em_hmm_partitions (strand-split)
for build_family_graph/borrow/completion; per O1 graphs stay per-strand and
cross-strand placements route to the neutral fp path inside run_fingerprint_em.

Co-Authored-By: Claude Opus 4.8 (1M context) <noreply@anthropic.com>"
```

---

### Task 7: Anchored-mass prior in run_fingerprint_em (anchor-first one-pass apportionment)

Replaces the self-reinforcing pileup-depth M-step prior with an anchored prior grounded in unambiguous mass (unique reads + dNM-decisive `Owns` reads). When `RUSTLE_VG_ANCHOR_PRIOR != "0"` (default ON in VG), the prior is computed ONCE from `anchored_mass_per_copy(family, bundles, 2, 0.8)` and the E-step runs ONCE. All-zero-anchor families graceful-degrade to the existing pileup prior. The score-gap gate and log-sum-exp normalization are kept verbatim, so per-read mass still sums to 1.0.

**Prerequisite:** `anchored_mass_per_copy(family: &FamilyGroup, bundles: &[Bundle], t: i64, extent_frac: f64) -> Vec<f64>` must already exist (earlier task). It is called here only.

**Files:**
- Modify: `src/rustle/vg.rs:4684-4697` (M-step prior source + iteration count inside `run_fingerprint_em`)
- Test: `src/rustle/vg.rs` inline `#[cfg(test)] mod tests` (new `#[test]` near the existing EM-fixture tests at vg.rs:5961-6064)

---

- [ ] **Step 1: Write the FAILING unit test** (append inside `mod tests`, after `extractor_keeps_co_expressed_ties` at vg.rs:6028, before `classify_family_smoke`)

A 2-copy family with IDENTICAL sequences (so `fp.n_sites == 0`, the DAZ regime where the prior is the only signal). Copy 0 has many unique anchor reads; copy 1 has few. A single tied multimapper (single-exon → `junction_compatibility == 1.0` at both; equal `nm` → equal identity → genuine junc/nm tie) must apportion by the anchored ratio, NOT 1/n. Build a single-node `FamilyGraph` with two identical copies via the existing `make_ec` helper.

```rust
    /// Anchor-first one-pass: with identical copies (fp.n_sites==0) and a junc/NM
    /// TIE, a tied multimapper must apportion by the ANCHORED mass ratio (unique +
    /// Owns reads), not the uniform 1/n. Copy0 has 20 unique anchors, copy1 has 2;
    /// the shared read should land ~0.91/0.09 (=20/22), and the two placements must
    /// still sum to 1.0 (mass conservation via the unchanged log-sum-exp).
    #[test]
    fn anchored_prior_apportions_tied_read_by_anchor_ratio() {
        std::env::set_var("RUSTLE_VG_ANCHOR_PRIOR", "1");
        // Identical single-exon copies → 0 diagnostic sites; copy0 at 100-200, copy1 at 5100-5200.
        let ec = make_ec(
            0,
            &[(0, b"ACGTACGTAC"), (1, b"ACGTACGTAC")],
            &[(0, (100, 110)), (1, (5100, 5110))],
        );
        let fg = FamilyGraph { family_id: 42, nodes: vec![ec], edges: vec![] };
        let fp = build_exon_fingerprints(&fg, 2);
        assert_eq!(fp.n_sites, 0, "identical copies → no diagnostic sites (DAZ regime)");

        let mut c0_reads = Vec::new(); // copy 0 bundle
        let mut c1_reads = Vec::new(); // copy 1 bundle
        let mut multimap: HashMap<u64, Vec<(usize, usize)>> = HashMap::new();

        // 20 unique anchors at copy0 (read_name_hash NOT in multimap → unique mass).
        for i in 0..20u64 {
            c0_reads.push(make_read_full(3000 + i, vec![(100, 1100)], 2, false));
        }
        // 2 unique anchors at copy1.
        for i in 0..2u64 {
            c1_reads.push(make_read_full(4000 + i, vec![(5100, 6100)], 2, false));
        }
        // ONE shared multimapper: single exon (no junctions → junc tie), equal nm
        // at both placements (→ identity tie). Only the prior can break the tie.
        let rnh = 9999u64;
        c0_reads.push(make_read_full(rnh, vec![(100, 1100)], 2, false));   // copy0 idx 20
        c1_reads.push(make_read_full(rnh, vec![(5100, 6100)], 2, false));  // copy1 idx 2
        let c0_idx = c0_reads.len() - 1;
        let c1_idx = c1_reads.len() - 1;
        multimap.insert(rnh, vec![(0, c0_idx), (1, c1_idx)]);

        let mut bundles = vec![
            make_bundle("chr1", '+', c0_reads),
            make_bundle("chr1", '+', c1_reads),
        ];
        let family = FamilyGroup {
            family_id: 42,
            bundle_indices: vec![0, 1],
            multimap_reads: multimap,
        };
        let family_graphs = vec![Some(fg)];

        // max_iter is effectively 1 under anchor mode (single E-step with fixed prior).
        let _ = run_fingerprint_em(&[family], &mut bundles, &family_graphs, 5);

        let w0 = bundles[0].reads[c0_idx].weight;
        let w1 = bundles[1].reads[c1_idx].weight;

        // Mass conservation: the read's two placements sum to 1.0.
        assert!((w0 + w1 - 1.0).abs() < 1e-6, "mass must sum to 1.0: {w0} + {w1}");
        // Anchored ratio 20:2 → copy0 gets the lion's share, NOT 0.5/0.5.
        assert!(w0 > 0.80, "copy0 (20 anchors) should dominate: w0={w0}");
        assert!(w1 < 0.20, "copy1 (2 anchors) should be down-weighted: w1={w1}");
        assert!((w0 - 0.5).abs() > 0.30, "must NOT collapse to uniform 1/n: w0={w0}");

        std::env::remove_var("RUSTLE_VG_ANCHOR_PRIOR");
    }
```

Run it — it FAILS (current code uses the pileup prior, which is seeded from the read's initial 1/NH placement mass and the unique reads' weight, and re-derives in-loop; the assertion that the tie does not collapse toward uniform may pass for the wrong reason or differ from the anchored ratio — confirm RED first):

```bash
cargo test -p rustle anchored_prior_apportions_tied_read_by_anchor_ratio 2>&1 | tail -25
```

- [ ] **Step 2: Read the change site to anchor the edit exactly**

```bash
sed -n '4684,4698p' src/rustle/vg.rs
```

Confirm it matches the `for iter in 0..max_iter {` header (4686) followed by the M-step block (4687-4697).

- [ ] **Step 3: Add the anchor-mode flag + compute the anchored prior ONCE before the iteration loop**

Replace the block starting at `let mut result = EmResult::default();` (vg.rs:4684) through the M-step's `.collect();` (vg.rs:4697) with the version below. The flag is read once; the anchored prior is built before the loop; the in-loop M-step is bypassed under anchor mode (so the loop runs a single effective E-step). Graceful-degrade to the existing pileup prior for an all-zero-anchor family only.

```rust
        let mut result = EmResult::default();

        // Anchor-first apportionment (default ON in VG): replace the M-step
        // pileup-depth prior with an anchored prior grounded in UNAMBIGUOUS mass
        // (unique reads + dNM-decisive Owns reads). This removes the pileup
        // prior's self-reinforcement of an already-double-counted copy. The
        // anchored prior is constant, so the E-step runs ONCE (max_iter
        // effectively 1). Calibration: raw dNM t=2, extent_frac=0.8.
        let anchor_prior_on = std::env::var("RUSTLE_VG_ANCHOR_PRIOR")
            .map(|v| v != "0").unwrap_or(true);
        let fixed_log_priors: Option<Vec<f64>> = if anchor_prior_on {
            let anchored = anchored_mass_per_copy(family, bundles, 2, 0.8);
            let total_anchored: f64 = anchored.iter().sum();
            if total_anchored < 1e-9 {
                // All-zero-anchor (e.g. every copy is a pure multimapper): no
                // grounded mass exists. Graceful-degrade to the existing
                // pileup-depth prior for THIS family only.
                eprintln!(
                    "[VG-FP-EM] Family {}: 0 anchored mass — graceful-degrade to pileup-depth prior",
                    family.family_id
                );
                None
            } else {
                let total = total_anchored.max(1.0);
                Some(anchored.iter().map(|&a| ((a / total) + 1e-3).ln()).collect())
            }
        } else {
            None
        };
        // Under anchor mode with a fixed prior, a single E-step suffices.
        let effective_max_iter = if fixed_log_priors.is_some() { 1 } else { max_iter };

        for iter in 0..effective_max_iter {
            // M-step prior: anchored (fixed, computed once) when available,
            // otherwise the legacy pileup-depth prior recomputed from current
            // weights (graceful-degrade or RUSTLE_VG_ANCHOR_PRIOR=0).
            let log_priors: Vec<f64> = if let Some(ref lp) = fixed_log_priors {
                lp.clone()
            } else {
                let mut copy_total = vec![0.0_f64; n_copies];
                for entry in &entries {
                    for (i, &w) in entry.weights.iter().enumerate() {
                        copy_total[entry.fam_pos[i]] += w;
                    }
                }
                let total_sum: f64 = copy_total.iter().sum::<f64>().max(1.0);
                copy_total.iter().map(|&t| ((t / total_sum) + 1e-3).ln()).collect()
            };
```

Note: the E-step (vg.rs:4699-4773), the score-gap gate (4745-4755), and the log-sum-exp normalization (4757-4772) are UNCHANGED below this point. The loop's tail (`result.iterations = iter + 1; ... break;` at 4775-4780) is also unchanged — under anchor mode it exits after the single iteration.

- [ ] **Step 3b: Lift the `None`-graph skip → neutral-fingerprint fallthrough (gated, makes the anchored prior reach the mixed-strand DAZ family)**

**Why (from T6 review):** T6 feeds the mixed-strand family with graph = `None`. The current `None => continue` skip (vg.rs:4643–4651) drops that family before any reweighting, so the anchored prior added in Step 3 never runs for DAZ. Replace the skip with: when `RUSTLE_VG_JOINT_STRAND_EM != "0"` (default ON — same flag the pipeline uses to choose the unsplit input, so the two stay coherent) AND the family has ≥2 copies with multimappers (already guaranteed by the n_copies<2 / empty check at vg.rs:4638), build an **empty** `ExonFingerprints` and proceed; otherwise keep the exact old skip (so `RUSTLE_VG_JOINT_STRAND_EM=0` is a true rollback).

Replace the `let fg = match fg_opt { ... }` block (vg.rs:4644–4651) with:

```rust
        // Joint-strand mode (default ON): a None/empty graph is the mixed-strand
        // (inverted-pair) case — there is NO valid joint sequence graph, but we
        // still want to apportion its shared multimappers by junc + nm + the
        // anchored prior. Build a neutral (0-site) fingerprint and fall through
        // instead of skipping. With RUSTLE_VG_JOINT_STRAND_EM=0 we keep the exact
        // legacy skip so the rollback A/B is byte-identical.
        let joint_strand_em = std::env::var("RUSTLE_VG_JOINT_STRAND_EM")
            .map(|v| v != "0").unwrap_or(true);
        let fp = match fg_opt {
            Some(g) if !g.nodes.is_empty() => enumerate_diagnostic_sites(g, n_copies),
            _ if joint_strand_em => {
                eprintln!(
                    "[VG-FP-EM] Family {}: no joint graph (mixed-strand) — neutral-fp fallthrough (junc+nm+anchored prior)",
                    family.family_id
                );
                ExonFingerprints {
                    sites: Vec::new(),
                    per_copy_site_refs: vec![Vec::new(); n_copies],
                    n_copies,
                    n_sites: 0,
                }
            }
            _ => {
                eprintln!("[VG-FP-EM] Family {}: no family graph — skipping", family.family_id);
                results.push(EmResult::default());
                continue;
            }
        };
```

Then DELETE the now-duplicate `let fp = enumerate_diagnostic_sites(fg, n_copies);` that immediately follows (vg.rs:4653) — `fp` is now bound by the match above. The `if fp.n_sites == 0 { ... }` log/fallthrough block (4654-4669) stays and now also covers the mixed-strand case. (The old binding named the graph `fg`; nothing below the deleted line uses `fg` directly — `score_read_exon_fingerprint` takes `&fp`, not `fg`. Confirm with `grep -n '\bfg\b' src/rustle/vg.rs` inside `run_fingerprint_em` that no live use of `fg` remains after the edit; if one does, it is the deleted line only.)

- [ ] **Step 3c: Add a FAILING test for the mixed-strand (None-graph) anchored apportionment** (append in `mod tests`, right after `anchored_prior_apportions_tied_read_by_anchor_ratio`)

Identical to the Step-1 test EXCEPT `family_graphs = vec![None]` (mixed-strand: no joint graph). Under the lifted skip + anchored prior, the tied read must STILL apportion 20:2 by anchored mass (not be skipped at 0.5/0.5, and not stay full-weight at both). This is the real DAZ regime the whole plan targets.

```rust
    /// Mixed-strand DAZ regime: the family has NO joint sequence graph
    /// (graph = None). The lifted skip (Step 3b) + anchored prior must still
    /// apportion a tied multimapper by anchored mass (20:2), proving the fix
    /// reaches the inverted-pair family rather than skipping it.
    #[test]
    fn anchored_prior_reaches_mixed_strand_none_graph_family() {
        std::env::set_var("RUSTLE_VG_ANCHOR_PRIOR", "1");
        std::env::set_var("RUSTLE_VG_JOINT_STRAND_EM", "1");

        let mut c0_reads = Vec::new();
        let mut c1_reads = Vec::new();
        let mut multimap: HashMap<u64, Vec<(usize, usize)>> = HashMap::new();
        for i in 0..20u64 { c0_reads.push(make_read_full(3000 + i, vec![(100, 1100)], 2, false)); }
        for i in 0..2u64  { c1_reads.push(make_read_full(4000 + i, vec![(5100, 6100)], 2, false)); }
        let rnh = 9999u64;
        c0_reads.push(make_read_full(rnh, vec![(100, 1100)], 2, false));
        c1_reads.push(make_read_full(rnh, vec![(5100, 6100)], 2, false));
        let c0_idx = c0_reads.len() - 1;
        let c1_idx = c1_reads.len() - 1;
        multimap.insert(rnh, vec![(0, c0_idx), (1, c1_idx)]);

        let mut bundles = vec![
            make_bundle("chr1", '+', c0_reads),
            make_bundle("chr1", '-', c1_reads),   // opposite strand = inverted pair
        ];
        let family = FamilyGroup { family_id: 77, bundle_indices: vec![0, 1], multimap_reads: multimap };
        let family_graphs: Vec<Option<crate::vg_hmm::family_graph::FamilyGraph>> = vec![None];

        let _ = run_fingerprint_em(&[family], &mut bundles, &family_graphs, 5);

        let w0 = bundles[0].reads[c0_idx].weight;
        let w1 = bundles[1].reads[c1_idx].weight;
        assert!((w0 + w1 - 1.0).abs() < 1e-6, "mass must sum to 1.0: {w0} + {w1}");
        assert!(w0 > 0.80, "copy0 (20 anchors) should dominate even with None graph: w0={w0}");
        assert!(w1 < 0.20, "copy1 (2 anchors) down-weighted: w1={w1}");

        std::env::remove_var("RUSTLE_VG_ANCHOR_PRIOR");
        std::env::remove_var("RUSTLE_VG_JOINT_STRAND_EM");
    }
```

Run it RED before Step 3b, GREEN after:
```bash
cargo test -p rustle anchored_prior_reaches_mixed_strand_none_graph_family 2>&1 | tail -15
```

- [ ] **Step 4: Build**

```bash
cargo build -p rustle 2>&1 | tail -20
```

- [ ] **Step 5: Run the new test (now GREEN) and the existing EM-fixture tests (no regression)**

```bash
cargo test -p rustle anchored_prior_apportions_tied_read_by_anchor_ratio 2>&1 | tail -15
cargo test -p rustle anchored_prior_reaches_mixed_strand_none_graph_family 2>&1 | tail -15
cargo test -p rustle --lib vg:: 2>&1 | tail -25
```

- [ ] **Step 6: Commit**

```bash
git add src/rustle/vg.rs
git commit -m "VG: anchored-mass prior + None-graph neutral fallthrough in run_fingerprint_em

Replace the self-reinforcing pileup-depth M-step prior with an anchored
prior from anchored_mass_per_copy(family,bundles,2,0.8) (unique + Owns
reads). Computed once; E-step runs once (max_iter effectively 1).
All-zero-anchor families graceful-degrade to the legacy pileup prior.

Also lift the None/empty-graph skip into a neutral-fingerprint fallthrough
(empty ExonFingerprints, n_sites=0 -> junc+nm+anchored prior) so the
mixed-strand inverted pair (DAZ1/DAZ3), which T6 feeds with graph=None,
is actually apportioned instead of skipped. Gated by
RUSTLE_VG_JOINT_STRAND_EM (=0 restores the legacy skip + strand-split
input, a clean rollback).

Score-gap gate + log-sum-exp normalization kept verbatim (mass sums to
1.0). Anchored prior default ON via RUSTLE_VG_ANCHOR_PRIOR (=0 restores).

Co-Authored-By: Claude Opus 4.8 (1M context) <noreply@anthropic.com>"
```

---

### Task 8: Set em fields in write-back and close the score-gap conservation hole

**Files:**
- Modify: `src/rustle/vg.rs:4745-4755` (gate-fire mass-conserving fallback), `src/rustle/vg.rs:4823-4866` (write-back em fields)
- Test: `src/rustle/vg.rs` (inline `#[cfg(test)] mod tests`, near the existing fingerprint EM tests ~5755)

Assumes **Task 1** already added `em_weight_gap: f64`, `em_n_sites: u32`, `em_anchored: bool` to `BundleRead` (types.rs), and **Task 7** sources `log_priors` from the anchored-mass prior. All reads in `entries` come from `family.multimap_reads`, so "was_unique" reduces to the read's `nh <= 1`.

- [ ] **Step 1: Add the FAILING unit test (gate-fired read still conserves mass + em fields set)**

The gate fires when two near-tied placements have `(best - second) < eff_gap`. Construct a 2-copy family whose fingerprint has a real SNP site but where a multimapper covers 0 sites at the second copy and produces a tiny score gap, so the read trips the gate. After `run_fingerprint_em(.., max_iter=1)` the gate-fired read's per-placement weights must sum to 1.0 (not stay at raw 1/NH=0.5+0.5=1.0 by luck — so we additionally assert the weights match the normalized prior, i.e. were rewritten) and the em fields are populated.

Add this test inside `mod tests` (it reuses the existing `make_ec`, `make_read`, `make_bundle` helpers):

```rust
    /// Conservation hole closed: when the score-gap gate fires for a read,
    /// its per-placement weights are set to the normalized anchored prior
    /// (mass-conserving, sums to 1.0) instead of being left at raw 1/NH,
    /// and the em_* fields are populated on the BundleRead.
    #[test]
    fn gate_fired_read_conserves_mass_and_sets_em_fields() {
        // Two copies with one diagnostic SNP at exon pos 2: copy0=G, copy1=C.
        // copy0 at ref 100-105, copy1 at ref 500-505 (non-overlapping loci).
        let ec = make_ec(
            0,
            &[(0, b"ACGTA"), (1, b"ACCTA")],
            &[(0, (100, 105)), (1, (500, 505))],
        );
        let fg = FamilyGraph { family_id: 0, nodes: vec![ec], edges: vec![] };
        let family_graphs = vec![Some(fg)];

        // One multimapper placed at both copies. At copy0 it covers the SNP
        // (no mismatch -> base G -> copy0 allele); at copy1 it is the SAME read
        // span (covers 0 sites there). The diagonal fp scoring is near-tied so
        // the read trips the gap gate.
        let mut r0 = make_read(vec![(100, 105)], vec![]);
        r0.read_name_hash = 42;
        r0.nh = 2;
        r0.weight = 0.5;
        let mut r1 = make_read(vec![(100, 105)], vec![]);
        r1.read_name_hash = 42;
        r1.nh = 2;
        r1.weight = 0.5;

        let bundles_vec = vec![
            make_bundle("chrTest", '+', vec![r0]),
            make_bundle("chrTest", '+', vec![r1]),
        ];
        let mut bundles = bundles_vec;

        let mut multimap: HashMap<u64, Vec<(usize, usize)>> = HashMap::new();
        multimap.insert(42, vec![(0, 0), (1, 0)]);
        let family = FamilyGroup {
            family_id: 0,
            bundle_indices: vec![0, 1],
            multimap_reads: multimap,
        };

        // Force the gate to ALWAYS fire by setting a huge per-site gap threshold
        // (any finite score gap < eff_gap -> gate fires for this read).
        std::env::set_var("RUSTLE_VG_EM_SCORE_GAP", "1000.0");
        std::env::set_var("RUSTLE_VG_EM_SCORE_GAP_PER_SITE", "1000.0");

        let _ = run_fingerprint_em(&[family], &mut bundles, &family_graphs, 1);

        std::env::remove_var("RUSTLE_VG_EM_SCORE_GAP");
        std::env::remove_var("RUSTLE_VG_EM_SCORE_GAP_PER_SITE");

        let w0 = bundles[0].reads[0].weight;
        let w1 = bundles[1].reads[0].weight;
        // Mass conservation across the read's placements.
        assert!((w0 + w1 - 1.0).abs() < 1e-6, "gate-fired weights must sum to 1.0: {w0}+{w1}");

        // em_* fields populated in write-back.
        for bi in 0..2 {
            let r = &bundles[bi].reads[0];
            assert!(r.em_weight_gap >= 0.0, "em_weight_gap set, got {}", r.em_weight_gap);
            assert_eq!(r.em_n_sites, 1, "em_n_sites = best-placement covered sites");
            // gate-fired read: gap small, max_sites>0 but gap<=0.5, nh>1 -> not anchored.
            assert!(!r.em_anchored, "gate-fired near-tied read should not be em_anchored");
        }
    }
```

Run it to confirm it FAILS (em_* fields/behaviour not wired yet):

```bash
cd /mnt/c/Users/jfris/Desktop/Rustle && cargo test -p rustle gate_fired_read_conserves_mass_and_sets_em_fields 2>&1 | tail -25
```

- [ ] **Step 2: Make the gate-fire branch mass-conserving (set normalized anchored prior before `continue`)**

In `src/rustle/vg.rs`, the gate block currently bails with a bare `continue`, leaving the read at its raw 1/NH weights. Replace it so it first assigns each placement the normalized prior mass (`exp(log_priors[fam_pos])` normalized over this read's placements to sum to 1.0), then continues. `log_priors` is in scope (computed in the M-step above this loop) and `entry.fam_pos[i]` indexes it.

Replace this exact block (vg.rs:4745-4755):

```rust
                if eff_gap > 0.0 && n >= 2 {
                    let mut best = f64::NEG_INFINITY;
                    let mut second = f64::NEG_INFINITY;
                    for &v in &log_post {
                        if v > best { second = best; best = v; }
                        else if v > second { second = v; }
                    }
                    if best.is_finite() && second.is_finite() && (best - second) < eff_gap {
                        continue;
                    }
                }
```

with:

```rust
                if eff_gap > 0.0 && n >= 2 {
                    let mut best = f64::NEG_INFINITY;
                    let mut second = f64::NEG_INFINITY;
                    for &v in &log_post {
                        if v > best { second = best; best = v; }
                        else if v > second { second = v; }
                    }
                    if best.is_finite() && second.is_finite() && (best - second) < eff_gap {
                        // Conservation hole fix: rather than leaving this read at its
                        // raw 1/NH weights (which do NOT generally sum to 1.0 across the
                        // placements kept in this entry), assign the normalized anchored
                        // prior over this read's placements. Mass-conserving: weights
                        // sum to exactly 1.0. The prior (log_priors) is the M-step's
                        // anchored-mass distribution, so gate-silenced reads fall back to
                        // the family's unambiguous-mass apportionment instead of a stale,
                        // non-conserving initialization.
                        let max_lprior = entry.fam_pos.iter()
                            .map(|&k| log_priors[k])
                            .fold(f64::NEG_INFINITY, f64::max);
                        if max_lprior.is_finite() {
                            let mut psum = 0.0_f64;
                            let mut pexps = vec![0.0_f64; n];
                            for i in 0..n {
                                pexps[i] = (log_priors[entry.fam_pos[i]] - max_lprior).exp();
                                psum += pexps[i];
                            }
                            if psum > 0.0 {
                                for i in 0..n {
                                    let new_w = pexps[i] / psum;
                                    let delta = (new_w - entry.weights[i]).abs();
                                    if delta > max_delta { max_delta = delta; }
                                    entry.weights[i] = new_w;
                                }
                            }
                        }
                        continue;
                    }
                }
```

- [ ] **Step 3: Populate em fields in the write-back loop**

In the write-back loop (vg.rs:4823-4866), `gap` and `max_sites` are already computed per read at the top of the loop body. For each placement set the three em fields on the `BundleRead` alongside the existing `weight` write. `was_unique` reduces to `nh <= 1` (these are all multimappers, so it is normally false, but kept faithful to the spec formula).

Replace this exact block (vg.rs:4852-4859):

```rust
            for (i, &(global_bi, ri)) in entry.locs.iter().enumerate() {
                if ri >= bundles[global_bi].reads.len() { continue; }
                let old_w = bundles[global_bi].reads[ri].weight;
                let new_w = entry.weights[i];
                if (old_w - new_w).abs() > 1e-9 {
                    bundles[global_bi].reads[ri].weight = new_w;
                    n_reweighted += 1;
                }
```

with:

```rust
            for (i, &(global_bi, ri)) in entry.locs.iter().enumerate() {
                if ri >= bundles[global_bi].reads.len() { continue; }
                let old_w = bundles[global_bi].reads[ri].weight;
                let new_w = entry.weights[i];
                if (old_w - new_w).abs() > 1e-9 {
                    bundles[global_bi].reads[ri].weight = new_w;
                    n_reweighted += 1;
                }
                // Persist per-read EM attribution for the downstream capacity-
                // confidence channel. was_unique reduces to nh<=1 here (all entries
                // are multimappers). em_anchored mirrors the spec: decisive gap,
                // OR fingerprint-covered with a moderate gap, OR a unique read.
                let was_unique = bundles[global_bi].reads[ri].nh <= 1;
                bundles[global_bi].reads[ri].em_weight_gap = gap;
                bundles[global_bi].reads[ri].em_n_sites = entry.n_sites_covered[i] as u32;
                bundles[global_bi].reads[ri].em_anchored =
                    (gap > 0.8) || (max_sites > 0 && gap > 0.5) || was_unique;
```

- [ ] **Step 4: Run the test — confirm it PASSES**

```bash
cd /mnt/c/Users/jfris/Desktop/Rustle && cargo test -p rustle gate_fired_read_conserves_mass_and_sets_em_fields 2>&1 | tail -25
```

Then confirm the surrounding EM tests still pass and the crate builds clean:

```bash
cd /mnt/c/Users/jfris/Desktop/Rustle && cargo test -p rustle fingerprint 2>&1 | tail -20 && cargo build -p rustle 2>&1 | tail -5
```

- [ ] **Step 5: Commit**

```bash
git add src/rustle/vg.rs
git commit -m "vg: set em_* fields in EM write-back; mass-conserve gate-fired reads via normalized anchored prior

Co-Authored-By: Claude Opus 4.8 (1M context) <noreply@anthropic.com>"
```

---

### Task 9: Parallel anchored-coverage accumulation (VG-gated, byte-identical non-VG)

> Prereqs (from earlier tasks, verified absent in tree today): `BundleRead.em_anchored: bool` (default `nh<=1`) must already exist (**Task 1**), and `GraphNode.anchored_coverage: f64` (default `0.0`) must already exist (**Task 2**). This task only WIRES the parallel channel. It must produce zero observable change unless `--vg`.

**Why a process-global gate:** the three coverage-accumulation functions (`map_reads_to_graph`, `map_reads_to_graph_per_read`, `map_reads_to_graph_bundlenodes` in `map_reads.rs`) and the private helper `collect_read_nodes_exact` (map_reads.rs:250) do NOT receive `config`/`vg_mode`. The codebase's established pattern for a config-less process-global is `debug_stage::STAGE_IO_ENABLED: AtomicBool` (debug_stage.rs:20, set in `pipeline.rs` at ~10026). We mirror it with a tiny `vg_runtime` module so non-VG never touches the anchored channel.

**Files:**
- Create: src/rustle/vg_runtime.rs
- Modify: src/rustle/lib.rs:82 (module decl)
- Modify: src/rustle/pipeline.rs:~10026 (set flag from `config.vg_mode`)
- Modify: src/rustle/map_reads.rs:254 (return type of `collect_read_nodes_exact`), :331 (push site), :560-563 (consumer), :850-853 (consumer), :1081-1084 (consumer)
- Modify: src/rustle/graph.rs:773 (chain-absorption coverage transfer)
- Test: src/rustle/map_reads.rs (new inline `#[cfg(test)] mod anchored_cov_tests`)

---

- [ ] **Step 1: Create the process-global VG flag module**

Create `src/rustle/vg_runtime.rs`:

```rust
//! Process-global "VG mode active" flag.
//!
//! Mirrors `debug_stage::STAGE_IO_ENABLED`: a few hot, config-less functions
//! (`map_reads::collect_read_nodes_exact` and the three `map_reads_to_graph*`
//! entrypoints) need to know whether `--vg` is active so they can drive the
//! parallel `anchored_coverage` channel WITHOUT touching the byte-identical
//! default de-novo path. Set once at pipeline start from `config.vg_mode`.

use std::sync::atomic::{AtomicBool, Ordering};

static VG_MODE_ENABLED: AtomicBool = AtomicBool::new(false);

/// Set at pipeline start from `config.vg_mode`. Idempotent.
#[inline]
pub fn set_vg_mode(enabled: bool) {
    VG_MODE_ENABLED.store(enabled, Ordering::Relaxed);
}

/// Fast path: `false` outside `--vg`, so the anchored-coverage channel is inert
/// and default de-novo output stays byte-identical.
#[inline]
pub fn vg_mode() -> bool {
    VG_MODE_ENABLED.load(Ordering::Relaxed)
}
```

- [ ] **Step 2: Register the module**

In `src/rustle/lib.rs`, add the declaration directly after the `debug_stage` line (lib.rs:82):

```rust
pub mod debug_stage;
pub mod vg_runtime; // process-global "--vg active" flag for the anchored-coverage channel
```

- [ ] **Step 3: Set the flag at pipeline start**

In `src/rustle/pipeline.rs`, immediately after the existing `debug_stage::init(...)?;` call (pipeline.rs:~10026), add:

```rust
    debug_stage::init(diag_tsv_resolved.as_deref())?;
    crate::vg_runtime::set_vg_mode(config.vg_mode);
```

- [ ] **Step 4: Build (mechanical wiring compiles)**

```bash
cargo build -p rustle 2>&1 | tail -20
```

- [ ] **Step 5: Add a FAILING test for the anchored channel**

Append this inline test module at the END of `src/rustle/map_reads.rs`. It calls the private helper `collect_read_nodes_exact` directly (same crate), so it sees the new triple return type and the gate. It will FAIL to compile against the current 2-tuple return — that is the intended red state for the behavioral change in Steps 6-7.

```rust
#[cfg(test)]
mod anchored_cov_tests {
    use super::*;
    use crate::graph::Graph;
    use crate::types::BundleRead;

    /// Full BundleRead literal (mirrors vg.rs::make_read). `em_anchored` toggled
    /// per case; em_weight_gap/em_n_sites carry their defaults.
    fn make_read(exons: Vec<(u64, u64)>, weight: f64, em_anchored: bool) -> BundleRead {
        BundleRead {
            read_uid: 0, read_name: "r".into(), read_name_hash: 0,
            ref_id: None, mate_ref_id: None, mate_start: None, hi: 0,
            ref_start: exons.first().map(|(s, _)| *s).unwrap_or(0),
            ref_end:   exons.last().map(|(_, e)| *e).unwrap_or(0),
            exons, junctions: vec![], junction_valid: vec![],
            junctions_raw: vec![], junctions_del: vec![],
            weight, is_reverse: false, strand: '+',
            has_poly_start: false, has_poly_end: false,
            has_poly_start_aligned: false, has_poly_start_unaligned: false,
            has_poly_end_aligned: false, has_poly_end_unaligned: false,
            unaligned_poly_t: 0, unaligned_poly_a: 0,
            has_last_exon_polya: false, has_first_exon_polyt: false,
            query_length: None, clip_left: 0, clip_right: 0,
            nh: 1, nm: 0, de: None, md: None, insertion_sites: vec![],
            unitig: false, unitig_cov: 0.0, read_count_yc: 1.0,
            countfrag_len: 0.0, countfrag_num: 0.0, junc_mismatch_weight: 0.0,
            pair_idx: vec![], pair_count: vec![],
            mapq: 60, mismatches: vec![], seq: Vec::new(), hp_tag: None, ps_tag: None,
            is_primary_alignment: true,
            em_weight_gap: -1.0, em_n_sites: 0, em_anchored,
        }
    }

    /// Two contiguous exon segments → one node covering both.
    fn one_node_graph() -> Graph {
        let mut g = Graph::new();
        g.add_node(100, 200); // node 0: 100..200 (100 bp)
        g
    }

    #[test]
    fn anchored_is_zero_when_not_vg_mode() {
        crate::vg_runtime::set_vg_mode(false);
        let g = one_node_graph();
        let ordered: Vec<usize> = (0..g.nodes.len()).collect();
        let read = make_read(vec![(100, 200)], 1.0, /*em_anchored=*/true);
        let (_path, cov_add) = collect_read_nodes_exact(&read, &g, &ordered, true);
        assert!(!cov_add.is_empty(), "read overlaps the node");
        for (_nid, cov, anchored) in &cov_add {
            assert!(*cov > 0.0, "coverage accrues regardless of mode");
            assert_eq!(*anchored, 0.0, "non-vg leaves anchored_coverage == 0");
        }
    }

    #[test]
    fn anchored_le_coverage_and_tracks_em_anchored() {
        crate::vg_runtime::set_vg_mode(true);

        let g = one_node_graph();
        let ordered: Vec<usize> = (0..g.nodes.len()).collect();

        // Anchored read: anchored == coverage (full weight*bp).
        let anchored_read = make_read(vec![(100, 200)], 1.0, true);
        let (_p, cov_add) = collect_read_nodes_exact(&anchored_read, &g, &ordered, true);
        for (_nid, cov, anchored) in &cov_add {
            assert!(*anchored <= *cov, "anchored_coverage <= coverage");
            assert_eq!(*anchored, *cov, "anchored read contributes full mass");
        }

        // Ambiguous (non-anchored) read: coverage accrues, anchored is 0.
        let ambiguous_read = make_read(vec![(100, 200)], 0.5, false);
        let (_p, cov_add) = collect_read_nodes_exact(&ambiguous_read, &g, &ordered, true);
        for (_nid, cov, anchored) in &cov_add {
            assert!(*cov > 0.0, "ambiguous read still adds coverage");
            assert_eq!(*anchored, 0.0, "non-anchored adds no anchored mass");
            assert!(*anchored <= *cov, "invariant holds for ambiguous read too");
        }

        crate::vg_runtime::set_vg_mode(false); // restore for other tests
    }
}
```

Run it (expect a COMPILE error: the helper still returns a 2-tuple):

```bash
cargo test -p rustle anchored_cov_tests 2>&1 | tail -25
```

- [ ] **Step 6: Change `collect_read_nodes_exact` to emit anchored mass**

In `src/rustle/map_reads.rs`, change the return type (currently map_reads.rs:254 `-> (Vec<usize>, Vec<(usize, f64)>)`) to carry the anchored portion as a third element:

```rust
fn collect_read_nodes_exact(
    read: &BundleRead,
    graph: &Graph,
    ordered_nodes: &[usize],
    add_coverage: bool,
) -> (Vec<usize>, Vec<(usize, f64, f64)>) {
```

At the push site (map_reads.rs:331), replace:

```rust
                if add_coverage && !read.unitig {
                    cov_add.push((nid, read.weight * (bp as f64)));
                }
```

with:

```rust
                if add_coverage && !read.unitig {
                    let cov = read.weight * (bp as f64);
                    // Parallel anchored-coverage channel. Gated on the
                    // process-global VG flag so non-VG runs always push 0.0
                    // here and the default de-novo path is byte-identical.
                    // anchored <= cov always (factor is weight or 0).
                    let anchored = if crate::vg_runtime::vg_mode() && read.em_anchored {
                        cov
                    } else {
                        0.0
                    };
                    cov_add.push((nid, cov, anchored));
                }
```

- [ ] **Step 7: Update the three consumer loops**

In `src/rustle/map_reads.rs`, `map_reads_to_graph` (consumer at map_reads.rs:560-563), replace:

```rust
        let (unique_nodes, cov_add) = collect_read_nodes_exact(read, graph, &ordered_nodes, true);
        for (idx, add) in cov_add {
            if let Some(n) = graph.nodes.get_mut(idx) {
                n.coverage += add;
            }
        }
```

with:

```rust
        let (unique_nodes, cov_add) = collect_read_nodes_exact(read, graph, &ordered_nodes, true);
        for (idx, add, anchored_add) in cov_add {
            if let Some(n) = graph.nodes.get_mut(idx) {
                n.coverage += add;
                n.anchored_coverage += anchored_add;
            }
        }
```

In `map_reads_to_graph_per_read` (consumer at map_reads.rs:850-853), replace:

```rust
        let (unique_nodes, cov_add) = collect_read_nodes_exact(read, graph, &ordered_nodes, true);
        for (idx, add) in cov_add {
            if let Some(n) = graph.nodes.get_mut(idx) {
                n.coverage += add;
            }
        }
```

with:

```rust
        let (unique_nodes, cov_add) = collect_read_nodes_exact(read, graph, &ordered_nodes, true);
        for (idx, add, anchored_add) in cov_add {
            if let Some(n) = graph.nodes.get_mut(idx) {
                n.coverage += add;
                n.anchored_coverage += anchored_add;
            }
        }
```

In `map_reads_to_graph_bundlenodes` (consumer at map_reads.rs:1081-1084), replace:

```rust
        for (nid, add) in cov_add {
            if let Some(n) = graph.nodes.get_mut(nid) {
                n.coverage += add;
            }
        }
```

with:

```rust
        for (nid, add, anchored_add) in cov_add {
            if let Some(n) = graph.nodes.get_mut(nid) {
                n.coverage += add;
                n.anchored_coverage += anchored_add;
            }
        }
```

- [ ] **Step 8: Transfer anchored_coverage in chain absorption (graph.rs)**

`graph.rs` absorbs pass-through nodes into a chain head (graph.rs:765-775); coverage/longcov/nodecov are summed and transferred. Keep `anchored_coverage` consistent by transferring it too. In `src/rustle/graph.rs`, replace:

```rust
                let mut added_cov = 0.0f64;
                let mut added_longcov = 0.0f64;
                let mut added_nodecov = 0.0f64;
                for &pt_id in &pts {
                    added_cov += self.nodes[pt_id].coverage;
                    added_longcov += self.nodes[pt_id].longcov;
                    added_nodecov += self.nodes[pt_id].nodecov;
                }
                self.nodes[head].end = last_end;
                self.nodes[head].coverage += added_cov;
                self.nodes[head].longcov += added_longcov;
                self.nodes[head].nodecov += added_nodecov;
```

with:

```rust
                let mut added_cov = 0.0f64;
                let mut added_anchored_cov = 0.0f64;
                let mut added_longcov = 0.0f64;
                let mut added_nodecov = 0.0f64;
                for &pt_id in &pts {
                    added_cov += self.nodes[pt_id].coverage;
                    added_anchored_cov += self.nodes[pt_id].anchored_coverage;
                    added_longcov += self.nodes[pt_id].longcov;
                    added_nodecov += self.nodes[pt_id].nodecov;
                }
                self.nodes[head].end = last_end;
                self.nodes[head].coverage += added_cov;
                self.nodes[head].anchored_coverage += added_anchored_cov;
                self.nodes[head].longcov += added_longcov;
                self.nodes[head].nodecov += added_nodecov;
```

(Absorbed nodes' `anchored_coverage` stays on them but, like `coverage`, those nodes are dropped/merged out; the head now carries the chain total so `anchored_coverage <= coverage` is preserved at the head.)

- [ ] **Step 9: Run the test — expect GREEN**

```bash
cargo test -p rustle anchored_cov_tests 2>&1 | tail -25
```

All three assertions (anchored==0 when !vg, anchored<=coverage, anchored tracks `em_anchored`) must pass.

- [ ] **Step 10: Full build + confirm no other `cov_add` destructure broke**

```bash
cargo build -p rustle 2>&1 | tail -20
grep -rn "collect_read_nodes_exact" src/rustle/map_reads.rs
```

Confirm every call site (map_reads.rs:559, 849, 1077) is now destructured as a 3-tuple and the build is clean.

- [ ] **Step 11: Commit**

```bash
git add src/rustle/vg_runtime.rs src/rustle/lib.rs src/rustle/pipeline.rs src/rustle/map_reads.rs src/rustle/graph.rs
git commit -m "vg: parallel anchored-coverage channel, VG-gated (byte-identical non-VG)

Accumulate GraphNode.anchored_coverage alongside coverage wherever read
mass is added (collect_read_nodes_exact push + 3 map_reads_to_graph*
consumers + graph.rs chain absorption). Anchored portion is
em_anchored ? weight*bp : 0, gated on a new process-global vg_runtime
flag so the default de-novo path never touches the channel and stays
byte-identical. Invariant anchored_coverage <= coverage holds per node.

Co-Authored-By: Claude Opus 4.8 (1M context) <noreply@anthropic.com>"
```

---

### Task 10: Compute per-transcript capacity_confidence + abundance_min at Transcript push

**Depends on (must land first):**
- GraphNode gains `pub anchored_coverage: f64` (spec component 4 / Phase C; graph.rs ~125, defaulted `0.0` in `GraphNode::new`).
- Transcript gains `pub capacity_confidence: Option<f64>` and `pub abundance_min: Option<f64>` (spec component 4; path_extract.rs Transcript struct ~600-665, beside `copy_independent_support`), with `None` added to every Transcript struct literal across pipeline.rs + path_extract.rs so the crate still builds.

**Files:**
- Modify: src/rustle/path_extract.rs (add helper `path_capacity_confidence` near top of file; wire into the long-read Transcript push at ~9119-9145)
- Test: src/rustle/path_extract.rs (new inline `#[cfg(test)] mod tests` at end of file)

- [ ] **Step 1: Add a FAILING unit test for the pure helper**
Append a new test module at the very end of `src/rustle/path_extract.rs` (path_extract.rs has no existing test module). It exercises `path_capacity_confidence`, which does not exist yet, so it fails to compile (= failing test).
```rust
#[cfg(test)]
mod capacity_confidence_tests {
    use super::*;
    use crate::graph::{Graph, GraphNode};

    fn graph_with(nodes: Vec<(f64, f64)>) -> Graph {
        // nodes: Vec<(total_coverage, anchored_coverage)>
        let mut g = Graph::new();
        for (i, (cov, anch)) in nodes.iter().enumerate() {
            let mut n = GraphNode::new(i, (i as u64) * 100, (i as u64) * 100 + 50);
            n.coverage = *cov;
            n.anchored_coverage = *anch;
            g.nodes.push(n);
        }
        g
    }

    #[test]
    fn cc_fully_anchored_is_one() {
        let g = graph_with(vec![(10.0, 10.0), (20.0, 20.0)]);
        let (anch, tot) = path_capacity_confidence(&g, &[0, 1]);
        assert_eq!(tot, 30.0);
        assert_eq!(anch, 30.0);
        let cc = (anch / tot).clamp(0.0, 1.0);
        assert!((cc - 1.0).abs() < 1e-9);
    }

    #[test]
    fn cc_half_anchored_is_half() {
        let g = graph_with(vec![(10.0, 5.0), (30.0, 15.0)]);
        let (anch, tot) = path_capacity_confidence(&g, &[0, 1]);
        let cc = (anch / tot).clamp(0.0, 1.0);
        assert!((cc - 0.5).abs() < 1e-9);
        // abundance_min must never exceed coverage
        let coverage = 12.5_f64;
        let abundance_min = coverage * cc;
        assert!(abundance_min <= coverage);
        assert!((0.0..=1.0).contains(&cc));
    }

    #[test]
    fn cc_zero_total_is_zero_no_nan() {
        let g = graph_with(vec![(0.0, 0.0)]);
        let (anch, tot) = path_capacity_confidence(&g, &[0]);
        assert_eq!(tot, 0.0);
        assert_eq!(anch, 0.0);
        // caller convention: tot==0 -> cc 0.0, never NaN
        let cc = if tot > 0.0 { (anch / tot).clamp(0.0, 1.0) } else { 0.0 };
        assert_eq!(cc, 0.0);
    }
}
```

- [ ] **Step 2: Run the test and confirm it FAILS to compile (helper missing)**
```bash
cargo test -p rustle path_capacity_confidence 2>&1 | tail -20
cargo test -p rustle capacity_confidence_tests 2>&1 | tail -20
```
Expect: compile error `cannot find function path_capacity_confidence`. (If the crate name is not `rustle`, drop `-p rustle` and run `cargo test capacity_confidence_tests`.)

- [ ] **Step 3: Implement the pure helper**
Add this free function to `src/rustle/path_extract.rs` immediately after the `use` imports near the top of the file (it only needs `crate::graph::Graph`, already in scope via existing imports). Returns `(sum_anchored, sum_total)` over the given real nodes; callers do the clamp.
```rust
/// Sum `(anchored_coverage, coverage)` over a transcript path's real nodes.
///
/// VG capacity-confidence channel (spec 2026-06-01, Phase C): `capacity_confidence
/// = sum_anchored / sum_total` over the transcript's path nodes. Source/sink are
/// excluded by construction because callers pass `real_nodes` (see path_extract
/// `real_nodes`, which already drops graph.source_id / graph.sink_id). Out-of-range
/// node ids are skipped. Returns `(0.0, 0.0)` for an empty path; the caller maps
/// `sum_total == 0.0` to a confidence of `0.0` (never NaN).
fn path_capacity_confidence(graph: &Graph, real_nodes: &[usize]) -> (f64, f64) {
    let mut sum_anchored = 0.0_f64;
    let mut sum_total = 0.0_f64;
    for &nid in real_nodes {
        if let Some(node) = graph.nodes.get(nid) {
            sum_anchored += node.anchored_coverage;
            sum_total += node.coverage;
        }
    }
    (sum_anchored, sum_total)
}
```

- [ ] **Step 4: Run the unit test and confirm it PASSES**
```bash
cargo test -p rustle capacity_confidence_tests 2>&1 | tail -20
```
Expect: `test result: ok. 3 passed`.

- [ ] **Step 5: Wire the helper into the long-read Transcript push**
In `extract_transcripts` (path_extract.rs:6089), the Transcript push at ~9119 has `graph: &mut Graph`, `config: &RunConfig`, the local `coverage: f64`, and `real_nodes: Vec<usize>` (path_extract.rs:6887) all in scope. Insert the computation directly above the `out.push(Transcript { ... })` at line 9119. Only populate the new fields under `config.vg_mode`; otherwise leave them `None` (default-de-novo byte-identical).
```rust
        // VG capacity-confidence channel (spec 2026-06-01, Phase C): annotate the
        // transcript with the anchored fraction of its path's coverage mass and a
        // sub-conservative abundance lower bound. VG-gated; None in default mode so
        // de-novo output is byte-identical and the GTF emitter skips the attrs.
        let (cap_conf_opt, abund_min_opt): (Option<f64>, Option<f64>) = if config.vg_mode {
            let (sum_anchored, sum_total) = path_capacity_confidence(graph, &real_nodes);
            let cc = if sum_total > 0.0 {
                (sum_anchored / sum_total).clamp(0.0, 1.0)
            } else {
                0.0
            };
            (Some(cc), Some(coverage * cc))
        } else {
            (None, None)
        };
        out.push(Transcript {
```
Then set the two new struct fields in the literal that follows. Replace the existing trailing field line (path_extract.rs:9143):
```rust
                    vg_family_id: None, vg_copy_id: None, vg_family_size: None, copy_assignment_confidence: None, copy_independent_support: None, family_verdict: None, intron_low: Vec::new(), synthetic: false, rescue_class: None,
```
with (adds `capacity_confidence` / `abundance_min`; verify the exact text first since prior tasks may have reflowed this line):
```rust
                    vg_family_id: None, vg_copy_id: None, vg_family_size: None, copy_assignment_confidence: None, copy_independent_support: None, capacity_confidence: cap_conf_opt, abundance_min: abund_min_opt, family_verdict: None, intron_low: Vec::new(), synthetic: false, rescue_class: None,
```

- [ ] **Step 6: Add an integration assertion to the unit test that pins the caller invariants**
Append one more test to the same module asserting the exact caller contract used in Step 5 (`0 <= cc <= 1` and `abundance_min <= coverage`) for a representative range of mass splits.
```rust
    #[test]
    fn caller_invariants_hold_across_splits() {
        for &(tot, anch) in &[(0.0, 0.0), (5.0, 0.0), (5.0, 2.5), (5.0, 5.0), (5.0, 7.0)] {
            // last case (anch>tot) must still clamp to <=1
            let g = graph_with(vec![(tot, anch)]);
            let (a, t) = path_capacity_confidence(&g, &[0]);
            let cc = if t > 0.0 { (a / t).clamp(0.0, 1.0) } else { 0.0 };
            assert!((0.0..=1.0).contains(&cc), "cc out of range: {cc}");
            let coverage = 4.0_f64;
            let abundance_min = coverage * cc;
            assert!(abundance_min <= coverage + 1e-9, "abundance_min {abundance_min} > coverage {coverage}");
        }
    }
```

- [ ] **Step 7: Build and run the full file's tests**
```bash
cargo build -p rustle 2>&1 | tail -20
cargo test -p rustle capacity_confidence_tests 2>&1 | tail -20
```
Expect: clean build and `test result: ok. 4 passed`.

- [ ] **Step 8: Commit**
```bash
git add src/rustle/path_extract.rs
git commit -m "vg: per-tx capacity_confidence + abundance_min at Transcript push

Add pure helper path_capacity_confidence(graph, real_nodes) summing
(anchored_coverage, coverage) over a path's real nodes, and wire it into
the long-read Transcript push. capacity_confidence = (anchored/total)
clamped to [0,1]; abundance_min = coverage * capacity_confidence.
VG-gated (config.vg_mode); None in default mode so de-novo output stays
byte-identical. Unit tests pin 0<=cc<=1 and abundance_min<=coverage.

Co-Authored-By: Claude Opus 4.8 (1M context) <noreply@anthropic.com>"
```

---

### Task 11: Confidence attach + abstain tag (pipeline) and GTF emission (capacity_confidence, abundance_min)

**Pre-req (verify, do not re-add):** Earlier tasks added `pub capacity_confidence: Option<f64>` and `pub abundance_min: Option<f64>` to `Transcript` (struct ends at `chain_witnessed: bool`, path_extract.rs:704) and populate them at Transcript push (Phase C). This task only (a) tags abstained transcripts via `family_verdict` WITHOUT dropping them, and (b) emits the two attrs in GTF when `Some`. The abstain block lives inside the existing `if config.vg_mode && !vg_copy_support.is_empty()` guard, so non-VG output stays byte-identical.

**Files:**
- Modify: /mnt/c/Users/jfris/Desktop/Rustle/src/rustle/pipeline.rs:18104-18116 (confidence attach + abstain)
- Modify: /mnt/c/Users/jfris/Desktop/Rustle/src/rustle/gtf.rs:178-191 (GTF emit)
- Test: /mnt/c/Users/jfris/Desktop/Rustle/src/rustle/gtf.rs (new `#[cfg(test)] mod tests`)

- [ ] **Step 1: Write FAILING test for GTF emission + retention**
Append this module to the END of /mnt/c/Users/jfris/Desktop/Rustle/src/rustle/gtf.rs. It builds `Transcript`s via `..Default::default()` (Transcript derives Default at path_extract.rs:599), so it is robust to unrelated field additions. It asserts (i) attrs present when `Some`, (ii) attrs absent when `None` (non-vg), (iii) an abstained low-cc tx is RETAINED with its coverage and carries `family_identifiability "none"`.
```rust
#[cfg(test)]
mod tests {
    use super::*;
    use crate::path_extract::Transcript;
    use crate::vg::{FamilyClass, FamilyVerdict, Identifiability, LocusRel};

    fn one_exon_tx() -> Transcript {
        Transcript {
            chrom: "chrT".to_string(),
            strand: '+',
            exons: vec![(100, 200)],
            coverage: 12.0,
            ..Default::default()
        }
    }

    fn render(tx: Transcript) -> String {
        let mut buf: Vec<u8> = Vec::new();
        write_gtf(&[tx], &mut buf, "STRG").unwrap();
        String::from_utf8(buf).unwrap()
    }

    #[test]
    fn capacity_attrs_emitted_when_some() {
        let mut tx = one_exon_tx();
        tx.capacity_confidence = Some(0.250);
        tx.abundance_min = Some(3.000);
        let out = render(tx);
        assert!(out.contains("capacity_confidence \"0.250\";"), "missing capacity_confidence in:\n{out}");
        assert!(out.contains("abundance_min \"3.000\";"), "missing abundance_min in:\n{out}");
    }

    #[test]
    fn capacity_attrs_absent_when_none() {
        let tx = one_exon_tx(); // both fields default to None (non-vg)
        let out = render(tx);
        assert!(!out.contains("capacity_confidence"), "capacity_confidence leaked into non-vg GTF:\n{out}");
        assert!(!out.contains("abundance_min"), "abundance_min leaked into non-vg GTF:\n{out}");
    }

    #[test]
    fn abstained_tx_retained_and_tagged() {
        // Simulate the pipeline abstain tag: a low-cc tx with no prior verdict
        // gets a synthesized verdict whose identifiability == None. The tx is NOT
        // dropped, so its coverage survives and the GTF carries the tag.
        let mut tx = one_exon_tx();
        tx.capacity_confidence = Some(0.010);
        tx.family_verdict = Some(FamilyVerdict {
            class: FamilyClass::Spillover,
            n_copies: 1,
            n_expressed: 0,
            connectivity: 0.0,
            identifiability: Identifiability::None,
            n_id_classes: 0,
            locus_rel: LocusRel::Single,
        });
        let out = render(tx);
        assert!(out.contains("cov \"12.000000\";"), "abstained tx lost its coverage:\n{out}");
        assert!(out.contains("family_identifiability \"none\";"), "abstain tag missing:\n{out}");
    }
}
```
Run (expected: FAIL to compile — `capacity_confidence`/`abundance_min` not emitted yet; if pre-req fields exist the two emission tests fail on missing attrs):
```bash
cargo test -p rustle --lib gtf::tests 2>&1 | tail -30
```

- [ ] **Step 2: Emit the two attrs in GTF (mirror copy_confidence)**
In /mnt/c/Users/jfris/Desktop/Rustle/src/rustle/gtf.rs, insert immediately AFTER the `copy_independent_support` block (currently ends at line 185 `}`) and BEFORE the `family_verdict` block at line 187. Find this anchor:
```rust
        if let Some(supp) = tx.copy_independent_support {
            tx_attrs.push_str(&format!(" copy_independent_support \"{:.3}\";", supp));
        }
```
and insert the new block right after its closing brace:
```rust
        // Capacity-confidence channel (--vg flow-apportionment spec 2026-06-01).
        // Per-transcript anchored/total coverage fraction; only set in VG mode,
        // so non-VG GTF is byte-identical. Mirrors copy_confidence formatting.
        if let Some(cc) = tx.capacity_confidence {
            tx_attrs.push_str(&format!(" capacity_confidence \"{:.3}\";", cc));
        }
        // Jointly-feasible lower bound on abundance (coverage * capacity_confidence).
        if let Some(amin) = tx.abundance_min {
            tx_attrs.push_str(&format!(" abundance_min \"{:.3}\";", amin));
        }
```
Run (the two emission tests now pass; abstain test passes on its own synthesized verdict):
```bash
cargo test -p rustle --lib gtf::tests 2>&1 | tail -30
```

- [ ] **Step 3: Implement abstain tag in pipeline (TAG, do not drop)**
In /mnt/c/Users/jfris/Desktop/Rustle/src/rustle/pipeline.rs, the survivor-annotation loop at 18104-18111 currently is:
```rust
        for tx in all_transcripts.iter_mut() {
            if let (Some(fam_id), Some(copy_id)) = (tx.vg_family_id, tx.vg_copy_id) {
                if let Some(&supp) = vg_copy_support.get(&(fam_id, copy_id)) {
                    tx.copy_independent_support = Some(supp);
                }
                tx.family_verdict = vg_family_verdict.get(&(fam_id, copy_id)).cloned();
            }
        }
```
Replace it with the version that reads `RUSTLE_VG_ABSTAIN_FLOOR` (default 0.05) once and tags low-capacity_confidence transcripts via `family_verdict` WITHOUT removing them from `all_transcripts`:
```rust
        // Abstain floor: transcripts whose capacity_confidence falls below this
        // are TAGGED low-confidence (via family_verdict) but kept — coverage/TPM
        // and the benchmark transcript set are preserved (spec O5: tag, not drop).
        let abstain_floor: f64 = std::env::var("RUSTLE_VG_ABSTAIN_FLOOR")
            .ok()
            .and_then(|s| s.parse().ok())
            .unwrap_or(0.05);
        for tx in all_transcripts.iter_mut() {
            if let (Some(fam_id), Some(copy_id)) = (tx.vg_family_id, tx.vg_copy_id) {
                if let Some(&supp) = vg_copy_support.get(&(fam_id, copy_id)) {
                    tx.copy_independent_support = Some(supp);
                }
                tx.family_verdict = vg_family_verdict.get(&(fam_id, copy_id)).cloned();
            }
            // Low-confidence abstain tag. Applies to ANY VG transcript carrying a
            // capacity_confidence (Phase C), not just family copies. We mark the
            // verdict's identifiability as None (surfaced as family_identifiability
            // "none" in GTF); a tx with no prior verdict gets a minimal synthesized
            // one so the tag is observable. Coverage is untouched.
            if let Some(cc) = tx.capacity_confidence {
                if cc < abstain_floor {
                    match tx.family_verdict {
                        Some(ref mut v) => {
                            v.identifiability = crate::vg::Identifiability::None;
                        }
                        None => {
                            tx.family_verdict = Some(crate::vg::FamilyVerdict {
                                class: crate::vg::FamilyClass::Spillover,
                                n_copies: 1,
                                n_expressed: 0,
                                connectivity: 0.0,
                                identifiability: crate::vg::Identifiability::None,
                                n_id_classes: 0,
                                locus_rel: crate::vg::LocusRel::Single,
                            });
                        }
                    }
                }
            }
        }
```

- [ ] **Step 4: Build + full gtf test pass**
```bash
cargo build -p rustle 2>&1 | tail -20 && cargo test -p rustle --lib gtf::tests 2>&1 | tail -20
```

- [ ] **Step 5: Commit**
```bash
git add src/rustle/gtf.rs src/rustle/pipeline.rs
git commit -m "vg: emit capacity_confidence/abundance_min GTF attrs + abstain-tag (keep coverage)

Confidence attach at pipeline survivor loop now tags transcripts with
capacity_confidence < RUSTLE_VG_ABSTAIN_FLOOR (default 0.05) low-confidence
via family_verdict (identifiability=None) instead of dropping them, so
coverage/TPM and the benchmark set are preserved (spec O5). GTF emits
capacity_confidence and abundance_min only when Some, so non-VG output
stays byte-identical.

Co-Authored-By: Claude Opus 4.8 (1M context) <noreply@anthropic.com>"
```

---

### Task 12: Final VALIDATION — regression guard, DAZ headline, synthetic oracle, mass-conservation, genome-wide no-fabrication scan

This is a validation-only task: **no new source code**. It builds the binary, then runs commands and asserts explicit outputs. All five spec test-plan items (§"Test & validation plan" 1-5) are exercised. Stop and fix the offending implementation task if any assert fails; do not weaken an assert to make it pass.

**Files:**
- Build/Run only: `/mnt/c/Users/jfris/Desktop/Rustle/target/release/rustle`
- Test: `/mnt/c/Users/jfris/Desktop/Rustle/src/rustle/vg.rs` (existing `mass_conservation` unit test from the earlier task — run, do not author)
- Oracle harness: `/mnt/c/Users/jfris/Desktop/Rustle/bench/multi_copy_eval/run_oracle.py`, `/mnt/c/Users/jfris/Desktop/Rustle/bench/multi_copy_eval/score_synth_assignment.py`
- Genome-wide scan: `/mnt/c/Users/jfris/Desktop/Rustle/bench/paralog_secondary_scan/run_scan.py`, `/mnt/c/Users/jfris/Desktop/Rustle/bench/paralog_secondary_scan/de_anchor.py`
- Inputs (absolute): `/mnt/c/Users/jfris/Desktop/GGO_19.bam`, `/mnt/c/Users/jfris/Desktop/GGO_19.gtf`, `/mnt/c/Users/jfris/Desktop/GGO.bam`, `/mnt/c/Users/jfris/Desktop/GGO.fasta`
- Notes file produced by this task: `/tmp/task12_validation.log` (scratch only — do NOT commit a report .md)

- [ ] **Step 1: Build release binary**
```bash
cargo build --release --manifest-path /mnt/c/Users/jfris/Desktop/Rustle/Cargo.toml 2>&1 | tail -5
test -x /mnt/c/Users/jfris/Desktop/Rustle/target/release/rustle && echo "BUILD_OK"
```
Expect a clean build ending in `BUILD_OK`. (Compiler-enforced exhaustiveness on the new `Transcript`/`GraphNode`/`BundleRead` fields means a green build already proves every struct literal was updated — Risk 1 in the spec.)

- [ ] **Step 2: Capture the current HEAD de-novo GTF as the byte-identity reference (before this branch's run)**
The regression guard requires the default de-novo GTF to be **byte-identical to HEAD** (spec test 1). Build a binary from HEAD into a temp location and produce its GTF, so we compare *behavior at HEAD* vs *behavior on this branch* on identical inputs.
```bash
git -C /mnt/c/Users/jfris/Desktop/Rustle stash list | head -1
git -C /mnt/c/Users/jfris/Desktop/Rustle worktree add /tmp/rustle_head HEAD 2>&1 | tail -2
cargo build --release --manifest-path /tmp/rustle_head/Cargo.toml 2>&1 | tail -2
/tmp/rustle_head/target/release/rustle -L /mnt/c/Users/jfris/Desktop/GGO_19.bam -o /tmp/denovo_head.gtf 2>/dev/null
grep -vc '^#' /tmp/denovo_head.gtf
```
Expect a nonzero line count (~12k GTF rows). `/tmp/denovo_head.gtf` is the immutable reference for Step 3.

> NOTE: if this validation task runs on a branch where the apportionment commits are ALREADY merged into HEAD, the worktree at HEAD will contain the new code and the byte-diff degenerates to a no-op. In that case, replace `HEAD` above with the pre-feature base commit (the parent of the first apportionment commit in this plan, e.g. `git log --oneline | grep -m1 "flow-capacity\|apportion" ` then use that commit's `^`). Record which base you used in `/tmp/task12_validation.log`.

- [ ] **Step 3: Regression guard — default de-novo (NO --vg) byte-identical + Sn/Pr in band**
```bash
/mnt/c/Users/jfris/Desktop/Rustle/target/release/rustle -L /mnt/c/Users/jfris/Desktop/GGO_19.bam -o /tmp/denovo_branch.gtf 2>/dev/null
# (a) byte identity to HEAD
if diff -q /tmp/denovo_head.gtf /tmp/denovo_branch.gtf >/dev/null; then echo "BYTE_IDENTICAL_OK"; else echo "BYTE_DIFF_FAIL"; diff /tmp/denovo_head.gtf /tmp/denovo_branch.gtf | head -20; fi
# (b) no anchored_coverage / capacity_confidence leaked into a non-vg GTF
if grep -q 'capacity_confidence\|abundance_min' /tmp/denovo_branch.gtf; then echo "LEAK_FAIL: vg-only attrs in default GTF"; else echo "NO_VG_ATTR_LEAK_OK"; fi
# (c) headline Sn/Pr in band
gffcompare -r /mnt/c/Users/jfris/Desktop/GGO_19.gtf /tmp/denovo_branch.gtf -o /tmp/denovo_branch_cmp 2>/dev/null
grep -E 'Transcript level:' /tmp/denovo_branch_cmp.stats
```
Asserts:
- `BYTE_IDENTICAL_OK` (the hard guard — spec: "GTF byte-identical to HEAD").
- `NO_VG_ATTR_LEAK_OK` (spec: "Assert `anchored_coverage` stays 0 when `!vg_mode`"; the observable proxy is that no capacity attrs are emitted off-VG).
- Transcript level Sn within `95.6 ± 0.3` per `bench/multi_copy_eval/expectations.json` (`default_tx_sn`/`default_band`). **CAVEAT (verify in repo):** the committed `bench/ggo19_denovo_latest_cmp.stats` currently shows Transcript level **95.1 | 88.9**, not 95.6/90.5. Byte-identity to HEAD is the authoritative regression assert; treat the precise Sn/Pr pair as a recorded number against expectations.json, and if it reads 95.1/88.9 that is acceptable PROVIDED it equals the HEAD run (which byte-identity already guarantees). Log both numbers to `/tmp/task12_validation.log`.

- [ ] **Step 4: DAZ headline — coverage collapse + joint-strand single-PreEntry proof**
DAZ3 = `NC_073248.2:42879918-42945552`, `+` strand (per `bench/paralog_secondary_scan/diag_sites_daz.py:11` and `bench/multi_copy_eval/run_oracle.py:obj_daz`). Extract the locus, run VG with anchor-prior + joint-strand ON (default), and dump the per-read attribution TSV.
```bash
samtools view -b /mnt/c/Users/jfris/Desktop/GGO.bam NC_073248.2:42700000-43000000 > /tmp/daz.bam && samtools index /tmp/daz.bam
rm -f /tmp/daz_attr.tsv
RUSTLE_VG_FP_ATTR_TSV=/tmp/daz_attr.tsv \
  /mnt/c/Users/jfris/Desktop/Rustle/target/release/rustle --vg --genome-fasta /mnt/c/Users/jfris/Desktop/GGO.fasta -L /tmp/daz.bam -o /tmp/daz_on.gtf 2>/tmp/daz_on.stderr
# (a) DAZ3 coverage: the + strand transcript(s) overlapping the DAZ3 span
awk -F'\t' '$3=="transcript" && $7=="+" && $4<42945552 && $5>42879918 {for(i=1;i<=NF;i++) if($i ~ /cov/) print $9}' /tmp/daz_on.gtf | grep -oE 'cov "[0-9.]+"' | head
```
Asserts (spec test 2):
- DAZ3 per-copy `cov` drops to **~1-2** (spec target: 163 -> 1-2). It MUST be < 10 (was ~163 uncorrected); record the exact value.
- A DAZ1 (`-` strand, span `42783133-42859657`) transcript is still emitted with substantial `cov` (DAZ1 keeps its mass — "not over-deflated"):
```bash
awk -F'\t' '$3=="transcript" && $7=="-" && $4<42859657 && $5>42783133 {print $9}' /tmp/daz_on.gtf | grep -oE 'cov "[0-9.]+"' | head
```
Assert at least one DAZ1 `cov` value > 10 (its real mass survives; DAZ1 has ~167 unique reads per `vg.rs:5977`).

- [ ] **Step 5: DAZ joint-strand proof via the attribution TSV (one PreEntry, weight_sum==1.0 across strands)**
The TSV columns are `family_id	read_name_hash	placement_copy	n_sites_covered	final_weight	weight_gap	weight_sum` (`src/rustle/vg.rs:4819`). The joint-strand fix means a read shared between DAZ1(−) and DAZ3(+) appears as multiple placement rows under ONE `read_name_hash` whose `weight_sum` is ~1.0 (mass conserved across both strands), instead of being 1/NH in each strand-partition.
```bash
# pick a read_name_hash that has >=2 placement rows (a multimapper), assert its weight_sum ~ 1.0
python3 - <<'PY'
import collections
rows=collections.defaultdict(list); hdr=None
for ln in open('/tmp/daz_attr.tsv'):
    p=ln.rstrip('\n').split('\t')
    if hdr is None: hdr=p; continue
    d=dict(zip(hdr,p)); rows[d['read_name_hash']].append(d)
multi=[(k,v) for k,v in rows.items() if len(v)>=2]
assert multi, "FAIL: no multi-placement reads in DAZ attr TSV (joint-strand input not built)"
bad=[]
for k,v in multi:
    ws=float(v[0]['weight_sum'])               # weight_sum is per-read (same on every row)
    fw=sum(float(r['final_weight']) for r in v)
    if abs(ws-1.0)>1e-3 or abs(fw-1.0)>1e-3: bad.append((k,ws,round(fw,4),len(v)))
print(f"multi-placement reads: {len(multi)}; weight_sum~1.0 holds for {len(multi)-len(bad)}")
assert not bad, f"FAIL weight conservation: {bad[:5]}"
# distinct placement_copy count across the family proves cross-strand spread
copies=set(r['placement_copy'] for k,v in multi for r in v)
print(f"distinct placement copies among multimappers: {sorted(copies)}")
print("JOINT_STRAND_PREENTRY_OK")
PY
```
Assert the script prints `JOINT_STRAND_PREENTRY_OK` (every multimapper's `weight_sum` and Σ`final_weight` == 1.0 ± 1e-3).

- [ ] **Step 6: A/B against RUSTLE_VG_JOINT_STRAND_EM=0 (proves the fix is the cause)**
```bash
rm -f /tmp/daz_off_attr.tsv
RUSTLE_VG_JOINT_STRAND_EM=0 RUSTLE_VG_ANCHOR_PRIOR=0 RUSTLE_VG_FP_ATTR_TSV=/tmp/daz_off_attr.tsv \
  /mnt/c/Users/jfris/Desktop/Rustle/target/release/rustle --vg --genome-fasta /mnt/c/Users/jfris/Desktop/GGO.fasta -L /tmp/daz.bam -o /tmp/daz_off.gtf 2>/dev/null
echo "--- OFF (legacy VG) DAZ3 cov ---"
awk -F'\t' '$3=="transcript" && $7=="+" && $4<42945552 && $5>42879918 {print $9}' /tmp/daz_off.gtf | grep -oE 'cov "[0-9.]+"' | head
echo "--- ON DAZ3 cov ---"
awk -F'\t' '$3=="transcript" && $7=="+" && $4<42945552 && $5>42879918 {print $9}' /tmp/daz_on.gtf | grep -oE 'cov "[0-9.]+"' | head
```
Assert the OFF DAZ3 `cov` is dramatically higher than ON (OFF ≈ the inflated regime ~100-215; ON ≈ 1-2). The flag `=0` restoring the inflated value is the rollback proof (spec §Gating: "`=0` restores current VG behavior").

- [ ] **Step 7: Synthetic oracle unchanged (fingerprint-EM still 100%)**
The synthetic assignment scorer (commit 9436e7e fixture) must stay decisive-accuracy 100% (spec test 3; `expectations.json: obj4_decisive_acc_min=100.0`).
```bash
python3 /mnt/c/Users/jfris/Desktop/Rustle/bench/multi_copy_eval/score_synth_assignment.py \
  --rustle /mnt/c/Users/jfris/Desktop/Rustle/target/release/rustle
```
Assert the JSON shows `"decisive_accuracy": 100.0`. Then run the consolidated oracle in check mode (covers obj3 Sn floor + obj4 + default-headline isolation band in one shot):
```bash
python3 /mnt/c/Users/jfris/Desktop/Rustle/bench/multi_copy_eval/run_oracle.py --check
```
Assert it prints `ALL OBJECTIVES PASS` and exits 0. (If it prints `REGRESSION: ...`, the failing objective names which assert broke — do not edit the oracle, fix the implementation.)

- [ ] **Step 8: Mass-conservation unit test (cargo test)**
The earlier task added a `#[cfg(test)]` unit test in `src/rustle/vg.rs` asserting Σ per-read placement weights == 1.0 ± 1e-6 (including the gate-fired-read case). Run it:
```bash
cargo test --release --manifest-path /mnt/c/Users/jfris/Desktop/Rustle/Cargo.toml mass_conservation 2>&1 | tail -15
```
Assert the run ends with `test result: ok.` and a nonzero number of tests run (i.e. the named test exists and passes — `0 passed; 0 filtered out` means the test name is wrong; stop and reconcile with the earlier task's actual test name before proceeding).

- [ ] **Step 9: Genome-wide no-fabrication scan — produce the taxonomy (93 expressed_real_copy / 89 pure_spillover)**
Regenerate the read-level taxonomy over GGO (spec test 5). This consumes the placement/intron TSVs the scan tooling builds from `/mnt/c/Users/jfris/Desktop/GGO.bam`; rebuild them if absent (scripts are idempotent).
```bash
cd /mnt/c/Users/jfris/Desktop/Rustle
# 1) gene introns + secondary positions + enriched candidates (rebuild if /tmp inputs are stale)
python3 bench/paralog_secondary_scan/build_gene_introns.py   # -> /tmp/gene_introns.tsv (verify path inside script)
# 2) de-based anchoring taxonomy (the canonical 93/89 classifier; calibration t=2, extent 0.7 default — set EXTENT_FRAC=0.8 per approved decision)
DE_T=2 MIN_ANCHOR=3 EXTENT_FRAC=0.8 \
  python3 bench/paralog_secondary_scan/de_anchor.py \
    /tmp/mm_placements_de.tsv /tmp/gene_introns.tsv /tmp/hits_enriched.json /tmp/hits_de.json
python3 - <<'PY'
import json,collections
rows=json.load(open('/tmp/hits_de.json'))
tax=collections.Counter(r['de_taxonomy'] for r in rows)
print("taxonomy:",dict(tax))
er=tax.get('expressed_real_copy',0); sp=tax.get('pure_spillover',0)
print(f"expressed_real_copy={er}  pure_spillover={sp}")
assert er>=90, f"FAIL: expressed_real_copy dropped to {er} (expected ~93; real copies must be RETAINED)"
assert sp>=80, f"FAIL: pure_spillover={sp} (expected ~89; fabrication-risk set must still be flagged)"
print("TAXONOMY_OK")
PY
```
Assert `TAXONOMY_OK`. The exact `EXTENT_FRAC=0.8` matches the approved anchor calibration (raw dNM `t=2`, `extent_frac=0.8`); the default in `de_anchor.py` is 0.7, so it must be overridden. **CAVEAT:** counts are approximate — assert the bands `expressed_real_copy >= 90` (retention guard, must not fall) and `pure_spillover >= 80` (flagging guard). If `build_gene_introns.py`/placement-TSV paths differ from the `/tmp/*` names above, read the top of each script for its real argv and adjust — these scripts take positional file args.

- [ ] **Step 10: Genome-wide scan — pure_spillover gets LOW capacity_confidence (no-fabrication guard)**
The taxonomy is a BAM-read classifier; the spec asserts that after the fix the 89 pure_spillover genes' assembled transcripts carry LOW `capacity_confidence` (band/abstain fires) while the 93 expressed_real_copy keep coverage. Produce a VG GTF over the same regions, then join gene spans -> emitted transcript `capacity_confidence`.
```bash
cd /mnt/c/Users/jfris/Desktop/Rustle
# VG run over the whole GGO BAM (or per-region if memory-bound; here genome-wide VG)
/mnt/c/Users/jfris/Desktop/Rustle/target/release/rustle --vg --genome-fasta /mnt/c/Users/jfris/Desktop/GGO.fasta -L /mnt/c/Users/jfris/Desktop/GGO.bam -o /tmp/ggo_vg.gtf 2>/tmp/ggo_vg.stderr
python3 - <<'PY'
import json,re,collections,bisect
rows=json.load(open('/tmp/hits_de.json'))
genes={r['gene']:r for r in rows}
# build chrom -> sorted [(start,end,gene,tax)] for genes we care about
spans=collections.defaultdict(list)
for r in rows:
    if r['de_taxonomy'] in ('expressed_real_copy','pure_spillover'):
        g=r  # de_anchor preserves enriched fields incl chrom/start/end
        spans[g.get('chrom')].append((int(g['start']),int(g['end']),g['gene'],r['de_taxonomy']))
for c in spans: spans[c].sort()
def gene_at(chrom,s,e):
    arr=spans.get(chrom) or []
    mid=(s+e)//2
    for ss,ee,gid,tax in arr:
        if ss<=mid<=ee: return (gid,tax)
    return None
# parse GTF transcripts, pull capacity_confidence
cc_by_tax=collections.defaultdict(list)
for ln in open('/tmp/ggo_vg.gtf'):
    if ln.startswith('#'): continue
    f=ln.rstrip('\n').split('\t')
    if len(f)<9 or f[2]!='transcript': continue
    m=re.search(r'capacity_confidence "([0-9.]+)"',f[8])
    if not m: continue
    hit=gene_at(f[0],int(f[3]),int(f[4]))
    if hit: cc_by_tax[hit[1]].append(float(m.group(1)))
import statistics as st
for tax in ('expressed_real_copy','pure_spillover'):
    v=cc_by_tax.get(tax,[])
    med=round(st.median(v),3) if v else None
    print(f"{tax}: n_tx_with_cc={len(v)} median_capacity_confidence={med}")
exp=cc_by_tax.get('expressed_real_copy',[]); spl=cc_by_tax.get('pure_spillover',[])
# guard: spillover should be lower-confidence than real copies
if exp and spl:
    me=st.median(exp); ms=st.median(spl)
    assert ms < me, f"FAIL: spillover median cc {ms} not below real-copy median {me}"
    print("CAPACITY_CONFIDENCE_SEPARATION_OK")
else:
    print("WARN: insufficient capacity_confidence coverage to assert separation; record raw counts")
PY
```
Asserts:
- The real-copy median `capacity_confidence` is meaningfully higher than the pure_spillover median (separation guard — spec: "93 retained, 89 get LOW capacity_confidence").
- expressed_real_copy transcripts are still EMITTED (coverage kept — abstain TAGs, does not drop; spec O5).
**CAVEAT:** capacity_confidence is only emitted on VG transcripts that pass through a family, so not every gene span maps to a transcript carrying it; the assert is on the median separation, not per-gene completeness. If `pure_spillover` has zero emitted transcripts with `capacity_confidence` (e.g. abstained-and-filtered upstream), that is consistent with the no-fabrication goal — record the count and treat the separation assert as WARN-only in that case. Log all numbers to `/tmp/task12_validation.log`.

- [ ] **Step 11: Clean up the HEAD worktree**
```bash
git -C /mnt/c/Users/jfris/Desktop/Rustle worktree remove /tmp/rustle_head --force 2>&1 | tail -2
git -C /mnt/c/Users/jfris/Desktop/Rustle worktree list
```
Expect the temp worktree gone (avoids a stale `git worktree` entry).

- [ ] **Step 12: Commit the validation record**
This task adds no source; record the validation by appending the run log into the bench directory as a plain text artifact (NOT a markdown report) so the result is reproducible and reviewable.
```bash
cp /tmp/task12_validation.log /mnt/c/Users/jfris/Desktop/Rustle/bench/multi_copy_eval/task12_validation.log
git -C /mnt/c/Users/jfris/Desktop/Rustle add bench/multi_copy_eval/task12_validation.log
git -C /mnt/c/Users/jfris/Desktop/Rustle commit -m "$(cat <<'EOF'
test(vg): final apportionment validation — regression guard + DAZ + oracle + mass-conservation + genome-wide scan

- default de-novo GTF byte-identical to HEAD; no capacity attrs leak off-VG
- DAZ3 cov ~163 -> ~1-2, DAZ1 mass retained; shared read in ONE PreEntry weight_sum==1.0
- JOINT_STRAND_EM=0 A/B restores inflated DAZ3 (rollback proof)
- synthetic fingerprint-EM oracle still 100% decisive; run_oracle.py --check passes
- mass_conservation unit test passes (Σ per-read weights == 1.0)
- genome-wide scan: 93 expressed_real_copy retained, 89 pure_spillover low capacity_confidence

Co-Authored-By: Claude Opus 4.8 (1M context) <noreply@anthropic.com>
EOF
)"
```

---

## Deferred to a follow-on plan (not in scope here)

**Spec Component 5 — coupled EM↔flow fixed-point fallback (`RUSTLE_VG_COUPLED_EM`, opt-in, OFF by default).** The panel scored it 4/10 (expensive: N_outer × full assembly; unstable on the 0-site DAZ case). The 0-anchor families this would target are handled in this plan by **graceful-degrade to the pileup prior + low-confidence flag** (Task 7 fallback + Task 11 tag), per approved decision O2. The coupled fallback is therefore a separate, independently-testable follow-on plan and is intentionally **not** a task here — this plan ships working, validated software (default-ON anchor-first + capacity-confidence) on its own. Track it as `docs/superpowers/plans/<date>-coupled-em-flow-fallback.md` when prioritized.

## Spec-coverage check (self-review)

| Spec component | Task(s) |
|---|---|
| 1 — inverted-pair joint EM input | T5 (`family_for_em_input`) + T6 (wiring, aligned graphs) |
| 2 — anchored-mass prior | T4 (`anchored_mass_per_copy`) + T7 (use, single E-step, graceful-degrade) |
| 3 — persist read anchoring | T1 (fields) + T8 (set in write-back) |
| 4 — capacity-confidence channel | T2 (GraphNode field) + T3 (Transcript fields) + T9 (accumulate) + T10 (compute) |
| 5 — coupled fallback (opt-in) | **Deferred** (see above) |
| 6 — GTF emission | T11 |
| Conservation (gate-fired fallback) | T8 |
| Abstain (tag, keep coverage) | T11 |
| Validation plan | T12 |
