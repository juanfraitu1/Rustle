# ExonClass Two-Threshold Borrowing Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Split the single `min_jaccard` parameter in `build_family_graph` into separate merge and refine thresholds so cross-copy ExonClasses survive into the assembly loop, activating junction propagation and dark-exon completion for sparse VG-mode paralogs.

**Architecture:** `build_family_graph` currently uses one threshold (0.30) for both `merge_singletons_by_sequence` (needs to be low: ~0.05 for GOLGA8-like families) and `refine_by_minimizer_jaccard` (needs to be high: 0.30, but 0.0 = no-op in non-HMM path where HMM profiles are not fitted). The main-EM call site in `pipeline.rs` is updated to pass `merge=0.05, refine=0.0` (non-HMM) or `merge=0.05, refine=0.30` (HMM). Filter/scoring call sites keep `(0.30, 0.30)` — they only need FG structure, not borrowing.

**Tech Stack:** Rust, rayon (existing); no new dependencies.

---

## Root Cause Context

The YAG (chrY GOLGA8) benchmark shows `max_sib=0.00` for every ExonClass in `RUSTLE_VG_COMPLETION_TRACE` output, meaning:
- `merge_singletons_by_sequence` with threshold 0.30 fails to merge singleton exons from different copies (GOLGA8 family k-mer Jaccard is 0.071–0.183, below 0.30).
- `refine_by_minimizer_jaccard` with threshold 0.30 splits the 3 position-clustered multi-copy ExonClasses back into per-copy singletons.
- Net: 12 copy-specific ExonClasses, 0 shared → zero synthetic injections.

`refine_by_minimizer_jaccard` with `min_jaccard=0.0` already behaves as a no-op (all pairs have Jaccard ≥ 0.0, so all stay merged). No change to its body is needed.

---

## File Map

| File | Change |
|------|--------|
| `src/rustle/vg_hmm/family_graph.rs` | Replace `min_jaccard: f64` with `merge_min_jaccard: f64, refine_min_jaccard: f64` in `build_family_graph`; thread to internal calls |
| `src/rustle/pipeline.rs` (line ~10212) | Pass `0.05, if do_hmm { 0.30 } else { 0.0 }` — main EM setup, borrowing-critical |
| `src/rustle/pipeline.rs` (line ~145) | Pass `0.30, 0.30` — kmer-Jaccard family filter, unchanged |
| `src/rustle/vg.rs` (line ~459) | Pass `0.30, 0.30` — kmer scoring, unchanged |
| `src/rustle/vg.rs` (line ~623) | Pass `0.30, 0.30` — POA scoring, unchanged |

---

## Task 1: Add unit test for the existing `refine_by_minimizer_jaccard` no-op behavior

**Files:**
- Modify: `src/rustle/vg_hmm/family_graph.rs`

- [ ] **Step 1: Write a failing test (it won't compile yet — signature change is in Task 2)**

Add this test module at the bottom of `src/rustle/vg_hmm/family_graph.rs`:

```rust
#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn refine_zero_threshold_keeps_all_merged() {
        // Two sequences with ~10% Jaccard similarity (representative of GOLGA8-level divergence).
        // A 200-nt sequence and a shifted version: most minimizers differ.
        let seq_a: Vec<u8> = b"ATCGATCGATCGATCGATCGATCGATCGATCG"
            .iter().cycle().take(200).copied().collect();
        let seq_b: Vec<u8> = b"GCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTA"
            .iter().cycle().take(200).copied().collect();

        let cluster = vec![(0usize, seq_a.clone()), (1usize, seq_b.clone())];

        // With threshold=0.0, refine is a no-op — all stay in one group.
        let groups_zero = refine_by_minimizer_jaccard(&cluster, 0.0, 15, 10);
        assert_eq!(groups_zero.len(), 1,
            "threshold=0.0 should keep all sequences in one ExonClass");
        assert_eq!(groups_zero[0].len(), 2);

        // With threshold=0.30, divergent sequences are split.
        let groups_high = refine_by_minimizer_jaccard(&cluster, 0.30, 15, 10);
        assert_eq!(groups_high.len(), 2,
            "threshold=0.30 should split GOLGA8-divergence sequences");
    }

    #[test]
    fn refine_identical_sequences_always_merged() {
        let seq: Vec<u8> = b"ATCGATCGATCGATCGATCG".iter().cycle().take(200).copied().collect();
        let cluster = vec![(0usize, seq.clone()), (1usize, seq.clone())];

        for threshold in [0.0, 0.05, 0.30, 0.99] {
            let groups = refine_by_minimizer_jaccard(&cluster, threshold, 15, 10);
            assert_eq!(groups.len(), 1,
                "identical sequences should always stay merged (threshold={threshold})");
        }
    }
}
```

- [ ] **Step 2: Run to confirm it compiles and passes with current code**

```bash
cd /mnt/c/Users/jfris/Desktop/Rustle && cargo test -p rustle refine_zero_threshold 2>&1 | tail -10
```

Expected: `test vg_hmm::family_graph::tests::refine_zero_threshold_keeps_all_merged ... ok`

The test should already pass — the no-op behavior at threshold=0.0 is a property of the existing union-find logic, not a new feature. This confirms the invariant before we change the signature.

---

## Task 2: Split `min_jaccard` into `merge_min_jaccard` and `refine_min_jaccard` in `build_family_graph`

**Files:**
- Modify: `src/rustle/vg_hmm/family_graph.rs` (lines 410–416, 456–469, 486–488)

- [ ] **Step 1: Update the function signature**

In `src/rustle/vg_hmm/family_graph.rs`, change:

```rust
pub fn build_family_graph(
    family: &FamilyGroup,
    bundles: &[Bundle],
    genome: Option<&GenomeIndex>,
    min_pos_recip: f64,
    min_jaccard: f64,
) -> Result<FamilyGraph> {
```

to:

```rust
pub fn build_family_graph(
    family: &FamilyGroup,
    bundles: &[Bundle],
    genome: Option<&GenomeIndex>,
    min_pos_recip: f64,
    merge_min_jaccard: f64,
    refine_min_jaccard: f64,
) -> Result<FamilyGraph> {
```

- [ ] **Step 2: Update the `merge_singletons_by_sequence` call**

Inside `build_family_graph`, find the call to `merge_singletons_by_sequence` (around line 458) and change:

```rust
let merged = merge_singletons_by_sequence(pos_clusters, &copies, g, min_jaccard);
```

to:

```rust
let merged = merge_singletons_by_sequence(pos_clusters, &copies, g, merge_min_jaccard);
```

- [ ] **Step 3: Update the `refine_by_minimizer_jaccard` call**

Inside the `for cluster in &effective_clusters` loop (around line 488), change:

```rust
refine_by_minimizer_jaccard(&with_seq, min_jaccard, 15, 10)
```

to:

```rust
refine_by_minimizer_jaccard(&with_seq, refine_min_jaccard, 15, 10)
```

- [ ] **Step 4: Verify compilation fails at call sites (expected)**

```bash
cd /mnt/c/Users/jfris/Desktop/Rustle && cargo check 2>&1 | grep "build_family_graph"
```

Expected: 4 errors about wrong number of arguments — one per call site. This confirms the signature change propagated.

---

## Task 3: Update the four call sites

**Files:**
- Modify: `src/rustle/pipeline.rs` (two call sites)
- Modify: `src/rustle/vg.rs` (two call sites)

- [ ] **Step 1: Update the main-EM call site in `pipeline.rs` (~line 10212)**

Find:
```rust
match crate::vg_hmm::family_graph::build_family_graph(
    fam, &bundles, genome_ref, 0.30, 0.30,
)
```

Replace with:
```rust
match crate::vg_hmm::family_graph::build_family_graph(
    fam, &bundles, genome_ref, 0.30, 0.05, if do_hmm { 0.30 } else { 0.0 },
)
```

This is the borrowing-critical path: low merge threshold (0.05) allows GOLGA8-level paralogs to share ExonClasses; zero refine threshold in non-HMM mode keeps them shared (HMM mode keeps 0.30 for profile accuracy).

- [ ] **Step 2: Update the kmer-Jaccard filter call site in `pipeline.rs` (~line 145)**

Find:
```rust
let fg = match crate::vg_hmm::family_graph::build_family_graph(
    &largest, bundles, Some(g), 0.30, 0.30,
) { Ok(g) => g, Err(_) => return None };
```

Replace with:
```rust
let fg = match crate::vg_hmm::family_graph::build_family_graph(
    &largest, bundles, Some(g), 0.30, 0.30, 0.30,
) { Ok(g) => g, Err(_) => return None };
```

This is only used to compute k-mer Jaccard scores for family filtering — borrowing is irrelevant, keep existing thresholds.

- [ ] **Step 3: Update the kmer-scoring call site in `vg.rs` (~line 459)**

Find:
```rust
let fg = match crate::vg_hmm::family_graph::build_family_graph(
    &largest, bundles, Some(genome), 0.30, 0.30,
) {
```

Replace with:
```rust
let fg = match crate::vg_hmm::family_graph::build_family_graph(
    &largest, bundles, Some(genome), 0.30, 0.30, 0.30,
) {
```

- [ ] **Step 4: Update the POA-scoring call site in `vg.rs` (~line 623)**

Find:
```rust
let fg = match crate::vg_hmm::family_graph::build_family_graph(
    &largest, bundles, Some(genome), 0.30, 0.30,
) {
```

Replace with:
```rust
let fg = match crate::vg_hmm::family_graph::build_family_graph(
    &largest, bundles, Some(genome), 0.30, 0.30, 0.30,
) {
```

- [ ] **Step 5: Verify all call sites compile**

```bash
cd /mnt/c/Users/jfris/Desktop/Rustle && cargo check 2>&1 | grep -E "error|warning.*build_family_graph"
```

Expected: no errors.

---

## Task 4: Build and run existing tests

**Files:** none

- [ ] **Step 1: Build release binary**

```bash
cd /mnt/c/Users/jfris/Desktop/Rustle && cargo build --release 2>&1 | tail -5
```

Expected: `Finished release [optimized] target(s)`.

- [ ] **Step 2: Run unit tests**

```bash
cd /mnt/c/Users/jfris/Desktop/Rustle && cargo test --release 2>&1 | tail -20
```

Expected: all tests pass. The two new tests from Task 1 should appear as `ok`.

---

## Task 5: Verify borrowing is now active on the YAG benchmark

**Files:** none (diagnostic run only)

- [ ] **Step 1: Run YAG with completion trace**

```bash
cd /mnt/c/Users/jfris/Desktop/Rustle && \
RUSTLE_VG_COMPLETION_TRACE=1 ./target/release/rustle \
  bench/multi_copy_eval/ggo_Y.bam \
  --genome-fasta /mnt/c/Users/jfris/Desktop/GGO.fasta \
  --vg --vg-solver em --vg-snp --vg-no-hmm \
  -o /tmp/yag_borrow_v2.gtf 2>/tmp/yag_borrow_v2.stderr
```

- [ ] **Step 2: Check that multi-copy ExonClasses now appear**

```bash
grep -E '^\[FG\].*multi-copy' /tmp/yag_borrow_v2.stderr
```

Expected: lines like `[FG]   after merge: clusters=N multi-copy=M` where M > 0, and the subsequent `[VG] Graph completion: X dark exon(s) across Y bundle(s)` or `[VG] Junction propagation: X synthetic junction(s)` lines appear.

- [ ] **Step 3: Check that `max_sib` is no longer always 0.00**

```bash
grep '\[VG-COMPL\]' /tmp/yag_borrow_v2.stderr | grep -v 'max_sib=0.00' | head -5
```

Expected: at least some lines with `max_sib > 0.00`.

- [ ] **Step 4: Check transcript count for no regression**

```bash
grep -c $'\ttranscript\t' /tmp/yag_borrow_v2.gtf
```

Expected: 327 (unchanged from stored baseline). Graph completion and junction propagation add structure but should not drop existing transcripts since they only add evidence, not filter it.

- [ ] **Step 5: Commit**

```bash
cd /mnt/c/Users/jfris/Desktop/Rustle && \
git add src/rustle/vg_hmm/family_graph.rs src/rustle/pipeline.rs src/rustle/vg.rs && \
git commit -m "$(cat <<'EOF'
vg: split ExonClass Jaccard into separate merge/refine thresholds

build_family_graph now takes merge_min_jaccard (default 0.05 for the
borrowing-critical EM path) and refine_min_jaccard (0.0 in non-HMM mode
to skip splitting; 0.30 in HMM mode for profile accuracy).

GOLGA8 paralogs have family-level k-mer Jaccard of 0.071–0.183. The
previous hardcoded 0.30 caused merge_singletons_by_sequence to fail
(no singleton merges) and refine_by_minimizer_jaccard to split all
position-clustered multi-copy ExonClasses back into per-copy singletons.
Result: max_sib=0 everywhere → zero synthetic injections.

With merge=0.05 / refine=0.0, shared ExonClasses survive into the
assembly loop, activating junction propagation and dark-exon completion
for under-covered VG-mode paralogs.

Co-Authored-By: Claude Sonnet 4.6 <noreply@anthropic.com>
EOF
)"
```

---

## Self-Review

**Spec coverage:**
- Two-threshold split: covered in Tasks 2–3
- merge_min_jaccard=0.05 for main-EM: covered in Task 3 Step 1
- refine_min_jaccard=0.0 for non-HMM: covered in Task 3 Step 1
- refine_min_jaccard=0.30 for HMM: covered in Task 3 Step 1
- Filter/scoring call sites keep 0.30: covered in Task 3 Steps 2–4
- Borrowing verification on YAG: covered in Task 5

**Placeholder scan:** None found.

**Type consistency:** `merge_min_jaccard: f64` and `refine_min_jaccard: f64` used consistently across all tasks. `refine_by_minimizer_jaccard` function signature unchanged (only the call site argument changes).
