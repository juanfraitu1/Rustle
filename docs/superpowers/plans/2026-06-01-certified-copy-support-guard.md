# Certified Copy-Support Guard — Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: superpowers:subagent-driven-development. Steps use `- [ ]`.

**Goal:** Suppress phantom VG copies (a copy whose reads all fit a sibling better, e.g. DAZ3) and emit `copy_independent_support` as a calibrated certificate on every VG transcript — without touching the default de-novo path.

**Architecture:** A pure metric fn (unit-tested) computes per-copy independent support from cross-locus NM; a thin extractor pulls the data from the family/bundles; a post-assembly pipeline pass suppresses below TAU and annotates survivors.

**Spec:** `docs/superpowers/specs/2026-06-01-certified-copy-support-guard-design.md`.

**Conventions:** build `cargo build --release 2>&1 | tail -3` (ends "Finished"); commits end with `Co-Authored-By: Claude Opus 4.8 (1M context) <noreply@anthropic.com>`; never stage `tools/stringtie`; whole-suite test uses `cargo test --release --lib` (external test files have a pre-existing `junction_pair_stats` compile break unrelated to this work).

---

### Task 1: Pure support-fraction metric (TDD)
**Files:** `src/rustle/vg.rs` (add `copy_support_fraction`); in-module `#[cfg(test)]` tests.

- [ ] **Step 1: failing tests.** Add an in-module test module. The contract:
```rust
// rate pairs are (rate_C, rate_min_sibling); n_unique reads always support C.
// A multimapper supports C iff rate_C <= rate_min_sibling + margin.
#[test] fn daz3_phantom_zero_support() {
    // 0 unique; all multimappers fit a sibling far better (rate_C 0.07 vs sib 0.005)
    let mm = vec![(0.07, 0.005); 30];
    assert!(copy_support_fraction(0, &mm, 0.01) < 0.05);
}
#[test] fn daz1_real_high_support() {
    // many unique reads + multimappers that fit C best
    let mm = vec![(0.005, 0.07); 14];
    assert!(copy_support_fraction(167, &mm, 0.01) > 0.95);
}
#[test] fn co_expressed_tie_supports_both() {
    // NM-tie multimappers (within margin) count as supporting C (no false suppression)
    let mm = vec![(0.004, 0.004); 14];
    assert!(copy_support_fraction(9, &mm, 0.01) > 0.95);
}
#[test] fn margin_boundary() {
    assert_eq!(copy_support_fraction(0, &vec![(0.015, 0.004)], 0.01), 0.0); // 0.015 > 0.004+0.01 -> belongs elsewhere
    assert_eq!(copy_support_fraction(0, &vec![(0.013, 0.004)], 0.01), 1.0); // within margin -> supports
}
#[test] fn no_reads_is_zero() { assert_eq!(copy_support_fraction(0, &[], 0.01), 0.0); }
```
- [ ] **Step 2:** run → fail (`copy_support_fraction` undefined).
- [ ] **Step 3: implement** (pure, ~10 lines):
```rust
/// Fraction of a copy's reads that fit IT at least as well as any sibling.
/// n_unique reads always support C; a multimapper (rate_C, rate_min_sibling)
/// supports C iff rate_C <= rate_min_sibling + margin. Returns 0.0 if no reads.
pub fn copy_support_fraction(n_unique: usize, multimappers: &[(f64, f64)], margin: f64) -> f64 {
    let total = n_unique + multimappers.len();
    if total == 0 { return 0.0; }
    let mm_support = multimappers.iter().filter(|&&(rc, rs)| rc <= rs + margin).count();
    (n_unique + mm_support) as f64 / total as f64
}
```
- [ ] **Step 4:** run → all pass. `cargo build --release` Finished.
- [ ] **Step 5: commit** `vg: copy_support_fraction — per-copy independent-support metric (pure, tested)`.

### Task 2: Extractor + pipeline guard wiring
**Files:** `src/rustle/vg.rs` (add `compute_copy_independent_support`), `src/rustle/pipeline.rs` (call + suppress + annotate), `src/rustle/path_extract.rs`/types + `src/rustle/gtf.rs` (carry + emit `copy_independent_support`).

- [ ] **Step 1: `compute_copy_independent_support(family, bundles, margin) -> HashMap<usize,f64>`.** For each copy (fam_pos) in `family.bundle_indices`: count C-unique reads (in C's bundle, read_name_hash NOT in `family.multimap_reads`); for each C-multimapper (entries of `family.multimap_reads` with a placement at C), compute `rate_C = nm/aligned_len` at C and `rate_min_sib = min over the read's OTHER placements`. Aligned length: use the read's exon span sum (or `query_length`). Missing `nm` → treat read as supporting C (no sibling evidence; never penalize for missing data). Return `copy -> copy_support_fraction(...)`.
- [ ] **Step 2:** verify against a hand-built fixture if Bundle construction is feasible; otherwise rely on Task 1's unit tests for the metric and validate end-to-end in Task 3 (note which).
- [ ] **Step 3: wire into the pipeline** after transcripts carry `vg_family_id`/`vg_copy_id`, before GTF write. Group VG transcripts by `(vg_family_id, vg_copy_id)`; per family call `compute_copy_independent_support`; **drop** transcripts whose copy support `< TAU`; set `copy_independent_support` on survivors. Read `TAU` from `RUSTLE_VG_MIN_INDEP_SUPPORT` (default 0.10), `margin` from `RUSTLE_VG_SUPPORT_NM_MARGIN` (default 0.01). Log `[VG] copy-support guard: suppressed N copies (independent_support < TAU)`.
- [ ] **Step 4: emit the attribute.** Add `copy_independent_support: Option<f64>` to the transcript struct (or reuse the existing copy-confidence carrier), and emit `copy_independent_support "<v>"` in `gtf.rs` for VG transcripts. Non-VG transcripts: no attribute (untouched).
- [ ] **Step 5:** `cargo build --release` Finished; `cargo test --release --lib` green.
- [ ] **Step 6: commit** `vg: certified copy-support guard — suppress phantom copies + emit copy_independent_support`.

### Task 3: End-to-end validation gate (controller-evaluated)
- [ ] **DAZ default `--vg`** (`rustle --vg --genome-fasta ../GGO.fasta -L /tmp/daz.bam`): DAZ3 (+strand 42.879-42.946M) **suppressed** (0 tx, or independent_support≈0 if flagged); DAZ1 (−) **kept**. Record numbers (was DAZ3=5tx/cov163).
- [ ] **fam 175** (`/tmp/f175.bam`, `/tmp/chr234.fa`): **both** copies kept; covB>covA preserved; each `copy_independent_support` high.
- [ ] **fam 214** (`/tmp/f214.bam`, `/tmp/chr224.fa`): both copies kept.
- [ ] **synthetic oracle**: `rm -f /tmp/synth_assign.attr.tsv; python3 bench/multi_copy_eval/run_oracle.py --fast --check` → ALL OBJECTIVES PASS (Obj-4 100%; both synthetic copies kept).
- [ ] **no-false-suppression**: confirm no genuine copy dropped in fam175/214/synthetic.
- [ ] **default de-novo isolation**: a non-`--vg` GGO_19 run → output unchanged (guard is a no-op without `vg_copy_id`; the change is VG-gated). Confirm.
- [ ] If DAZ3 is NOT suppressed, trace `independent_support(DAZ3)` — it must be ≈0; if not, the extractor's cross-locus NM lookup for cross-strand siblings is wrong (DAZ1/DAZ3 must be linked in `family.multimap_reads`). Fix before proceeding.

### Task 4: Oracle check + expectations
**Files:** `bench/multi_copy_eval/score_copy_support.py` (new), `expectations.json`, `run_oracle.py`.
- [ ] Write `score_copy_support.py`: run DAZ + fam175 + fam214, assert DAZ3 suppressed (0 tx in the +strand window) and fam175/214 both copies present; emit JSON.
- [ ] Add `copy_guard_daz3_suppressed: true`, `copy_guard_fam175_both_kept: true` to `expectations.json`; wire into `run_oracle.py --full --check` (SKIPPED when GGO data absent).
- [ ] **commit** `oracle: copy-support guard check (DAZ3 suppressed, fam175/214 kept)`.

### Task 5: Docs + final review
- [ ] Update `docs/experiments/DAZ_vg_verification.md`: DAZ3 is now correctly suppressed (the false positive is gone on the default path); record the before/after.
- [ ] Final code-review pass over the diff; confirm default de-novo unchanged + the gate green.

---

## Self-Review
- **Spec coverage:** §3 metric → Task 1 + Task 2 Step 1; §3 suppress+certificate → Task 2 Steps 3–4; §4 location → Task 2 Step 3; §5 gate → Task 3 + Task 4; §6 isolation → Task 3 default-de-novo step. Covered.
- **Placeholders:** the metric fn is complete code; the extractor + wiring give exact inputs/fields (the implementer locates the precise emission point — acceptable, it's integration).
- **Type consistency:** `copy_support_fraction(usize, &[(f64,f64)], f64) -> f64` used identically in Task 1 tests and Task 2 extractor; `copy_independent_support: Option<f64>` consistent across struct + gtf.
