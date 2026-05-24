# VG Enhanced Coverage/Junction Sharing Design

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task.

**Goal:** Replace the conservative structural-prior borrow (25%/cap-4.0) with full per-copy-expected sharing, add a borrow-OFF gate for controlled benchmarking, and run a three-condition GOLGA8 benchmark to measure impact.

**Architecture:** Keep per-bundle assembly unchanged. Change the coverage and junction sharing *policy* inside the existing borrow mechanism — no coordinate remapping, no new graph structures. Add `RUSTLE_VG_NO_BORROW=1` to cleanly disable borrowing for the OFF condition.

**Tech Stack:** Rust (src/rustle/pipeline.rs, src/rustle/vg.rs), bash benchmark script, gffcompare.

---

## Background and Motivation

In `--vg` mode each gene-family copy assembles independently from its own reads plus a small borrow from siblings. For highly accurate long reads (IsoSeq / ONT Q20+) a single read can fully characterize a transcript. The 7 unresolved GOLGA8N/NL copies fail not because their sequences are indistinguishable but because each copy has too few reads for `path_extract` to emit a path. Shared exons (identical across all GOLGA8 paralogs) should accumulate full family coverage; the copy-specific junction — seen in even 1–2 long reads — then becomes the anchor signal.

The current structural prior (25% of sibling, capped at 4.0) was calibrated for short-read noise avoidance. For long reads where per-copy expected coverage is 10–20 reads this cap cuts signal by 4–5×. The enhanced policy removes the cap and uses the mathematically correct expected value.

---

## Design

### 1. Borrow-OFF Gate

**File:** `src/rustle/pipeline.rs`, around lines 10287–10291.

Add one env-var check before populating `bundle_borrow_cov` and `bundle_borrow_junctions`:

```rust
let no_borrow = std::env::var_os("RUSTLE_VG_NO_BORROW").is_some();
if !no_borrow {
    bundle_borrow_cov = crate::vg::build_bundle_borrow_coverage(&vg_families, &family_graphs);
    bundle_borrow_junctions = crate::vg::build_bundle_borrow_junctions(
        &vg_families, &family_graphs, &bundles,
    );
}
```

When `RUSTLE_VG_NO_BORROW=1` both maps stay empty, so neither the structural prior nor the junction injection fires downstream. This is the clean OFF condition.

**Scope:** `RUSTLE_VG_NO_BORROW` disables only the ExonClass-based coverage/junction borrow. It does NOT affect `RUSTLE_VG_FAMILY_BOOST` (primary-read boost, gated separately) or EM reweighting. Those remain active so the OFF condition isolates the structural borrow signal specifically.

### 2. Enhanced Coverage Sharing Policy

**File:** `src/rustle/pipeline.rs`, around lines 12580–12611 (inside the per-bundle assembly loop).

The block currently has two paths: `RUSTLE_VG_BORROW_FLOOR` (legacy blind floor) and the default structural prior (25%/cap-4.0). Add a third path selected by `RUSTLE_VG_BORROW_LEGACY=1`, and make the new enhanced formula the default:

```rust
if use_legacy_floor {
    // RUSTLE_VG_BORROW_FLOOR=1: original blind 5% floor (unchanged)
    if node.coverage == 0.0 && total_fam_cov > 0.0 {
        node.coverage = 1.0_f64.min(total_fam_cov * 0.05);
    }
} else if std::env::var_os("RUSTLE_VG_BORROW_LEGACY").is_some() {
    // RUSTLE_VG_BORROW_LEGACY=1: old 25%/cap-4.0 structural prior
    let expected = if n_copies_total > 0 { total_fam_cov / n_copies_total as f64 } else { total_fam_cov };
    if this_cov < expected * 0.5 && max_sibling_cov > 1.0 {
        let prior = (max_sibling_cov * 0.25).min(4.0);
        node.coverage = node.coverage.max(prior);
    }
} else {
    // DEFAULT (new): per-copy-expected floor for shared nodes
    let expected_per_copy = if n_copies_total > 0 {
        total_fam_cov / n_copies_total as f64
    } else {
        total_fam_cov
    };
    node.coverage = node.coverage.max(expected_per_copy);
}
```

This lets the benchmark run all three conditions (OFF / Legacy / Enhanced) with a single binary.

**Rationale:** `total_fam_cov / n_copies_total` is the unbiased expected per-copy contribution. It respects the copy's own reads (`.max()` keeps higher values) and adds nothing for copy-specific bubble nodes (`if copy_specific { continue; }` guard is already in place).

### 3. Enhanced Junction Sharing Policy

**File:** `src/rustle/pipeline.rs`, around lines 12469–12484 (inside the junction injection block).

Change the injected count from 25% (capped at 5.0) to 50% (no hard cap):

```rust
// OLD
let c = (sibling_count * 0.25).min(5.0);

// NEW
let c = sibling_count * 0.5;
```

**Rationale:** A shared junction seen in `N` reads across all siblings translates to an expected `N / n_copies` per copy. For a family of 13 copies with 20 family-wide supporting reads, the per-copy expected count is ~1.5 — and `0.5 * 20 = 10` is already conservative relative to "full sharing." The 5.0 hard cap was a noise floor for short reads; for long reads it cuts real signal.

---

## Benchmark Protocol

**Script:** `bench/vg_borrow_benchmark.sh`

Three conditions on the GOLGA8 locus (`NC_073240.2:23-36Mb`), using `GGO_19.bam`, `GGO_19.gtf`, and the gorilla genome FASTA at `/mnt/c/Users/jfris/Desktop/GGO.fasta`:

| Label | Env vars | Notes |
|-------|----------|-------|
| **OFF** | `RUSTLE_VG_NO_BORROW=1` | No ExonClass borrow at all |
| **Legacy** | `RUSTLE_VG_BORROW_LEGACY=1` | Old 25%/cap-4.0 structural prior |
| **Enhanced** | _(none)_ | New per-copy-expected floor (default) |

All three use the same binary compiled from this branch. All three use `--vg --genome-fasta <path> --family-manifest <manifest>`.

Each condition:
1. Run Rustle, output GTF to a temp file
2. Run `gffcompare -r GGO_19.gtf -o bench/vg_borrow_benchmark/<label>/cmp <output.gtf>`
3. Extract exact-match count from gffcompare `.stats` file
4. Print summary table

**Success criterion:** `Enhanced ≥ Current ≥ OFF` on exact matches. If `Enhanced ≈ OFF`, borrowing is not the limiting factor — this is equally useful information (redirects focus to EM or junction thresholds).

**Baseline to beat:** Current VG = 6 exact matches vs StringTie = 3.

---

## Testing

- **Regression:** GOLGA8 YAG benchmark must pass with all three conditions; existing 12-exact-match baseline must not regress under OFF or Enhanced conditions vs pre-change baseline.
- **Synthetic family test** (`test_data/synthetic_family/`): must pass with `RUSTLE_VG_NO_BORROW=1` set and unset — confirms EM pathway is unaffected.
- **No new unit tests** needed: the sharing formula is a one-line arithmetic change; the gate is a one-line env-var check. Coverage comes from the GOLGA8 benchmark.

---

## Files Touched

| File | Change |
|------|--------|
| `src/rustle/pipeline.rs` | Borrow-OFF gate (~line 10287); three-way coverage formula (~line 12597); enhanced junction count (~line 12472) |
| `bench/vg_borrow_benchmark.sh` | New benchmark script (three conditions, single binary) |

No other files need to change. `vg.rs` build functions (`build_bundle_borrow_coverage`, `build_bundle_borrow_junctions`) are unchanged — the policy change is entirely in how the caller uses their output.

## Environment Variable Reference

| Variable | Effect |
|----------|--------|
| `RUSTLE_VG_NO_BORROW=1` | Disable ExonClass coverage + junction borrow entirely (OFF condition) |
| `RUSTLE_VG_BORROW_LEGACY=1` | Use old 25%/cap-4.0 structural prior (Legacy condition) |
| `RUSTLE_VG_BORROW_FLOOR=1` | Use original blind 5% floor (pre-existing, unchanged) |
| `RUSTLE_VG_FAMILY_BOOST=1` | Primary-read boost for underpowered LOC-prefix copies (separate mechanism, unaffected) |
