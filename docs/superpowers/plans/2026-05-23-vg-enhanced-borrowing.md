# VG Enhanced Coverage/Junction Sharing Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Replace the conservative per-copy structural-prior borrow with a full per-copy-expected floor, add a borrow-OFF gate for controlled benchmarking, and run a three-condition GOLGA8 benchmark to measure impact on multi-copy transcript recovery.

**Architecture:** Three targeted edits to `src/rustle/pipeline.rs` (borrow-OFF gate, three-way coverage policy, junction injection formula) plus a new benchmark shell script. No new files, no new structs, no coordinate remapping — only policy changes inside the existing VG borrow mechanism.

**Tech Stack:** Rust (cargo), bash, gffcompare 0.12.10, samtools.

---

## File Structure

| File | Role |
|------|------|
| `src/rustle/pipeline.rs` | Three edits (see tasks 1–3) |
| `bench/vg_borrow_benchmark.sh` | New benchmark script — three conditions, single binary |

---

## Background (read before touching code)

In `--vg` mode Rustle assembles each gene-family copy independently. After EM reweighting, two "borrow" maps are built from the FamilyGraph:

- `bundle_borrow_cov`: per-bundle list of `(ec_start, ec_end, copy_specific, this_cov, total_fam_cov, max_sibling_cov, n_copies_total)` — used to boost shared-exon node coverage in the per-bundle graph.
- `bundle_borrow_junctions`: per-bundle list of `(donor, acceptor, strand, sibling_count)` — used to inject synthetic junction stats from well-covered siblings.

These maps are built once (around `pipeline.rs:10287`) and consumed in the per-bundle assembly loop (~`pipeline.rs:12466` and `12580`).

**Current coverage policy** (the "Legacy" condition in the benchmark):
```rust
if this_cov < expected * 0.5 && max_sibling_cov > 1.0 {
    let prior = (max_sibling_cov * 0.25).min(4.0);
    node.coverage = node.coverage.max(prior);
}
```
For long reads with 10–20 reads per copy the cap of 4.0 is 3–5× too conservative.

**Enhanced policy** (the new default):
```rust
let expected_per_copy = total_fam_cov / n_copies_total as f64;
node.coverage = node.coverage.max(expected_per_copy);
```
Every shared-exon node gets at least the per-copy-expected coverage from family evidence.

**Junction policy** (same for both Legacy and Enhanced):
Change from `(sibling_count * 0.25).min(5.0)` to `sibling_count * 0.5`.

---

## Task 1: Add borrow-OFF gate (`RUSTLE_VG_NO_BORROW=1`)

**Files:**
- Modify: `src/rustle/pipeline.rs:10286-10302`

- [ ] **Step 1: Locate the block**

Open `src/rustle/pipeline.rs` and find line 10286:
```rust
                    if build_graph {
                        bundle_borrow_cov = crate::vg::build_bundle_borrow_coverage(
                            &em_hmm_partitions,
                            &family_graphs,
                        );
                        bundle_borrow_junctions = crate::vg::build_bundle_borrow_junctions(
                            &em_hmm_partitions,
                            &family_graphs,
                            &bundles,
                        );
                        if std::env::var_os("RUSTLE_VG_COMPLETION_OFF").is_none() {
                            bundle_completion_nodes = crate::vg::build_bundle_completion_nodes(
                                &em_hmm_partitions,
                                &family_graphs,
                            );
                        }
                    }
```

- [ ] **Step 2: Replace with gated version**

Replace the entire `if build_graph { ... }` block with:
```rust
                    if build_graph {
                        if std::env::var_os("RUSTLE_VG_NO_BORROW").is_none() {
                            bundle_borrow_cov = crate::vg::build_bundle_borrow_coverage(
                                &em_hmm_partitions,
                                &family_graphs,
                            );
                            bundle_borrow_junctions = crate::vg::build_bundle_borrow_junctions(
                                &em_hmm_partitions,
                                &family_graphs,
                                &bundles,
                            );
                        }
                        if std::env::var_os("RUSTLE_VG_COMPLETION_OFF").is_none() {
                            bundle_completion_nodes = crate::vg::build_bundle_completion_nodes(
                                &em_hmm_partitions,
                                &family_graphs,
                            );
                        }
                    }
```

`bundle_completion_nodes` is intentionally NOT gated — it uses a separate env var (`RUSTLE_VG_COMPLETION_OFF`) and is not part of the coverage/junction borrow being benchmarked.

- [ ] **Step 3: Verify it compiles**

```bash
cargo build --release 2>&1 | tail -5
```
Expected: no errors, warnings only (if any).

---

## Task 2: Three-way coverage sharing policy

**Files:**
- Modify: `src/rustle/pipeline.rs:12590-12608`

- [ ] **Step 1: Locate the block**

Find lines 12590–12608 in `pipeline.rs`. The `if use_legacy_floor { ... } else { ... }` block looks like this:

```rust
                                if use_legacy_floor {
                                    // Old behaviour: blind 5% floor for dark nodes
                                    if node.coverage == 0.0 && total_fam_cov > 0.0 {
                                        node.coverage = 1.0_f64.min(total_fam_cov * 0.05);
                                    }
                                } else {
                                    // Structural prior: boost under-represented copies
                                    let expected = if n_copies_total > 0 {
                                        total_fam_cov / n_copies_total as f64
                                    } else {
                                        total_fam_cov
                                    };
                                    if this_cov < expected * 0.5
                                        && max_sibling_cov > 1.0
                                    {
                                        let prior = (max_sibling_cov * 0.25).min(4.0);
                                        node.coverage = node.coverage.max(prior);
                                    }
                                }
```

- [ ] **Step 2: Replace the `else` branch with three-way logic**

Replace the entire `if use_legacy_floor { ... } else { ... }` with:

```rust
                                if use_legacy_floor {
                                    // RUSTLE_VG_BORROW_FLOOR=1: original blind 5% floor
                                    if node.coverage == 0.0 && total_fam_cov > 0.0 {
                                        node.coverage = 1.0_f64.min(total_fam_cov * 0.05);
                                    }
                                } else if std::env::var_os("RUSTLE_VG_BORROW_LEGACY").is_some() {
                                    // RUSTLE_VG_BORROW_LEGACY=1: old 25%/cap-4.0 structural prior
                                    let expected = if n_copies_total > 0 {
                                        total_fam_cov / n_copies_total as f64
                                    } else {
                                        total_fam_cov
                                    };
                                    if this_cov < expected * 0.5 && max_sibling_cov > 1.0 {
                                        let prior = (max_sibling_cov * 0.25).min(4.0);
                                        node.coverage = node.coverage.max(prior);
                                    }
                                } else {
                                    // Enhanced (default): per-copy-expected floor for shared nodes
                                    let expected_per_copy = if n_copies_total > 0 {
                                        total_fam_cov / n_copies_total as f64
                                    } else {
                                        total_fam_cov
                                    };
                                    node.coverage = node.coverage.max(expected_per_copy);
                                }
```

- [ ] **Step 3: Update the comment block immediately above the borrow block**

Find and replace this exact comment block (lines ~12565–12578):

```rust
                    // Structural graph merging: boost shared-exon nodes whose
                    // coverage in this copy is far below what the family signal
                    // predicts. Uses ExonClass equivalences (structural signal)
                    // rather than per-read sequence features (SNPs). Only shared
                    // exons are eligible; copy-specific bubble branches are skipped.
                    //
                    // Condition: this copy's ExonClass-level coverage is below
                    //   50% of (total_fam_cov / n_copies_total) AND at least one
                    //   sibling has meaningful coverage (max_sibling_cov > 1.0).
                    // Amount: 25% of the best sibling's coverage, capped at 4.0.
                    //   This is enough to let path_extract emit paths without
                    //   inflating TPM to ST-comparable levels.
                    //
                    // Legacy crude floor (RUSTLE_VG_BORROW_FLOOR=1) still works
                    // as a fallback; structural prior is the default path.
```

Replace with:

```rust
                    // Coverage borrowing from sibling copies. Three modes:
                    //   RUSTLE_VG_BORROW_FLOOR=1  — blind 5% floor (pre-existing)
                    //   RUSTLE_VG_BORROW_LEGACY=1 — old 25%/cap-4.0 structural prior
                    //   default                   — per-copy-expected floor
                    //     (total_fam_cov / n_copies_total), unconditional for shared nodes
                    // Only non-copy-specific ExonClass nodes are eligible.
```

- [ ] **Step 4: Verify it compiles**

```bash
cargo build --release 2>&1 | tail -5
```
Expected: no errors.

---

## Task 3: Update junction injection formula

**Files:**
- Modify: `src/rustle/pipeline.rs:12472`

- [ ] **Step 1: Locate the line**

Find line 12472 inside the `for &(donor, acceptor, strand, sibling_count) in borrow_jcts` loop:

```rust
                                let c = (sibling_count * 0.25).min(5.0);
```

- [ ] **Step 2: Replace**

```rust
                                let c = sibling_count * 0.5;
```

No other change in this block.

- [ ] **Step 3: Compile and run tests**

```bash
cargo build --release 2>&1 | tail -3
cargo test 2>&1 | tail -20
```

Expected: build succeeds, all tests pass (including the synthetic family VG test).

- [ ] **Step 4: Commit**

```bash
git add src/rustle/pipeline.rs
git commit -m "vg: enhanced coverage/junction sharing policy + borrow-OFF gate

- RUSTLE_VG_NO_BORROW=1: skip bundle_borrow_cov/junctions population (OFF condition)
- RUSTLE_VG_BORROW_LEGACY=1: restore old 25%/cap-4.0 structural prior (Legacy condition)
- Default enhanced: per-copy-expected floor (total_fam_cov / n_copies_total)
- Junction injection: 0.25*min(5) → 0.5 (no cap)

Co-Authored-By: Claude Sonnet 4.6 <noreply@anthropic.com>"
```

---

## Task 4: Write the benchmark script

**Files:**
- Create: `bench/vg_borrow_benchmark.sh`

- [ ] **Step 1: Create the script**

```bash
cat > bench/vg_borrow_benchmark.sh << 'SCRIPT'
#!/usr/bin/env bash
# Three-condition VG borrow benchmark on GOLGA8 region.
# Usage: bash bench/vg_borrow_benchmark.sh
# Requires: target/release/rustle built, gffcompare on PATH
set -euo pipefail

RUSTLE=target/release/rustle
BAM=bench/multi_copy_eval/golga8_region.bam
REF=bench/multi_copy_eval/ref_GOLGA8.gff
FASTA=/mnt/c/Users/jfris/Desktop/GGO.fasta
OUTDIR=bench/vg_borrow_benchmark

mkdir -p "$OUTDIR"/{off,legacy,enhanced}

run_condition() {
    local label=$1
    local extra_var=${2:-}          # optional VAR=value string
    echo "=== Running: $label ==="
    if [[ -n "$extra_var" ]]; then
        env "$extra_var" \
            "$RUSTLE" "$BAM" \
            --vg --vg-solver em --vg-snp \
            --genome-fasta "$FASTA" \
            -o "$OUTDIR/$label/out.gtf" \
            2>"$OUTDIR/$label/rustle.stderr"
    else
        "$RUSTLE" "$BAM" \
            --vg --vg-solver em --vg-snp \
            --genome-fasta "$FASTA" \
            -o "$OUTDIR/$label/out.gtf" \
            2>"$OUTDIR/$label/rustle.stderr"
    fi
    gffcompare -r "$REF" "$OUTDIR/$label/out.gtf" \
        -o "$OUTDIR/$label/cmp" \
        2>/dev/null
    local matching
    matching=$(grep "Matching transcripts:" "$OUTDIR/$label/cmp.stats" \
               | awk '{print $NF}')
    echo "  Matching transcripts: $matching"
    echo "$matching" > "$OUTDIR/$label/matching_transcripts.txt"
}

run_condition "off"      "RUSTLE_VG_NO_BORROW=1"
run_condition "legacy"   "RUSTLE_VG_BORROW_LEGACY=1"
run_condition "enhanced"

OFF=$(cat "$OUTDIR/off/matching_transcripts.txt")
LEG=$(cat "$OUTDIR/legacy/matching_transcripts.txt")
ENH=$(cat "$OUTDIR/enhanced/matching_transcripts.txt")

echo ""
echo "=== GOLGA8 borrow benchmark summary ==="
echo "Baseline (StringTie 3.0): 6 exact matches (historical reference)"
printf "%-12s %s\n" "Condition" "Matching transcripts"
printf "%-12s %s\n" "---------" "--------------------"
printf "%-12s %s\n" "OFF"       "$OFF"
printf "%-12s %s\n" "Legacy"    "$LEG"
printf "%-12s %s\n" "Enhanced"  "$ENH"
SCRIPT
chmod +x bench/vg_borrow_benchmark.sh
```

- [ ] **Step 2: Verify script syntax**

```bash
bash -n bench/vg_borrow_benchmark.sh && echo "syntax OK"
```
Expected: `syntax OK`

---

## Task 5: Run the benchmark and document results

- [ ] **Step 1: Run the benchmark**

```bash
bash bench/vg_borrow_benchmark.sh 2>&1 | tee bench/vg_borrow_benchmark/run.log
```

Expected output (exact numbers TBD — success criterion is Enhanced ≥ Legacy ≥ OFF):
```
=== GOLGA8 borrow benchmark summary ===
Baseline (StringTie 3.0): 6 exact matches (historical reference)
Condition    Matching transcripts
---------    --------------------
OFF          ??
Legacy       ??
Enhanced     ??
```

- [ ] **Step 2: Check for regression vs existing baseline**

The existing `golga8_fp_cmp.stats` has 12 matching transcripts. The Enhanced condition should be ≥ 12:

```bash
ENH=$(cat bench/vg_borrow_benchmark/enhanced/matching_transcripts.txt)
BASELINE=12
if [ "$ENH" -ge "$BASELINE" ]; then
    echo "PASS: Enhanced ($ENH) >= baseline ($BASELINE)"
else
    echo "REGRESSION: Enhanced ($ENH) < baseline ($BASELINE)"
fi
```

- [ ] **Step 3: Document results in RESULTS.md**

Append a new section to `bench/multi_copy_eval/RESULTS.md`:

```markdown
## Borrow-ON/OFF controlled benchmark (2026-05-23)

Three conditions, single binary, GOLGA8 region (`golga8_region.bam`),
reference `ref_GOLGA8.gff`. Metric: matching intron chains (gffcompare).

| Condition | Env var | Matching transcripts |
|-----------|---------|---------------------|
| OFF | `RUSTLE_VG_NO_BORROW=1` | ?? |
| Legacy | `RUSTLE_VG_BORROW_LEGACY=1` | ?? |
| Enhanced | _(default)_ | ?? |
| StringTie 3.0 (historical) | — | 6 |

**Interpretation:** [fill in after run]
```
Replace `??` with actual numbers from Step 1.

- [ ] **Step 4: Commit** (fill in actual OFF/Legacy/Enhanced numbers from Step 1 before committing)

```bash
git add bench/vg_borrow_benchmark.sh bench/multi_copy_eval/RESULTS.md \
        bench/vg_borrow_benchmark/
# Replace OFF, LEGACY, ENH with actual numbers from the run log
git commit -m "bench: VG borrow-ON/OFF three-condition GOLGA8 benchmark

Results: OFF=<n> Legacy=<n> Enhanced=<n>

Co-Authored-By: Claude Sonnet 4.6 <noreply@anthropic.com>"
```

---

## Interpretation Guide

After the benchmark runs, use this decision table:

| Result | Interpretation | Next step |
|--------|---------------|-----------|
| Enhanced > Legacy > OFF | Borrowing works; enhanced formula is better | Investigate which copies recovered; consider raising junction multiplier further |
| Enhanced ≈ Legacy > OFF | Coverage formula doesn't matter much; borrowing itself helps | Check if junction injection (0.5×) is the main driver |
| Enhanced ≈ Legacy ≈ OFF | Borrowing has no effect on GOLGA8 exact matches | Investigate a different bottleneck (EM assignment, junction support thresholds, guide mode) |
| Enhanced < Legacy | Enhanced formula introduces noise (over-boosting) | Tune the per-copy-expected formula; consider a cap like `min(expected_per_copy, max_sibling_cov)` |

## Environment Variable Summary

| Variable | Effect |
|----------|--------|
| `RUSTLE_VG_NO_BORROW=1` | Skip coverage + junction borrow maps entirely |
| `RUSTLE_VG_BORROW_LEGACY=1` | Use old 25%/cap-4.0 structural prior |
| `RUSTLE_VG_BORROW_FLOOR=1` | Use original blind 5% floor (pre-existing) |
| _(none of the above)_ | Enhanced per-copy-expected floor (new default) |
