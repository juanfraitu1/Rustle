# VG over-collapse recovery — clean re-bundle + primary-gated union — Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Make VG mode recover the ~109 baseline-recoverable annotated isoforms it currently drops at secondary-bearing family bundles (over-collapse), via a clean re-bundle of primary reads feeding the existing strictly-additive union plus a primary-support gate — opt-in, default byte-identical.

**Architecture:** Replace the broken clone-and-patch sub-bundle production (pipeline.rs:12185-12225) with a *build-fresh* `Bundle` whose strand is derived from its own reads (not inherited from the polluted parent mega-bundle), routed through the existing assembly par_iter and the existing holdout→union stage (pipeline.rs:19337). Add a `longcov ≥ N` gate to the union.

**Tech Stack:** Rust (the `rustle` crate). Tests via `cargo test --lib`. Spec: `docs/superpowers/specs/2026-06-09-vg-union-baseline-rebundle-design.md`.

---

## Context the engineer needs

- **Working dir:** `/mnt/c/Users/jfris/Desktop/Rustle`. Branch `vg/flow-capacity-apportionment` (do not switch).
- **The defect:** `RUSTLE_VG_UNION_BASELINE` (pipeline.rs:12150-12235) re-splits a secondary-bearing bundle's primary reads, but builds each sub-bundle via `sub = b.clone()` (line 12198) — inheriting the parent mega-bundle's `strand` (and `hp_tag`/`ps_tag`). A re-split sub-gene on the OPPOSITE strand keeps the parent's strand → mis-assembly. Recovers only ~6/13.
- **`Bundle` struct (src/rustle/types.rs:83), 15 fields:** `chrom, start, end, strand, reads, junction_stats, junction_pair_stats, bundlenodes, read_bnodes, bnode_colors, synthetic, rescue_class, vg_family_id, hp_tag, ps_tag`.
- **`BundleRead`** has `strand: char` and `weight: f64`.
- **Helper:** `compute_initial_junction_stats_for_reads(reads: &[BundleRead], bundle_start: u64, bundle_end: u64, config: &RunConfig) -> (JunctionStats, HashMap<(Junction,u64),u32>)` (pipeline.rs:1370).
- **Union stage (pipeline.rs:19343-19365):** drains `union_baseline_holdout`, unions a tx iff multi-exon AND its intron chain is absent from `all_transcripts`. `Transcript` has `longcov: f64` (path_extract.rs:641).
- **Run cargo tests with the OOM watchdog** (per `project_genomewide_oom_protocol`): the runaway test `dp_scales_to_many_reads` is already `#[ignore]`; suite is ~277 tests. Use `cargo test --lib`.
- **OOM-safe validation:** per-locus slices only; whole-chromosome `--vg` is ~2.2GB serial (fine), never concurrent.

## File structure

- Modify `src/rustle/pipeline.rs`: add `majority_strand` + `build_fresh_baseline_subbundle` helpers; replace the clone block; add the `longcov` gate.
- Tests: Rust `#[cfg(test)]` unit tests in `pipeline.rs` (or its test module) for the pure helpers; integration validation via the harness in `bench/flow_recall_phase0/`.

---

## Task 1: `majority_strand` pure helper (the residual-bug fix)

**Files:** Modify `src/rustle/pipeline.rs` (add fn + test).

- [ ] **Step 1: Write the failing test**

Add to the test module in `src/rustle/pipeline.rs` (find `#[cfg(test)] mod tests` or add one):
```rust
#[test]
fn majority_strand_picks_read_weighted_majority() {
    // minus-weight dominates -> '-', regardless of any parent bundle strand
    assert_eq!(super::majority_strand([('+', 1.0), ('-', 3.0), ('-', 2.0)].into_iter()), '-');
    // plus dominates -> '+'
    assert_eq!(super::majority_strand([('+', 5.0), ('-', 1.0)].into_iter()), '+');
    // tie / empty -> '+' (deterministic default)
    assert_eq!(super::majority_strand([('+', 1.0), ('-', 1.0)].into_iter()), '+');
    assert_eq!(super::majority_strand(std::iter::empty()), '+');
}
```

- [ ] **Step 2: Run test to verify it fails**

Run: `cargo test --lib majority_strand_picks_read_weighted_majority 2>&1 | tail -15`
Expected: compile error / FAIL — `majority_strand` not found.

- [ ] **Step 3: Write minimal implementation**

Add near the other free functions in `src/rustle/pipeline.rs` (e.g. just above `stamp_union_baseline_rescue_class` at line 713):
```rust
/// Read-weight majority strand for a fresh sub-bundle. The over-collapse re-bundle must derive
/// strand from the sub-bundle's OWN reads — NOT inherit the parent mega-bundle's strand, which is
/// the residual mis-assembly the clone path could not fix. Ties / empty default to '+'.
fn majority_strand(strand_weights: impl Iterator<Item = (char, f64)>) -> char {
    let (mut plus, mut minus) = (0.0f64, 0.0f64);
    for (s, w) in strand_weights {
        match s {
            '+' => plus += w,
            '-' => minus += w,
            _ => {}
        }
    }
    if minus > plus { '-' } else { '+' }
}
```

- [ ] **Step 4: Run test to verify it passes**

Run: `cargo test --lib majority_strand_picks_read_weighted_majority 2>&1 | tail -8`
Expected: PASS (test result: ok. 1 passed).

- [ ] **Step 5: Commit**

```bash
git add src/rustle/pipeline.rs
git commit -m "feat(vg): majority_strand helper for clean re-bundle (read-derived strand)"
```

---

## Task 2: `build_fresh_baseline_subbundle` + wire it into the block

**Files:** Modify `src/rustle/pipeline.rs` (add fn; replace clone block 12185-12225).

- [ ] **Step 1: Add the build-fresh helper**

Add just below `majority_strand`:
```rust
/// Build a FRESH primary-only sub-bundle from re-split reads — every field set from the reads,
/// nothing inherited from the (polluted) parent except `chrom`. This replaces the clone-and-patch
/// path whose inherited `strand`/`hp_tag`/`ps_tag` caused ~7/13 over-collapse sub-bundles to
/// mis-assemble. Tagged `UnionBaseline` so the end-stage union holds it out and unions it back by
/// novel intron chain (strictly additive).
fn build_fresh_baseline_subbundle(
    parent: &crate::types::Bundle,
    reads: Vec<crate::types::BundleRead>,
    config: &RunConfig,
) -> crate::types::Bundle {
    let start = reads
        .iter()
        .filter_map(|r| r.exons.first().map(|e| e.0))
        .min()
        .unwrap_or(parent.start);
    let end = reads
        .iter()
        .filter_map(|r| r.exons.last().map(|e| e.1))
        .max()
        .unwrap_or(parent.end);
    let strand = majority_strand(reads.iter().map(|r| (r.strand, r.weight)));
    let (junction_stats, junction_pair_stats) =
        compute_initial_junction_stats_for_reads(&reads, start, end, config);
    crate::types::Bundle {
        chrom: parent.chrom.clone(),
        start,
        end,
        strand,
        reads,
        junction_stats,
        junction_pair_stats,
        bundlenodes: None,
        read_bnodes: None,
        bnode_colors: None,
        synthetic: false,
        rescue_class: Some(crate::vg_family::diagnostic::RescueClass::UnionBaseline),
        vg_family_id: None,
        hp_tag: None,
        ps_tag: None,
    }
}
```

- [ ] **Step 2: Replace the clone-and-patch loop body**

In `src/rustle/pipeline.rs`, replace the block that currently reads (lines ~12185-12225):
```rust
            for grp in crate::bundle::split_spans_by_runoff(&spans, runoff) {
                let reads: Vec<crate::types::BundleRead> =
                    grp.iter().map(|&i| prim[i].clone()).collect();
                let start = reads
                    .iter()
                    .filter_map(|r| r.exons.first().map(|e| e.0))
                    .min()
                    .unwrap_or(b.start);
                let end = reads
                    .iter()
                    .filter_map(|r| r.exons.last().map(|e| e.1))
                    .max()
                    .unwrap_or(b.end);
                let mut sub = b.clone();
                sub.reads = reads;
                sub.start = start;
                sub.end = end;
                sub.bundlenodes = None;
                sub.read_bnodes = None;
                sub.bnode_colors = None;
                sub.synthetic = false;
                sub.vg_family_id = None;
                sub.rescue_class =
                    Some(crate::vg_family::diagnostic::RescueClass::UnionBaseline);
                let (js, jps) = compute_initial_junction_stats_for_reads(
                    &sub.reads,
                    sub.start,
                    sub.end,
                    &config,
                );
                sub.junction_stats = js;
                sub.junction_pair_stats = jps;
                clones.push((next_idx, sub));
                next_idx += 1;
            }
```
with:
```rust
            for grp in crate::bundle::split_spans_by_runoff(&spans, runoff) {
                let reads: Vec<crate::types::BundleRead> =
                    grp.iter().map(|&i| prim[i].clone()).collect();
                let sub = build_fresh_baseline_subbundle(b, reads, &config);
                clones.push((next_idx, sub));
                next_idx += 1;
            }
```
NOTE: read the actual lines first (`sed -n '12184,12226p' src/rustle/pipeline.rs`) and match the exact current text — line numbers may have shifted by Task 1's insertion. The semantic target is: the `for grp in split_spans_by_runoff` loop body that builds `sub` and pushes to `clones`.

- [ ] **Step 3: Verify it compiles**

Run: `cargo build --release 2>&1 | tail -5`
Expected: `Finished` (no errors). If `prim[i]` / `b` / `spans` / `runoff` are not in scope, you replaced the wrong span — re-read the block.

- [ ] **Step 4: Confirm default still byte-identical (no flag)**

Run: `cargo test --lib 2>&1 | tail -5`
Expected: suite passes (no behavioral change with the flag unset; this block is inside `if ... RUSTLE_VG_UNION_BASELINE`).

- [ ] **Step 5: Commit**

```bash
git add src/rustle/pipeline.rs
git commit -m "fix(vg): build-fresh re-bundle (derive strand) replaces clone-and-patch"
```

---

## Task 3: primary-support gate on the union

**Files:** Modify `src/rustle/pipeline.rs` (union stage ~19343-19365).

- [ ] **Step 1: Read the current union block**

Run: `sed -n '19343,19372p' src/rustle/pipeline.rs`
Confirm the `for mut t in std::mem::take(&mut union_baseline_holdout)` loop with the `if t.exons.len() < 2 { continue; }` and `if vg_chains.contains(&chain) || !seen.insert(chain) { continue; }` guards.

- [ ] **Step 2: Add the longcov gate**

Immediately before the `for mut t in std::mem::take(&mut union_baseline_holdout) {` line, add:
```rust
        let min_longcov: f64 = std::env::var("RUSTLE_VG_UNION_BASELINE_MIN_LONGCOV")
            .ok()
            .and_then(|s| s.parse().ok())
            .unwrap_or(2.0);
```
Then inside the loop, immediately after the existing `if t.exons.len() < 2 { continue; }` guard, add:
```rust
            if t.longcov < min_longcov {
                continue; // primary-support gate: only union real primary-backed dropped isoforms
            }
```

- [ ] **Step 3: Verify compile + suite**

Run: `cargo build --release 2>&1 | tail -3 && cargo test --lib 2>&1 | tail -4`
Expected: builds; suite passes.

- [ ] **Step 4: Commit**

```bash
git add src/rustle/pipeline.rs
git commit -m "feat(vg): primary-support gate (longcov>=N) on union-baseline recovery"
```

---

## Task 4: Integration — recovers the real regression, default byte-identical

**Files:** none (validation commands).

- [ ] **Step 1: Confirm default (no flag) is byte-identical on a regression locus**

```bash
PC=/mnt/c/Users/jfris/Desktop/Rustle/bench/copy_recovery_eval/results_genomewide/perchrom
F=/home/juanfra/winloci_scratch/GGO.fasta
R=target/release/rustle
samtools view -b $PC/NC_073224.2/c.bam NC_073224.2:39996000-40059000 > /tmp/un.bam && samtools index /tmp/un.bam
RUSTLE_VG_TANDEM=1 RUSTLE_VG_TANDEM_PRIMARY_JUNCTIONS=1 RUSTLE_VG_DECISIVE_GATE=1 RUSTLE_VG_DECISIVE_GATE_MIN_PRIM=4 $R --vg --vg-snp --genome-fasta $F -L /tmp/un.bam -o /tmp/un_off.gtf 2>/dev/null
```
Expected: runs clean. (This is the no-union VG output.)

- [ ] **Step 2: Confirm the flag recovers XM_063708549.1's chain**

```bash
RUSTLE_VG_UNION_BASELINE=1 RUSTLE_VG_TANDEM=1 RUSTLE_VG_TANDEM_PRIMARY_JUNCTIONS=1 RUSTLE_VG_DECISIVE_GATE=1 RUSTLE_VG_DECISIVE_GATE_MIN_PRIM=4 $R --vg --vg-snp --genome-fasta $F -L /tmp/un.bam -o /tmp/un_on.gtf 2>/tmp/un.log
grep VG-UNION-BASELINE /tmp/un.log
python3 -c "
import re,sys
sys.path.insert(0,'bench/flow_recall_phase0'); import lib
ref=lib.ref_chain('NC_073224.2','rna-XM_063708549.1')[2]
def has(p):
    tx={}
    for line in open(p):
        f=line.split(chr(9))
        if len(f)<9 or f[2]!='exon': continue
        m=re.search(r'transcript_id \"([^\"]+)\"',f[8]); tx.setdefault(m.group(1),[]).append((int(f[3]),int(f[4])))
    for ex in tx.values():
        ex.sort()
        if tuple((ex[i][1]+1,ex[i+1][0]-1) for i in range(len(ex)-1))==tuple(ref): return True
    return False
print('OFF has chain:', has('/tmp/un_off.gtf'), ' ON has chain:', has('/tmp/un_on.gtf'))
"
```
Expected: `OFF has chain: False  ON has chain: True` — the flag recovers the dropped isoform.
If ON is still False: the sub-bundle still mis-assembles — inspect with `RUSTLE_VG_DISCOVERY_TRACE`/the `[VG-UNION-BASELINE]` logs to see whether the sub-bundle was created and what strand it got (do NOT add more fixes blindly — return to root-cause analysis per systematic-debugging).

- [ ] **Step 3: Commit a short result note**

```bash
mkdir -p bench/flow_recall_phase0 && printf "OFF/ON recovery for NC_073224.2 XM_063708549.1: see un_off/un_on\n" >> bench/flow_recall_phase0/UNION_REBUILD_NOTES.md
git add bench/flow_recall_phase0/UNION_REBUILD_NOTES.md
git commit -m "test(vg): integration confirm union-baseline recovers over-collapse locus"
```

---

## Task 5: Genome-wide validation + record

**Files:** `bench/flow_recall_phase0/UNION_REBUILD_NOTES.md`.

- [ ] **Step 1: Run ±RUSTLE_VG_UNION_BASELINE across the high-regression chromosomes**

For each of NC_073247.2, NC_073229.2, NC_073242.2, NC_073224.2, NC_073228.2 (the top-regression contigs): run `--vg` full-config with and without `RUSTLE_VG_UNION_BASELINE=1` on `perchrom/<C>/c.bam`, gffcompare both vs `awk '$1==C' bench/copy_recovery_eval/results/ref.gtf`, and tally FSM + total tx. Serial, OOM-safe (~2.2GB each). Use a background script with a `free -m` watchdog killing any `rustle` proc >9GB (per `project_genomewide_oom_protocol`).

- [ ] **Step 2: Measure recall gain vs FP cost**

For each chrom: recovered regressions (from `bench/flow_recall_phase0/vg_regressions.json`) now in the ON FSM set, and the added non-annotated tx (FP). Compute the genome-wide recall:FP ratio. **Gate:** ratio should be clearly > 1.0 (the validated direction was 1.27); if the `longcov` default of 2 yields a worse ratio, sweep `RUSTLE_VG_UNION_BASELINE_MIN_LONGCOV` ∈ {2,3,5} and pick the knee.

- [ ] **Step 3: Record results + recommendation**

Write the per-chrom table, the chosen `MIN_LONGCOV`, recovered-regression count, FP cost, and ratio to `bench/flow_recall_phase0/UNION_REBUILD_NOTES.md`. State whether it meets the net-positive bar and is a default-flip candidate (a SEPARATE decision, not done here). Update memory `project_vg_baseline_regression.md` with the outcome.

```bash
git add bench/flow_recall_phase0/UNION_REBUILD_NOTES.md
git commit -m "docs(vg): genome-wide validation of union-baseline rebundle recovery"
```

---

## Phase 2 (NOT in this plan)

Flipping `RUSTLE_VG_UNION_BASELINE` to default-on is a separate post-validation decision (mirrors the strand-bundling flip), contingent on Task 5 showing a net-positive genome-wide ratio.

## Self-review notes

- Spec coverage: build-fresh re-bundle (T1+T2), primary-gated union (T3), env gating (default-off, inside the existing `RUSTLE_VG_UNION_BASELINE` guard — T2/T3), integration + byte-identical guard (T4), genome-wide validation (T5), default-flip out of scope. All spec sections mapped.
- The unit-testable unit is `majority_strand` (the residual-bug fix); the build-fresh wiring and gate are validated by integration (T4) + suite byte-identical (T2/T3), since they are pipeline wiring, not pure functions.
- Type consistency: `build_fresh_baseline_subbundle(parent, reads, config)` and `majority_strand(iter of (char,f64))` are used consistently; `RUSTLE_VG_UNION_BASELINE_MIN_LONGCOV` default 2.0 matches the spec.
