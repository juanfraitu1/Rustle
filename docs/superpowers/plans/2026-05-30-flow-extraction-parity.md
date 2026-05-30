# Flow-Extraction Parity Port (sub-project 2) — Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Under `RUSTLE_FLOW_ST=1` (default OFF, composes with `RUSTLE_PREDCLUSTER_ST`), make Rustle's long-read flow extraction reproduce StringTie's candidate set and per-path coverage, so that with the sub-project-1 selection the final output approaches ST's (F1 → 95.27 at the 2a milestone, toward 96.3 at 2b). Flip both flags to default if combined F1 > baseline 93.78.

**Architecture:** Activate the existing `parse_trflong_st.rs` `_st` variants (`long_max_flow_st`, `parse_trflong_st`, `back_to_source_fast_long_st`, `fwd_to_sink_fast_long_st` — currently wired in COMPARISON/diff mode only) as a real production dispatch behind `RUSTLE_FLOW_ST`. Phase 2a = coverage-allocation parity (`long_max_flow_st` depletion); phase 2b = candidate-enumeration parity (seed order + path tracing). Parity-diffs are the gates.

**Tech Stack:** Rust (`cargo build --release`), instrumented ST `./tools/stringtie/stringtie`, Python 3 (`bench/gtf_chain_diff.py`, `bench/selection_oracle.py`, a new `bench/path_extracted_diff.py`), the parity JSONL (`path_extracted` carries source/cov/longcov/introns).

**Reference:** spec `docs/superpowers/specs/2026-05-30-flow-extraction-parity-design.md`; sub-project 1 spec/plan (`2026-05-30-predcluster-selection-parity*`). ST: `long_max_flow` rlink.cpp:8471 (depletion tail 8627-8665), `parse_trflong` 9758, `back_to_source_fast_long`/`fwd_to_sink_fast_long` (ST equivalents). Rustle: `src/rustle/parse_trflong_st.rs` (long_max_flow_st:950, parse_trflong_st:1036, back/fwd_st), `max_flow.rs:1085` (long_max_flow_direct), `path_extract.rs` (parse_trflong:5863, main loop ~6507, back_to_source:4913, fwd_to_sink:3943, the _st comparison dispatch ~7077-7400).

**Key facts:**
- Baseline (deterministic): TP=1750/FN=64/FP=168/F1=93.78 vs `../GGO_19.gtf` (1814 chains). ST: 1754/60/67 (96.3 Pr). Sub-project-1 flag-ON (`RUSTLE_PREDCLUSTER_ST=1`): 1734/80/155, F1 93.65 (capped by coverage).
- chain_counts: multi-intron (strand, intron-tuple), intron=(exon_end,next_start).
- `path_extracted` cov: ST is tlen-scaled (e.g. 55974) vs Rustle per-bp — compare WITHIN-LOCUS RATIOS (cov / locus-dominant-cov), not absolute.
- The `_st` variants exist but run in comparison mode (`emit_diff_if_diverges`); no `RUSTLE_FLOW_ST` flag yet. `st_flow()` predicate to be TDD'd in Task 2.
- Flags default-OFF → default path byte-unaffected (verify each step). Transient regression of flag-ON accepted (parity per phase, F1 at phase ends).

---

### Task 1: Oracle — confirm the decomposition reaches ST (abort gate, NO production code)

**Files:** Modify `bench/selection_oracle.py` (add a sub-project-2 mode) or create `bench/flow_oracle.py`.

**Context:** Verify that extraction-parity (ST candidates+cov) + sub-project-1 selection = ST's output, with no leftover divergence outside scope (graph boundaries, emit). Mechanism: feed ST's candidate set into Rustle's `select_predictions_st` is hard cross-process; instead measure two bounds.

- [ ] **Step 1: ST's final F1 (the absolute ceiling)** — `gtf_chain_diff.py` ST vs ref: confirms ceiling = ST 1754/60/67 (96.3 Pr). Already known; reconfirm.
- [ ] **Step 2: Selection-on-ST-candidates bound** — extend `selection_oracle.py`: instead of "keep ST winners on matching clusters," compute the set = ST's `path_extracted` candidates that PASS sub-project-1's selection criteria as ALREADY measured (use ST's final GTF as the proxy for "ST candidates that survive ST selection"). Compare to ST final. Then compute: of ST's final chains, how many are reachable if Rustle extracted ST's candidate set (i.e. ST_final ∩ everything) → this is trivially ST_final. The MEANINGFUL check: are there ST final chains that even ST's own candidates+selection wouldn't give Rustle because of a graph/emit difference? Quantify ST_final chains whose introns are NOT all present in Rustle's accepted junctions (un-buildable in Rustle's graph) — these are outside sub-project-2 scope.
- [ ] **Step 3: Abort gate** — PROCEED if ≥~95% of ST's final chains are buildable in Rustle's graph (junctions present) — confirming extraction+selection is the complete decomposition. ABORT/re-scope if a large fraction of ST final chains are un-buildable (graph-boundary parity needed first). Record counts.
- [ ] **Step 4: Commit** (`feat(flow): sub-project-2 oracle — decomposition reachability check` with the numbers + verdict).

---

### Task 2: `RUSTLE_FLOW_ST` flag + assess `_st` scaffold + production dispatch (default unchanged)

**Files:** Modify `src/rustle/stringtie_parity.rs`, `src/rustle/path_extract.rs`, `src/rustle/parse_trflong_st.rs`.

- [ ] **Step 1: Failing test for the flag** — in stringtie_parity.rs tests:
```rust
    #[test]
    fn st_flow_default_off() {
        use super::st_flow_from;
        assert!(!st_flow_from(None));
        assert!(st_flow_from(Some("1")));
        assert!(!st_flow_from(Some("0")));
    }
```
- [ ] **Step 2: Run → fail.** `cargo test -p rustle st_flow_default_off 2>&1 | tail -10`.
- [ ] **Step 3: Implement the predicate** in stringtie_parity.rs:
```rust
#[inline]
pub fn st_flow_from(v: Option<&str>) -> bool { matches!(v, Some(s) if s != "0") }
/// ST-faithful flow extraction (default OFF). Enable RUSTLE_FLOW_ST=1. Composes with
/// st_predcluster(). See docs/superpowers/specs/2026-05-30-flow-extraction-parity-design.md.
#[inline]
pub fn st_flow() -> bool { st_flow_from(std::env::var("RUSTLE_FLOW_ST").ok().as_deref()) }
```
- [ ] **Step 4: Run → pass.**
- [ ] **Step 5: Assess + document the `_st` scaffold** — read `parse_trflong_st.rs` (long_max_flow_st:950, parse_trflong_st:1036, back/fwd_st) and the comparison dispatch (path_extract.rs:7077-7400). Write a 1-paragraph note (in the commit message or a code comment) on what's implemented vs stubbed: does `long_max_flow_st` do real depletion, or just diff against the Rustle path? Does `parse_trflong_st` produce a full extraction or only compare back/fwd tracing? This determines how much Phase 2a/2b must add.
- [ ] **Step 6: Wire production dispatch** — where the comparison-mode `_st` calls live (path_extract.rs:7252+ `comparison_active()`), add a `st_flow()` branch that uses the `_st` extraction RESULT in production (not just diffs it). With the `_st` variants in whatever state they're in, flag-ON may differ from baseline — that's expected; the point is the wiring + default-OFF safety.
- [ ] **Step 7: Build + default-unchanged guard** (cargo build; flag-OFF Intron chain ~96.5/91.2). Report flag-ON behavior (may diverge; that's the starting point for 2a/2b).
- [ ] **Step 8: Commit** (`feat(flow): RUSTLE_FLOW_ST flag + activate _st extraction dispatch (default off)`).

---

### Task 3 (Phase 2a): coverage-allocation parity in `long_max_flow_st`

**Files:** Modify `src/rustle/parse_trflong_st.rs` (long_max_flow_st) / `src/rustle/max_flow.rs`; Create `bench/path_extracted_diff.py`.

**Context:** Port ST's `long_max_flow` abundance-subtraction depletion (rlink.cpp:8471, tail 8627-8665): after extracting a path, subtract allocated flow from each participating transfrag's abundance (clamp 0 below DBL_ERROR), so competing-path residual coverage matches ST. Target: the within-locus cov ratio of the 42 contested minorities (ST ~40/dom vs Rustle ~24.6/dom). NOT a cov formula change (those were no-ops) — the depletion algorithm.

- [ ] **Step 1: Build `bench/path_extracted_diff.py`** — join Rustle vs ST `path_extracted` by (strand, intron-chain); on matched chains report the WITHIN-LOCUS cov ratio (cov / locus-max-cov) distribution and the count of contested minorities where Rustle's ratio << ST's (the 42). This is the 2a gate.
- [ ] **Step 2: Complete `long_max_flow_st` depletion** to faithfully mirror rlink.cpp:8627-8665 (per-node, per-transfrag flow subtraction, DBL_ERROR clamp), using Rustle's flow structures (the Task-2 assessment says how much is already there). Gate under `st_flow()`.
- [ ] **Step 3: Build + default guard.**
- [ ] **Step 4: Gate (cov parity)** — `RUSTLE_FLOW_ST=1` path_extracted vs ST via `path_extracted_diff.py`: the contested-minority cov ratios move toward ST's. Record.
- [ ] **Step 5: Combined F1 (the 2a milestone)** — `RUSTLE_FLOW_ST=1 RUSTLE_PREDCLUSTER_ST=1` chain-level TP/FN/FP/F1 (3 runs). EXPECTED: the ~42 previously over-killed TPs are recovered (selection no longer over-kills them now their cov clears ST-thresholds) → F1 toward 95.27. Report vs baseline 93.78 and sub-project-1-alone 93.65.
- [ ] **Step 6: Commit** (`feat(flow): ST-faithful long_max_flow_st depletion — coverage parity, +<>pp combined F1`).

---

### Task 4 (Phase 2b-i): re-evaluate ASC seed order composed with selection

**Files:** none (measurement) or a small dispatch tweak.

**Context:** ASC order (`RUSTLE_PARSE_TRFLONG_ST_ORDER=1`) failed standalone (−33 TP) because Rustle's selection picked wrong winners under minority-first. Now that selection is ST-faithful (`RUSTLE_PREDCLUSTER_ST`) AND coverage is ST-faithful (Task 3), re-test ASC composed.

- [ ] **Step 1: Measure** `RUSTLE_FLOW_ST=1 RUSTLE_PREDCLUSTER_ST=1 RUSTLE_PARSE_TRFLONG_ST_ORDER=1` chain F1 vs without the order flag (3 runs each). If ASC now HELPS (recovers minority isoforms without the prior TP loss), fold it into `st_flow()` (make `st_flow()` imply ASC order). If it still hurts, leave it off and record why (the composition still doesn't fix it).
- [ ] **Step 2: Commit** the decision (`feat/docs(flow): ASC seed order composed with ST selection+coverage — <result>`).

---

### Task 5 (Phase 2b-ii): path-tracing parity (back_to_source / fwd_to_sink)

**Files:** Modify `src/rustle/parse_trflong_st.rs` (back_to_source_fast_long_st / fwd_to_sink_fast_long_st).

**Context:** Close the 70 path-enum FP (ST never extracts) + 26 FN by matching ST's source/sink path tracing so seeds aren't over-extended into the dominant TSS (the alt-TSS collapse). ST: back_to_source/fwd_to_sink linking in parse_trflong.

- [ ] **Step 1: Candidate-set diff** — `path_extracted_diff.py` extended to chain-SET mode: under `RUSTLE_FLOW_ST=1`, count Rustle-only candidates (the 70) + ST-only (the 26) vs ST path_extracted. This is the 2b gate.
- [ ] **Step 2: Port the tracing parity** in back_to_source_fast_long_st / fwd_to_sink_fast_long_st to match ST's linking (the diff data from the comparison mode + the alt-TSS finding guide the exact divergence). Gate under `st_flow()`.
- [ ] **Step 3: Build + default guard.**
- [ ] **Step 4: Gate** — Rustle-only / ST-only candidate counts drop toward 0. Record.
- [ ] **Step 5: Combined F1** — `RUSTLE_FLOW_ST=1 RUSTLE_PREDCLUSTER_ST=1` chain F1 (3 runs) → toward 96.3. Report.
- [ ] **Step 6: Commit** (`feat(flow): ST-faithful path tracing — candidate-set parity, combined F1 <>`).

---

### Task 6: Combined validation + default decision

**Files:** Modify `docs/STRINGTIE_PARITY_FINDINGS.md`.

- [ ] **Step 1: Final combined measurement** — `RUSTLE_FLOW_ST=1 RUSTLE_PREDCLUSTER_ST=1` chain TP/FN/FP/Sn/Pr/F1 (3 runs, determinism) + `gtf_chain_diff.py` vs ST (final-chain parity). Confirm default-OFF unchanged.
- [ ] **Step 2: Default decision** — if combined F1 > baseline 93.78 with Sn not hurt, flip BOTH `st_flow()` and `st_predcluster()` default-on (opt-out). Else keep opt-in; record the achieved combined F1 and the residual.
- [ ] **Step 3: Record** — append "Flow-extraction parity (sub-project 2) — DONE" to `docs/STRINGTIE_PARITY_FINDINGS.md`: 2a coverage-parity result + 2a-milestone F1, 2b candidate-parity result + combined F1, the default decision, residual gap vs ST's 96.3. Commit.

---

## Self-Review
- **Spec coverage:** §3 architecture (RUSTLE_FLOW_ST, activate _st variants, compose with sub-project 1) → Task 2. §4 2a coverage (`long_max_flow_st` depletion) → Task 3; 2b seed order → Task 4, path tracing → Task 5. §6 validation (cov-ratio + candidate-set parity gates, F1 at phase ends) → Tasks 3/5/6. §7 oracle-first + abort → Task 1; default flip → Task 6. §2 non-goals (graph-boundary minority, guided) honored — Task 1 abort-checks graph buildability, doesn't fix it.
- **Placeholder note:** the porting tasks (3,5) are specified by ST source-of-truth + the parity gate + the discipline because the `_st` variants are a research port read-from-ST; the oracle/flag/scaffold/gates/commands are concrete. Task 2's scaffold-assessment is explicitly an investigation step that informs how much 3/5 add.
- **Type consistency:** `st_flow()`/`st_flow_from()` + `RUSTLE_FLOW_ST` consistent; composes with the existing `st_predcluster()`/`RUSTLE_PREDCLUSTER_ST` and `RUSTLE_PARSE_TRFLONG_ST_ORDER`; gates are the established parity-diff pattern (path_extracted cov-ratio + chain-set, gtf_chain_diff).
- **Risk:** Task 1 abort can end early (graph-boundary parity needed first). Task 3 (depletion) is the hard phase; Task 4 (ASC) may stay off. The +1.49pp lands at the Task-3 (2a) milestone even if 2b is partial.
