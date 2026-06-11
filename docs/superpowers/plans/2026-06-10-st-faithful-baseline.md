# ST-Faithful Baseline — Implementation Plan (Milestone 1: converge mini3)

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Flip rustle's default to ST-faithful and converge the mini3 fixture to `rustle-vs-StringTie = 0/0` (currently 3 rustle-only / 16 ST-only), with today's behavior preserved under `RUSTLE_PRECISE=1`.

**Architecture:** A master `precise_mode()` gate. Default → ST-faithful (the converging baseline). `RUSTLE_PRECISE=1` → today's strict-early behavior, held byte-identical to commit `4705ab1`. Divergences are ported in coherent stage-groups (G4 extraction → G1 junctions → G3 filters), each tracer-driven on the 0.1s mini3 loop. Metric is rustle-vs-StringTie, not vs annotation.

**Tech Stack:** Rust (`src/rustle/`), the parity harness (`/tmp` JSONL + `bench/gtf_chain_diff.py`), `bench/mini3/` fixture.

**Design:** `docs/superpowers/specs/2026-06-10-st-faithful-baseline-design.md`. **Revert point:** git `4705ab1`.

**Scope note:** This is Milestone 1 (mini3 → 0/0). Full-BAM and genome-wide convergence are follow-on plans. Each port task is *discovery-driven*: the failing test is "specific mini3 chain(s) don't match ST"; the implementation is found by tracing that chain's first divergence, then porting that exact decision behind the gate. Code for the gate/tests is given in full; code for each ported decision is found during the task (the spec forbids guessing it).

---

## File Structure

- `src/rustle/stringtie_parity.rs` — add `precise_mode()` master gate (mirrors existing `st_shadow()`). Cache the env read.
- `src/rustle/graph_build.rs`, `junction_graph*.rs`, `killed_junctions.rs` — G1 junction-acceptance forks.
- `src/rustle/parse_trflong_st.rs`, `path_extract.rs` — G4 extraction forks.
- `src/rustle/transcript_filter.rs` — G3 filtering forks.
- `bench/mini3/` — fixture + `check.sh` (exists) + new `expected_precise.gtf` (frozen `4705ab1` behavior) + `check_precise.sh` (escape-hatch invariant).
- `tests/` (or `#[cfg(test)]`) — the convergence + invariant integration assertions.

Each precise-vs-ST-faithful fork uses the exact pattern:
```rust
if crate::stringtie_parity::precise_mode() {
    // today's strict-early behavior (unchanged)
} else {
    // ST-faithful behavior (the new default)
}
```

---

### Task 1: `precise_mode()` gate + escape-hatch reference

**Files:**
- Modify: `src/rustle/stringtie_parity.rs`
- Create: `bench/mini3/expected_precise.gtf`, `bench/mini3/check_precise.sh`

- [ ] **Step 1: Freeze today's behavior as the precise reference**

The current default (HEAD = `4705ab1`) is today's strict-early behavior. Capture it as the `RUSTLE_PRECISE=1` target BEFORE any flip:
```bash
target/release/rustle -L bench/mini3/mini3.bam -o bench/mini3/expected_precise.gtf 2>/dev/null
```

- [ ] **Step 2: Write the failing test (escape-hatch invariant)**

Create `bench/mini3/check_precise.sh`:
```bash
#!/usr/bin/env bash
# Invariant: RUSTLE_PRECISE=1 must byte-match the frozen 4705ab1 behavior.
set -euo pipefail
cd "$(git rev-parse --show-toplevel)"
OUT=$(mktemp --suffix=.gtf)
RUSTLE_PRECISE=1 target/release/rustle -L bench/mini3/mini3.bam -o "$OUT" 2>/dev/null
if diff <(grep -v '^#' bench/mini3/expected_precise.gtf) <(grep -v '^#' "$OUT") >/dev/null; then
  echo "ESCAPE-HATCH OK: RUSTLE_PRECISE=1 byte-matches 4705ab1"
else
  echo "ESCAPE-HATCH DRIFT — RUSTLE_PRECISE=1 diverged from 4705ab1"; diff <(grep -v '^#' bench/mini3/expected_precise.gtf) <(grep -v '^#' "$OUT") | head; exit 1
fi
rm -f "$OUT"
```
Make it executable. Run it now: it passes trivially (gate not yet added, default == precise). This is the regression guard for every later task.

- [ ] **Step 3: Add the `precise_mode()` gate**

In `src/rustle/stringtie_parity.rs`, mirroring the existing `st_shadow()` (~line 193):
```rust
/// Master gate for the ST-faithful baseline flip. DEFAULT (unset) = ST-faithful
/// (the converging baseline). `RUSTLE_PRECISE=1` = today's strict-early behavior
/// (the escape hatch; must stay byte-identical to commit 4705ab1).
pub fn precise_mode() -> bool {
    use std::sync::OnceLock;
    static P: OnceLock<bool> = OnceLock::new();
    *P.get_or_init(|| std::env::var_os("RUSTLE_PRECISE").is_some())
}
```

- [ ] **Step 4: Build + verify both invariants hold (no forks yet)**

```bash
cargo build --release 2>&1 | tail -3
bash bench/mini3/check_precise.sh      # ESCAPE-HATCH OK
bash bench/mini3/check.sh               # still 3 / 16 (no behavior change yet)
```
Expected: gate compiles, both checks unchanged (the gate is dormant until forks reference it).

- [ ] **Step 5: Commit**

```bash
git add src/rustle/stringtie_parity.rs bench/mini3/expected_precise.gtf bench/mini3/check_precise.sh
git commit -m "st-faithful: add precise_mode() gate + escape-hatch invariant (mini3)"
```

---

### Task 2: G4 extraction — converge the 10 "transfrag-exists-not-extracted" ST-only chains

**Files:** `src/rustle/parse_trflong_st.rs`, `src/rustle/path_extract.rs` (the `_st` extension / canonical path).

These 10 chains have a rustle transfrag but rustle's extraction is stricter than StringTie's and doesn't emit them. The known shape (§6q–s): rustle's back/fwd extension rejects low-cov ~1.0 sub-paths ST admits.

- [ ] **Step 1: Identify the failing chains (the test)**

```bash
STRINGTIE_PARITY_LOG=/tmp/m3_st.jsonl tools/stringtie/stringtie -L bench/mini3/mini3.bam -o /tmp/m3_st.gtf 2>/dev/null
RUSTLE_PARITY_LOG=/tmp/m3_ru.jsonl target/release/rustle -L bench/mini3/mini3.bam -o /tmp/m3_ru.gtf 2>/dev/null
python3 bench/gtf_chain_diff.py /tmp/m3_ru.gtf /tmp/m3_st.gtf   # 3 / 16 baseline
```
Pick one G4 ST-only chain (transfrag exists, not extracted) — locus B (36012544-36069321) has most. Record its intron chain.

- [ ] **Step 2: Trace its first divergence**

Run the first-divergence tracer (the §6w worked example: junction_accept → transfrag_collapse → backfwd_extension/extension_scan_step → node_flux → path_extracted) on that chain from `/tmp/m3_{ru,st}.jsonl`, scoped to locus B. Identify the exact extension step where ST commits a candidate rustle skips (the `extension_scan_step` `skip`/`lower_cov_reject` vs ST `select`). This names the precise decision (e.g. the canonical free-sink / low-cov sub-path admission, `parse_trflong_st.rs:375` region).

- [ ] **Step 3: Port that decision behind the gate**

Wrap the divergent extension predicate: `if precise_mode() { strict } else { ST-faithful (admit the candidate ST admits) }`. Use the trace to pin the exact condition; do not guess. (This is likely turning the canonical `_st` path / the relaxed admission ON by default — i.e. making `!precise_mode()` take the ST-faithful branch the canonical mode already implements.)

- [ ] **Step 4: Build + verify mini3 improves AND escape-hatch holds**

```bash
cargo build --release 2>&1 | tail -3
bash bench/mini3/check_precise.sh      # ESCAPE-HATCH OK (RUSTLE_PRECISE=1 unchanged)
bash bench/mini3/check.sh               # ST-only should DROP (toward 0); rustle-only may transiently rise
```
Expected: ST-only decreases. If rustle-only rises (the permissive extraction over-produces), that is accepted transient regression — it will be reclaimed by the G3 filter task. Record the new counts.

- [ ] **Step 5: Commit**

```bash
git add src/rustle/parse_trflong_st.rs src/rustle/path_extract.rs
git commit -m "st-faithful(G4): admit ST's low-cov extension candidates by default (gate strict behind RUSTLE_PRECISE)"
```

- [ ] **Step 6: Repeat Steps 1–5 for the remaining G4 chains until the 10 are extracted**

Each iteration: pick a still-missing G4 chain, trace, port the next divergent predicate behind the gate, verify mini3 ST-only drops + escape-hatch holds, commit. Stop when `transfrag-exists-not-extracted` ST-only = 0 (re-run the Task-categorizer from the design session to recount).

---

### Task 3: G1 — converge the 6 never-built chains (junctions / read-through)

**Files:** `src/rustle/graph_build.rs`, `junction_graph*.rs`, `killed_junctions.rs`.

The 4 ST-only never-built (locus C read-through across LOC109023803−→ANXA2R+) + 2 rustle-only ST-never-built. Baseline rustle must *reproduce* ST's read-through artifact (the strict-early junction/extension avoidance moves behind the gate).

- [ ] **Step 1: Failing test** — confirm the 4 locus-C ST chains are absent from rustle (`check.sh`), and the 2 rustle-only chains ST lacks are present.

- [ ] **Step 2: Trace** — for a locus-C ST chain, trace junction_accept (does rustle accept the junctions ST uses to bridge the two genes?) and the extension. For a rustle-only chain, trace why ST never builds it (a junction ST rejects that rustle accepts, or an extension ST stops). Identify the exact junction-acceptance / extension-termination decision.

- [ ] **Step 3: Port behind the gate** — make rustle accept ST's junction set (the read-through-enabling junctions) and stop where ST stops, by default; gate the strict-early rejection behind `precise_mode()`. ⚠ Expect this to be where the default regresses most (permissive junctions); it is reclaimed by G3.

- [ ] **Step 4: Build + verify** — `check_precise.sh` (OK), `check.sh` (locus-C ST chains now built; the 2 rustle-only resolved). Accept transient rustle-only rise.

- [ ] **Step 5: Commit** — `git commit -m "st-faithful(G1): match ST junction acceptance + read-through by default"`

- [ ] **Step 6: Repeat** until the 6 never-built/extra are resolved.

---

### Task 4: G3 — converge the filter divergences + reclaim transient rustle-only

**Files:** `src/rustle/transcript_filter.rs` (print_predcluster analog).

The 2 ST-only built-then-killed + 1 rustle-only ST-killed + any rustle-only over-production accumulated from G1/G4. Match ST's `pred_filter` cascade (`included_drop`, `retained_intron`, `isofrac`, readthr) exactly, on the now-ST-faithful coverage.

- [ ] **Step 1: Failing test** — `check.sh` rustle-only count (should now be the accumulated over-production from G1/G4) + the 2 ST-only rustle kills wrongly.

- [ ] **Step 2: Trace the kills** — `pred_kill` events both tools at the divergent chains (the worked example: `included_drop` killer_cov vs `retained_intron` victim_cov). Identify where rustle's kill rule/threshold differs from ST's.

- [ ] **Step 3: Port the filter decision behind the gate** — make rustle's pairwise/RI/isofrac/readthr decisions match ST's exact logic + thresholds by default; gate rustle's stricter/leaner variant behind `precise_mode()`.

- [ ] **Step 4: Build + verify** — `check_precise.sh` (OK), `check.sh` → **0 / 0** (rustle-only reclaimed, ST-only kills fixed).

- [ ] **Step 5: Commit** — `git commit -m "st-faithful(G3): match ST pred_filter cascade by default"`

---

### Task 5: Milestone gate — mini3 = 0/0, escape-hatch intact

- [ ] **Step 1:** `bash bench/mini3/check.sh` → `Rustle-only: 0  ST-only: 0`.
- [ ] **Step 2:** `bash bench/mini3/check_precise.sh` → ESCAPE-HATCH OK.
- [ ] **Step 3:** Run the existing rustle test suite under `RUSTLE_PRECISE=1` (precise behavior must be unbroken): `RUSTLE_PRECISE=1 cargo test --release 2>&1 | tail -20`. Fix any test that asserted precise behavior to set `RUSTLE_PRECISE=1`.
- [ ] **Step 4:** Quick full-BAM smoke (not full convergence — that is Milestone 2): `target/release/rustle -L GGO_19.bam -o /tmp/full_ru.gtf 2>/dev/null; python3 bench/gtf_chain_diff.py /tmp/full_ru.gtf /tmp/st_all.gtf | tail -1` — record the full-BAM rustle-vs-ST counts (expected: improved vs 187/104 but not yet 0; Milestone 2 closes it).
- [ ] **Step 5: Commit** the milestone record + update the spec's status line.

```bash
git commit -am "st-faithful: Milestone 1 — mini3 converged to 0/0; full-BAM = <counts>"
```

---

## Notes for the executor

- **Every port task gates strict behind `precise_mode()`** and must keep `check_precise.sh` green — that is the non-negotiable invariant.
- **Transient regression is expected and accepted** (G1 especially). Do NOT add compensating knobs to mask it (`feedback_full_stringtie_mimicry`); converge it in the coupled group (G3).
- **Measure rustle-vs-StringTie, never vs annotation**, for this work.
- ⚠ Never `pkill -f rustle` (cwd path contains "rustle"); use `pgrep -af target/release/rustle`.
- The 0.1s mini3 loop is the workhorse; only smoke-test full-BAM at the milestone gate.
