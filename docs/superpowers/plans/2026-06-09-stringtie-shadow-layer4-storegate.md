# StringTie Shadow Mode — Layer 4 (flow STORE GATE) Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Under `RUSTLE_ST_SHADOW=1`, port StringTie's flow STORE GATE (reject post-depletion low-flux/independent-abundance sibling seed-paths) + ST-restrictive checktrf store, driving the HONEST final-output gate (`bench/gtf_chain_diff.py`) rustle-only chains from ~1506 toward the ~113 filter floor — with a pre-registered go/no-go fork.

**Architecture:** A `crate::stringtie_parity::st_shadow()`-gated change in rustle's long-read seed store decision (`src/rustle/path_extract.rs` parse_trflong main loop). Research-first: Task 1 names the exact store-predicate divergence (rustle stores vs ST rejects), Tasks 2-3 port it, Task 4 validates at the final-output gate and applies the fork. Default (shadow off) byte-identical throughout.

**Tech Stack:** Rust (`cargo build --release`), instrumented ST fork (`./tools/stringtie/stringtie`), Python 3 (`bench/gtf_chain_diff.py` — the gate), `gffcompare` (default regression guard). Parity dataset: `GGO_19.bam` (repo root) vs `../GGO_19.gtf` (the established chr19 shadow baseline).

**Spec:** `docs/superpowers/specs/2026-06-09-stringtie-shadow-layer4-storegate-design.md`.

---

## Context the engineer needs

- **Working dir:** `/mnt/c/Users/jfris/Desktop/Rustle`. Branch `vg/flow-capacity-apportionment` (do NOT switch).
- **What's already true:** `RUSTLE_ST_SHADOW=1` (`crate::stringtie_parity::st_shadow()`) enables Layers 1+2 (junction-acceptance parity). The prior Layer-4 attempt ported ST's transfrag-abundance *depletion* and found it **gate-neutral** (1506→1505) — rustle already has depletion; the surviving siblings carry **independent abundance** the *store gate* should reject. So this layer targets the STORE decision, not depletion.
- **Rustle store sites:** `path_extract.rs` parse_trflong main loop (the seed store-vs-defer decision); markers already in code: zero-flux candidates deferred to checktrf (path_extract.rs:577, :8880), "strict behavior: failed direct long-rec seeds deferred" (:6355). The `StoreReason` doc-enum is near :569-577. The exact store predicate is what Task 1 pins.
- **ST source-of-truth:** `tools/stringtie/rlink.cpp` `parse_trflong` flow loop — the store gate (after `long_max_flow`, the condition that stores an extracted path vs discards it) and the checktrf store (`!mixedMode` long-read mode rarely stores). NOTE: memory's line numbers (9807/9917/9975) are ~11 days stale; locate by content (`store`, `nodeflux`, `checktrf`).
- **The gate** (HONEST, final output): `python3 bench/gtf_chain_diff.py <rustle.gtf> <st.gtf>` → rustle-only / ST-only intron-chain counts. NEVER use the pre-collapse transfrag proxy.
- **Default regression guard:** shadow-OFF intron-chain must stay `96.5 / 90.7`±0.1 (±0.1 = run-to-run nondeterminism, not a leak).
- **Run cargo tests with the OOM watchdog** (per `project_genomewide_oom_protocol`); `dp_scales_to_many_reads` is `#[ignore]`d.

## File structure
- Modify `src/rustle/path_extract.rs` (the store-gate branch under `st_shadow()`; checktrf gate).
- Modify `docs/STRINGTIE_PARITY_FINDINGS.md` (Task 1 §6j characterization, Task 4 §6k result).
- No new files (a pure-predicate unit test goes in path_extract.rs's test module if extractable).

---

## Task 1: Characterize the store-gate divergence (research — output grounds Task 2)

**Files:** Modify `docs/STRINGTIE_PARITY_FINDINGS.md` (add §6j).

- [ ] **Step 1: Re-verify the gate on CURRENT code.**
```bash
RUSTLE_ST_SHADOW=1 ./target/release/rustle GGO_19.bam -L -o /tmp/ru_final.gtf 2>/dev/null
./tools/stringtie/stringtie GGO_19.bam -L -o /tmp/st_final.gtf 2>/dev/null
python3 bench/gtf_chain_diff.py /tmp/ru_final.gtf /tmp/st_final.gtf
```
Record rustle-only and ST-only. Expected ≈ rustle-only 1506 / ST-only ~322 (if materially different, the 11-day-old baseline shifted — note the current numbers; they become this layer's baseline).

- [ ] **Step 2: Confirm the flow-divergence share.** Re-run with parity logging and split rustle-only chains into flow-divergence (ST never extracts the path) vs filter-divergence (ST extracts then kills):
```bash
RUSTLE_ST_SHADOW=1 RUSTLE_PARITY_LOG=/tmp/ru.jsonl RUSTLE_PARITY_FILTER_STEPS=path_extracted,pred_kill ./target/release/rustle GGO_19.bam -L -o /tmp/ru_final.gtf 2>/dev/null
STRINGTIE_PARITY_LOG=/tmp/st.jsonl STRINGTIE_PARITY_FILTER_STEPS=path_extracted,pred_kill ./tools/stringtie/stringtie GGO_19.bam -L -o /tmp/st_final.gtf 2>/dev/null
```
For each rustle-only final chain: present in ST `path_extracted`? If not → flow-divergence; if yes but ST `pred_kill`s it → filter-divergence. Record the split (expected ~92% flow). Pick 2-3 flow-divergence chains as trace cases.

- [ ] **Step 3: Read both store decisions.** In `src/rustle/path_extract.rs`, find the parse_trflong main-loop store decision (where an extracted seed-path is stored as a transcript vs deferred/discarded) — the actual condition on flux/abundance after `long_max_flow`. In `tools/stringtie/rlink.cpp` parse_trflong, find ST's corresponding store gate (after `long_max_flow`; the condition on `nodeflux`/path flux that stores vs discards) and the checktrf store. For ONE trace-case sibling chain, determine: what flux/abundance value does rustle compute for it (it stores), and what is ST's store condition that rejects it.

- [ ] **Step 4: Write §6j.** Append to `docs/STRINGTIE_PARITY_FINDINGS.md` a "§6j — Layer 4 store-gate divergence" note: the exact rustle store site (file:line), the exact ST store condition (rlink.cpp:line), and the named divergence (e.g. "rustle stores a path if flux>0; ST additionally requires path flux ≥ X / the path's transfrags weren't already consumed — siblings with independent abundance pass rustle's gate but fail ST's"). This is the concrete contract Task 2 implements.

- [ ] **Step 5: Commit.**
```bash
git add docs/STRINGTIE_PARITY_FINDINGS.md
git commit -m "docs(shadow): Layer 4 store-gate divergence characterization (§6j)"
```

---

## Task 2: Port ST's store gate under `st_shadow()`

**Files:** Modify `src/rustle/path_extract.rs` (the store site named in Task 1 §6j; + test module if predicate is extractable).

**Context:** Implement exactly the ST store condition characterized in §6j, gated on `crate::stringtie_parity::st_shadow()`. The non-shadow path must be byte-identical.

- [ ] **Step 1 (TDD, if predicate extractable): Write the failing test.** If the store decision can be expressed as a pure predicate `fn st_store_path(path_flux: f64, ...) -> bool` (inputs named in §6j), add to path_extract.rs's `#[cfg(test)] mod tests`:
```rust
#[test]
fn st_store_gate_rejects_low_flux_sibling() {
    // ST stores only when the path clears its store condition (from §6j); a sibling whose
    // post-depletion flux is below ST's threshold is rejected.
    assert!(!super::st_store_path(/* low-flux sibling args from §6j */));
    assert!(super::st_store_path(/* dominant-path args from §6j */));
}
```
(If the decision is too entangled to extract purely, SKIP this step and rely on the integration gate in Task 4 — note that in the commit. Do NOT fabricate a pure function that doesn't match the real code.)

- [ ] **Step 2: Run the test, expect FAIL** (`cargo test --lib st_store_gate_rejects_low_flux_sibling 2>&1 | tail -8`) — function not defined.

- [ ] **Step 3: Implement the shadow store gate.** At the store site (§6j), add:
```rust
if crate::stringtie_parity::st_shadow() {
    // ST store gate (rlink.cpp store condition from §6j): reject this extracted path unless it
    // clears ST's flux/abundance store threshold, mirroring ST discarding flux<=epsilon siblings
    // after the dominant path consumes shared transfrags. Use the real flux/abundance values
    // computed at this site (named in §6j).
    if !st_store_path(/* the §6j inputs */) {
        // defer to checktrf / discard exactly as the non-store branch already does here
        // (do NOT store this seed-path as a transcript)
    }
}
```
Use the real rustle variable/struct names from Task 1. Keep the non-shadow code path unchanged. If a pure `st_store_path` was added in Step 1, call it; otherwise inline the §6j condition.

- [ ] **Step 4: Build.** `cargo build --release 2>&1 | grep -E "^error" | head; echo "exit=${PIPESTATUS[0]}"` → exit=0. If a Step-1 test exists: `cargo test --lib st_store_gate_rejects_low_flux_sibling 2>&1 | tail -6` → PASS.

- [ ] **Step 5: Default regression guard (shadow OFF must be unchanged).**
```bash
./target/release/rustle GGO_19.bam -L -o /tmp/shadow_off.gtf 2>/dev/null
gffcompare -r ../GGO_19.gtf -R -Q -o /tmp/shadow_off /tmp/shadow_off.gtf 2>/dev/null
grep "Intron chain level" /tmp/shadow_off.stats
```
Expected: `96.5 | 90.7` (or `90.8`). Any other value = the `st_shadow()` gate leaked into the default path → fix before proceeding.

- [ ] **Step 6: Drive the gate.**
```bash
RUSTLE_ST_SHADOW=1 ./target/release/rustle GGO_19.bam -L -o /tmp/ru_final.gtf 2>/dev/null
python3 bench/gtf_chain_diff.py /tmp/ru_final.gtf /tmp/st_final.gtf
```
Record rustle-only / ST-only vs Task 1's baseline. **Watch ST-only:** a sharp rise above ~322 = over-rejection (store gate too strict) → revisit the §6j condition. rustle-only should DROP toward the ~113 filter floor.

- [ ] **Step 7: Commit** (record the gtf_chain_diff before/after in the message):
```bash
git add src/rustle/path_extract.rs
git commit -m "feat(shadow): Layer 4 - ST-faithful flow store gate under st_shadow() [gate: 1506->NNN]"
```

---

## Task 3: ST-restrictive checktrf store under shadow

**Files:** Modify `src/rustle/path_extract.rs` (checktrf store gates).

**Context:** ST in `!mixedMode` long-read mode rarely stores checktrf rescues. Rustle's checktrf store gates are default-off → extra rustle-only chains. Under shadow, apply ST's restrictive condition.

- [ ] **Step 1: Read** the checktrf store path in `path_extract.rs` (search `checktrf` + the deferred-to-checktrf store) and ST's checktrf store in `rlink.cpp` (`!mixedMode` long mode). Identify ST's store condition there.

- [ ] **Step 2: Gate under shadow.** Wrap the checktrf store in `if crate::stringtie_parity::st_shadow() { /* apply ST's restrictive store condition */ } else { /* existing default behavior */ }`, so under shadow rustle only stores checktrf rescues ST would. Default unchanged. (Concrete condition from Step 1.)

- [ ] **Step 3: Build + default guard** (Task 2 Steps 4-5). Default must stay `96.5/90.7`±0.1.

- [ ] **Step 4: Gate run** (Task 2 Step 6). Record rustle-only / ST-only; should drop further toward the filter floor.

- [ ] **Step 5: Commit** `feat(shadow): Layer 4 - ST-restrictive checktrf store under st_shadow() [gate: NNN->MMM]`.

---

## Task 4: Validate + apply the go/no-go fork

**Files:** Modify `docs/STRINGTIE_PARITY_FINDINGS.md` (§6k).

- [ ] **Step 1: Final gate run** (shadow ON), record rustle-only / ST-only vs the Task-1 baseline. Also `gffcompare` shadow-ON vs `../GGO_19.gtf` for the precision number (informational — should rise from ~62.9 as spurious chains are suppressed; NOT the gate).

- [ ] **Step 2: Confirm default** (shadow OFF) still `96.5/90.7`±0.1.

- [ ] **Step 3: Apply the pre-registered fork and record §6k.** Append "§6k — Layer 4 store-gate RESULT":
  - If rustle-only dropped materially (ST-only stable): **gate MOVED** — the store gate was the lever; record the new rustle-only and recommend the combined Layers 4+5+6 port.
  - If rustle-only ≈ unchanged (store gate ported cleanly, ST-only stable): **gate-neutral** — decisive evidence the store gate alone is insufficient; recommend STOP and re-scope as the combined full-downstream unit. Update `project_st_shadow_mode` memory with the outcome.

- [ ] **Step 4: Commit** `docs(shadow): Layer 4 store-gate result (§6k) + fork decision`.

---

## Self-review notes
- **Spec coverage:** characterize (Component 1) → Task 1; store gate (Component 2) → Task 2; checktrf (Component 3) → Task 3; validate + fork (Component 4) → Task 4. Default-OFF guard + ST-only alarm + gate = gtf_chain_diff on every code task. Layers 5/6 out of scope (noted).
- **Research-port dependency (not a placeholder):** Task 2's exact predicate is the concrete OUTPUT of Task 1 §6j (named divergence + ST condition + rustle store site:line). The ST source-of-truth and the gate/commands/guard are all concrete. This mirrors the validated 2026-05-29 Layer-4 plan pattern.
- **Type consistency:** `st_store_path` is the only introduced symbol, conditional on extractability (Task 2 Step 1), used consistently in Steps 1/3; if not extractable, no symbol is introduced and the condition is inlined.
