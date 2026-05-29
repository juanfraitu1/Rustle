# StringTie Shadow Mode — Layer 4 (Flow Depletion Parity) Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans. Steps use checkbox (`- [ ]`) syntax.

**Goal:** Under `RUSTLE_ST_SHADOW=1`, make Rustle's per-seed flow extraction deplete transfrag abundance like StringTie's `long_max_flow`, so that after a dominant path is extracted, sibling seeds over shared nodes return flux≈0 and are NOT stored — driving the final-output `bench/gtf_chain_diff.py` Rustle-only chains from 1506 toward 0 (92% = 1387 are this flow divergence).

**Architecture:** Dispatch on `crate::stringtie_parity::st_shadow()` (default OFF — default 96.5/90.8 untouched). The change is in Rustle's long-read seed-extraction loop (`src/rustle/path_extract.rs` main loop ~6507, with `src/rustle/max_flow.rs` / `src/rustle/global_flow.rs`): replace/strengthen the per-seed depletion to mirror ST `long_max_flow` (rlink.cpp:8471, abundance-subtraction tail 8627-8665), plus make the checktrf-rescue store gate ST-restrictive.

**Tech Stack:** Rust (`cargo build --release`), instrumented ST fork (`./tools/stringtie/stringtie`), Python 3 (`bench/gtf_chain_diff.py` — the HONEST gate; `bench/transfrag_parity_diff.py` is a pre-collapse proxy, do NOT use it as the gate), `gffcompare` for the default regression guard.

**Reference:** findings `docs/STRINGTIE_PARITY_FINDINGS.md` §6i (Layer 4 localization) + §6h (why final-output is the gate). ST: `long_max_flow` rlink.cpp:8471 (depletion tail 8627-8665), `parse_trflong` 9693, store gate 9807/9917, checktrf 9975. Rustle: `parse_trflong` path_extract.rs:5863, main extraction loop ~6507, checktrf store 10248, checktrf gates 10146/10162.

**Key facts:**
- The default pipeline path is the long-read seed loop in `path_extract.rs` (parse_trflong → per-seed flow). `RUSTLE_ST_SHADOW=1` enables shadow; Layers 1+2 already active under it.
- HONEST gate baseline (shadow ON, before this layer): `python3 bench/gtf_chain_diff.py /tmp/ru_final.gtf /tmp/st_final.gtf` → Rustle-only **1506** (1387 flow-div, 113 filter-div). Default (shadow OFF): 219 Rustle-only, 96.5/90.8 vs truth.
- ST depletion essence (rlink.cpp:8627-8665): for each node on the extracted path, for each transfrag starting there with `flow[n1][n2]>0`: subtract `min(flow, abundance)` from `transfrag->abundance` via `update_capacity`, set `abundance=0` if `<DBL_ERROR`. Net effect: a path's transfrags are consumed, so siblings sharing those transfrags get flux 0 next.
- Default Intron-chain figure is run-to-run NONDETERMINISTIC ±0.1pp (90.7↔90.8) — not a regression.

---

### Task 1: Characterize Rustle's current per-seed depletion vs ST's

**Files:** none (investigation; output is a findings note)

- [ ] **Step 1:** Read Rustle's long-read seed extraction loop: `path_extract.rs` ~6507-6700 (the `local_nodecov` / `original_abundances` depletion), `max_flow.rs`, `global_flow.rs`. Identify exactly HOW Rustle reduces transfrag/node abundance after extracting a seed's path, and the store decision (what flux/cov threshold causes a seed to be stored vs demoted to checktrf).
- [ ] **Step 2:** Read ST `long_max_flow` (rlink.cpp:8471-8674), focusing on the abundance-subtraction tail (8627-8665) and `update_capacity`. Identify the precise difference: does Rustle deplete per-junction edgecov (weaker) instead of per-transfrag abundance-to-zero? Does Rustle's store gate accept flux that ST would have zeroed?
- [ ] **Step 3:** Write the difference (the specific depletion-strength and store-gate divergence, file:line both tools) to a `docs/STRINGTIE_PARITY_FINDINGS.md` note (§6j) and commit (`docs(shadow): Layer 4 depletion divergence characterization`). This grounds Task 2's exact change. Confirmed-known: `RUSTLE_NO_EDGECOV_DEPL=1` only moved 1506→1505, so the lever is the depletion magnitude / store gate, NOT the edgecov toggle.

---

### Task 2: Port ST's transfrag-abundance depletion under `st_shadow()`

**Files:** Modify `src/rustle/path_extract.rs` (and/or `max_flow.rs`/`global_flow.rs` per Task 1).

**Context:** After a seed's path is extracted, deplete every participating transfrag's abundance by the allocated flow (clamp to 0 below the long-read epsilon), mirroring ST rlink.cpp:8627-8665, so sibling seeds over shared transfrags get flux≈0. Gate on `crate::stringtie_parity::st_shadow()`; default path unchanged.

- [ ] **Step 1:** At the depletion site found in Task 1, add a shadow branch that, after extracting a path, subtracts the path's allocated flow from each participating transfrag's `abundance` (using the existing flow/abundance structures), setting `abundance = 0` when below the long-read epsilon (Rustle's analog of `DBL_ERROR`). Mirror ST's `min(flow, abundance)` two-case logic (8639-8654). Use the real Rustle variable/struct names from Task 1; keep the non-shadow path byte-identical.
- [ ] **Step 2:** Build: `cargo build --release 2>&1 | grep -E "^error" | head; echo "exit=${PIPESTATUS[0]}"` → exit=0.
- [ ] **Step 3:** Default regression guard (±0.1pp nondeterministic): `./target/release/rustle GGO_19.bam -L -o /tmp/shadow_off.gtf 2>/dev/null; gffcompare -r ../GGO_19.gtf -R -Q -o /tmp/shadow_off /tmp/shadow_off.gtf 2>/dev/null; grep "Intron chain level" /tmp/shadow_off.stats` → `96.5 | 90.7` or `90.8`.
- [ ] **Step 4:** Drive the HONEST gate:
```bash
RUSTLE_ST_SHADOW=1 ./target/release/rustle GGO_19.bam -L -o /tmp/ru_final.gtf 2>/dev/null
./tools/stringtie/stringtie GGO_19.bam -L -o /tmp/st_final.gtf 2>/dev/null
python3 bench/gtf_chain_diff.py /tmp/ru_final.gtf /tmp/st_final.gtf
```
Expected: Rustle-only DROPS substantially from 1506 (target: toward the ~113 filter-div floor). Record before/after. Also note ST-only (should stay near 322; a large rise means over-depletion).
- [ ] **Step 5:** Commit (`feat(shadow): Layer 4 - ST-faithful transfrag-abundance flow depletion under st_shadow()`), recording the gtf_chain_diff before/after in the message.

---

### Task 3: Make the checktrf-rescue store gate ST-restrictive under shadow

**Files:** Modify `src/rustle/path_extract.rs` (checktrf store ~10248, gates ~10146/10162).

**Context:** ST in `!mixedMode` long-read mode rarely stores checktrf rescues (rlink.cpp:9975). Rustle's checktrf gates (`RUSTLE_CHECKTRF_JUNC_GATE`:10146, `RUSTLE_CHECKTRF_RI_SUPPRESS`:10162) are default-OFF → +310 ru-only. Under shadow, apply the restrictive store gate.

- [ ] **Step 1:** Read the checktrf store path (path_extract.rs ~10100-10260) and ST's checktrf store (rlink.cpp:9975-9980). Identify ST's store condition in `!mixedMode` long-read mode.
- [ ] **Step 2:** Under `crate::stringtie_parity::st_shadow()`, gate the checktrf store to ST's restrictive condition (only store source/sink-linked rescues per ST). Keep default unchanged.
- [ ] **Step 3:** Build + default regression guard (as Task 2 Steps 2-3).
- [ ] **Step 4:** Re-run `bench/gtf_chain_diff.py` — Rustle-only should drop further (toward the filter-div floor). Record.
- [ ] **Step 5:** Commit (`feat(shadow): Layer 4 - ST-restrictive checktrf store gate under st_shadow()`).

---

### Task 4: Validate Layer 4 and document

**Files:** Modify `docs/STRINGTIE_PARITY_FINDINGS.md`.

- [ ] **Step 1:** Final gate run (shadow ON) with gtf_chain_diff; record Rustle-only / ST-only vs the 1506 baseline. Also run `gffcompare` shadow-ON vs truth to see whether precision recovered from 62.9 (it should rise as spurious chains are suppressed — informational, NOT the gate).
- [ ] **Step 2:** Confirm default (shadow OFF) still 96.5/90.7-or-.8.
- [ ] **Step 3:** Append a "Shadow Layer 4 — DONE" note to findings (§6k): gtf_chain_diff Rustle-only 1506 → after; residual breakdown (the ~113 filter-div remain for Layer 5; any flow residual); shadow-ON precision before/after. Commit.

---

## Self-Review
- **Spec coverage:** Layer 4 flow depletion (the 1387 flow-div, 92%) → Tasks 1-2; checktrf store (+310) → Task 3; validation on the HONEST final-output gate → Task 4. Filter-div (113) is explicitly deferred to a future Layer 5. Coverage/abundance divergence deferred to Layer 6.
- **Placeholder note:** Task 2's exact Rust is parameterized on Task 1's characterization (the depletion site/structures are read first) — a research-port dependency, not a placeholder; the ST source-of-truth (rlink.cpp:8627-8665) and the gate/commands are concrete.
- **Gate discipline:** the gate is `bench/gtf_chain_diff.py` (final output), NOT transfrag_pre_depl (pre-collapse proxy). Default-OFF regression guard on every code task.
- **Risk:** over-depletion would raise ST-only (Rustle missing ST chains) — watch ST-only at each step; if it rises sharply, the depletion is too strong (mismatch vs ST's two-case min(flow,abundance)).
