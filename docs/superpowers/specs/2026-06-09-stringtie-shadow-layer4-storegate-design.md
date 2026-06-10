# StringTie Shadow Mode — Layer 4 (flow STORE GATE) design

**Date:** 2026-06-09
**Status:** Design (not yet implemented)
**Builds on:** `docs/superpowers/plans/2026-05-29-stringtie-shadow-layer4-flow.md` (depletion port,
attempted + gate-neutral), `docs/STRINGTIE_PARITY_FINDINGS.md`, `project_st_shadow_mode` memory.

## Problem & framing

Under `RUSTLE_ST_SHADOW=1` (Layers 1+2 active: junction-acceptance parity), rustle's final output
diverges from StringTie by ~1,506 rustle-only intron chains, **92% (1,387) = FLOW divergence** —
sibling seed-paths that ST's flow never extracts but rustle stores. The prior Layer-4 attempt
ported ST's transfrag-abundance **depletion** and found it **gate-neutral** (1,506 → 1,505):
rustle already has the depletion; the siblings survive because they carry **independent
abundance** the store decision should reject. So the un-done lever is ST's **store gate**
(`rlink.cpp:9807/9917`) — the decision that *rejects* a flux≈0/low-abundance sibling after the
dominant path is consumed — plus the ST-restrictive `checktrf` store (`rlink.cpp:9975`).

**This is the documented "downstream / parse_trflong path-selection" lever every analysis
converges on.** The decisive prior conclusion is that shadow layers **do not converge
incrementally** (each ST-faithful change is correct but gate-neutral until enough of the
downstream is faithful together; precision is catastrophic ~62.9 throughout). Therefore this
layer is scoped as a **decisive, pre-registered experiment**, not an assumed win.

## Goal

Under `st_shadow()` only (default byte-identical), port ST's flow **store gate** + restrictive
`checktrf` store, driving the HONEST final-output gate (`bench/gtf_chain_diff.py`) rustle-only
chains from ~1,506 toward the ~113 filter-divergence floor.

## Pre-registered go/no-go fork

- **Gate moves** (rustle-only drops materially, ST-only stable) → the store gate WAS the lever;
  the combined Layers 4+5+6 full-downstream port is justified — continue.
- **Gate-neutral** (store gate ports cleanly, ST-only stable, but rustle-only stays ~1,506) →
  decisive evidence the store gate alone is insufficient, confirming the "everything-at-once"
  verdict → STOP this layer and re-scope as the combined full-downstream unit.

Both outcomes are informative. Defining the fork up front is the point: this contained layer
tells us which world we're in before committing to the open-ended port.

## Architecture & components

A single `crate::stringtie_parity::st_shadow()`-gated change in rustle's long-read
seed-extraction store decision, validated only at the final-output gate. Non-shadow path
byte-identical throughout.

### Component 1 — Characterize the store divergence (research, FIRST)
- Re-verify the gate on CURRENT code (the 1,506/1,387 are 11 days old): shadow-ON rustle + ST,
  `python3 bench/gtf_chain_diff.py /tmp/ru_final.gtf /tmp/st_final.gtf`; confirm rustle-only ≈1,506
  and the flow-divergence share via the path_extracted/pred_kill split.
- Read rustle's store decision (`path_extract.rs` main loop ~6507 — where an extracted seed-path
  is stored vs demoted to checktrf) against ST's store gate (`rlink.cpp:9807/9917`). Name exactly
  WHY a flux≈0/independent-abundance sibling ST rejects, rustle stores (the predicate divergence).
- Output: `docs/STRINGTIE_PARITY_FINDINGS.md` §6j note with both file:lines. No porting before the
  divergence is named.

### Component 2 — Port ST's store gate under `st_shadow()`
- At the store site, add the shadow branch applying ST's store condition (reject a seed-path whose
  post-depletion flux/abundance is below ST's threshold), using rustle's real flow/abundance
  structures (named in Component 1). Default path unchanged.
- Where the predicate is extractable as a pure function (path flux/abundance + threshold -> store
  bool), add a Rust unit test pinning ST's two-case condition (pattern: `st_shadow_from`).

### Component 3 — ST-restrictive checktrf store
- Under `st_shadow()`, gate the checktrf rescue to ST's `!mixedMode` long-read store condition
  (the default-off gates `path_extract.rs` ~10146/10162 become shadow-implied). Default unchanged.

### Component 4 — Validate at the gate + apply the fork
- Re-run `gtf_chain_diff.py`; record rustle-only before/after AND ST-only. Apply the Section
  go/no-go fork. Append a findings note (§6k).

## Validation discipline (hard-won)
- **Gate = `bench/gtf_chain_diff.py`** (final-output rustle-only chains), NEVER the pre-collapse
  transfrag proxy (3.8x inflated).
- **F1 NOT measured mid-stack** — catastrophic by design, no signal until the downstream is
  faithful.
- **Default-OFF regression guard every code change:** shadow-off intron-chain = `96.5 / 90.7`±0.1
  (±0.1 is run-to-run nondeterminism, not a leak). Real default movement = gating leaked -> fix.
- **ST-only = over-rejection alarm:** porting a *reject* gate risks dropping chains ST keeps.
  ST-only rising sharply from baseline (~322) = store gate too strict -> back off.

## Risk & abort
The pre-registered fork IS the risk management — gate-neutral is a valid decisive result, not a
failure; we stop rather than pile on isolated layers. Concrete abort: store gate ports cleanly
(ST-only stable) but rustle-only doesn't move from ~1,506 -> gate-neutral -> re-scope combined.

## Out of scope
- Filter-divergence (~113, Layer 5) and abundance/coverage (Layer 6) — deferred unless the gate
  moves.
- Default-mode behavior — untouched (all changes `st_shadow()`-gated, default-off).
- Genome-enabled parity (consensus wiring) — separate track.
