# Design: StringTie-faithful color + mm_negative alignment (bad-junction handling)

Status: DESIGN (approved 2026-05-30). Next: implementation plan (writing-plans).
Grounding: this-session color/CGroup audit (`project_color_cgroup_parity` memory; findings §6j);
the mm_negative floor (`project_jfp_missr_characterization`; findings §6h); exact-ST philosophy
(`feedback_full_stringtie_mimicry`).

## 1. Problem & insight

Rustle does NOT produce StringTie's colors / bundle segmentation. The color audit (findings §6j)
located the root precisely, and it is NOT mm_negative as first hypothesized:

- **Half 1 — color-break test.** Production call `pipeline.rs:14819` passes `junction_stats=None`,
  so `color_should_break_left` (`bundle_builder.rs:137`) breaks on `!has_good_left` (junction not
  in `good_junctions`) instead of ST's `strand==0` (`rlink.cpp:1069`). `has_good_left` is strictly
  narrower than `strand!=0`, so Rustle over-breaks color: **2510 of 2536 RU break sites (99%) land
  on a junction ST keeps with strand ±1** → pervasive interrupted colors ST never produces.
- **Half 2 — mm_negative edge handling.** Rustle's mm_negative reject (`graph_build.rs:845`)
  filters the junction out of the graph EDGE list but never sets `JunctionStat.strand=0` (verified:
  all 3467 RU mm_negative junctions keep jstrand=±1; zero RU junctions have jstrand==0). ST's
  good_junc instead sets `strand=0` on the demoted junction and keeps it as a color-break signal.
  So Rustle is MORE aggressive on the edge (removes it) but, unlike ST, never breaks color on it.

**The entanglement.** Flipping ONLY half 1 (`RUSTLE_SUBBUNDLE_USE_STATS=1`, ST's strand==0 break)
REGRESSES parity: bundles 3405→2569 (now UNDER-segmenting vs ST 3430), mismatches 125→155. The
`None`/`!has_good_left` default is a TUNED compensation — Rustle over-breaks color (broad) AND
over-filters mm_negative edges; the two partially cancel. True faithfulness requires aligning BOTH
halves together.

**The risk.** Aligning mm_negative handling to ST lands in the known mm_negative floor (findings
§6h), which is globally F1-catastrophic. So this experiment must be ORACLE-FIRST: bound the prize
before building, and abort if the prize is structurally negligible.

**Insight (tempering the prize).** Union-find re-merges most color breaks — only 11 of 125
boundary-mismatch bundles are TRUE splits, and only ~3 FPs (RSTL.103.1 class p, RSTL.331.2 class j;
a third split benign) were traced. The reachable prize is plausibly ~3–8 net FP/FN (~0.15–0.4pp).

## 2. Goal & non-goals

- **Goal:** measure — oracle-first, default-OFF — whether aligning both halves of bad-junction
  handling to ST converges Rustle's junction set, colors, and bundle segmentation to ST's, and what
  that does to F1. Decide ship/shelve from the measured prize BEFORE building Phase 1.
- **Decision rule:** oracle-bound first, then decide (user-chosen). Hard abort gate between Phase 0
  and Phase 1.
- **Success criterion (if Phase 1 reached):** bundle-segmentation convergence to ST
  (`graphnode_diff` mismatches 125→0) + junction divergent-treatment count → 0. F1 measured and
  recorded; transient regression accepted per exact-ST philosophy.
- **Non-goals:** the ~10174 ST-only low-support accepted junctions (Layer-3 support-accounting gap,
  out of scope); flow/selection parity (sub-projects 1/2); guided/nascent modes; the committed
  `get_min_start` tie-break fix (independent, already shipped).

## 3. Architecture — four phases, hard abort gate

- **Phase 0a — Junction-set parity (foundational).** Definitively answer "do we have exactly ST's
  junctions?" at three layers; classify only-RU / only-ST / shared-same / shared-divergent.
- **Phase 0b — Prize-bound at divergent loci.** Locus-level final-GTF attribution of FP/FN to
  segmentation across the 125 mismatch bundles; counterfactual force-merge on true-split loci to
  confirm causation. Produces the net reachable FP+FN number.
- **Abort gate.** If net reachable prize < 5 net FP/FN (~0.25pp F1): STOP, record, shelve Phase 1.
  (The ~3 traced FPs suggest this is the likely outcome — and a successful, bounding result.)
- **Phase 1 (only if gate clears) — two-halves faithful alignment** behind one default-OFF flag.

## 4. Components

### Phase 0a — Junction-set parity (`bench/junction_set_diff.py`, new)
- **Inputs:** `junction_accept` JSONL from both tools (exists; we know RU 7303 ⊆ ST 17477, 0 RU-only
  at this layer), plus two new instrumentation-only emit points to see all three layers:
  1. **Raw observed** splice sites (pre-filter) — Rustle at the junction-collection site; ST at the
     matching `rlink.cpp` junction-vector build.
  2. **`good_junctions` membership** (the set `has_good_left` queries) — Rustle where `good_junctions`
     is finalized; ST where good_junc sets/keeps strand (the strand-zeroing site).
- These are env-gated parity emits (`RUSTLE_PARITY_FILTER_STEPS` / `STRINGTIE_PARITY_*`), NO logic
  change.
- **Output:** classification table — only-RU / only-ST / shared-same / shared-divergent-treatment
  (e.g. RU mm_negative-edge-filtered-but-jstrand±1 vs ST strand=0), bucketed by read-support,
  strand, mm_negative status. Feeds Phase 1 expectation (how many junctions should change treatment).

### Phase 0b — Prize-bound (`bench/segmentation_prize.py`, new)
- **Inputs:** `bench/graphnode_diff.py` output (the 125 mismatch bundles) + `bench/gtf_chain_diff.py`
  (final RU vs ST vs reference).
- **Logic:** per mismatch bundle, attribute each FP/FN to segmentation — a RU FP inside a RU-only
  split sub-bundle, or an FN where a RU split dropped a chain ST kept. Sum net FP+FN.
- **Confirmation:** counterfactual force-merge (throwaway hack) on the ≤11 true-split loci — merge
  RU's split bundles, re-measure, confirm the attributed FP disappears (causation, not co-location).
- **Output:** net reachable FP+FN → the abort number.

### Phase 1 — two-halves alignment (only if gate clears), flag `RUSTLE_ST_BADJUNC` (default OFF)
- **Half 1:** `pipeline.rs:14819` — pass real `junction_stats` so `color_should_break_left` uses
  ST's `strand==0` test (the existing `RUSTLE_SUBBUNDLE_USE_STATS` path), under the combined flag.
- **Half 2:** `graph_build.rs:845` — instead of edge-filtering mm_negative, set
  `JunctionStat.strand=0` (matching ST good_junc) so it is BOTH a color-break signal AND
  edge-handled like ST. New helper, flag-gated.
- The two halves compose under the single flag (neither alone — half 1 alone is known to regress).

## 5. Data flow (Phase 1, flag ON)

reads → junctions (mm_negative → strand=0, ST-faithful) → colors (break on strand==0, ST-faithful)
→ bundles (converge to ST segmentation) → graph → flow → selection → final GTF.

## 6. Validation (parity gates; F1 at phase ends only)

- **Phase 0a:** diagnostic — the classification IS the deliverable (no pass/fail).
- **Phase 0b:** net reachable FP+FN → **abort if < 5 (~0.25pp)**.
- **Phase 1 (if reached):** `graphnode_diff` mismatches 125→0 (or accounted) + junction
  divergent-treatment → 0 = success criterion. F1 measured/recorded; regression accepted.
- **Default-unchanged guard** at every step (flags default OFF; baseline 95.6/90.5 verified).

## 7. Safety & exit

- All instrumentation env-gated, default silent. Phase 1 behind `RUSTLE_ST_BADJUNC` (default OFF) →
  zero risk to the shipped 95.6/90.5 default; verified each step.
- **Abort at 0b (likely):** record prize ceiling + junction-set picture in findings/memory; shelve
  Phase 1. A successful bounding result — closes the color-divergence question.
- **Proceed past 0b → Phase 1:** measure convergence + F1. If F1 holds/improves → consider default
  flip; if it regresses (the floor) → keep as documented ST-faithful opt-in, record the regression
  magnitude as the cost of exact-ST color parity.
- Multi-step but small; Phase 0 is pure analysis (no production risk). Transient regressions
  accepted in Phase 1 per exact-ST philosophy (parity per stage, F1 at phase ends).
