# Design: prediction-population parity — donor-snap (2A) + chimeric-suffix-rescue suppression (2B)

Status: DESIGN (approved 2026-05-31). Next: implementation plan (writing-plans).
Grounding: needy-loci decision trace + the lever-2 mechanism verification
(`project_needy_loci_decision_trace`); the retained-intron shelve lesson (same memory); the
flow/coverage-depletion floor (`project_coverage_metrics_deviation`).

## 1. Problem & insight

The needy-loci trace named 3 over-enumeration levers; lever #1 (retained-intron) was SHELVED because it
tried a downstream filter and ST's killers don't survive in RU's final set — the divergence is the
upstream prediction POPULATION. Lever #2 ("snap/fold at transfrag-collapse") was mechanism-verified and
turned out to be TWO independent sub-levers, NEITHER in transfrag-collapse, BOTH acting on the
prediction population (the right locus):

- **2A — near-duplicate DONOR snap.** ST's `higherr` block (`rlink.cpp:14233-14373`, apply `:14535-14640`)
  moves reads on a weak donor onto a stronger same-strand donor within 25bp BEFORE any transfrag is
  built (e.g. 30800990 → 30800998, 8bp apart; ST builds 0 seeds on 30800990, RU builds 9 → emits FP
  isoforms RSTL.398.16/.11/.5). RU HAS the port (`correct_bundle_junctions_higherr`,
  `junction_correction.rs:151`, wired `junctions.rs:626` → `pipeline.rs:11443`) but the `is_bad` gate
  (`junction_correction.rs:178-189`) was narrowed to `intron>100kb AND nreads_good<10`, excluding this
  case. **Scope gap in an existing mechanism.**
- **2B — minority/partial-read "fold."** NOT a transfrag fold: RU's FP (RSTL.69.9) is created by a
  DOWNSTREAM RU-only path rescue, `chimeric_suffix_rescue` (`path_extract.rs:2992-3034`, call site
  `:9485`), which REVIVES a 5'-truncated suffix transfrag as a standalone prediction. ST has no analog —
  the suffix's coverage stays folded in the full-length flow. **Extra RU-only mechanism.** Existing
  toggle: `RUSTLE_DISABLE_CSR` (`path_extract.rs:9485`).

**Insight:** both act on the prediction population (junction stage / path-extraction), where every prior
experiment located the real divergence. They are independent (different code paths, fix types, risks).

**The hard lesson (retained-intron):** an oracle PRIZE is a ceiling, not proof a fix realizes it. Each
oracle must verify the FP disappears AND the real isoform the read belongs to is not lost. TP-cost is
the dominant abort trigger. Mechanism mis-attribution happened 3× before; this design is pinned to the
verified file:line mechanisms.

## 2. Goal & non-goals

- **Goal:** measure — oracle-first, default-OFF — whether 2B (chimeric-suffix-rescue suppression) and 2A
  (junction donor-snap) remove over-enum FPs at acceptable TP-cost; decide ship/shelve EACH from its
  oracle. Sequenced 2B-first (it has a free probe).
- **Decision rule (user-chosen):** oracle-bound first; each sub-lever independently gated; TP-cost is the
  dominant trigger.
- **Non-goals:** lever #3 (checktrf rescue-vs-drop); the flow/coverage-depletion floor; transfrag-collapse
  / `compatible_long` (verified NOT the mechanism for either sub-lever).

## 3. Architecture — two sequenced, independently-gated sub-experiments

### 2B (first — has free probe)
- **Phase 0 free probe:** `RUSTLE_DISABLE_CSR=1` → gffcompare + `gtf_chain_diff` vs reference & baseline.
  Net = full-rescue FP-removed − real-alt-TSS-TP-lost. If strongly net-negative (rescue earns its keep),
  only the contained subset is suppressible.
- **Phase 0b refine** (`bench/csr_classify.py`): classify each `chimeric_suffix_rescue` prediction as
  CONTAINED (node-chain is a strict suffix of a kept full-length path — ST folds these) vs genuine
  alt-TSS; FP-vs-TP vs reference. Targeted prize = contained-FP; targeted cost = contained-TP.
- **Abort gate:** abort 2B if targeted net (contained-FP − contained-TP) ≤ 0.
- **Phase 1 (gated):** suppress the revival when the suffix is contained in a kept full-length path,
  behind `RUSTLE_CSR_FOLD` (default OFF), at `is_chimeric_suffix_rescue` (`path_extract.rs:2992`) / call
  site (`:9485`).

### 2A (second)
- **Phase 0 oracle** (`bench/donor_snap_prize.py`): from `junction_accept` (both tools) + reference +
  RU final, find RU junctions where a stronger same-strand donor exists within 25bp (the `is_bad`-excluded
  snap candidates). Prize = RU-only-FP final chains using the weak donor; cost = RU chains whose weak-donor
  intron IS a real reference intron (genuine alt-donor the snap would destroy = small-exon regression).
  Net = prize − cost.
- **Abort gate:** abort 2A if cost ≥ prize or net < ~3.
- **Phase 1 (gated):** a NARROWED re-enable of the `is_bad` gate (`junction_correction.rs:178-189`) behind
  `RUSTLE_HIGHERR_SNAP` (default OFF) — add a support-floor / good_junc check so it snaps noise donors but
  spares well-supported real alt-donors (NOT a blanket `nm>=nreads`, which "killed small exons" — comment
  at `junction_correction.rs:164-172`). Redirect→read-repair plumbing already wired.

## 4. Components

- 2B Phase 0: existing `RUSTLE_DISABLE_CSR` + gffcompare + `bench/gtf_chain_diff.py`.
- 2B Phase 0b: `bench/csr_classify.py` (new) — reads RU parity (`path_emit`/rescue events carrying
  `source=checktrf_rescue`, `rescue_class=chimeric_suffix_rescue`; confirm the emit at `path_extract.rs:9485`)
  + RU final + reference; CONTAINED-vs-altTSS + FP/TP.
- 2B Phase 1: flag `RUSTLE_CSR_FOLD` predicate (`stringtie_parity.rs`) + suppression guard at
  `path_extract.rs:2992/9485`.
- 2A Phase 0: `bench/donor_snap_prize.py` (new) — near-dup-donor prize/cost from `junction_accept` + ref.
- 2A Phase 1: flag `RUSTLE_HIGHERR_SNAP` predicate + narrowed `is_bad` gate at `junction_correction.rs:178`.

## 5. Data flow (flags ON)

reads → **[2A: junction donor-snap]** → junctions → graph → transfrags → flow →
**[2B: suppress contained-suffix rescue]** → selection → final. Both change the prediction POPULATION.

## 6. Validation (gates)

- **2B:** free-probe full net; decisive gate = Phase-0b targeted net (contained-FP − contained-TP) > 0.
- **2A:** prize − cost; abort if cost ≥ prize or net < ~3.
- **Phase 1 (each):** genome-wide j-FP reduction + Sn vs 95.6/90.5; ship only on net F1 gain; else opt-in.
- **Realizability cross-check (the lesson):** before believing a Phase-1 gain, confirm the removed FPs
  are the oracle-predicted chains AND the prize chains' reads fold into a surviving TP (no orphaned
  coverage / new artifact).
- **Default-unchanged guard** at every step (flags default OFF).

## 7. Safety & exit

- 2B Phase 0 = existing default-OFF toggle; 2A Phase 0 = analysis-only; Phase-1 flags default OFF; default
  95.6/90.5 verified each step.
- **Dominant risks:** 2B — suppressing genuine alt-TSS isoforms (measured by contained-TP); 2A — snapping
  real NAGNAG-style alt-donors within 25bp (measured by weak-donor-real-intron count). Both surfaced by
  the oracle before code.
- Each sub-lever shelves/ships independently; 2B's free probe may short-circuit.
- **Honest expectation:** these are ~5-FP levers on 3 loci; retained-intron's +25 ceiling proved
  unrealizable. Genome-wide net may be small and/or TP-offset. An abort with a clear TP-cost reason is a
  successful result — it shows the rescue/snap is doing real work, and further localizes the floor.
