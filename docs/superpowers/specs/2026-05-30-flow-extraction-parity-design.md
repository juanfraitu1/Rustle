# Design: StringTie-faithful flow-extraction parity (sub-project 2)

Status: DESIGN (approved 2026-05-30). Next: implementation plan (writing-plans).
Builds on: sub-project 1 (`RUSTLE_PREDCLUSTER_ST`, committed opt-in selection port,
spec `2026-05-30-predcluster-selection-parity-design.md`).
Grounding: this-session predcluster map + the sub-project-1 result (memory
`project_jfp_missr_characterization`); coverage findings (`project_coverage_metrics_deviation`).

## 1. Problem & insight

Sub-project 1 ported StringTie's prediction SELECTION (`RUSTLE_PREDCLUSTER_ST`) faithfully — Pr rose
(91.6→91.8) and FP fell (168→155), but it was **F1-neutral standalone** (93.65 vs baseline 93.78)
because it runs on Rustle's flow outputs, which diverge from ST in two ways:

- **Coverage allocation:** on contested minority chains Rustle's `long_max_flow` gives lower cov than
  ST (e.g. ST 40 vs Rustle 24.6). ST-faithful isofrac/RI thresholds therefore correctly over-kill ~42
  TPs in Rustle (their cov falls below threshold) that ST keeps (its cov clears it). This is the cap on
  sub-project 1.
- **Candidate enumeration:** 70 FP are paths Rustle's flow extracts that ST never builds; 26 FN are
  paths ST extracts that Rustle doesn't.

The graph/bundle structure is ~98% identical (3351/3405 byte-identical spans), so the divergence is
NOT graph construction — it's the **flow extraction** (which paths, and what coverage) on top of the
(mostly shared) graph. The oracle (force ST winners on candidate-matching clusters) measured the
selection ceiling at **F1 95.27 (+1.49pp)**, reachable via coverage parity alone; closing the 70 FP /
26 FN goes beyond, toward ST's 96.3 Pr.

**Insight:** the back half of the pipeline only reproduces ST's output when BOTH selection (sub-project
1) AND extraction (this sub-project) are ST-faithful. Sub-project 1 is the selector; sub-project 2
gives it ST's inputs (candidate set + coverage). Together → ST's output.

## 2. Goal & non-goals

- **Goal:** under `RUSTLE_FLOW_ST=1` (default OFF, composes with `RUSTLE_PREDCLUSTER_ST`), Rustle's
  long-read flow extraction reproduces ST's candidate set and per-path coverage, so that with the
  sub-project-1 selection the final output approaches ST's (target Pr → ~96, F1 → toward 95.3 at the
  2a milestone and toward 96.3 at 2b). Flip to default if the combined `RUSTLE_FLOW_ST` +
  `RUSTLE_PREDCLUSTER_ST` beats baseline F1 (93.78).
- **Phasing (internal):** 2a coverage-allocation parity first (cashes the +1.49pp ceiling with the
  sub-project-1 selection); 2b candidate-enumeration parity second (toward 96.3).
- **Non-goals:** graph/bundle-boundary parity for the ~54/79 divergent bundles (small; deferred unless
  it blocks 2a/2b); guided/nascent modes; rewriting Rustle's flow (keep the ~40-isoform recovery).

## 3. Architecture — complete the `_st` flow variants, composed with sub-project 1

A `RUSTLE_FLOW_ST=1` flag (default OFF) routes long-read extraction through ST-faithful variants:
`long_max_flow_st` (coverage allocation) and `parse_trflong_st` (seed order + path tracing). These
`_st` scaffolds already exist in Rustle (`long_max_flow_st`, `parse_trflong_st`, and the
`RUSTLE_PARSE_TRFLONG_ST_ORDER` order flag); this sub-project completes/corrects them to match ST.
Composes with `RUSTLE_PREDCLUSTER_ST` (sub-project 1) so the back half — extraction → selection — is
ST-faithful end-to-end. Mirrors the established `_st` dispatch pattern.

## 4. Components

### Phase 2a — coverage-allocation parity (`long_max_flow_st`)
Port ST's `long_max_flow` abundance-subtraction depletion faithfully (rlink.cpp:8471, tail 8627–8665):
after extracting a path, subtract the allocated flow from each participating transfrag's abundance
(clamp to 0 below DBL_ERROR), so the residual coverage of competing paths matches ST's. The target is
the cov *ratio* of contested minorities to their dominant (ST 40/dom vs Rustle 24.6/dom) — raise
Rustle's to ST's so the sub-project-1 isofrac/RI thresholds stop over-killing the ~42 TPs. NOTE: prior
*formula*-level coverage experiments (longcov/tlen, raw_flow) were all no-ops — the divergence is the
**depletion/allocation algorithm among competing paths**, not the cov formula. This phase targets that.

### Phase 2b — candidate-enumeration parity (`parse_trflong_st`)
1. **Seed order:** ASC/reverse iteration (`RUSTLE_PARSE_TRFLONG_ST_ORDER`, exists). It failed standalone
   (−33 TP) *because* Rustle's selection picked wrong winners under minority-first; with sub-project-1's
   ST-faithful selection now available, re-evaluate it composed (it may now help).
2. **Path tracing:** `back_to_source_fast_long` / `fwd_to_sink_fast_long` (path_extract.rs:4913/3943) —
   match ST's source/sink linking so seeds are not over-extended into the dominant TSS (the alt-TSS
   collapse), eliminating the divergent candidates (70 FP / 26 FN).

## 5. Data flow

reads → graph (≈identical) → **[`RUSTLE_FLOW_ST`: `long_max_flow_st` + `parse_trflong_st`]** →
candidate set + per-path coverage → **[`RUSTLE_PREDCLUSTER_ST`: `select_predictions_st`]** → final.

## 6. Validation (parity gates; F1 only at phase ends)

- **Coverage parity (2a):** capture per-path cov (`path_extracted`) both tools; on matched chains —
  especially the 42 contested minorities — drive the cov ratio → 1 (or the cov-to-dominant ratio → ST's).
  Then re-run with `RUSTLE_PREDCLUSTER_ST` ON and confirm the 42 TPs are no longer over-killed → F1
  toward 95.27.
- **Candidate-set parity (2b):** `path_extracted` chain-set diff (Rustle vs ST) → 70 FP / 26 FN → 0
  (or accounted). Then combined F1 → toward 96.3.
- **Default-unchanged guard** at every step (flags default OFF).
- **Oracle-first:** feed ST's candidate set + coverage into Rustle (override), run with the
  sub-project-1 selection, confirm combined F1 reaches ~96.3 — bounding the full prize before building.
  Abort criterion: if forcing ST's candidates+coverage + ST-selection does NOT reach ~ST's F1, a
  divergence remains outside this scope (e.g. the graph-boundary minority, or post-selection emit) →
  re-scope.

## 7. Build order & increments

1. **Oracle** (de-risk): force ST candidates+coverage + `RUSTLE_PREDCLUSTER_ST`; measure combined F1
   ceiling vs ST's 96.3. Abort/re-scope if it falls short.
2. **Phase 2a:** complete `long_max_flow_st` depletion → coverage parity on the 42 contested minorities
   → with selection ON, F1 toward 95.27. Milestone: the +1.49pp cashed.
3. **Phase 2b-i:** re-evaluate ASC seed order composed with `RUSTLE_PREDCLUSTER_ST`.
4. **Phase 2b-ii:** `back_to_source`/`fwd_to_sink` tracing parity → candidate-set parity → F1 toward 96.3.
5. **Combined validation + default decision:** if `RUSTLE_FLOW_ST`+`RUSTLE_PREDCLUSTER_ST` beats baseline
   F1 with Sn not hurt, flip both to default-on (opt-out). Record.

## 8. Safety & exit

- Flags default-OFF → zero risk to the shipped 96.5/91.2 baseline; verified each step.
- Oracle abort at step 1 if the combined ceiling falls short of ST (a divergence outside scope).
- 2a is the high-value, better-defined phase (cashes the measured +1.49pp); 2b is the riskier
  enumeration phase (no intrinsic discriminator standalone — relies on composition with sub-project 1).
- Multi-session; each phase lands with its own parity gate. Transient regressions accepted (parity per
  phase, F1 at phase ends).
- This sub-project + sub-project 1 together are the "match StringTie's output" goal; if both land and
  beat baseline, the combined default flip realizes the full prize.
