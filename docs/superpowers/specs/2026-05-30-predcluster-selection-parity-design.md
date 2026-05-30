# Design: StringTie-faithful predcluster + merge SELECTION parity (sub-project 1)

Status: DESIGN (approved 2026-05-30). Next: implementation plan (writing-plans).
Grounding: the 2026-05-29 predcluster divergence map (memory `project_jfp_missr_characterization`,
`docs/STRINGTIE_PARITY_FINDINGS.md`). Diff tools: `bench/gtf_chain_diff.py`, the parity-decisions
JSONL (`pred_kill`, `path_extracted`).

## 1. Problem & insight

Rustle's de-novo `-L` output is 96.5 Sn / 91.2 Pr (1750/64/169) vs the reference `GGO_19.gtf`. The
*current* StringTie binary scores 96.7 / 96.3 (1754/60/67) against the same reference — so Rustle's
real gap is **+102 false positives** (and +4 FN). The 2026-05-29 map showed **147 divergent-survivor
clusters account for 162/169 FP + 61/64 FN**, and split the error budget into two distinct subsystems:

- **SELECTION** (this spec): of the FP, **81 are chains ST extracts then KILLS** (Type A — predcluster
  decision); of the FN, **29 are ref chains Rustle extracts but silently drops** (Type B — merge/selection).
  Total **110 selection errors**.
- **CANDIDATE EXTRACTION** (separate sub-project, OUT of scope): 70 FP ST never extracts + 26 FN ST
  extracts but Rustle doesn't (~96 errors). ST's selection cannot kill a candidate it never had.

**Insight (why a coherent port, not piecemeal):** the selection rules — retained-intron containment,
isofrac, merge-collapse — are **load-bearing and Sn↔Pr-coupled**. The map proved every isolated
threshold change trades ~1:1 (stricter RI: −18 FP / −29 TP; RI off: +276 FP). The `PER_MAXINT` isofrac
windowing toggle is a literal no-op. The rules only produce ST's decisions when run **together, in ST's
order, on ST-equivalent state** — exactly as ST does. So we port the whole selection stage as one unit.

## 2. Goal & non-goals

- **Goal:** under `RUSTLE_PREDCLUSTER_ST=1`, Rustle's per-cluster winner-selection reproduces the current
  StringTie binary's selection decisions, driving the **110 selection errors → 0** on clusters where the
  candidate sets already match. Expected partial F1 gain ≈ +3pp Pr (the selection half). If it beats
  baseline F1, it becomes the default.
- **Non-goals:** candidate-extraction parity (the 70 path-enum FP + 26 FN — sub-project 2); guided/nascent
  modes; bit-exact byte reproduction (we target the *decision*, i.e. which chains survive, not internal
  byte state). Full output parity / the remaining ~+2pp needs sub-project 2.

## 3. Architecture — parallel `_st` selection function (dispatch, not piecemeal edits)

New module `src/rustle/predcluster_st.rs` exposing
`select_predictions_st(candidates, graph, …) -> kept_set`. It takes the **candidate predictions Rustle's
flow already produces** (the same input as today's predcluster) and applies StringTie's selection
sub-stages **in ST's order**. Gated on `RUSTLE_PREDCLUSTER_ST=1` (default OFF during development; flip to
default if it strictly beats baseline F1). Mirrors Rustle's established `_st` convention
(`parse_trflong_st`, `junction_graph_st`, `long_max_flow_st`). Same input/output contract as the existing
predcluster → a pure dispatch swap at the predcluster call site (transcript_filter.rs
`print_predcluster_with_summary_multi`).

## 4. Components

1. **Input adapter** — map Rustle's candidate `Transcript`/prediction to the fields ST's `CPrediction`
   selection needs: `cov`, `longcov`, intron-pattern bitvec (over graph nodes), node list, `longstart`/
   `longend`, terminal `hardstart`/`hardend`, `guide` flag. All already present in Rustle (the shipped
   `compatible_long` uses this convention). The adapter is the boundary; the `_st` functions consume only it.
2. **`pairwise_containment_st`** — port ST rlink.cpp:7363–7404 + `has_retained_intron` (7191–7201):
   for prediction pair (i⊆j by intron bitvec): if `has_retained_intron(j,i)` and `cov_i<cov_j` → kill i;
   else (contained, no RI) → kill i UNCONDITIONALLY; plus `included_drop` (rlink.cpp:18438 DROP logic).
   Emits `pred_kill` reason=`retained_intron` / `included_drop`.
3. **`isofrac_st`** — port ST rlink.cpp:18734–18794 per-maxint `longunder` loop: iterate maximal-coverage
   intervals, seed `usedcov[s]` from the single dominant prediction through each interval, kill
   `cov < isofraclong·usedcov[s]` (relaxed `isofraclong` for high-cov via CHI_WIN), accumulate
   `usedcov[s] += cov` per interval. Emits `pred_kill` reason=`isofrac`.
4. **`merge_transfrags_st`** — port ST's compatible-prediction collapse (`merge_transfrags`). Rustle's
   current merge emits no decision event; add a `pred_merge` parity event (kept-rep + absorbed members +
   cov) so it is validatable, mirroring the shipped `transfrag_collapse` event.
5. **`select_predictions_st()` orchestrator** — runs (4) → (2) → (3) in ST's order on the candidate set,
   returns the kept winners. Dispatched when the flag is on.

## 5. Data flow

flow extraction → candidate predictions → **[flag ON: `select_predictions_st`]** / [flag OFF: existing
predcluster] → kept set → final emit. Identical contract; the only change is the dispatch swap. The
adapter isolates ST-faithful logic from Rustle's `Transcript` internals.

## 6. Validation (per-stage parity gates; not unit tests for the integration — graph/selection behavior)

Validate on **candidate-matching clusters only** (where both tools extracted the same candidate set), so
selection divergence is isolated from extraction divergence:
- After each sub-stage: `pred_kill` parity diff (Rustle vs ST, by reason + intron-chain) → that reason's
  kills match ST.
- Merge stage: the new `pred_merge` event diff → same reps/absorbed members as ST.
- Combined gate: `bench/gtf_chain_diff.py(Rustle, ST)` restricted to the selection-divergent clusters →
  the 110 selection errors → 0 (or accounted).
- F1 vs reference measured at the END (and only partially closeable here — sub-project 2 closes the rest).
The `RUSTLE_PREDCLUSTER_ST` predicate itself is unit-testable (default-off); the sub-functions are
unit-testable against ST's documented logic on small candidate sets.

## 7. Build order & increments (oracle-first)

1. **Oracle (de-risk first).** Feed ST's *actual kept-set* (from ST's final GTF) into Rustle as the
   selection result on candidate-matching clusters; measure the F1/selection-cluster ceiling. **Abort
   criterion:** if forcing ST's winners does NOT close the 110 selection errors on candidate-matching
   clusters, the errors are not purely selection (e.g. coverage/merge artifacts upstream) → replan before
   building. (Proven methodology — killed multiple false starts this project.)
2. Input adapter + `predcluster_st.rs` skeleton + `RUSTLE_PREDCLUSTER_ST` flag + dispatch (pass-through
   first: with the flag on but empty logic, fall back to keeping all candidates; verify default unchanged).
3. Port `pairwise_containment_st` (RI + included_drop) → `pred_kill` RI/included parity gate.
4. Port `isofrac_st` (per-maxint) → `pred_kill` isofrac parity gate.
5. Port `merge_transfrags_st` → `pred_merge` parity gate.
6. Combined: selection-cluster `gtf_chain_diff` → 0; measure F1 vs reference; if F1 > baseline (93.78) with
   Pr up and Sn not hurt, flip `RUSTLE_PREDCLUSTER_ST` to default-on (opt-out). Else keep opt-in + report.

## 8. Safety & exit

- Default-OFF flag → zero risk to the shipped 96.5/91.24 baseline; verified by a default-unchanged check
  at each step (as in every prior layer).
- Abort/re-scope at step 1 (oracle) if the decomposition is wrong.
- The candidate-set mismatch CAPS achievable parity — success is selection-cluster convergence, not full
  output parity. The remaining gap is explicitly sub-project 2 (candidate-extraction parity).
- Multi-session effort; each sub-stage lands with its own parity gate.
