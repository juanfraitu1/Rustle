# Copy-assignment significance gate (IsoCon-style) — design

**Date:** 2026-06-26
**Status:** approved design, pre-implementation
**Scope:** one focused change to `src/rustle/vg_family/copy_assign.rs` + tests + evaluation scripts.

## Goal

Replace the fixed log-LR margin τ as the ASSIGN/ABSTAIN decision in copy-assignment with IsoCon's
**real-vs-error hypothesis test**, producing a **calibrated per-read p-value certificate** and a
power-aware identifiability condition. This is the most IsoCon-faithful improvement
([[reference_isocon_sahlin]]), is Canzar-aligned (a standard significance level, not a hand-set margin),
and closes the open L17 "abstain certificate" loose end.

Non-goals (explicitly deferred): iterative template refinement (IsoCon's other step), the Clair3-RNA
RNA-editing PSV filter, any new module/file.

## Background — current model (unchanged parts)

`assign_read(read, copies, params)` computes, per copy, a per-base-quality log-likelihood over the family's
PSV columns (match → `ln(1-e)`, mismatch → `ln(e/3)`, `e` = the read's per-column quality or a flat default)
plus a junction log-odds term (`±junction_weight`). It picks `best = argmax`, computes
`margin = logl[best] − logl[second]`, and decides: `Tied` if `n_decisive < 1`, else `Assigned` if
`margin ≥ τ`, else `Ambiguous`. τ = `ln((1−p)/p)` is a principled operating point but a *fixed* margin —
it only approximates a target per-read misassignment rate and ignores how much distinguishing *power* the
read actually has.

**Kept as-is:** the per-copy log-likelihood, the `best` selection, and the reported `log_lr_margin`.
**Changed:** the decision (gate) and the identifiability condition.

## The statistical model

For the chosen `best` copy *b* and each competitor copy *c* ≠ *b*, restrict to the **distinguishing
observations** the read spans:

- **Distinguishing PSV columns**: columns `j` where `b.alleles[j]` and `c.alleles[j]` are both defined and
  differ, and the read spans `j` (`psv_obs[j] = Some`).
- **Distinguishing junctions**: junctions in the **read's own junction set** that are present in exactly
  one of `b.junctions` / `c.junctions` (within `boundary_tol`). Restricting to junctions the read *has*
  matches the existing machinery (which iterates `read.junctions`) and avoids having to assert a read
  "lacks" a junction it never spanned. A read-junction matching `b` but not `c` is a success; matching `c`
  but not `b` is a non-success (it supports the competitor).

Each distinguishing observation *j* is a Bernoulli trial:

- **success** = "the read shows *b*'s allele / junction-pattern at *j*". Observed count `k` = number of
  successes over the `n` distinguishing observations the read spans for this pair.
- **success probability under H₀ "read came from *c*"** = εⱼ:
  - PSV column: εⱼ = eⱼ/3 (the read's per-column error landing on *b*'s specific allele; eⱼ from quality
    or the flat default), capped to `[0,1)`.
  - junction: εⱼ = `junction_err` (a small fixed splice/mapping error probability).

```
p(b,c) = P( PoissonBinomial({εⱼ : j ∈ distinguishing(b,c) ∩ read}) ≥ k )
p_read = max over c≠b of p(b,c)            # IsoCon's least-significant-pair rule
```

`p(b,c)` is the probability of seeing **this much b-support by error alone** if the read were truly *c*.
Small `p_read` ⇒ the read's evidence rejects every competitor ⇒ confident assignment.

**Empty-set conventions (these define the Tied edge case).** A competitor *c* with **no** distinguishing
observations the read spans has 0 trials: `poisson_binomial_upper_tail(0, [])  = 1.0` and the empty product
`Π ∅ = 1.0`. So such a competitor forces `p(b,c) = 1.0` and contributes `1.0` to `min_p_value` — i.e. the
read is indistinguishable from *c* and the decision routes to `Tied` (since `Tied` is checked first, before
`Assigned`). `poisson_binomial_upper_tail(k>0, [])  = 0.0` (0 trials can't reach k>0). For `k = n` the tail
equals `Π εⱼ`, matching the `min_p_value` definition.

### Identifiability certificate (falls out for free)

The *minimum attainable* p against competitor *c* (read supports *b* at **every** distinguishing position,
`k = n`) is `Π εⱼ`. Define `min_p_value = max_c (Π_{j ∈ distinguishing(b,c)} εⱼ)` (the hardest competitor).

- **Resolvable** iff `min_p_value < α`: there exists an evidence configuration that could reject every
  competitor at level α.
- **Unresolvable** iff `min_p_value ≥ α`: some competitor has too few / too low-quality distinguishing
  positions; **no** evidence configuration can ever reach significance.

This replaces the heuristic `n_decisive ≥ 1` with a **power-calibrated** condition: not "≥1 position
differs" but "the distinguishing positions carry enough combined certainty (`Π εⱼ < α`) to *possibly*
resolve." (`n_decisive` is still computed and reported, for continuity.)

### Decision (with the Bonferroni family-wide correction)

`best` is the **argmax over n copies**, so a per-competitor p-value at level α only bounds the *pairwise*
error; the winner's-curse of selecting `best` inflates the family-wide misassignment to the union bound
`(n−1)·α`. To control the **family-wide** rate at α, each per-competitor test must clear `α/(n−1)`
(`thr = α / (n−1)`; for `n=2` this is just α, so the common two-copy case is unchanged):

```
thr    = α / (n − 1)                            # Bonferroni over the n−1 competitors `best` was chosen against
status = Tied      if min_p_value ≥ thr         # identifiability-limited (unresolvable even at full support)
       = Assigned  if p_read < thr AND margin>0 # rejects every competitor at the corrected level, and best is the strict MLE
       = Ambiguous otherwise                    # resolvable in principle, but this read's evidence is short
```

The `margin > 0` guard makes exact-LLR-tie / balanced-conflict reads abstain (the one-sided test alone would
assign the index-tiebreak winner). With the Bonferroni `thr`, the realized family-wide misassignment among
ASSIGNED reads is bounded by α — a theorem, not just the empirical 2-copy result.

## API changes (`copy_assign.rs`)

- `Assignment` gains:
  - `p_value: f64` — `p_read`, the certificate (least-significant competitor p).
  - `min_p_value: f64` — the identifiability bound (best attainable p over the hardest competitor).
  - keeps `best_copy`, `log_lr_margin`, `n_decisive`, `resolvable`, `status`. `resolvable = min_p_value < α`.
- `AssignParams` gains:
  - `alpha: f64` — the **family-wide** significance level (target misassignment rate over ASSIGNED reads,
    via the Bonferroni `α/(n−1)` correction). **Default 1e-3** (matches the production operating point τ=6.9 ≡
    p=1e-3 / Eichler-AS≥10; minimises surprise vs shipped numbers).
  - `junction_err: f64` — εⱼ for a distinguishing junction. **Default 1e-4** (strictly below the default α so
    a *single* copy-specific junction is decisive, matching the old `junction_weight`/τ behaviour; HiFi splice
    junctions are sharp). Tunable — a more conservative value (e.g. 1e-3) makes a single junction abstain,
    appropriate if RtSwitch / template-switch mis-splice false-presence is a concern.
  - `use_margin_gate: bool` — default `false`. When `true`, reverts to the legacy τ-margin decision (for
    reproducing legacy numbers and the A/B comparison). The significance test is otherwise the default gate.
  - new constructor `AssignParams::for_alpha(α)`.
- New helper `poisson_binomial_upper_tail(k: usize, probs: &[f64]) -> f64` — exact O(n²) DP convolution of
  the per-trial success probabilities, returning `P(Σ ≥ k)`. Distinguishing positions are few (<~100), so
  exact is cheap and avoids normal-approximation error.

`tau_from_p` / `p_from_tau` / `for_target_misassignment` are retained (the LLR margin is still reported and
the legacy gate still selectable).

## Integration

The change is confined to `assign_read`. Downstream callers (`copy_assign_pipeline::assign_family_detailed`,
the assigned/ambiguous/tied aggregators) branch on the `AssignStatus` enum, whose three variants are
unchanged, so they keep working — the counts simply reflect the new gate. No call-site signature changes.

## Evaluation

1. **sim5x labeled-truth ladder** (`bench/sim_reads.py`, K=0..8 PSV ladder — the load-bearing identifiability
   benchmark). Significance gate vs legacy τ gate, at matched recall:
   - accuracy of assigned reads, recall (% assigned), abstain breakdown;
   - **calibration**: as α sweeps {1e-2, 1e-3, 1e-4}, does the *realized* misassignment rate on the labeled
     reads track α? A calibrated p-value (realized error ≤ α) is the headline result — the property the
     fixed τ-margin could only approximate;
   - the certificate must flag the K=0 (exonically-identical) reads as `Tied` via `min_p_value ≥ α`.
2. **Real GGO**: re-run copy-assignment on the conflict-catalog substrate with the significance gate; report
   assigned/ambiguous/tied % and silver agreement vs the τ-gate baseline.

## Testing

- `poisson_binomial_upper_tail`: equal probs ⇒ binomial upper tail; `k=0 ⇒ 1.0`; all-zero probs with `k>0`
  ⇒ `0.0`; monotone non-increasing in `k`.
- Certificate behavior: 1 high-Q distinguishing PSV supporting `best` ⇒ small `p`; 2 ⇒ `Assigned` at
  α=1e-3; 0 distinguishing ⇒ `Tied` (`min_p_value ≥ α`); conflicting support (matches `best` at some,
  competitor at others) ⇒ `Ambiguous`.
- Calibration unit test: simulated reads at a known error rate ⇒ realized misassignment ≤ α.
- Update the handful of existing `assign_read` tests whose expected `status`/`margin` assumed the τ gate —
  move them to the significance gate, or pin `use_margin_gate = true` where the test is specifically about
  the margin.

## Risks / decisions

- **Default behavior changes** (τ-margin → significance). Accepted: significance is the new default gate;
  legacy reproducible via `use_margin_gate`. Headline assignment numbers to be recomputed on GGO.
- **Uniform-error assumption** (εⱼ = eⱼ/3): treats the 3 non-true bases as equally likely error targets.
  Conservative *in expectation*, but NOT per-substitution worst-case: a biased substitution makes it
  anti-conservative (up to 3×). The canonical case is **RNA editing (A→I read as A→G)** — at an A/G PSV the
  true-A read's rate of showing G is the editing rate, which can be ≫ eⱼ/3. The deferred Clair3-RNA editing
  filter ([[reference_clair3_rna]]) addresses exactly this; until then an unfiltered editing site at a PSV
  can be anti-conservative. A per-substitution error model is a future refinement.
- **junction_err** is a single fixed value (not per-junction). Adequate; junctions are sharp. Tunable.

## Updates applied during implementation (final-review reconciliation)

- **Bonferroni `α/(n−1)`** added to the gate so the family-wide misassignment rate is controlled at α (the
  un-corrected gate was calibrated per-competitor only; with `n` copies the winner's-curse inflates realized
  error to `(n−1)·α` — verified in `sig_gate_bonferroni_controls_familywide_error`). For `n=2` `thr=α`, so
  the validated two-copy + sim5x numbers are unchanged.
- **`margin > 0` guard** added: exact-LLR-tie / balanced-conflict reads abstain (`Ambiguous`) instead of
  going to the index-tiebreak argmax.
- **`junction_err` default = 1e-4** (the spec's earlier "~0.01" would have made a single junction abstain).
- **L17 CLOSED:** `p_value` + `min_p_value` are emitted in `<out>.assignments.tsv` (the per-read certificate).
