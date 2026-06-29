# Kill-or-validate: what actually beats the production vote-rule on copy assignment

The workflow's first recommendation was "wire the LLR (log-likelihood) engine into production." A
labeled-data decomposition (`psv_llr_vs_votes_killtest.py`, ground-truth sim copies) shows that the
**likelihood engine *per se* buys nothing** — the real lever is the **admission gate**. Three effects,
isolated on identical reads:

| effect | comparison | result |
|---|---|---|
| **SCORING** (votes → likelihood) | `VOTE_gate1` vs `LLR_flat` | **identical in 16/16 configs.** With flat per-base error the LLR margin = vote_margin·ln(3(1−e)/e), so likelihood ≡ vote-count. Switching scorers alone = 0 gain. |
| **GATE** (`min_psv=3` → `n_decisive≥1`) | `VOTE_prod` vs `VOTE_gate1` | **the dominant lever.** A family with only 1–2 PSVs can never reach 3 votes, so the production rule **discards 100% of the near-identical tail** (recall 0.000). The `n_decisive≥1` gate resolves it: +1151–2901 reads recovered at **94.7–100% precision**. |
| **QV-INFO** (flat vs calibrated per-base QV) | `LLR_flat` vs `LLR_cal` | the **only** thing votes structurally cannot do. Real but **modest and high-error-only**: at normal HiFi error it removes ~0 misassignments; under stress (err×3) it ~halves them (e.g. C5K1 97→38 wrong, prec 0.947→0.978) at a small recall cost. A no-op at realistic QV. |

## The reframed lever

The headroom over Eichler's AS≥10 gate (measured: ~28% of discarded GGO reads recoverable) AND over
Rustle's own production rule comes from the **same place** — the **1–2 PSV near-identical tail** that
both gates throw away:

- Eichler discards it because the *aggregate* alignment-score margin is <10 (signal diluted over the
  full molecule).
- Rustle-production discards it because `min_psv=3` (a 3-concordant-PSV hard floor) can't be met.

The fix is one principled gate: **`n_decisive≥1` (identifiability) + a single calibrated decisive
margin τ = ln((1−p)/p)** that sets the precision/recall operating point by a *stated* misassignment
bound P(misassign)≤p — the Canzar-clean restatement of Eichler's conservative AS≥10, with the magic
integers (`min_psv=3`, `margin=1`) collapsed into one interpretable knob.

The LLR engine (`copy_assign.rs::assign_read`) is the right **vehicle** for τ (it makes τ a probability
bound and can fold in QV + junction columns), but the **mechanism that recovers reads is the gate, not
the likelihood**. Sell it that way.

## What this does NOT cross

The `n_decisive=0` floor is untouched: reads spanning no column where candidate copies differ
(exonically identical silent copies; single-shared-exon reads) stay discarded under every arm. Every
lever here re-adjudicates the **resolvable** tail more correctly; none manufactures a decisive feature.

## Operating-point evidence (the one knob)

Resolving the 1–2 PSV tail is a precision/recall trade, governed by τ:
- K=2 families: recovered at ~99–100% precision (cheap, near-free win).
- K=1 families: recovered at ~95–98% precision (one PSV → one error can flip it; this is where τ and QV
  earn their keep). Production's `min_psv=3` is just the τ→∞ corner (discard the whole tail).

Artifacts: `psv_llr_vs_votes_killtest.py` · `psv_llr_vs_votes_killtest.png` · `.json`.
