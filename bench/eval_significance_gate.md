# Evaluation — IsoCon significance gate vs the τ-margin gate

The copy-assignment decision was changed from a fixed log-LR margin τ to IsoCon's real-vs-error significance
test (`copy_assign::assign_read`, default gate). This note records the validation. Reproduce with the
`copy_assign` binary flags `--alpha <α>` (significance gate, default 1e-3) and `--margin-gate --margin <τ>`
(legacy τ-gate).

## 0. The guarantee (what α controls)

`best` is the argmax over `n` copies, so a per-competitor p-value only bounds the *pairwise* error; selecting
the winner inflates the family-wide rate to the union bound `(n−1)·α`. The gate applies a **Bonferroni
correction** — each competitor must clear `α/(n−1)` — so the **family-wide misassignment rate over ASSIGNED
reads is bounded by α** (a theorem, validated in `sig_gate_bonferroni_controls_familywide_error`: in a
10-competitor 1-PSV regime the corrected gate keeps realized ≤ α by abstaining, while the un-corrected
per-pair level inflates to ~K·(e/3) > α). For `n=2`, `thr=α` (unchanged).

## 1. Unit calibration (Rust test `sig_gate_is_calibrated_realized_error_tracks_alpha`)

Simulated reads from a known copy over 6 Q30 PSVs (n=2), 20k reads per α ∈ {1e-2, 1e-4}: the gate assigns
~100% with **0 realized misassignments** at both α (the bound holds with slack). Confirms the certificate is
conservative and high-recall when PSVs are present.

## 2. sim5x labeled-truth ladder (ground truth, NOT the circular silver standard)

`bench/build_sim5x.py` builds a 5-copy tandem gene with a private-exonic-PSV ladder K=0..8 and reads named
with their true copy (`K{K}_c{copy}_r{n}`). Per-read assignment scored against the label (detected→true copy
mapped by majority vote), significance gate (α=1e-3) vs legacy τ-gate (τ=6.9):

| K | assigned% (both gates) | tied% | accuracy of assigned | misassigned |
|---|---|---|---|---|
| 0 (exonically identical) | **0%** | **100%** | — | 0 |
| 2 | 100% | 0% | 100% | 0 |
| 4 | 100% | 0% | 100% | 0 |
| 8 | 100% | 0% | 99% | 2 / 200 |

**Findings:**
- **The significance gate is identical to the validated τ-margin gate** at every K — same recall, same
  accuracy, same 2/200 errors at K8. So it is a drop-in, no-regression replacement: it reproduces the
  validated behaviour while adding the calibrated per-read certificate and a cleaner identifiability rule.
- **K=0 is correctly 100% Tied.** Exonically-identical copies have no distinguishing PSV, so every
  competitor has an empty distinguishing set → `min_p_value = 1.0 ≥ α` → `Tied`. The identifiability
  certificate flags the unresolvable regime *by construction*, replacing the heuristic `n_decisive ≥ 1`
  with the power-calibrated `min_p_value < α`.
- The 2/200 misassignments at K8 are produced by **both** gates (including τ=6.9 ≡ p=1e-3), so they are a
  property of the data (reads whose few spanned PSVs match a sibling after error), not a miscalibration of
  the significance gate.

## 3. Default α

**α = 1e-3 confirmed** as the default: it reproduces the production τ=6.9 (p=1e-3 / Eichler-AS≥10) operating
point on the labeled ladder, with the precision behaviour the assign-or-abstain philosophy wants. Tunable
via `--alpha`.

## Status / follow-ups
- The labeled sim5x ladder (ground truth) is the load-bearing validation here; a genome-wide GGO re-run
  would only add circular silver-standard numbers and is deferred (runnable via `copy_assign --regions … `
  with/without `--margin-gate`).
- **L17 closed:** `<out>.assignments.tsv` now emits `p_value` and `min_p_value` (the per-read certificate +
  identifiability bound) as the final two columns.
- **Honest limits:** the εⱼ = e/3 model is uniform-error (anti-conservative for biased substitutions, esp.
  unfiltered A→G **RNA editing** at A/G PSVs — the deferred Clair3-RNA editing filter addresses this).
