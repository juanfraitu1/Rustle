# EM copy-assignment with an identifiability-gated soft posterior — Design

**Date:** 2026-07-08  **Substrate:** gorilla (GGO) HiFi Iso-Seq; existing O2 copy-assignment machinery.

## Goal

Replace the perception that per-read copy-assignment is "wishful thinking" with a **provably-consistent
estimator plus the empirical curve that confirms it**. Build an EM (expectation–maximization) copy-assignment:
the existing per-read PSV/junction/divergence likelihoods become the E-step, copy abundances the M-step,
iterated to convergence, with a per-read soft posterior labeled by the `min_p` identifiability bound. The
headline deliverable is a **consistency theorem** (assignment → truth as coverage → ∞ in the identifiable
regime) and a **coverage sweep** on planted ground-truth that demonstrates it — i.e. "it works with enough
data," stated and shown, not asserted.

## Motivation

The shipped O2 assignment is a one-shot, per-read significance gate (assign if the PSV p-value certifies, else
abstain). It is statistically sound per read but never demonstrates that the **joint** solution is
self-consistent or that it recovers a known truth. The standard, accepted approach in multimapping / transcript
quantification (RSEM, kallisto, Salmon; the LP/facility-location cousin is Canzar's own line) is **iterative
EM**: co-estimate the latent copy origins and the copy abundances, alternating to convergence. EM supplies the
three things a skeptical reader accepts as evidence: (1) a monotone objective that converges; (2) errors and
artifacts wash out because they do not cohere with any copy's model across iterations; (3) individually-
ambiguous reads are resolved by borrowing strength through the global abundance estimate.

## Non-goals / scope

- **Copy-model refinement (stage 2) is DEFERRED.** This spec covers the abundance-EM core only: copy PSV/
  junction models `θ` are FIXED from the catalog; the M-step re-estimates only the abundances `π`. Re-estimating
  each copy's consensus/PSV profile from its assigned reads (IsoCon/SDA-like self-correction) is a documented
  follow-up, not built here. A tight, provable core is the stronger first deliverable.
- **Additive, non-destructive.** The EM path is a new `--em` mode. It does NOT replace or alter the validated
  one-shot gate; both remain runnable. Default behavior is byte-identical to today.
- No new likelihood model: the E-step reuses the existing quality-weighted PSV base likelihoods, junction
  compatibility, and divergence terms (`copy_assign.rs`, `psv_linkage.rs`). The editing / allele / spanning
  filters of the current gate cascade still decide which columns contribute to `L_rk`.

## The model (per co-located family)

A family has `K` copies. Each read `r` has a hidden copy-of-origin `z_r ∈ {1…K}`. Parameters:
- **`π`** = copy abundances (mixing weights), `Σ_k π_k = 1`.
- **`θ`** = copy models (per-copy allele at each PSV column + junction profile), FIXED from the catalog.

Per-read likelihood `L_rk = P(read r | copy k, θ_k)` = the product over the PSV columns the read spans of the
quality-weighted per-base allele-match probability, times junction compatibility, times the divergence term —
**the existing per-read likelihood already computed in the gate**. Reads spanning no distinguishing feature have
`L_rk` equal across the copies they are compatible with (this is what makes them non-identifiable).

## The EM loop

Initialize `π` from the confidently-placed reads (the current gate's certified assignments, or unique-mapping
reads); if none, uniform. Then iterate:

- **E-step (responsibilities):** `γ_rk = π_k · L_rk / Σ_j π_j · L_rj`. This is the soft posterior over copies
  for read `r`. `γ` is an `N × K` matrix; each row sums to 1.
- **M-step (abundances):** `π_k = (Σ_r γ_rk) / N`.
- **Convergence:** compute the observed-data log-likelihood `ℓ = Σ_r log(Σ_k π_k L_rk)`; stop when
  `ℓ^(t) − ℓ^(t−1) < ε` (default `ε = 1e-6·|ℓ|`) or a max-iteration cap. `ℓ` is **non-decreasing** every
  iteration (standard EM guarantee) — this monotonicity is an assertion in the tests, and the convergence
  certificate replaces any tuning threshold.

## Output: soft posterior + identifiability label

For every read emit its posterior `γ_r` and a label derived from the `min_p` identifiability bound (already in
`copy_assign.rs`):
- **`Certified`** — the read spans ≥1 distinguishing feature AND its top posterior copy clears the identifiability
  threshold (`min_p < α`, Bonferroni-corrected family-wide). A validatable hard call.
- **`SoftZone`** — the posterior stays spread over the consistent zone (the copies the read cannot be
  distinguished among). Honest; reported as a soft posterior over that zone, **never** collapsed to a hard 1/k
  call.

`<out>.em.tsv`: `read_name  family_id  argmax_copy  label  posterior(k1:p1;k2:p2;…)  n_iter`.
`<out>.em_abundance.tsv`: `family_id  copy_id  pi_hat  n_reads_soft`.

## The consistency theorem (`bench/em_consistency.md`)

**Claim.** Let a family's copies be *identifiable* iff every pair is separated by ≥1 PSV column at which the
per-read likelihood certificate is achievable (`min_p < α`). In the identifiable regime the finite mixture over
copies is identifiable, and the EM/MLE estimator is **consistent**: as per-copy coverage `n → ∞`,
`π̂ → π*` and the MAP assignment `ẑ_r → z*_r` (almost surely). Non-identifiable copies (K=0: no distinguishing
PSV) form an equivalence class on which the posterior is provably invariant to the truth and stays at the
abundance-prior — correctly surfaced as `SoftZone`, never as a confident copy.

**Basis.** Finite-mixture identifiability + MLE consistency (Redner–Walker 1984), specialized to the discrete
PSV emission model; the novelty is that the identifiability partition is *exactly* the `min_p` per-read
certificate. This is the **soft/EM relaxation of the MWCA facility-location objective** already proven in the
theory chapter (`bench/copy_assignment_theory.md`), so it slots into the existing lattice rather than standing
alone. The write-up states the assumptions, the theorem, the proof sketch, and the identifiability partition,
and points to the coverage sweep as the empirical confirmation.

## Validation harness (the demonstration)

1. **Coverage sweep on planted sim** (`bench/em_coverage_sweep.py`, extends `bench/sim_genome.py` /
   `sim_reads.py`): plant a family with known copy abundances `π*` and known per-read origins `z*`. Simulate at
   coverage `{1,2,5,10,20,50,100}×`. At each coverage run the EM and record: (a) **assignment accuracy** vs
   `z*` on identifiable reads; (b) **abundance L1 error** `‖π̂ − π*‖₁`; (c) the fraction of reads certified vs
   soft-zone. Expected/asserted: accuracy → 100% and abundance L1 → 0 as coverage grows in the identifiable
   regime; K=0 families stay at the identifiability floor (soft-zone), never wrongly forced. The resulting curve
   IS the theorem, confirmed. Emits `bench/EM_COVERAGE_SWEEP.md` + a plot-ready TSV.
2. **Real-data cross-checks** (reuse existing runs): EM assignments vs the **silver standard** (reads that also
   map uniquely — independent aligner placement) and **held-out-PSV** confirmation on a real gorilla family
   (e.g. the 8,461-read `o2_chk` family). Expected: EM certified calls agree with silver at the ~100% already
   observed for the gate.
3. **Head-to-head vs the one-shot gate:** on the same family, EM `Certified` calls match the gate's assignments
   where identifiable, and EM additionally reports soft-zone posteriors where the gate abstained. Report the
   agreement rate and the added-coverage of soft-zone localizations.

## Files

- **Create** `src/rustle/vg_family/em_copy_assign.rs`: `em_assign(reads, copies, likelihoods, alpha, eps,
  max_iter) -> EmResult` with `EmResult { posteriors: Vec<Vec<f64>>, abundances: Vec<f64>, labels: Vec<EmLabel>,
  n_iter, loglik_trace: Vec<f64> }` and `enum EmLabel { Certified, SoftZone }`. Pure functions: `e_step`,
  `m_step`, `loglik`, `label_read` (consumes the existing `min_p`). Reuses per-read `L_rk` from the current
  likelihood code — no new emission model.
- **Modify** `src/bin/copy_assign.rs`: add `--em` (+ `--em-max-iter`, `--em-eps`) to run the EM path and emit
  `<out>.em.tsv` / `<out>.em_abundance.tsv`. Off by default (existing output unchanged).
- **Modify** `src/rustle/vg_family/mod.rs`: register the module.
- **Create** `bench/em_coverage_sweep.py`, `bench/EM_COVERAGE_SWEEP.md`, `bench/em_consistency.md`.
- **Test** `src/rustle/vg_family/em_copy_assign.rs` (`#[cfg(test)]`) + `tests/em_copy_assign_integration.rs`.

## Testing (TDD)

Unit (in `em_copy_assign.rs`):
- `e_step` responsibilities sum to 1 per read and equal `π_k L_rk / Σ`.
- `m_step` abundance = mean responsibility; recovers planted `π` on a clean 2-copy toy.
- `loglik` is **non-decreasing** across iterations (monotonicity guarantee) — assert over a multi-iter run.
- `label_read`: a read spanning a decisive PSV with concentrated posterior → `Certified`; a read spanning no
  distinguishing feature → `SoftZone` with posterior ≈ abundance-proportional (K=0 stays soft).
- Convergence: EM stops when Δℓ < ε; deterministic given a fixed seed/order.

Integration (`tests/em_copy_assign_integration.rs`):
- On a small planted 3-copy family, EM abundances match `π*` within tolerance and certified accuracy = 100%.
- `--em` emits both TSVs with the documented columns; default (no `--em`) output byte-identical to today.

## Reproduce

```
cargo test --lib em_copy_assign
copy_assign --bam GGO_mm.bam --fasta GGO.fasta --regions <family> --em --out ca_em
python bench/em_coverage_sweep.py     # -> bench/EM_COVERAGE_SWEEP.md (accuracy/abundance vs coverage)
```
