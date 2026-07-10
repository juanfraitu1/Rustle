# EM copy-assignment = soft SDA PSV-clustering on the PSV-aware VG — Design

**Date:** 2026-07-08 (rev. 2 — grounded in SDA/Vollger + the PSV-aware VG per advisor steer).
**Substrate:** gorilla (GGO) HiFi Iso-Seq; existing PSV-aware VG (`psv_linkage.rs`/`layer2.rs`) + O2 gate.

## Goal

Make copy-assignment a mechanism the advisor already accepts, **derived from the data and anchored in the
paper he sent** — not a fitting procedure we picked. Concretely: express read→copy assignment as the
**maximum-likelihood soft relaxation of SDA's PSV correlation-clustering** (Vollger et al., Nat Methods 2019),
run on the thesis's existing **PSV-aware variation graph** (copies = paths, PSVs = bubbles, assignment =
threading). The EM adds copy abundances and soft responsibilities; the identifiability is the thesis's already-
proven **MCC = χ(H) / Strong-Separation** theorem (the provable-conditions layer *under* SDA's heuristic). The
headline deliverable is that this derivation is **consistent**, on two distinct axes: per-read certified
accuracy on the *identifiable* set is a `δ`/error-rate property that is already high at every coverage (NOT
coverage-driven), while abundance error `‖π̂−π*‖₁ → 0` as coverage grows (MLE consistency of `π̂`) — stated as
a theorem and shown by a coverage sweep. (The identifiable-set accuracy curve explains, but is not numerically
equal to, SDA's overall 91–93% floor, which mixes in its hard-called unidentifiable residue — see §6.)

## Prior-art grounding (why he'll believe it)

- **Vollger SDA 2019** (`reference_sda_vollger`) — *the advisor sent this as prior art (2026-06-21).* SDA builds
  a PSV graph: nodes = PSVs, reads = edges; a read carrying two PSVs on one molecule is an **attraction** edge
  (→ same copy), mutually-exclusive PSVs a **repulsion** edge (→ different copies); **correlation-clustering**
  assigns PSVs→copies and WhatsHap partitions reads. "Reads carrying PSVs on the same molecule" *is* SDA's
  attraction edge. SDA hits the K=0 floor ("virtually identical duplications need >100 kb reads").
- **Our `copy_assignment_theory.md`** — the column/allele view of the VG (columns = bubbles, allele-vectors =
  paths, MCC = χ(H) = min path-cover, Strong-Separation = bubbles distinguish paths). Memory records it as
  *"the provable-conditions theory UNDER SDA's heuristic."* This is where the consistency theorem lives.
- **PSV-aware VG** (`project_psv_aware_vg`, production Rust) — copies = paths, PSVs = bubbles, assign =
  threading; genome-wide it already resolves 90.3% of 145 families at 95.4% single-copy agreement. The EM is the
  probabilistic layer this VG was missing.
- **IsoCon 2018** (per-position real-vs-error test → which columns are true PSVs vs error), **Sudmant 2010 SUN**
  (single-position private markers), **Clair3-RNA** (A→I editing filter) — the per-column filters that decide
  which bubbles enter the graph. Already in the gate cascade; the EM reuses them unchanged.

## The object: one PSV-aware VG per family

A family is one variation graph: a shared backbone with **parallel paths = copies** and **bubbles = PSV
columns**. Copy `k` is the path carrying a specific allele at each bubble — its `CopyProfile.alleles` vector
(already built by `psv_linkage.rs` from read-supported PSV bubbles, i.e. *derived from the data*, SDA-style, not
a handed-down catalog). A read is a **partial path**: its `psv_obs` vector = the PSV alleles it carries on one
molecule (SDA's co-occurring PSVs). Threading = matching the read's carried alleles to a copy-path.

## The model (SDA correlation → likelihood)

Each read `r` has a hidden copy-of-origin `z_r ∈ {1…K}`. Parameters: copy **abundances** `π` (path usage) and
the copy **paths** `θ` (PSV-bubble allele vectors, from the VG). The per-read likelihood
`L_rk = P(read r threads copy k)` = the product over the bubbles the read spans of the quality-weighted allele-
match probability × junction compatibility × divergence — the existing `read_copy_evidence` likelihood
(Task 1). This likelihood **is** SDA's attraction/repulsion made continuous: a read that carries copy `k`'s
private allele at a bubble contributes `log(1−e)` to `k` and `log(e/3)` to the others — the soft version of an
attraction edge to `k` and repulsion from the rest. A read spanning no distinguishing bubble has `L_rk` equal
across the copies it is compatible with (SDA's un-attractable read = our K-frontier).

## The EM loop (soft correlation clustering)

Initialize `π` uniform (documented; SDA uses 15 random inits — abundance-from-certified-reads is a follow-up).
- **E-step:** `γ_rk = π_k·L_rk / Σ_j π_j·L_rj` — soft assignment of the read (partial path) to a copy-path. This
  is SDA's read-partition made soft (a fractional WhatsHap).
- **M-step:** `π_k = (Σ_r γ_rk)/N` — re-estimate path usage (copy abundance).
- **Convergence:** observed-data log-likelihood `ℓ = Σ_r log Σ_k π_k L_rk` is **non-decreasing** each iteration
  (tested invariant); stop at `Δℓ < ε`.

*(Copy-path refinement — re-estimating `θ_k` from γ-weighted reads, the direct EM analog of SDA re-clustering —
is DEFERRED to a follow-up. This spec fixes `θ` from the VG and estimates `π` + soft assignments only.)*

## Output: soft posterior + identifiability label

Per read: the posterior `γ_r` plus a label from the `min_p` bound (= the K-frontier test) — **`Certified`**
(spans a distinguishing bubble and clears `min_p < alpha/(K−1)`; a validatable hard call) vs **`SoftZone`**
(spread over the consistent zone; honest, never a hard 1/k). Files: `<out>.em.tsv`
(`read_name  family_id  argmax_copy  label  posterior  n_iter`), `<out>.em_abundance.tsv`
(`family_id  copy_id  pi_hat  n_reads_soft`).

`<out>.em_abundance.tsv` is the convergent, gate-likelihood EM (the estimator the consistency theorem below
describes): it uses `params.error_rate` and a log-likelihood convergence gate. This is a *different* estimator
from the legacy `<out>.quant.tsv` (`copy_assign_pipeline::soft_quantify_em`, fixed error rate 0.01, fixed 100
iterations, no convergence check) that this feature does not touch, so the two files' `pi_hat` for the same
family can legitimately differ — do not treat them as redundant. Prefer `.em_abundance.tsv` when citing the
theorem's estimator.

## The consistency theorem (`bench/em_consistency.md`)

**Claim.** The EM above is the ML soft relaxation of SDA's PSV correlation-clustering. In the identifiable
regime — every copy pair separated by a bubble where `min_p < alpha/(K−1)` is achievable (Strong-Separation) —
the mixture over copy-paths is identifiable and the EM/MLE is **consistent**: as per-copy coverage `n → ∞`,
`π̂ → π*` and the MAP assignment `ẑ_r → z*_r`. This is exactly the regime `copy_assignment_theory.md` proves
Strong-Separation ⟹ unique minimum path-cover (MCC = χ(H)); EM converges to that cover. Non-identifiable copies
(K=0: no distinguishing bubble) are SDA's ">100 kb" floor — a posterior-invariant class surfaced as `SoftZone`.

**Basis.** Finite-mixture identifiability + MLE consistency (Redner–Walker 1984) specialized to the discrete PSV
emission, with the identifiability partition = the `min_p` per-read certificate = SDA's attraction/repulsion
separability. The theorem *explains SDA's empirical accuracy floor* and predicts the coverage sweep.

## Validation harness (the demonstration)

1. **Coverage sweep on planted sim** (`bench/em_coverage_sweep.py`, extends `bench/sim_genome.py`): plant copies
   with known paths `θ*`, abundances `π*`, per-read origins `z*`. At coverage `{1,2,5,10,20,50,100}×`, run
   `copy_assign --em`; record assignment accuracy vs `z*` on identifiable reads, `‖π̂−π*‖₁`, certified fraction.
   Two distinct axes, not one: expected/asserted per-read accuracy on the *identifiable* set is **high at every
   coverage** (a `δ`/error-rate property, not coverage-driven — this is the curve that *explains* SDA's 91–93%
   overall floor without reproducing that number, since SDA's floor also mixes in its hard-called
   unidentifiable residue); abundance-L1 **→ 0 as coverage grows** (the genuine coverage-driven MLE-consistency
   effect); K=0 families stay soft-zone at every coverage (a boundary, not a coverage artifact). Emits
   `bench/EM_COVERAGE_SWEEP.md`.
2. **Real-data cross-checks:** EM vs **unique-mapper agreement** (independent aligner placement) + **held-out-PSV**
   on a real gorilla family; EM `Certified` calls should match the gate/VG at the ~100% already observed.
3. **Head-to-head vs the one-shot gate + the VG scan:** EM matches where identifiable, adds soft-zone posteriors
   where each abstained.

## Non-goals / scope

- Copy-path refinement (M-step on `θ`) deferred; `θ` fixed from the VG here.
- Additive & non-destructive: `--em` is a new mode; default output byte-identical; the one-shot gate and the
  VG layer are unchanged.
- No new likelihood/emission model and no new PSV caller: reuse `read_copy_evidence` (Task 1) and the existing
  IsoCon/Clair3/SUN column filters verbatim.

## Files

- **Create** `src/rustle/vg_family/em_copy_assign.rs`: `e_step`, `m_step`, `loglik`, `label_read`, `em_assign`
  (+ `EmResult`, `EmLabel`), consuming Task 1's `ReadEvidence`. No new emission model.
- **Modify** `src/bin/copy_assign.rs`: `--em` (+ `--em-max-iter`, `--em-eps`); build `Vec<ReadEvidence>` from the
  same family reads/copy-paths the VG/gate path constructs; emit the two TSVs. Off by default.
- **Modify** `src/rustle/vg_family/mod.rs`: register the module.
- **Create** `bench/em_consistency.md` (SDA-derivation + theorem), `bench/em_coverage_sweep.py`,
  `bench/EM_COVERAGE_SWEEP.md`.
- **Test** `src/rustle/vg_family/em_copy_assign.rs` (`#[cfg(test)]`) + `tests/em_copy_assign_integration.rs`.

## Testing (TDD)

Unit: `e_step` responsibilities sum to 1 and apply the abundance prior; `m_step` = mean responsibility; `loglik`
non-decreasing across iterations; `label_read` Certified/SoftZone via `min_p`; K=0 stays soft (posterior ≈ π).
Integration: planted 3-copy family → abundances match `π*` within tol, certified accuracy 100%; K=0 → all
SoftZone; `--em` emits both TSVs; default output byte-identical.

## Reproduce

```
cargo test --lib em_copy_assign
copy_assign --bam GGO_mm.bam --fasta GGO.fasta --regions <family> --em --out ca_em
python bench/em_coverage_sweep.py   # -> bench/EM_COVERAGE_SWEEP.md (SDA 91-93% -> 100% as coverage grows)
```
