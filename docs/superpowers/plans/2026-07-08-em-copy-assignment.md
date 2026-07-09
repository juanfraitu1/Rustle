# EM copy-assignment — Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** An EM copy-assignment whose E-step is the existing per-read PSV/junction likelihood and whose M-step is copy abundance, iterated to convergence, emitting a per-read soft posterior labeled by the `min_p` identifiability bound — plus a consistency theorem and a coverage sweep that demonstrates convergence to truth as coverage grows.

**Architecture:** New `em_copy_assign.rs` module of pure EM functions that consume per-read per-copy log-likelihoods extracted (DRY) from `copy_assign::assign_read_editing`. A `--em` mode in the `copy_assign` binary builds the evidence from a family's reads/copies (reusing the existing integration), runs the EM, and emits two TSVs. Validation is a planted-ground-truth coverage sweep + a theorem write-up. Additive: default output unchanged.

**Tech Stack:** Rust (existing `vg_family` crate), Python (bench validation), reusing `copy_assign.rs`/`psv_linkage.rs`/`sim_genome.py`.

## Global Constraints

- **Additive & non-destructive.** The `--em` path must not change any existing output. With `--em` absent, the binary's `.assignments.tsv`/`.families.tsv` are byte-identical to today. The one-shot gate (`assign_read_editing`) keeps its exact current behavior.
- **DRY likelihood.** The EM must consume the SAME per-copy log-likelihood the gate computes — extracted into one shared function, not re-derived. No second emission model.
- **No arbitrary thresholds.** Convergence is `Δℓ < eps` on the monotone log-likelihood; the identifiability label uses the existing Bonferroni `alpha/(K-1)` `min_p` rule. No new tuning constants beyond `eps` (numeric tolerance) and the reused `alpha`.
- **Monotonicity is a tested invariant.** The observed-data log-likelihood is non-decreasing every iteration; a test asserts it.
- Reuse `AssignParams` (`error_rate`, `alpha`, `boundary_tol`, `junction_weight`) verbatim; do not introduce a parallel params struct.

---

## File Structure

- `src/rustle/vg_family/copy_assign.rs` — MODIFY: extract `read_copy_evidence` (Task 1).
- `src/rustle/vg_family/em_copy_assign.rs` — CREATE: EM core + driver (Tasks 2–3).
- `src/rustle/vg_family/mod.rs` — MODIFY: register module (Task 2).
- `src/bin/copy_assign.rs` — MODIFY: `--em` wiring + emit (Task 4).
- `bench/em_consistency.md` — CREATE: theorem (Task 5).
- `bench/em_coverage_sweep.py`, `bench/EM_COVERAGE_SWEEP.md` — CREATE: coverage sweep (Task 6).
- `tests/em_copy_assign_integration.rs` — CREATE: end-to-end (Task 7).

---

### Task 1: Extract per-read per-copy evidence (DRY the likelihood)

**Files:**
- Modify: `src/rustle/vg_family/copy_assign.rs`
- Test: `src/rustle/vg_family/copy_assign.rs` (`#[cfg(test)]`)

**Interfaces:**
- Produces: `pub(crate) struct ReadEvidence { pub logl: Vec<f64>, pub min_p: f64, pub n_decisive: usize }` and `pub(crate) fn read_copy_evidence(read: &ReadFeatures, copies: &[CopyProfile], p: &AssignParams, editing_cols: &[bool]) -> ReadEvidence`. `logl[ci]` = the per-copy log-likelihood currently accumulated in `assign_read_editing`'s `logl` vector (PSV term + junction term); `min_p` = the identifiability bound currently computed there; `n_decisive` = the decisive-feature count.

- [ ] **Step 1: Write the failing test** — pin that the extracted evidence equals the gate's internals.

```rust
#[test]
fn read_copy_evidence_matches_assignment_internals() {
    // two copies differ at column 0 (A vs C); read observes A -> favors copy 0.
    let copies = vec![
        CopyProfile { copy_id: 0, alleles: vec![Some(b'A')], junctions: vec![] },
        CopyProfile { copy_id: 1, alleles: vec![Some(b'C')], junctions: vec![] },
    ];
    let read = ReadFeatures { psv_obs: vec![Some(b'A')], psv_qual: vec![], junctions: vec![] };
    let p = AssignParams::for_alpha(1e-3);
    let ev = read_copy_evidence(&read, &copies, &p, &[false]);
    // logl favors copy 0; min_p equals the Assignment's identifiability bound; one decisive column.
    assert!(ev.logl[0] > ev.logl[1]);
    assert_eq!(ev.n_decisive, 1);
    let a = assign_read_editing(&read, &copies, &p, &[false]).unwrap();
    assert!((ev.min_p - a.min_p_value).abs() < 1e-12);
    // posterior consistency: softmax(logl) == Assignment.posterior
    let m = ev.logl.iter().cloned().fold(f64::NEG_INFINITY, f64::max);
    let mut e: Vec<f64> = ev.logl.iter().map(|&l| (l - m).exp()).collect();
    let z: f64 = e.iter().sum(); for x in &mut e { *x /= z; }
    for (i, &pv) in a.posterior.iter().enumerate() { assert!((e[i] - pv).abs() < 1e-9); }
}
```

- [ ] **Step 2: Run it, verify it fails** — `cargo test --lib read_copy_evidence_matches` → FAIL (function not defined).
- [ ] **Step 3: Extract the function.** Move the `logl` accumulation (the PSV term loop and the junction term loop, current `copy_assign.rs:265–322`) and the `min_p`/`n_decisive` computation into `read_copy_evidence`, returning `ReadEvidence`. Refactor `assign_read_editing` to call `read_copy_evidence` and then do its argmax/margin/gate/posterior on `ev.logl`/`ev.min_p`/`ev.n_decisive` — behavior unchanged. (Note: `min_p` in the gate depends on `best`; keep `min_p` computation as-is by returning it from the same pass — `read_copy_evidence` computes `logl`+`n_decisive`, then argmax `best`, then the `copy_pair_significance` loop for `min_p`, exactly as today, and returns all three. `assign_read_editing` reuses `best` via recomputing argmax on `ev.logl`, which is O(K) and identical.)
- [ ] **Step 4: Run tests** — `cargo test --lib copy_assign` → all PASS (new test + every pre-existing `copy_assign` test unchanged).
- [ ] **Step 5: Commit** — `git add src/rustle/vg_family/copy_assign.rs && git commit -m "refactor(copy_assign): extract read_copy_evidence (per-copy logl + min_p) for EM reuse"`

---

### Task 2: EM core pure functions

**Files:**
- Create: `src/rustle/vg_family/em_copy_assign.rs`
- Modify: `src/rustle/vg_family/mod.rs` (add `pub mod em_copy_assign;`)
- Test: `src/rustle/vg_family/em_copy_assign.rs` (`#[cfg(test)]`)

**Interfaces:**
- Consumes: `Vec<Vec<f64>>` per-read per-copy log-likelihoods (from Task 1's `ReadEvidence.logl`).
- Produces: `pub fn e_step(logl: &[Vec<f64>], pi: &[f64]) -> Vec<Vec<f64>>`, `pub fn m_step(gamma: &[Vec<f64>]) -> Vec<f64>`, `pub fn loglik(logl: &[Vec<f64>], pi: &[f64]) -> f64`.

**Math:** `γ_rk = softmax_k(logl_rk + ln π_k)`; `π_k = (Σ_r γ_rk)/N`; `ℓ = Σ_r logsumexp_k(logl_rk + ln π_k)`. Use log-sum-exp for numeric stability. Treat `π_k = 0` as `ln π_k = -inf` (copy drops out).

- [ ] **Step 1: Write failing tests.**

```rust
use super::em_copy_assign::{e_step, m_step, loglik};

#[test] fn e_step_rows_are_normalized_posteriors() {
    let logl = vec![vec![0.0_f64, -10.0], vec![-10.0, 0.0]];
    let pi = vec![0.5, 0.5];
    let g = e_step(&logl, &pi);
    for row in &g { assert!((row.iter().sum::<f64>() - 1.0).abs() < 1e-12); }
    assert!(g[0][0] > 0.99 && g[1][1] > 0.99);           // each read favors its copy
}
#[test] fn e_step_applies_the_abundance_prior() {
    // read equidistant in likelihood; prior 0.9/0.1 must tilt the posterior to 0.9/0.1.
    let g = e_step(&vec![vec![0.0, 0.0]], &vec![0.9, 0.1]);
    assert!((g[0][0] - 0.9).abs() < 1e-9 && (g[0][1] - 0.1).abs() < 1e-9);
}
#[test] fn m_step_is_mean_responsibility() {
    let g = vec![vec![1.0, 0.0], vec![0.0, 1.0], vec![1.0, 0.0]];
    let pi = m_step(&g);
    assert!((pi[0] - 2.0/3.0).abs() < 1e-12 && (pi[1] - 1.0/3.0).abs() < 1e-12);
}
#[test] fn loglik_is_nondecreasing_under_em() {
    // 3 reads, 2 copies, well-separated; one full E/M sweep must not decrease ℓ.
    let logl = vec![vec![0.0, -8.0], vec![0.0, -8.0], vec![-8.0, 0.0]];
    let mut pi = vec![0.5, 0.5];
    let l0 = loglik(&logl, &pi);
    let g = e_step(&logl, &pi); pi = m_step(&g);
    let l1 = loglik(&logl, &pi);
    assert!(l1 >= l0 - 1e-12, "EM decreased loglik: {l0} -> {l1}");
}
```

- [ ] **Step 2: Run, verify fail** — `cargo test --lib em_copy_assign` → FAIL (module missing).
- [ ] **Step 3: Implement** the three functions with log-sum-exp (register the module in `mod.rs`).
- [ ] **Step 4: Run** — `cargo test --lib em_copy_assign` → PASS.
- [ ] **Step 5: Commit** — `git commit -m "feat(em_copy_assign): E-step/M-step/loglik pure functions"`

---

### Task 3: EM driver + identifiability label

**Files:**
- Modify: `src/rustle/vg_family/em_copy_assign.rs`
- Test: same file (`#[cfg(test)]`)

**Interfaces:**
- Consumes: `&[ReadEvidence]` (Task 1), K, `alpha`, `eps`, `max_iter`.
- Produces: `pub enum EmLabel { Certified, SoftZone }`; `pub struct EmResult { pub posteriors: Vec<Vec<f64>>, pub abundances: Vec<f64>, pub labels: Vec<EmLabel>, pub n_iter: usize, pub loglik_trace: Vec<f64> }`; `pub fn em_assign(evidence: &[super::copy_assign::ReadEvidence], k: usize, alpha: f64, eps: f64, max_iter: usize) -> EmResult`; `pub fn label_read(min_p: f64, alpha: f64, k: usize) -> EmLabel` (`Certified` iff `min_p < alpha/((k-1).max(1))`, else `SoftZone`).

**Driver:** init `pi` uniform (documented; abundance-from-certified-reads is a follow-up); loop E→M, push `loglik` to `loglik_trace`, stop when `Δℓ < eps*|ℓ|` or `max_iter`; final `posteriors` = last E-step; `labels[r] = label_read(evidence[r].min_p, alpha, k)`.

- [ ] **Step 1: Write failing tests.**

```rust
#[test] fn em_recovers_planted_abundance_3copy() {
    // 3 copies, one decisive column each (one-hot logl), planted 5:3:2 reads -> pi ~ 0.5,0.3,0.2
    let ev = |best: usize| super::copy_assign::ReadEvidence {
        logl: (0..3).map(|c| if c==best {0.0} else {-12.0}).collect(), min_p: 1e-6, n_decisive: 1 };
    let mut reads = vec![]; for _ in 0..5 {reads.push(ev(0));} for _ in 0..3 {reads.push(ev(1));} for _ in 0..2 {reads.push(ev(2));}
    let r = em_assign(&reads, 3, 1e-3, 1e-9, 200);
    assert!((r.abundances[0]-0.5).abs()<0.02 && (r.abundances[1]-0.3).abs()<0.02 && (r.abundances[2]-0.2).abs()<0.02);
    for w in r.loglik_trace.windows(2) { assert!(w[1] >= w[0]-1e-9); }   // monotone
    assert!(r.labels.iter().all(|l| matches!(l, EmLabel::Certified)));   // all identifiable
}
#[test] fn em_k0_reads_stay_soft_and_abundance_proportional() {
    // no distinguishing feature: logl equal across copies, min_p >= alpha -> SoftZone, posterior == pi.
    let flat = super::copy_assign::ReadEvidence { logl: vec![0.0, 0.0], min_p: 1.0, n_decisive: 0 };
    let reads = vec![flat.clone(), flat.clone(), flat];
    let r = em_assign(&reads, 2, 1e-3, 1e-9, 50);
    assert!(r.labels.iter().all(|l| matches!(l, EmLabel::SoftZone)));
    for row in &r.posteriors { assert!((row[0]-r.abundances[0]).abs() < 1e-6); } // soft == prior, never a hard 1/k call
}
#[test] fn label_read_uses_bonferroni_min_p() {
    assert!(matches!(label_read(1e-6, 1e-3, 3), EmLabel::Certified));  // 1e-6 < 1e-3/2
    assert!(matches!(label_read(0.4,  1e-3, 3), EmLabel::SoftZone));
}
```

- [ ] **Step 2: Run, verify fail.** `cargo test --lib em_copy_assign` → FAIL.
- [ ] **Step 3: Implement** `EmLabel`, `EmResult`, `label_read`, `em_assign`.
- [ ] **Step 4: Run** — PASS.
- [ ] **Step 5: Commit** — `git commit -m "feat(em_copy_assign): em_assign driver + identifiability label (Certified/SoftZone)"`

---

### Task 4: `--em` binary wiring + emit

**Files:**
- Modify: `src/bin/copy_assign.rs`
- Test: `src/bin/copy_assign.rs` (`#[cfg(test)]`) or `tests/em_copy_assign_integration.rs` (Task 7 covers end-to-end; here add a focused gating test)

**Interfaces:** add CLI `--em` (bool), `--em-max-iter` (default 500), `--em-eps` (default 1e-6). When `--em`: for each family, build `Vec<ReadEvidence>` by calling `read_copy_evidence` on the same `ReadFeatures`/`CopyProfile` the gate path already constructs; run `em_assign`; write `<out>.em.tsv` (`read_name  family_id  argmax_copy  label  posterior  n_iter`) and `<out>.em_abundance.tsv` (`family_id  copy_id  pi_hat  n_reads_soft`). `posterior` field = `k:p` pairs joined by `;`. `n_reads_soft` = `Σ_r γ_rk`.

- [ ] **Step 1: Write failing test** — assert `--em` produces the two files with correct headers, and that WITHOUT `--em` neither file is written and the existing outputs are unchanged. (Fixture: reuse the smallest existing `copy_assign` test BAM/region, or a synthetic 2-copy `ReadFeatures` set invoked through the family loop.)
- [ ] **Step 2: Run, verify fail.**
- [ ] **Step 3: Implement** the flag + emit, gated so the default path is untouched.
- [ ] **Step 4: Run** — the new test + the existing binary tests PASS; confirm default `.assignments.tsv` byte-identical (diff against a pre-change run of the fixture).
- [ ] **Step 5: Commit** — `git commit -m "feat(copy_assign): --em mode emits soft-posterior + abundance TSVs (default unchanged)"`

---

### Task 5: Consistency theorem write-up

**Files:** Create `bench/em_consistency.md`.

Doc deliverable (no code). Contents: (1) the mixture model + assumptions; (2) **Theorem** — in the identifiable regime (every copy pair separated by a PSV where `min_p < alpha/(K-1)` is achievable) the mixture is identifiable and the EM/MLE is consistent (`π̂→π*`, MAP `ẑ_r→z*_r` as coverage `n→∞`); (3) proof sketch citing finite-mixture identifiability + MLE consistency (Redner–Walker 1984) specialized to the discrete PSV emission; (4) the **identifiability partition = the `min_p` per-read certificate**, K=0 copies form a posterior-invariant class surfaced as `SoftZone`; (5) the tie to the MWCA facility-location objective in `bench/copy_assignment_theory.md` (EM = its soft/LP relaxation); (6) pointer to the coverage sweep (Task 6) as empirical confirmation.

- [ ] **Step 1:** Write `bench/em_consistency.md` per the outline above.
- [ ] **Step 2: Commit** — `git commit -m "docs(em): consistency theorem (MLE consistency in the identifiable=min_p regime)"`

---

### Task 6: Coverage sweep (the demonstration)

**Files:** Create `bench/em_coverage_sweep.py`, `bench/EM_COVERAGE_SWEEP.md`.

**Interfaces:** extend `bench/sim_genome.py`/`sim_reads.py` to plant a family with known abundances `π*` and per-read origins `z*`. For coverage in `{1,2,5,10,20,50,100}`: simulate reads, run `copy_assign --em`, parse `<out>.em.tsv`/`.em_abundance.tsv`, compute (a) certified-read assignment accuracy vs `z*`, (b) `‖π̂−π*‖₁`, (c) certified fraction. Write a plot-ready TSV + `EM_COVERAGE_SWEEP.md` with the table and the interpretation (accuracy→100%, abundance-L1→0 in the identifiable regime; K=0 planted family stays soft-zone).

- [ ] **Step 1:** Implement `em_coverage_sweep.py` (plant π*/z*, sweep, run binary, score).
- [ ] **Step 2: Run** it end-to-end; confirm accuracy rises and abundance-L1 falls monotonically-ish with coverage in the identifiable family, and the K=0 family stays soft-zone.
- [ ] **Step 3:** Write `bench/EM_COVERAGE_SWEEP.md` with the results table + interpretation.
- [ ] **Step 4: Commit** — `git commit -m "bench(em): coverage sweep — accuracy->100%, abundance-L1->0 as coverage grows"`

---

### Task 7: End-to-end integration test

**Files:** Create `tests/em_copy_assign_integration.rs`.

- [ ] **Step 1: Write the test** — planted 3-copy family (in-memory `ReadFeatures`/`CopyProfile`, or a tiny fixture): `em_assign` abundances within 2% of `π*`, all identifiable reads `Certified` with 100% argmax accuracy vs truth, `loglik_trace` monotone. Plus: a K=0 family → all `SoftZone`.
- [ ] **Step 2: Run, verify it exercises the real API.**
- [ ] **Step 3: Run** — `cargo test --test em_copy_assign_integration` → PASS; then `cargo test --lib` full suite green.
- [ ] **Step 4: Commit** — `git commit -m "test(em): end-to-end integration — planted family recovery + K=0 soft-zone"`

---

## Self-Review

- **Spec coverage:** model+E/M (Tasks 2–3), output soft-posterior+label (Task 3–4), theorem (Task 5), coverage sweep (Task 6), additive `--em` (Task 4), DRY likelihood (Task 1) — all covered.
- **Type consistency:** `ReadEvidence{logl,min_p,n_decisive}` produced in Task 1, consumed in Tasks 3/4; `e_step/m_step/loglik` signatures consistent Task 2→3; `EmResult`/`EmLabel` Task 3→4/7.
- **Placeholder scan:** none — every code step carries real code or a concrete content outline (Tasks 5/6 are doc/validation deliverables with explicit contents).
