# Copy-assignment significance gate (IsoCon-style) Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Replace the fixed log-LR margin τ with IsoCon's real-vs-error hypothesis test as the default ASSIGN/ABSTAIN gate in copy-assignment, producing a calibrated per-read p-value certificate and a power-aware identifiability condition.

**Architecture:** One focused change to `src/rustle/vg_family/copy_assign.rs`: a new exact Poisson-binomial helper, two new `AssignParams` knobs (`alpha`, `junction_err`) + a legacy escape hatch (`use_margin_gate`), two new `Assignment` fields (`p_value`, `min_p_value`), and the significance logic inside `assign_read`. The per-copy log-likelihood ranking and the reported `log_lr_margin` are unchanged. Downstream callers branch only on the unchanged `AssignStatus` enum, so no call sites change.

**Tech Stack:** Rust, `cargo test --lib`. Statistics: exact Poisson-binomial DP (no external crates).

## Global Constraints

- Edit only `src/rustle/vg_family/copy_assign.rs` for the engine; the calibration test lives in its `#[cfg(test)] mod tests`. Evaluation scripts live under `bench/`.
- Per-base error from the read's own quality when present (`super::copy_split::phred_err(q)`), else the flat `AssignParams::error_rate`. PSV error-to-a-specific-allele probability εⱼ = `eⱼ/3` (uniform 3-way error), clamped to `[0,1)`.
- Default `alpha = 1e-3` (matches the production τ=6.9 ≡ p=1e-3 / Eichler-AS≥10 operating point). Default `junction_err = 0.01`. Default `use_margin_gate = false` (significance is the default gate).
- Decision order is **Tied first**: `Tied` if `min_p_value ≥ alpha`, else `Assigned` if `p_value < alpha`, else `Ambiguous`.
- Empty-set conventions: `poisson_binomial_upper_tail(0, []) = 1.0`; `poisson_binomial_upper_tail(k>0, []) = 0.0`; empty product `Π ∅ = 1.0`.
- The repo is **not git-initialized**; the `Commit` steps are optional. The real checkpoint each task is `cargo test --lib` passing. If git is later initialized, run the shown commit.

---

### Task 1: Exact Poisson-binomial upper tail

**Files:**
- Modify: `src/rustle/vg_family/copy_assign.rs` (add a free `pub(crate)` fn near the top, after `p_from_tau`)
- Test: same file, `#[cfg(test)] mod tests`

**Interfaces:**
- Produces: `pub(crate) fn poisson_binomial_upper_tail(k: usize, probs: &[f64]) -> f64` — returns `P(Σ Bernoulli(probs[i]) ≥ k)`.

- [ ] **Step 1: Write the failing tests**

Add to `mod tests`:

```rust
#[test]
fn poisson_binomial_matches_binomial_for_equal_probs() {
    // n=4, p=0.5: P(X>=2) = (6+4+1)/16 = 11/16.
    let probs = [0.5_f64; 4];
    let got = poisson_binomial_upper_tail(2, &probs);
    assert!((got - 11.0 / 16.0).abs() < 1e-12, "got {got}");
}

#[test]
fn poisson_binomial_edge_cases() {
    assert_eq!(poisson_binomial_upper_tail(0, &[]), 1.0); // P(>=0) = 1
    assert_eq!(poisson_binomial_upper_tail(1, &[]), 0.0); // 0 trials cannot reach 1
    assert_eq!(poisson_binomial_upper_tail(1, &[0.0, 0.0]), 0.0); // no possible success
    assert!((poisson_binomial_upper_tail(1, &[1.0]) - 1.0).abs() < 1e-12); // certain success
    // k == n: tail equals the product of probs (all must succeed)
    let probs = [0.1_f64, 0.2, 0.05];
    let prod: f64 = probs.iter().product();
    assert!((poisson_binomial_upper_tail(3, &probs) - prod).abs() < 1e-12);
}

#[test]
fn poisson_binomial_monotone_in_k() {
    let probs = [0.3_f64, 0.4, 0.2, 0.6];
    let mut prev = poisson_binomial_upper_tail(0, &probs);
    for k in 1..=probs.len() {
        let cur = poisson_binomial_upper_tail(k, &probs);
        assert!(cur <= prev + 1e-12, "P(>={k}) must not exceed P(>={})", k - 1);
        prev = cur;
    }
}
```

- [ ] **Step 2: Run tests to verify they fail**

Run: `cargo test --lib poisson_binomial 2>&1 | tail -5`
Expected: FAIL — `cannot find function poisson_binomial_upper_tail`.

- [ ] **Step 3: Write the implementation**

Add after `p_from_tau` (around line 78):

```rust
/// Exact upper tail of a Poisson-binomial: `P(Σ Bernoulli(probs[i]) >= k)`. O(n^2) DP convolution of the
/// per-trial success probabilities (the number of trials = distinguishing positions, typically < ~100, so
/// exact is cheap and avoids any normal-approximation error). Conventions: `k == 0 -> 1.0`; `k > probs.len()
/// -> 0.0` (also `k > 0` with no trials). Probabilities are clamped to `[0, 1]`.
pub(crate) fn poisson_binomial_upper_tail(k: usize, probs: &[f64]) -> f64 {
    if k == 0 {
        return 1.0;
    }
    if k > probs.len() {
        return 0.0;
    }
    // dp[s] = P(exactly s successes) after the processed trials. Iterate s downward so dp[s-1] is still the
    // pre-trial value when read.
    let mut dp = vec![0.0f64; probs.len() + 1];
    dp[0] = 1.0;
    for &p in probs {
        let p = p.clamp(0.0, 1.0);
        for s in (0..dp.len()).rev() {
            let from_prev = if s > 0 { dp[s - 1] * p } else { 0.0 };
            dp[s] = dp[s] * (1.0 - p) + from_prev;
        }
    }
    dp[k..].iter().sum()
}
```

- [ ] **Step 4: Run tests to verify they pass**

Run: `cargo test --lib poisson_binomial 2>&1 | tail -5`
Expected: PASS — 3 tests ok.

- [ ] **Step 5: Commit (optional — repo not git-initialized)**

```bash
git add src/rustle/vg_family/copy_assign.rs && git commit -m "feat(copy_assign): exact Poisson-binomial upper tail"
```

---

### Task 2: Struct fields (AssignParams + Assignment), behavior unchanged

**Files:**
- Modify: `src/rustle/vg_family/copy_assign.rs` — `AssignParams` (around 82-110), `Assignment` (around 51-59), the `Assignment { .. }` construction in `assign_read` (around line 202)
- Test: same file, `mod tests`

**Interfaces:**
- Consumes: nothing new.
- Produces: `AssignParams { error_rate, junction_weight, boundary_tol, margin, alpha: f64, junction_err: f64, use_margin_gate: bool }`; `AssignParams::for_alpha(alpha: f64) -> Self`; `Assignment { best_copy, log_lr_margin, n_decisive, resolvable, status, p_value: f64, min_p_value: f64 }`.

This task only ADDS fields with neutral defaults and keeps the existing τ-margin gate, so all existing tests stay green. The significance computation is Task 3.

- [ ] **Step 1: Write the failing test**

```rust
#[test]
fn assign_params_alpha_defaults_and_constructor() {
    let d = AssignParams::default();
    assert_eq!(d.alpha, 1e-3);
    assert_eq!(d.junction_err, 0.01);
    assert!(!d.use_margin_gate);
    let a = AssignParams::for_alpha(0.05);
    assert_eq!(a.alpha, 0.05);
    assert_eq!(a.error_rate, d.error_rate); // other fields inherit the default
}
```

- [ ] **Step 2: Run test to verify it fails**

Run: `cargo test --lib assign_params_alpha 2>&1 | tail -5`
Expected: FAIL — no field `alpha` on `AssignParams`.

- [ ] **Step 3: Add the fields + constructor + Assignment fields**

In `AssignParams` struct, add after `pub margin: f64,`:

```rust
    /// Significance level / target per-read false-assignment rate for the IsoCon gate. Default 1e-3.
    pub alpha: f64,
    /// Per-distinguishing-junction error probability ε used in the significance test (junctions are sharp).
    pub junction_err: f64,
    /// Revert to the legacy τ-margin gate (for reproducing legacy numbers / the A/B comparison). Default false.
    pub use_margin_gate: bool,
```

In `impl Default for AssignParams`, change the constructed value to:

```rust
        AssignParams {
            error_rate: 0.003,
            junction_weight: 5.0,
            boundary_tol: 4,
            margin: 2.0,
            alpha: 1e-3,
            junction_err: 0.01,
            use_margin_gate: false,
        }
```

In `impl AssignParams`, add:

```rust
    /// Construct with the significance level set directly (the IsoCon-gate knob): `alpha = p`.
    pub fn for_alpha(alpha: f64) -> Self {
        AssignParams { alpha, ..Self::default() }
    }
```

In `Assignment` struct, add after `pub status: AssignStatus,`:

```rust
    /// IsoCon certificate: P(this assignment by error if the read were the least-distinguishable competitor).
    pub p_value: f64,
    /// Identifiability bound: the best attainable `p_value` (read supports `best` at every distinguishing
    /// position) against the hardest competitor. `>= alpha` ⇒ the read is unresolvable (`Tied`).
    pub min_p_value: f64,
```

In `assign_read`, change the final `Some(Assignment { .. })` to include the new fields with neutral
placeholders (real values come in Task 3):

```rust
    Some(Assignment {
        best_copy: copies[best].copy_id,
        log_lr_margin: margin,
        n_decisive,
        resolvable,
        status,
        p_value: 0.0,
        min_p_value: 0.0,
    })
```

- [ ] **Step 4: Run tests to verify everything passes**

Run: `cargo test --lib copy_assign 2>&1 | tail -5`
Expected: PASS — the new test + all existing `copy_assign` tests (still on the τ-margin gate).

- [ ] **Step 5: Commit (optional)**

```bash
git add src/rustle/vg_family/copy_assign.rs && git commit -m "feat(copy_assign): add alpha/junction_err params + p_value fields"
```

---

### Task 3: Significance computation + IsoCon gate in `assign_read`

**Files:**
- Modify: `src/rustle/vg_family/copy_assign.rs` — body of `assign_read` (the status/decision block around lines 180-203), plus update the existing status-asserting tests.
- Test: same file, `mod tests`

**Interfaces:**
- Consumes: `poisson_binomial_upper_tail` (Task 1); `AssignParams.alpha`, `.junction_err`, `.use_margin_gate`; `Assignment.p_value`, `.min_p_value` (Task 2); existing `boundary_present`, `super::copy_split::phred_err`.
- Produces: `assign_read` now sets `p_value`/`min_p_value` and gates on significance by default.

- [ ] **Step 1: Write the failing certificate tests**

Add to `mod tests` (helpers `cp`/`feat` build a `CopyProfile`/`ReadFeatures`; mirror the style of the existing `assign_read_assigns` test — copies with `alleles: Vec<Option<u8>>`, reads with `psv_obs`/`psv_qual`/`junctions`):

```rust
// two copies differing at columns 0 and 1; helper to build them.
fn two_copies_2psv() -> Vec<CopyProfile> {
    vec![
        CopyProfile { copy_id: 0, alleles: vec![Some(b'A'), Some(b'C')], junctions: vec![] },
        CopyProfile { copy_id: 1, alleles: vec![Some(b'G'), Some(b'T')], junctions: vec![] },
    ]
}

#[test]
fn sig_two_highq_psv_supporting_best_is_assigned() {
    // read matches copy 0 at both columns, Q40 (e=1e-4, eps≈3.3e-5) -> p≈1e-9 << 1e-3 -> Assigned.
    let r = ReadFeatures {
        psv_obs: vec![Some(b'A'), Some(b'C')],
        psv_qual: vec![Some(40), Some(40)],
        junctions: vec![],
    };
    let a = assign_read(&r, &two_copies_2psv(), &AssignParams::default()).unwrap();
    assert_eq!(a.best_copy, 0);
    assert_eq!(a.status, AssignStatus::Assigned);
    assert!(a.p_value < 1e-3, "p_value {} should clear alpha", a.p_value);
}

#[test]
fn sig_no_distinguishing_columns_is_tied() {
    // read spans neither distinguishing column -> min_p_value = 1.0 -> Tied.
    let r = ReadFeatures { psv_obs: vec![None, None], psv_qual: vec![None, None], junctions: vec![] };
    let a = assign_read(&r, &two_copies_2psv(), &AssignParams::default()).unwrap();
    assert_eq!(a.status, AssignStatus::Tied);
    assert!(a.min_p_value >= 1e-3);
}

#[test]
fn sig_conflicting_support_is_ambiguous() {
    // read matches copy 0 at col0 but copy 1 at col1 (k=1 of n=2, high-Q) -> not significant -> Ambiguous,
    // but resolvable in principle (min_p_value < alpha).
    let r = ReadFeatures {
        psv_obs: vec![Some(b'A'), Some(b'T')],
        psv_qual: vec![Some(40), Some(40)],
        junctions: vec![],
    };
    let a = assign_read(&r, &two_copies_2psv(), &AssignParams::default()).unwrap();
    assert_eq!(a.status, AssignStatus::Ambiguous);
    assert!(a.min_p_value < 1e-3, "two Q40 columns CAN resolve in principle");
    assert!(a.p_value >= 1e-3, "but this read's split evidence is not significant");
}

#[test]
fn sig_legacy_margin_gate_still_selectable() {
    // a single Q40 PSV under the legacy margin gate (margin=2.0) is Assigned, as before.
    let copies = vec![
        CopyProfile { copy_id: 0, alleles: vec![Some(b'A')], junctions: vec![] },
        CopyProfile { copy_id: 1, alleles: vec![Some(b'G')], junctions: vec![] },
    ];
    let r = ReadFeatures { psv_obs: vec![Some(b'A')], psv_qual: vec![Some(40)], junctions: vec![] };
    let p = AssignParams { use_margin_gate: true, ..AssignParams::default() };
    let a = assign_read(&r, &copies, &p).unwrap();
    assert_eq!(a.status, AssignStatus::Assigned);
}
```

- [ ] **Step 2: Run tests to verify they fail**

Run: `cargo test --lib sig_ 2>&1 | tail -8`
Expected: FAIL — statuses are computed by the old τ gate; `p_value`/`min_p_value` are still 0.0.

- [ ] **Step 3: Implement the significance certificate + gate**

Replace the decision block in `assign_read` (the part from `let margin = ...` through the final `Some(Assignment { .. })`) with:

```rust
    let margin = if n > 1 { logl[best] - second } else { f64::INFINITY };

    // --- IsoCon significance certificate: test `best` vs each competitor on their DISTINGUISHING obs ---
    // p_read = least-significant (max) pairwise p; min_p = identifiability bound (best attainable p).
    let (mut p_read, mut min_p) = (0.0f64, 0.0f64);
    for c in 0..n {
        if c == best {
            continue;
        }
        let mut eps: Vec<f64> = Vec::new();
        let mut k = 0usize;
        // distinguishing PSV columns the read spans
        for j in 0..n_cols {
            let obs = match read.psv_obs[j] {
                Some(b) => b,
                None => continue,
            };
            let ba = copies[best].alleles.get(j).copied().flatten();
            let ca = copies[c].alleles.get(j).copied().flatten();
            if let (Some(ba), Some(ca)) = (ba, ca) {
                if ba != ca {
                    let e = match read.psv_qual.get(j).copied().flatten() {
                        Some(q) => super::copy_split::phred_err(q),
                        None => p.error_rate,
                    };
                    eps.push((e / 3.0).clamp(0.0, 1.0));
                    if obs == ba {
                        k += 1;
                    }
                }
            }
        }
        // distinguishing junctions (from the read's own junction set)
        for &jb in &read.junctions {
            let in_best = boundary_present(jb, &copies[best].junctions, p.boundary_tol);
            let in_c = boundary_present(jb, &copies[c].junctions, p.boundary_tol);
            if in_best != in_c {
                eps.push(p.junction_err.clamp(0.0, 1.0));
                if in_best {
                    k += 1; // read carries a junction `best` has and `c` lacks -> supports best
                }
            }
        }
        let pbc = poisson_binomial_upper_tail(k, &eps);
        if pbc > p_read {
            p_read = pbc;
        }
        // best attainable p for this pair = Π ε  (empty distinguishing set -> 1.0 -> forces Tied)
        let attain = if eps.is_empty() { 1.0 } else { eps.iter().product::<f64>() };
        if attain > min_p {
            min_p = attain;
        }
    }
    if n <= 1 {
        // single copy: trivially resolved/assigned.
        p_read = 0.0;
        min_p = 0.0;
    }

    let (resolvable, status) = if p.use_margin_gate {
        let resolvable = n_decisive >= 1;
        let status = if !resolvable {
            AssignStatus::Tied
        } else if margin >= p.margin {
            AssignStatus::Assigned
        } else {
            AssignStatus::Ambiguous
        };
        (resolvable, status)
    } else {
        let resolvable = min_p < p.alpha;
        let status = if !resolvable {
            AssignStatus::Tied
        } else if p_read < p.alpha {
            AssignStatus::Assigned
        } else {
            AssignStatus::Ambiguous
        };
        (resolvable, status)
    };

    Some(Assignment {
        best_copy: copies[best].copy_id,
        log_lr_margin: margin,
        n_decisive,
        resolvable,
        status,
        p_value: p_read,
        min_p_value: min_p,
    })
```

- [ ] **Step 4: Run the new certificate tests**

Run: `cargo test --lib sig_ 2>&1 | tail -8`
Expected: PASS — the 4 `sig_` tests.

- [ ] **Step 5: Fix the existing status-asserting tests broken by the default gate change**

Run the full module to find the fallout: `cargo test --lib copy_assign 2>&1 | tail -20`.
For each FAILING pre-existing test that asserts `status == Assigned` / `Tied` under default params (e.g. the
single-PSV `assign_read_assigns`-style test): make its intent explicit by EITHER (a) pinning the legacy gate
`AssignParams { use_margin_gate: true, ..Default::default() }`, OR (b) giving it two high-Q PSVs so it clears
α=1e-3. Tests that only assert on `log_lr_margin` (the quality-weighting tests) are unaffected — `margin` is
still computed. Do not weaken a test to pass; preserve its original intent under the appropriate gate.

- [ ] **Step 6: Run the whole suite green**

Run: `cargo test --lib 2>&1 | tail -3`
Expected: PASS — `test result: ok.` (no failures).

- [ ] **Step 7: Commit (optional)**

```bash
git add src/rustle/vg_family/copy_assign.rs && git commit -m "feat(copy_assign): IsoCon significance gate as default decision"
```

---

### Task 4: Calibration integration test (the headline validation)

**Files:**
- Test: `src/rustle/vg_family/copy_assign.rs` `mod tests`

**Interfaces:**
- Consumes: `assign_read`, `AssignParams::for_alpha`, `AssignStatus`.

This proves the certificate is *calibrated*: among reads the gate ASSIGNS, the realized misassignment rate stays at/below α — the property the τ-margin only approximated.

- [ ] **Step 1: Write the calibration test**

```rust
#[test]
fn sig_gate_is_calibrated_realized_error_tracks_alpha() {
    // Two copies differing at 6 PSV columns. Simulate reads from a known true copy, inject errors at a
    // realistic HiFi rate, run the gate, and check the realized misassignment among ASSIGNED reads.
    let m = 6usize;
    let copy0: Vec<Option<u8>> = (0..m).map(|_| Some(b'A')).collect();
    let copy1: Vec<Option<u8>> = (0..m).map(|_| Some(b'C')).collect();
    let copies = vec![
        CopyProfile { copy_id: 0, alleles: copy0.clone(), junctions: vec![] },
        CopyProfile { copy_id: 1, alleles: copy1.clone(), junctions: vec![] },
    ];
    // deterministic LCG
    let mut state = 0x2545F4914F6CDD1Du64;
    let mut next = || {
        state ^= state << 13;
        state ^= state >> 7;
        state ^= state << 17;
        state
    };
    let bases = [b'A', b'C', b'G', b'T'];
    let q = 30u8; // e = 1e-3
    let e = 10f64.powf(-(q as f64) / 10.0); // inline (phred_err's path differs inside mod tests)
    let run = |alpha: f64, next: &mut dyn FnMut() -> u64| -> (usize, usize) {
        let p = AssignParams::for_alpha(alpha);
        let (mut assigned, mut wrong) = (0usize, 0usize);
        for _ in 0..20_000 {
            let truth = (next() % 2) as usize;
            let template = if truth == 0 { &copy0 } else { &copy1 };
            // each column: emit the true base, or with prob e a random base
            let psv_obs: Vec<Option<u8>> = template
                .iter()
                .map(|t| {
                    let tb = t.unwrap();
                    if (next() as f64 / u64::MAX as f64) < e {
                        Some(bases[(next() % 4) as usize])
                    } else {
                        Some(tb)
                    }
                })
                .collect();
            let r = ReadFeatures { psv_obs, psv_qual: vec![Some(q); m], junctions: vec![] };
            let a = assign_read(&r, &copies, &p).unwrap();
            if a.status == AssignStatus::Assigned {
                assigned += 1;
                if a.best_copy != truth {
                    wrong += 1;
                }
            }
        }
        (assigned, wrong)
    };
    let (a_hi, w_hi) = run(1e-2, &mut next);
    let (a_lo, w_lo) = run(1e-4, &mut next);
    assert!(a_hi > 1000 && a_lo > 1000, "should assign a substantial fraction ({a_hi}, {a_lo})");
    // realized misassignment among assigned must respect the bound (conservative p-value), with finite-sample slack
    let rate_hi = w_hi as f64 / a_hi as f64;
    let rate_lo = w_lo as f64 / a_lo as f64;
    assert!(rate_hi <= 3e-2, "alpha=1e-2 realized {rate_hi}");
    assert!(rate_lo <= 3e-4, "alpha=1e-4 realized {rate_lo}");
    // calibration is monotone: a stricter alpha yields a lower (or equal) realized error
    assert!(rate_lo <= rate_hi + 1e-9, "stricter alpha must not increase error ({rate_lo} vs {rate_hi})");
}
```

- [ ] **Step 2: Run it to verify it fails (until Task 3 is in) / passes**

Run: `cargo test --lib sig_gate_is_calibrated -- --nocapture 2>&1 | tail -8`
Expected: PASS (Task 3 implemented). If it FAILS on the rate bounds, the certificate math is wrong — debug `assign_read`, do not loosen the bounds.

- [ ] **Step 3: Commit (optional)**

```bash
git add src/rustle/vg_family/copy_assign.rs && git commit -m "test(copy_assign): calibration test for the significance gate"
```

---

### Task 5: Evaluation on sim5x + real GGO (compare gates, pick the operating point)

**Files:**
- Create: `bench/eval_significance_gate.md` (a short results note)
- Uses: the existing `copy_assign` binary + sim5x harness (`bench/sim_reads.py`) + the GGO conflict-catalog substrate.

This task is an execution/measurement step (no new unit test). It produces the numbers that confirm the gate in practice and choose the default α.

- [ ] **Step 1: Expose the gate knobs on the binary, then build**

Inspect `src/bin/copy_assign.rs` for where it constructs `AssignParams`. If it does not already expose them,
add two clap flags wired into that construction (≈5 lines): `--alpha <f64>` → `AssignParams.alpha`, and
`--margin-gate` (bool) → `AssignParams.use_margin_gate` (for the legacy A/B). Default α unchanged (1e-3).
Then: `cargo build --release --bin copy_assign 2>&1 | tail -2` → Expected `Finished`.

- [ ] **Step 2: sim5x labeled-truth comparison**

Locate the sim5x runner (`bench/sim_reads.py` generates the K=0..8 PSV ladder; the `copy_assign` CLI assigns).
Run the assignment on the sim5x labeled reads under: (a) the new significance gate at α ∈ {1e-2, 1e-3, 1e-4},
and (b) the legacy gate (`RUSTLE_*`/CLI flag for `use_margin_gate`, or `for_target_misassignment`). Record per
α: recall (% Assigned), accuracy of Assigned reads vs labels, and the realized-vs-α calibration. Verify K=0
(exonically-identical) reads are `Tied` via `min_p_value ≥ α`. Write the table into `bench/eval_significance_gate.md`.

- [ ] **Step 3: Real GGO comparison**

Re-run copy-assignment on the GGO conflict-catalog substrate (the per-region driver in
`/home/juanfra/winloci_scratch/ca_gw/`) with the significance gate; record assigned/ambiguous/tied % and
unique-mapper agreement, vs the legacy τ-gate baseline. Append to `bench/eval_significance_gate.md`.

- [ ] **Step 4: Decide the default α**

If the sim5x calibration confirms α=1e-3 gives the intended precision at acceptable recall, keep the default.
Otherwise adjust `AssignParams::default().alpha` and note the rationale in the eval doc. Re-run `cargo test --lib`.

- [ ] **Step 5: Commit (optional)**

```bash
git add bench/eval_significance_gate.md src/rustle/vg_family/copy_assign.rs && git commit -m "eval: significance gate vs tau-margin on sim5x + GGO"
```

---

## Notes for the implementer

- `phred_err(q)` converts a Phred byte to an error probability (`10^(-q/10)`); it already exists in `super::copy_split`.
- The per-copy log-likelihood loop and `n_decisive` counting are UNCHANGED — only the decision block at the end of `assign_read` is replaced. Keep them.
- The significance loop reuses the same per-column quality lookup as the likelihood loop; keep the εⱼ = `e/3` convention consistent between them.
- Do not add external crates. The Poisson-binomial DP is the only statistics needed.
