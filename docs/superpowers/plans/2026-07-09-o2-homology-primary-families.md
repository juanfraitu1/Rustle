# O1↔O2 harmony: E_r family membership in `copy_assign` — Implementation Plan

> REQUIRED SUB-SKILL: superpowers:subagent-driven-development.

**Goal:** Opt-in `copy_assign --homology-primary` in which family **membership** is the E_r homology component
(not the E_c conflict graph), so uniquely-mappable copies stop being dropped and returning spurious `tied`.
Measure the delta; no default flip.

**Architecture:** Swap the **edge oracle** only. `reps → {conflict_edges | homology_edges_all_reps + gamma_quasi_clique_partition} → SplitFamily → colocated_families → copies → assign`. Everything downstream (PSV discovery, assignment, EM, `chi_H`, `famcn_readonly`) is untouched. Conflict/PSV/χ(H) remain *within* the family.

## Global Constraints
- `--homology-primary` OFF ⟹ the E_c path is literally untouched ⟹ byte-identical (except the separately-declared `--min-copies` default change).
- **Pure E_r** when on — conflict edges are NOT unioned into membership.
- **No silent degradation:** `homology_edges_all_reps` returns `Result` (shells minimap2). On failure the opt-in path must **abort loudly** with a message naming `RUSTLE_MINIMAP2` — never fall back to E_c.
- γ = **0.20** (the value O1's `gamma_quasi_clique_partition` tests use). Pass it, don't hard-code a new constant elsewhere.
- Reuse: `homology_edges_all_reps` (denovo_pipeline.rs:1739, `pub(crate)`), `homology_refine_params` (:1446), `gamma_quasi_clique_partition` (family_split.rs:477), `conflict_to_split_families` (denovo_pipeline.rs).

## File Structure
- `src/rustle/vg_family/denovo_pipeline.rs` — MODIFY (T1 core extraction, T2 oracle branch).
- `src/bin/copy_assign.rs` — MODIFY (T3 flag + `min_copies` default).
- `bench/homology_primary_delta.py`, `bench/HOMOLOGY_PRIMARY_DELTA.md` — CREATE (T4).

---

### Task 1: Generalize the components→SplitFamily conversion (DRY, no behavior change)

**Files:** Modify `src/rustle/vg_family/denovo_pipeline.rs`. Test: in-crate.

Today `conflict_to_split_families(families: &[Vec<usize>], c_edges: &[(usize,usize,usize)], p: &SplitParams) -> Vec<SplitFamily>` converts `usize` edge weights to `f64` internally, then runs `community_stats` + `classify` + a deterministic sort.

**Extract the core:** `fn to_split_families(families: &[Vec<usize>], edges: &[(usize,usize,f64)], p: &SplitParams) -> Vec<SplitFamily>` containing the existing body from the `float_edges` line onward. Re-implement `conflict_to_split_families` as: convert `c_edges` to `f64` then delegate. The homology path (T2) calls `to_split_families` directly with weight `1.0` edges.

- [ ] **Test `to_split_families_matches_conflict_wrapper`:** build 4 reps, 3 conflict edges with weights `(0,1,5),(1,2,3),(3,?)`; assert `conflict_to_split_families(fams, &c_edges, &p)` equals `to_split_families(fams, &c_edges_as_f64, &p)` field-for-field (members, class). Write RED (fn missing), implement, GREEN.
- [ ] Confirm every pre-existing `denovo_pipeline` test still passes (behavior-preserving refactor).
- [ ] Commit: `refactor(denovo_pipeline): extract to_split_families(f64 edges) core for the homology path`

### Task 2: Branch the edge oracle in `detect_and_assign`

**Files:** Modify `src/rustle/vg_family/denovo_pipeline.rs`. Test: in-crate.

**Add** `pub homology_primary: bool` to `DenovoConfig` (default `false` in its `Default` impl).

In `detect_and_assign` (~line 881, "Conflict-graph: AUTHORITATIVE family criterion"): keep `build_read_placements` (its `placements` are used downstream, e.g. mapq0 support), then branch how `(families, edges_f64)` are produced:

- `cfg.homology_primary == false` → exactly as today: `conflict_edges` → `conflict_families`, plus the existing `de⊆AS` audit `eprintln!`. Weight-map to `f64` for `to_split_families`.
- `cfg.homology_primary == true` →
  ```
  let refine = homology_refine_params(None, <threads from cfg or the existing threads value>);
  let e2 = homology_edges_all_reps(&reps, &refine)
      .expect("--homology-primary: homology_edges_all_reps failed — is minimap2 on PATH or RUSTLE_MINIMAP2 set?");
  let edges_f64: Vec<(usize,usize,f64)> = e2.iter().map(|&(a,b)| (a,b,1.0)).collect();
  let families = gamma_quasi_clique_partition(reps.len(), &edges_f64, 0.20);
  eprintln!("[detect_and_assign] homology (E_r): {} edges -> {} families", e2.len(), families.len());
  ```
  (The `.expect` is the loud abort; do NOT fall back to E_c.)

Then, both branches: `let split_fams = to_split_families(&families, &edges_f64, &split_params);` and continue into `colocated_families` unchanged.

- [ ] **Test `homology_primary_admits_uniquely_mappable_rep`:** construct reps where rep0 and rep1 are homologous (their exon-sum sequences are near-identical) but generate BAM reads such that **no conflict edge** forms between them (reads map uniquely). Assert: with `homology_primary=false` the two reps land in **separate** families (or rep1 is absent); with `homology_primary=true` they are in **one** family. If building real reads is heavy, factor the branch so the *oracle* is testable directly: assert `gamma_quasi_clique_partition` over `homology_edges_all_reps`' output unions the pair while `conflict_families` over an empty edge set does not.
- [ ] **Test `homology_primary_off_is_byte_identical`:** with `homology_primary=false`, the produced `SplitFamily` vector is identical to the pre-change code path on a fixture.
- [ ] Commit: `feat(denovo_pipeline): --homology-primary edge oracle (E_r membership) in detect_and_assign`

### Task 3: CLI — `--homology-primary` and `--min-copies` default 3→2

**Files:** Modify `src/bin/copy_assign.rs`. Test: in-crate/CLI.

- Add `#[arg(long)] homology_primary: bool` → sets `cfg.homology_primary`.
  Doc-comment must state: membership becomes E_r homology; conflict/PSV/χ(H) stay within-family; **assignments will shift** (larger copy set ⟹ stricter Bonferroni `α/(K−1)`); requires minimap2.
- Change `min_copies` default from `3` to **`2`**. Doc-comment: this is an **independent, deliberate** change (2-copy homologous families are the majority and were invisible to O2); it alters default family detection on its own, separate from `--homology-primary`.

- [ ] **Test:** `--homology-primary` off ⟹ no behavior change beyond `min_copies`; the flag is wired to `cfg.homology_primary`. Confirm `cargo build --release --bin copy_assign` compiles.
- [ ] Commit: `feat(copy_assign): --homology-primary flag; --min-copies default 3->2 (declared)`

### Task 4: Measure the delta (the deliverable)

**Files:** Create `bench/homology_primary_delta.py`, `bench/HOMOLOGY_PRIMARY_DELTA.md`.

⚠ **Crash rule (this WSL2 box):** run `copy_assign` **FOREGROUND, serial (no `--region-threads`), small batches**, outputs under `/home/juanfra/winloci_scratch/` (never `/tmp` — wiped on session restart). No `nohup`, no background waiters, no `pkill -f`.

1. **Decisive sim** — `bench/sim_k0_flank.py` already plants copy **B** (2%-diverged, MAPQ 60, currently dropped). Run `copy_assign` on its region twice (`--min-copies 2`, with and without `--homology-primary`) and report:
   - copies detected (expect **3 → 4**, B admitted);
   - B's reads: `tied` count (expect **98 → ~0**, now assigned);
   - total tied mass before/after.
2. **Real families** — 2–3 regions from `/home/juanfra/winloci_scratch/silver/probe6.txt`, same two modes. Report per region: families, copies, reads, tied fraction, and how many assignments **changed** (expected, by design).
3. Write `bench/HOMOLOGY_PRIMARY_DELTA.md`: both tables; state plainly that assignments shift because the copy set grew (stricter `α/(K−1)`); recommend (or not) a default flip based on the numbers. If the sim does **not** recover B, say so — a negative result is the finding.

- [ ] Commit: `bench(o2): homology-primary delta — copies recovered, spurious tied mass removed`

## Self-Review
- Spec coverage: oracle swap (T2), DRY core (T1), CLI + declared default change (T3), delta measurement (T4).
- Type consistency: `to_split_families(&[Vec<usize>], &[(usize,usize,f64)], &SplitParams) -> Vec<SplitFamily>` produced T1, consumed T2; `DenovoConfig.homology_primary` T2→T3.
- Placeholders: γ = 0.20 is stated explicitly; `threads` source in T2 marked "from cfg or the existing threads value" — the implementer must locate it, not invent a parameter.
