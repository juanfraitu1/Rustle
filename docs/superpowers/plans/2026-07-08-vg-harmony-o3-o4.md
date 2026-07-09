# VG harmony: O4 absent-copies → EM/copy-number, and O3 junctions into the EM — Plan

> REQUIRED SUB-SKILL: superpowers:subagent-driven-development.

**Goal:** Make the EM (O2) + reference-free copy number use the VG's ability to represent copies NOT in the
linear genome. (H1) Verify + TEST that O4 reference-absent copies flow into the EM assignment and `chi_H`/
`famcn_readonly` — so copies not in the genome are assigned and counted. (H2) Thread O3 copy-specific junctions
into the EM likelihood (currently PSV-only), so a copy distinguished by a NOVEL junction (a splice not in the
reference) is used for assignment and surfaces as a copy — harmonizing O2↔O3.

**Architecture:** Both reuse the existing objects. O4 `absent_copy::admit_candidate` already adds admitted
copies to the family copy set (`denovo_pipeline.rs:626`, `--absent-copies`) → they reach `copy_psv_alleles` →
the EM + `chi_h`. `read_copy_evidence` already computes a JUNCTION term; the EM just isn't fed junctions
(`em_assign_family` passes `junctions: vec![]`). H2 exposes per-copy + per-read junctions on `FamilyAssignment`
(parallel to `copy_psv_alleles`/`read_psv_obs`) and threads them.

## Global Constraints
- Additive & non-destructive: existing outputs byte-identical; new behavior only under existing opt-in flags
  (`--absent-copies`, `--em`).
- Reuse the existing junction likelihood in `read_copy_evidence` (no new junction model).

---

### Task H1: O4 → EM → copy-number harmony (verify + test)

**Files:** in-crate `#[cfg(test)]` (in `em_copy_assign.rs` or `readonly_copy_number.rs`); read
`src/rustle/vg_family/absent_copy.rs` + `denovo_pipeline.rs:620-645` to confirm the flow.

**Verify (report in the test's doc-comment):** trace that an `Admission::Copy(t)` from `absent_copy::
admit_candidate` (a divergent copy that does NOT remap to its host ≥98%) is pushed into the family's copies →
`build_family_profiles` → `copy_psv_alleles`. So it is a first-class copy for the EM and `chi_h`.

**Test `absent_copy_is_assigned_and_counted`:** simulate the *result* of O4 admission — a family whose
`copy_psv_alleles` includes a copy whose distinguishing alleles are NOT shared with the others (the admitted
absent copy). Generate reads carrying that copy's alleles. Assert: (a) `chi_h(&copy_alleles)` counts it (the
absent copy adds a color); (b) `em_assign_family(...)` assigns those reads to it (`Certified`, argmax = the
absent copy) and its abundance > 0. This proves a copy not in the reference, once admitted by O4, is assigned
by the EM and counted by the reference-free copy number. Commit.

*(If a full genome sim with a planted absent copy is cheap to add via `bench/sim_genome.py`, add an end-to-end
variant; otherwise the flow test above + the verified wiring is the deliverable — note the choice in the report.)*

### Task H2: O3 copy-specific junctions into the EM

**Files:** Modify `src/rustle/vg_family/denovo_pipeline.rs` (FamilyAssignment fields), `copy_assign_pipeline.rs`
(populate), `src/rustle/vg_family/em_copy_assign.rs` (`em_assign_family` signature), `src/bin/copy_assign.rs`
(pass them). Test in-crate.

**Interfaces:**
- Add to `FamilyAssignment`: `pub copy_junctions: Vec<Vec<i64>>` (per copy, parallel to `copy_psv_alleles`) and
  `pub read_junctions: Vec<Vec<i64>>` (per read, parallel to `read_psv_obs`). Populate them in the assign path
  from the already-computed per-copy `CopyProfile.junctions` (`copy_boundaries`) and per-read junctions
  (`read_features`/`read_introns`). Default empty when not computed (byte-identical to today for consumers that
  ignore them).
- Change `em_assign_family(read_obs, copy_alleles, params, eps, max_iter)` →
  `em_assign_family(read_obs, copy_alleles, read_junctions: &[Vec<i64>], copy_junctions: &[Vec<i64>], params,
  eps, max_iter)`. Build `CopyProfile { alleles, junctions: copy_junctions[k] }` and `ReadFeatures { psv_obs,
  junctions: read_junctions[i] }` (empty vecs if the junction slices are empty → identical to current behavior).
  `read_copy_evidence` already folds junctions into `logl`, `n_decisive`, and `min_p`, so no likelihood change.
- The binary's `--em` call passes `&fa.read_junctions`, `&fa.copy_junctions`.

**Test `junction_only_copy_is_resolved_by_em`:** two copies with IDENTICAL PSV alleles (PSV-only → K=0/SoftZone)
but DIFFERENT junction vectors (one has a copy-specific intron boundary). Reads carrying the distinguishing
junction. Assert: with empty junctions, `em_assign_family` labels them `SoftZone` (PSV-only can't separate); with
the real `copy_junctions`/`read_junctions` threaded, the same reads become `Certified` and argmax the correct
copy. This proves a copy defined by a NOVEL junction (not a PSV) is now resolvable — the O2↔O3 harmony.

**Steps:** TDD the two-mode test first (RED on the second mode). Wire fields + populate + thread + binary. Verify
`--em` still additive (default outputs byte-identical; the two new FamilyAssignment fields are additive).

**Verify (both tasks):** `PATH=/home/juanfra/miniforge3/bin:$PATH cargo test --lib` green; the O4 flow test +
the junction-resolution test pass; `cargo build --release --bin copy_assign` compiles; default outputs unchanged.
Commit each task separately.

### Task H3: Doc

Append a "VG harmony (O3/O4)" section to `bench/READONLY_COPY_NUMBER.md` (or a new `bench/VG_HARMONY.md`):
copies not in the linear genome are found two ways — (1) COLLAPSED via `depth_cn` (already), (2) DIVERGENT via
O4 admission (non-remapping) → EM + `chi_H`; O3 junctions now let novel-junction copies be resolved by the EM.
State the remaining frontier honestly: reads that never map to the family region need genome-wide VG-alignment
(vg-giraffe/GraphAligner) — the linear-alignment step is the last reference-bias source, scoped as future work.
