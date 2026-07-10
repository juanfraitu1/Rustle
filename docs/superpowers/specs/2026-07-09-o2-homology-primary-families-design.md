# O1↔O2 harmony: E_r homology-primary family membership in `copy_assign` — Design

**Date:** 2026-07-09. **Substrate:** gorilla (GGO) HiFi Iso-Seq; `denovo_pipeline::detect_and_assign` (O2).

## The defect

`copy_assign` (O2) forms families with **`read_conflict::conflict_edges` / `conflict_families`** — the de-tie
conflict graph (**E_c**), which links two copies only when reads map ambiguously between them. The **E_r
homology-primary** definition (`homology_edges_all_reps`, `gamma_quasi_clique_partition`) — built precisely to
fix this — is reachable **only from `gw_family_catalog`**; `copy_assign` never references it.

**Consequence:** a copy whose reads map uniquely (MAPQ high) forms no conflict edge and is **dropped from the
family**. Its reads are still pulled into the family region and come back **`tied`** — not because they are
unassignable, but because *their true copy was never admitted to the copy set*. Demonstrated in
`bench/sim_k0_flank.py`: copy **B** (2%-diverged, MAPQ 60) is dropped; its **98 reads are 100% tied**.

**A family is defined by homology, not by read ambiguity.** Easily-assignable members are still members.

This makes a **third cause of abstention** (documented in `bench/K0_FLANK_EXPERIMENT.md`):

| cause | mechanism | fixable |
|---|---|---|
| K = 0 | no distinguishing column exists | No — proven wall |
| Coverage | columns exist, read doesn't span them | Yes — longer reads |
| **Missing copy** | true copy never admitted (E_c drops uniquely-mappable members) | **Yes — this spec** |

## Goal

Give `copy_assign` an opt-in **`--homology-primary`** mode in which family **membership** is the E_r homology
component, so O1's fixed definition finally flows into O2. Measure the delta; flip the default only afterwards.

## Design — swap the edge oracle, nothing else

```
today (E_c):     reps ─► conflict_edges ─► conflict_families ────────┐
--homology-primary: reps ─► homology_edges_all_reps (E_r)            ├─► SplitFamily ─► colocated_families ─► copies ─► assign
                          ─► gamma_quasi_clique_partition(γ)         ┘
```

- **Membership only changes.** `colocated_families` (co-location, per-`(chrom,strand)` split, `win`,
  `min_copies`), PSV discovery, assignment, the EM, `chi_H`, `famcn_readonly` are **untouched**.
- **Conflict/PSV/χ(H) remain WITHIN the family** — copy-number evidence and assignment, never a membership
  oracle. This is exactly the recorded E_c→E_r reframe.
- **Pure E_r when on.** Conflict edges are NOT unioned in. (E_c ⊆ E_r operationally; the ~3 "shared-exon leak"
  edges E_c finds but homology rejects are, by this definition, not families.)

### Interfaces (all exist today)
- `homology_edges_all_reps(reps: &[DenovoTranscript], params: &RefineParams) -> Result<Vec<(usize, usize)>>`
  (`pub(crate)`, denovo_pipeline.rs:1739) — exon-sum nt homology (asm20 ∪ sensitive tier).
- `homology_refine_params(min_identity: Option<f64>, threads: usize) -> RefineParams` (denovo_pipeline.rs:1446).
- `gamma_quasi_clique_partition(n: usize, edges: &[(usize, usize, f64)], gamma: f64) -> Vec<Vec<usize>>`
  (family_split.rs:477) — γ = 0.20, the same operator O1 uses (locate the existing constant; do not invent one).
- `conflict_to_split_families(families, c_edges: &[(usize,usize,usize)], p: &SplitParams) -> Vec<SplitFamily>`
  (denovo_pipeline.rs) — **generalize**: extract a core taking `&[(usize,usize,f64)]` so both the conflict path
  (weights = read counts) and the homology path (weights = 1.0) share it. DRY; no duplicated `community_stats`/
  `classify` logic.

### Config / CLI
- `DenovoConfig.homology_primary: bool` (default `false`).
- `copy_assign --homology-primary` sets it. Reuse `RefineParams` defaults via `homology_refine_params`
  (`RUSTLE_MINIMAP2` honored).
- **`--min-copies` default changes `3 → 2`.** This is an **independent, deliberate change** (2-copy homologous
  families are the majority and are currently invisible to O2). It alters default family detection on its own,
  separate from the E_r flag — stated here so it is not mistaken for a side-effect.

### Failure behavior (no silent degradation)
`homology_edges_all_reps` shells minimap2 and returns `Result`. If it fails (binary absent, alignment error),
`--homology-primary` must **abort with a clear error**, never silently fall back to E_c. (Precedent: the famCN
`unwrap_or_default()` silent-degradation bug.)

## Consequence, stated plainly

Admitting a dropped copy **enlarges the copy set** → more PSV columns → a stricter Bonferroni `α/(K−1)`
certificate. **Existing assignments will shift.** This is intended: a more complete copy set is more correct.
It is why the mode is opt-in and measured before any default flip. **`--homology-primary` off ⟹ the E_c path is
literally untouched ⟹ byte-identical** (modulo the separate `min_copies` default change above).

## Validation — the delta IS the deliverable

1. **Decisive sim** (`bench/sim_k0_flank.py`, already planted with a uniquely-mappable copy B):
   - *without* the flag: 3 copies, **B dropped**, B's **98 reads 100% `tied`**.
   - *with* the flag: **4 copies, B admitted**, B's reads **assigned** (not tied); tied mass drops by ~98.
   This isolates cause-3 (missing copy) from causes 1–2 (K=0, coverage) — the sim already contains the fixture.
2. **Real families** — a few probe regions (`winloci_scratch/silver/probe6.txt`), **foreground, serial, small
   batches, outputs under `winloci_scratch`** (background sweeps crash this WSL2 box; see the crash rule).
   Report: copies recovered, families gained, tied-mass reduction, and any assignment that *changed* (expected).
3. Emit `bench/HOMOLOGY_PRIMARY_DELTA.md` with both tables + the honest note that assignments shift by design.

## Non-goals

- No default flip in this spec (decide after the delta).
- No union with E_c; no change to `colocated_families`, PSV discovery, assignment, EM, `chi_H`, or `famcn_readonly`.
- Not genome-wide (that is `gw_family_catalog`'s job); this is per-region within `copy_assign`.

## Files

- **Modify** `src/rustle/vg_family/denovo_pipeline.rs`: `DenovoConfig.homology_primary`; branch the edge oracle
  in `detect_and_assign`; generalize `conflict_to_split_families` to an f64-edge core.
- **Modify** `src/bin/copy_assign.rs`: `--homology-primary`; `--min-copies` default `3 → 2`.
- **Create** `bench/homology_primary_delta.py` (runs both modes on the sim + probe regions),
  `bench/HOMOLOGY_PRIMARY_DELTA.md`.
- **Test** (in-crate): a uniquely-mappable homologous rep that `conflict_edges` misses IS admitted by the
  homology oracle; flag-off leaves the E_c families byte-identical; `homology_edges_all_reps` failure aborts.

## Reproduce

```
copy_assign --bam k0.bam --fasta k0.fasta --region <r> --min-copies 2 --skip-poa-diagnostic --out e_c
copy_assign --bam k0.bam --fasta k0.fasta --region <r> --min-copies 2 --homology-primary --skip-poa-diagnostic --out e_r
python bench/homology_primary_delta.py     # copies recovered, tied-mass reduction
```
