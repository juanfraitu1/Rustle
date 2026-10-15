# Handoff — laptop → main machine (2026-10-15)

Work done on the second machine (WSL laptop, 6 cores / 9 GB, no `/mnt/linuxdisk`) lives on the branch
**`laptop/o3-psv-multimappers`**, on top of `main` = `49e62be3`. It is NOT merged. `main` was not touched from the laptop.

## What is on the branch

| commit | what |
|---|---|
| `797f55f8` | `pixi.toml` / `pixi.lock` (minimap2 2.30, samtools 1.22.1, bench and IsoCon envs); `.gitignore` exceptions for the two fixture sets that were never committed (below) |
| `5022350c` | O3 positive control re-measured (`docs/PREREG_missing_copy_posctl_2026-10-14.md`, PASS); few-read false flags diagnosed (`docs/PREREG_o3_fewread_flags_2026-10-14.md`: F = 0 on simulation); exon-structure projection in `candidate_copies` (`structure.tsv`, `contigs.gtf`; `docs/PREREG_o3_structure_projection_2026-10-14.md`, PASS) |
| `efc84916` | O3 per the advisor's definition (cluster the reads O2 places on a family by their PSVs; a cluster no copy explains = a copy absent from the genome). Prototype `bench/o3_unexplained/psvclust.py`; `docs/PREREG_o3_unexplained_2026-10-15.md`, PASS |
| `fdebaec0` | Rust port: **`copy_assign --psv-clusters`** (`rustle::missing_copy::psv_clusters`, OPT-IN), v2 with the linkage test against RNA editing and shared errors; `docs/PREREG_o3_psv_clusters_v2_2026-10-15.md`, all bars PASS (425/426 absent copies, 0/152 controls flagged, 16/16 runs in parity with the prototype) |
| this file | the steps below |

Design notes: `docs/superpowers/specs/2026-10-14-o3-structure-projection-design.md`,
`docs/superpowers/specs/2026-10-15-o3-unexplained-multimappers-design.md`. `REPRODUCE.md`, `bench/README.md` and the O3 row of
`docs/THESIS_OBJECTIVES.md` are updated on the branch.

## Steps on the main machine, in order

Use the canonical repo's usual build target (`CARGO_TARGET_DIR=...` on `/mnt/linuxdisk`) and the usual locks.

### 1. Fetch the branch

```bash
cd /mnt/c/Users/jfris/Desktop/Rustle        # the canonical clone
git fetch github laptop/o3-psv-multimappers   # or `origin`, whichever remote points at juanfraitu1/Rustle
git switch -c laptop/o3-psv-multimappers --track github/laptop/o3-psv-multimappers
```

If `main` moved on the main machine since `49e62be3`, rebase the branch on it (or merge `main` into it) before step 3.

### 2. Commit the test fixtures that exist only on the main machine

`.gitignore` ignores `*.gtf` and `*.fa`, so these were never committed. A fresh clone fails 10 tests without them. The branch
already has the `.gitignore` exceptions (`!/tests/fixtures/bridge_units/*.gtf`, `!/tests/fixtures/o2_origin_resolve/*.fa`).

```bash
git add tests/fixtures/bridge_units/plain.gtf \
        tests/fixtures/bridge_units/expected.units.gtf \
        tests/fixtures/bridge_units/s_f0.5.plain.gtf \
        tests/fixtures/bridge_units/s_f0.5.expected.units.gtf \
        tests/fixtures/o2_origin_resolve/copies.fa \
        tests/fixtures/o2_origin_resolve/genome.fa
git status --short tests/fixtures      # anything else the tests read that is still untracked? add it too
```

Check sizes first (`du -sh`); all should be small synthetic files.

### 3. Run the full suite

```bash
cargo build --release --bins
cargo test --release 2>&1 | tee /tmp/rustle_test_laptop_branch.log
```

Expected on the laptop: 1,299 passed, 13 failed. The 13 are known:
- 10 are the missing fixtures (step 2). They should pass here.
- 3 are stale tests in `tests/pipeline_driver_dispatch.rs`. The driver's `all)` branch now calls `run_all`, which the test does not
  stub, so it exits 127. Fix the test (stub `run_all`) or record it as known.

Anything else failing is new and must be looked at before merging. `copy_assign` without `--psv-clusters` is byte-identical to
before on the simulation runs, but that was not checked on real data.

### 4. Merge to `main` and push

Only after step 3 is clean (or only the 3 dispatch tests remain and are recorded):

```bash
git switch main && git merge --no-ff laptop/o3-psv-multimappers
git push github main
```

Add the laptop results to `docs/NEGATIVE_RESULTS_REGISTER.md` (the posctl, few-read, structure-projection, v1 and v2 preregs above).
They were not given row numbers on the laptop to avoid clashing with rows added on the main machine.

### 5. M3 foreign-species check (pre-registered, not run)

`docs/PREREG_o3_candidates_foreign_check_2026-10-14.md`. Inputs are already on `/mnt/linuxdisk/tmp/rna_allele/{a14a15,a15h}/`.
It asks whether the remaining `candidate_copies` false flags best match another ape. Treat any result as evidence about the reads,
not about how they entered the library. That is not established.

### 6. O3 on real data: pre-register first, then run

The simulation level is done. The next level, from the v2 prereg's reading, is:

- Pick real multi-copy families from the gorilla (KB3781) catalog. **Remove copies from the GENOME** (not just the catalog, which
  is what the 2026-09 excision experiments did, so reads stayed on the removed locus), re-map the real reads with the exact O2
  settings, run `copy_assign --families ... --copies-fa ... --no-as-tied-only --flag-missing-copies --psv-clusters`. The removed
  copy's reads should form an `unexplained` cluster on its nearest paralog. Then add the copy back: the cluster must disappear.
- **Single-copy-gene control**: run the same on single-copy genes. Any `unexplained` cluster labelled `copy` there is the false
  rate the real heterozygous allele background gives. Simulation cannot measure this.
- Write the prereg (bars, families, held-out split) and commit it before running.

Run with `--no-as-tied-only`. With exact AS ties only (the shipped gate), O3 sees only identical-twin cases. The absorbed reads of
a missing copy are mostly unique mappers.

### 7. Smaller open items

- `o3_export.py` does not yet read `candidate_copies`' `contigs.gtf` (secondary under the advisor's framing).
- HG002 plan: `docs/PENDING_2026-10-11.md`.

## Laptop notes (if work continues there)

- pixi envs: `pixi run <cmd>` (default), `pixi run -e bench`, `pixi run -e isocon`. Build with `-j 4`
  (`.cargo/config.toml`'s `jobs = 48` is for the main machine).
- Data: `~/rustle_data/` (CHM13 chr20 + RefSeq annotation, the simulation worlds). Real BAMs will come from the NAS.
- The laptop keeps committing to `laptop/*` branches, never to `main`.
