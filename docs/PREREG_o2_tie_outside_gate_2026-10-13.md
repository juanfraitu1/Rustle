# PREREG — the tie-outside gate depends on region sweep order: the fix (2026-10-13)

**Found by** the exact-tie fix's held-out check (`docs/PREREG_o2_exact_tie_2026-10-13.md`, F3; register 1274). Two
human MCL23 molecules were `assigned` by the binary before the fix and `tied` after it, with an unchanged copy and real
margins (6.9 / 13.8). The same before binary on 1 region thread also gives `tied`.

**Root cause.** The §6gz gate says a molecule with an AS-tied placement outside every supplied family unit is never
`Assigned`.
- The registry of such molecules (`register_tie_outside`) is filled by each region's worker as it sweeps its records.
- The gate reads that registry inside the worker that certifies the molecule's family (family.rs, the read-star verdict).

So when the region holding the outside placement is swept after the family's region (serially in region order, or
concurrently under `--region-threads`), the gate misses. The call stays `assigned`, while the write-time column
`tie_outside_catalog` says 1. Which calls leak depends on order and thread timing. Assigned rows with
`tie_outside_catalog` = 1 in the shipped runs: gorilla_KB3781 1 / 624, chimp_PTR 1 / 54; human anchor run 1 / 675.

**Fix.**
- The same demotion as the worker's gate (`tie_row`: `Tied`, not resolvable, uniform posterior), plus the
  `families.tsv` recount (`recount_after_union`), is re-applied to every region's families once ALL regions have been
  computed: `denovo_pipeline::apply_tie_outside_gate`, called in `copy_assign` before the cross-family pass and the
  drain.
- Every emitted status therefore sees the complete registry: assignments, posterior, psv_reads, phased reads, the
  placement fallback.
- Quant is unaffected: abundance comes from the per-family EM over PSV observations, and hard counts ignore status.
- Skipped under `--union-certificate`, whose pass scores outside placements itself.

**Not fixed (documented residual).** With `--families --gtf`, GTF routing by O2 status runs inside the worker. It can
still see a leaked `Assigned`. The pipeline never uses this combination: its assemble stage is `--assemble-only`.

**Bars.**
- **G1 (unit).** `apply_tie_outside_gate` demotes exactly the `assigned` rows whose read is outside, sets the worker
  gate's fields, recounts the family, and leaves other rows and families untouched. Red first.
- **G2 (determinism, real data).** The human anchor batch b4 (MCL23's batch) with the fixed binary on 1 and on 4
  region threads gives byte-identical `assignments.tsv`, and both MCL23 rows are `tied`. Before the fix these
  differed.
- **G3 (no other change, real data).** gorilla_KB3781 full assign, fixed binary against the exact-tie binary, same
  args (`--region-threads 4`): `assignments.tsv` differs only in rows that go assigned -> tied and have
  `tie_outside_catalog` = 1. Their count is reported (the exact-tie run had 0 such rows at 4 threads, so 0 changes
  are expected there). Then the fixed binary on 1 region thread (the pipeline default) must give 0 assigned rows
  with `tie_outside_catalog` = 1. The shipped run had 1.
- **G4.** `cargo test --release --no-fail-fast` passes, except the 3 known `pipeline_driver_dispatch` failures.

Outputs `/mnt/linuxdisk/tmp/o2_tie/`.

## Amendment 1 (2026-10-13, after G2's first run, before G3)

**G2's first run.** `assignments.tsv` was identical on 1 and 4 threads, but `families.tsv` and `quant.tsv` were not,
in MCL23 only. The statement "quant is unaffected" above was **wrong**, for two reasons:
- On the read-star path the abundance is `star_soft_abundance`, the normalised sum of the molecules' posteriors.
  The worker computes it after its own gate, so an in-time demotion (uniform posterior) moves it.
- `uniq` / `uniq_agree` are counted in the worker from the then-current status, and the recount does not rebuild them.

**The fix now also:**
- records each row's unique-mapper contribution, `FamilyAssignment::uniq_contrib`, and subtracts it when a row is
  demoted late;
- recomputes the read-star soft abundance for families with a demoted row. This uses the same function
  (`star_soft_abundance_of`, which `star_soft_abundance` now calls), the same row order and the same width, so it is
  bit-identical to what the worker computes when its gate fires in time. It is skipped under the VG re-align leg,
  whose own EM replaces the abundance.

**Bars, extended.**
- **G2:** every output identical on 1 and 4 threads.
- **G3 adds:** the fixed binary on the full gorilla run gives byte-identical `assignments` / `families` / `quant` on
  1 and 4 threads.

## Results

(pending)
