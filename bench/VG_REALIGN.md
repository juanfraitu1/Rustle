# VG re-align supplement

**Date:** 2026-07-08. A bounded, per-family supplement that re-aligns the reads a linear aligner mis-placed or
dropped to O1's family **copy-paths**, and reports — under a significance certificate — which reads should be
re-assigned and which suggest a copy not in the linear reference. Opt-in (`--vg-realign`), default output
byte-identical.

## Where the non-linear reference comes from — O1

The copy-paths are O1's: a family's copies (their spliced consensuses, `DenovoTranscript::seq`) are the paths of
the `FamilyGraph`. The supplement aligns a candidate read to each copy-path (not just its one linear locus), so
a read the linear reference biased onto the wrong copy — or failed to place — can find its true copy-path. O2
(the EM) and the reference-free copy number consume the result. This is **not** a genome-wide graph aligner:
only candidate reads, only the family's own copy-paths.

## Mechanism (all in `src/rustle/vg_family/vg_realign.rs`)

1. **Candidate** (`is_candidate`): a read (among the family's region-mapped reads) whose linear fit is poor —
   low MAPQ, high soft-clip fraction, or high per-base divergence (`BamRead::de`) to its locus. *(Unmapped-read
   routing — `unmapped_reads_from_bam` + `route_unmapped` minimizer pre-screen — is BUILT and unit-tested but
   NOT yet called from the pipeline; the wired path only sees region-mapped reads. See follow-up (a).)*
2. **Re-align** (`realign_to_paths`): score the read against every copy-path with `bridge_detector::aln_id`
   (edlib infix identity, tries both strands); `best_copy` = argmax, `id_best`, and `id_linear` = fit to the
   read's current linear copy.
3. **Accept** (`accept_realignment`): re-assign only if the re-alignment is **significant** — the same `ε^δ`
   certificate the gate/EM use, here `min_p = (error_rate/3)^n_decisive < alpha`, with `n_decisive =
   round((id_best − id_linear)·read_len)` = the positions where the read matches the best copy-path but not the
   linear locus. No tuned margin. `Reassign(best)` vs `Reject`.
4. **Novel candidates:** a candidate read that fits NO existing copy-path (`realign_to_paths` → `None`) is
   flagged **per-read** as `novel-candidate` in the report. *(Pooling those unfit reads by mutual identity —
   `pool_novel`, ≥ `min_reads` = a candidate reference-absent copy — is BUILT and unit-tested but NOT yet wired;
   see follow-up (b).)*

`copy_assign --vg-realign` emits `<out>.vg_realign.tsv`
(`read_name  family_id  action  target_copy  id_best  linear_copy`), `action ∈ {reassigned, rejected,
novel-candidate}`.

## What's shipped vs. what's next (honest scope)

- **Shipped END-TO-END (833 tests):** under `--vg-realign`, corrections and admissions now **change the
  outputs**, not just the report:
  - **Corrections:** a `reassigned` read's PSV evidence is re-extracted at the corrected copy-path via a
    strand-oriented **traceback aligner** (`align_traceback` + `path_obs_at`; there is no edlib, so the
    alignment path is hand-rolled and validated against `hw_distance`), its full `Assignment` is **re-derived**
    (`assign_read_editing` over the corrected obs — status/posterior/p_value all consistent), and the family's
    `copy_abundance` is recomputed by the EM. So `.assignments.tsv`/`.em_abundance.tsv`/`.posterior.tsv` reflect
    the correction.
  - **Admissions:** `pool_novel` clusters of unfit reads are turned into `CollapsedCandidate`s and admitted via
    `absent_copy::admit_candidate` (the O4 remap gate); an admitted copy is appended to the family copy set and
    counted by `chi_H`/`famcn_readonly`, with `copy_abundance` widened + recomputed over the new roster.
  - **Additive & reviewed:** with `--vg-realign` off, nothing runs and every existing output is byte-identical
    (verified by a real-fixture binary diff + the `vg_realign_off_is_byte_identical` test). A whole-feature
    review caught + fixed a minus-strand `path_obs` corruption, an abundance-length desync on admission, stale
    `Assignment` fields on correction, and a `.posterior.tsv` width issue.
  - **`psv_positions_for`** derives each family PSV column's offset in a copy's spliced consensus by inverting
    `build_family_profiles`' `copy_gpos` through `gen2off` (the consensus used for alignment IS the one the
    offsets index).
- **Follow-up (not wired):** (a) routing UNMAPPED reads into families (`unmapped_reads_from_bam` +
  `route_unmapped` — the helpers exist, unit-tested) so out-of-region reads become candidates; (b) the
  end-to-end genome sim (`bench/sim_vg_realign.py`, plant a divergent mis-mapping + a genome-absent copy and
  assert corrected assignments + admitted-and-counted copy) — validation currently rests on the unit tests, the
  real-fixture additivity diff, and the whole-feature review; (c) admission yield on real GGO is expected to be
  data-limited (the O4-divergent frontier); the mechanism is wired, real yield reported not assumed.
- **Frontier (still future):** genome-wide VG-alignment (vg-giraffe/GraphAligner) for reads that never map near
  any family — the last linear-alignment reference-bias source.

Relates to `bench/VG_HARMONY.md`, `bench/READONLY_COPY_NUMBER.md`, `absent_copy.rs`, `family_graph.rs`.
