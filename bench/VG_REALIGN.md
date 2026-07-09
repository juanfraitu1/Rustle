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

1. **Candidate** (`is_candidate`): a read whose linear fit is poor or absent — low MAPQ, high soft-clip
   fraction, or high per-base divergence (`BamRead::de`) to its locus. (Unmapped reads: `unmapped_reads_from_bam`
   + `route_unmapped` minimizer pre-screen route them to candidate families.)
2. **Re-align** (`realign_to_paths`): score the read against every copy-path with `bridge_detector::aln_id`
   (edlib infix identity, tries both strands); `best_copy` = argmax, `id_best`, and `id_linear` = fit to the
   read's current linear copy.
3. **Accept** (`accept_realignment`): re-assign only if the re-alignment is **significant** — the same `ε^δ`
   certificate the gate/EM use, here `min_p = (error_rate/3)^n_decisive < alpha`, with `n_decisive =
   round((id_best − id_linear)·read_len)` = the positions where the read matches the best copy-path but not the
   linear locus. No tuned margin. `Reassign(best)` vs `Reject`.
4. **Novel candidates** (`pool_novel`): reads that fit NO existing copy-path are pooled by mutual identity;
   a pool of ≥ `min_reads` is a candidate reference-absent copy.

`copy_assign --vg-realign` emits `<out>.vg_realign.tsv`
(`read_name  family_id  action  target_copy  id_best  linear_copy`), `action ∈ {reassigned, rejected,
novel-candidate}`.

## What's shipped vs. what's next (honest scope)

- **Shipped (T1–T5, 822 tests):** the full candidate → re-align → significance-accept → novel-pool logic, and
  the pipeline wiring that RUNS it per family and REPORTS the decisions in `.vg_realign.tsv`. Additive: with
  `--vg-realign` off, nothing runs and every existing output is byte-identical (fixture-verified). The unit
  tests demonstrate the correction (a read from copy 1 linearly mis-placed on copy 0, low MAPQ → `reassigned`
  to copy 1) and the candidate gating.
- **Follow-up (deliberately not wired):** (a) routing accepted `reassigned` reads back into the EM assignment /
  `chi_H` / `famcn_readonly` (so a correction changes the counts, not just the report) — invasive re-entry into
  the assignment, kept out to protect the byte-identical guarantee; (b) admitting `novel-candidate` pools as
  reference-absent copies via `pool_novel` + `absent_copy::admit_candidate` (the O4 gate) so a discovered copy
  enters the copy set end-to-end; (c) refining the identity-based accept to per-PSV `path_obs` +
  `read_copy_evidence` (an edlib alignment-path walk `aln_id` doesn't currently expose).
- **Frontier (still future):** genome-wide VG-alignment (vg-giraffe/GraphAligner) for reads that never map near
  any family — the last linear-alignment reference-bias source.

Relates to `bench/VG_HARMONY.md`, `bench/READONLY_COPY_NUMBER.md`, `absent_copy.rs`, `family_graph.rs`.
