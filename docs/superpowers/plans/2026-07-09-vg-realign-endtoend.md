# VG re-align end-to-end wiring — Implementation Plan

> REQUIRED SUB-SKILL: superpowers:subagent-driven-development.

**Goal:** Under `--vg-realign`, make corrections and admissions change the outputs: a corrected read's PSV
evidence is re-extracted at the right copy-path (new traceback aligner) and its assignment/abundance updated;
a novel-candidate pool is admitted as a reference-absent copy (O4) and counted. Default (flag off) byte-identical.

**Architecture:** `assign_family_detailed` derives read PSV evidence from the LINEAR alignment, so corrected/
admitted reads (whose linear alignment is wrong for their true copy) cannot be fixed by a plain re-run — the
`path_obs` must OVERRIDE the linear evidence. Implement as an in-place **patch** of the `FamilyAssignment`
(corrected reads' `read_psv_obs` + `assignments`; appended admitted copies + their reads) followed by an EM
abundance recompute — the bounded family-local update of the spec's §4. Reuses the new traceback aligner,
`vg_realign::{realign_to_paths, accept_realignment, pool_novel}`, `absent_copy::admit_candidate`, `em_assign_family`.

## Global Constraints
- Additive: `cfg.vg_realign` false ⟹ no patch runs ⟹ `.assignments.tsv`/`.families.tsv`/`.quant.tsv`/`.em*.tsv`/
  `.famcn_readonly.tsv` byte-identical. All changes live inside the `if cfg.vg_realign` block.
- No new admission model: novel copies go through `absent_copy::admit_candidate` (the O4 remap gate) verbatim.
- The traceback aligner MUST agree with `bridge_detector::hw_distance` on edit distance (tested).

## File Structure
- `src/rustle/vg_family/vg_realign.rs` — MODIFY: `align_traceback`, `path_obs_at`, `apply_realign`.
- `src/rustle/vg_family/denovo_pipeline.rs` — MODIFY: patch `FamilyAssignment` under `cfg.vg_realign`.
- `bench/sim_vg_realign.py` — CREATE; `bench/VG_REALIGN.md` — UPDATE.

---

### Task 1: Traceback aligner + path_obs

**Files:** `src/rustle/vg_family/vg_realign.rs`. Test in-crate.

**Interfaces:**
- `pub fn align_traceback(query: &[u8], target: &[u8]) -> Vec<(Option<usize>, Option<usize>)>` — infix/HW
  alignment (free end-gaps on `target`) via an O(|query|·|target|) DP with a backtrack matrix; returns the
  ordered aligned pairs `(query_idx, target_idx)` (`None` on the gapped side). Choose the best-scoring end
  column (min edit distance), backtrack to the first query base.
- `pub fn path_obs_at(align_map: &[(Option<usize>, Option<usize>)], psv_positions_in_consensus: &[usize], query: &[u8]) -> Vec<Option<u8>>`
  — for each PSV position `t` (index into `target`), find the pair whose `target_idx == Some(t)`; return
  `Some(query[qi])` if that pair's `query_idx == Some(qi)`, else `None` (read gaps / doesn't span it). One entry
  per `psv_positions_in_consensus`, in order.

- [ ] **Test `traceback_edit_distance_matches_hw`:** on the same pairs `hw_distance` is tested on (exact infix, substitutions, indels, empty), the count of non-matching aligned columns + gaps in `align_traceback`'s output equals `hw_distance(query, target)`. (Mirror `bridge_detector`'s `hw_distance_corners` corpus.)
- [ ] **Test `path_obs_reads_bases_at_psv_columns`:** `query = target` with a 1-base substitution at consensus position 5 → `path_obs_at(map, &[5, 10], query)` returns `[Some(query[5]), Some(query[10])]`; a query that ends before position 20 → `path_obs_at(map, &[20], query)` returns `[None]`.
- [ ] Implement; GREEN; commit `feat(vg_realign): align_traceback + path_obs_at (read PSV bases at a copy-path)`.

### Task 2: `apply_realign` orchestrator

**Files:** `src/rustle/vg_family/vg_realign.rs`. Test in-crate.

**Interface:**
```
pub struct RealignApply {
    pub corrected: std::collections::HashMap<usize, (usize, Vec<Option<u8>>)>, // read_idx -> (new_copy, new_read_psv_obs)
    pub admitted: Vec<DenovoTranscript>,                                        // novel reference-absent copies
    pub records: Vec<RealignRecord>,                                           // the report rows (as today)
}
pub fn apply_realign(
    bam_reads: &[BamRead], copies: &[DenovoTranscript], copy_seqs: &[Vec<u8>],
    psv_pos_per_copy: &[Vec<usize>],           // family PSV positions in each copy's consensus coords
    linear_copy_of: &[Option<usize>],          // per read, its linear copy (best_overlap_copy)
    genome: &GenomeIndex, fasta_path: &str,
    rp: &RealignParams, error_rate: f64, alpha: f64,
) -> RealignApply
```
Logic: for each candidate read (Task-1 supplement's `is_candidate`): `realign_to_paths` → if `Some(hit)`,
`accept_realignment` → on `Reassign(best)`, compute `path_obs_at(align_traceback(read_seq, copy_seqs[best]),
psv_pos_per_copy[best], read_seq)` and insert `corrected[read_idx] = (best, path_obs)` + a `reassigned` record;
on `Reject`, a `rejected` record. Reads with `realign_to_paths == None` collect into the unfit pool. Then
`pool_novel(unfit, min_id, rp.min_reads)`; for each pool, build a `CollapsedCandidate` from the pool consensus
(reuse `copy_split`'s allele-vector construction / `recover_collapsed_candidates`' shape) and call
`absent_copy::admit_candidate(cand, host, genome, fasta_path, &AbsentCopyParams::default())`; push
`Admission::Copy(t)` into `admitted` + a `novel-copy` record (else `novel-candidate`/rejected).

- [ ] **Tests:** `apply_corrects_mismapped_read` — a read carrying copy 1's alleles but `linear_copy_of = Some(0)` with low MAPQ → `corrected` contains `read_idx -> (1, obs)` where `obs` matches copy 1's alleles. `apply_admits_novel_pool` — ≥`min_reads` mutually-similar unfit reads that don't remap (mock/params so `admit_candidate` returns `Copy`) → one `admitted`. `apply_no_candidates_noop` — clean reads → empty everything. Commit `feat(vg_realign): apply_realign (corrections via path_obs + O4 admissions)`.

### Task 3: Patch the FamilyAssignment under `--vg-realign`

**Files:** `src/rustle/vg_family/denovo_pipeline.rs`. Test additivity in-crate.

In `detect_and_assign`, inside the existing `if cfg.vg_realign` block (where `run_family_realign` was called):
replace the report-only call with `apply_realign(...)` and PATCH the family's `FamilyAssignment fa`:
1. **Corrections:** for each `(read_idx, (new_copy, obs))` in `corrected`: set `fa.read_psv_obs[i] = obs` and
   update `fa.assignments[i]`'s assigned copy to `new_copy` (find the position `i` whose `assignments[i].0 ==
   read_idx`).
2. **Admissions:** append each `admitted` copy's PSV allele vector to `fa.copy_psv_alleles`, its junctions to
   `fa.copy_junctions`, its tid to `fa.copy_tids`; bump `fa.n_copies`; assign the pool's reads to the new copy
   index (patch their `assignments`/`read_psv_obs` like a correction).
3. **Recompute abundance:** `fa.copy_abundance = em_assign_family(&fa.read_psv_obs, &fa.copy_psv_alleles,
   &fa.read_junctions, &fa.copy_junctions, &params, eps, max_iter).abundances` (the EM now sees corrected +
   admitted state). Keep `fa.realign_records = records`.
`chi_H`/`famcn_readonly` are recomputed in the binary from `fa.copy_psv_alleles`/`fa.n_reads`, so appending an
admitted copy raises `chi_H` automatically. When `cfg.vg_realign` is false NONE of this runs.

- [ ] **Test additivity:** a fixture/synthetic family with `cfg.vg_realign = false` ⟹ `fa` unchanged vs the pre-patch build (all fields byte-identical). With a planted mis-mapped read + `true` ⟹ that read's `assignments` entry shows the corrected copy. Commit `feat(vg_realign): patch FamilyAssignment (corrections+admissions) under --vg-realign; recompute EM`.

### Task 4: End-to-end sim + doc

**Files:** `bench/sim_vg_realign.py` (extends `bench/sim_genome.py`), update `bench/VG_REALIGN.md`.

Plant a family with (i) a ~6–10% divergent copy some of whose reads mis-map to a sibling copy's locus, and (ii)
a genome-absent divergent copy (in reads, removed from the reference FASTA). Run `copy_assign` with and without
`--vg-realign`. Assert: (a) without the flag the mis-mapped reads land on the wrong copy in `.assignments.tsv`;
with it, on the correct copy. (b) without the flag the genome-absent copy is absent from the catalog; with it,
it appears as an admitted copy and `chi_H`/`famcn_readonly` count it. Write results into `bench/VG_REALIGN.md`
(move corrections/admissions from "follow-up" to "shipped"; keep the honest real-GGO caveat). Commit.

## Self-Review
- Spec coverage: path_obs (T1), corrections+admissions orchestration (T2), state patch + EM recompute (T3),
  sim (T4). Additivity in every wiring task.
- Type consistency: `align_traceback`/`path_obs_at` (T1) consumed by `apply_realign` (T2); `RealignApply`
  consumed by the patch (T3); reuses `RealignRecord`, `BamRead`, `DenovoTranscript`, `CollapsedCandidate`,
  `absent_copy::admit_candidate`, `em_assign_family`, `GenomeIndex`.
- Placeholder scan: the pool-consensus → `CollapsedCandidate` construction (T2) says "reuse
  `recover_collapsed_candidates`' shape" — the implementer must locate and mirror that existing construction,
  not invent a candidate format.
