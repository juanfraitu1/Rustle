# VG re-align end-to-end wiring — Design

**Date:** 2026-07-09. Follows the report-only VG re-align supplement (`bench/VG_REALIGN.md`,
`src/rustle/vg_family/vg_realign.rs`). Substrate: gorilla (GGO) HiFi Iso-Seq; the EM, `chi_H`/`famcn_readonly`,
O4 `absent_copy`.

## Goal

Turn the report-only supplement into state changes: under `--vg-realign`, a **corrected** read's PSV evidence
is re-extracted at the right copy-path and fed to the EM, and a **novel-candidate** pool is admitted as a
reference-absent copy via O4 — so `.assignments.tsv`, `.em*.tsv`, `chi_H`, and `famcn_readonly` all *reflect*
the corrections and discoveries. With `--vg-realign` off, nothing runs and every output is byte-identical.

## 1. `path_obs` via a traceback aligner (the one new algorithmic piece)

There is NO edlib crate — `bridge_detector::hw_distance` is a hand-rolled 2-row DP that yields edit *distance*
only, no alignment path. So `path_obs` needs a traceback alignment:
- `pub fn align_traceback(query: &[u8], target: &[u8]) -> Vec<(Option<usize>, Option<usize>)>` — an infix/HW
  alignment (free end-gaps on `target`) computed with a backtrack matrix, returning the column↔row map (aligned
  pairs; `None` on a gap side). Sequences are kb-scale (read vs copy consensus), so O(mn) with a byte backtrack
  matrix is acceptable per candidate read.
- **Correctness anchor:** the number of mismatches+indels in `align_traceback`'s path MUST equal
  `hw_distance(query, target)` — assert this in tests so the traceback can't silently diverge from the shipped
  distance function.
- `pub fn path_obs_at(align_map: &[(Option<usize>, Option<usize>)], psv_positions_in_consensus: &[usize], query: &[u8]) -> Vec<Option<u8>>`:
  for each family PSV column's position in the copy consensus (`target`), return the query (read) base aligned to
  it, or `None` if the read gaps/doesn't span it.

The family PSV columns' positions in the best copy's consensus coordinates come from the existing PSV discovery
(`copy_assign_pipeline::discover_psvs` / the per-copy PSV allele layout) — the caller supplies
`psv_positions_in_consensus` for the chosen copy.

## 2. Corrections → EM

For each read whose `accept_realignment` returned `Reassign(best)`: compute `path_obs` at `best`'s consensus
(§1) → this replaces the read's `read_psv_obs` entry. When the family assignment/EM re-runs (§4) with the
corrected `read_psv_obs`, the read's PSV evidence now matches its true copy, so the EM posterior + the hard
assignment place it there. `chi_H` is unaffected by a re-assignment (it counts copy hap-vectors, not reads);
`famcn_readonly`'s `depth_cn` is unaffected (total family depth is unchanged) — corrections move a read between
existing copies and are reflected in `.assignments.tsv`/`.em*.tsv`.

## 3. Admissions → copy set

`pool_novel` the reads that fit no existing copy-path → build a pool consensus (POA of the pool, or the longest
read as a seed with the others aligned) → construct a `CollapsedCandidate` (its `psv_pos` + `allele_vector` from
the consensus at the family PSV columns) → `absent_copy::admit_candidate(cand, host, p, remap)` with a **real
remap closure**: align the pool consensus to the genome and return its best identity to any existing locus;
`admit_candidate` admits it iff it is min_p-distinct, strand-symmetric, and does NOT remap ≥98% (genuinely
absent). Reuse the exact remap the O4 path in `denovo_pipeline` already constructs. Admitted synthetic copies are
appended to the family's copy set; the pool's reads then assign to them, and `chi_H` (+1 color) and
`famcn_readonly` count the new copy.

## 4. Apply + re-run (behind `--vg-realign`, family-local)

In `detect_and_assign`, when `cfg.vg_realign`, after the initial assignment: (a) collect `Reassign`/`NovelCopy`
decisions from the realign pass; (b) patch the corrected reads' `read_psv_obs`; (c) admit novel copies and
extend the family's copy set; (d) **re-run the family-local assignment** (the same `assign_family_detailed` path,
bounded to this family — NOT a global re-run) so `FamilyAssignment` (`assignments`, `copy_psv_alleles`,
`copy_abundance`, `realign_records`) and the downstream TSVs reflect (b)+(c). Flag off ⟹ none of (a)-(d) runs ⟹
byte-identical.

## 5. Validation (the end-to-end sim)

`bench/sim_vg_realign.py` (extends `bench/sim_genome.py`): plant a family with (i) a copy divergent enough
(~6–10%) that some of its reads mis-map to a sibling copy's locus, and (ii) a genome-absent divergent copy
(present in reads, removed from the reference). Run `copy_assign --vg-realign --lambda-global <λ>`; assert:
- **Correction:** the mis-mapped reads' `.assignments.tsv` copy = the correct copy (vs the wrong copy without
  `--vg-realign`).
- **Admission:** the genome-absent copy appears as an admitted copy; `chi_H`/`famcn_readonly` count it (higher
  than without `--vg-realign`).
Emit results into an updated `bench/VG_REALIGN.md`. Honest real-GGO caveat retained (admission may be
data-limited like O4).

## 6. Reuse / non-goals

Reuses the new `align_traceback`/`path_obs_at`, `pool_novel`, `absent_copy::admit_candidate` + its remap, the
`assign_family_detailed` re-run, the EM/`chi_H`/`famcn_readonly`. **Non-goals:** genome-wide VG-align (frontier);
unmapped-read routing (stays a separate opt-in); changing default (flag-off) behavior.

## Files

- **Modify** `src/rustle/vg_family/vg_realign.rs`: `align_traceback`, `path_obs_at`, a `apply_realign`
  orchestrator returning `(corrected_read_psv_obs: HashMap<usize, Vec<Option<u8>>>, admitted_copies:
  Vec<DenovoTranscript>)` from the family's candidates/copies/genome.
- **Modify** `src/rustle/vg_family/denovo_pipeline.rs`: in `detect_and_assign` under `cfg.vg_realign`, apply
  corrections + admissions and re-run `assign_family_detailed`.
- **Create** `bench/sim_vg_realign.py`; **update** `bench/VG_REALIGN.md` (move corrections/admissions from
  "follow-up" to "shipped"; add sim results).
- **Test:** in-crate `#[cfg(test)]` for `align_traceback` (path mismatch-count == `hw_distance`), `path_obs_at`
  (bases at PSV columns; gaps → `None`), and `apply_realign` (a mis-mapped read → corrected `read_psv_obs` that
  now matches the right copy; a novel pool → an admitted copy via a mock remap). Integration/sim for end-to-end.

## Testing (TDD)

Unit: `align_traceback` path mismatch+indel count equals `hw_distance` on random + adversarial pairs (mirror
the existing `hw_distance` test corpus); `path_obs_at` returns the read's base at each PSV position and `None`
on gaps; `apply_realign` corrects a planted mis-mapped read's `read_psv_obs` and admits a planted novel pool
(mock remap → not-remapping). Additivity: `cfg.vg_realign` false ⟹ `apply_realign` not called ⟹ byte-identical.

## Reproduce

```
copy_assign --bam GGO_mm.bam --fasta GGO.fasta --regions <fam> --vg-realign --lambda-global <λ> --out ca
python bench/sim_vg_realign.py   # planted mis-map + genome-absent copy -> corrected assignments + admitted+counted copy
```
