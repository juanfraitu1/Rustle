# VG re-align supplement — Implementation Plan

> REQUIRED SUB-SKILL: superpowers:subagent-driven-development.

**Goal:** Re-align poor-fit/unmapped candidate reads to O1's per-family copy-paths (the `FamilyGraph`'s paths);
accept a re-alignment only when it is significantly better than the linear one (`min_p` certificate) — correcting
mis-assignment and admitting reference-absent copies (via O4). Opt-in `--vg-realign`, default byte-identical.

**Architecture:** The family's copy-paths (`recover_paralog_path` → copy consensus sequences) ARE the VG paths.
A candidate read is aligned to each copy-path (reusing the existing edlib/`aln_id` machinery), the best path's
PSV-space observations are scored with `read_copy_evidence`, and a re-assignment is accepted iff significant AND
better than the read's linear-locus evidence. Reads that thread a consistent NOVEL path are pooled and admitted
as reference-absent copies via `absent_copy::admit_candidate` (O4 discipline). Corrected assignments + admitted
copies feed the existing EM / `chi_H` / `famcn_readonly`. (POA-graph shared-exon alignment is an optimization,
not required — aligning to the copy-paths achieves the reference-bias fix; note it as a follow-up.)

## Global Constraints
- Additive & non-destructive: `--vg-realign` off ⟹ no `.vg_realign.tsv`, all other outputs byte-identical.
- No arbitrary thresholds for ADMISSION: acceptance uses the `min_p`/Poisson-binomial certificate
  (`read_copy_evidence`) and O4's `admit_candidate`. Candidate-SELECTION thresholds are cheap pre-screens only.
- Reuse: `recover_paralog_path`, `read_copy_evidence`/`min_p`, `absent_copy::admit_candidate`, the EM/copy-number.

## File Structure
- `src/rustle/vg_family/vg_realign.rs` — CREATE: candidate selection, re-align-to-copy-paths, accept, pooling.
- `src/rustle/vg_family/denovo_assemble.rs` — MODIFY: add an unmapped-read fetch (separable stage T2).
- assign path (`copy_assign_pipeline.rs`/`denovo_pipeline.rs`) — MODIFY: invoke behind `vg_realign` flag.
- `src/bin/copy_assign.rs` — MODIFY: `--vg-realign` + emit.
- `bench/sim_vg_realign.py`, `bench/VG_REALIGN.md` — CREATE.

---

### Task 1: Candidate selection (mis-mapped)

**Files:** Create `src/rustle/vg_family/vg_realign.rs` (register in `mod.rs`). Test in-crate.

**Interface:** `pub struct RealignParams { pub max_mapq: u8, pub min_clip_frac: f64, pub min_div: f64, pub min_reads: usize }` (defaults e.g. `max_mapq: 20, min_clip_frac: 0.20, min_div: 0.05, min_reads: 3`).
`pub fn is_candidate(read: &BamRead, div: f64, clip_frac: f64, p: &RealignParams) -> bool` = `read.mapq <= p.max_mapq || clip_frac >= p.min_clip_frac || div >= p.min_div`. (Divergence + clip fraction are computed by the caller from the CIGAR/NM — provide small helpers `clip_fraction(cigar)`/`read_divergence` or reuse existing ones in `denovo_assemble`/`copy_split`; confirm which exist.)

- [ ] Test `is_candidate_flags_poor_fit`: a low-MAPQ read, a high-clip read, and a high-div read each return true; a MAPQ-60 / no-clip / low-div read returns false. Write RED, implement, GREEN. Commit.

### Task 2: Unmapped-read fetch (separable stage)

**Files:** Modify `src/rustle/vg_family/denovo_assemble.rs`. Test in-crate.

`aligned_reads_from_bam` skips unmapped records (`flags.is_unmapped()`). Add `pub fn unmapped_reads_from_bam(bam_path, threads) -> Result<Vec<(String, Vec<u8>)>>` returning `(read_name, seq)` for `is_unmapped()` primary records (skip secondary/supplementary). This is the input to the minimizer pre-screen (Task 3). Keep it independent — if this task is deferred, Task 3's pre-screen simply has no unmapped input and the mis-mapped path still works.

- [ ] Test on a tiny fixture BAM with ≥1 unmapped record (or construct a noodles record) that it returns the unmapped read's name+seq and excludes mapped ones. Commit.

### Task 3: Re-align to copy-paths + minimizer routing

**Files:** `src/rustle/vg_family/vg_realign.rs`. Test in-crate.

**Interfaces:**
- `pub fn route_unmapped(seq: &[u8], family_consensuses: &[(usize, Vec<u8>)], min_shared: usize) -> Vec<usize>` — families whose consensus shares ≥ `min_shared` canonical minimizers (`minimizers`, MINIMIZER_K/W) with `seq` (the pre-screen for unmapped reads; mis-mapped reads skip this — their family is their locus).
- `pub fn realign_to_paths(read_seq: &[u8], copy_seqs: &[Vec<u8>], psv_positions_in_path_space) -> Option<RealignHit>` where `RealignHit { best_copy: usize, path_obs: Vec<Option<u8>>, aln_id: f64 }` — align the read to each copy-path consensus (reuse the existing edlib/`aln_id` alignment used elsewhere — confirm the helper), pick the best `aln_id`, and extract the read's observed bases at the family PSV columns along that alignment (so the significance gate scores it in the same PSV/allele space as the linear evidence). Return `None` if no path aligns above a minimal `aln_id` floor.

- [ ] Test `realign_routes_and_scores`: a read identical to copy 1's consensus → `best_copy == 1`, high `aln_id`, `path_obs` matches copy 1's alleles; `route_unmapped` returns the family whose consensus shares the read's minimizers. Commit.

### Task 4: Significance-gated accept + novel-copy pooling

**Files:** `src/rustle/vg_family/vg_realign.rs`. Test in-crate.

**Interfaces:**
- `pub enum RealignAction { Reassign(usize), NovelCopy, Reject }`.
- `pub fn accept_realignment(hit: &RealignHit, copies: &[CopyProfile], linear_copy: Option<usize>, params: &AssignParams) -> RealignAction`: build `ReadFeatures{ psv_obs: hit.path_obs, .. }`, call `read_copy_evidence(&rf, copies, params, &[])`. If the best copy's evidence is significant (`min_p < params.alpha/((K-1).max(1))`) AND (no `linear_copy`, or the VG path's evidence beats the linear copy's evidence at the read's decisive columns) → `Reassign(hit.best_copy)` if `best_copy` is an existing copy, else `NovelCopy`. Else `Reject`.
- `pub fn admit_novel_copies(novel_reads_by_path, host: &DenovoTranscript, remap: F, params: &AbsentCopyParams) -> Vec<DenovoTranscript>`: pool `NovelCopy` reads by their threaded path; for each pool with `≥ min_reads`, build a `CollapsedCandidate` and call `absent_copy::admit_candidate` (reuse O4's gates: min_p-distinct, strand-symmetric, non-remap ≥98%). Return admitted synthetic copies.

- [ ] Tests: `accept_significant_path_reassigns` (significant + beats linear → Reassign); `accept_insignificant_rejects` (tied path → Reject); `novel_path_below_min_reads_not_admitted`; `novel_path_admitted_with_o4_gate` (≥min_reads + strand-sym + non-remap via a mock remap closure → admitted). Commit.

### Task 5: Pipeline wiring + emit

**Files:** Modify the assign path + `src/bin/copy_assign.rs`. Test additivity.

Behind `--vg-realign` (+ `--vg-realign-max-mapq`, `--vg-realign-min-clip`, `--vg-realign-min-div`, `--vg-realign-min-reads`): per family, collect candidates (Task 1) [+ unmapped via Task 2/3 routing], re-align (Task 3), accept (Task 4); route `Reassign` into the family's read→copy assignment (so the EM/`chi_H` see the corrected placement) and admitted `NovelCopy` copies into the copy set (like O4). Emit `<out>.vg_realign.tsv` (`read_name  family_id  action  target_copy  path  min_p  linear_locus`). Off ⟹ nothing runs, outputs byte-identical.

- [ ] Test: `--vg-realign` off ⟹ `.assignments.tsv`/`.families.tsv`/`.famcn_readonly.tsv` byte-identical + no `.vg_realign.tsv`; on ⟹ the file appears. Confirm on a fixture or by gated-code inspection. Commit.

### Task 6: Sim validation + doc

**Files:** Create `bench/sim_vg_realign.py`, `bench/VG_REALIGN.md`.

`sim_vg_realign.py` (extends `bench/sim_genome.py`): plant a family with one copy divergent enough (~6–10%) that its reads mis-map to a sibling locus or fail to map on the linear reference. Run `copy_assign --vg-realign`; assert (a) the mis-mapped reads are re-assigned to the correct copy-path; (b) the divergent copy is admitted as a novel copy and counted by `chi_H`/`famcn_readonly`. Write `bench/VG_REALIGN.md` with the sim result + the honest real-GGO caveat (admission may be data-limited like O4; the correction leg should still move reads). Commit.

## Self-Review
- Spec coverage: candidates (T1 mis-mapped + T2/T3 unmapped), re-align (T3), significance accept + O4 admission (T4), wiring/emit (T5), validation (T6) — all covered.
- Type consistency: `RealignParams`/`RealignHit`/`RealignAction` produced T1/T3/T4, consumed T5; reuses `BamRead`, `CopyProfile`, `ReadFeatures`, `read_copy_evidence`, `AssignParams`, `absent_copy::admit_candidate`, `DenovoTranscript`.
- Placeholder scan: alignment helper (`aln_id`/edlib) and divergence/clip helpers are marked "confirm which exists" — the implementer must locate the existing one and reuse it (DRY), not invent a new aligner.
