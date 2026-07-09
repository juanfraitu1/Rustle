# Reference-free copy number in the pipeline + O1↔O2 harmony — Plan

> REQUIRED SUB-SKILL: superpowers:subagent-driven-development.

**Goal:** Emit a genome-free copy-number estimate per family in the shipped `copy_assign` pipeline —
`chi_H` (distinguishable copies = O1 conflict count) + `depth_cn` (recovers identical/collapsed copies from
read depth, incl. unassignable reads) + `regime` + `famCN_readonly = max(chi_H, depth_cn)` — and pin the
O1↔O2 coupling with a test. Validated: `depth_cn` = 66% within 25% of `asm_hapCN`, recovering Tier-3 identical
copies no per-read assignment can (advisor's point); honest caveat = it counts EXPRESSED copies.

**Architecture:** The reads-only copy number is powered by the SAME family object O1 builds (copies → PSV
bubbles → EM). `chi_H` = # distinct pairwise-conflicting copy hap-vectors (from `copy_psv_alleles`). `depth_cn`
= `E_fam / λ_global`, where `λ_global` = median `n_reads` over single-copy (non-family) de-novo transcripts —
a reference-free RNA expression floor (Sudmant 2010 parCN, ported). The depth leg sums ALL family reads incl.
the EM's SoftZone/unassignable mass — that's why unassignable reads indicate copy number.

## Global Constraints
- Additive & non-destructive: new `<out>.famcn_readonly.tsv` only; existing outputs byte-identical.
- Reference-free: no genome projection / no assembly used — only RNA reads + the RNA single-copy expression floor.
- Match the reference logic: `chi_H` per `bench/family_copy_number.py` (copyonly conflict colors); `depth_cn`/
  `λ_global` per `bench/rna_copy_number_depth.py` (`E_fam/λ_global`, median over single-copy transcripts).

---

### Task R1: Reads-only copy number emit in the binary

**Files:** Modify `src/bin/copy_assign.rs` (+ a helper in `src/rustle/vg_family/em_copy_assign.rs` or a new
`readonly_copy_number.rs` if cleaner). Test: in-crate `#[cfg(test)]`.

**Interfaces / logic:**
- `pub fn chi_h(copy_alleles: &[Vec<Option<u8>>]) -> usize` = number of colors of the conflict graph H over
  copies, where copies i,j CONFLICT iff they differ at ≥1 column where both are `Some` (distinguishable), and
  non-conflicting copies (identical hap-vector on shared columns) share a color. Implement as: greedy grouping —
  two copies in the same group iff they do NOT conflict; `chi_h` = number of groups. (Matches
  `family_copy_number.py` copyonly_K = # distinct hap-vectors.)
- `λ_global` (computed once per run): median `n_reads` over de-novo transcripts that are SINGLE-COPY (not in any
  multi-copy family) and multi-exon. The binary already has all detected families/transcripts; take the
  `n_reads` of single-copy ones, median them.
- Per family: `E_fam = fa.n_reads`; `depth_cn = E_fam / λ_global` (f64; `NaN`/blank if `λ_global<=0`); `chi_H =
  chi_h(&fa.copy_psv_alleles)`; `regime = if fa.n_copies as f64 >= chi_H { "reference_resolved" } else {
  "reference_collapsed" }` (or per family_copy_number.py's `n_loci < chi_H` rule — use `fa.n_copies` as n_loci);
  `famcn_readonly = max(chi_H as f64, depth_cn)` (fallback to `chi_H` if `depth_cn` is NaN).
- Emit `<out>.famcn_readonly.tsv`, header `family_id\tchrom\tn_copies\tchi_H\tdepth_cn\tregime\tfamcn_readonly`
  — one row per family. Always written (it needs no extra flag; it's cheap and genome-free), OR behind `--em`
  if the controller prefers to keep it opt-in — DECIDE: write it always (additive new file, no behavior change).

**Steps:** TDD `chi_h` first (identical copies → 1; K pairwise-distinct → K; a mixed case). Then wire λ_global +
the per-family loop + emit. Verify `<out>.famcn_readonly.tsv` appears with correct columns on a fixture and
existing outputs are byte-identical. Commit.

### Task R2: O1↔O2 harmony test

**Files:** in-crate `#[cfg(test)]` (near `em_assign_family` or the pipeline).

Pin that the EM (O2) and the reads-only copy number (O1) consume the SAME copy object, so improving O1 flows
into O2: a test that builds a small family's `copy_psv_alleles`, asserts `chi_h(&copy_alleles)` equals the
number of distinct copies the EM sums abundance over (`em_assign_family(...).abundances.len() == copy_alleles
.len()`), and that collapsing two copies to identical hap-vectors drops BOTH `chi_h` by 1 AND leaves the EM
unable to separate them (both SoftZone). This encodes "O1 defines K/copies; O2 assigns within it; a better O1
(more distinguishing bubbles) raises chi_h AND lets the EM certify more." Commit.

### Task R3: Validation doc

**Files:** Create `bench/READONLY_COPY_NUMBER.md`.

Write up the validation (already computed): per-gene aggregated, `depth_cn` = 66% within 25% of `asm_hapCN`
(corr 0.52), `chi_H` = 49%; the advisor's-point wins where depth recovers identical copies no assignment can
(LOC115930538 `chi_H=1 → depth 11.4 ≈ asm 12`; LOC109025447 `1 → 15.8`; LOC129526550 `5 → 11.8`); the honest
EXPRESSED-copy caveat (under-counts silent copies: 13/59 genes, e.g. LOC129531752 sees 2 of 22). State the
O1↔O2 harmony (one family object; chi_H = O1 conflict count = the EM's K; depth uses all reads incl. the EM
SoftZone mass). Reference `bench/family_copy_number.py`, `rna_copy_number_depth.py`, `em_consistency.md`. Commit.
