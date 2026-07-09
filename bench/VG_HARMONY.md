# VG harmony: finding copies not in the linear genome, tied to the EM and copy number

**Date:** 2026-07-08. The variation graph exists to escape **linear-reference bias**: reads from copies that
are divergent from — or absent in — the linear assembly get mis-mapped, clipped, or dropped, so those copies
are invisible to a linear pipeline. This note records how the EM copy-assignment (O2) and the reference-free
copy number are harmonized with the reference-absent (O4) and junction (O3) machinery so that copies **not in
the linear genome** are found, assigned, and counted — and where the honest frontier still is.

## Three ways a copy not in the linear reference is recovered

1. **COLLAPSED copies (present in the assembly but merged) → `depth_cn`.** The reference-free copy number's
   depth leg (`E_fam/λ_global`) recovers copies the linear assembly collapsed onto one locus, from read depth
   alone — no per-read assignment needed. Validated: LOC115930538 `chi_H=1 → depth 11.4 ≈ asm 12`
   (`bench/READONLY_COPY_NUMBER.md`). This is a direct linear-reference-bias win.

2. **DIVERGENT reference-absent copies → O4 admission → EM + `chi_H`.** `absent_copy::admit_candidate` admits a
   candidate as a synthetic copy iff it does **not** remap to its host locus at ≥98% identity (genuinely
   divergent, i.e. not in the reference). Admitted copies are pushed into the family copy set
   (`denovo_pipeline.rs`, `--absent-copies`) → `build_family_profiles` → `copy_psv_alleles`. **Pinned by test
   `absent_copy_is_assigned_and_counted`:** such a copy raises `chi_H` by a color (counted by the reference-free
   copy number) AND the EM assigns its reads `Certified` with argmax on the absent copy (abundance > 0). So a
   copy not in the reference, once admitted, is a first-class copy for both O2 and the copy number.

3. **Copies defined by a NOVEL JUNCTION (O3) → threaded into the EM.** A splice not in the reference can be the
   only thing distinguishing a copy (identical PSVs). The EM was PSV-only; O3 copy-specific junctions are now
   threaded into `em_assign_family` (via `read_junctions`/`copy_junctions` on `FamilyAssignment`, folded into
   the existing `read_copy_evidence` junction term). **Pinned by test `junction_only_copy_is_resolved_by_em`:**
   two copies with identical PSVs but distinct junctions are `SoftZone` PSV-only and become `Certified` once the
   junctions are threaded — the O2↔O3 harmony.

## O1 ↔ O2 ↔ O3 ↔ O4 — one object

All of these consume the same family object: O1 defines copies (incl. O4-admitted absent ones) and their PSV
bubbles + junctions; O2 (the EM) assigns reads to those copy-paths using PSVs **and** junctions; the
reference-free copy number counts them (`chi_H` for distinguishable, `depth_cn` for collapsed). Improving O1
(cleaner copies, more bubbles) or O4 (more admitted absent copies) or O3 (junction resolution) each raises the
copies the EM can certify and the copy number can count — they move together (pinned by
`o1_o2_share_one_copy_object` + the two tests above).

## Honest frontier (scoped as future work)

The one reference-bias source that remains: reads are still **aligned to the linear genome first** (minimap2 →
BAM). A copy whose reads never map to the family region at all is invisible *before* the VG threading step —
O4's collapsed/divergent recovery operates on reads that DID map near the host. Eliminating this fully requires
**genome-wide VG-alignment** (align reads to a variation graph with all copies as paths — vg-giraffe /
GraphAligner) so absent-copy reads align to their true path. That is a separate, larger build; the current
harmony captures collapsed copies (depth), divergent admitted copies (O4), and junction-defined copies (O3),
which is most of the reachable gain without a graph aligner.

Relates to `bench/READONLY_COPY_NUMBER.md`, `bench/em_consistency.md`, `absent_copy.rs`, `allele_specific_junctions.rs`.
