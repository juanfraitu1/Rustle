# Design: Certified copy-support guard (fabrication abstention)

Status: DESIGN (approved 2026-06-01; phantom handling = suppress-clear-phantoms + always emit certificate).
Next: implementation plan.

Grounding: this session's findings — the DAZ3 "recovery" is a FALSE POSITIVE (its reads fit DAZ1 ~15× better,
NM 6 vs 88; DAZ3 has ~0 genuinely-originating reads); the phantom persists under BOTH the HMM and the
phasing-EM because the cross-strand DAZ1↔DAZ3 decision is never phased (partition_and_remap splits by strand)
and DAZ3 is ASSEMBLED from preserved cross-mapped secondary reads. This is the gap-map item #4 (no
assembly-time fabrication guard). Memory: `project_vg_wiring`.

## 1. Problem & insight
StringTie/the current `--vg` pipeline will emit a gene copy assembled entirely from reads that actually
belong to a *sibling* copy (a chimeric/phantom copy). The boundary theorem says: a copy is identifiable only
if some reads distinguish it from its siblings. **Operationalized: a copy is real only if some reads fit IT
at least as well as any sibling.** DAZ3 fails this completely. The guard computes that fraction, abstains on
copies that fail (suppress), and emits the fraction as a calibrated confidence certificate on every copy.

## 2. Goal & non-goals
- **Goal:** before emitting a VG copy's transcripts, compute its **independent support** (fraction of its
  reads that fit it at least as well as any sibling) and (a) **suppress** copies below a conservative
  threshold (the phantom case — DAZ3), (b) **always emit** `copy_independent_support` as a GTF attribute.
- **Decision rule / invariants:** never suppress a genuinely-supported copy; default de-novo (non-`--vg`)
  output is byte-identical; the guard is conservative (only clear phantoms dropped).
- **Non-goals:** cross-strand joint phasing (option A — separate, for genuinely-expressed inverted copies);
  changing the EM; reference-absent discovery.

## 3. The metric
For each family F and each copy C (= a bundle assembled with `vg_copy_id = C`):
- **C-unique reads**: reads in C's bundle that do not multi-map to any sibling in F. They can only fit C →
  they always support C (a copy's anchor of independence).
- **C-multimappers**: reads in `F.multimap_reads` with a placement at C. For read r, let
  `rate_C = NM_C(r)/aligned_len(r)` and `rate_min_sib = min over siblings C' of NM_{C'}(r)/aligned_len`.
  r **supports C** iff `rate_C ≤ rate_min_sib + MARGIN`; else it **belongs elsewhere**.
  (NM and aligned length are on `BundleRead`; sibling NM comes from the read's placements in F's other
  bundles — already enumerated by family discovery. Works for cross-strand siblings, e.g. DAZ1/DAZ3, because
  the multimap index links them.)
- **independent_support(C)** = (#C-unique + #C-multimappers-that-support-C) / (#C-unique + #C-multimappers).

Behavior:
- **Suppress** copy C (drop the transcripts tagged `vg_copy_id = C`) iff `independent_support(C) < TAU`.
- **Emit** `copy_independent_support "<value>"` on every kept VG transcript (per-copy value).
- Edge cases: a copy with **no multimappers and ≥1 unique read** → support = 1.0 (kept). A copy with **0
  reads** → not assembled anyway. Non-VG transcripts (no `vg_copy_id`) → untouched, no attribute.

Defaults (tunable via env, calibrated in the plan on DAZ/fam175):
- `MARGIN` = 0.01 (NM-rate; ~1% — a sibling must fit *clearly* better to count as belongs-elsewhere) —
  `RUSTLE_VG_SUPPORT_NM_MARGIN`.
- `TAU` = 0.10 (suppress only copies with <10% independently-supporting reads — conservative) —
  `RUSTLE_VG_MIN_INDEP_SUPPORT`.

**Why this respects the identifiability floor:** NM-tie reads (genuinely ambiguous between near-identical
co-expressed copies) satisfy `rate_C ≤ rate_min_sib + MARGIN` for *both* copies, so they count as supporting
*each* → neither copy is falsely suppressed. We only suppress when reads fit a sibling *decisively* better
(the DAZ3 regime), or when a copy has no anchoring unique reads.

## 4. Where in the pipeline
A post-assembly pass inside the `--vg` flow (after transcripts are tagged with `vg_family_id`/`vg_copy_id`,
before final GTF write). Inputs already in hand: the per-family `multimap_reads` index (read → [(copy, ri)]),
`bundles[*].reads[ri].{nm, query_length/aligned span}`, and the assembled transcripts' `vg_copy_id`.
New function (e.g. `vg.rs::compute_copy_independent_support(family, bundles) -> HashMap<copy, f64>`), called
from the pipeline; transcripts of suppressed copies are filtered, survivors annotated. No new alignment, no
EM change.

## 5. Validation gate (every item; encoded in the plan + a new oracle check)
| Check | Target |
|---|---|
| DAZ default `--vg` | DAZ3 `independent_support ≈ 0` → **suppressed**; DAZ1 kept (support ≈ 1) |
| fam 175 | both copies `independent_support` high → **both kept**; correct B>A ratio |
| fam 214 | both copies kept |
| synthetic fixture | Obj-4 still 100%; both copies kept |
| no-false-suppression | on a panel of genuinely co-expressed copies, **zero** real copies dropped |
| default de-novo GGO_19 | **unchanged** (no `vg_copy_id` → guard is a no-op) |
| certificate emitted | every VG transcript carries `copy_independent_support` |

A new oracle check (`score_copy_support` in `bench/multi_copy_eval/`) asserts DAZ3 suppressed + fam175/214
kept, wired into `expectations.json`.

## 6. Risks
- **False suppression of a real lowly-expressed copy** (highest risk): mitigated by (a) the conservative TAU
  (0.10) + MARGIN, (b) unique reads always counting as support, (c) the no-false-suppression validation on a
  co-expressed panel. If a real copy is ever dropped, lower TAU / raise the bar — tunable.
- **Aligned-length/NM availability**: `nm` may be absent on some records → treat missing NM as "no sibling
  evidence" (read supports C; never penalize a copy for missing data).
- **Isolation**: guard runs only on transcripts with `vg_copy_id`; the default de-novo path produces none, so
  output is unchanged (verified in the plan).

## 7. What "done" looks like
`rustle --vg --genome-fasta` on DAZ no longer emits DAZ3 (it's suppressed, support ≈ 0), and every VG
transcript carries a calibrated `copy_independent_support`; fam 175/214/synthetic copies are all retained;
default de-novo is byte-identical. The tool now *abstains from fabricating* a copy whose reads belong to a
sibling — the certified-abstention contribution, demonstrated on the showcase locus.
