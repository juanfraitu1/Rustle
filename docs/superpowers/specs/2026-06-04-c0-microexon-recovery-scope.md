# Scope: recover c0's micro-exon chain in tandem mode (the post-gate residual)

**Date:** 2026-06-04. **Status:** root-caused + fix proxy-validated; principled implementation scoped
(NOT built). Follows the over-enumeration gate work (`2026-06-04-over-enumeration-research.md`).

## Problem
After the primary-junction gate removes the phantom over-enum, the dominant residual on the RBMY
mega-bundle is c0 (LOC129530243): VG retains a **4-micro-exon block** (89/111/111/115 bp) as one exon
instead of splicing it, producing a class-`m` transcript that doesn't match c0's annotated chain.
Baseline splices it correctly. c0 IS genuinely expressed (4–6 full-length primary reads, chain matches
RefSeq XM_055377108 within ±1 bp).

## Root cause (grounded)
The tandem sub-bundle assembly diverges from baseline's bundle assembly for the SAME primary reads:
- It is **not** the `synthetic` flag (a post-hoc label, affects downstream filters, not
  `extract_transcripts`), nor sub-bundle boundary truncation (the c0 sub-bundle covers the gene).
- It **is read WEIGHT**: c0's primary reads are under-weighted in the tandem path. The existing tandem
  weight-floor (`RUSTLE_VG_TANDEM_WEIGHT_FLOOR`, raises primary reads to a floor) is **support-gated**
  (`RUSTLE_VG_TANDEM_SUPPORT_FLOOR=0.75`, only floors "certified" copies) and **excludes c0** — so
  c0's reads keep low EM-apportioned weight → low flow → `path_extract` prefers the retained-intron
  path. Baseline assembles c0's reads at weight 1.0 and splices the micro-exons.

## Fix — proxy-validated
Flooring c0's primary-read weight recovers it. Measured (`gate + EM-off + RUSTLE_VG_TANDEM_WEIGHT_FLOOR=0.5
RUSTLE_VG_TANDEM_SUPPORT_FLOOR=0.0`):

| | Tx Sn/Pr | Locus Sn/Pr | c0 class | exact matches |
|---|---|---|---|---|
| baseline | 65.0/81.2 | 83.3/100 | — | 13/20 |
| **this lever** | **70.0**/73.7 | 83.3/83.3 | **`=` (10 exons)** | **14/20** |

c0 becomes an **exact match**, and Tx Sn **70 > baseline 65** — the first config where VG *beats*
baseline on sensitivity, with c0 genuinely recovered (the +1 is c0).

## The honesty boundary — why `support-floor=0` is a PROXY, not the fix
Blanket `support-floor=0` floors **every** tandem copy, including **c6** (LOC129530242, 2 genuine
reads, held-out concordance 0.09 = non-identifiable). With the blanket floor, c6 also turns class-`=`
— a structural match assembled largely from **borrowed/secondary structure**, not its 2 own reads.
That is the DAZ3 fabrication pattern (confident structure on a non-identifiable copy). So the blanket
floor recovers c0 (genuine) AND fabricates c6 (borderline) — the c0 win is real but the lever is too
coarse.

## Principled implementation (the actual scope)
**Evidence-gate the weight floor by IDENTIFIABILITY, not by `independent_support`.** Floor a tandem
copy's primary reads only when the copy is genuinely identifiable — reuse the per-read `em_ev_decisive`
signal (fix #3) / the divergence the held-out PSV concordance measures: c0 (concordance 1.0,
em_ev_decisive high) gets floored → recovered; c6 (concordance ~0, non-identifiable) does NOT → no
fabrication. This gives the c0 win without the c6 borderline.

Work items:
1. **Compute per-copy identifiability** at the weight-floor pass (mean `em_ev_decisive` of the copy's
   reads, already computed for `copy_confidence`) and gate the floor on it (replace/augment the
   `independent_support >= support_floor` gate). Behind the existing `RUSTLE_VG_TANDEM` + a new
   sub-flag; default-off until validated.
2. **EM-off-for-tandem — DONE (item #2, the scoped anchor-prior).** Grounded: the EM's two roles
   are SEPARATE writes — `read.weight` apportionment (vg.rs:5630, → assembly flow) vs the attribution
   fields (5648-5656). And the attribution (`em_ev_decisive`) is the PRE-PRIOR evidence margin (fix
   #3), so it is BYTE-IDENTICAL with the legacy prior (verified RBMY: c0 1.000, c4 0.026 under EM-on
   AND EM-off). So the anchor-prior contributes NOTHING to attribution — it is purely the
   sink-collapse assembly harm. FIX: `anchor_prior_on` (vg.rs:5369) is now ` && !is_tandem_family` —
   tandem families use the legacy pileup prior (no sink-collapse → c0 assembles), dispersed families
   (DAZ) keep the anchor-prior (DAZ3 phantom stays suppressed). Result: RBMY vs-RefSeq 60→70 BY
   DEFAULT (no EM-off flags), c0 still `=`, attribution preserved; DAZ3 4.04, suite 222/0, headline
   95.6/90.5, obj5 1.0/abstain, tier1 100. Dual-reference verdict: VG >= baseline-real AND finds more.
3. **Clean up the c0 residual partial** (the extra class-`m` 2-exon tx, Pr 73.7 vs baseline 81.2) — a
   retained-intron partial alongside the exact match; likely the same micro-exon coverage issue.

## Risks
- **Fabrication** (the c6 pattern) — mitigated by the identifiability gate; MUST validate c6 stays
  non-`=` and that no non-identifiable copy gains a borrowed-structure match.
- **Over-enum** — flooring raises coverage; a too-high floor over-enumerates (wf=2.0 → c0 3 tx).
  Keep wf ≤ 1.0.
- **Flag sprawl** — the proxy uses 4 env knobs; the shipped version should fold the principled gate
  in so it's not knob-tuning.

## Validation plan (before any default change)
RBMY (c0 `=`, c6 NOT `=`, Sn ≥ baseline) · synthetic novelty preserved (hidden/inverted) · ≥2 seeds ·
DAZ3 byte-identical (scoped to TandemCopy) · default headline 95.6/90.5 · suite 222/0 · obj5/tier1.
Do NOT promote a borrowed-structure c6 as recovery.
