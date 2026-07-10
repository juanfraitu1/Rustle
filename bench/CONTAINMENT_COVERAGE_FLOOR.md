# The Containment "defect" is a coverage floor, not a prunable artifact

**Date:** 2026-07-10. Investigated after RFPL and r4 emitted `Containment`-flagged overlapping copies. A
7-agent workflow characterised the flagged copies, proposed three independent discriminators, and adversarially
tested each against the two tests that protect genuine adjacent paralogs. **Verdict: do not add a pruning rule.**

## What the flagged copies actually are

- **RFPL** (`NC_086018.1:30203643-30381055`, 104 primary reads over 177 kb). The raw overlapping transcripts
  `{30314371-30368310, 30320521-30368310}` **share the exact junction `(30366265, 30368092)`**, so
  `prune_same_locus` clause (a) already collapses them — they are 5′-truncated isoforms of one sparsely
  expressed, unannotated minus-strand unit (100% identical over their overlap). The residual `Containment` flag
  is a *pooling* artifact: the 707-read pooled consensus sprawls across the region and staggers against a
  28-read fragment. RFPL2 itself forms no family.
- **r4** (`NC_073228.2:182473722-182663103`, gene-dense, ~50–200 reads). The flagged pair is a **convergent
  opposite-strand** overlap: a `+` lncRNA-derived transcript whose exons interleave into the *introns* of the
  `−` FAM153B gene, overlapping only 8127 bp of intronic span — landing **585 bp below** the antisense clause's
  bar (8127 < 0.5·17424 = 8712). A low-coverage chimeric fragment, not a paralog. The r4 family core
  (FAM153B− / FAM153B+) is a *legitimate* annotated paralog pair; the extra "copies" are nested same-gene
  fragments plus this one chimera.

## Why no pruning rule can remove them

Every candidate discriminator was defeated by the same counterexample, **CAFAM0** — and it defeats them by
sitting in the *identical feature cell* as the genuine adjacent tandem paralogs the test suite pins:

| feature | protected tests T1/T2 | CAFAM0 (must-drop) | verdict |
|---|---|---|---|
| reciprocal overlap (ov / shorter) | **0.40** | 0.27 | CAFAM0 overlaps *less* — a containment threshold cannot cut between them |
| minority-copy reads | 9 | **28** | a "drop if reads < K" gate needs K ≤ 9 to keep the tests and K ≥ 29 to drop CAFAM0 — contradiction |
| strand | same | same | antisense clause does not apply |
| POA homology | low (random seq) | passed the bar | homology clause does not fire |

The full drop-target set straddles the protected pair on containment fraction:
`0.168 < 0.273 < (T1/T2 = 0.40) < 0.466`. There is no coordinate, read-count, or homology property in
`prune_same_locus` or `colocated_families` that removes the artifacts without also deleting real overlapping
paralogs. The one "principled" separator — overlap-window sequence identity — is degenerate on real data:
reference-aligned transcripts at shared coordinates are ~100% identical *by construction*, so that rule reduces
to "collapse any pair overlapping by ≥ 1 bp," which would destroy genuine tandem and inverted-duplication
paralogs. (The tests survive it only because they fabricate divergent random sequence at shared coordinates —
biologically impossible for aligned reads.)

**This is a region-level coverage floor** (104 reads / 177 kb), invisible to the pruner by design. The clean
families are all span-**disjoint** — DAZ (~40 kb gap), GSTM, MAGEA, RBMY — so they never enter any overlap
clause and are untouched whatever the region coverage.

## What was done instead

Nothing in the pruner. The `catalog_overlaps` warning already flags all three pairs; it was only **describing
them wrong** — it told the `DuplicateLocus` "min_p == 1 masquerade" story for `Containment` pairs too. The
warning is now kind-accurate:

- **DuplicateLocus** (recip ≈ 1): one locus admitted twice; reads abstain at `min_p == 1` (not the K=0 wall).
- **Containment** (recip ≪ 1): a fragment/readthrough nested in a real copy, **inflating the copy count** on
  low-coverage regions — reported, not pruned, because it shares its feature cell with real overlapping
  paralogs.

## The honest risk statement

The residual `Containment` flags remain in the output, and a reader could mistake them for real extra copies —
the kind-accurate warning is the mitigation. The **far larger** risk is the opposite path: any threshold tuned
to catch these (containment nudge to 0.466, read gate to K ≥ 29, overlap-homology) silently collapses genuine
overlapping paralogs and breaks the adjacent-paralog tests — trading a cosmetic low-coverage flag for real
biological false negatives, and introducing exactly the arbitrary threshold the advisor rejects.

The genuine upstream cause of CAFAM0 (a multi-exon readthrough bridge pooling a 707-read consensus over the
region) lives in the readthrough filter / E_r oracle, not in copy-set pruning. If RFPL is ever a priority, the
lever is more reads or a multi-exon-readthrough discriminator — not a containment threshold.

Related: `bench/YAG_VS_ISOCON.md`, `bench/DAZ2_RECOVERED.md`, `denovo_pipeline::catalog_overlaps`.
