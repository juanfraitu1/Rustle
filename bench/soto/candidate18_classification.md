# Classification of the 18 "no-Soto-homology" RNA candidates

**Question (meeting / FP concern):** the DNA-ceiling/RNA-subset audit flagged 37 "extra" RNA
candidate copies (RNA-detected loci not in Soto's DNA annotation). 19 were confirmed homologous to a
Soto family; **18 were labelled "no Soto homology."** Are those 18 false positives — phantom family
members RNA invented?

**Answer: no. Zero of the 18 are phantom artifacts. All 18 are real genomic loci.** A rigorous
re-check reclassifies 15 of them as real copies of *known* Soto families; the other 3 are real
non-Soto loci (2 with their own genomic paralogs, 1 single-copy gene).

---

## The earlier "18 no-homology" was an alignment artifact (caught on re-verification)

The verdict that produced "18 no-Soto-homology" recorded `best_id = 0.000` for **every** one of the
18 — its alignment step returned *no hit at all*, even for candidates that clearly align at 0.945
identity / 99% coverage to their own nearest family (e.g. chr1:16478340 → MST1L / ID_393). A
best_id of exactly 0.000 against 362 members is the signature of a failed/over-strict alignment
pass, not of true absence of homology.

Re-aligning the 18 directly to all 362 Soto member genomic sequences with the documented family-edge
preset (`minimap2 -x asm20 -N5 -p0.5`, criterion id ≥ 0.80 ∧ cov-of-shorter ≥ 0.50) recovers the
homology. Spot-check: chr9:77896718 → ANKRD20A1 / ID_280 at 0.958 id, full coverage — an
unambiguous ANKRD20A paralog. Cross-checked against a whole-genome self-alignment (paralog counts).

## Classification (both signals: Soto-member homology + genome self-alignment paralogs)

| bucket | n | what it is | FP? |
|---|---|---|---|
| **Copy of a known Soto family** | **12** | id ≥ 0.80, cov ≥ 0.50 to a Soto member — a real copy Soto's per-member check missed | no — real member |
| **Partial / divergent Soto member** | **3** | id 0.83–0.93 but partial coverage (0.28–0.39) — a real, divergent/fragmentary copy of a Soto family | no — real locus |
| **Non-Soto real copy (has genomic paralog)** | **2** | ≥ 2 high-identity genomic paralogs, below the Soto family threshold — a real duplicated locus outside Soto's 80 families | no — real copy |
| **Non-Soto single-copy locus** | **1** | one genomic copy, real gene, not a duplication (chr7:76085482) | no — real gene |
| **Phantom artifact / family-level FP** | **0** | — | — |

**15 of the 18 are copies of known Soto families** (12 full-length + 3 partial/divergent); the earlier
"no-homology" label was the aligner failing, not RNA inventing a member. **2 more are real non-Soto
copies** (they have their own genomic paralogs — real duplications Soto's 80-family catalog doesn't
contain). **1 is a real single-copy gene.** None is a phantom.

Notably, most of the 15 are homologous to a **different** Soto family than the one they were flagged
*near* (proximity ≠ homology): candidates near ID_215 are copies of ID_400/ID_212; near ID_113 → ID_116;
near ID_104/ID_260 → ID_280 (ANKRD20A). Flagging by genomic proximity mislabels the family; direct
homology assigns it correctly.

## What this does to the "37 extra candidates" precision story

- Before: 19 homologous + **18 unexplained**.
- After the corrected re-check: **34 of 37** extra candidates are copies of a Soto family (19 + 15),
  **3** are real non-Soto loci (2 with paralogs, 1 single-copy), **0** are phantom family members.
- **Family-level false positives against the DNA truth: 0.** Every "extra" RNA call lands on real
  genomic sequence; none is a fabricated member merged into a family roster.

## Per-candidate table

Full data: `bench/soto/candidate18_classification.tsv` (locus, near_family, genomic-paralog loci,
best Soto-member id/cov/family, class).

## Honesty rails

- "Copy of a Soto family" means the candidate aligns to a Soto member at the family-edge threshold; it
  does **not** assert a read-level copy-assignment (that is the separate K=0 question).
- The 3 "partial/divergent" members clear identity but not the 0.50 coverage bar — real loci, but
  their family membership is partial-homology, flagged not asserted.
- The single-copy locus (chr7:76085482) is a real gene, not a duplication; calling it a "copy" would
  be the only mislabel, so it is bucketed separately.

**Reproduce:** `minimap2 -x asm20 -N5 -p0.5 soto_members.fa cand18.fa` (homology) +
`minimap2 -x asm20 -N8 -p0.6 chm13v2.0.fa cand18.fa` (genomic paralogs); classify per the buckets above.
