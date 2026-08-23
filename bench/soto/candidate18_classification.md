# SD classification of the 18 "no-Soto-homology" RNA candidates

**Question (meeting / FP concern):** the DNA-ceiling/RNA-subset audit flagged 37 "extra" RNA
candidate copies (RNA-detected loci not in Soto's DNA annotation). 19 were confirmed homologous to a
Soto SD cluster; **18 were labelled "no Soto homology."** Are those 18 false positives — phantom loci
RNA invented?

**Answer at the locus level: no. Zero of the 18 are phantom artifacts. All 18 are real genomic
loci.** This does not by itself establish gene-family membership: Soto groups segmental duplications,
and an SD cluster may contain several overlapping gene families or no intact gene.

---

## The earlier "18 no-homology" was an alignment artifact (caught on re-verification)

The verdict that produced "18 no-Soto-homology" recorded `best_id = 0.000` for **every** one of the
18 — its alignment step returned *no hit at all*, even for candidates that clearly align at 0.945
identity / 99% coverage to their own nearest family (e.g. chr1:16478340 → MST1L / ID_393). A
best_id of exactly 0.000 against 362 members is the signature of a failed/over-strict alignment
pass, not of true absence of homology.

Re-aligning the 18 directly to all 362 Soto member genomic sequences with the documented SD-edge
preset (`minimap2 -x asm20 -N5 -p0.5`, criterion id ≥ 0.80 ∧ cov-of-shorter ≥ 0.50) recovers the
homology. Spot-check: chr9:77896718 → ANKRD20A1 / ID_280 at 0.958 id, full coverage — an
ANKRD20A-containing SD. Cross-checked against a whole-genome self-alignment (paralog counts).

One legacy row also used an invalid coverage value greater than one. Re-alignment of
`chr1:148601552-148620972` with the current single-record axis rule gives identity 0.992 and coverage
1.000 to the ID_400 block (NBPF26/NOTCH2NLB), not 0.700/1.320. The TSV contains the corrected values.

## SD classification (Soto-member homology + genome self-alignment paralogs)

| bucket | n | what it is | FP? |
|---|---|---|---|
| **Full Soto SD-cluster match** | **13** | id ≥ 0.80, cov ≥ 0.50 to a Soto member — real homologous sequence Soto's per-member check missed | no — real locus; gene type unresolved here |
| **Partial / divergent Soto SD match** | **3** | id 0.83–0.93 but partial coverage (0.28–0.39) | no — real locus; cluster membership flagged only |
| **Non-Soto duplicated locus** | **1** | ≥ 2 high-identity genomic paralogs, below the Soto SD threshold | no — real duplicated locus |
| **Non-Soto single-copy locus** | **1** | one genomic copy, real gene, not a duplication (chr7:76085482) | no — real gene |
| **Phantom locus** | **0** | — | — |

**16 of the 18 match known Soto SD clusters** (13 full-length + 3 partial/divergent); the earlier
"no-homology" label was the aligner failing, not RNA inventing sequence. One more is a duplicated
non-Soto locus and one is single-copy. None is a phantom locus.

Notably, most are homologous to a **different** Soto SD cluster than the one they were flagged
*near* (proximity ≠ homology): candidates near ID_215 match ID_400/ID_212; near ID_113 → ID_116;
near ID_104/ID_260 → ID_280 (ANKRD20A). Flagging by genomic proximity mislabels the family; direct
homology assigns it correctly.

## What this does to the "37 extra candidates" precision story

- Before: 19 homologous + **18 unexplained**.
- After the corrected re-check: **35 of 37** extra candidates match a Soto SD cluster (19 + 16),
  one is a duplicated non-Soto locus, one is single-copy, and **0** are phantom loci.
- This is a locus/SD result, not a gene-family precision result. The typed audit in
  `bench/o1_gene_family_audit/` applies the additional gene-family gate.

## Gene-family typing

The typed audit requires SD homology, at least two genomic loci, and independent same-family gene
annotation. It supports four unique Soto-missing gene copies in this set: **NBPF8, NBPF19, GOLGA6L1,
and the putative ANKRD20A2 locus LOC128966611**. The HERC2-like `LOC102723564` record is retained as a
noncoding family locus, not counted as a gene copy. Three unannotated loci remain novel-copy
candidates rather than confirmed copies. The CHRFAM7A interval is explicitly rejected as an ULK4P
gene-family member even though it matches the ID_179 SD block.

## Per-candidate table

Full data: `bench/soto/candidate18_classification.tsv` (locus, near_family, genomic-paralog loci,
best Soto-member id/cov/family, class).

## Honesty rails

- "Soto SD-cluster match" means the candidate aligns to a Soto member at the edge threshold; it does
  **not** assert gene-family membership or read-level copy assignment (the latter is the separate K=0
  question).
- The 3 "partial/divergent" members clear identity but not the 0.50 coverage bar — real loci, but
  their family membership is partial-homology, flagged not asserted.
- The single-copy locus (chr7:76085482) is a real gene, not a duplication; calling it a "copy" would
  be the only mislabel, so it is bucketed separately.

**Reproduce:** `minimap2 -x asm20 -N5 -p0.5 soto_members.fa cand18.fa` (homology) +
`minimap2 -x asm20 -N8 -p0.6 chm13v2.0.fa cand18.fa` (genomic paralogs); classify per the buckets above.
