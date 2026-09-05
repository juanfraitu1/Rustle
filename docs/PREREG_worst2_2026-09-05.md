# PREREG — adjudication of the two worst sweep-v2 families (2026-09-05)

Families: fam_MCL32_073244 (ZSCAN5A/5C ×3 units; 27 unit reads, 24 ambiguous) and fam_MCL38_073244
(2 lncRNA + ZNF875 + ZNF569-like; MAPQ-60 placement agreement 119/310 — 190 reads placed by the aligner at
the ZNF569-like locus with NM 6/3,955 are assigned by O2 to ZNF875 at p 1e-133..1e-177).

Diagnosis before testing (from the assignments, the BAM records and the unit table):
- Both families' units are FRAGMENTS of the transcripts the reads carry. ZNF569-like: the GFF gene spans
  48660619-48689258 and its 3.3-kb 3' exon (48660619-48663874) is annotated, but core refinement trimmed the
  locus to the SD hull 48675744-48689016 and the read chain is bounded by the hull ⟹ 809-bp unit = the 5'
  fragment; the reads' 3.5 kb lie outside every unit. ZSCAN5A: the reads splice 86 kb into exons at
  68.42 Mb (annotated, same gene) that the hull-bounded chain excludes.
- Every ambiguous ZSCAN5A molecule has 6 BAM records (1 primary MAPQ 60 + 5 secondaries at the paralogues).
  O2 scores every record; the molecule contract demotes a molecule whose Assigned records name ≥2 copies.
- For ZNF569-like the molecule representative = the record with max n_decisive; the primary record overlaps
  its own 809-bp unit by ~150 bp, a secondary record over the ZNF875 zone covers that unit fully ⟹ the
  secondary is the representative ⟹ "assigned to ZNF875".

Predictions (fail rule: a failed prediction means the mechanism above is wrong; report it, do not patch):
- P1 (primary-only BAM, MCL32): ambiguous 24 → ≤4; assigned ≥20, every assigned read on unit 2 (ZSCAN5A);
  MAPQ-60 agreement 100%.
- P2 (primary-only BAM, MCL38): the 190 disagreeing reads are NOT assigned to ZNF875 any more (disagreement
  ≤10); they become tied/ambiguous, i.e. assigned drops by ≥150. Agreement among assigned ≥ 95%.
- P3 (units built with --no-core-refine, secondaries kept, MCL38): the ZNF569-like unit contains
  48660619-48663874; MAPQ-60 agreement ≥ 90% with the ORIGINAL BAM.
Decision rule: P2 ∧ P3 ⟹ the defect is the truncated node, and the fix belongs in O1 (unit extent must
cover the read-supported chain, not be clipped to the hull) — NOT a MAPQ/secondary filter in O2.
P1 ∧ ¬P3 ⟹ the record-reduction rule is the defect. Nothing is shipped from this adjudication.

## Addendum (after E1/E2, before T1–T3) — the read-support PSV filter
Found: `read_supported_columns` (PSV_MIN_ALLELE_READS 2, PSV_MIN_JUDGE_COV 4) keeps a copy-vs-copy column only
if ≥2 alleles each reach ≥2 reads. In a family with ONE expressed copy every true PSV is monomorphic in the
reads and is dropped; the survivors are the expressed copy's own polymorphic sites. MCL32: 4 columns survive,
copies 0/1 identical at all 4, the 11 wrong reads carry TGCC (2 columns = the ZSCAN5C allele) ⟹ p 6.7e-9.
Predictions with RUSTLE_PSV_READFILTER=0:
- T1 MCL32 (units v2, primary-only): ≥100 columns; 0 reads assigned to copy 0/1; agreement 100% of assigned.
- T2 MCL38 (units v2, original BAM): the 190 disagreements persist (mechanism = truncated unit, D1), ±20.
- T3 sweep v2 (35 families, original BAM): MAPQ-60 agreement does not fall below 95.3% and the number of
  assigned reads does not fall by more than 10%. Fail ⟹ the filter guards something real; keep it, report.
