# PREREG — root-cause fix for MCL "fragmentation" (row 694): fold within clusters + exon-less span (2026-09-05)
Diagnosis (frag_diag.txt): the 9 "fragmented" Soto members and 19 of the 43 annotated misses are records that
overlap a DIFFERENT family's record on exon bases; fold-first + representative-only discards them. 20 Soto members
exist only as RefSeq pseudogene records without exon children (no exonic denominator, no edge possible).
Arms on the same slice/scorer as §6ev (min_size 2): (a) --fold-within-clusters; (b) --exonless-span; (c) both.
Reference: fold-first representative-only, min_size 2: detection 272/362 = 0.751, band-[0.90,1) precision
506/533 = 0.949, recall|both 0.942, family exact 51/67. Record-level (no fold): detection 0.843, precision 0.848.
Predictions (≥50 % floor):
- P1 (a): detection ≥ 0.80 (recovers the folded members) AND band precision ≥ 0.94 (nested duplicate records of
  one gene are still folded, so record-level's 107 FP do not return). Family exact ≥ 51.
- P2 (b): detection rises by ≥ 5 members over the reference with band precision ≥ 0.94.
- P3 (c): detection ≥ 0.82, band precision ≥ 0.94, recall|both ≥ 0.94; the Soto-silent band's asserted count
  does not grow by more than 20 %.
- P4 gorilla 3-contig with (c): NPIP stays one family of 29 loci (±2) with LCR16u separate (the block certificate
  survives), and every anchored cluster (tandem 60+58 array, the 3 artefacts) persists.
Fail rule: ¬P1 ⟹ the fold order is not the mechanism; report and stop. P4 fails ⟹ the fix trades the block
certificate for recall and does not ship.
