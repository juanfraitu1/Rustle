# PREREG EXON-SPACE PLACEMENT, MULTI-ISOFORM REFERENCE (register 1643, 2026-11-12)

Why: r1642 failed with one spliced model per copy; 25 of its 44 losses to the primary were reads the single model could
not hold (no record or cover < 0.90), 19 were sequence. The user asked for the untested arm: several sequences per copy.

Change (only the reference; set, truth, rule, scorer `bench/o2_contested/exon_place.py` and minimap2 line as r1642):
`bench/o2_contested/multi_iso_ref.py` writes, per catalog copy of `kb.fam.copies.tsv`:
- `cat`: the catalog's own spliced sequence (r1642's reference),
- `tx`: every `kb.gtf` transcript (the sample's own de novo assembly) overlapping the copy's span and no copy at another
  locus, spliced; a transcript over two loci is left out (counted),
- `pre`: the copy's genomic span unspliced (unspliced and intron-retaining reads).
A copy's score = the best AS over all its sequences.

Expected, stated before the run: bar 2 still fails. At most the 25 model-limited losses can turn (15 vs 19 at best, bar
needs 2 to 1); bar 4 should fail again because the reference still holds no sequence from outside the catalog.

Bars: r1642's five, unchanged (reach >= 30%; arm right / primary wrong >= 2 x the reverse and >= 10; same ratio on
held-out alone and minus the largest family, >= 5; negative molecules assigned <= 10%; wall <= 10 min on 8 threads).
Readout: r1642's, plus the r1642 three-way split of the losses and which sequence kind carried the winning record.
Verdicts as r1642. A second FAIL closes exon-space re-alignment for placement whatever the reference.

Outputs: `/anvme/workspace/v138da13-rustle/r1643_multiiso/`.

## RESULT r1643 (2026-11-12, login node, 9 s with the reference build; `r1643_multiiso/result.txt`)
Reference: 5,518 sequences (812 cat, 812 pre, 3,664 transcripts placed, 55 over two loci left out).
- Bar 1 PASS: assigns 1,276 of 1,570 (0.813; r1642 0.699); truth copy without a record 21 (r1642 172).
- Bar 2 FAIL: arm right / primary wrong 31, the reverse 36. Assigned molecules: arm 1,192 / 84, primary 1,197 / 79.
- Bar 3 mixed: held-out alone 22 vs 8 (meets the ratio); minus MCL115 31 vs 36 (fails). Panel 9 vs 28.
- Bar 4 FAIL: 833 of 1,419 negative molecules assigned (0.587; r1642 0.456, O2 0.162).
- Bar 5 PASS.
**FAIL as registered**: with several sequences per copy the arm equals the primary (and O2: 494 vs 493 right on 529
shared), it does not beat it, and it assigns more of the reads whose origin is outside the catalog.
Readout: margin >= 10 gives 25 vs 27; cover gate off 42 vs 48. Unspliced primaries 12 vs 3, spliced primaries 19 vs 33.
MCL9 alone is 0 vs 15 (r1642: 0 vs 0), MCL5 0 vs 8. Top record kind per read: tx 1,495, pre 743, cat 702. The
three-way split of the losses was not rerun.
Exon-space re-alignment is closed for placement with either reference. Not followed up (post hoc if ever): the
held-out 22 vs 8 and the unspliced-primary class.
