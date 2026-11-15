# PREREG r1620: tie-group loci by correlation clustering of the records (2026-11-10)

User (2026-11-10): try correlation clustering (Strainy-like) on the tie groups. The open defect (r1612) is records -> loci:
the span rule under-splits (a long seed swallows what lies in its intron), the block rule against one seed over-splits
(interleaved shapes open a locus each time). Loci -> groups (>= 2 shared tied molecules) is not open and is unchanged.
Written before the rule was run on either sample.

## Rule (one change, `tie_groups.py build ... cc`)
Records with the same chromosome, strand and aligned blocks are one shape (weight = its records). Two shapes ATTRACT when their
blocks overlap by >= 50% of the smaller block total (the r1612 criterion), REPEL when their spans overlap and they do not attract.
Shapes are taken heaviest first (ties: position); a shape joins the locus with the largest positive vote (sum over ALL its
member shapes of +weight for attract, -weight for repel) and opens a new locus when no vote is positive. Deterministic.
Recipe otherwise frozen: `build PREFIX as98 0 1 RMSK_BED:0.5 stranded`, scorer `full`, SEDEF co-duplicated.
Self-checks (`test_tie_groups.py`): the long-intron case of r1612 still passes; 2 shapes alternating 10 times give 2 loci.

## Defect measures (new; the old bars passed r1612 while it was broken)
- D1 duplicate models: models whose blocks overlap another same-strand model's blocks by >= 50% of the smaller, share of all models.
- D2 swallowed loci: loci whose model span is > 5x the median span of its records, share of all loci.
Both are printed by `build` for all three rules (span / block / cc) on the same inputs, before the repeat filter.
Audit references: block D1 15,791 of 20,008; span D2 489 of 5,541.

## Bars (both samples; A119b `r1603/recs.slim.tsv`, gorilla OR6737 `r1611/recs.slim.tsv`)
- D1 <= 0.05 and D2 <= 0.01.
- The r1607/r1611 bars unchanged: B1g grouped reach (scorer's bar), B2s >= 0.90, B3 TAF11L <= 5 kb >= 20 (A119b), B5r, B4 <= 1 h.
- PASS on both: cc is the recipe's locus rule and its figures replace r1607/r1611. FAIL: reported; locus rule stays open.
Gorilla is only partly unseen: its r1612 block output was looked at.

## Outcome (2026-11-10, job 4458996, `r1620/`): FAIL as registered
| sample, locus rule | loci | D1 duplicate models (<= 0.05) | D2 swallowed loci (<= 0.01) | B1g | B2s | B3 |
|---|---|---|---|---|---|---|
| A119b cc | 12,256 | **837 = 0.068 FAIL** | 41 = 0.003 | 69 of 72 | 0.999 | 30 of 31 |
| A119b span | 5,541 | 319 = 0.058 | 8 = 0.001 | | | |
| A119b block | 27,919 | 17,877 = 0.640 | 19 = 0.001 | | | |
| gorilla cc | 1,579 | 25 = 0.016 | 4 = 0.003 | 127 of 128 | 1.000 | - |
| gorilla span | 1,209 | 31 = 0.026 | 2 = 0.002 | | | |
| gorilla block | 2,065 | 543 = 0.263 | 0 = 0.000 | | | |
- Gorilla passes every bar; A119b misses D1 (0.068), so cc is not adopted and the locus rule stays open. No rule passes D1 on
  A119b: span 0.058, block 0.640. cc removes block's over-splitting (0.640 -> 0.068) and keeps 2.2x the loci of span.
- B4: A119b build 3,290 s (bar 1 h, passes by 310 s); gorilla 223 s. B5r passes on both (largest group on 1 chromosome).
- Self-check note: under cc the long read of the fixture joins the copy its second block covers by half, so the expectation in
  `bench/famcmp/test_tie_groups.py` differs from block's; that is the rule working as written.
