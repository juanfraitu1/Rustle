# PREREG row 1754 (2026-11-15): find a skipped exon of <= 15 bp by searching the intron for the read's extra bases

Row 1751 (register 1753) failed on reach: an exon dropped in every read has no junction in the locus to realign to
(reach 0.26-0.52 of row 1746's exons). This row searches the intron itself. Written and committed before the search was
run on any sample (one small gorilla contig was run for mechanics only, output not read). Rows 1746 and 1751 results on
all four samples were already seen, so no sample is unseen here; nothing in the rule is tuned on them.

## Rule (`bench/famcmp/skip_search.py`)

- A primary read's intron (d, a) with canonical ends and L = 3..15 extra read bases in the junction window (row 1751's
  window, 25 bp each side, clamped to the read's flanking exons).
- Sites: every L bp stretch of the intron with the intron's own splice signals around it (AG before, GT after, on the
  intron's strand) and >= 50 bp of intron on each side.
- The search is used only when sites x 4^-L <= 0.01 (expected chance matches per read). Short exons in long introns are
  therefore not searched.
- A site is taken only when it is the single one whose sequence is in the read, with <= 1 edit left in the window and
  none of the L extra bases unexplained.

## Step 0 (read level, four samples)

- **Gate (decoy):** the same search on the reversed intron. Reads placed by the decoy / reads placed by the search
  <= 0.05 on each sample, else INVALID.
- **Precision:** found units whose outer junctions (d, a) flank an annotated exon of <= 15 bp: the found exon is the
  annotated one in >= 0.95 on each sample with >= 20 such units.
- **Reach:** row 1746's exons (>= half the skip reads carry the exon; 26 / 27 / 34 / 27) reached by this search or
  by row 1751's rule: >= 0.75 on >= 3 of the 4 samples. The search alone and the lengths of the exons not reached are
  reported.
- Reported, no bar: found units by annotation; "neither" = unannotated short exons.

## Step 1 (assembly; read only if step 0 passes all three)

As row 1751: panel = reads within 100 kb of a fired unit, arms `base` and `fix` (`skip_fix.py fix` on the units of both
rules), default assemble, same binary. S1 skip chains at exon-only annotated fired units fall by >= 50%; S3 panel
precision of `fix` >= `base` - 0.0005. S2 is changed from row 1751 on what its losses turned out to be (all but 2 of
15 lost correct chains were annotated skip forms at units where both forms are annotated, their reads moved to the exon
form): lost correct chains that are not such skip forms <= new correct chains. PASS = S1-S3 on >= 3 of 4 samples. A PASS
makes the two rules together a candidate for an opt-in stage (user decision). No default changes here.

## Limits

Exons of 3-5 bp are almost never searchable (the chance bound), so reach is expected to miss them. Only an exact (d, a)
skip intron with the bases inserted next to it is reached. The rewrite's read test (2-edit margin) is looser than the
search, so it also moves reads of the same intron that carry the exon with an error.

## OUTCOME (2026-11-16, job 4473405, `r1754_skipsearch/{or,kb,a,hg}.{search,prec,score}.txt`): FAIL on precision

| | OR6737 | KB3781 | A119b | HG002 |
|---|---|---|---|---|
| gate: decoy / found reads (<= 0.05) | 0.0010 | 0.0001 | 0.0018 | 0.0003 |
| precision vs the annotated exon (>= 0.95) | 24 / 25 (0.960) | 15 / 19 (under the 20 floor) | 114 / 116 (0.983) | 79 / 87 (0.908) |
| reach, search / row 1751 rule / either (>= 0.75) | 14 / 9 / 23 of 26 (0.885) | 10 / 14 / 24 of 27 (0.889) | 16 / 9 / 25 of 34 (0.735) | 12 / 12 / 24 of 27 (0.889) |
| units found (neither form annotated) | 252 (152) | 300 (202) | 955 (670) | 712 (458) |

Gate valid. Reach passes (3 of 4); the two rules reach different exons, and those left are 3-6 bp. Precision fails on
HG002, so the row fails as registered and step 1 decides nothing. Rewrite check: differing records = logged, all four.

Post hoc:

- The 12 disagreeing finds (KB3781 4, HG002 8): 10 share a boundary with the annotated exon or sit 2-3 bp inside it;
  8 of the 12 rest on one read. Two are a pair of near-identical 6 bp exons 516 bp apart (TCACAA / TCACAG) in the same
  gene of both species; on KB3781 1,253 reads carry the unannotated one. A floor of 2 reads per unit is the obvious
  fix; not tested, and no sample is unseen for it.
- No exon over 10 bp was found on any sample. Checked on OR6737 units row 1751 fixed: the reads that drop an 11-15 bp
  exon carry 2 or more errors in the window and fail the 1-edit limit; the code finds the site. Row 1751's rule covers
  that length range.
- Step 1 on the panel:

| | OR6737 | KB3781 | A119b | HG002 |
|---|---|---|---|---|
| reads rewritten | 18,524 | 76,924 | 30,563 | 59,344 |
| skip chains at exon-only units, base -> fix | 72 -> 0 | 57 -> 0 | 210 -> 20 | 112 -> 13 |
| chains with the exon at those units | 45 -> 107 | 66 -> 116 | 140 -> 326 | 104 -> 200 |
| new correct / lost correct chains | 14 / 20 | 14 / 14 | 96 / 18 | 73 / 11 |
| panel precision, base -> fix | 0.3279 -> 0.3263 | 0.3442 -> 0.3436 | 0.1686 -> 0.1703 | 0.2358 -> 0.2375 |

Which of the lost correct chains are annotated skip forms (the changed S2) was not computed. No stage built. The panel
BAMs were deleted; `*.fix.bam.fixes.tsv` keeps every rewritten read.
