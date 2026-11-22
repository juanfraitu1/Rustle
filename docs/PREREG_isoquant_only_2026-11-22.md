# Row 1896: which annotated chains does IsoQuant report that Rustle lacks, on twenty datasets? (step 0, pre-registered 2026-11-22)

User (2026-11-22): "is there anything from isoquant that we could try?" Record: IsoQuant was read as a source of rules four
times (rows 1765, 1838, 1853, 1854). Its rules were mapped one by one against Rustle's (row 1838's table); its precision is an
operating point (Rustle cut to its correct count is 6-16 points more precise on real libraries); the two differences in kind
were tested and failed (a contained chain kept by its own read-end peak: rows 1840, 1841, 1843; the 2-exon MAPQ rule: row
1854, which would drop 4-6% of the correct 2-exon chains of duplicated genes on apes). Those rows asked "why is IsoQuant more
precise". None asked the recall-side question on all datasets, and IsoQuant outputs now exist for twenty (row 1880).

## Measure (existing files, no run; `r1896_isoquant_only/iq_only.py`)
Per dataset: IsoQuant's correct chains; how many are missing from Rustle's `final.gtf`; how many are missing from the whole
ranked list (`final.gtf` + one-read tier). The ones missing from the whole list, by the first class that fits, against the
chains of Rustle's list:
1. PIECE: a contiguous part of a list chain (Rustle has the longer form only).
2. LONGER: a list chain is a contiguous part of it (Rustle has the shorter form only).
3. NEAR: equal to a list chain except one intron within 10 bp.
4. SHARES INTRONS: every intron in some list chain, combination absent.
5. NEW INTRON: an intron in no list chain.

## Reading, fixed now
- **NOTHING TO PORT ON THE RECALL SIDE** when on at least 17 of 20 datasets the correct IsoQuant chains missing from Rustle's
  whole list are under 2% of IsoQuant's correct chains.
- Otherwise a class holding at least half of the missing ones on at least 5 datasets in both groups is a LEAD, to be read
  against IsoQuant's source for the rule that makes those chains; a port would be its own row.
- Either way the answer to the user's question is given with this table and row 1838's rule mapping.

## Result (2026-11-22): a LEAD by the reading, and it is the sub-chain collapse again; a direct check finds nothing to port

`r1896_isoquant_only/q.NAME.txt` (`iq_only.py`), existing files of row 1880.
| dataset | group | IsoQuant correct | missing from `final.gtf` | missing from the whole list | PIECE | LONGER | NEAR | SHARES INTRONS | NEW INTRON |
|---|---|---|---|---|---|---|---|---|---|
| fly simulated | N | 6,592 | 190 (0.029) | 120 (0.0182) | 47 | 22 | 3 | 4 | 44 |
| mouse simulated | N | 24,098 | 1905 (0.079) | 1653 (0.0686) | 893 | 117 | 29 | 311 | 303 |
| zebrafish simulated | N | 11,503 | 513 (0.045) | 317 (0.0276) | 71 | 62 | 47 | 9 | 128 |
| mouse ES | N | 13,443 | 795 (0.059) | 765 (0.0569) | 568 | 45 | 10 | 102 | 40 |
| mouse heart | N | 11,755 | 438 (0.037) | 357 (0.0304) | 236 | 41 | 6 | 39 | 35 |
| mouse adrenal | N | 12,880 | 594 (0.046) | 509 (0.0395) | 336 | 47 | 11 | 74 | 41 |
| mouse cortex | N | 13,416 | 664 (0.049) | 578 (0.0431) | 393 | 65 | 8 | 65 | 47 |
| WTC11 r2 | P | 16,609 | 621 (0.037) | 437 (0.0263) | 219 | 96 | 21 | 38 | 63 |
| A119b | P | 18,903 | 1018 (0.054) | 662 (0.0350) | 370 | 185 | 13 | 17 | 77 |
| KB3781 gorilla fibroblast | P | 15,785 | 401 (0.025) | 321 (0.0203) | 117 | 117 | 10 | 21 | 56 |
| HG002 | P | 16,087 | 876 (0.054) | 668 (0.0415) | 270 | 240 | 6 | 27 | 125 |
| OR6737 gorilla testis | P | 16,677 | 442 (0.027) | 258 (0.0155) | 107 | 81 | 12 | 19 | 39 |
| Sumatran fibroblast | P | 15,667 | 409 (0.026) | 272 (0.0174) | 90 | 99 | 7 | 18 | 58 |
| Sumatran testis | P | 12,866 | 477 (0.037) | 215 (0.0167) | 78 | 95 | 8 | 8 | 26 |
| bonobo testis | P | 13,109 | 366 (0.028) | 185 (0.0141) | 80 | 62 | 5 | 13 | 25 |
| bonobo fibroblast | P | 13,139 | 581 (0.044) | 369 (0.0281) | 149 | 169 | 4 | 6 | 41 |
| Bornean fibroblast | P | 14,803 | 472 (0.032) | 312 (0.0211) | 105 | 130 | 7 | 14 | 56 |
| Bornean testis | P | 14,094 | 525 (0.037) | 280 (0.0199) | 119 | 105 | 5 | 16 | 35 |
| chimpanzee cell line | P | 16,579 | 536 (0.032) | 366 (0.0221) | 104 | 177 | 8 | 16 | 61 |
| chimpanzee testis | P | 15,761 | 600 (0.038) | 268 (0.0170) | 70 | 116 | 11 | 33 | 38 |

- **Not "nothing to port" by the count:** the correct IsoQuant chains missing from Rustle's whole list are under 2% on 7 of 20
  datasets (bar 17); they are 1.4% to 6.9% (mouse simulated), 215 to 1,653 chains.
- **The LEAD is PIECE:** at least half of the missing chains on 7 datasets in both groups (mouse simulated 54%, mouse ES 74%,
  heart 66%, adrenal 66%, cortex 68%, WTC11 50%, A119b 56%). On the non-human apes LONGER is as large or larger (the pattern of
  row 1883: computed annotations hold longer forms). NEAR is 0.5-15%, NEW INTRON 5-40%.
- A PIECE is an annotated chain that is a contiguous part of a chain Rustle reports: Rustle's sub-chain collapse removed it,
  IsoQuant kept it. This is the lever rows 1839-1843 closed (read ratio, own end peak, end peak against run-through reads,
  3'-side exemption), where IsoQuant's own criterion, the read-end peak, was rebuilt in Rustle and did not separate.

**Direct check, post hoc (`pieces.py`; not registered).** Those rows tested Rustle's rebuild of the criterion, not IsoQuant's
actual choices. So: among the chains the collapse removes from `final.gtf` (arm `ism0` minus `base`, row 1864's runs), are the
ones IsoQuant also reports more often annotated?
| library | collapsed chains | reported by IsoQuant too | not reported by IsoQuant | the same number taken by most reads | library precision |
|---|---|---|---|---|---|
| mouse ES | 1,193, 24.2% annotated | 586, 29.0% | 607, 19.6% | 586, 34.5% | 41.0% |
| WTC11 r2 | 554, 18.8% | 403, 19.6% | 151, 16.6% | 403, 22.3% | 63.9% |
IsoQuant's kept contained chains are a little more often annotated than the ones it drops, less often than the same number
picked by read count alone, and far below either library's precision. IsoQuant gains these chains by keeping them, not by
knowing which to keep. Rustle's switch for the same trade exists (`RUSTLE_POLISH_ISM_RATIO`, 0 = keep all).

**Answer to the question.** Nothing is left in IsoQuant to port: its rules have Rustle counterparts (row 1838's table), its
precision is an operating point (rows 1765, 1838, 1853), its 2-exon rule hurts duplicated genes (row 1854), and its recall
edge over Rustle's list is contained chains kept without a better selection than a read count (this row). IsoQuant is closed
as a source of rules by five rows.
