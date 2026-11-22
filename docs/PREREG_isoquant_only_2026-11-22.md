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
