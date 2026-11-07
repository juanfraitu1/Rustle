# PREREG row 1572: locus pieces confined to loci with family evidence (follow-up to row 1571)

## Why
Row 1571 cut 3,328 loci: 27 of the 71 absorbed pairs recovered, but 154 pairs A hit were lost (152 family splits, net
-48). Designed post hoc from 1571's output on the same sample (development evidence, disclosed). Step 0 on 1571's arm B:
154 loci keep >= 2 pieces as family copies, 174 keep exactly 1; the other ~3,000 cut loci have no piece in any family,
so cutting them can only cost.

## Arm C (two-pass)
Pass 1 = row 1571 arm B (`/anvme/workspace/v138da13-rustle/r1571/Bm.fam.copies.tsv`). S = original gene_ids with >= 1
piece among Bm's copies (328). Pass 2: `a1.families.gtf` with only the loci in S replaced by their pieces (the same
locus_pieces.py cuts, crossers dropped), then the families stage (HEAD mcl_families, minimap2-gate), then piece_merge.py.
Every locus outside S is exactly as in arm A.

## Bars (C vs row 1571 arm A), as row 1571
1. Of A's 71 absorbed pairs, <= 35 still missed.  2. Precision >= A - 0.01.  3. Recall >= A.  4. Copies <= 1.05 x A.

## Prediction
Bar 1 FAIL (1571 recovered 27; C cannot reach more), bars 2-4 PASS with a small net hit gain.

## Result (2026-11-07, job 4452544, `/anvme/workspace/v138da13-rustle/r1571/C*`)
**FAIL (bars 1, 3).** 328 loci cut; copies 8,286 (PASS <= 8,385); precision 0.806 (PASS); recall 0.724 < 0.740 (FAIL);
A's 71 absorbed pairs: 27 recovered, 44 still missed (FAIL). 151 pairs A hit are lost (150 split): NPIPB 54, GOLGA 40,
NPIPA 36. 137 of A's 1,565 families are broken in C and 77 of those hold no cut locus: re-running MCL on a graph with a few
extra piece nodes redraws the large NPIP/GOLGA families. The collateral is MCL's, not the cuts'.
