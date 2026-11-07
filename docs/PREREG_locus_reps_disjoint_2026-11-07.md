# PREREG row 1575: multiple representatives per locus with exon-disjoint groups (fix of row 1574's grouping)

## Why
Row 1574 FAILED on a mechanism bug: readthroughs kept in one group made the groups' exon unions overlap, and the
families stage folds such loci into one node, so the secondary representatives vanished. Fifth variant on A119b
(development evidence, designed after seeing 1571-1574; held-out = option 1, next).

## Arm F
`locus_reps.py --disjoint`: as row 1574, but a transcript sharing bases with >= 2 representatives is dropped, then any
transcript touching another group's exon union, until the groups' unions are disjoint (families input only). Families
stage (HEAD mcl_families, minimap2-gate), piece_merge.py, piece_attach.py onto row 1571 arm A's unchanged families.
Descriptive: F's groups attached onto row 1573's D table (union).

## Bars (F vs A), as row 1571
1. Of A's 71 absorbed pairs, <= 35 still missed.  2. Precision >= A - 0.01.  3. Recall >= A.  4. Copies <= 1.05 x A.

## Prediction
Bar 1 FAIL (17 genes reachable at step 0, but attach needs a family hit; ~20 pairs), bar 2 borderline, 3-4 PASS.

## Result (2026-11-07, job 4452574, `/anvme/workspace/v138da13-rustle/r1571/F*`, `Fa`, `DF`)
**FAIL (bars 1, 2).** Fa: 0.787 / 0.772 (1,686 hits, +71; 0 pairs A hit lost), copies 8,251 (PASS); A's 71 absorbed:
9 recovered, 62 still missed. Union DF (D + F attached): 0.779 / 0.776, 16 recovered, copies 8,302, 0 lost.
Trace of the 25 genes: 17 now have a node of their own (vs ~0 in 1574: the fold is fixed); 10 reach a family in F's own
clustering; 4 survive the attach. Lost at two identified steps: (a) piece_merge/attach keeps only the most-read piece when
two pieces of one locus land in one family without aligning to each other (NPIPB7 loses to the overlapping lncRNA group,
TCP11X1/X2, FAM25G, OR4M2, ZNF324B: 6 genes); (b) 8 nodes get no family edge at all (UGT2B17, ZNF585B, ZNF528, ZNF880,
KLRC2, TBC1D28, C1QTNF9B, AC092329.3), consistent with the overlapping-span / self-hit rejection between pieces of one
original locus (not verified). Prediction right on bar 1, bar 2 fails by 0.005.
