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
