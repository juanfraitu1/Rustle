# PREREG row 1574: more than one representative per locus (option 2 for the 71 absorbed pairs)

## Why
Rows 1571-1573: positional cuts cannot reach overlapping or nested genes (35 of the 71 absorbed pairs). Step 0 (A119b,
same 25 genes, development evidence): exon-disjoint secondary representatives reach 17/25 (12/12 fused, 5/13 overlapping;
positional cuts reached <= 16) and touch 2,558 of 76,161 loci; allowing 10% exon overlap adds no gene (not used). The 8
unreached overlap the primary through a readthrough or shared exons (AMY2B/RNPC3, ZC3H11A/ZBED6, NXF2B ...).

## Arm E
`bench/denovo_eval/locus_reps.py` on `a1.families.gtf` (greedy exon-disjoint representatives, >= 2 reads; groups
GENE.pK), families stage on the whole modified GTF (HEAD mcl_families, minimap2-gate), piece_merge.py, then
piece_attach.py onto row 1571 arm A's families (unchanged; no re-clustering; row 1573's rule).

## Bars (E vs A), as row 1571
1. Of A's 71 absorbed pairs, <= 35 still missed.  2. Precision >= A - 0.01.  3. Recall >= A.  4. Copies <= 1.05 x A.
Descriptive: E with row 1573's D pieces attached as well (union), same four numbers.

## Prediction
Bar 1 PASS (17 genes reachable cover > half of the pairs), bars 2-4 PASS. Held-out next (option 1).

## Result (2026-11-07, job 4452547, `/anvme/workspace/v138da13-rustle/r1571/E*`, `Ea`, `DE`)
**FAIL (bars 1, 2).** Ea (E attached to A): 0.788 / 0.756 (1,650 hits), copies 8,199; A's 71 absorbed pairs: 8 recovered,
63 still missed; 4 pairs A hit lost. Union DE (D + E attached): 0.779 / 0.775, 15 recovered, copies 8,282.
Cause (mechanism, not the idea): transcripts joining two representatives (readthroughs) were kept in one group, so the
groups' exon unions overlap (UGT2B17 714 bp, ZNF585B 363 bp, ZC3H11A 6,475 bp) and the families stage folds them into one
node ("fold-within-clusters"): most secondary groups never become a node of their own. Prediction wrong on bars 1 and 2.
