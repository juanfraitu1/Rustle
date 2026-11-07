# PREREG row 1571: locus pieces for the families input (genes with no locus of their own)

## Why (step 0, A119b a1, `bench/denovo_eval/pairdiag.py`)
71 missed Compara pairs come from 25 expressed genes that have no locus of their own: their transcripts sit in a locus
whose representative maps to another gene (12 fused neighbours: ZNF585A/B, ZNF443/ZNF709, UGT2B15/17, TCP11X1/NXF2 ...;
13 overlapping or nested genes: ZBED6 in ZC3H11A, KLRC1/2, AMY2B/RNPC3 ...). Step 0 on read structure alone: a
single-transcript articulation rule reaches 3/25 (fused loci have parallel readthroughs); a positional cut reaches 9/25
(16/25 without a read-minority condition) but 71-79% of its cuts split ONE annotated gene, and no read feature separates
the two (best: crosser share < 0.1, 41% two-gene; odd chromosomes). r1516/1517 (tandem cuts with a self-homology gate)
reached 5-6 of 69. New here: cut permissively, then let the family stage's own sequence evidence decide.

## Arms (one job, mcl_families built at HEAD, the families stage command of tools/rustle_pipeline.sh)
- A: `a1.families.gtf` unchanged.
- B: `bench/denovo_eval/locus_pieces.py` on `a1.families.gtf`: in each gene_id, a cut at the midpoint of every gap between
  merged exon blocks where both sides hold a transcript of >= 2 reads lying wholly on that side; transcripts crossing a
  cut are dropped; pieces become gene_id GENE.pK (families input only; the assembly GTFs are untouched). Then the families
  stage, then `bench/denovo_eval/piece_merge.py`: pieces of one original gene_id in the same family stay separate copies
  only when their loci share a gated-PAF record of >= 300 bp at >= 0.70 identity (single linkage); otherwise each
  linked group keeps one row, its piece with the most reads.

## Bars (B vs A, pairscore simiiformes, same-symbol ignored; pairdiag classes)
1. Absorbed missed pairs (pairdiag `absorbed`) B <= 0.5 x A.
2. Pair precision B >= A - 0.01.
3. Pair recall B >= A.
4. Copies B <= 1.05 x A (over-splitting guard).
PASS = all four -> proposed as an opt-in families-input layer, held-out on a second sample before any default change.

## Disclosed limits
Development evidence on A119b (the 25 genes were seen). Nested and overlapping genes (13 of 25) cannot be separated by a
positional cut by construction. O2 consumes the copy table: B's extra copies would change assignment; not measured here.

## Prediction
Bar 1 FAIL (reach on the 25 is at most 16, and merged-back pieces lose some), bars 2-4 PASS.

## Result (2026-11-07, jobs 4452474 arm A / 4452502 arm B, `/anvme/workspace/v138da13-rustle/r1571/`)
**FAIL (bars 1, 3).** Arm A (HEAD mcl_families, minimap2-gate) reproduces the base run's score exactly (0.802 / 0.740,
1,615 hits; 7,986 copies vs 7,987 in October: today's duplicate-locus fix). Arm B: 3,328 loci cut into 6,693 pieces,
54,541 crossers dropped; piece_merge rejoined 3 rows (8 family x locus groups with >= 2 pieces).
1. A's 71 absorbed pairs: 27 recovered, 44 still missed (21 absorbed, 14 split, 9 unlinked). FAIL (bar <= 35).
   (The pairdiag absorbed class alone falls 71 -> 21; strict reading reported, since a merged-away piece leaves its
   gene "unlinked", not recovered.)
2. Precision 0.800 vs 0.802: PASS.  3. Recall 0.718 vs 0.740: FAIL.  4. Copies 8,347 vs <= 8,385: PASS.
Collateral: 154 pairs hit by A are missed by B, 152 of them `split`: replacing full-length representatives by partial
pieces cut family edges. Net -48 hits. Prediction right on bar 1, wrong on bar 3.
