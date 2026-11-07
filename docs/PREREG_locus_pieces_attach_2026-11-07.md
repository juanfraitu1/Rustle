# PREREG row 1573: attach locus pieces to the unchanged families (no re-clustering)

## Why
Rows 1571/1572 recover 27 of the 71 absorbed pairs but lose ~150 pairs because MCL, re-run with extra piece nodes, redraws
the large NPIP/GOLGA families (77 of 137 broken families hold no cut locus). Third variant on A119b, designed post hoc
from 1571/1572 (development evidence; held-out needed before any opt-in).

## Arm D (`bench/denovo_eval/piece_attach.py`)
Start from arm A's copy table (families unchanged). For every piece copy of row 1572's arm C after piece_merge (Cm; loci
in S only): target = the A family sharing the most original gene_ids with the piece's C family, excluding the piece's own
locus (ties: lower family id; no shared member -> the piece is not attached). In each target family the parent locus's A
row, if present, is removed and the attached pieces of that locus are added; two pieces of one locus in one target family
stay separate only when they share a gated-PAF record >= 300 bp at >= 0.70 identity in C's PAF (piece_merge rule), else
the most-read piece. No family is created or split.

## Bars (D vs A), as row 1571
1. Of A's 71 absorbed pairs, <= 35 still missed.  2. Precision >= A - 0.01.  3. Recall >= A.  4. Copies <= 1.05 x A.

## Prediction
Bar 1 FAIL (reach <= 27), bars 2-4 PASS (net gain ~+20 hits).

## Result (2026-11-07, `/anvme/workspace/v138da13-rustle/r1571/D.fam.copies.tsv`)
**FAIL (bars 1, 2).** A 0.802 / 0.740 -> D 0.789 / 0.774 (1,615 -> 1,689 hits, +74; false 398 -> 451); copies 8,184.
1. A's 71 absorbed pairs: 16 recovered, 55 still missed: FAIL.  2. Precision 0.789 < 0.792: FAIL (by 0.003).
3. Recall +0.034: PASS.  4. Copies <= 8,385: PASS. Collateral gone: 3 pairs A hit are lost (AGAP4/6/9, NXF2/NXF2B) vs
~150 in rows 1571/1572. Most of the +74 are not absorbed pairs: second genes of fused loci that were unlinked or split.
Across rows 1571-1573, 44 of the 71 absorbed pairs are never recovered: 35 involve overlapping or nested genes (a
positional cut cannot separate them by construction), 9 fused neighbours. Prediction wrong on bar 2.
Worth a held-out (another sample) as an opt-in families post-step: the recall gain is large, the precision cost 0.013.
