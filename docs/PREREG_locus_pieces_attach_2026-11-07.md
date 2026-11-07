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
