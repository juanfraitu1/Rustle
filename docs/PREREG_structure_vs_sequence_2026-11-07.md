# PREREG row 1570: which family definition agrees better with Compara, sequence (O1) or exon architecture?

## Why
Row 1569: O1 families share exon architecture far above a null, but structure alone looked weak below 95% identity.
User: test which is the better DEFINITION against Ensembl Compara paralogues. Scoring is the standing readout
(`bench/denovo_eval/pairscore.py`, level simiiformes, `PAIRSCORE_SAME_SYMBOL=ignore`; baseline reproduced today:
O1 full 0.802 / 0.740 on 2,183 truth pairs).

## Universe (same loci for every arm)
A119b base run `~/runs/human_A119b/a1.*`. U2 = loci whose representative has >= 4 exons (n_int >= 2 internal exons;
one internal exon collides by chance too often to define anything). Representative exons: `a1.fam.copies.tsv` for
family copies, `a1.fam.loci.gff3` otherwise. pairscore `--gtf` = `a1.families.gtf` restricted to U2 gene_ids, so
truth pairs are Compara pairs among protein-coding genes reached by U2 loci.

## Arms (copies-format TSVs, one row per U2 locus in a family of >= 2)
- SEQ: O1 families restricted to U2.
- STRUCT_EXACT: U2 loci with identical internal-exon length tuples (5'->3').
- STRUCT_TOL: same n_int and every internal exon within 3 bp, single linkage (union-find).
- AND: each O1 family (in U2) split into connected components of row-1569 ALIKE edges (LCS within 3 bp, S >= 0.5).

## Rule
An arm BEATS SEQ if precision is higher and recall is not lower by more than 0.01, or recall is higher and precision
is not lower by more than 0.01. SEQ BEATS an arm by the mirror rule. Otherwise: trade-off, no winner.
Headline question: does either STRUCT arm beat SEQ? (structure as the better definition)

## Disclosed limits
Compara builds trees from protein sequence: closer to SEQ's evidence than to STRUCT's, a residual tilt toward SEQ.
U2 leaves out 1-3 exon genes (processed pseudogenes, fragments), where structure has nothing to say.
Exon boundaries in every arm come from the same read alignments.

## Prediction
SEQ beats STRUCT_EXACT and STRUCT_TOL (lower precision: unrelated genes share architecture; lower recall: diverged
paralogs gain/lose exons). AND vs SEQ: higher precision, recall lower by > 0.01, so a trade-off.

## Result (2026-11-07, `bench/structfam/structdefs.py` + pairscore; `~/runs/structfam_r1570/pairscore.tsv`)
**SEQ (O1) is the better definition; SEQ beats both STRUCT arms; AND vs SEQ: no winner (identical within 0.001).**
U2 = 15,251 loci (1,062 O1 members); truth 1,468 Compara pairs (universe 12,292).

| arm | families | loci | pred pairs | hit | false (related) | precision | recall |
|---|---|---|---|---|---|---|---|
| SEQ | 260 | 889 | 1,334 | 1,180 | 154 (140) | 0.885 | 0.804 |
| STRUCT_EXACT | 441 | 1,180 | 2,101 | 852 | 1,249 (1,198) | 0.406 | 0.580 |
| STRUCT_TOL | 633 | 3,449 | 177,624 | 962 | 176,662 (6,394) | 0.005 | 0.655 |
| AND | 246 | 848 | 1,332 | 1,179 | 153 (140) | 0.885 | 0.803 |

Prediction right on both STRUCT arms; wrong on AND (no precision gain: O1 families in U2 are already structurally
coherent, splitting them removed 2 pairs). STRUCT_TOL single linkage chains into one 1,136-locus component.
Post hoc reading: 1,198 of STRUCT_EXACT's 1,249 false pairs are relatives (Compara at an older level, same HGNC group
or symbol), so shared architecture tracks deep homology, older than the simiiformes truth level; it does not separate
young families. Structure stays a corroborating view (row 1569), not a definition.
