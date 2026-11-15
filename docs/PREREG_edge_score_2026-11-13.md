# Row 1675 — one weighted score over the row 1674 features, human fit, gorilla test (2026-11-13)

User choice after row 1674's STOP: "Weighted score, gorilla test". Families layer (O1), opt-in only if found, not for
`advisor_share`. Registered before any score is computed. Seen already: row 1674's tables and its in-sample ceiling for
threshold conjunctions (20 true pairs cost precision 0.59-0.69), so the odds are low and that is stated here.

**Model (frozen, `bench/denovo_eval/edgescore.py`).** L2 logistic regression (IRLS, lambda 1.0, 30 steps) on the 20
row 1674 features, standardized; aln_bp, n_records, exon_bp, rank_max, degree_min, degree_max enter as log1p. No
interactions, no feature selection, no class weights. Population as row 1674.

**Fit.** A119b + HG002 pooled (`r1674_edgerule/feat.{A,H}.tsv`). Threshold = the lowest score at which the 5-fold
out-of-fold ranking of the training pool has precision >= 0.75 with >= 20 true pairs (folds by CRC32 of the gene pair).
**Go/no-go:** no such threshold -> STOP, gorilla not scored.

**Test (frozen weights and threshold, no refit).** OR6737 + KB3781 (`feat.{OR,KB}.tsv`, protein-hit truth, same-human-
gene pairs set aside): pooled new true pairs >= 20 at precision >= 0.60. Void if pooled reachable true < 50 (it is 210).

**If it passes:** row 1674's step 2 (family level through `--extra-edges`, same bars) decides whether a rule is found.
**If it stops or fails:** the row closes; no further model on these features.

## Result: STOP at the go/no-go (2026-11-13, `r1674_edgerule/score.out`)

Training pool 370,077 pairs, 277 true. No out-of-fold threshold reaches precision 0.75 with 20 true pairs.
Out-of-fold ranking: 20 true with 13 false (0.606), 40 with 55 (0.421), 80 with 147 (0.352). Gorilla not scored.
The row closes; no further model on these features.

**Descriptive, after the verdict (`r1674_edgerule/score_related.out`): what the false pairs are.** Of the
out-of-fold false pairs, those that are older Compara paralogs or share an HGNC gene group ("related", pairscore's
definition):

| top of the ranking holding | false | related | unrelated | precision as scored | with related set aside |
|---|---|---|---|---|---|
| 20 true | 13 | 9 | 4 | 0.606 | 0.833 |
| 40 true | 55 | 45 | 10 | 0.421 | 0.800 |
| 80 true | 147 | 126 | 21 | 0.352 | 0.792 |
| 120 true | 273 | 225 | 48 | 0.305 | 0.714 |

The score separates homologous gene pairs from unrelated ones (about 4 in 5 of its picks beyond the truth are older
relatives); what it cannot do is date the duplication, which is what the simiiformes truth asks. Unrelated examples
at the 40-true point: NBPF x MSTRG.359, ZBED6 x ZC3H11B, CCDC144A x PRPSAP2, FOXO3B x ZNF286A.
This does not change the verdict: under the registered scoring the pairs are false.
