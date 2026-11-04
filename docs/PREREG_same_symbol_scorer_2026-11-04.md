# PREREG r1530: pairscore's "false" pairs include same-symbol CAT duplicate ids; correct the scorer, re-score r1528 (2026-11-04)

Finding (r1529, after the numbers): 213 of the base catalog's 611 false coding gene pairs are SAME-SYMBOL pairs (AMY1A~AMY1A,
GAGE12C~GAGE12C, CDY1~CDY1): CAT/Liftoff gives identical copies duplicate gene ids under one symbol; Compara has no pair to name,
so pairscore counts them as false. This row corrects the measurement; it moves no verdict of r1528 (FAIL as registered stands).

## Rule (fixed before any number is computed)
- `PAIRSCORE_SAME_SYMBOL=ignore`: a predicted coding pair whose two genes carry the same symbol is IGNORED (neither hit nor
  false), reported in a new column `ignored_same_symbol`. Everything else unchanged.
- Adoption check, first: of the base's 213 same-symbol false pairs, the share whose two CAT genes map to the SAME Ensembl source
  gene id (G[x]['ens'] intersects G[y]['ens']). **Bar: >= 0.90.** Below it, the pairs are not provably the same gene; the rule is
  NOT adopted and the row closes with the share recorded.
- If adopted: re-score base (0.60), f0.45, f0.30 from /anvme/workspace/v138da13-rustle/fam1528 with the corrected scorer.
  Bars for an arm, mirroring r1528: corrected precision >= corrected base - 0.01; recall, covered members, largest family, wall
  as already measured in r1528 (recall >= 0.771, covered >= 711, largest <= 330, wall <= 900 s).
- Prediction written down: corrected base precision ~ 1617 / (2228 - 213) = 0.80; the arms stay below base - 0.01 because their
  extra 197 / 259 false pairs were two-thirds domain-level relatives (r1528 result table), not same-symbol duplicates.
  If an arm passes every bar under the corrected scorer, that is a NEW result for a new row (the gate is not changed here).

## RESULT (~/r1529/r1530_rescore.out)
Adoption check: 213 of 213 same-symbol false pairs share an Ensembl source id = 1.000 (bar >= 0.90) -> rule adopted.
| catalog | pred pairs | hit | false | precision as registered | corrected precision | recall | ignored_same_symbol |
|---|---|---|---|---|---|---|---|
| 0.60 (base) | 2228 -> 2013 | 1617 -> 1615 | 611 -> 398 | 0.726 | **0.802** | 0.740 | 215 |
| 0.45 | 2451 -> 2236 | 1643 -> 1641 | 808 -> 595 | 0.670 | 0.734 FAIL (< 0.792) | 0.752 | 215 |
| 0.30 | 2571 -> 2356 | 1701 -> 1699 | 870 -> 657 | 0.662 | 0.721 FAIL (< 0.792) | 0.778 | 215 |
215 pairs ignored (213 false + 2 truth pairs that Compara did name under one symbol; hit 1617 -> 1615, recall 0.741 -> 0.740).
The arms' extra false pairs are untouched by the correction (0 extra same-symbol pairs), as predicted. **r1528 verdict stands.**
Scorer: `PAIRSCORE_SAME_SYMBOL=ignore` is the figure to report; base catalog pair precision 0.802 / recall 0.740 (primates).
