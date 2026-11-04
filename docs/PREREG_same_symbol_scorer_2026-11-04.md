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
