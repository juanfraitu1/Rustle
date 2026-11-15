# PREREG row 1663: a wider truth level for family precision (2026-11-13, user: "lets register the wider truth level")

**Why.** Row 1662 closed every edge lever under the simiiformes truth. Its last step (protein hit plus proximity)
failed on precision, and the false pairs shown were real but older paralogs (ATAD3A/B/C, FBXO2/6/44, AKR7A2/A3).
The simiiformes level scores a family that correctly holds an older paralog as wrong.

**What is registered.** A second track, W, reported beside the standing figure and never replacing it:

- recall: unchanged (Compara pairs at simiiformes among expressed coding genes, PAIRSCORE_SAME_SYMBOL=ignore);
- `precision_wide`: a predicted pair counts as correct when it is a Compara paralog pair at **eutheria** or younger
  (`famtruth.LEVELS['eutheria']`, the family level defined 2026-10-30). `pairscore.py --wide-level eutheria` adds
  the columns `false_older` and `precision_wide`. HGNC groups and Compara levels older than Eutheria do NOT count
  (pairscore's `of_which_related` does count them; that column stays descriptive);
- wide members: genes placed in a family beside at least one partner that is correct on track W.

Eutheria is fixed here, before any count at that level, and is not to be moved after.

**Already seen (so nothing below is blind on these):** A119b / HG002 counts of false pairs related at ANY level
(row 1662 audit: 110 of 165 on the A119b base table), and six example names from step 0f.

**Readings and bars.**

1. Descriptive, no bar: base vs default tables on A119b (`asg_def`), KB3781 (`kb_all`), on track W.
2. The shipped default re-judged: PASS on track W if `precision_wide` of the default table >= that of the base table
   on A119b AND on KB3781 (protein-hit truth). I expect A119b to fall slightly, from the audit above.
3. Step 0f again on track W (protein hit, same chromosome within 0.2 / 1 / 5 Mb, identity floors 50 / 70, A119b,
   `port1593` tables). Go on to a prototype only if some cell gives >= 30 genes newly placed beside a wide-correct
   partner at new-pair `precision_wide` >= 0.75. Otherwise STOP. A go is a prototype on A119b with gorilla as the
   transfer test, behind a flag.

Standing figure for the thesis stays the simiiformes one. Family layer only; nothing here touches read assignment.

## Results (same day)

| table | precision | recall | false | of which eutheria or younger | precision_wide |
|---|---|---|---|---|---|
| A119b base | 0.912 | 0.782 | 165 | 10 | 0.917 |
| A119b default | 0.902 | 0.808 | 192 | 10 | 0.907 |
| KB3781 base (protein-hit) | 0.205 | 0.408 | 379 | 2 | 0.210 |
| KB3781 default | 0.219 | 0.450 | 385 | 2 | 0.223 |

1. The wide level changes almost nothing: of the 110 "related" false pairs on A119b only 10 are Compara pairs at
   eutheria or younger. The rest are older than Eutheria or related only through an HGNC group.
2. Shipped default on track W: **FAIL** (A119b 0.917 -> 0.907; KB3781 passes, 0.210 -> 0.223). Same verdict as row 1654.
3. Step 0f on track W: **STOP.** Best `precision_wide` 0.166 (5 Mb, identity >= 70: 56 correct / 282 false, 75 genes
   newly placed); the member bar is met in every cell (49 to 114 genes) but precision is 0.12 to 0.17 against 0.75.

Track W is closed as a way to rescue edge levers. The pairs these rules add are mostly older than Eutheria.
