# Row 1890: is there a better order for the one-read tier, judged on all twenty datasets? (pre-registered 2026-11-22)

User (2026-11-22, asked about the lost StringTie comparison on Sumatran orangutan testis): "Keep looking". Rows 1883 and 1889
ruled out extension and the seeding pool. What is left inside the list is the ORDER of its second tier. On Sumatran testis the
frozen one-read ranking (`oneread_score`, row 1806) is no better than chance: the best 17,747 one-read chains hold 743 annotated
ones (4.2%), the whole tier 4.1%; on chimpanzee testis the same cut holds 14% against 7.6%. A different order chosen on that one
sample would be fitting it, so this row judges one alternative on all twenty datasets of row 1880 at once.

**Alternative, fixed now:** order tier 2 by the transcript's `cov` attribute (read coverage of the model, which counts the
partial reads lying on it), most first; ties by `oneread_score`. No fit, no parameter. Nothing else in the list changes.

## Measure (existing runs of row 1880, no new run; `r1890_tier2_order/order.py`)
For each dataset and each rival larger than `final.gtf`: correct chains of the list at the rival's size, under the current order
and under `cov`. Also the best 10% of the tier (the `list_sensitive` cut) under both.

## Reading, fixed now
- **ADOPTABLE** = at StringTie's size the `cov` order gives at least as many correct chains as the current order on at least 17 of
  20 datasets, on at least 6 of the 7 non-primate ones, loses no AHEAD of row 1880, and turns Sumatran testis to AHEAD.
- **BETTER BUT NOT ENOUGH** = the first three hold and Sumatran testis stays BEHIND: reported; changing the order is then the
  user's call on its general merit alone.
- Anything else: the current order stays, and the Sumatran loss stands; no further order is tried (one alternative was
  registered, not a search).

## Result (2026-11-22): the coverage order is worse everywhere; the current order stays, the Sumatran loss stands

`r1890_tier2_order/order.NAME.txt`, existing runs of row 1880. At StringTie's size (17 datasets where StringTie is larger than
`final.gtf`; on fly, mouse and zebrafish simulated it is smaller, so tier 2 is not reached):
- The `cov` order gives fewer correct chains than `oneread_score` on 17 of 17, by 23 to 1,189 chains (mouse ES -230, heart -26,
  WTC11 -538, chimpanzee testis -757, bonobo fibroblast -1,189).
- Sumatran testis -686 -> -1,070; Bornean testis would turn from AHEAD (+666) to BEHIND (-506).
- **Verdict: neither ADOPTABLE nor BETTER.** The frozen ranking of row 1806 is the better order on every real library,
  including the one where it looked like chance. As registered, no further order is tried.

**State of "keep looking" on the Sumatran sample.** Three registered attempts, none closes it against RefSeq: transfrag
extension (row 1882: closes at StringTie's size, fails complete truth, mouse and human in row 1883), the seeding pool of `all`
with the repair (row 1889: -686 -> -582), the order of the one-read tier (this row). What the sample needs is chains longer
than its reads, which is what the complete truths say Rustle should not build.
