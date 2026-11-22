# Row 1895: StringTie's selection as an order of the one-read tier (locus share), twenty datasets (pre-registered 2026-11-22)

User (2026-11-22): try StringTie's coverage-based strategies. Row 1893: in long-read mode StringTie does not trim and builds
almost nothing; it SELECTS among read chains, keeping a low-count chain when it is a large enough share of its locus (isoform
fraction 0.01) at a coverage of one read. About half of its lead on Sumatran orangutan testis is that selection. Rustle's
counterpart is the order of the one-read tier of the ranked list. Row 1892 called further variants on that sample closed; this
row re-opens the order once, on the user's request, with a rule taken from the rival rather than searched for.

**Order, fixed now.** For a one-read chain, TOP = the largest `reads` among the chains of `final.gtf` that overlap its span on
the same strand of the contig (0 when none overlaps). Tier 2 is ordered by 1 / TOP, largest first: chains in loci with nothing
assembled first, then chains beside weakly expressed transcripts, last the ones beside the most expressed. Ties by
`oneread_score`. Nothing else changes. No parameter, no fit.

## Measure (existing runs of row 1880, no new run; `r1894_locus_share_order/share.py`)
Per dataset and rival larger than `final.gtf`: correct chains of the list at the rival's size under the current order and
under the locus-share order; the best 10% of the tier under both.

## Reading, fixed now (the bars of row 1890)
- **ADOPTABLE** = at StringTie's size the new order gives at least as many correct chains as the current one on at least 15 of
  the 17 datasets where StringTie is larger than `final.gtf`, on at least 3 of the 4 non-primate ones among them, loses no
  AHEAD of row 1880, and turns Sumatran testis to AHEAD.
- **BETTER, NOT ENOUGH** = the first three hold, Sumatran testis stays BEHIND: a general improvement of the order, the user's
  call on that merit alone.
- Otherwise the frozen order stays. This is the second order tried (row 1890 the first); a third is not registered.

(Registered as row 1894; another session took that number at the same time, so this is row 1895. Workspace
`r1894_locus_share_order`, alias `r1895_locus_share_order`. No change of content.)

## Result (2026-11-22): better on the four mouse libraries, worse on all thirteen primate ones; the frozen order stays

`r1894_locus_share_order/s.NAME.txt`, existing runs of row 1880. Correct chains of the list at StringTie's size, current order
-> locus-share order (17 datasets where StringTie is larger than `final.gtf`):
- **Non-primate, 4 of 4 better:** mouse ES +1,003 -> +1,048 ahead; heart +544 -> +727; adrenal +647 -> +782; cortex +1,272 -> +1,496.
- **Primate, 0 of 13 better:** from -9 (Sumatran testis) to -945 (Bornean testis); WTC11 -230, HG002 -450, A119b -697,
  chimpanzee testis -673, bonobo fibroblast -874.
- **Sumatran testis:** -686 -> -695. **Bornean testis turns from AHEAD (+666) to BEHIND (-279).**
- **Verdict: not ADOPTABLE, not BETTER** (4 of 17, an AHEAD lost). The frozen `oneread_score` order stays; no third order is
  registered.

**Reading (mine).** StringTie's selection rule, moved onto Rustle's one-read chains, helps exactly where the libraries are
shallow and well annotated (mouse tissue, 4-6 M reads) and hurts on every primate library, human ones included. So it is not
what gives StringTie its lead on Sumatran testis either: there the two orders are 9 chains apart. Together with row 1893 this
exhausts the StringTie strategies that can be named: trimming is inert in long-read mode, the built chains are wrong on complete
truth, and its selection rule does not transfer to Rustle's one-read tier. The split by species group is the reverse of a
primate fit and is worth knowing: a per-library choice between the two orders is NOT proposed (no reference-free statistic
has separated library types in this project, rows 1834 / 1883).
