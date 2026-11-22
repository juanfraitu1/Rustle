# Row 1894: StringTie's selection as an order of the one-read tier (locus share), twenty datasets (pre-registered 2026-11-22)

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
