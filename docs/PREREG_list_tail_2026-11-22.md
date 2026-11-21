# Row 1884: a third tier at the end of the ranked list, against FLAIR on HG002 (pre-registered 2026-11-22)

Row 1880: the ranked list (final.gtf by reads, then the one-read tier by score) is behind FLAIR on one dataset of fifteen,
HG002: FLAIR holds 38,489 correct chains in 728,951; Rustle's whole list holds 32,318 in 318,884 and simply ends. FLAIR has no
isoform fraction, no sub-chain collapse and no retained-intron rule; the chains those rules remove are in neither tier of the
list. This row appends them.

**Tier 3, fixed now:** the chains of an assembly run with a one-read floor and the four dropping rules off
(`RUSTLE_MIN_CHAIN_READS=1 RUSTLE_GATE_MIN_READS=1 RUSTLE_POLISH_ISM_RATIO=0 RUSTLE_POLISH_RETAINED_RATIO=0
RUSTLE_POLISH_SPLICED_SHARE=0`, `--polish-isoform-fraction 0`; the `f1all` arm of row 1864 plus the gate bypass), its `x.gtf`,
that are in neither earlier tier, by `reads`, most first. It is appended AFTER tier 2, so every comparison of row 1880 at a
size inside the first two tiers is unchanged by construction; only rivals longer than the list can move.

## Runs
HG002 (the miss) and A119b (the other dataset where FLAIR is longer than the list by more than 50,000 chains and human).
Driver snapshot and binaries of row 1880. Workspace `r1884_list_tail`; scorer = `rank.py` of row 1880 with the third tier.

## Reading, fixed now
- **CLOSED** = with tier 3 the list cut at FLAIR's size has at least FLAIR's correct chains on HG002, and A119b stays AHEAD.
- Otherwise the miss stands; the table gives how far tier 3 goes.
- Reported: tier 3's size and precision. It is a recall tail (expected precision under 0.1, as rows 1837 / 1864 measured for
  these chains); it is not proposed for any default output. If CLOSED, the candidate "ranked list as a product" of row 1880
  gains a third, clearly labelled tier, the user's call like the rest.

## Limits
- Tier 3 costs a third assembly. Its chains are the ones four registered rules remove on purpose.
- Two datasets, both human; the other thirteen FLAIR comparisons are already AHEAD without it.
