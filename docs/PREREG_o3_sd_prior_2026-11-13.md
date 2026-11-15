# Row 1684 — step 0: duplicated sequence on the draft as a prior for absent-copy piles (O3) (2026-11-13)

Goal (user, 2026-11-13): "explore what other ways can we include the segmental duplications so they help the
pipeline"; user's choice "O3 prior from SD regions". O3 measurement, opt-in only, nothing wired, not for
`advisor_share` (the user has said O3 is not yet integrated).

**Background.** Row 1539's blind divergence-pile detector on the Kamilah draft (KB3781 panel reads) found 24 piles:
11 true (reads of a copy absent from the draft), 13 false, all 13 with no labelled read. Precision 0.46 against a 0.5 bar.

**No SD map of the draft exists** (SEDEF covers CHM13, HG002 paternal and the gorilla T2T only; a full build costs
17-40 h). The step 0 therefore asks the SD question locally: each pile region is aligned to the whole draft
(`minimap2 -x asm20 -c -N 50 -p 0.5`).

**Rule (frozen).** A pile is "duplicated on the draft" when it has a hit outside its own interval with an aligned
block >= 1 kb at identity >= 0.90 (the usual SD definition) whose query segment is < 50% soft-masked.
**Stated before the run:** theory cuts both ways. A copy absent from the draft because the draft collapsed a recent
duplication leaves an absorbing locus that is single-copy on the draft, so true piles may be the ones WITHOUT a
duplicate. Both directions are read.

**Bars (the better of the two groups, duplicated / not duplicated):** precision >= 0.60 within the group AND >= 9 of
the 11 true piles in it; a group holding <= 2 or >= 22 of the 24 piles is vacuous. Otherwise FAIL.
**Limits stated now:** 24 piles, one sample, labels cover panel reads only; a pass would be a lead for a registered
test on another draft, not a rule.

## Result (2026-11-13, job 4470757, workspace `r1684_o3sd`) — FAIL

| group | piles | true | false | precision |
|---|---|---|---|---|
| duplicated on the draft | 9 | 6 | 3 | 0.667 |
| not duplicated | 15 | 5 | 10 | 0.333 |

Neither group meets the bar (>= 9 of the 11 true piles at precision >= 0.60): the duplicated group has the precision
but holds only 6 true piles. The direction favours the prior (0.667 against 0.333) but on 24 piles that difference is
within chance (Fisher exact, two-sided, about 0.2), and using it as a filter would lose 5 of 11 true piles. Both
readings stated before the run occur: 6 true piles sit on sequence the draft holds twice, 5 on sequence it holds once
(the collapsed case). Row closed as a rule; as a ranking hint it would need another draft and more piles.
Per-pile table: `r1684_o3sd/run.4470757.out`.
