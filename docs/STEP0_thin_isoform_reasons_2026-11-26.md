# Step 0: why are annotated isoforms with reads of their own not assembled? (row 1966, 2026-11-26)

Row 1965: on a gorilla sample about 10,000 annotated isoforms have 1-4 full-length reads of their own (primary, mapping
quality >= 1, exactly the annotated intron chain) and are not in the assembly; relatives' reads reach 0.2% of them. The
user: start exploring the recovery of thin isoforms. Before any rule: what removes them today. A description, no change.

## Known before looking (memory `assembly`, do not rediscover)
- The read floor is 2 on distinct coordinate records, by the user's choice (precision, clean SIRV set). An isoform with
  one read is below it by design and lives in the one-read tier of `all` (`PREFIX.oneread.gtf`, `PREFIX.ranked.gtf`).
- Closed, not to be re-offered without a new mechanism: fixed or depth-scaled read floors (row 1852), sub-chain collapse
  variants (rows 1839-1843; `RUSTLE_POLISH_ISM_RATIO` is the user's dial), one-read junction / end proofs (rows 1805,
  1806), retained-intron relaxation, chain-score refits.
- Gorilla RefSeq holds computed models; here every isoform counted has a read with exactly its chain, so the usual
  "longer than any read" bias does not apply to this count.

## Method (`bench/o2_truth/thin_reasons.py`)
Today's assembly rebuilt with the dropped-chain record (`pool_bam.py as:0.98` + `--seed-pool all` = today's pool, row
1942; `PREFIX.polish_reasons.tsv` names the rule that removed each candidate). OR6737 and KB3781.
For every annotated isoform with own reads that is not assembled, by own-read count (1, 2, 3, 4, 5-9, 10+):
`removed by RULE` (a candidate with exactly the chain was built and that rule removed it), or `no candidate` (none was
built). For `no candidate`: the distinct coordinate records among its own reads (the floor counts those).
## Predictions (mine)
Most of the missed thin isoforms have exactly 1 own read (the floor). Of those with 2-4, the largest classes are the
sub-chain collapse and the isoform fraction; a visible share has 2 reads that are 1 distinct record.
What would be worth a test afterwards: only a class that is neither the floor nor a closed dial.

## Result (jobs 4510059 / 4510060; the rebuilt assembly equals the stored one: KB3781 61,772 of 61,772 chains, OR6737 73,128 of 73,173)
Annotated isoforms with own full-length reads that are not assembled:

| own reads | OR6737: not assembled | of them: no candidate / isoform fraction / sub-chain collapse / other rules | KB3781: not assembled | no candidate / isoform fraction / sub-chain collapse / other rules |
|---|---|---|---|---|
| 1 | 8,813 | 8,800 / 10 / 0 / 3 | 6,578 | 6,564 / 12 / 0 / 2 |
| 2 | 1,169 | 176 / 824 / 105 / 64 | 1,433 | 159 / 1,066 / 114 / 94 |
| 3 | 535 | 17 / 412 / 56 / 50 | 776 | 20 / 615 / 85 / 56 |
| 4 | 342 | 8 / 252 / 43 / 39 | 487 | 4 / 404 / 45 / 34 |
| 5-9 | 601 | 12 / 451 / 72 / 66 | 1,098 | 12 / 858 / 115 / 113 |
| 10 or more | 404 | 45 / 228 / 77 / 54 | 1,199 | 57 / 859 / 143 / 140 |
| all | 11,864 | | 11,571 | |
- **One read: 74% (OR6737) and 57% (KB3781) of everything missed.** No candidate is built: the read floor of 2, as designed.
- **With two or more reads the isoform fraction removes 71% and 76%** of what is missed (2,167 and 3,802 isoforms); the
  sub-chain collapse 12% and 10%.
- **No candidate with two reads:** 150 and 145 of 176 and 159 are two reads with identical coordinates (one distinct record).
- **No candidate with five or more reads in two or more records:** 56 and 68 isoforms; 48 and 62 of them have a
  non-canonical intron in the annotation (the strict junction gate, by design); 15 and 18 are assembled with a junction
  within 10 bases (junction snap).
What each rule removes, and how much of it is annotated (spliced candidates):

| rule | OR6737: removed / annotated | KB3781: removed / annotated |
|---|---|---|
| sub-chain collapse | 52,780 / 0.7% | 64,525 / 0.8% |
| isoform fraction | 26,782 / 8.1% | 53,309 / 7.2% |
| retained intron | 2,375 / 4.1% | 3,528 / 2.9% |
| spliced share | 1,134 / 6.3% | 3,364 / 5.2% |
Isoform fraction by the candidate's reads (OR6737): 2 reads 13,613 removed, 6.1% annotated; 3-4: 7,819, 8.6%; 5-9: 3,898,
11.6%; 10-49: 1,408, 15.5%; 50 or more: 44, 25%. KB3781: 4.7%, 6.6%, 9.1%, 14.0%, 22.6% (499 at 50 or more).

## Reading
- My predictions: "most have one read" right; "sub-chain collapse and isoform fraction dominate the rest" half right
  (the fraction alone is three quarters; the collapse a tenth).
- **Nothing here is a defect.** Every class is a rule set on purpose: the floor (user's choice, row 1850), the isoform
  fraction (rows 859, 1083, 1423, 1837, 1849, 1876: a trade, 5-26% of what it removes is right; here 7-8%), the junction
  gate. The pool of "thin isoforms" is large because these rules are strict, and at every read count more than 75% of
  what the fraction removes is not annotated.
- Both classes already have a home outside the main list: one-read chains in the ranked list's second tier
  (`PREFIX.oneread.ranked.gtf`), rule-dropped chains in its opt-in third tier (`RUSTLE_RANKED_TAIL=1`, row 1891).
- So a recovery of thin isoforms is a question of RANKING inside those tiers (which of the 8,813 one-read and 2,167
  fraction-removed annotated isoforms come first), not of a new assembly rule. What is measurable next without
  reopening a closed dial: how many of them the existing tiers place near the top on this sample.

## Step 1 (description, registered before the run; the user agreed 2026-11-26): where the existing tiers put them
OR6737 and KB3781: stage `oneread` with `RUSTLE_RANKED_TAIL=1` on the rebuilt assembly gives `PREFIX.ranked.gtf` (tier 1
= the final list by reads, tier 2 = one-read chains by the shipped score `r1899`, tier 3 = rule-dropped chains by
reads). `bench/o2_truth/ranked_thin.py`: per tier, its size and annotated share by tenth of the tier's order, and how
many of the annotated isoforms with own reads that the main assembly lacks sit in each tenth.
Predictions (mine): tier 2 holds at least 90% of the 8,813 / 6,578 one-read annotated isoforms, and its first tenth
holds about 30% of them at about 35% annotated; tier 3 holds most of the fraction-removed ones at 5-10% annotated with
little concentration at the top. No claim, no default involved.

### Step 1 result (jobs 4510073 / 4510074)
| | OR6737 | KB3781 |
|---|---|---|
| tier 1 (main list): spliced chains / annotated | 58,557 / 41.7% | 53,509 / 43.8% |
| tier 2 (one-read chains): chains / annotated / annotated ones the main assembly lacks | 123,751 / 5.0% / 6,044 | 68,558 / 4.9% / 3,338 |
| tier 2, first tenth: chains / annotated / share of the tier's new annotated | 12,375 / 11.5% / 23.0% | 6,855 / 13.8% / 28.1% |
| tier 2, first half: share of the tier's new annotated | 76.3% | 82.3% |
| tier 3 (rule-dropped chains): chains / annotated / annotated ones the main assembly lacks | 262,580 / 2.2% / 5,646 | 357,682 / 2.2% / 7,885 |
| tier 3, first tenth: annotated / share of the tier's new annotated | 3.1% / 13.5% | 4.9% / 21.9% |
- The tiers DO hold the missed isoforms: 6,044 + 5,646 = 11,690 on OR6737 against 11,864 missed with own reads; 11,223
  against 11,571 on KB3781.
- They hold them diluted: one annotated isoform in 20 chains of tier 2 and one in 45 of tier 3, against one in 2.4 of
  the main list. The shipped score lifts the first tenth of tier 2 to 11.5-13.8% (2.3-2.8 times the tier's average);
  the read order of tier 3 hardly concentrates anything.
- My predictions were too optimistic: I expected the first tenth of tier 2 near 35% annotated; it is 12-14%.
- So on these samples: taking the first tenth of the one-read tier adds 1,389 / 937 annotated isoforms (+5.6% / +4.0%
  of the main assembly's) together with about 11,000 / 5,900 chains that are not annotated.

## Conclusion of the exploration
The thin isoforms are not lost: they are in the labelled tiers. They cannot be lifted into the main list by any rule
found here, because at one or two reads an annotated isoform is outnumbered 7 to 45 times by chains that look the same.
The levers that would separate them have been tried and are closed (one-read proofs and scores, rows 1805, 1806, 1890,
1895, 1899; isoform fraction, rows 1837, 1876), and evidence from a second individual is excluded by the standing rule
that no rule may need a second sample. Relatives' reads (row 1965) reach 0.2% of them. Nothing new to offer.

## Step 2 (description, registered before the run; user 2026-11-26: could clustering or StringTie's nascent feature help)
Clustering (IsoCon, Strainy) is not run: it separates sequences and needs depth, the thin class has one or two reads,
and its ports are closed (rows 1375-1377, 1463, 1823; IsoCon and `isoseq cluster2` as assemblers).
The nascent idea, as a check of the tiers only: `RUSTLE_POLISH_NASCENT=tag` (register 1635, StringTie's `mark_if_nascent`:
a transcript whose introns are the first introns of a better-supported one from the 5' end and which runs on into its
next intron) on the two assemblies the tiers are made from (one-read floor; one-read floor with the dropping rules off).
`bench/o2_truth/ranked_thin.py --tags`: per tier, the tagged share among the annotated chains the main assembly lacks and
among the unannotated ones, and what the first tenth looks like once tagged chains are put last.
Prediction (mine): the tag covers under 10% of the unannotated chains of either tier (the retained-intron rule and the
sub-chain collapse already take the common forms) and about 1% of the annotated ones; no useful reordering.

### Step 2 result (jobs 4510109 / 4510110)
| | OR6737 | KB3781 |
|---|---|---|
| tier 2: tagged among the annotated new / among the unannotated | 31 of 6,044 (0.5%) / 1,356 of 117,577 (1.2%) | 19 of 3,338 (0.6%) / 921 of 65,193 (1.4%) |
| tier 2, first tenth annotated: as shipped / tagged chains last | 11.5% / 11.5% | 13.8% / 13.8% |
| tier 3: tagged among the annotated new / among the unannotated | 307 of 5,646 (5.4%) / 21,851 of 256,798 (8.5%) | 488 of 7,885 (6.2%) / 37,408 of 349,692 (10.7%) |
| tier 3, first tenth annotated: as shipped / tagged chains last | 3.1% / 3.0% | 4.9% / 5.1% |
As predicted: the tag covers 1% of the one-read tier and 9-11% of the rule-dropped tier, and it marks annotated chains
at more than half the rate of unannotated ones (5-6% against 9-11%), so it separates nothing. Moving tagged chains last
changes the first tenth by 0.2 points at most. The nascent idea is closed for the tiers; clustering was not run (see above).
