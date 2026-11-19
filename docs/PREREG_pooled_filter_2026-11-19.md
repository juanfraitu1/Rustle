# Row 1818: a chain score fitted on primate and non-primate data together, judged leave-one-species-out (pre-registered 2026-11-19)

User (2026-11-19): a quick study on wrong chains to find a filter that works across the board. Rows 1814 / 1815: single read
features give one rule (5' truncation). Row 1813: the stable score r1813 was fitted on 11 primate samples and filters little on
other data. The untried step: the same score fitted with non-primate data in the fit. Registered before any fit. Nothing here
changes a default.

**Data.** Feature tables of `PREFIX.gtf` at today's defaults for nine datasets in six species groups: fly simulated; zebrafish
simulated; mouse (simulated, ES, heart, adrenal); gorilla OR6737; human A119b; chimpanzee testis (tables of rows 1813-1815).
**Model.** Design `frozen9` exactly as r1813 (same 23 features and extras, pooled mean / sd / deciles stored, same fit routine).
Each dataset enters a fit as at most 25,000 chains drawn at random (seed 1818), so the deep libraries do not outweigh the rest.
**Judging.** Leave one species group out: fit on the other groups, score the held-out datasets. On each held-out dataset the new
score keeps its best chains until it has as many correct chains as today's filter keeps (r1813 at 0.018 + 5' truncation rule);
compared: precision of the kept chains.

**Bar.** ACROSS THE BOARD = precision at equal correct chains is higher than today's filter on all nine held-out datasets, and by
at least 0.01 on at least six. Otherwise the datasets where it is not higher are named, and whether they are the simulated ones,
the real non-primate ones or the apes.
**If it passes:** one model fitted on all nine is written and put to the user; a default needs a check on a species in none of
these tables.

**Limits.** Labels differ in quality (complete truth on simulated reads, GENCODE on mouse, thinner RefSeq on the apes), so the fit
mixes three meanings of "wrong". Four of nine datasets are mouse. Leave-one-species-out for mouse removes all four at once.
Equal correct chains compares one point of the curve.

## Result (2026-11-19): NOT ACROSS THE BOARD

Job 4482948, table `r1818_pooled_filter/loso.txt`. Precision of the kept chains at the number of correct chains today's filter keeps:

| held-out dataset | today's filter | pooled model, species held out | difference |
|---|---|---|---|
| fly simulated | 0.944 | 0.935 | -0.009 |
| zebrafish simulated | 0.845 | 0.841 | -0.003 |
| mouse simulated | 0.924 | 0.930 | +0.006 |
| mouse ES | 0.344 | 0.346 | +0.002 |
| mouse heart | 0.565 | 0.562 | -0.002 |
| mouse adrenal | 0.560 | 0.560 | 0.000 |
| OR6737 | 0.364 | 0.349 | -0.015 |
| A119b | 0.212 | 0.209 | -0.004 |
| chimpanzee testis | 0.397 | 0.381 | -0.016 |

Higher on 2 of 9, by 0.01 or more on none. Putting non-primate data into the fit does not make the score better on non-primate
data, and a fit without the held-out ape is worse on that ape than the primate-fitted r1813 (which saw it). No model written.

**Shown without a bar: what is still wrong in `PREFIX.final.gtf` where the truth is complete** (`left.txt`; runs of row 1817):

| class of wrong chain | fly simulated (261) | mouse simulated (1,327) | zebrafish simulated (1,514) | SIRV (28) |
|---|---|---|---|---|
| piece of a true chain | 70 | 525 | 547 | 0 |
| new intron within 10 bp of a true one | 65 | 361 | 279 | 22 |
| new intron elsewhere, some true introns | 50 | 237 | 287 | 1 |
| no true intron | 58 | 79 | 259 | 4 |
| true introns, combination not in truth | 18 | 125 | 142 | 1 |

More than half of each class has over 5 reads. The near-miss test the features already hold (`wobble`: an intron within 10 bp of
one with 5 times its reads; `wobble.txt`) removes 3% to 14% of the wrong chains but costs 0.6% to 2.9% of the correct ones
(A119b 983 correct for 5,887 wrong), outside the 1% used in rows 1814 / 1815.

**Reading.** With these 23 read features the limit is reached: one rule (row 1815), one stable score (row 1813), and neither a
pooled fit nor a second single rule adds to them. What is left needs evidence the features do not hold.
