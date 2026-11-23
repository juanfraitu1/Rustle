# PREREG: alignment scores inside the EM of the `o2` labels (row 1921, 2026-11-23)

The user asked to try "the EM over alignment scores". What the `o2` stage runs today is a counts-only EM
(`tools/eq_classes.py em`): the scores decide only which loci a read is EXACTLY tied between (minimap2 `ms`); a read
whose best score is at one locus alone is placed there whatever the margin, and tied reads are shared by abundance.
The scores themselves are not in the likelihood. This row tests the variant where they are.

## The variant (`bench/o2_families/em_filter.py --weighted T`, opt-in, written before any run)
Every alignment within 10 `ms` of the read's best is a candidate with weight exp(-gap / T); the EM shares the read by
abundance x weight; a read with one candidate stays placed. The `expressed` label keeps its meaning (a read whose best
score is at this locus alone), so only the EM count, and with it the `homology_only` verdict (EM < 2), can change.
T -> 0 is today's test (checked on one planted set before registration: the labels are the same on 135 of 136 loci and
the EM counts within 0.2; no score of the weighted arms was looked at beyond the one summary line each run prints).
Two values, fixed now: T = 1 (a substitution, about 5 `ms`, is a 150-fold likelihood ratio: the scale of HiFi error
rates) and T = 5 (a soft arm, one substitution = factor 2.7).

Expected mechanism: a read of an expressed copy with a sequencing error at a site where a silent relative differs gets
its best score at the silent copy and is placed there today; the weighted EM lets the abundant copy claim it back. The
risk is the reverse: a weakly expressed copy next to an abundant relative loses its own reads.

## Data (all already planted with ragged ends; reuse is the user's rule, a fresh seed follows if the claims hold)
Human named families seeds 1920-1924, scenarios E and S (`r1918_rerun_ragged/lab`); Y-ampliconic set E, D, S
(`r1918_rerun_ragged/yag`; D = one copy 10x, the others a third: the risk case); gorilla own families E, S
(`r1917_ragged_reads/g`). Input = the same `rg.multicopy.gtf`, BAM and score table as the default labels; scorer
`label_score.py` (units = copies with identical planted isoforms).

## Claims (T = 1 against the default labels)
- W1: over the five human S seeds, silent units labelled `expressed` or `shared` (summed) fall.
- W2: over all sets, expressed units labelled `homology_only` (summed) do not rise by more than 1.
- W3: gorilla S: silent units labelled `expressed` or `shared` do not rise.
W1 + W2 + W3 = confirm on a fresh seed (1927), then propose as the default of the labels; the default is the user's call.
T = 5 is reported beside it and is not a candidate unless it meets the same three.

## Predictions (mine)
Small effect. W2 and W3 hold. W1 is uncertain: most silent units that are not `homology_only` today are `shared`
(exact ties over the whole read, which no weight changes), so I expect at most a handful of units to move, and W1 may
hold only by one or two units or not at all. T = 5 costs expressed units in scenario D.

## Result on the registered sets (job 4496780, 2026-11-23)
Summed over the sets (default / T = 1 / T = 5):

| | silent units labelled `expressed` or `shared` | of which `expressed` | expressed units labelled `homology_only` |
|---|---|---|---|
| human, five seeds (E + S) | 23 / 21 / 20 | 6 / 4 / 3 | 11 / 11 / 92 |
| Y-ampliconic set (E, D, S) | 3 / 3 / 3 | 1 / 1 / 1 | 2 / 2 / 29 |
| gorilla (E, S) | 3 / 3 / 3 | 1 / 1 / 0 | 5 / 5 / 13 |

- W1 HOLDS by two units (23 -> 21: seeds 1921 and 1922 each lose one silent unit labelled `expressed`). W2 HOLDS (no
  expressed unit changes label anywhere, scenario D included). W3 HOLDS (3 = 3).
- T = 5 is not a candidate: it demotes 8-15 expressed units per set to `homology_only` (92 against 11 on the human
  seeds), as predicted.
- Prediction met: a handful of units, and only on the `expressed`-labelled silent ones; the `shared` ones do not move.
- By the registered rule the next step is the confirmation on unseen seeds. Amendment written before that run: three
  fresh seeds (1927, 1928, 1929; scenario S, human named families) instead of one, because the effect is two units in
  five seeds and one seed cannot show it. Confirmed = T = 1 has no more silent units labelled `expressed` or `shared`
  than the default on any of the three, fewer in the sum, and no more expressed units labelled `homology_only` in the sum.

## Confirmation on three unseen seeds (job 4496782, scenario S, human named families)

| seed | silent units labelled `expressed` or `shared`: default / T = 1 / T = 5 | expressed units labelled `homology_only` |
|---|---|---|
| 1927 | 2 / 1 / 1 | 1 / 1 / 6 |
| 1928 | 6 / 5 / 5 | 0 / 0 / 3 |
| 1929 | 4 / 3 / 3 | 1 / 1 / 6 |
| sum | 12 / 9 / 9 | 2 / 2 / 15 |

**CONFIRMED for T = 1:** no worse on any seed, three fewer in the sum (one false `expressed` label removed per seed), no
expressed unit changes label. Over all eight human seeds: 35 -> 30 silent units mislabelled, 0 expressed units lost.
T = 5 again costs expressed units (15 against 2). Size: about one unit per run; the units left are `shared` (exact ties
over the read) and a few `expressed` ones the weights do not reach. Not run on real samples (no truth for the labels
there). The flag is opt-in (`em_filter.py --weighted 1`); making it the default of stage `o2` is the user's call.

**Default since 2026-11-23 (user: "yes, lets make it default").** Stage `o2` calls `em_filter.py --weighted 1`;
`RUSTLE_O2_EM_WEIGHT=0` gives the counts-only EM of before, any other value is T. The benchmark runner follows
(`O2_EM_WEIGHT`, default 1). Planted results quoted before this date used the counts-only EM. Not validated on real
samples (no truth for the labels there); on real gorilla OR6737 the two EMs are compared in workspace `r1921_weighted_em/real`.

## First look on a real sample, and an amendment (gorilla OR6737, driver stage `o2` end to end, job 4496783; no truth)
Both settings ran through the driver (exit 0; 20,120 loci). Labels, counts-only -> weighted: expressed 18,868 -> 18,844,
shared 46 -> 46, homology_only 1,206 -> 1,230. 24 loci go from `expressed` to `homology_only` (12 of them in a multi-copy
family), 2 from `shared` to `homology_only`, 2 the other way. All 24 are thin: 1-4 reads with their best score there
(21 of them 1 or 2 reads).
**Artifact found there and fixed (amendment after the registered runs):** 6 of the 24 had an EM count that prints as
2.00: two reads each leaking about 1e-4 of a read to a rival ten `ms` away fell under the cut of exactly 2. The cut is
now taken at the precision the table prints (`round(em, 2) >= --min-em`). Every planted set of this row (the 15
registered ones and the 3 confirmation seeds) was rerun with the fix: all 18 score identically (job 4496785), so the
results above stand. With the fix 18 loci change on OR6737.
**Known limit, not tuned further:** a locus with exactly two reads still changes label when a relative of similar
abundance sits one substitution away (EM 1.95-1.99 in 5 of the 18). Planted reads do not show this (no expressed unit
lost in 18 sets) because planted copies have more reads; real libraries have many two-read loci. Whether those 18 are
errors corrected or thin real loci demoted is not known: there is no truth. `RUSTLE_O2_EM_WEIGHT=0` is the way back.
