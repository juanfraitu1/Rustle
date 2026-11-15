# PREREG r1621: families with negative edges (correlation clustering inside MCL families) (2026-11-10)

User (2026-11-10): try correlation clustering on the families. MCL sees only positive edges; pairs of loci that align on
exons but share too few exons (gate-0 graph edge with SHARE < 0.5) are dropped, not used as evidence against a family.
Written before any of it was computed. Known: the default's scores (A119b damped .918/.782, 153 false pairs of which 98
related; HG002 damped .703/.355), so A119b has at most ~0.08 precision to give.

## Step 0 (stop rule, A119b, default damped partition `r1593/b_T0.5_I2.0.fam.copies.tsv`)
For every scored within-family gene pair (pairscore's coding pairs): is there a NEGATIVE edge between loci of the two genes
(`r1592/g0.tsv` pair with `share.tsv` < 0.5)? GO only if negative edges are on >= 20% of the false pairs AND at >= 2x the
rate of the hit pairs. Otherwise STOP: negative edges carry no signal about false pairs.

## Rule (if GO; `bench/denovo_eval/negsplit.py`)
Inside each default family, agglomerative correlation clustering on its signed graph: positive weight = the edge's weight
when SHARE >= 0.5, negative = -LAMBDA x weight when SHARE < 0.5 (gate-0 weights as dumped, undamped). Start from single
loci; merge the two clusters with the largest net weight while it is > 0. Families are only ever split. Deterministic.
LAMBDA = 0 gives back the default family whenever it is connected: that is the control (algorithm swaps alone moved nothing
before; only the gap to the control counts). Loci without a node in the dumped graph stay with their largest sub-cluster.
Output = the default copy table with `family_id` relabelled; route check: LAMBDA = 0 must reproduce the default's pred and
hit counts within the control's own difference, reported.

## Pricing and bars
Baseline = the default's own damped SHARE curve (I 1.6 / 1.8 / 2.0 / 2.4; `r1593/dev.score.tsv`, `r1593/hg002/hg.score.tsv`),
linear interpolation at the candidate's recall (r1593 pricing).
- Dev A119b: LAMBDA in {0.5, 1, 2, 4}; eligible = recall >= .772; choose the largest gain; go to HG002 only if gain >= +0.02
  and above the control's gain.
- Held out HG002 (once, the chosen LAMBDA, default `r1593/hg002/b_T0.5_I2.0`): PASS if gain >= +0.03 with recall >= .345.
- PASS: proposed as an opt-in refinement after MCL. FAIL or STOP: registered; MCL stays as is.

## Outcome (2026-11-10, `r1621/step0.txt`): STOP at step 0
| within-family gene pairs (A119b default, damped) | direct positive edge | only a negative edge | no edge |
|---|---|---|---|
| true (1,783) | 1,334 (0.748) | 110 (0.062) | 339 (0.190) |
| false (162) | 90 (0.556) | 25 (0.154) | 47 (0.290) |
- Negative edges sit on 15.4% of the false pairs (bar 20%), 2.5x the rate of the true pairs (bar 2x met). **STOP as
  registered**; the split rule was written but not run or scored, HG002 untouched.
- In counts the signal is the wrong way round for a split: 110 true pairs carry only a negative edge against 25 false ones,
  so cutting on negative edges would cost about four true pairs per false pair removed. More than half of the false pairs
  (90 of 162) are joined by a direct positive edge that passes the exon-share rule; those are not a clustering problem.
- The step-0 pair counts (1,783 / 162) differ slightly from pairscore's line for this partition (1,707 hit / 153 false):
  step 0 does not apply pairscore's expressed-universe filter to predicted pairs. Not chased, the verdict does not turn on it.
