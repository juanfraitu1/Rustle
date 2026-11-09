# PREREG row 1593: exon sharing (SHARE >= 0.5) at low inflation, against baselines that reach its recall; held out on HG002

Written 2026-11-09, after row 1592 (STOP on both tracks; post hoc its T 0.5 I 2.0 arms, .869/.783 plain and .918/.782
damped, lay above both baseline curves' recall ranges) and before any run of this row. **Dev is partly seen:** the
A119b T 0.5 I 2.0 arms are row 1592's, known when this was written. HG002 is the decisive test and is unseen.

## Rule (unchanged from row 1592)

Gate-0 graph (`mcl_families --min-shared-exon-frac 0.0 --dump-graph`), edges kept when the locus pair's SHARE >= 0.5
(`exonshare.py --locus-scores` on that run's `fam.loci.gated.paf`), fed back with `--extra-edges-only` at gate 0.0;
track (b) adds `--sd-damp-bed BED --sd-damp 0.25`. HG002: the route check of row 1592 amendment 1 (g60 fed back gives
the g60 partition and pairscore line) runs before any HG002 arm is scored.

Candidates per track: I in {1.6, 1.8, 2.0, 2.4}.

## Baselines (extended so they can reach the candidates' recall)

- (a) plain inflation curve I {1.2, 1.4, 1.6, 1.8, 2.0, 2.4, 2.8, 3.2, 4.0} (A119b: `r1571/clu`, g60 for 2.8; HG002:
  `r1580/hg002/plain_I*` plus new I 1.2, 1.4).
- (b) damped curve (`--sd-damp 0.25`) I {1.2, 1.4, 1.6, 2.0, 2.4, 2.8, 3.2, 3.6} (A119b: r1578span / r1581dev plus new
  I 1.2, 1.4; HG002: `r1580/hg002/sd_I*` plus new I 1.2, 1.4).

## Pricing (fixed now)

Gain of a candidate at recall r over a baseline curve (points sorted by recall; exact hit counts):
- r inside the curve's recall range: precision minus the curve's precision at r (linear interpolation between the two
  bracketing points; ties in recall: the higher precision).
- r above the curve's highest recall: precision minus the precision of the curve's highest-recall point (the candidate
  then has more recall than every baseline point; it must beat that point's precision).
- r below the curve's lowest recall: not priced, FAIL.

## Dev (A119b) and held-out (HG002)

Dev: per track, eligible candidates have recall >= the track's shipped default (a .740, b .718); choose the largest
gain. A track goes to HG002 only if that gain is >= +0.04. The chosen I is appended here before any HG002 arm is scored
(HG002 arms may be computed in parallel; pairscore runs on them only after the append).

**HG002 PASS** (per track, on the dev-chosen I): gain >= +0.03 by the pricing above against the HG002 curve of that
track, AND recall >= (a) .305 (plain default .345 - 0.04) / (b) .321 (damped default). Reported, no bar: all four I.

## What follows

- A track PASSES: an exon-sharing gate in `mcl_families` (port byte-identical on membership to this route) proposed
  opt-in at the chosen I; defaults are the user's call.
- Both FAIL: register; exon sharing stays an edge-level finding (row 1591).

## Setting chosen (2026-11-09, A119b dev, job 4453933; appended before any HG002 arm is scored)

Damped baseline extended: I 1.2 .659/.776, I 1.4 .862/.741 (plain I 1.2 .472/.771, I 1.4 .627/.781). `r1593/select.py`
on exact counts (truth 2,183):

| arm | precision / recall | priced | gain |
|---|---|---|---|
| a I1.6 / 1.8 / **2.0** / 2.4 | .826/.799, .862/.786, **.869/.783**, .878/.732 | above, above, above, below floor | +.199, +.235, **+.241**, (+.059 ineligible) |
| b I1.6 / 1.8 / **2.0** / 2.4 | .902/.788, .910/.783, **.918/.782**, .919/.730 | above, above, above, inside | +.243, +.251, **+.259**, +.022 |

**Chosen: I 2.0 on both tracks** (gains >= +0.04, recall >= floors). Both gains are priced against the baseline's
highest-recall point (plain I 1.4 .627/.781; damped I 1.2 .659/.776). HG002 route check (job 4453934): identical
partition (644 families) and pairscore line (.466/.345 = `r1580/hg002/plain_I2.8`); 0 new nodes in every arm.

## Outcome (2026-11-09, HG002 held-out, job 4453934, scored after the append above): PASS on both tracks

Truth 414 pairs. HG002 curves extended (plain I 1.2 .122/.498, 1.4 .233/.483; damped I 1.2 .246/.507, 1.4 .473/.374).

| verdict arm | pred / hit / false | precision / recall | baseline at recall | gain (bar +0.03) | recall bar |
|---|---|---|---|---|---|
| **(a) plain, T0.5 I2.0** | 304 / 164 / 140 | **.539 / .396** | plain .223 (inside) | **+0.317** | .396 >= .305 |
| **(b) damped, T0.5 I2.0** | 209 / 147 / 62 | **.703 / .355** | damped .558 (inside) | **+0.145** | .355 >= .321 |

**(a) PASS, (b) PASS.** Against the shipped defaults: plain .466/.345 (307 / 143 / 164) -> (a) +21 true pairs and 24
fewer false; damped .719/.321 (185 / 133 / 52) -> (b) +14 true pairs for 10 more false. Reported, no bar: (a) I 1.6 / 1.8
/ 2.4 .402/.401, .481/.396, .575/.389; (b) .668/.384, .691/.367, .718/.314. Unlike A119b dev, both HG002 verdict arms lie
INSIDE their curves (no above-range pricing).

What follows (as registered): an exon-sharing edge gate inside `mcl_families`, byte-identical on membership to this route,
proposed opt-in at I 2.0; defaults are the user's call.
