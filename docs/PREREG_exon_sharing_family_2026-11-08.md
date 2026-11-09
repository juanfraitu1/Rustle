# PREREG row 1592: exon sharing as the families-stage edge rule (family level, pairscore), held out on HG002

Written 2026-11-08, after row 1591 (edge level, PASS: SHARE beats exon coverage at equal recall, HG002 +0.368) and
before any family-level exon-sharing run. Setting chosen on A119b and appended here before any held-out run.

## Rule

Replace the shipped edge gate (best record covers >= 60% of the smaller member's exonic length, `--min-shared-exon-frac
0.60`) by **SHARE >= T**: the share of the fewer-exon member's exons that have a >= 50% reciprocal partner exon through
the pair's best record (`bench/denovo_eval/exonshare.py --locus-scores`, max over a locus pair's records).

Route (no binary change):
1. `mcl_families ... --min-shared-exon-frac 0.0 --dump-graph g0.tsv` (gate off, no damping): every candidate edge with
   the binary's own weight. Also `--min-shared-exon-frac 0.60 --dump-graph g60.tsv` (the shipped graph).
2. Edges of g0 whose locus pair has SHARE >= T -> `--extra-edges FILE --extra-edges-only` (weights unchanged; nodes
   unchanged; `--sd-damp` applies after, as for any edge).
3. **Route check before any arm is scored:** g60 fed back the same way at I = 2.8 must give `PREFIX.fam.copies.tsv`
   byte-identical to the shipped run's (the g60 dump run itself; the route run uses `--min-shared-exon-frac 0.0` exactly as the arms do). Else
   stop and fix the route.

## Setting (dev A119b; fixed before any held-out run)

Two tracks, each with candidates T in {0.5, 0.6, 0.7} x I in {2.0, 2.8, 3.6}:
- **(a) plain:** no damping; baseline = the binary plain inflation curve (row 1586: I 1.6 .689/.761, 1.8 .733/.753,
  2.0 .750/.750, 2.4 .781/.743, 2.8 .802/.740, 3.2 .837/.723, 4.0 .863/.704). Eligible: recall >= .740 - 0.04.
- **(b) damped:** `--sd-damp 0.25` (CHM13 SEDEF); baseline = the damped curve (row 1581: I 1.6 .883/.736, 2.0
  .909/.716, 2.4 .914/.722, 2.8 .925/.718, 3.2 .943/.716, 3.6 .943/.688). Eligible: recall >= .718.

Per track, choose the eligible candidate with the largest precision gain over its baseline curve at its recall (linear
interpolation between the two bracketing points, points sorted by recall). A track goes to held-out only if that dev
gain is >= +0.04 and its recall lies inside the baseline's recall range.

## Held-out (HG002) and bars

The chosen (T, I) of each track at I in {1.6, 2.0, 2.8, 3.6, 5.0} is reported; the verdict arm is the chosen (T, I).
- **(a) PASS** iff precision >= HG002 plain curve (`r1580/hg002/plain_I*`) at its recall + 0.03 AND recall >= HG002
  plain default (.345) - 0.04.
- **(b) PASS** iff precision >= HG002 damped curve (`r1580/hg002/sd_I*`) at its recall + 0.03 AND recall >= HG002
  damped default (.321).
Recall outside the baseline's range: FAIL. Each track is judged on its own. Gorilla is not run (no exon-sharing dev
there and 88 truth pairs).

## What follows

- A track PASSES: an exon-sharing gate inside `mcl_families` (binary port, byte-identical to this route) proposed
  opt-in; defaults are the user's call.
- Both FAIL: register; the edge-level gain of row 1591 does not survive MCL.

## Amendment 1 (2026-11-08, after the route check, before any arm is run or scored)

Job 4453793 (dev, gate-minimap2): the route run gives `route.fam.copies.tsv` not byte-identical to `g60.fam.copies.tsv`,
but the same 1,565 families with identical member sets (0 differing). The bytes differ only because
`--extra-edges-only` clears the idents (`max_family_identity` = NA, known before this row) and copy order follows them.
The check is restated as: **identical partition (family member sets) and identical `pairscore.py` line**; it passes on
the partition, and the pairscore line is printed at the top of the arms job (`r1592/arms.sh`) before any arm is scored.
Also fixed before any arm: SHARE is computed on `g0.fam.loci.gated.paf` (the gate-0 run whose edges are filtered), not
on g60's, so that edges only g0 holds can be scored. A binary port of a passing track would still be checked
byte-identical on membership.

## Outcome (2026-11-08/09, A119b dev, jobs 4453793 / 4453848): STOP on both tracks, no held-out run

Route check (amended): `route.fam.copies.tsv` 2013 / 1615 hit, .802/.740 = g60 exactly. Exon-sharing scores on
`g0.fam.loci.gated.paf`: 4,139,552 records, 2,876,406 locus pairs. Gate-0 graph 1,154,423 edges (shipped g60: 525,890);
kept at T 0.5 / 0.6 / 0.7: 435,224 / 130,884 / 130,211. Gains below use exact counts (truth 2,183; the plain curve
re-scored from `r1571/clu/I*` + g60 reproduces the registered values).

| arm | precision | recall | baseline at recall | gain |
|---|---|---|---|---|
| a T0.5 I2.8 | .8949 | .7100 | plain .8550 | **+.0399** |
| a T0.5 I2.0 | .8685 | .7833 | above plain range (max .7613) | not priced |
| b T0.5 I2.8 | .9213 | .7293 | damped .8983 | **+.0230** |
| b T0.5 I2.0 | .9177 | .7820 | above damped range (max .7361) | not priced |
| all T0.6 / T0.7, T0.5 I3.6 | .906-.956 | .602-.708 | | below the recall floors |

(a) best eligible gain +0.0399 < +0.04: **STOP** (by 0.0001 on exact counts; +0.040 on the rounded pairscore columns).
(b) +0.023: **STOP**. Both as registered; no HG002 run.

Post hoc, no claim: the T 0.5 I 2.0 arms reach the highest recall of any families run (.783 / .782) at .869 / .918, where
the rule has no baseline point. Plain inflation reaches that recall only at I 1.4 (.627/.781, `r1571/clu`, outside the
registered curve): +0.24 at equal recall descriptively. T >= 0.6 cuts the graph to a quarter of the shipped edges and
recall to .60-.65 at every inflation (2-exon loci share 1 of 2 = 0.5). Testing T 0.5 at low inflation needs a new
registered row with a baseline curve that reaches recall .78.
