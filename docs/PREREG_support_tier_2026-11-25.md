# PREREG: two tiers of the multi-copy mode, by whether a primary alignment carries the transcript (row 1952, 2026-11-25)

Diagnosis (workspace `false_isoforms_2026-11-25`, seeds 1964 / 1965, a description): at 60 reads per copy the tied
secondary alignments add no exact isoform to the mode (197 = 197) and add 50-56 chains that are no planted isoform of
their copy; 53-78 of the mode's non-exact chains per set are absent from the primaries-only assembly. FLAIR avoids this by
construction: its isoforms come from primary alignments and tied reads are shared among isoforms that exist. The user
agreed (2026-11-25) to test a TIER, not a filter: nothing is dropped.

## The tier (`bench/o2_families/support_tier.py`, from `tools/read_provenance.py`'s transcript table)
`primary` = at least one primary alignment has the transcript's intron chain (exactly, or after the junction replacement
the provenance tool infers), or for a single-exon transcript lies inside it. `tied_only` = none.

## Data
Planted human named families, FRESH seeds 1967, 1968, 1969 (ragged ends, today's defaults): all expressed (E) and half
silent (S) at 60 reads per copy; all expressed at 3 reads per copy. Real: gorilla OR6737, the mode of row 1937.

## Claims
- **C1, precision:** in each of the six sets at 60 reads the `primary` tier's isoform precision is at least 0.10 above the
  whole mode's, with isoform recall not lower by more than 0.01.
- **C2, false units:** in each half-silent set the `primary` tier has fewer false units than the whole mode.
- **C3, the other tier earns its place:** at 3 reads per copy, in each seed at least 25% of the mode's exact isoforms are in
  the `tied_only` tier (so dropping that tier would cost the low-coverage recall).
- **R, real sample:** on OR6737 the `primary` tier's junction precision at the named families is not lower than the whole
  mode's. The annotated copies of each tier are reported, no claim.
- Reported for every set: transcripts, exact isoforms, isoform recall and precision, units and false units of both tiers.

## Predictions (mine)
C1 holds (I expect about 0.80 against 0.65). C2 holds. C3 holds, with 40-50% of the exact isoforms at 3 reads in the
`tied_only` tier. R holds. If C3 fails, the tier is simply a filter and `RUSTLE_SEED_SPLICED=support` would do.

## Result (job 4504820; seeds 1967 / 1968 / 1969)

| | whole mode | `primary` tier | `tied_only` tier | Rustle primaries only | FLAIR | StringTie |
|---|---|---|---|---|---|---|
| isoform precision, all expressed, 60 reads | 0.654, 0.654, 0.659 | 0.811, 0.824, 0.804 | 0 of 58, 62, 54 exact | 0.831, 0.828, 0.812 | 0.769, 0.736, 0.729 | 0.769, 0.777, 0.766 |
| isoform recall, all expressed, 60 reads | 0.824, 0.824, 0.824 | 0.824, 0.824, 0.824 | 0 | 0.820, 0.824, 0.816 | 0.808, 0.803, 0.820 | 0.669, 0.669, 0.657 |
| isoform precision, half silent | 0.452, 0.417, 0.488 | 0.609, 0.621, 0.630 | 0 of 59, 83, 48 exact | 0.646, 0.620, 0.620 | 0.585, 0.597, 0.636 | 0.621, 0.590, 0.601 |
| false units, half silent (of 47, 49, 46) | 29, 34, 23 | 13, 13, 9 | 22, 28, 19 | 8, 12, 7 | 14, 16, 9 | 8, 15, 8 |
| exact isoforms at 3 reads (of 239) | 166, 162, 168 | 142, 140, 146 | 24, 22, 22 | 83, 70, 72 | 37, 49, 41 | 79, 64, 69 |
| isoform precision at 3 reads | 0.755, 0.726, 0.740 | 0.835, 0.819, 0.811 | 0.480, 0.423, 0.468 | 0.954, 0.921, 0.935 | 0.841, 0.961, 0.872 | 0.760, 0.744, 0.750 |
| expressed units found at 3 reads (of 108) | 103, 103, 103 | 94, 97, 101 | 34, 25, 28 | 66, 65, 64 | 40, 42, 42 | 84, 75, 81 |

Gorilla OR6737 (named families: annotated copies with a spliced transcript / junction precision / spliced transcripts):
whole mode 37 / 0.67 / 176; `primary` tier 35 / 0.68 / 125; `tied_only` tier 19 / 0.64 / 51. Genome-wide the tier split is
70,718 transcripts with a primary alignment and 2,455 without (3.4%).

- **C1 HOLDS in all six sets:** precision +0.14 to +0.20 at identical recall. At 60 reads the `tied_only` tier holds NOT
  ONE exact isoform (0 of 364 transcripts over the six sets): everything the tied secondaries add there without a primary
  is a chain the copy does not express.
- **C2 HOLDS:** 13, 13, 9 false units against 29, 34, 23.
- **C3 FAILS:** at 3 reads the `tied_only` tier holds 24 of 166, 22 of 162 and 22 of 168 exact isoforms (13-15%, bar 25%),
  at precision 0.42-0.48. My prediction (40-50%) was wrong. The mode's low-coverage recall does NOT come mainly from
  transcripts built on secondaries alone: the `primary` tier by itself has 140-146 exact isoforms against 70-83 for the
  primaries-only assembly. Tied alignments help mostly by adding SUPPORT to chains a primary alignment carries (they lift
  them over the read floor), which is the way FLAIR's paralog option uses them.
- **R HOLDS:** 0.68 against 0.67. But this measure cannot see the problem: a paralog's isoform laid onto a copy has that
  copy's annotated junctions, so the `tied_only` tier scores 0.64 on real data while it is 0% exact on planted reads.
  Junction precision at annotated copies is not evidence that a tied-only transcript is expressed.
- **What it gives.** The `primary` tier alone: isoform precision 0.80-0.82 at full recall at 60 reads (FLAIR 0.73-0.77,
  StringTie 0.77), and at 3 reads 0.59-0.61 recall at 0.81-0.84 precision (a FLAIR with two supporting reads and its
  paralog option, row 1948: 0.46-0.51 at 0.83-0.86). The `tied_only` tier costs 13-15% of the low-coverage exact isoforms
  and 2-9 units if left out.
- Not decided, not wired: the tier is a tag made after the run by `support_tier.py` from the provenance table. Whether the
  tag becomes part of stage `o2`, and whether any product should show only the `primary` tier, is the user's call.
