# PREREG: FLAIR with its own multimapper option, `--allow_paralogs` (row 1946, 2026-11-25)

`docs/MULTIMAPPERS_ACROSS_ASSEMBLERS_2026-11-25.md` found that FLAIR has a switch for reads that fit two paralogs' isoforms
equally (`--allow_paralogs`: the read counts for every tied isoform instead of the one whose name sorts first). It was never
an arm of the family benchmark. The user asked to test it (2026-11-25): it is the fairest "another tool with its own
multimapper option" comparison.

## Arms (planted human named families, reads with ragged ends, today's defaults)
- FLAIR as in the benchmark (`flair collapse --trust_ends`) and the same command with `--allow_paralogs`, on the same read
  table and reads.
- Rustle on primaries, Rustle's multi-copy mode, the mode with the EM labels; StringTie for reference.

## Data
Fresh seeds 1964, 1965, 1966. (a) scenarios E (every copy expressed) and S (half the copies silent) at 60 reads per copy;
(b) scenario E at 3 reads per copy (`plant.py --reads 3`), the coverage where rows 1917-1918 found tied secondaries to help.

## Measures
`units_score.py`: expressed units found, false units (silent units with a transcript). `score.py`: isoform recall and
precision. For the labelled mode: silent units labelled `expressed` or `shared` (`label_score.py`).

## Claims
- **P1, the option helps FLAIR:** summed over the three seeds at 3 reads, `--allow_paralogs` finds >= 3 more expressed units
  than FLAIR's default, AND summed over the three S sets at 60 reads it has at most 3 more false units.
- **P2, Rustle still differs:** at 3 reads the multi-copy mode finds more expressed units than FLAIR `--allow_paralogs` in
  each of the three seeds; AND in each S set the labelled mode has fewer silent units labelled `expressed` or `shared`
  than FLAIR `--allow_paralogs` has false units.
- Reported whatever the outcome: isoform recall and precision of every arm, and the number of transcripts the option adds.

## Predictions (mine)
P1 fails: the option acts after FLAIR's first-pass isoforms were built from primary alignments only, so it can add reads
to isoforms that exist but cannot create a copy that received no primary; I expect at most one unit gained. P2 holds.
On the one set looked at for the study (seed 1951, S, 60 reads; not part of this test) the option added 13 transcripts,
one false copy and 0.008 isoform recall.

## Result (job 4504621, seeds 1964 / 1965 / 1966)

| | FLAIR | FLAIR `--allow_paralogs` | Rustle primaries | Rustle mode | mode, labelled | StringTie |
|---|---|---|---|---|---|---|
| expressed units found at 3 reads (of 108) | 44, 40, 40 | 54, 48, 47 | 71, 62, 76 | 104, 102, 104 | - | 74, 76, 79 |
| isoform recall at 3 reads | 0.197, 0.172, 0.197 | 0.280, 0.264, 0.272 | 0.331, 0.318, 0.356 | 0.649, 0.690, 0.699 | - | 0.293, 0.272, 0.322 |
| isoform precision at 3 reads | 0.855, 0.872, 0.940 | 0.882, 0.887, 0.929 | 0.908, 0.916, 0.895 | 0.783, 0.717, 0.795 | - | 0.769, 0.722, 0.811 |
| false units, half silent, 60 reads (of 48, 47, 47) | 18, 13, 15 | 19, 15, 17 | 16, 15, 12 | 36, 29, 34 | 4, 3, 4 labelled `expressed` or `shared` | 17, 15, 14 |
| isoform recall, all expressed, 60 reads | 0.820, 0.816, 0.816 | 0.824, 0.824, 0.824 | 0.824, 0.833, 0.816 | 0.824, 0.824, 0.820 | - | 0.661, 0.674, 0.657 |
| isoform precision, all expressed, 60 reads | 0.742, 0.753, 0.741 | 0.701, 0.719, 0.704 | 0.811, 0.829, 0.823 | 0.644, 0.657, 0.651 | - | 0.763, 0.797, 0.766 |
| isoform recall / precision, half silent | 0.828 / 0.587, 0.810 / 0.580, 0.802 / 0.574 | 0.869 / 0.561, 0.825 / 0.547, 0.810 / 0.545 | 0.852 / 0.605, 0.817 / 0.566, 0.810 / 0.593 | 0.877 / 0.433, 0.825 / 0.466, 0.817 / 0.424 | - | 0.697 / 0.607, 0.667 / 0.596, 0.651 / 0.607 |

The option adds 15-26 transcripts per set (for example 283 -> 306, 55 -> 80).

- **P1 FAILS as registered, on its second half.** The option does help FLAIR at 3 reads: +10, +8, +7 = 25 more expressed
  units (bar 3). But it adds 1 + 2 + 2 = 5 false units over the three half-silent sets (bar: at most 3).
- **P2 HOLDS.** At 3 reads the mode finds 104, 102, 104 units against 54, 48, 47; in the half-silent sets the labelled mode
  marks 4, 3, 4 silent units `expressed` or `shared` against 19, 15, 17 false units.
- **My prediction was wrong on P1's first half:** I expected at most one unit gained; it gained 25. The option does more
  than add reads to existing isoforms: with tied reads counted at every tied isoform, isoforms pass FLAIR's support
  threshold that did not before.
- **Fairness caveat for the 3-read rows:** FLAIR requires 3 supporting reads per isoform (`--support`, default 3), so at 3
  reads per copy it sits at its own floor; Rustle's floor is 2 reads and its family rescue 1. FLAIR's low recall there is
  mostly that threshold, not its treatment of multimappers. The comparison that isolates the option is FLAIR against FLAIR.
- At 60 reads the option changes little: isoform recall +0.004 to +0.041, precision -0.03 to -0.04, 1-2 more false units.
- Unlabelled, Rustle's mode has the MOST false units (36, 29, 34); the labels are what bring it to 4, 3, 4. Without the
  labels the mode must not be called more precise than any tool here.
