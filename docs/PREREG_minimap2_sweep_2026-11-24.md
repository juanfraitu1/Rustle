# PREREG: can a minimap2 parameter change remove or decide the near-tied reads (row 1940, 2026-11-24)

## Outcome (2026-11-24, `r1940_minimap2_sweep/cells/*.tsv`)
NO CELL AVOIDS, NO CELL DECIDES. B0 reproduces 553 / 47 and 18.
- Run: 36 of 37 cells scored. The first job stood still behind k 11 (repeat cutoff 76,192); resubmitted with a 40 minute
  limit per cell (deviation from the registration); k 11 did not finish in 40 minutes (base 7.5).
- Four cells were no test: `-g 2000` and `-z 200` are the preset's own values, `-C 0` and `-C 9` gave byte-identical
  alignments (cause not confirmed). So 32 real cells.
- Avoid: reads with an exactly tied rival 3,513 (annotated junctions) to 5,046 over all cells, base 3,684; bar 1,842.
- Decide: best cell on the tuning half is w 10 (20 wrong of 300 against 29, S1 guard kept); on the judging half it has
  20 wrong of 300 against the baseline 18 (bar 12). All reads: 40 against 47 (11 fixed, 4 broken).
- Range of wrong near ties of 600 over the other real cells: 43 to 54; bandwidth 200 / 1,000 give 353 / 282 and 3,865 /
  3,539 wrong clear-lead reads.
Prediction (no cell avoids, none decides; seeding moves counts by about 10% at most): right, w 10 moved it by 15% on all reads.
Limit: one-at-a-time cells on one sample; combinations and other aligner versions not run.
Page: section "32 minimap2 settings: the ties stay" of artifact 688f0fe1-7489-4d2f-92bc-fd94b55522b2.

## Current state
Run: no cell avoids, no cell decides.
