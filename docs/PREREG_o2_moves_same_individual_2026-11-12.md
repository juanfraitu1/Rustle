# PREREG O2 MOVES OFF THE PRIMARY, SAME-INDIVIDUAL SAMPLES (register 1652, 2026-11-12)

Why: r1651's side finding (registered as an O2 error): on HG002 reads against CHM13, when `copy_assign --tie-ms 10`
assigns a read to a copy off its primary's locus it is right 2 times and wrong 121; kept on the locus 168 / 54. Before
designing a guard: is that a cross-individual effect, or does O2 move reads wrongly everywhere?

Runs: none new. `bench/o2_contested/o2_moves.py` on existing slack-0 dumps and primary record spans:
- KB3781 on its own assembly: `r1642_exonplace/pof_{panel,heldout}.dump.tsv`, `r1648_oarfish/kb.spans.tsv`.
- HG002 on its own paternal haplotype: `r1644_hg002_unspliced/dump.tsv`, `r1648_oarfish/hg.spans.tsv`.
- HG002 on CHM13 (r1651) rerun through the same script as the parity check (must print 2 / 121 and 168 / 54).
Rows = decidable assigned rows; primary judged by record overlap; moved = O2's copy span does not overlap the primary
record.

Decision rule (per same-individual sample, on moved rows; needs >= 20 moved rows, else UNDERPOWERED):
- right share >= 0.80 on both samples: moves are safe on the reference individual; the guard to design is
  cross-individual only.
- right share < 0.50 on either: O2's moves are unsafe in general; the guard is designed for every sample.
- otherwise: UNDECIDED, reported with the per-family split.
Readout: kept rows, the six-way table, families holding the wrong rows.

Outputs: `/anvme/workspace/v138da13-rustle/r1652_o2_moves/`.

## RESULT r1652 (2026-11-12, no new run; `r1652_o2_moves/{kb,hg,hg13}.txt`)
Parity: HG002 on CHM13 reprints 2 / 121 moved and 168 / 54 kept.
| sample | kept: right / wrong | moved: right / wrong | moved share of assigned rows |
|---|---|---|---|
| KB3781 on its own assembly | 520 / 46 | 2 / 2 | 4 of 570 (0.7%) |
| HG002 on its own paternal haplotype | 479 / 23 | 1 / 22 | 23 of 525 (4.4%) |
| HG002 on CHM13 (r1651) | 168 / 54 | 2 / 121 | 123 of 345 (36%) |
By the rule: KB3781 UNDERPOWERED (4 moved rows); HG002-on-pat right share 0.043 < 0.50, so **O2's moves off the primary
are unsafe in general and the guard is to be designed for every sample.**
Reading:
- O2 almost never moves a read when reads and reference are one individual (0.7% and 4.4% of its assignments); across
  individuals it moves 36%, and that is where the errors are.
- On every sample, a move is wrong far more often than right: 5 right / 145 wrong pooled.
- One family dominates each time: MCL338 on HG002-on-pat (20 of the 22 wrong moves), MCL189 on CHM13.
- Kept rows: 8% wrong on KB3781, 5% on HG002-on-pat, 24% on CHM13: the primary and O2 wrong together.
What O2 adds, on this evidence, is confidence in the primary (assign vs abstain), not relocation. A guard "never assign
off the primary's locus" would remove 150 rows and lose 5 right ones on these three sets. That guard is its own row
(panel identity, held-out, and the count of rows removed on A119b), not adopted here.
Limit: 23 moved rows in one sample decide the rule, 20 of them one family.
