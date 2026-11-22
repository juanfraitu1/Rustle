# PREREG: do the copies added by secondaries have reads of their own (row 1863, 2026-11-21)

User (2026-11-21): "lets measure requirement 2" of the multi-copy mode (row 1861): a copy seeded only by near-tie secondaries
should have evidence of its own before it is kept. Written before the measurement.

## Measurement
Members of the named families (as `toolcompare.py`: human MAINSTAYS + YAGS on A119b, GGO_MAINSTAYS on OR6737). A copy is
`primary` when row 1858's `primary/x.gtf` has a spliced transcript mapped to it, `added` when only `good/x.gtf` has one,
`none` otherwise. For each copy, every spliced record (primary or secondary, no supplementaries) that `map_copy` places on
it on its strand, classed by ALIGNMENT SCORE (AS; `as_table` best / second of the molecule):
- `own`: this record holds the molecule's best AS and the second best is lower (lead >= 1); `own_clear` = lead > 10.
- `tie`: this record's AS = best = second.
- `behind`: this record's AS is below the molecule's best.
Per copy: molecules per class. Per family: copies by status and by verdict.

## Verdict per copy (2 molecules = the floor used throughout)
- OWN: >= 2 `own` molecules. TIE ONLY: fewer, but >= 2 `own` + `tie`. BORROWED: otherwise (its reads score better elsewhere).

## How it is read
- The filter of requirement 2 is useful if, among `added` copies on A119b, BORROWED is >= 25% (there is something to
  remove) AND among `primary` copies BORROWED is <= 5% (it does not remove what primaries support).
- TIE ONLY copies cannot be called expressed or not from scores: they are reported as groups, not decided.

## Predictions (mine)
Few `added` copies are OWN (a copy with own reads would have been seeded by primaries); most TSPY additions are TIE ONLY;
BORROWED between a quarter and a half.

## Limits
AS only (not ms, not copy-private positions); the sample's copies differ from CHM13's, so a read of copy A can score equal
or better on copy B; annotated copies of fifteen families, one testis sample per species.

## Result (2026-11-21): a BORROWED filter removes little; the added copies are tie groups
Job 4491330, tables `r1863_copy_evidence/{a1,or}.tsv` (`ev.py`).

| A119b copies | OWN | TIE ONLY | BORROWED | no reads |
|---|---|---|---|---|
| `primary` (93) | 78 | 15 | 0 | |
| `added` (30) | 3 | 23 | 4 | |
| `none` (10) | 4 | 0 | 5 | 1 |

Registered reading: BORROWED among `added` is 4 of 30 (13%, bar >= 25%): NOT USEFUL as a filter on its own; among `primary`
0% (bar <= 5%, holds). Predictions: few OWN and TSPY mostly TIE ONLY held; BORROWED a quarter to a half did not.
- `added` OWN: TBC1D3I (11 own molecules), CT47A7, one BPY2B. BORROWED: 3 TSPY, VCY.
- TSPY: every copy carries the same ~200 molecules (own + tie + behind = 200 at nearly each); 5 copies OWN, 5 `primary` and
  21 `added` TIE ONLY, 3 BORROWED. Ties come in sets of 7, 17, 27-29, 46-52 molecules.
- `primary` TIE ONLY: CT47A2-6, HSFY1 / HSFY2 (1,554 tied molecules), 5 TSPY, 3 XAGE: StringTie's and Rustle's primary copies
  there rest on tied reads whose primary flag is arbitrary.
- `none` yet OWN (>= 2 own spliced molecules, no Rustle transcript in either arm): A119b PRY and PRY2, one TSPY, one XAGE;
  OR6737 6 SPTBN5 and 4 JMJD7-PLA2G4B. A loss that secondaries do not explain.
OR6737: `primary` 25 OWN; `added` 1 OWN (JMJD7-PLA2G4B LOC134758743), 1 BORROWED (LRPAP1 LOC134757218, 4,394 behind).
