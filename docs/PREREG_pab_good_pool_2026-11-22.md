# Row 1889: the Sumatran testis comparison with the seeding pool of `all` (pre-registered 2026-11-22)

Row 1880's one lost comparison: StringTie on Sumatran orangutan testis (`PAB_mm`), -686 correct chains at StringTie's size
against RefSeq. That row ran `assemble` on its own, which seeds with primary alignments (row 1812); the full pipeline `all`
seeds with the "good" pool (tied secondary alignments at AS >= 0.98 x best), and this library carries about three secondary
records per read. StringTie reads every record of the BAM. This row runs the list once with the pool `all` uses and with
today's default repair. No new lever: both are existing defaults of one entry point or the other.

## Runs
`assemble --seed-pool good` + `oneread` on `PAB_mm` (driver of row 1886, repair on by default), scorer `rank.py` of row 1880.
Controls with the same command: `PPY_mm` (the other short-read testis library) and mouse cortex. Workspace `r1888_pab_good_pool` (alias `r1889_pab_good_pool`; the number 1888 was taken by another session while this was being registered).

## Reading, fixed now
- **CLOSED** = on `PAB_mm` the list is AHEAD at StringTie's size and no rival of row 1880 turns to BEHIND on the three datasets.
- **NOT CLOSED** = anything else; then the loss is stated as it is and nothing further is tried on this sample (extension is
  ruled out by row 1883; choosing a ranking on one sample would be fitting it).

## Result (2026-11-22): NOT CLOSED; the loss stands

Jobs 4494843-4494845, every stage exit 0; `r1888_pab_good_pool/score.NAME.txt`. Difference in correct chains, Rustle's list
minus StringTie, row 1880 (primary pool, no repair) -> this row (good pool, repair on):
| dataset | at StringTie's size | top 25% | top 50% | top 75% | other rivals |
|---|---|---|---|---|---|
| Sumatran testis | -686 -> **-582** | -673 -> -753 | -460 -> -538 | -62 -> -103 | IsoQuant +641, FLAIR +1,878, bambu +2,387: all still AHEAD |
| Bornean testis | +666 -> +731 | -235 -> -304 | +326 -> +291 | +842 -> +865 | all AHEAD |
| mouse cortex | +1,272 -> +1,306 | +349 -> +359 | +664 -> +672 | +882 -> +903 | IsoQuant +1,277 |

- **NOT CLOSED as registered.** With the pool of `all` and the repair, `final.gtf` on Sumatran testis is 55,988 / 15,700 / 0.2804
  (row 1880: 54,748 / 15,552 / 0.2841): 148 more correct chains, not 686.
- No rival turns on the three datasets.
- Per the reading fixed above nothing further is tried on this sample. State of the comparison against RefSeq: the ranked list
  is ahead of every other assembler at its own output size on every dataset tested, with two exceptions: StringTie on Sumatran
  orangutan testis (-582 to -686; ahead by 237 against the CAT annotation, row 1880 post hoc), and FLAIR on HG002 unless the
  third tier of row 1884 is counted (the user chose not to write that tier).
