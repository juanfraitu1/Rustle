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
