# PREREG: should `members` default to 0.9 of the best score with the repeat filter on? (row 1967, 2026-11-26)

The user: check whether that change of defaults is beneficial. Today: `RUSTLE_MEMBERS_MIN_RATIO` 0.3, no repeat filter
(`RUSTLE_MEMBERS_REPEATS` unset). The stage's purpose (row 1910): the relatives of every assembled locus through ALL of
its reads' alignments, "where relatives are, expressed or not". Rows 1956 / 1957 / 1960 (A119b): relatives are 84.6% of
the rows at 0.9 and above with the filter and 66.9% at 0.7-0.9. Row 1955 (planted): silent pseudogenes are reached at
0.3-0.9, almost never at 0.9 and above (5 of 417 against 173 at 0.3). So the change buys precision and may cost the
stage its purpose. Both sides are measured, on the sample the numbers came from and on one that was not used.

## Data
A119b (development: `r1915_members_blocks/a1.members.tsv`) and WTC11 replicate 3 (human, CHM13, not used for `members` or
for rows 1955-1960; assembled and run here with `RUSTLE_MEMBERS_MAX_REPEAT=2`, the column only).

## Grid, from one unfiltered table per sample (`bench/o2_truth/members_grid.py`)
Threshold 0.3 (today), 0.5, 0.7, 0.9 x repeat filter off / on (footprint under half repeat).
- PRECISION: of the (place, source) rows with an annotated source, the share at a relative (classes 1-3 of row 1956, or
  another class with the similarity test of row 1957 positive).
- REACH: of the listed relatives of the sources that have any row (relation table), the share that some kept row lands
  on, by kind: recent paralogue, older paralogue, pseudogene, copy of the same source gene.
- size: rows, distinct places, places on no assembled locus.

## Claims, all on WTC11 (the unused sample); A119b is reported beside it
- **M1:** at 0.9 with the filter, precision is at least 0.80.
- **M2:** that is at least 0.30 above today's default.
- **M3:** at 0.9 with the filter the reach of recent paralogues is at least 0.8 x today's.
- **M4:** at 0.9 with the filter the reach of pseudogenes is at least 0.5 x today's.
The change is beneficial as a default only if all four hold. Also read off the grid: the setting with the filter on and
the threshold unchanged.
## Predictions (mine)
M1 and M2 hold, M3 holds, M4 FAILS (pseudogenes sit below 0.9). So I expect to recommend the filter as a default and
the threshold left where it is, with the score ratio kept as a column to cut on.
A practical limit whatever the result: the filter needs a repeat BED for the genome, which the pipeline does not ship;
a default can only mean "on when `FASTA.repeats.bed` exists", like the exclusion list.

## Result (jobs 4510133 / 4510134). Precision = share of rows at a relative; reach in brackets = share of the listed relatives.
### WTC11 replicate 3 (the unused sample; 16,256 rows with an annotated source at today's default)
| threshold | filter | rows | distinct places | at a relative | recent paralogues (628) | older paralogues (24,330) | pseudogenes (5,391) | copies of one source gene (545) |
|---|---|---|---|---|---|---|---|---|
| 0.3 (today) | off | 16,256 | 9,425 | 0.436 | 488 (0.777) | 376 | 3,595 (0.667) | 76 |
| 0.3 | on | 9,147 | 7,017 | 0.717 | 459 (0.731) | 352 | 3,523 (0.653) | 57 |
| 0.5 | on | 7,517 | 5,738 | 0.742 | 421 (0.670) | 227 | 2,937 (0.545) | 57 |
| 0.7 | on | 4,797 | 3,584 | 0.774 | 341 (0.543) | 106 | 1,879 (0.349) | 57 |
| 0.9 | off | 2,499 | 1,723 | 0.742 | 253 (0.403) | 41 | 746 (0.138) | 62 |
| 0.9 | on | 2,002 | 1,457 | 0.814 | 233 (0.371) | 41 | 728 (0.135) | 47 |
### A119b (development; 59,643 rows at today's default)
| threshold | filter | rows | distinct places | at a relative | recent paralogues (2,736) | older paralogues (49,965) | pseudogenes (6,541) | copies of one source gene (2,590) |
|---|---|---|---|---|---|---|---|---|
| 0.3 (today) | off | 59,643 | 21,996 | 0.329 | 2,468 (0.902) | 1,458 | 4,256 (0.651) | 1,574 |
| 0.3 | on | 19,287 | 9,034 | 0.673 | 2,142 (0.783) | 871 | 3,872 (0.592) | 852 |
| 0.5 | on | 15,819 | 7,313 | 0.710 | 2,118 (0.774) | 677 | 3,274 (0.501) | 852 |
| 0.7 | on | 11,096 | 4,753 | 0.767 | 1,979 (0.723) | 316 | 2,208 (0.338) | 846 |
| 0.9 | off | 11,518 | 3,500 | 0.707 | 1,913 (0.699) | 386 | 1,029 (0.157) | 1,252 |
| 0.9 | on | 6,127 | 2,088 | 0.846 | 1,668 (0.610) | 123 | 919 (0.140) | 700 |

### Claims (WTC11)
- **M1 HOLDS** (0.814). **M2 HOLDS** (+0.378 over today's 0.436).
- **M3 FAILS:** recent paralogues reached 233 against 488 (0.48 x; bar 0.8). I predicted it would hold: wrong.
- **M4 FAILS:** pseudogenes reached 728 against 3,595 (0.20 x; bar 0.5). A119b the same: 0.68 x and 0.22 x.
- **So the proposed default (0.9 with the filter) is NOT beneficial:** it keeps a fifth of the pseudogenes and half of
  the recent paralogues the stage finds today, for a table that is 81-85% relatives.

### The other setting on the grid: the filter on, the threshold left at 0.3 (a description, not a registered claim)
| | WTC11 | A119b |
|---|---|---|
| at a relative | 0.436 -> 0.717 | 0.329 -> 0.673 |
| rows | 16,256 -> 9,147 | 59,643 -> 19,287 |
| recent paralogues reached | 0.94 x | 0.87 x |
| pseudogenes reached | 0.98 x | 0.91 x |
| copies of one source gene reached | 0.75 x | 0.54 x |
The filter alone gives most of the precision (+0.28 and +0.34) and keeps 87-98% of the paralogues and pseudogenes. Its
cost is on copies of one source gene, half of them on A119b: those are tandem arrays, which the repeat track partly
annotates as satellite or simple repeat (not tested: the filter without those two classes).
Recommendation: threshold stays 0.3; the repeat filter is the candidate default, as "on when a repeat BED sits beside the
genome"; the score ratio and `repeat_share` stay columns to cut on. The default itself is the user's decision.
