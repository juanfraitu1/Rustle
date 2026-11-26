# PREREG: a repeat filter for `members` that keeps tandem arrays (row 1968, 2026-11-26)

Row 1967: the repeat filter alone (threshold 0.3) lifts the share of places at a relative from 0.44 to 0.72 and keeps
87-98% of the paralogues and pseudogenes, but only 75% (WTC11 r3) and 54% (A119b) of the copies of one source gene:
tandem arrays. The user: register the tandem arrays in a way that keeps them.

## Two ways, both judged
- **V1, narrower repeat classes:** the filter counts only interspersed repeats (SINE, LINE, LTR, DNA, Retroposon), not
  satellite, simple repeat or low complexity, which is how the repeat track labels much of an array.
- **V2, a neighbour exemption:** the filter as in row 1967 (all eight classes), but a place on the same contig within
  1 Mb of its source locus is kept whatever its repeat share: a tandem copy is by definition next to its source. Needs
  no annotation beyond the repeat BED.
`bench/o2_truth/place_classes.py` gets the class list from `PLACE_CLASSES_REP` and writes the neighbour flag;
`members_grid.py` adds the exemption as a third filter column.

## Data
WTC11 replicate 2 (human, CHM13; not used for `members` or rows 1955-1967): the claims. WTC11 replicate 3 and A119b
(both used in row 1967) are reported beside it.

## Claims, at threshold 0.3, each variant against today's unfiltered table on WTC11 replicate 2
- **K1:** copies of one source gene reached are at least 0.9 x today's.
- **K2:** the share of rows at a relative is at least today's + 0.20.
- **K3:** recent paralogues and pseudogenes reached are each at least 0.9 x today's.
A variant that meets all three is the candidate default (as "on when a repeat BED sits beside the genome").
## Predictions (mine)
V2 meets all three. V1 meets K2 and K3 and misses K1 (about 0.8 x): array genes also carry interspersed repeats.
Not established beforehand: how many "copies of one source gene" are neighbours at all; the run reports it.

## Result (jobs 4510198-4510200; threshold 0.3). Share at a relative; reached, with the ratio to today's in brackets.
### WTC11 replicate 2 (the unused sample)
| filter | rows | at a relative | recent paralogues | pseudogenes | copies of one source gene |
|---|---|---|---|---|---|
| off (today) | 14,208 | 0.497 | 446 | 3,554 | 58 |
| all eight classes (row 1967's filter) | 9,303 | 0.726 | 429 (0.96) | 3,495 (0.98) | 51 (0.88) |
| V2: all eight, neighbours kept | 9,381 | 0.724 | 437 (0.98) | 3,498 (0.98) | 52 (0.90) |
| V1: interspersed classes only | 9,341 | 0.726 | 439 (0.98) | 3,501 (0.99) | 51 (0.88) |
| V1 + V2 | 9,413 | 0.724 | 444 (1.00) | 3,504 (0.99) | 52 (0.90) |
### A119b (used before) and WTC11 replicate 3 (used before)
| filter | A119b: at a relative / recent / pseudogenes / copies of one source | WTC11 r3: the same |
|---|---|---|
| off (today) | 0.329 / 2,468 / 4,256 / 1,574 | 0.436 / 488 / 3,595 / 76 |
| all eight classes | 0.673 / 0.87 / 0.91 / 0.54 | 0.717 / 0.94 / 0.98 / 0.75 |
| V2 | 0.669 / 0.93 / 0.92 / 0.60 | 0.715 / 0.97 / 0.98 / 0.80 |
| V1 | 0.660 / 0.90 / 0.93 / 0.72 | 0.717 / 0.96 / 0.98 / 0.76 |
| V1 + V2 | 0.657 / 0.96 / 0.93 / 0.74 | 0.715 / 0.98 / 0.98 / 0.82 |

### Claims (WTC11 replicate 2)
- **K2 HOLDS** for both variants (+0.227 and +0.229). **K3 HOLDS** for both (0.98-0.99).
- **K1 FAILS for both, by one copy:** 51 of 58 (0.879) for V1, 52 of 58 (0.897) for V2 and for the two together; the bar
  was 0.9 (53 copies). My prediction that V2 would meet it was wrong by that one copy.
- **My premise was wrong:** copies of one source gene are mostly NOT tandem neighbours. Reached through a place within
  1 Mb on the same contig: 11 of 58 (WTC11 r2), 14 of 76 (r3), 478 of 1,574 (A119b).
- **What is still lost with both variants:** on A119b 403 copies, by biotype lncRNA 143, unprocessed pseudogene 114,
  processed pseudogene 89, and no protein-coding gene among the five largest classes (224 protein-coding copies are
  kept); the sources with most lost copies are repeat-embedded loci (CR381670.1, LINC01666, FRG1BP, AC018638.4); the
  median lost pair is 58% interspersed repeat at its best place. On WTC11 r2 six copies are lost. These are copies that
  are themselves mostly transposable-element sequence: the filter does to them what it is meant to do.

### Reading
The two together (interspersed classes only, neighbours kept) are the best setting measured: on the unused sample the
share at a relative rises from 0.50 to 0.72 with 100% of the recent paralogues, 99% of the pseudogenes and 52 of 58
copies of one source gene kept. It misses the registered bar for those copies by one. Not a candidate default by the
rule written above; the closest thing to one. Not built into `tools/member_discovery.py` (the neighbour exemption would
be a small addition; the class list is only a matter of which BED is given). The user's decision.

## Built (user 2026-11-26: build the neighbour exemption, filter on when a repeat file sits beside the genome)
- `tools/member_discovery.py`: `--neighbour-bp` (1,000,000; 0 = none): with `--repeats`, a place on its source locus's
  contig within that distance is kept whatever its repeat share.
- Driver, stage `members`: the filter is on when `FASTA.repeats.bed` exists; `RUSTLE_MEMBERS_REPEATS=FILE` names another
  BED, set empty = off; a named file that is missing stops the stage (exit 2), a missing sidecar does not;
  `RUSTLE_MEMBERS_NEIGHBOUR_BP`, `RUSTLE_MEMBERS_MAX_REPEAT` as before. The stage says when the filter is on.
- The BED to use: interspersed repeats only. For CHM13: `r1956_place_classes/rmsk/interspersed.bed` (4,642,787
  intervals of the public RepeatMasker track). NOT placed beside the shared genome: that switches the filter on for
  every session and is the user's step.
Checked through the driver on WTC11 replicate 2 (job 4510202), a private genome folder holding the sidecar:

| run | result |
|---|---|
| sidecar present | filter on: 10,160 of 24,380 rows kept, every one a row of the unfiltered run; 102 kept above the limit as neighbours |
| sidecar present, `RUSTLE_MEMBERS_REPEATS=` | 24,380 rows, 13 columns: the unfiltered table |
| no sidecar | the same 24,380 rows |
| sidecar, `RUSTLE_MEMBERS_NEIGHBOUR_BP=0` | 10,058 rows, highest repeat share 0.499 |
| a named file that does not exist | exit 2 with the message |
Whole `cargo test`: 1,591 passed, 0 failed. (The unfiltered table here has 5 rows of one locus that differ from the table
of the grid run above: that run's prefix held a second GTF for the stage; within this job every run reads the same one.)
