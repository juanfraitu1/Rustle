# PREREG: a tuned FLAIR, and FLAIR's paralog option on a real sample (row 1948, 2026-11-25)

The two items left open by `docs/MULTIMAPPERS_ACROSS_ASSEMBLERS_2026-11-25.md` (user: "can we do the tuned FLAIR and the
paralog option for a real sample?").

## A. Tuned FLAIR on planted reads
Row 1946's caveat: FLAIR needs 3 supporting reads per isoform, so at 3 reads per copy it sat at its own floor. Here FLAIR
gets the floor Rustle has: `--support 2`, alone and with `--allow_paralogs`, on the nine planted sets of row 1946 (seeds
1964-1966; these two arms have not been run on them). Compared with the arms of row 1946 (already computed).
- **T1:** at 3 reads per copy Rustle's multi-copy mode finds more expressed units than FLAIR `--support 2
  --allow_paralogs` in each of the three seeds.
- **T2:** in each half-silent set at 60 reads, FLAIR `--support 2 --allow_paralogs` has more false units than the labelled
  mode has silent units labelled `expressed` or `shared`.
- Reported: units, false units, isoform recall and precision of both tuned arms on all nine sets.
Predictions (mine): the lower floor helps FLAIR a lot at 3 reads (I expect 75-95 units of 108 with both switches) and T1
still holds; T2 holds, and FLAIR's false units rise with the lower floor.

## B. FLAIR's paralog option on a real sample
Gorilla testis OR6737 (reference mGorGor1). FLAIR as in the benchmark (`flair collapse --trust_ends` on its own read table
of the same BAM) and the same with `--allow_paralogs`. Scored with `toolcompare.py --species gorilla` (the named families'
annotated copies with a spliced transcript, junction precision) beside Rustle's multi-copy mode (run of row 1937) and
StringTie (the run of 2026-11-21).
- **R1:** `--allow_paralogs` models at least as many annotated family copies as FLAIR's default, with junction precision not
  lower by more than 0.02.
- **R2:** Rustle's mode models more annotated family copies than FLAIR `--allow_paralogs`.
Predictions (mine): R1 holds with a gain of 0-2 copies; R2 holds (the mode had 37 on this sample against StringTie's 28).
One real sample only; no truth for expression there (copies MODELED, not shown expressed).

## Result (jobs 4504691 and 4504692; the paralog arm of part A was rerun: my first job script put the reads inside FLAIR's temporary folder, which the first FLAIR run empties)

### A. Tuned FLAIR on planted reads (seeds 1964 / 1965 / 1966)
| | FLAIR (row 1946) | FLAIR `--allow_paralogs` (row 1946) | FLAIR `--support 2` | FLAIR `--support 2 --allow_paralogs` | Rustle mode | mode, labelled |
|---|---|---|---|---|---|---|
| expressed units found at 3 reads (of 108) | 44, 40, 40 | 54, 48, 47 | 95, 99, 90 | 98, 101, 96 | 104, 102, 104 | - |
| isoform recall at 3 reads | 0.20, 0.17, 0.20 | 0.28, 0.26, 0.27 | 0.41, 0.44, 0.41 | 0.46, 0.51, 0.49 | 0.65, 0.69, 0.70 | - |
| isoform precision at 3 reads | 0.86, 0.87, 0.94 | 0.88, 0.89, 0.93 | 0.83, 0.83, 0.87 | 0.84, 0.83, 0.86 | 0.78, 0.72, 0.80 | - |
| false units, half silent, 60 reads | 18, 13, 15 | 19, 15, 17 | 20, 14, 16 | 20, 17, 17 | 36, 29, 34 | 4, 3, 4 |
| isoform recall / precision, all expressed, 60 reads | 0.82 / 0.74 | 0.82 / 0.70 | 0.82 / 0.72-0.74 | 0.82-0.83 / 0.67-0.70 | 0.82 / 0.64-0.66 | - |

- **T1 HOLDS, narrowly:** 104 > 98, 102 > 101, 104 > 96. In seed 1965 the margin is ONE unit.
- **T2 HOLDS:** 20, 17, 17 false units against 4, 3, 4 silent units the labelled mode marks `expressed` or `shared`.
- **What this changes.** Most of the gap of row 1946 at 3 reads was FLAIR's support threshold, as its caveat said: with a
  floor of 2 reads FLAIR finds 90-99 units, with its paralog option 96-101, against the mode's 102-104. At the level of
  UNITS FOUND the mode's lead over a FLAIR given the same floor is 1 to 8 units of 108. What stays clearly different at 3
  reads is ISOFORM recall (0.65-0.70 against 0.46-0.51, at 0.05-0.11 lower precision), and the labels.
- My prediction for the tuned FLAIR was 75-95 units; it reached 96-101.

### B. FLAIR's paralog option on gorilla OR6737 (named families; annotated copies with a spliced transcript / junction precision / spliced transcripts)
| | FLAIR | FLAIR `--allow_paralogs` | Rustle mode | StringTie |
|---|---|---|---|---|
| all four families (49 annotated copies) | 23 / 0.55 / 116 | 24 / 0.48 / 172 | 37 / 0.67 / 176 | 28 / 0.59 / 155 |
| LRPAP1 (11) | 10 / 0.46 | 10 / 0.40 | 11 / 0.67 | 10 / 0.67 |
| MAPKBP1 (9) | 1 / 0.78 | 1 / 0.72 | 2 / 0.92 | 1 / 0.82 |
| JMJD7-PLA2G4B (19) | 12 / 0.57 | 13 / 0.48 | 16 / 0.52 | 16 / 0.53 |
| SPTBN5 (10) | 0 | 0 | 8 / 0.86 | 1 / 1.00 |
| transcripts genome-wide | 153,402 | 208,122 | - | - |

- **R1 FAILS:** one more annotated copy (24 against 23) but junction precision falls from 0.55 to 0.48 (bar: at most 0.02
  lower), and the option adds 36% transcripts genome-wide.
- **R2 HOLDS:** 37 annotated copies against 24.
- One real sample, FLAIR at its default support of 3; copies are MODELED, not shown expressed.
