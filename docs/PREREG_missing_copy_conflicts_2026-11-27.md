# PREREG: the missing-copy detector counting conflicting sequences, on fresh seeds (register 1986)

2026-11-27. User: "lets take the next step" (after row 1985). Written before any set of this row was planted.

## What is judged
`bench/o2_families/locus_sequences.py --count conflicts` (file frozen: `r1986_conflict_count/frozen.md5`): as row 1983's
detector, but per cell (locus, isoform) the count is the largest set of sequences that pairwise disagree at a site both
cover, not the number of sequences. flag = a cell's count exceeds 2 x the loci of its class.
Row 1985 found that row 1983's false flags at single loci come from two isoforms pooled in one cell, and wrote this
count after reading row 1983's seeds; on those seeds it gave 2 of 60 and 12 of 171 false flags at 53 of 68 groups.

## Test data
Fresh: `plant.py` seeds 1986, 1987, 1988 (`--private 0.001 --het 0.001`, 60 reads per copy), `hide_copy.py --seed`
2986, 2987, 2988. Scenario E carries the claims, three sets pooled, planted error rate (5 per kb).
Nulls as row 1983: groups with no read of a hidden copy and >= 20 reads on the reduced reference; all groups with
>= 20 reads on the whole reference. Row 1983's count of sequences is run on the same sets for the paired comparison.

## Claims (the bars of row 1983, unchanged)
- D1: groups with >= 10 hidden-copy reads flagged >= 0.60.
- D2: false flags <= 0.10 on the null groups AND <= 0.10 on the whole reference.
- D3: hidden copies >= 6 edits from the nearest remaining copy found >= 0.70.

## Description
The count of sequences on the same sets; single-locus against tied groups; scenario D; near-HiFi errors
(`--sub 0.0005 --ins 0.0004 --dele 0.0004`).

## Reading
All three hold: an opt-in candidate; next = real samples with haplotype truth (HG002, KB3781); defaults are the user's call.
D2 holds and D1 or D3 fails: the cleaner count costs too much recall; report the trade.

## Prediction
D1 holds (about 0.75). D2 holds (0.03-0.08). D3 holds narrowly (about 0.75) at planted errors and may fall under 0.70
in the near-HiFi arm (0.696 on row 1983's seeds).

## Outcome (2026-11-27, job 4512972, `r1986_conflict_count/pooled.txt`; scripts unchanged, `frozen.md5`)
Scenario E, three fresh sets pooled:

| | bar | conflicts, planted errors | sequences, planted errors | conflicts, near-HiFi | sequences, near-HiFi |
|---|---|---|---|---|---|
| D1 groups with a hidden copy flagged | >= 0.60 | 42 of 64 (0.656) HOLDS | 47 of 64 (0.734) | 38 of 67 (0.567) | 44 of 67 (0.657) |
| D2 false flags, null groups | <= 0.10 | 3 of 66 (0.045) HOLDS | 8 of 66 (0.121) | 2 of 74 (0.027) | 6 of 74 (0.081) |
| D2 false flags, whole reference | <= 0.10 | 10 of 172 (0.058) HOLDS | 26 of 172 (0.151) | 4 of 181 (0.022) | 21 of 181 (0.116) |
| whole reference: single locus / tied | | 1 of 101 / 9 of 71 | 5 of 101 / 21 of 71 | 1 of 108 / 3 of 73 | 5 of 108 / 16 of 73 |
| D3 hidden copies >= 6 edits away found | >= 0.70 | 26 of 45 (0.578) FAILS | 32 of 45 (0.711) | 25 of 45 (0.556) | 29 of 45 (0.644) |
| hidden copies >= 20 edits away found | | 18 of 34 | 23 of 34 | 16 of 34 | 19 of 34 |
| hidden copies found, tied / placed uniquely | | 35 of 52 / 33 of 53 | 37 of 52 / 37 of 53 | 23 of 51 / 35 of 54 | 26 of 51 / 38 of 54 |

Scenario D: conflicts 35 of 64 (planted) and 41 of 66 (near-HiFi) groups, false flags 3 of 65 and 3 of 75; sequences
40 of 64 and 46 of 66, false flags 8 of 65 and 6 of 75.
Predictions: D1 held (0.656, predicted 0.75), D2 held as predicted, D3 predicted to hold narrowly and FAILED (0.578).
Reading as registered: D2 holds and D3 fails, so the cleaner count costs recall. Flag rate where a copy is hidden over
where none is (whole reference): sequences 4.9 x, conflicts 11 x; near-HiFi 5.7 x and 26 x.
Seed dependence: the count of sequences itself gave 0.824 / 0.804 on row 1983's seeds and 0.734 / 0.711 here; the
hidden copies differ. Quote ranges over both rows, not one of them.
