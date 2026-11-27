# PREREG: a missing-copy detector that reads every read of a locus (register 1983)

2026-11-26. User goal: "Judge a missing-copy detector on fresh planted sets: every read of a locus, hidden copies vs a
two-allele null, recall and false flags reported in a registered row." Follows row 1982 (the planted sample; the
existing sequence test sees a hidden copy only when its reads tie). Written before any fresh set was planted.

## Detector (`bench/o2_families/locus_sequences.py`, frozen at this commit of the file)
IsoCon's test as ported in `recurrence_test.py` (row 1926), applied per locus and isoform to EVERY read there:
primaries plus secondaries exactly tied on ms. A sequence counts with >= 3 reads. Sequences sharing half of the
smaller one's reads are one (a tied read is read at each of its loci); cells (locus, isoform) holding one sequence
form a class. flag = a class has more sequences than 2 x the loci it is on (two alleles per reference copy).
One detector only; `copy_phasing` is not run (it needs the legacy assignment inputs).

## Development (not evidence)
Row 1982's three sets (seeds 1939-1941). First rule (slots = 2 x loci of the whole tie group): 10 of 35 hidden copies
found on one set; changed to slots per class after reading the misses: 26 / 28 / 29 of 35, groups with >= 10
hidden-copy reads flagged 15 / 15 / 18 of 22, false flags 1-3 of 23 null groups, 7-9 of 57 on the whole reference.

## Test data
Fresh planted sets: `plant.py` seeds 1983, 1984, 1985, `--private 0.001 --het 0.001`, 60 reads per copy; hidden copies
by `hide_copy.py --seed` 2983, 2984, 2985 (a different draw per set). Scenario E carries the claims, pooled over the
three sets. Null = groups with no read of a hidden copy and >= 20 reads (every planted copy has two alleles), and the
same reads on the whole reference, where every flag is false.

## Claims (scenario E, planted error rate, pooled)
- D1: groups with >= 10 hidden-copy reads flagged: >= 0.60.
- D2: false flags <= 0.10 of null groups on the reduced reference AND <= 0.10 of groups on the whole reference.
- D3: hidden copies >= 6 edits from the nearest remaining copy found: >= 0.70.

## Description
By edits class, tied against uniquely placed, single-locus groups; scenario D; a second arm with errors near real HiFi
(`--sub 0.0005 --ins 0.0004 --dele 0.0004`, 1.3 per kb) on the same seeds and hidden copies.

## Reading
D1-D3 hold: an opt-in candidate, to be run on real samples with haplotype truth (HG002, KB3781) before anything else.
D2 fails: the count of sequences over-calls at two-allele loci; say at which kind of group.

## Prediction
D1 holds (about 0.72), D3 holds (about 0.8), D2 fails on the whole reference (about 0.14, from groups of tied loci)
and holds on the reduced reference's null.

## Caveats carried
"Found" credits a copy when its group is flagged: a large chained group (TSPY, 33 loci) credits every hidden copy in it;
the group-level D1 is the fair unit. Credit uses read names (truth). 60 reads per copy; N block instead of a collapse;
hidden copies are annotated reference copies. A real three-allele locus does not exist in the planted sets.

## Outcome (2026-11-26, job 4512729, `r1983_missing_copy_detector/pooled.txt`; scripts unchanged since registration, `frozen.md5`)
Scenario E, three fresh sets pooled:

| | planted errors (5 per kb) | near-HiFi errors (1.3 per kb) |
|---|---|---|
| D1 groups with >= 10 hidden-copy reads flagged | 56 of 68 (0.824) HOLDS | 56 of 71 (0.789) |
| D3 hidden copies >= 6 edits away found | 37 of 46 (0.804) HOLDS | 37 of 46 (0.804) |
| D2 false flags, null groups of the reduced reference | 10 of 60 (0.167) FAILS | 9 of 64 (0.141) |
| D2 false flags, whole reference | 35 of 171 (0.205) FAILS | 24 of 180 (0.133) |
| whole reference, single-locus groups | 11 of 98 (0.112) | 6 of 105 (0.057) |
| hidden copies found, all | 85 of 104 | 90 of 104 |
| twins / 1-5 / 6-19 / >= 20 edits | 26 of 32 / 21 of 24 / 12 of 12 / 25 of 34 | 30 of 32 / 23 of 24 / 12 of 12 / 25 of 34 |
| tied / placed uniquely | 38 of 48 / 47 of 56 | 43 of 48 / 47 of 56 |

Scenario D (uneven expression): groups 50 of 67 and 49 of 68; copies >= 20 edits away 18 of 34; false flags 8 of 58 and 6 of 69.
Predictions: D1 and D3 held (both under-predicted); D2 predicted to fail on the whole reference (right) and to hold on
the reduced reference's null (wrong, 0.167). False flags on the whole reference: groups of tied loci 24 of 73, single-locus 11 of 98.
Reading as registered: D2 fails, so the count of sequences over-calls at two-allele loci, most in groups of tied loci.
Not done: why a two-allele single locus shows a third sequence (recurring planted errors, or the second isoform's
partial reads); real samples with haplotype truth; a locus with three real alleles.
