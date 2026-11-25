# PREREG r1951: the seeding pool on exact `ms` ties only (2026-11-25)

Why: rows 1949 / 1950 measured, on HG002 reads of known origin, what a secondary record of a near-tied read is worth: a
record that ties the read's best `ms` exactly is truly an origin of the read for 2,529 of 2,586 (0.98); a record 1-10
below the best for 792 of 1,212 (0.65), and nothing tested says which. User, 2026-11-25: can that guide assembly too.
Today's pool admits a secondary alignment at AS >= 0.98 x the read's best AS. Row 1942 tested `ms` within 10 (M2): more
precise on four real samples, fewer copies on HG002 (34 against 37): a trade, not adopted. The exact-tie pool was not tested.

## How (row 1942's method and files, one new arm; the assembler is not changed)
- **M0 = `msgap:0`**: `bench/o2_families/pool_bam.py` keeps only the secondaries at the read's best `ms`; assembled with
  `--seed-pool all`, then `families`, `multicopy`, labels on the ORIGINAL BAM and score table (`r1942_ms_seed_pool/arms.sh`).
- Beside it, already on disk from row 1942: A0 (today) and M2 (within 10) for every dataset.
- Data: planted seeds 1951-1953, scenarios E and S (sample = reference); OR6737, KB3781 (gorilla, `toolcompare.py
  --species gorilla`); A119b, HG002 on CHM13 (`--species human`). All looked at before.

## Claims (row 1942's rule, now over four real samples)
M0 is BETTER than A0 when, on planted reads, expressed units found are not lower in any set, silent units labelled
`expressed` or `shared` are not higher by more than 1 per set, and the mode's isoform precision or recall is higher by
>= 0.01 in the sum without the other falling by more than 0.005; AND on all four real samples annotated copies are not
lower and junction precision is not lower by more than 0.01. EQUIVALENT when all within the same tolerances both ways.
Otherwise a TRADE (state what is gained and lost) or WORSE (nothing gained).
Printed, no bar: secondaries admitted, rescued loci, transcripts, M0 against M2.

Prediction (mine): a TRADE. Planted sets unchanged in units and labels (sample = reference, so the true copies tie
exactly) with isoform precision above M2's; junction precision highest of the three arms on the real samples; annotated
copies lower than A0 on at least HG002 and rescued loci well below M2's.

Workspace `r1951_exact_tie_pool`.

## RESULT (2026-11-25, jobs 4504788-4504792; `r1951_exact_tie_pool/`)
Secondaries admitted by M0: planted 0.10-0.13 of all secondary alignments; OR6737 0.002, KB3781 0.005, A119b 0.010, HG002 0.008.

Planted, the mode (A0 today / M2 within 10 / M0 exact), seeds 1951, 1952, 1953:
| | scenario E (all expressed) | scenario S (half silent) |
|---|---|---|
| expressed copies found | 129 of 133 in every arm | 67, 67, 66 of 70 in every arm |
| isoform recall | equal in every arm | equal in every arm |
| isoform precision A0 | 0.639, 0.645, 0.644 | 0.432, 0.433, 0.421 |
| isoform precision M2 | 0.664, 0.676, 0.679 | 0.445, 0.409, 0.413 |
| isoform precision M0 | 0.744, 0.717, 0.746 | 0.505, 0.488, 0.492 |
| false copies at silent places A0 / M2 / M0 | 0 | 48 / 46 / 34; 42 / 42 / 33; 45 / 43 / 33 |
Labels of expressed units identical in every arm; silent units labelled `expressed` or `shared` not higher.

Real, the named families (annotated copies with a spliced transcript / junction precision / spliced transcripts there):
| | A0 | M2 | M0 | rescued loci A0 / M2 / M0 |
|---|---|---|---|---|
| OR6737 | 37 / 0.67 / 176 | 37 / 0.70 / 135 | 35 / 0.70 / 121 | 1,688 / 1,303 / 961 |
| KB3781 | 28 / 0.63 / 108 | 29 / 0.68 / 98 | 28 / 0.63 / 105 | 514 / 375 / 286 |
| A119b | 127 / 0.49 / 1,096 | 127 / 0.51 / 1,007 | 123 / 0.54 / 865 | 4,398 / 3,402 / 2,340 |
| HG002 | 37 / 0.45 / 442 | 34 / 0.47 / 328 | 31 / 0.46 / 262 | 1,584 / 1,209 / 787 |

**A TRADE, as predicted; not BETTER.** The planted condition holds with room (isoform precision +0.07 to +0.10 in all six
sets, nothing lower), the real one fails: annotated copies are lower on three of four samples (-2, 0, -4, -6).
Junction precision: +0.03, 0, +0.05, +0.01. Wrong in my prediction: M0 is NOT the most precise arm on the real samples
(KB3781 0.63 against M2's 0.68, HG002 0.46 against 0.47); only on A119b.
Reading: on planted reads, where the sample is the reference, true copies tie exactly and the exact pool removes a quarter
of the false copies at no cost. On real samples a copy of the individual often differs from the reference copy by a base
or two under the read, so its reads tie its neighbours within 10 but not exactly, and those copies are lost. The truth of
rows 1949 / 1950 (records 1-10 below the best: 0.65 truly an origin) says the same from the other side: two thirds of
those records are right. The within-10 pool (row 1942) stays the better compromise on real data; neither is adopted.
