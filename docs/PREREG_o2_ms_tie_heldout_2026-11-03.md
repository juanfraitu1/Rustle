# PREREG O2 MS-TIE HELD-OUT (register 1510, 2026-11-03)

Why: r1509 failed as registered because its truth (haplotype hits within NM + 2 of the best) cannot score near ties on one
individual's near-identical copies; its post hoc slack-0 re-read (unique best haplotype hit) gave O2 284 right / 2 wrong on
286 reads against the primary's 280 / 6. That number was read after the fact on the same panel. This row fixes the truth
BEFORE a new run and scores reads r1509 never aligned.

Truth (fixed now): `bench/o2_contested/strict_truth.py` at g = 0: a read's truth = the catalog copies overlapping its
haplotype hits with NM equal to its best, given a near-perfect best (NM <= max(2, 0.003 x aligned)); decidable when one
locus; rows whose best-NM hits sit on >= 2 loci are UNVERIFIABLE (identical loci over the read) and are reported as their
own class, never as successes. The primary's verdict: its genome-wide position inside a truth copy.

Run: KB3781, ALL 494 families, `copy_assign --tie-ms 10` (same binary and flags as r1509 otherwise), haplotype_truth on
every contested row (SAMPLE raised so abstained rows are all aligned: `RUSTLE_TRUTH_SAMPLE`), then strict_truth.
Held-out set = the molecules of the 454 families outside r1509's 40-family panel plus, inside the panel, molecules not in
r1509's MS10 table (none expected). Reported beside it, labelled, the panel again.

Decision rule (pre-declared):
1. Held-out decidable assigned molecules >= 30, else UNDERPOWERED (report only).
2. On them, O2 wrong rate <= half the primary's wrong rate on the same molecules, AND O2 assigns >= 20% of the held-out
   decidable molecules.
3. Report, no bar: unverifiable fraction of assignments; abstained-but-decidable count (coverage O2 leaves); molecules
   (not rows) contested / assigned; wall and run time.
Counting unit everywhere: distinct molecules (a read x family row is folded to its read's best row: assigned > abstained).

## RESULT (2026-11-03, job 4421539; `~/o2contested/r1510_strict.txt`, per-row dumps workspace o2mstie/*.strict_dump.tsv)
Run: all 494 families, 116,915 contested rows, 5,804 assigned, assign wall 462 s; truth on every contested row (970 s).
Held-out (205 families with contested rows outside the panel): 764 assigned rows with hits, 304 decidable at g = 0:
O2 143 right / 161 wrong; the primary 147 / 157 on the same rows. Bar 2 FAILS as registered (O2 wrong rate is not
half the primary's; it is the same). Assigned 304 of 8,424 decidable held-out molecules (3.6%): bar 2's 20% coverage
also FAILS. ❌ FAIL.
Panel (the 40 r1509 families, rerun inside the all-families job): 290 decidable, O2 288 / 2, primary 284 / 6 (r1509's
slack-0 numbers reproduced).
Post hoc (flagged as such): 131 of the 161 held-out wrong rows are one family, MCL13 (3 single-exon copies of ~760 bp on
three chromosomes + a 3-exon copy on chrX), all 131 assigned to copy 0 with the primary also on copy 0, while the reads'
unique best haplotype hit is copy 2's locus. Without MCL13: O2 143 / 30, primary 143 / 26: O2 is no better than the
primary on the held-out families either way. The gain seen on the panel (288 / 2 vs 284 / 6) does not transfer.
