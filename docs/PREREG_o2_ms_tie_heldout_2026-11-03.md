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
