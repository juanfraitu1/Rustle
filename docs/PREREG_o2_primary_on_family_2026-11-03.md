# PREREG O2 PRIMARY-ON-FAMILY GUARD (register 1514, 2026-11-03)

Why: r1510 held-out FAILED and its post hoc dissection found the cause before any fix was written: 3,465 of the 4,529
held-out assigned rows (77%) are molecules whose PRIMARY record lies on no catalog copy at all (MCL148 3,283: spliced
reads of an uncatalogued parent gene at NC_073231.2:47.0 Mb assigned to an 850-bp single-exon retrocopy; MCL13 148: the
parent at NC_073230.2:52.6 Mb, MAPQ 60, assigned to copy 0 with p 1e-11). They entered the family through a secondary
record within 10 ms of the primary, and the PSV model, which only knows catalog copies, picked the nearest one with
confidence. The r1509 panel had 0 such molecules, which is why it passed. `tie_outside_catalog` (§6gz) marks AS ties
outside the catalog and is 0 on all of them.

Rule: `copy_assign --primary-on-family` (opt-in for this row): a molecule enters the tie gate only if one of its records
in the family's windows is a PRIMARY (not secondary, not supplementary). A molecule seen only through secondaries has its
best placement elsewhere, outside this family: O2 cannot test that origin and must not assign the read. Such molecules
are left out of the gate exactly as unique mappers are (no row). Without the flag: byte-identical.

Run: KB3781, all 494 families, same binary/flags as r1510 plus the flag; truth = r1510's `truth_all10` hits (the new
contested set is a subset of r1510's), `strict_truth` with copy loci requiring a near-perfect best per haplotype
(commit 54b1b4c4). Amendment (2026-11-03, before the run finished): the r1510 loci definition (best + 2) is NOT
reported beside it; strict_truth now has only the near-perfect rule, and bar 2's reference numbers (124 / 25) are from
that rule.

Decision rule (pre-declared), held-out 205 families (the 454 non-panel families with contested rows):
1. Assigned rows whose truth is OUTSIDE every catalog copy fall from 326 to <= 5% of assigned rows.
2. Decidable assigned rows: O2 right / wrong within +-5 of r1510's 124 / 25 (no loss of what was right).
3. Panel: the 1,275 assigned rows unchanged (all had a primary on the family) -> the panel table byte-identical.
Pass all three: the guard becomes part of `--tie-ms` (not a separate flag). Report: held-out assigned count, the
statuses of the removed molecules' rows in r1510, wall.
