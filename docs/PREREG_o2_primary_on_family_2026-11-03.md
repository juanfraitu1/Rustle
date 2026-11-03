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

## RESULT (2026-11-03, job 4421555; `~/o2contested/r1514_strict.txt`, dumps workspace o2mstie/pof_*.strict3_dump.tsv)
Run: 494 families, wall 341 s (r1510: 462 s); rows 54,733 (116,915), assigned 2,543 (5,804).
Held-out (205 families): assigned 4,529 -> 1,109; the primary of every assigned molecule now lies within 50 kb of a copy
of its own family (1,105) or of another family (4); the 3,465 no-copy assignments of r1510 are gone (their r1510 rows:
3,471 assigned, 45,155 ambiguous, 117 tied molecules removed in all).
Bar 1 (truth OUTSIDE every copy among assigned rows with hits <= 5%): 188 of 601 = 31% under the registered near-perfect
truth: FAIL. Under haplotype_truth's original loci rule (best + 2, `STRICT_LOCI=best2`, comparison only) it is 27 of 601
= 4.5% (r1510: 36 of 764); the registered truth calls a read "outside" whenever its copy has no near-perfect allele on
that haplotype, so most of the 188 are the truth's gap, not O2's; as registered the bar still FAILS.
Bar 2 (decidable right / wrong within +-5 of 124 / 25): 129 / 22: PASS. The primary on the same rows: 125 / 26.
Cross-tab on the 151 decidable held-out molecules: both right 124, O2 right where the primary is wrong 5, the reverse 1,
both wrong 21. O2's gain over the aligner on the ms near-tie set is 5 vs 1 molecules; the errors are shared.
Bar 3 (panel byte-identical): FAIL. All 1,275 r1510 assignments kept and 159 gained (MCL9 107, MCL3 26, MCL111 26):
removing the outside molecules from a family's pool changes its PSV model and p-values, so more same-family molecules pass.
Panel decidable: 393 / 26 (primary 390 / 29); r1510's 302 / 16 on the smaller set.
Verdict as registered: ❌ FAIL (bars 1 and 3). What the row established: the outside-origin assignment class is removed
entirely at no cost to what was right, and the panel gains assignments. The flag stays opt-in; a re-registration with
bars stated on the near-perfect truth's decidable class and "panel superset" instead of identity is the honest next row,
then default. Follow-up (not registered here): the MCL148 class is also a shape mismatch (spliced read, single-exon copy);
a separate guard, measured separately.
