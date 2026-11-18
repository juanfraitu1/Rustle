# PREREG r1781: the wider gate (`--tie-ms 10`) as an arm on the neutral truth, HG002 (2026-11-17)

Written before the run. Workspace `r1781_tie_ms_arm`, scorer `bench/o2_truth/assign_arm.py` (unchanged from r1778).

## Question
r1778: the default assign stage (exact AS ties only) makes 7 calls on the 710 decidable near-tie reads. Of the 57
primary-wrong reads whose two loci are copies of one family, 33 have no row; all 33 have a best-over-runner-up AS lead of
1-10, so the exact-tie gate excludes them by design. Does the existing opt-in gate `--tie-ms 10` (runner-up within 10 on
`ms`; primary-on-family guard on, as shipped) recover them?

Correction recorded here: r1778's second diagnostic called 7 of the 33 "exact AS ties". That compared the best locus with
the ORIGIN locus. In those 7 the origin is the best-AS locus and the primary alignment sits 1-5 AS lower (minimap2 ranks
on ms). None of the 33 is an exact tie of the best two loci; nothing is dropped by a defect.

## Prior rows checked
r1509 (`--tie-ms` panel: assigned 1,265 vs 78; post hoc truth 284 / 2 vs primary 280 / 6), r1514 (held out with the guard
129 / 22 vs 125 / 26), r1542 (guard on under the flag), r1778. Those truths used the gorilla individual's haplotypes with
the stage's catalogue loci; this is the catalogue-free truth on HG002.

## Design
`copy_assign` exactly as the driver's assign stage ran it in r1778 (same BAM, families table, copies FASTA, regions, skip
list, `--tied-prefilter`), plus `--tie-ms 10`. Same table, same scorer, same arms.

## Decision rule (as r1777 / r1778; primary 639 right / 71 wrong on the 710)
- PASS: assign-or-primary wrong <= 47 AND own calls wrong <= 5% AND guard S1 assign-or-primary wrong <= primary + 3.
- FAIL: assign-or-primary wrong >= 71, or own calls wrong > 10% (read at >= 30 calls).
- Otherwise PARTIAL.
- Reported: how many of the 33 get a row and what the row says; calls on the outside / tie / no-hit classes.

## Prediction
FAIL. Most of the 33 get a row; calls on the 710 rise to between 20 and 80; assign-or-primary wrong stays between 62 and
80; the newly gated rows mostly end ambiguous with margin 0, as the 24 did.

## Limits
As r1778: regional families, one sample. `--tie-ms` is opt-in; a PASS here would be evidence for the user's decision on
the default, not a default change.
