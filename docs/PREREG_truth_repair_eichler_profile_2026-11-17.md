# PREREG r1785 (repair of the neutral truth) and r1786 (how far the adapted Eichler filter carries) (2026-11-17)

> Naming corrected 2026-11-18: "Eichler's rule" here always meant our adaptation to one reference of the Eichler lab's 10-point filter (Guitart 2024); the paper's own procedure is different and was not tested. Definition: `docs/copy_assignment_definition.md`. Rules, bars and predictions are unchanged.

Written before the repaired truth was run and before any per-lead count was looked at. Workspaces `r1785_truth_repair`,
`r1786_eichler_profile`; code `bench/o2_truth/neutral_truth.py` (switch `NT_REPAIR=1`), `bench/o2_truth/eichler_profile.py`.

## r1785: the repair
Defect (docs/NEUTRAL_TRUTH_PLACEMENT_2026-11-17.md): the origin test was a span overlap without strand or exon check, and
a placement with no near-perfect haplotype hit was dropped as a candidate. Known before this run: the 8 chr1 reads of
r1783. Nothing else about the repaired table is known.

Under `NT_REPAIR=1` (the old behaviour stays the default until the user decides, because another session is scoring with
this script today):
1. A placement is an origin when one of its haplotype hits and the read's best hit are on the same haplotype, contig and
   strand and share at least one exon base (aligned blocks, introns excluded).
2. Competitors of a read = its placements within 10 AS of the best, plus the primary's placement. If a competitor has no
   near-perfect haplotype hit, the read's status is `competitor unlocatable` and it is not scored.
Everything else is unchanged (slack, near-perfect bound, strata, caps, hash order).

Re-scored on the repaired tables, with the bars of their own rows and no new tuning: primary, both adaptations of the Eichler filter,
Winnowmap (KB3781), r1777 allele graph, r1778 assign, r1781 `--tie-ms 10`, r1783 own alleles. Reported next to the old
counts. The old verdicts are not rewritten; each row gets a line with its repaired counts.

Reading of the repair itself:
- The 8 chr1 reads must leave the "primary right" cell (become primary-wrong or unlocatable): if they do not, the repair
  did not do what it was written for and nothing else is read.
- Reported: how many of the 710 (HG002) and 62 (KB3781) decidable near ties change status or origin.

Prediction: 5 to 25% of the 710 change; most changes are to `competitor unlocatable`; the room (primary wrong) stays
between 50 and 110; no rule passes its bar on the repaired table, r1783 comes closest.

## r1786: the adapted Eichler filter as a heuristic (user: "useful, but things are more complex; can we measure that?")
The rule: an alignment more than 10 AS above the rest gets the read; within 10 the read is undecided. Measured on the
repaired tables of both samples (HG002 on CHM13; KB3781 on its own T2T assembly), decidable reads, S1 sample weighted to
the full stratum.

M1 reach (no truth): share of family-gene molecules with one placement, a clear lead, a near tie.
M2 lead profile: by lead bin (0, 1-5, 6-10, 11-20, 21-50, 51-200, > 200): reads, share where the best-AS locus is the
   origin, share where the primary is the origin. A bin is read at n >= 30.
M3 score dependence (count without truth, accuracy with): reads with a lead > 10 whose best-AS locus is not the primary;
   who is right there.
M4 size dependence (no truth): near-tie share by the best alignment's AS (< 1000, 1000-1999, 2000-3999, >= 4000), a
   fixed 10 against reads of different length.

Reading rules, fixed now:
- "Above 10 is safe": holds if the best-AS locus is wrong in <= 1% of the bins 11-20 and 21-50 taken together; fails at > 5%.
- "10 is a threshold": holds if the best-AS wrong share in 6-10 is at least 3 times that in 11-20; otherwise the error
  falls gradually and 10 is not special.
- "Within 10 is undecided": holds if the primary is right in < 80% of the bin 6-10; if >= 90%, the rule discards
  placements that are mostly right.
- M3 and M4 are reported, not judged.

Predictions: above 10 is NOT safe as "move to the best AS" (the errors are the reads where best AS and primary disagree)
but is safe as "keep the primary"; the error falls gradually (no threshold at 10); the primary is right in >= 90% of the
bin 6-10; near ties are more frequent among short reads.

Limits: two individuals, one of them the reference's own individual (KB3781); family-gene molecules only; the truth is the
repaired one, itself new; exact ties (lead 0) have no defined best locus, the first in the table is used.
