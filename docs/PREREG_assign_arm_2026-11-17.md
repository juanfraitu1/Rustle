# PREREG r1778: our assign stage as an arm on the neutral read-origin truth, HG002 (2026-11-17)

Written before the families and assign stages were run on this input; no assign output exists. Workspace
`r1778_assign_arm`, scorer `bench/o2_truth/assign_arm.py`.

## Question
The advisor's bar for the thesis: show the method works at multi-copy gene families. The family gives the candidate copies,
the sites where the copies differ (PSVs) decide. Does the default pipeline's read assignment (`copy_assign --families` on the
families copy table) place near-tied reads better than minimap2's primary alignment, judged by the catalogue-free truth of
r1775 / r1776?

## Prior rows checked
r1776 (the room: 71 of 710 decidable near-tie reads have a primary alignment that is not the origin), r1777 (graphs with the
sample's own uncertified alleles FAIL), r1374 (own variant calls FAIL held-out), r1699 (PSV-assigned 0.56% of the A119b
near-tie pool), o2 notes (moved rows wrong 59-78% under the catalogue truths). New here: the assign stage has never been
scored by a truth that does not use its catalogue. HG002 is spent; nothing is tuned.

## Design
- Input: `r1777_allele_graph/region.bam` (every HG002 record on the 501 candidate regions) and the assembly already made
  from it with current defaults (`asm/x.*`, copied). Stages `families` then `assign` of `tools/rustle_pipeline.sh`,
  defaults, no switches.
- Table: the r1777 table (`allele_graph.table`): all near-tie reads (S2) + the first 5,000 clear-lead reads (S1, guard) of
  `r1776_hg002_room/hg/reads.tsv`.
- A read's assign call = its row with status `assigned`: the placement of the read on the same chromosome and strand that
  overlaps the assigned copy's locus. An assigned copy that overlaps none of the read's placements is a call "off the
  placements" and counts as wrong for a decidable read (its origin is one of the placements). Any other status or no row =
  abstain.
- Arms: primary; assign (its own calls); assign-or-primary (the call, else the primary alignment).

## Decision rule (the bars of r1777, on the 710 decidable near-tie reads; primary 639 right / 71 wrong)
- PASS: assign-or-primary wrong <= 47 AND assign's own calls wrong <= 5% AND guard S1: assign-or-primary wrong <= primary + 3.
- FAIL: assign-or-primary wrong >= 71, or assign's own calls wrong > 10%.
- Otherwise PARTIAL. Own-call share is read only at >= 30 calls; below that the verdict rests on the net count alone and
  is stated as "abstains".
- Reported, not judged: status of every table read (no row, tied, ambiguous, assigned); calls on the outside / tie /
  no-hit classes; calls off the placements.

## Prediction
FAIL by "no net gain": fewer than 100 calls on the 710, assign-or-primary wrong between 65 and 80, most table reads with no
row or not assigned. Guard unchanged within 3.

## Limits stated beforehand
- Families and copies come from the 501 regions only, not a whole-sample run (147 GB BAM does not fit the 1 h 30 limit);
  a whole-sample run can give different families and so different calls.
- The assign stage gates on its own tie definition (ms near ties on the families' copies), not on r1776's AS lead; reads it
  never considers are abstentions here.
- One sample, 25% hash sample of family-gene molecules, substitutions and junctions as the stage uses them.
