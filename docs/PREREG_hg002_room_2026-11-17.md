# PREREG r1776: room for graph assignment at near ties on a sample that is not its own reference (2026-11-17)

Written before any count. Workspace `r1776_hg002_room`, script `bench/o2_truth/neutral_truth.py` (the scorer of r1775, unchanged
but for a hash sample in `select`).

## Why
User and advisor, 2026-11-17: build graphs from the produced GTF, group loci whose graphs look alike into families, and assign
reads with near-tied secondary alignments on those graphs. Earlier graph rows (r1643, r1649, r1696) equalled minimap2's primary
on KB3781, and r1775 showed the room there is 5 reads (62 decidable near ties, the primary right on 57): KB3781 is the reference
animal. This row is step 0: count the room on HG002 aligned to CHM13, before any graph is built.

## Sample, reads, truth
- HG002 Kinnex, `HG002_kinnex.t2t.bam` (CHM13). Reads = molecules with a primary record on the exons and strand of a member of a
  truth family (Compara simiiformes, CHM13 CAT annotation), a 25% hash sample (`NT_SELECT_PCT=25`), realigned with minimap2 (our flags).
- Truth = both haplotypes of the HG002 T2T assembly v1.1 (maternal and paternal, downloaded from the human-pangenomics bucket),
  indexed with `-x splice:hq`. Rule exactly as r1775: near-perfect best hit, slack 0, a placement is the origin when its
  reference sequence lands where the read does; decidable = exactly one placement.
- Strata as r1775: S0 no other locus, S1 lead > 10 in AS, S2 lead <= 10 (near tie). Caps 5,000 / 30,000 / 30,000.

## Reading, fixed now
- ROOM = decidable S2 reads on which minimap2's primary is not the origin.
- GO for the graph arm (graphs from the produced GTF carrying the sample's own alleles, scored on these reads) when ROOM >= 30.
  Below 30: STOP, the graph arm is not built on this sample.
- Reported, not judged: the primary's right share on decidable S2 and S1 reads; the tie / outside / no-hit split; reads where
  the best-AS locus is not the primary.

## Predictions
ROOM >= 30 (GO). Primary right on 75-92% of decidable S2 reads, below KB3781's 0.919. Ties remain the largest S2 class.

## Limits stated beforehand
- CHM13 differs from HG002 everywhere, so the placement's reference sequence is not expected to be near-perfect on the
  haplotypes; the rule asks only that it lands at the read's locus (hits within NM + 2 of its own best, best near-perfect).
  If most placements have no near-perfect hit, "outside" will be inflated; the outside share is reported for that reason.
- One sample; HG002 has been used for family rules, not for this count.
- The room is an upper bound for any method, not a result for graphs.
