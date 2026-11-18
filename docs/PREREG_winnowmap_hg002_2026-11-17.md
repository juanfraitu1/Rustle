# PREREG r1782: Winnowmap on the near-tied reads of HG002 against CHM13 (2026-11-17)

Written before any count of this row. Workspace `r1782_winnowmap_hg002` (`wm.sbatch`, `wm_ties.py`), scorer
`bench/o2_truth/neutral_truth.py` unchanged.

## Why
User, 2026-11-17: re-align the near-tied reads with Winnowmap and see whether that resolves the ties. r1775 ran Winnowmap
as an arm on KB3781 (62 decidable near ties, primary 57 / 5, Winnowmap 58 / 4), but KB3781 is the reference animal and
left room for 5 reads. r1776 found room on HG002 against CHM13: 710 decidable near ties, minimap2's primary wrong on 71.
This row runs the Winnowmap arm there. Each read is mapped on its own, so mapping all selected reads and reading the
near-tie stratum equals mapping the near-tied reads alone.

## Sample, reads, truth
- The reads of r1776 (`r1776_hg002_room/hg.sel.fa`, 25% hash sample of HG002 Kinnex molecules with a primary on a member
  of a truth family) and its minimap2 BAM (`mm.sel.bam`). Strata and caps as r1775 / r1776, set by minimap2's records:
  S2 = lead <= 10 AS (near tie), S1 = lead > 10.
- Winnowmap 2.03 (the r1775 build), `-W` list from `meryl` k = 15 on CHM13 (`distinct=0.9998`), `-ax splice:hq` and the
  flags of our BAMs. Known from r1775: the splice presets set `SVaware = false`, so this is the weighted-minimizer
  method, not the Winnowmap2 method.
- Truth = both haplotypes of HG002 T2T v1.1 (`r1776_hg002_room/{mat,pat}.mmi`), rule exactly as r1775.
- The scorer keeps Winnowmap's pick as a placement, so the decidable set can differ from r1776's 710 by a few reads.
  Both arms are judged on the decidable set of this run.

## Readings, fixed now
1. GAIN when, on the decidable S2 reads (>= 200 of them), Winnowmap's primary is right on at least 1 percentage point
   more than minimap2's primary (reading 4 of r1775), AND on the decidable S1 reads its wrong count is at most the
   primary's + 3 (the guard of r1777). Otherwise NO GAIN as a spliced aligner.
2. Reported, not judged: primary x Winnowmap on the decidable S2 reads (fixed / broken / both wrong); the share of reads
   with the same primary locus; and "resolved on its own scores" = S2 reads that have a lead > 10 AS in Winnowmap's
   records, with Winnowmap's right / wrong among the decidable ones.

## Predictions
NO GAIN. Same primary locus on >= 97% of reads. Winnowmap fixes fewer than 10 of the primary-wrong reads and breaks about
as many. Fewer than 5% of S2 reads get a clear lead on Winnowmap's scores.

## Limits stated beforehand
- One sample; a 25% sample of the family reads. HG002 families are spent on earlier rows; this row tunes nothing.
- The Winnowmap2 method cannot be tested on spliced reads with this build.
- Limits of the truth as r1775 / r1776 (mismatch counts; "outside" inflated by reference differences).
