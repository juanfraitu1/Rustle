# PREREG r1791: what the "true ties" of r1790 are (2026-11-18)

Written before any count of this row. Workspace `r1791_tie_split` (`tie_split.py`, `split.sbatch`). Inputs: the read and
placement sequences of r1784 (`r1784_winnowmap2_spliced/hg/{reads,placements}.fa`, `reads.tsv`), HG002 haplotype indexes of
r1776. The scorer is not changed; its functions are imported and its rule is re-applied to the tied reads only.

## Why
User, 2026-11-18: are the true ties undecidable even with PSVs? In `neutral_truth.py` a near-tied read has status `tie`
when more than one candidate locus of the reference coincides with a best hit of the read on HG002's haplotypes. That
covers two different things, and only the first means "no sequence difference in the read can place it":
- SEVERAL COPIES: the read has equally good hits at two or more places of one HG002 haplotype, and the candidate loci
  land on different ones of them.
- ONE PLACE: on every haplotype all the read's origin loci land on one and the same hit of the read. The reference loci
  are not kept apart on HG002 (one copy fewer, or copies within the 2-mismatch slack of the locus test). This is a limit
  of the truth, not evidence about the read.

## Rule, fixed now
For a tied read: per haplotype, D = the number of distinct best read hits that an origin locus coincides with
(`same_place`, the old test, as r1784). SEVERAL COPIES when D >= 2 on at least one haplotype, else ONE PLACE.
Check: the status recomputed here must be `tie` for every read that r1784 called `tie` (the job fails otherwise).

## Readings, fixed now
1. The sentence "59% of the near ties in the 14 families cannot be placed by any method" STANDS when at least 80% of
   their true ties are SEVERAL COPIES. Otherwise it is replaced by the SEVERAL COPIES share, and the ONE PLACE reads are
   reported as "not judged by the truth".
2. Reported, not judged: the split per family (14 families and all), and for ONE PLACE reads whether the read has
   further best hits on that haplotype that no candidate locus reaches.

## Prediction
STANDS: 80 to 95% SEVERAL COPIES overall; EIF3C / EIF3CL, FAM156 and SMN above 90%.

## Limits stated beforehand
- SEVERAL COPIES means the read fits two HG002 copies with the same mismatch count; a PSV cancelled by a read error is
  counted here. HG002's two haplotypes are treated separately; a read that fits one copy on each haplotype is ONE PLACE.
- One sample, 25% of the family reads; the known defects of the truth rule apply.
