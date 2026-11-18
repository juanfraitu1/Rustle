# PREREG r1784: the Winnowmap2 method forced on for spliced reads, HG002 near ties (2026-11-17)

Written before any count of this row (a 2,000-read smoke run showed only that the patched build runs and reports
"SV-aware mode"; no placement was counted). Workspace `r1784_winnowmap2_spliced` (`wm2.sbatch`, patched source in
`Winnowmap/`), scorer `bench/o2_truth/neutral_truth.py` unchanged, tie table `r1782_winnowmap_hg002/wm_ties.py`.

## Why
User, 2026-11-17, after r1782 (Winnowmap worse than minimap2 on the near ties, 605 / 106 against 639 / 72): check the
part r1782 could not test. The Winnowmap2 method (minimal confidently alignable substrings: map substrings of the read,
keep the ones that map with MAPQ >= 5, build the alignment from their anchors) is the part meant to choose between
copies. The stock build cannot run it on our reads for two reasons read from the source: the `splice` presets set
`SVaware = false`, and the method only applies to reads >= 10,000 bases with substrings of 2,000 to 16,000 (no read of
the 136,074 is that long; mean 2,006).

## The patch (one hunk in `src/options.c`, splice presets)
`SVaware = true`, `SVawareMinReadLength = 1000`, substrings 500 / 1,014 / 2,056 bases sampled every 500 bases
(`minPrefixLength = suffixSampleOffset = 500`, `maxPrefixLength = stage2_max_gap = 4000`). These values are a single
choice made before the run (the map-pb values cut by 2, to fit reads of about 2 kb); no second setting is run.
Reads under 1,000 bases get the weighted-minimizer method of r1782.

## Sample, reads, truth
As r1782: reads `r1776_hg002_room/hg.sel.fa`, minimap2 BAM of r1776, `-W` list of r1782, the same flags, truth = both
HG002 T2T v1.1 haplotypes. Strata by minimap2's records.

## Readings, fixed now
0. Step 0: the run is judged only if >= 0.95 of reads map and the spliced share of primaries is within 0.05 of
   minimap2's (0.897). Otherwise INVALID: the method does not carry over to spliced reads with this patch.
1. GAIN when, on the decidable S2 reads (>= 200), the patched Winnowmap's primary is right on at least 1 percentage
   point more than minimap2's primary, AND on the decidable S1 reads its wrong count is at most the primary's + 3.
   Otherwise NO GAIN.
2. Reported, not judged: primary x Winnowmap on decidable S2 (fixed / broken); same-primary-locus share; S2 reads with
   a clear lead on its own scores; the same figures of r1782 beside them.

## Predictions
Step 0 passes. NO GAIN: right share on decidable S2 between 80% and 90% (minimap2 89.9%, r1782 85.1%); the guard fails.
Reason: a substring that maps confidently to one copy carries the same differences the whole-read score already counts.

## Limits stated beforehand
- The method is run outside the read lengths it was built for; a FAIL says nothing about it on genomic reads.
- One parameter choice, one sample, a 25% read sample; truth limits of r1775 / r1776.
