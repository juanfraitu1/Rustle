# PREREG r1783: copy differences from the individual's own reads, certified against the other copy, HG002 (2026-11-17)

Written before any count was looked at. Workspace `r1783_own_psv`, scorer `bench/o2_truth/own_psv.py`.

## Question (user, 2026-11-17)
The family gives the options; PSVs tell copies apart. Reference PSVs do not (r1778: 17 of 24 gated reads tie on them),
uncertified sample alleles mislead (r1777). Do differences taken from the individual's own confidently placed reads, used
only where the OTHER copy's own reads do not show the same base, place near-tied reads better than the primary alignment?

## Prior rows checked, and why this is not a reopening
- r1333 (sample's own copies by read clustering, HG002): FAIL at 4 wrong / 86 placed, 126 / 11 vs 48 / 7. Its named, never
  run, fix: teach only with reads resolved by a clear margin. Here the teaching reads are MAPQ-60 primaries only.
- r1374 (sample's variant calls): showcase PASS, held-out FAIL, errors in CNV families.
- r1777 (alleles forgiven per locus, no check against the other copy): FAIL.
- New: (1) the pairwise certificate below, (2) anchors restricted to MAPQ 60, (3) the catalogue-free truth of r1776.
HG002 is spent; nothing is tuned: the thresholds are r1777's (>= 3 reads, >= 0.2 of anchor coverage) and the stage's
"two deciding sites".

## Rule
- Anchors: primary alignments with MAPQ 60 in `r1777_allele_graph/region.bam`. Own alleles of a reference position = bases
  with >= 3 anchor reads and >= 0.2 of the anchor coverage there (the reference base counts only if the anchors show it).
- Competitors of a read: its placements within 10 AS of its best, plus the primary's placement. For a clear-lead read
  (guard) the best two placements.
- For two competitors X, Y, a read position aligned (match or mismatch column) in both is a CERTIFIED site when both
  positions have own alleles and the two allele sets are disjoint. The read votes X if its base is an own allele of X,
  Y if of Y.
- Call X when, against every other competitor Y, X has >= 2 votes and Y has 0. Otherwise abstain.
- Reported only: the looser call (votes X - votes Y >= 2 against every Y).
- Substitutions only; no pooling of reads across copies; no second sample.

## Decision rule (as r1777 / r1778 / r1781; 710 decidable near-tie reads, primary 639 right / 71 wrong)
- PASS: own-or-primary wrong <= 47 AND own calls wrong <= 5% AND guard S1 own-or-primary wrong <= primary + 3.
- FAIL: own-or-primary wrong >= 71, or own calls wrong > 10% (read at >= 30 calls).
- Otherwise PARTIAL.
- Reported: reads with at least one certified site; calls on the outside / tie / no-hit classes; moves off the primary.

## Prediction
PARTIAL. Calls on 40 to 150 of the 710 (copies this similar have few MAPQ-60 reads, so few certified sites); own calls
wrong between 5% and 10%; own-or-primary wrong between 55 and 70.

## Limits stated beforehand
- The truth uses HG002's haplotypes, so a rule reading HG002's own differences is close to what the truth measures; a PASS
  shows the reads carry the signal, not that it transfers to a sample without such truth.
- MAPQ-60 reads can themselves be misplaced; then the other copy's base enters the allele set and the site is dropped
  (conservative), or a wrong allele is certified (the failure to look for).
- One sample, 25% hash sample of family-gene molecules, candidate regions only.
