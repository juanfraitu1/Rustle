# PREREG: Winnowmap on its own target, as a positive control for the page (row 1938, 2026-11-24)

## Current state
Registered, not run.

User (2026-11-24): "make it clear why winnowmap works for certain conditions but not for multi-copy gene families".
Every Winnowmap number we have is on transcript reads at gene families (rows 1775, 1782, 1784). That it works on its own
target is so far taken from its paper. This row measures it here, so the page can show both conditions from our own runs.
Written before any read was simulated.

## Setup
- Reference CHM13 v2.0. Simulated DNA reads (no real DNA reads of these samples exist): a window of the reference, read
  from either strand, with 0.5% errors (half substitutions, a quarter 1-base insertions, a quarter 1-base deletions),
  seed 1938. 2,000 reads per class:
  - SAT10: 10 kb, from the 21 alpha-satellite HOR arrays longer than 1 Mb (`chm13v2.0_censat_v2.1.bed`).
  - FAM10: 10 kb, centred in the minimap2 primary locus of a near-tie read of the HG002 table (row 1785), one per read.
  - FAM2: 2 kb, same centres (a DNA read as short as a transcript read).
- Aligners: `minimap2 -ax map-hifi`; Winnowmap 2.03 (unpatched build of row 1775) `-W repetitive_k15.txt -ax map-pb`
  (CHM13 list of row 1782). Other options default.
- A read is RIGHT when its primary record lies on the true chromosome and overlaps the true window by at least half the
  read length; unmapped or elsewhere = WRONG.

## Reading
- WORKS ON ITS TARGET: on SAT10 minimap2 is wrong on >= 1% of reads and Winnowmap on at most half as many.
- NOT REPRODUCED: otherwise; then the page says the satellite claim is the paper's, not ours.
- Reported beside, no bar: FAM10 and FAM2 for both aligners, and MAPQ >= 30 share of right and wrong reads.

## Predictions (mine)
SAT10: minimap2 wrong a few percent, Winnowmap clearly fewer (WORKS ON ITS TARGET). FAM10: both >= 99% right, no
difference (a 10 kb DNA read carries intron and flank differences). FAM2: both worse than FAM10 and equal to each other.

## Limits
Simulated reads with uniform errors from the reference itself: the sample equals the reference, which favours every
aligner and says nothing about real satellite variation. A control for a diagram, not a benchmark of Winnowmap.
