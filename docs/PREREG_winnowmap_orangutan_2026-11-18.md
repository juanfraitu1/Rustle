# PREREG r1787: Winnowmap (stock and patched) on a second ape, with tandem copies read on their own (2026-11-18)

Written before any count of this row. Workspace `r1787_winnowmap_orangutan` (`orang.sbatch`, `wm_ties.py`,
`tandem_split.py`), scorer `bench/o2_truth/neutral_truth.py` as it stands today.

## Why
User, 2026-11-17 / 18: the pipeline excludes satellites, so tandem gene copies are the place Winnowmap could still help;
the HG002 split (r1784, post hoc) had 85 decidable tandem reads, too few to judge. Use another ape.

## Sample, reads, truth
- Sumatran orangutan AG06213 fibroblast Iso-Seq (`reads/PAB/PAB_AG06213_fib_mm.bam`), reference `genomes/PAB/PAB.fa`
  (RefSeq mPonAbe1). Chosen because both haplotype assemblies of this animal are public (GenomeArk
  `mPonAbe1.hap1.cur.20231205`, `hap2`) and orangutans are the most heterozygous great ape on hand. Not a tuning sample.
- Reads = ALL molecules (no hash sample) with a primary on a member of a truth family: Compara simiiformes pairs on the
  RefSeq annotation, genes mapped to human by best protein hit (the r1611 recipe, built in the job).
- Arms: minimap2 primary; stock Winnowmap 2.03 (r1782); patched Winnowmap with the substring method on (r1784 build,
  unchanged). `-W` list from meryl on PAB.fa. All three scored in one call, so on one decidable set.
- Truth rule as r1775. Strata by minimap2's records. Tandem = the primary's locus and its best-AS rival locus on the same
  chromosome and strand within 1 Mb (the r1784 definition, unchanged).

## Readings, fixed now
0. Each Winnowmap arm is judged only if >= 0.95 of reads map and its spliced share is within 0.05 of minimap2's.
1. All near ties, per arm: GAIN when right share on decidable S2 reads (>= 200) is >= 1 percentage point above the
   primary's AND wrong on decidable S1 reads <= primary's + 3. Otherwise NO GAIN.
2. Tandem near ties, per arm: GAIN when, with >= 200 decidable tandem S2 reads, the arm's right share is >= 1 point
   above the primary's. Under 200 decidable: NOT READABLE, counts reported with fixed / broken.
3. Reported: the class split, ties resolved on each arm's own scores, same-primary-locus share.

## Predictions
The reference is one haplotype of the same animal (as KB3781, where the primary was right on 57 of 62), so the primary
is right on >= 90% of decidable near ties and the room is small. NO GAIN for both arms; the patched arm breaks more
than it fixes. Tandem: fewer than 200 decidable (NOT READABLE) is likely; if readable, no gain.

## Limits stated beforehand
- Self-reference: reads of the reference haplotype match it exactly, so this sample is easier for any aligner than
  HG002 on CHM13. An ape that is not its own reference, with both haplotypes assembled, is not on hand.
- Fibroblast; one sample; the substring method is run outside the read lengths it was built for.
