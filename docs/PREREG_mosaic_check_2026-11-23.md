# PREREG: what are the sequence splits the haplotypes do not confirm; the mosaic check (row 1929, 2026-11-23)

Row 1928 (KB3781): where `recurrence_test.py` reports two or more sequences at a locus, the animal's own haplotypes
confirm 5 of 54 pairs. The separating variants are substitutions with a DNA-like spectrum (not homopolymer errors,
not editing), recurrent, in neither haplotype where the reads fit best. Suggested explanation to test: mosaic
molecules, one part from one copy and the rest from a related copy (PCR-mediated recombination, or template switching),
which match no genomic copy but tie exactly between the two parents.

## Measurement (`r1929_mosaic_check/mosaic.py`; the same BAM, no new alignment)
The reads are exactly tied between the loci of their group, and minimap2 `-Y` keeps the sequence in every record. So
for each read of a non-reference sequence the base-by-base state (match / mismatch) is read at the locus the tool read
it at (L1) and at every other locus of the group it has a tied record at. Per read, against the other locus that
explains most: A = read positions mismatching at L1 and matching at L2; B = matching at L1 and mismatching at L2;
C = mismatching at both. Sequences taken: every kept non-reference sequence (>= 3 reads) at a locus x isoform with
>= 2 kept sequences on KB3781 (the splits); as a comparison, the same on gorilla OR6737 (testis; row 1926's run).

## Claims and what each outcome would mean (written before any of it is computed)
- K1: >= 0.80 of the mismatches at L1 (A + C, pooled over the reads of the split sequences) are matches at another
  tied locus (A). No = the variants do not come from the tied relatives; the mosaic idea is wrong.
- K2: among reads with >= 2 A positions and >= 2 B positions, >= 0.80 have ONE switch (all A on one side of all B).
  Yes = two-part molecules; no (interleaved) = the sequence of a third copy that shares variants with both.
- K3: in at least half of the split sequences the reads have >= 3 distinct start-end pairs. Few distinct ends = copies
  of one amplified molecule (a PCR chimera multiplied; the fix is to count distinct molecules as support, as the
  assembler's read floor already does). Many distinct ends with one shared switch = independent molecules of one real
  recurrent sequence that the haplotype placement does not show (not an artifact of amplification).

## Predictions (mine)
K1 holds (the spectrum is that of related copies). K2 holds. K3: I expect few distinct ends, i.e. amplified chimeras,
but with low confidence: the KB3781 library's duplication rate is unknown to me.

## Result (workspace `r1929_mosaic_check`, `mosaic.py`)

| | KB3781 (51 split sequences, 495 reads) | OR6737 (63 split sequences, 685 reads) |
|---|---|---|
| K1 mismatches at L1 that match another tied locus | 686 of 2,430 (0.28) | 839 of 3,728 (0.23) |
| K2 reads with one switch | 68 of 146 | 3 of 114 |
| K3 sequences with >= 3 distinct start-end pairs | 42 of 51 | 61 of 63 |

- **K1 FAILS, K2 FAILS: the mosaic explanation is wrong.** Three quarters of the positions where these reads differ
  from the locus they were read at differ from the other tied loci too. They are not one part of copy A joined to a
  part of copy B; they carry a sequence that is equally far from all the tied loci at those positions.
- **K3 holds: they are independent molecules** (different ends), not copies of one amplified molecule. So not PCR
  duplicates of a chimera either. My prediction for K1 and K2 was wrong, and for K3 as well.

## Follow-up, post hoc (`in_dna.py`, KB3781): are the variants in the animal's DNA
For every read of a non-reference sequence, each of its variant positions was looked up in the same read's best hit on
the maternal and the paternal haplotype (match there = the variant is in that haplotype).

| non-reference sequences | variants in the DNA (>= 0.8 of read-sites) | mixed | NOT in the DNA (< 0.2) |
|---|---|---|---|
| at loci with >= 2 sequences (the splits) | 12 | 3 | 29 |
| at loci with one sequence | 35 | 11 | 18 |

- More splits are real than row 1928's pair test showed (12 of 44 carry variants that are in a haplotype), but two
  thirds are not in either haplotype, and neither are 18 of the 64 single non-reference sequences.
- So on this sample about half of the recurrent non-reference patterns (47 of 108) are confirmed in the animal's
  assembled DNA and about half are not: recurrent in independent molecules, substitutions with the spectrum of related
  copies, tied between reference loci, and absent from both haplotypes at the place the read fits best.
- **What they are is still open, with the mosaic and the duplicate explanations removed.** Remaining: reads of a
  sequence the assemblies do not hold where the read lands (a copy collapsed or missing in a haplotype assembly: the
  O3 question on the T2T individual itself); or the read's best haplotype hit is not its place of origin (the origin
  is another copy that minimap2 scores lower). The decisive data is the animal's DNA reads: are these variants in raw
  genomic HiFi reads of KB3781. Not done.
- For the opt-in column the reading of row 1928 stands: a non-reference sequence on real data is a recurrent pattern
  of independent molecules, confirmed in the DNA about half of the time on the one sample with truth.
