# PREREG: can a locus that only secondary alignments support be real? (row 1971, 2026-11-26)

The user: loci that only secondary alignments could discover and support; check whether any of them could be real.
Known: on chromosomes such loci were shadows by haplotype truth (rows 1767, 1783 era; KB3781, HG002), scaffold rDNA
copies were real; transcripts with no primary alignment are 0 of 361 exact at 60 planted reads per isoform and 22-24 of
about 50 at 3 reads (row 1952). Not measured: the question per LOCUS, on real samples, with independent evidence.

## Definition (`bench/o2_truth/secondary_only.py`)
A locus of the assembly (gene_id, its merged exons). PRIMARY-SUPPORTED: at least one primary alignment of any mapping
quality on >= 50 of its exon bases. SECONDARY-ONLY: none, and secondary alignments of at least 2 distinct reads.
A read whose primary was drawn elsewhere by the tie-break is the only way such a locus can be the true origin of a read.

## Evidence
1. **Planted truth** (row 1952's sets, seeds 1967-1969: E and S at 60 reads per isoform, E at 3): a secondary-only locus
   is on an EXPRESSED copy (real), a SILENT copy, or no planted copy.
2. **A second animal** (OR6737's loci looked up in KB3781's reads, same reference; and the reverse): the share of
   secondary-only loci that have a primary alignment with mapping quality >= 1 in the other animal, against the same
   share for the primary-supported loci of the same assembly. Validation only; no rule uses it.
3. **Read-origin truth** (HG002, family loci, row 1785's sample): secondary-only loci with a read whose true origin is
   the locus.
Also reported: annotated gene under the locus, spliced or not, reads per locus.

## Claims
- **S1:** planted, 60 reads per isoform: at most 2% of the secondary-only loci are on an expressed copy; at 3 reads per
  isoform at least 20% are.
- **S2:** in both gorillas the secondary-only loci have a primary alignment in the other animal at less than half the
  rate of the primary-supported loci.
- **S3:** HG002: under 5% of the secondary-only family loci have a read that originates there.
## Predictions (mine)
All three hold: real at very low depth (every read of a thinly expressed copy can be drawn away), almost never at depth.
The simple reason, to be shown with the numbers: for a copy tied with one other, the chance that all k of its reads are
drawn away is 2^-k.

## Result (jobs 4510229-4510232)
### Planted truth (seeds 1967 / 1968 / 1969): secondary-only loci, by what lies under them
| set | secondary-only loci | on an expressed copy | on a silent copy | on no planted copy |
|---|---|---|---|---|
| all expressed, 60 reads per isoform | 5 / 7 / 7 | 0 / 0 / 0 | 0 / 0 / 0 | 5 / 7 / 7 |
| half silent, 60 reads | 17 / 24 / 20 | 0 / 0 / 0 | 10 / 14 / 11 | 7 / 10 / 9 |
| all expressed, 3 reads | 14 / 13 / 13 | 2 / 1 / 0 | 0 / 0 / 0 | 12 / 12 / 13 |
### Real samples
| | OR6737 | KB3781 | HG002, family loci |
|---|---|---|---|
| loci: primary-supported / secondary-only | 20,558 / 732 | 17,683 / 351 | 1,780 / 328 |
| secondary-only: spliced model / median reads with a secondary there | 84% / 31 | 74% / 322 | 43% / 164 |
| secondary-only on a protein-coding gene / on no annotated exon | 43% / 41% | 57% / 21% | 11% / 53% |
| with a placed read (mapping quality >= 1) in the other gorilla: primary-supported / secondary-only | 79.9% / 12.3% | 95.9% / 43.6% | |
| the same, 5 or more placed reads: secondary-only | 9 of 732 | 38 of 351 | |
HG002 read-origin truth (a sample of family reads): of the 328 secondary-only family loci, 21 are touched by the sample,
and 4 of those have reads that ORIGINATE there (41, 29, 12 and 1 reads: chr1:144,556,557, chr17:37,368,245,
chr17:37,384,266 in the TBC1D3 region, chr19:23,994,803); 17 appear only as another candidate. Of the 720 sampled
primary-supported loci 414 have an origin read.

### Claims
- **S1 holds in its first half, FAILS in its second:** at 60 reads per isoform 0 of 80 secondary-only loci are on an
  expressed copy; at 3 reads 3 of 40 (7.5%; I had said at least 20%).
- **S2 HOLDS:** 12.3% against 79.9% (OR6737), 43.6% against 95.9% (KB3781; just under half).
- **S3 HOLDS as registered** (4 of 328 = 1.2%), but the denominator that can be judged is 21, and 4 of 21 is 19%.

### Reading
- In the setting where the sample IS the reference (planted), a secondary-only locus is never an expressed copy at
  depth and rarely at 3 reads: the 2^-k argument holds.
- In a real individual it CAN be real, and not by the draw: at three HG002 loci 12-41 sampled reads truly originate
  there and not one primary lands on them. A fair draw cannot lose 41 of 41; the reads fit ANOTHER reference copy
  better than the reference's copy at their own position, because the individual's copy is not the reference's (two of
  the three are in the TBC1D3 arrays, where copy content differs between people). The primary is then wrong
  systematically, and only the secondary alignment points at the true place.
- So: most secondary-only loci are relatives' shadows (the other animal does not place reads there: 12% against 80%),
  and a minority are real copies whose reads the reference pulls elsewhere. Read counts do not tell them apart (the real
  ones have 15 to 5,000 secondary reads). What tells them apart here is the individual's own assembly, which a de novo
  run does not have; the stage that looks for this without it is the missing-copy flag (O3).
- KB3781's 38 secondary-only loci with 5 or more placed reads in OR6737 are loci a gorilla does express; whether KB3781
  does is not shown.
