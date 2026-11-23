# PREREG: the sequence count on real reads, judged by the sample's own two haplotypes (row 1928, 2026-11-23)

Rows 1926 / 1927: IsoCon's test ported to Rustle's tied groups counts the distinct sequences among copies the reference
cannot separate (0.80 exact on planted homozygous copies; 0.53 with two alleles per copy; never tested on real errors).
User: test it where the DNA answer is known. KB3781 (gorilla fibroblast) is the assembled individual with maternal and
paternal haplotype assemblies (`o2contested/haps/KB3781.{mat,pat}.mmi`).

## Written before any haplotype-based number exists
The pipeline job (4496969: assemble with the good pool, families, `o2` with `RUSTLE_O2_SEQUENCES=1`) prints step-0 counts
only (sequences, loci x isoforms, how many with >= 2 sequences). Nothing below has been computed.

## Truth from the haplotypes (`r1928_real_haplotypes/hap_score.py`)
Every read of a kept sequence is aligned to both haplotypes (same minimap2 recipe, secondaries kept). Per haplotype:
cost = edit distance + clipped bases of its best hits, and the regions of all hits at that cost (regions = overlapping
hit spans merged). Reads whose regions overlap form one COPY of the sample (copies identical in both haplotypes fall
together: one sequence, correctly). A read prefers the maternal or paternal haplotype when its cost is lower there;
equal cost = not informative about the allele. A copy shows TWO sequences when >= 3 of its reads prefer each haplotype,
else one. A sequencing error costs the same on both haplotypes, so it cannot create a preference.

## Measurements and claims (kept sequences with >= 3 reads; loci x isoforms with a truth copy)
- R1 purity: a kept sequence is ONE real sequence when >= 0.9 of its reads belong to one copy and, if that copy shows two
  sequences, >= 0.9 of its allele-informative reads prefer one haplotype. Claim: >= 0.80 of the kept sequences.
- R2 count: the tool's number of sequences at a locus x isoform equals the number the haplotypes show. Claim: >= 0.60
  (planted: 0.80 homozygous, 0.53 heterozygous).
- R3 false splits: among pairs of kept sequences at one locus x isoform, pairs on the same copy and the same
  haplotype side (the haplotypes show one sequence where the tool reports two). Claim: <= 0.15 of the pairs.
- Descriptive, the point of the row: of the remaining pairs, how many are two ALLELES of one copy (same copy, opposite
  haplotypes) and how many are different COPIES. This is the label the RNA reads cannot give (row 1927).
All of R1-R3 = the count transfers to real HiFi reads on this sample; then HG002 on CHM13 (both haplotypes on disk,
sample differs from the reference) is run with the same script as the second sample.

## Predictions (mine)
R1 holds. R2 near 0.6, either side. R3 is the real risk: HiFi errors recur at homopolymers, IsoCon's error model
allows for that only through the read's own indel count, and I expect 10-30% false splits. Because KB3781 is the
reference individual most non-reference sequences should be the other haplotype, so alleles should outnumber copies
among the true pairs. Fibroblast expresses few multi-copy families: the step-0 count may be small, and under 30 loci x
isoforms with >= 2 sequences the claims are reported as counts without a verdict.

## What happened first (job 4496969)
The pipeline ran (assemble 1,354 s, families, `o2` to 3,665 s; 17,684 loci labelled, 225 families) but the sequence
test crashed: `OverflowError`, IsoCon's placement factor is an exact integer and exceeds any float when two candidates
differ at hundreds of sites (planted copies never do). The stage reported the failure and went on, as designed. Fixed
(p-value and factor combined in logarithms; self-test with 400 sites) and the sequence step rerun on the stage's files
(job 4497034, 272 s). No haplotype-based number existed before the fix.

## Result (KB3781; `hap_score.py`)
Step 0: 322 sequences at 272 loci x isoforms; 30 of them with >= 2 sequences (the registered minimum for a verdict);
214 sequences with >= 3 reads; 2,837 tied reads, all with a hit on the maternal haplotype, 1,812 on the paternal.

| measure | result | claim |
|---|---|---|
| R1 purity: kept sequences that are one real sequence | 175 of 178 (0.98) | >= 0.80: HOLDS |
| R2 count equals the haplotypes' sequences | 126 of 147 loci x isoforms (0.86); over 18, under 3 | >= 0.60: HOLDS |
| R3 pairs of sequences the haplotypes do not separate | 49 of 54 (0.91) | <= 0.15: **FAILS** |
| of the other pairs | two alleles of one copy 2, different copies 3 | descriptive |

- **Where the tool reports one sequence it agrees with the haplotypes (most loci); where it reports two or more, the
  haplotypes confirm only 5 of 54 pairs.** The count transfers for "one"; the splits do not.
- Post hoc, to see what the unconfirmed sequences are (not registered):
  - A second pair test by each sequence's majority placement on the haplotypes: 11 of 58 pairs at different places
    (copies), 47 at the same place and side. The first measure was not simply too strict.
  - The separating variants are NOT the error type I predicted: 954 substitutions against 46 insertions over the kept
    sequences, 70% transitions, A>G 19% on the transcript strand. Not homopolymer indel errors, not RNA editing: the
    spectrum of DNA differences between related sequences.
  - A crude cost check: at loci with >= 2 sequences, 44 of 51 non-reference patterns are not explained by either
    haplotype (median cost of their reads on the best haplotype hit >= their number of variants).
  So: recurrent, strongly supported (3-55 reads, p down to 1e-92) substitution patterns that look like paralogue
  differences and are in neither haplotype of the same animal at the place the reads fit best. Candidates, untested:
  mosaic molecules made in library preparation (PCR-mediated recombination between transcripts of related copies);
  reads of a copy that is neither among the tied loci nor well represented where they place; a truth that cannot see
  the difference (cost ties between places). I do not know which.
- Predictions: R1 right, R2 better than expected (0.86), R3 fails as feared but for a different reason than the one I
  gave (homopolymer errors).
- **HG002 is not run** (the registered condition was all three claims).
- **For the opt-in column:** on real data a locus with `tied_sequences` of 2 or more is unverified; do not present such
  splits as expressed copies or alleles. The next step is to find out what these sequences are (compare their variant
  alleles with the other members of the family at the homologous positions: a mosaic would match them site by site).

**Correction 2026-11-23 (row 1931):** the pair measure R3 (copies by overlapping hit regions, allele by cost) is too
coarse: on HG002 one locus with 42 sequences gives 861 of its 1,318 pairs, and a site-by-site lookup of the variants in
the haplotypes (rows 1929-1931) contradicts it. Do not quote "49 of 54 pairs unconfirmed"; quote 52 of 107 sequences
in the DNA (KB3781) and 517 of 689 (HG002).
