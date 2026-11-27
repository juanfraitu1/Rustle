# PREREG: recurring differences among tied reads as a warning for the cluster verdict (row 1978, 2026-11-26)

The cluster verdict (rows 1974 / 1975, `RUSTLE_O2_CLUSTERS=1`) is wrong for 2.7% of HG002's decided tied reads: their true
locus is in another cluster, in families where the individual's copy is not the reference's (NBPF, SRGAP2C). The user:
can errors be clustered as IsoCon does. IsoCon's test is already ported (`bench/o2_families/recurrence_test.py`,
register 1926): among the tied reads of a group of loci it finds the sequences the reads really show, a difference from
the reference that RECURS across reads being a real variant and one seen once an error.

## Idea
A group of tied reads that agree with each other on a difference from the reference carry a sequence the reference does
not have at that place: the individual's own version of a copy, or a copy from elsewhere. Such reads are where a verdict
taken on the reference is most likely wrong. So: flag the decided reads that belong to a NON-REFERENCE sequence (a
candidate of the test with at least one allele that differs from the reference).

## Data and method (`bench/o2_truth/cluster_sequences.py`)
HG002 family reads (the sorted BAM, tie table, clusters and verdicts of row 1975's built step, workspace
`r1975_cluster_filter_denovo/hg_built`), truth = row 1972's haplotype run, any assembled locus on the truth gene
accepted. `recurrence_test.py` at its defaults on the same BAM and GTF. Each decided read with a truth is: in a
non-reference sequence, in a reference sequence, or not tested (not in an exactly tied group, or a sequence with fewer
than 2 reads, or ambiguous between sequences).

## Claims
- **F1:** decided reads in a non-reference sequence are wrong at least 3 times as often as decided reads in a reference
  sequence.
- **F2:** flagging the reads of non-reference sequences catches at least 50% of the wrong decided reads that the test
  covers, and flags at most 20% of the right ones it covers.
- **F3 (reach):** the test covers at least 50% of the decided reads that have a truth.
## Predictions (mine)
F1 holds. F2 uncertain: heterozygous sites of a correctly placed copy also give non-reference sequences (alleles, not
copies), so right reads will be flagged too; I expect the 20% bound to be the one at risk. F3 fails: the test needs
exactly tied groups and depth, and many decided reads are near ties that are not exact.
Not tested here: assigning a sequence to a member of its cluster by matching it to that member's bases (the test
reports alleles against one locus only; it would need the members aligned to each other).

## Result (job output `r1978_cluster_sequences/run.*.out`)
The recurrence test on HG002's tied family reads: 59 groups, 168 sequences kept, 128 of them non-reference, 40 reference.

| decided reads with a truth (5,257; 140 wrong, 2.7%) | reads | wrong | share wrong |
|---|---|---|---|
| in a non-reference sequence | 833 | 21 | 2.5% |
| in a reference sequence | 134 | 2 | 1.5% |
| not tested | 4,290 | 117 | 2.7% |
Flagging the non-reference sequences: catches 21 of the 23 wrong reads the test covers (91%), and flags 812 of the 944
right ones (86%). By family (flagged reads, wrong among them): EIF3C 550 (0), FAM156A 140 (0), CCZ1B 66 (0), NPIPB5 17
(13), RGPD8 9 (0), GTF2IRD2B 7 (0), RASA4 6 (0), NPIPB3 6 (6).

### Claims
- **F1 FAILS:** 2.5% against 1.5% wrong, a ratio of 1.7 (bar 3). I predicted it would hold.
- **F2 FAILS:** it catches 91% of the covered wrong reads only because it flags 86% of the right ones (bar 20%).
- **F3 FAILS** (as predicted): the test covers 18.4% of the decided reads with a truth.

### Reading
- **A recurring difference from the reference is the normal case, not a warning.** Three quarters of the sequences the
  tied reads show are non-reference (128 of 168): a person's expressed copies differ from another person's reference by
  alleles and private variants almost everywhere. EIF3C alone has 550 flagged reads and none of them is wrong.
- **So the flag does not separate the errors.** The wrong verdicts are concentrated in one family (NPIP: 19 of the 21
  flagged wrong reads), not in "non-reference" as such. What marks them is the family, which a run cannot know in advance.
- **Reach is the second limit:** four decided reads in five are near ties that are not exact, or too few per sequence,
  and the test does not look at them.
- IsoCon's test keeps its validated use (counting the sequences expressed among copies the reference cannot separate,
  row 1926). As a warning for the cluster verdict it is not usable. Not tried, and the only variant that could still
  work: comparing each sequence with every member of its cluster, to tell "matches another member" (a wrong place) from
  "matches none" (an allele); it needs the members aligned to each other.
