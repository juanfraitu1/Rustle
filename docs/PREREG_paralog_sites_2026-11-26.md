# PREREG: are the recurring differences of a tied sequence the bases of another copy? (row 1979, 2026-11-26)

Row 1978: flagging every sequence that differs from the reference does not find the wrong cluster verdicts; a
difference from the reference is the normal case in another individual. Left untried there, and asked for by the user:
compare each sequence with the OTHER copies, to tell "these differences are the bases of another copy" (a sign the reads
belong elsewhere, or to a mixture of copies) from "these differences match no copy" (an allele or a private variant).

## Method (`bench/o2_truth/paralog_sites.py`)
For every non-reference sequence of the recurrence test (row 1978's `seq.tsv`) and each of its substitution sites: the
reads of the sequence carry a base there that differs from the locus they are read at. Each of those reads also has
other alignments in the BAM (`-N 50 -p 0.1`, `--eqx`). At an alignment to another locus Y, the same read base is either a
match or a mismatch. A site is EXPLAINED BY Y when most of the sequence's reads that reach Y match there. Loci of the
read's decided cluster are not counted as Y (they tie with the reads).
A sequence is `paralog-explained` when one locus Y explains every one of its substitution sites, `private` when no
locus explains any, `mixed` otherwise. Insertions and deletions are not used.

## Data
As row 1978: HG002 family reads, the built step's verdicts, truth = row 1972's haplotype run.

## Claims, over the decided reads with a truth that sit in non-reference sequences (833; 21 wrong)
- **H1:** reads of paralog-explained sequences are wrong at least 5 times as often as reads of private sequences.
- **H2:** flagging the paralog-explained sequences catches at least 50% of those 21 wrong reads and flags at most 20%
  of the 812 right ones.
- **H3:** for the flagged wrong reads, the explaining locus Y lies in the cluster of the truth locus for at least 50%.
## Predictions (mine)
H1 holds. H2 fails on its second half: copies of a family exchange sequence, so a correctly placed copy's own alleles
often coincide with another copy's bases. H3 uncertain. If H2 fails the idea is closed.

## Result (`r1979_paralog_sites/out.txt`; one pass over the family-read BAM, run in the session)
Non-reference sequences with substitution sites: 114; paralog-explained 74, mixed 13, private 27.

| sequence | decided reads with a truth | wrong | share wrong |
|---|---|---|---|
| paralog-explained (one other locus matches at every site) | 50 | 21 | 42.0% |
| mixed | 13 | 0 | 0% |
| private (no other locus matches at any site) | 676 | 0 | 0% |
Flagging the paralog-explained sequences: catches 21 of 21 wrong reads and flags 29 of the 718 right ones (4.0%).
For 11 of the 21 flagged wrong reads the explaining locus lies at the truth.

### Claims
- **H1 HOLDS:** 42.0% wrong against 0% (0 of 676).
- **H2 HOLDS,** against my prediction: 100% of the wrong reads caught, 4.0% of the right ones flagged (bars 50% and 20%).
- **H3 HOLDS narrowly:** 11 of 21 (52%; bar 50%).

### Reading
- **This separates what row 1978 could not.** A recurring difference is harmless when it matches no other copy (an
  allele, a private variant: 676 reads, none wrong) and a warning when another copy carries exactly those bases (50
  reads, 21 wrong). The reads then look like that other copy read at the wrong place, or a copy of the individual that
  is a mixture of two reference copies.
- **For half of the flagged wrong reads the explaining copy is where the read truly comes from,** so the flag also
  names a candidate for the right place.
- **Narrow, and to be said with the number:** 21 wrong reads, 19 of them in one family (NPIP), one individual. The test
  reaches only reads in exactly tied groups with enough depth and a substitution: 739 of the 5,257 decided reads with a
  truth, and 21 of the 140 wrong ones. It is a flag for a minority of the errors, not a repair of the verdict.
- Not confirmed on other data: KB3781 and the planted sets have almost no wrong verdicts to find (0.7%, 0-0.1%), so
  they can test the false-flag rate only. Insertions and deletions are not used.
