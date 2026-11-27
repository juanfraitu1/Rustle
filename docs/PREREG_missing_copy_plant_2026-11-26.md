# STEP 0 + PREREG: a planted sample with copies the reference lacks (register 1982)

2026-11-26. User: "lets focus on the planted missing copy sample" (open item since row 1927: `plant.py --private / --het`
plant a sample that differs from the reference, never one with a copy the reference lacks). Written before any read
was mapped to a reduced reference.

## The sample
`bench/o2_families/hide_copy.py`: the planted reads are unchanged; a quarter of the copies of every family (at least one,
never all) are replaced by N in a copy of the genome, and the reads are mapped to that genome. Coordinates do not move.
Planted sets: row 1927's three (`r1927_het_alleles/h1939|h1940|h1941/plant`, `--private 0.001 --het 0.001`, 60 reads
per copy), scenarios E (all copies expressed) and D (one copy 10 x, others a third). These seeds were looked at in row
1927 for another question; they are development data here. The same reads mapped to the whole genome exist in row 1927
(the control of every locus).

## What is measured (`bench/o2_families/missing_copy_census.py`)
Per hidden copy: where the primaries of its reads land (remaining copy of the family / elsewhere / unmapped), the share
on the top absorbing copy, differences per kb of its reads there against the absorber's own reads, soft-clipped share,
MAPQ 0 share; all split by the edits between the hidden copy and the nearest remaining copy (0 = an identical twin
remains, 1-5, 6-19, >= 20). Then the existing sequence test (`recurrence_test.py`, row 1926) on the reduced-reference run:
a hidden copy is SEEN when a called sequence has >= 5 reads and >= 80% of them are its own (credit by read names: a
use of the truth in scoring).

## Predictions
- P1: >= 90% of the hidden copies' reads have a primary on a remaining copy of their family.
- P2: for hidden copies >= 20 edits from the nearest remaining copy the top absorber takes >= 80% of the placed reads.
- P3: same copies: their reads carry >= 1.5 x the differences per kb of the absorber's own reads (median of ratios);
  for twins (0 edits) the ratio is under 1.2.
- P4: of hidden copies >= 6 edits from the nearest with >= 20 placed reads, the sequence test sees >= 50%;
  it sees none of them when the absorbing locus has no tied reads (the test reads tied groups only).

## Reading
This row makes the sample and says what a hidden copy looks like in the alignments. No detector is judged here; the
detector row is registered after, on fresh seeds, with a null (loci with two alleles and no hidden copy).

## Carried caveats ([[planted-benchmark-rules]])
Planted errors are 5 per kb (real HiFi 1.1-1.4); 60 reads per copy; an N block instead of a collapse; hidden copies
are annotated reference copies, so their divergence from the family is the reference's, not that of a young
individual-specific duplicate.

## Outcome (2026-11-26, job 4512616, `r1982_missing_copy_plant/h19*/{E,D}.census.{txt,tsv}`)
The same 35 of 133 copies are hidden in the three sets (one hide seed): the sets differ in reads and variants only.
Scenario E, three sets:
- P1 FAILS narrowly: 0.890 / 0.893 / 0.896 of the hidden copies' reads on a planted family copy; the others on relatives
  outside the planted set (non-coding members, e.g. chr15:82.2 Mb for GOLGA6L4, chrY 9.0-9.7 Mb for TSPY). Unmapped: 1 read.
- P2 FAILS: top absorber share for copies >= 20 edits away = 0.667 (median) in every set.
- P3: holds for >= 20 edits (1.86 / 1.74 / 1.85); twins 0.98 / 1.00 / 1.20 (not under 1.2 in one set).
- P4 FAILS: 6 / 6 / 7 of 16 hidden copies >= 6 edits away are seen by the sequence test. By placement, all classes:
  tied (MAPQ 0 for >= 20% of the reads) 43 seen of 47; placed uniquely 1 of 58.
What a hidden copy looks like: every read placed, full length, on relatives; when a relative is near-identical the reads tie
and the existing sequence test isolates them; when the copy is diverged its reads sit uniquely on one to three relatives
with about twice the differences per kb, and nothing in the pipeline reads them.
Not done: a detector with a null; other hide seeds; real error rate (planted 5 per kb).
