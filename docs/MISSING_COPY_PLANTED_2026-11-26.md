# A planted sample with copies the reference lacks, and a first detector (rows 1982, 1983)

## Current state (2026-11-26)
- The sample exists and is cheap to make (`hide_copy.py`). Every read of a hidden copy is placed, full length, on relatives.
- The existing sequence test (`recurrence_test.py`, tied reads only) sees a hidden copy when its reads tie (43 of 47) and
  not when they are placed uniquely (1 of 58).
- A detector that reads every read of a locus (`locus_sequences.py`) flags 0.82 of the groups holding a hidden copy on
  fresh sets, and 0.17-0.20 of the groups holding none: it fails the registered false-flag bar (0.10).
- Nothing is wired into the driver; nothing is a default; nothing has run on a real sample.
- Why a two-allele single locus shows a third sequence is known (row 1985): two isoforms folded into one cell, each
  allele counted once per isoform. A count of mutually conflicting sequences removes it on the same seeds (development);
  judged on fresh seeds in row 1986: false flags 0.05-0.06 (bar 0.10 met), groups flagged 0.66 (bar met), hidden
  copies >= 6 edits away found 0.58 (bar 0.70 missed).

Pre-registrations with their outcome tables: `PREREG_missing_copy_plant_2026-11-26.md` (row 1982),
`PREREG_missing_copy_detector_2026-11-26.md` (row 1983).

## Why
Until now every planted set had sample = reference, or a sample whose copies differ from the reference by their own
variants (`plant.py --private`, `--het`; rows 1924, 1927). The case that matters for copies the reference lacks or
collapses (O3) was missing: reads from a copy that has no place of its own on the reference.

## The sample
`bench/o2_families/hide_copy.py PLANT_DIR GENOME.fa OUT.fa [--seed N] [--share 0.25]`

- Input: a planted set (`plant.py`, then `units.py` for `units.tsv` and `dist.tsv`). The reads are not touched: every
  copy is expressed as planned and each read's copy of origin stays in `SC.truth.tsv`.
- Per family of two or more copies, a quarter of the copies (at least one, never all) are drawn with `--seed` among
  copies whose span overlaps no other planted copy. Their span, first to last planted exon, becomes N in `OUT.fa`.
- Coordinates do not move, so `loci.gtf`, `units.tsv` and the truth tables stay valid on the reduced reference.
- `OUT.fa.hidden.tsv`: family, copy, name, unit, span, `nearest` (the remaining copy with the fewest edits), `edits`,
  `twin` (1 = a copy of the same unit, identical in every planted isoform, remains).
- Then map and assemble as for any planted set: `run_tools.sh PLANT_DIR SCENARIO OUT.fa OUTDIR`.

Cost: one 3 GB genome copy per hidden draw; three sets x two scenarios mapped and assembled in 7 minutes on one node.

Classes used everywhere below, by the edits between a hidden copy's first isoform and the nearest remaining copy on the
REFERENCE: `0 twin` (0 edits), `1-5`, `6-19`, `>= 20`, `no alignment`. With `--private 0.001 --het 0.001` every copy
also carries about one variant of its own per kb and a second allele, so a twin is not identical in the sample.

## What a hidden copy looks like (row 1982, `missing_copy_census.py`)
Three sets of row 1927 (seeds 1939-1941), the same 35 of 133 copies hidden in each (one hide seed), scenario E:

| question | result |
|---|---|
| reads unmapped | 1 of 6,300 |
| reads soft-clipped | none, except one copy with no aligned relative (16% of the read) |
| primary on a planted family copy | 0.89; the rest on relatives outside the planted set (non-coding members) |
| share of the placed reads on the top absorbing copy | 0.67 (median, copies >= 20 edits away): the reads scatter |
| the top absorber is the nearest copy | 16-21 of 35 |
| differences per kb against the absorber's own reads | 1.74-1.86 x for copies >= 20 edits away; 0.98-1.20 x for twins |
| seen by the tied-reads sequence test | tied reads (MAPQ 0 for >= 20%): 43 of 47; placed uniquely: 1 of 58 |

On the fresh sets of row 1983 (other hidden copies) the share on a planted family copy was lower (0.75 in one set):
that figure depends on which copies are hidden and on what else of the family the genome holds.

## The detector (row 1983, `locus_sequences.py`)
`locus_sequences.py ASSEMBLY.gtf BAM AS_TABLE GENOME.fa OUT.tsv [--min-seq-reads 3] [--names]`

1. Reads of a locus = its primaries plus the secondaries exactly tied on ms with their read's best alignment.
2. Per locus and isoform (intron chain), the decision of `recurrence_test.py` (IsoCon's test, row 1926) gives the
   distinct sequences; one counts with >= 3 reads.
3. Sequences of different loci or isoforms that share at least half of the smaller one's reads are one sequence
   (a tied read is read at each of its loci). Cells (locus, isoform) holding one sequence form a class.
4. flag = a class has more sequences than 2 x the loci it is on: more than two alleles per reference copy.

Step 4 was changed once on development data. With slots = 2 x the loci of the whole tie group, 10 of 35 hidden copies
were found (groups chain through a few tied reads and the slots grow with them); with slots per class 26-29 of 35.
The rule was frozen after that (`r1983_missing_copy_detector/frozen.md5`) and judged on fresh sets only.

## Judged result (row 1983, `missing_copy_score.py`)
Fresh sets: plant seeds 1983-1985, hide seeds 2983-2985 (35 of 133 copies each, a different draw per set), 60 reads
per copy, scenario E, three sets pooled. A group is positive with >= 10 primaries of hidden copies, null with none and
>= 20 primaries. The same reads on the whole reference are a second null (every flag there is false).

| | bar | planted errors (5 per kb) | near-HiFi errors (1.3 per kb) |
|---|---|---|---|
| groups holding a hidden copy flagged | >= 0.60 | 56 of 68 (0.824) HOLDS | 56 of 71 (0.789) |
| hidden copies >= 6 edits away found | >= 0.70 | 37 of 46 (0.804) HOLDS | 37 of 46 (0.804) |
| false flags, null groups, reduced reference | <= 0.10 | 10 of 60 (0.167) FAILS | 9 of 64 (0.141) |
| false flags, whole reference | <= 0.10 | 35 of 171 (0.205) FAILS | 24 of 180 (0.133) |
| whole reference, single-locus groups | | 11 of 98 (0.112) | 6 of 105 (0.057) |
| whole reference, groups of tied loci | | 24 of 73 (0.33) | 18 of 75 (0.24) |
| hidden copies found: twins / 1-5 / 6-19 / >= 20 | | 26 of 32 / 21 of 24 / 12 of 12 / 25 of 34 | 30 of 32 / 23 of 24 / 12 of 12 / 25 of 34 |
| hidden copies found: tied / placed uniquely | | 38 of 48 / 47 of 56 | 43 of 48 / 47 of 56 |

Scenario D (one copy 10 x, the others a third): groups 50 of 67 and 49 of 68; copies >= 20 edits away 18 of 34;
false flags 8 of 58 and 6 of 69.

What may be said: a flag is four to five times more frequent where a copy is hidden than where none is. It is a
ranking signal, not a certificate.

## False flags (row 1985, description)
Every counted sequence of the 17 falsely flagged single-locus groups (six whole-reference runs, nothing hidden) was read
against the read names and the planted variant tables (`r1985_third_sequence/why.py`, `why.txt`, `by_isoform.txt`).

- 14 of 17: the locus has two planted isoforms, and the intron chain of the second is a contiguous piece of the first.
  `recurrence_test.by_isoform` treats such a chain as shorter reads of the longer isoform (IsoCon's end-invariance,
  which is right for truncated reads) and puts both in one cell. The second isoform covers sites the first does not,
  so its patterns do not fold into the first's, and each allele is counted once per isoform: four sequences from two
  alleles, or three when the second isoform holds no heterozygous site.
- 14 of the 20 single loci with two isoforms in one cell are flagged; 3 of the 183 with one isoform.
- The 14 are three genes recurring over seeds and error arms (CHM13_G0018251, CHM13_G0056284, CHM13_G0018345).
- Not the cause: reads of another copy (no foreign primary in any of the 17); recurring errors (3 of 17 carry a site
  absent from the planted variants, all at 5 errors per kb, in sequences of 3 reads).

A count that removes this, DEVELOPMENT ONLY (written after reading these seeds; `locus_sequences.py --count conflicts`,
opt-in): per cell, the largest set of sequences that pairwise disagree at a site both cover. Two sequences that agree
wherever both are covered can be one haplotype seen through two isoforms or two read lengths.

| same six runs | count of sequences (row 1983) | conflicting sequences |
|---|---|---|
| false flags, null groups, reduced reference | 10 of 60 | 2 of 60 |
| false flags, whole reference | 35 of 171 | 12 of 171 |
| whole reference, single-locus groups | 11 of 98 | 2 of 98 |
| whole reference, groups of tied loci | 24 of 73 | 10 of 73 |
| groups holding a hidden copy flagged | 56 of 68 | 53 of 68 |
| near-HiFi: false flags null / whole reference | 9 of 64 / 24 of 180 | 1 of 64 / 9 of 180 |
| near-HiFi: groups holding a hidden copy flagged | 56 of 71 | 51 of 71 |
| near-HiFi: hidden copies >= 20 edits away found | 25 of 34 | 20 of 34 |

These numbers are not a judged result: the same seeds were used to find the cause. A judged row needs fresh plant and
hide seeds with the bars of row 1983. Still unexplained: the false flags left in groups of tied loci (10 of 73).

## The conflict count judged on fresh seeds (row 1986)
Plant seeds 1986-1988, hide seeds 2986-2988, scenario E pooled, the bars of row 1983
(`PREREG_missing_copy_conflicts_2026-11-27.md`).

| | bar | conflicting sequences | count of sequences, same sets |
|---|---|---|---|
| groups holding a hidden copy flagged | >= 0.60 | 42 of 64 (0.656) HOLDS | 47 of 64 (0.734) |
| false flags, null groups | <= 0.10 | 3 of 66 (0.045) HOLDS | 8 of 66 (0.121) |
| false flags, whole reference | <= 0.10 | 10 of 172 (0.058) HOLDS | 26 of 172 (0.151) |
| hidden copies >= 6 edits away found | >= 0.70 | 26 of 45 (0.578) FAILS | 32 of 45 (0.711) |

Near-HiFi errors: conflicts 38 of 67 groups, false flags 2 of 74 and 4 of 181; sequences 44 of 67, 6 of 74 and 21 of 181.
The trade: the conflict count cuts false flags to about a third and loses 5-6 groups in 64-67. A flag is 11 x more
frequent where a copy is hidden than where none is (count of sequences: about 5 x); 26 x at near-HiFi errors.
Both counts miss the most diverged hidden copies most often (>= 20 edits: 18 and 23 of 34), and the count of sequences
itself scored lower here than on row 1983's seeds (0.734 against 0.824): which copies are hidden matters, so quote
the range over both rows.

## How to reproduce
`r1982_missing_copy_plant/run.sbatch` (sample and census) and `r1983_missing_copy_detector/run.sbatch` (fresh sets,
both error arms, whole-reference control, scoring; 55 minutes, most of it the other assemblers that `run_tools.sh`
also runs). Pooled counts: `r1983_missing_copy_detector/pooled.txt`; per set `sSEED.ARM/{E,D}.score.txt`.

## Limits
- "Found" credits a hidden copy when its group is flagged; one chained group (TSPY, 33 loci) credits every hidden
  copy in it. The group-level row is the fair unit.
- Credit uses read names, a use of the truth in scoring.
- 60 reads per copy; real family loci often have two.
- An N block is not a collapse; for the mapper both mean no sequence there.
- Hidden copies are annotated reference copies: their distance to the family is the reference's, not that of a young
  duplicate private to an individual.
- No planted locus has three real alleles, and no real sample has been run.
- One detector. `copy_phasing` was not run (it needs the inputs of the legacy assignment stage).
