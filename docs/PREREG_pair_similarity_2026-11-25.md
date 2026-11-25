# PREREG: are the member places on "another gene" unlisted paralogues? (row 1957, 2026-11-25)

Row 1956: at 0.9 or more of the best score, 18% of the A119b member rows (34% once repeat-dominated rows are dropped) are
on a gene that is no listed relative, no relative at one remove and not of the same named family. The user asked to
resolve this class first, then to add a repeat filter to `members`.

## Method (`bench/o2_truth/pair_similarity.py`)
For every gene in the rows at >= 0.9: its exons merged, spliced together, repeat bases (CHM13 RepeatMasker, the classes
of row 1956) replaced by N. All sequences against all (`minimap2 -c -X -x asm20 -k15 -w10 -N 200 -p 0.05`). A pair is
SIMILAR when the aligned, non-repeat bases cover at least 30% of the shorter sequence's non-repeat bases and at least
200 bases. A row is similar when the source gene is similar to any gene under the place. This uses the annotated gene
models only, not the reads.
Controls on the same sequences: POSITIVE = the rows of class 1 (listed relative) that are not repeat-dominated;
NEGATIVE = 5,000 random pairs of the same genes that are not listed relatives of each other.
A pair lost to the `-N 200` cap counts as not similar (conservative).

## Claims
- **D1, the test works:** at least 85% of the positive-control rows are similar, at most 2% of the random pairs.
- **D2, the unresolved class is mostly paralogues:** at least 70% of the class 5 rows are similar. Then relatives are at
  least 87% of the repeat-filtered rows ((3,901 + 0.7 x 2,078) / 6,127).

## Predictions (mine)
D1 holds. D2 holds at about 80%: the examples are clone-named genes beside NBPF and PRAMEF members.
If D1 fails the test is not usable and D2 is not judged.

## Result (job output `r1957_pair_similarity/run.*.out`; 8,267 genes, 972 with under 200 non-repeat exonic bases)

| rows at 0.9 and above | rows | similar | share |
|---|---|---|---|
| 1 listed relative, not repeat-dominated (positive control) | 2,746 | 2,371 | 86.3% |
| 1 listed relative, repeat-dominated | 1,189 | 1,028 | 86.5% |
| 2 relative at one remove | 693 + 307 | 613 + 282 | 88.5% / 91.9% |
| 3 same named family | 462 + 299 | 375 + 205 | 81.2% / 68.6% |
| 4 repeat | 3,596 | 1,165 | 32.4% |
| 5 another gene: protein-coding | 1,072 | 689 | 64.3% |
| 5 another gene: pseudogene | 919 | 561 | 61.0% |
| 5 another gene: long non-coding, other | 87 | 34 | 39.1% |
| **5 all** | **2,078** | **1,284** | **61.8%** |
| 6 nothing annotated | 148 | 0 | (no gene to compare) |
| random pairs, not listed relatives (negative control) | 5,000 | 3 | 0.06% |

- **D1 HOLDS** (86.3% and 0.06%).
- **D2 FAILS** (61.8%, bar 70%; I predicted about 80%). The test misses 14% of the listed relatives, so the true share
  is somewhat higher (61.8 / 0.863 = 72% if the miss rate is the same; an estimate, not a result).
- **Relatives among the repeat-filtered rows:** (3,901 + 1,284) / 6,127 = 84.6% measured.
- **Not expected:** 32% of the repeat-dominated rows of class 4, and 86% of the repeat-dominated listed relatives, are
  similar in their NON-repeat sequence. A place that is mostly repeat is often still at a real relative (the reads end
  in a repeat-rich untranslated region). Over all 11,518 rows: relatives (classes 1-3 or similar) 8,145 (70.7%) without
  a filter; with the repeat filter 5,185 of 6,127 (84.6%), and 2,960 relative rows (36%) are dropped.

## Step 2: the repeat filter in `members` (opt-in; job 4506082)
`tools/member_discovery.py --repeats BED [--max-repeat 0.5]`, driver `RUSTLE_MEMBERS_REPEATS=BED` (`RUSTLE_MEMBERS_MAX_REPEAT`):
a last column `repeat_share`, and rows at or above the limit left out of both files (2 = column only). Unset, the stage
is unchanged. Run through the driver on A119b with the eight repeat classes of the CHM13 RepeatMasker track
(`r1956_place_classes/rmsk/repeats.bed`): 97,978 rows -> 23,794 (74,184 left out), no row that the earlier table lacks;
by best score ratio 0.3-0.7: 52,115 -> 9,362; 0.7-0.9: 26,500 -> 6,338; 0.9 and above: 19,363 -> 8,094 (the number row
1956 gives for the rows under half repeat). 53 minutes, as before. Whole `cargo test` passes.
Known cost (above): at 0.9 and above the filter also leaves out 36% of the rows at real relatives.
