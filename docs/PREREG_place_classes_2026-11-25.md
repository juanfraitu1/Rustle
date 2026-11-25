# PREREG: what are the real-sample member places that are not on a listed relative? (row 1956, 2026-11-25)

Row 1955, addendum 2: on A119b only 35% of the member places at 0.9 or more of the read's best score are on a listed
relative of the source gene; 62% are on "another annotated gene". The user asked to classify those, starting at 0.9.
If most are relatives the table lacks, a second family layer at 0.9 is defensible; if not, it is not.

## Data
`r1915_members_blocks/a1.members.tsv` (A119b, stage `members`: 5 reads, 0.3 of the best AS), rows with best ratio >= 0.9.
Source locus -> annotated gene as in row 1955 (assembly `r1920_mono_from_primaries/real_a1/g.gtf`). Relations:
`annot_related/chm13.related.relations.tsv`. Gene groups: `~/data/denovo_eval/hgnc_complete_set.txt`. Repeats: the public
CHM13 RepeatMasker track (workspace `r1956_place_classes/rmsk/rmsk.bed`).

## Classes of a (place, source) row, first match wins (`bench/o2_truth/place_classes.py`)
1. **listed relative**: a gene under the place is a listed relative of the source gene.
2. **relative at one remove**: a gene under the place and the source share a listed relative (two pseudogenes of one
   parent; a pseudogene and its parent's paralogue).
3. **same named family**: they share an HGNC gene group, or the same symbol once a pseudogene suffix `P<number>` is removed.
4. **repeat**: half or more of the place's aligned footprint is an interspersed or simple repeat (SINE, LINE, LTR, DNA,
   Retroposon, Simple_repeat, Low_complexity, Satellite).
5. **another gene**, by the biotype of the gene under the place (protein-coding, long non-coding, pseudogene, other).
6. **nothing annotated** (split: repeat or not).
Rows whose source is on no annotated gene are counted apart (repeat or not). Counted as rows and as distinct places.
Repeat share is also reported inside classes 1-3.

## Claim
- **C1:** classes 1-3 hold at least 80% of the rows with an annotated source at >= 0.9. Then the layer at 0.9 can be
  called "relatives".
## Predictions (mine)
C1 fails narrowly: classes 1-3 about 60-75%, class 2 the largest of the new ones (ribosomal-protein pseudogenes);
repeats under 10% at this threshold; the rest long non-coding and overlapping genes.
The same table at 0.7-0.9 is reported as a description.

## Result (workspace `r1956_place_classes/classes.txt`)
Rows = (place, source). Shares are of the rows whose source is on an annotated gene.

| class | at 0.9 and above (11,518 rows) | mostly repeat | at 0.7-0.9 (15,867 rows) |
|---|---|---|---|
| 1 listed relative | 3,935 (34.2%) | 1,189 | 2,424 (15.3%) |
| 2 relative at one remove | 1,000 (8.7%) | 307 | 837 (5.3%) |
| 3 same named family | 761 (6.6%) | 299 | 420 (2.6%) |
| **classes 1-3** | **5,696 (49.5%)** | 1,795 | **3,681 (23.2%)** |
| 4 repeat | 3,596 (31.2%) | all | 9,594 (60.5%) |
| 5 another gene: protein-coding | 1,072 (9.3%) | 0 | 748 (4.7%) |
| 5 another gene: pseudogene | 919 (8.0%) | 0 | 1,319 (8.3%) |
| 5 another gene: long non-coding, other | 87 (0.8%) | 0 | 125 (0.8%) |
| 6 nothing annotated | 148 (1.3%) | 0 | 400 (2.5%) |
Rows whose source is on no annotated gene: 7,845 at 0.9 and above, 5,878 of them mostly repeat; 10,633 at 0.7-0.9,
9,266 mostly repeat.

- **C1 FAILS** (49.5%, bar 80%). Both predictions were wrong: classes 1-3 are below the 60-75% I expected, and repeats
  are 31%, not under 10%. A read that ends in or runs through a transposable element aligns almost as well at other
  copies of that element; a high score ratio does not exclude that.
- **Without the repeat-dominated rows** (description, not a registered claim): 6,127 rows remain at 0.9 and above, of
  which classes 1-3 are 3,901 (63.7%), another gene 2,078 (33.9%), nothing annotated 148 (2.4%).
- **Class 5 is not resolved.** Its examples are mostly genes named by clone (AL137798.1, AL669831.1, AC243961.1) beside
  NBPF, PRAMEF and similar families: probably paralogues and pseudogenes without an entry in the relation table, but
  that was not tested (a sequence comparison of the two genes would).
- Reading: at 0.9 the layer cannot be called "relatives" as it stands. A repeat filter is the one clear lever (it
  removes a third of the rows and four fifths of the rows from unannotated sources).
