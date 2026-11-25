# PREREG: a lower seeding threshold, and the relatives layer at 0.7-0.9 (rows 1959 and 1960, 2026-11-25)

The user: could the share of the best score be lowered, or should there be two modes, one for ties (seeding) and one for
lower-score relatives? My answer before any run: do not lower the seeding threshold (lower-score places are silent
relatives and repeats, rows 1955-1957), keep lower-score places as a list (`members`). The user asked for both to be measured.

## Row 1959: seeding pools at 0.95 and 0.90 of the best AS (workspace `r1959_seed_ratio_sweep`)
Arms by `bench/o2_families/pool_bam.py --rule as:0.95 | as:0.90` + `--seed-pool all` (the route of rows 1942 / 1951; `as:0.98`
reproduces today's pool exactly, row 1942), against today's assemblies:
- planted: row 1952's seeds 1967-1969, sets E and S (60 reads per isoform) and low E (3 reads), nine sets, with the layer
  split of `support_tier.py`; same environment as row 1952's assemblies;
- real: OR6737 as in row 1951's `real.sbatch` (same settings, `toolcompare.py`, A0 = `r1937_shared_floor1/or`).
Claims:
- **P1:** overall isoform precision at 0.95 and at 0.90 is not above today's in any of the nine planted sets.
- **P2:** exact transcripts rise by at most 2% in every planted set at both thresholds.
- **P3 (real):** annotated family copies on OR6737 rise by at most 1 at both thresholds, and junction precision does not rise.
Predictions (mine): all hold; precision falls a little at 0.95 and more at 0.90; the unspliced-secondary restriction of
the pool limits the effect.

## Row 1960: the member places at 0.7-0.9 with the repeat filter (workspace `r1960_relatives_07`)
Rows 1956 / 1957 repeated unchanged on the band 0.7-0.9 of the A119b table (`place_classes.py`, `pair_similarity.py` with
its band argument). Relatives = classes 1-3, or class 5 and similar; over the rows with an annotated source and a repeat
share under 0.5.
- **Q1:** relatives are at least 80% of those rows (the bar of row 1956). At 0.9 and above the figure was 84.6%.
Prediction (mine): fails, 60-70%.

## Row 1960 result (job 4506251)
Band 0.7-0.9: 15,867 rows with an annotated source, 10,898 of them repeat-dominated (69%); of the 4,969 others,
relatives (classes 1-3, or class 5 and similar) are 3,324 = **66.9%. Q1 FAILS** (bar 80%; predicted 60-70%). The same
computation at 0.9 and above gives 5,185 of 6,127 = 84.6%, as in row 1957. In this band the similarity test finds only
64.6% of the listed relatives (86.3% at 0.9 and above): more diverged pairs, so the 66.9% is an undercount of unknown
size. Random pairs: 1 of 5,000.

## Row 1959 result, planted (job 4506252). Overall = transcripts on planted copies; exact in brackets.
| set | today (0.98) | 0.95 | 0.90 | tied_only transcripts, today / 0.95 / 0.90 (exact) |
|---|---|---|---|---|
| 1967 E | 0.654 (197) | 0.595 (197) | 0.547 (197) | 58 / 84 / 113 (0 / 0 / 0) |
| 1968 E | 0.654 (197) | 0.610 (197) | 0.553 (197) | 62 / 83 / 115 (0 / 0 / 0) |
| 1969 E | 0.659 (197) | 0.606 (197) | 0.540 (197) | 54 / 78 / 117 (0 / 0 / 0) |
| 1967 S | 0.452 (103) | 0.410 (103) | 0.377 (103) | 59 / 85 / 105 (0 / 0 / 0) |
| 1968 S | 0.417 (105) | 0.383 (105) | 0.342 (105) | 83 / 105 / 138 (0 / 0 / 0) |
| 1969 S | 0.488 (104) | 0.441 (104) | 0.419 (104) | 48 / 66 / 77 (0 / 0 / 0) |
| 1967 low E | 0.755 (166) | 0.712 (168) | 0.697 (166) | 50 / 67 / 74 (24 / 25 / 26) |
| 1968 low E | 0.726 (162) | 0.752 (161) | 0.750 (159) | 52 / 49 / 52 (22 / 21 / 21) |
| 1969 low E | 0.740 (168) | 0.681 (173) | 0.671 (173) | 47 / 73 / 78 (22 / 23 / 23) |
False copies in the S sets: 45 / 47 / 50, 48 / 52 / 54, 40 / 45 / 47.
- **P1 FAILS** by one set: 1968 low E is more precise at both thresholds (0.752, 0.750 against 0.726); the other eight
  are less precise at both. **P2 FAILS** by one set: 1969 low E gains 5 exact transcripts (+3.0%); elsewhere -3 to +2.
- Reading: at 60 reads per isoform a wider pool adds only wrong transcripts (exact unchanged in all six sets, precision
  -4 to -6 points at 0.95 and -7 to -12 at 0.90, more silent copies modelled). At 3 reads it is noise around zero:
  exact +2, -1, +5 at 0.95 and 0, -3, +5 at 0.90, precision -0.04, +0.03, -0.06 and -0.06, +0.02, -0.07. The primary
  layer does not move in any set. No case for a lower threshold; the real-sample arm (P3) follows.

## Row 1959 result, real (OR6737, job 4506253; `toolcompare.py`, 49 annotated family copies, 47 modelled by any arm)
| | today (0.98) | 0.95 | 0.90 |
|---|---|---|---|
| annotated family copies with a spliced transcript | 37 | 44 | 45 |
| junction precision at them | 0.67 | 0.65 | 0.66 |
| spliced transcripts at them | 176 | 314 | 435 |
| transcripts in all, with the mode (rescued) | 73,173 (2,992) | 75,144 (3,673) | 77,212 (3,941) |
- **P3 FAILS, against my prediction:** 7 and 8 more annotated family copies get a transcript; junction precision does
  not rise (0.67 -> 0.65, 0.66). By the yardstick rows 1942 / 1951 used against the narrower pools (copies lost = a
  cost), this is a gain, at 1.8 and 2.5 times the transcripts on those copies.
- What it does not show: that the added copies are EXPRESSED. On planted reads the same widening models more silent
  copies and adds no true transcript at 60 reads. OR6737 has no expression truth, so the two readings are not separated
  here. The other three real samples follow (KB3781, A119b, HG002).

## Row 1959, the other real samples and what the added copies are (jobs 4506299-4506301; `bench/o2_truth/gained_copies.py`)
| annotated family copies with a spliced transcript / junction precision / spliced transcripts | today (0.98) | 0.95 | 0.90 |
|---|---|---|---|
| OR6737 (of 49) | 37 / 0.67 / 176 | 44 / 0.65 / 314 | 45 / 0.66 / 435 |
| KB3781 (of 49) | 28 / 0.63 / 108 | 34 / 0.60 / 198 | 45 / 0.67 / 276 |
| HG002 (of 63) | 37 / 0.45 / 442 | 39 / 0.44 / 520 | 44 / 0.41 / 615 |
P3 fails on all three (A119b not finished when this was written).

Description, after the result: for each ADDED copy, the primary alignments with mapping quality >= 1 on >= 50 of its exon
bases (reads placed there without a tie-break):

| | copies | 0 such reads | 1-4 | 5 or more |
|---|---|---|---|---|
| OR6737, modelled today | 37 | 1 | 8 | 28 |
| OR6737, added at 0.95 | 8 | 5 | 2 | 1 (35 reads) |
| OR6737, added at 0.90 | 10 | 5 | 2 | 3 (21, 23, 35) |
| KB3781, modelled today | 28 | 0 | 3 | 25 |
| KB3781, added at 0.95 | 7 | 3 | 0 | 4 (21, 26, 45, 76) |
| KB3781, added at 0.90 | 18 | 14 | 0 | 4 (the same four) |
| HG002, modelled today | 37 | 10 | 7 | 20 |
| HG002, added at 0.95 | 2 | 2 | 0 | 0 |
| HG002, added at 0.90 | 8 | 6 | 2 | 0 |
HG002 origin truth (row 1785's sample of family reads): 27 of today's 37 copies have a read whose origin is the copy; of
the copies added, 0 of 2 and 0 of 8 (5 of the 8 appear only as another candidate of reads from elsewhere). The KB3781
truth table has no read on these family copies (0 of 28 for today's too): not usable here.
- Reading: the added copies are of two kinds. A FEW have many reads of their own that today's pool fails to turn into a
  spliced transcript (4 on KB3781, 1-3 on OR6737): a real gain, and all of it is reached at 0.95. MOST have no read of
  their own (5 of 8, 3 of 7, 2 of 2 at 0.95; 5 of 10, 14 of 18, 6 of 8 at 0.90): transcripts built from relatives' reads.
  Going from 0.95 to 0.90 adds almost only the second kind (KB3781: +11 copies, all with 0 reads).
- Not established: why the copies with 21-76 placed reads get no spliced transcript today (worth a look on its own).

## Why copies with 21-76 reads of their own are "not modelled" today (user: explore why the 76 reads disagree; KB3781)
They are modelled. The four copies sit in tandem units of two adjacent annotated genes (a JMJD7-like gene of 9-10
exons and a cytosolic phospholipase A2 beta-like gene of 20 exons: the JMJD7-PLA2G4B read-through arrangement), and
their reads are READ-THROUGH transcripts over both genes:
- LOC129531983 (JMJD7-like, 9 exons): all 76 placed reads run beyond the gene into the neighbour LOC129531982; they have
  15-26 introns (median 22). Today's assembly has 7 transcripts of 20-26 exons over the pair.
- LOC134758743 (PLA2G4B-like, 20 exons): 40 of its 45 placed reads run beyond it into the JMJD7 gene next to it; today's
  assembly has two transcripts of 25 and 26 exons over the pair. LOC129534646 and LOC134758753: the same picture (three
  and two transcripts of 24-26 exons over the pair).
The benchmark (`toolcompare.py`, and `gained_copies.py` after it) gives a transcript to ONE gene, the one it shares the
most exonic bases with, so the other gene of the pair counts as not modelled although its reads are all in a transcript.
What makes the copy "modelled" at 0.95 are new SHORT transcripts inside one gene (10 and 19 exons at LOC134758743),
compatible with 1 and 2 of the placed reads and with 15 and 61 secondary alignments: built from relatives' reads.
**So the "real gain" at 0.95 is not one.** These copies are an artifact of the count, and the transcripts that satisfy
the count rest on secondary alignments. With that, no added copy in the three samples is shown to be a copy with reads
of its own that today's threshold misses (OR6737's one to three such copies are in the same family and were not opened
one by one).
Consequences: (1) rows 1942 and 1951 judged narrower pools by this same count ("copies lost"); at read-through pairs a
change in which short transcripts exist moves the count without any change in what is expressed, so those losses were
not checked for this either. (2) A count that credits every gene a transcript covers would be the fairer yardstick
(not built).

## Summary of rows 1959 and 1960
- Planted reads: a wider seeding pool adds only wrong transcripts at 60 reads per isoform and noise at 3.
- Real samples: more annotated copies get a transcript (37 / 44 / 45, 28 / 34 / 45, 37 / 39 / 44), but the added copies
  either have no read of their own, or are read-through partners already covered today. HG002 truth: 0 of the added
  copies have a read that originates there.
- Relatives list: 84.6% real relatives at 0.9 and above with the repeat filter, 66.9% at 0.7-0.9.
- Recommendation: seeding stays at 0.98; lower-score places belong to the `members` list at 0.9 with the repeat filter.
  A119b's arm was still running when this was written (job 4506300).
