# PREREG: can secondaries at 0.9 or more of the best score help seed real loci? AS and ms (row 1961, 2026-11-25)

The user's goal: see how the multimappers at more than 0.9 of the best alignment could help seed real loci, with both
scores. Row 1959 showed that taking ALL of them (AS) adds transcripts without reads of their own. The constructive
question: is there a part of that band that supports loci which are really expressed, and can it be picked out?

## Idea: an anchored pool
A secondary below today's pool (AS under 0.98 of the read's best) is admitted only where the sample has reads of its
own: at least 50 of its aligned bases lie where 3 or more primary alignments with mapping quality 1 or more are stacked
(`bench/o2_families/pool_bam.py --rule anchor:as:0.9:3 | anchor:ms:0.9:3`; today's pool is always kept). A place that only
tie-breaks or nothing puts reads at gets no seed. K = 3 and R = 0.9 are fixed here, not tuned.

## Step 0, the mechanism (`bench/o2_truth/band_census.py`; planted sets, truth known)
Every secondary alignment that is not at its read's true copy, by band on each score (>= 0.98, 0.9-0.98, under 0.9 of the
read's best AS; the same on ms): where it lands (an EXPRESSED copy, a SILENT copy, no planted copy), whether the place is
anchored, and whether its intron chain is part of an isoform that is truly expressed at the landing copy (then it adds
support to a real transcript) or not (then it can only build a wrong one).

## Step 1, assembly (the route of rows 1942 / 1951 / 1959: filtered BAM + `--seed-pool all`)
Arms: anchored on AS, anchored on ms, and `ms:0.90` unanchored (row 1959 has `as:0.90`). Nine planted sets (row 1952's
seeds 1967-1969: E, S at 60 reads per isoform, low E at 3) with the layer split; real: OR6737, KB3781, HG002 with
`toolcompare.py` and `gained_copies.py`, to which a second count is added first: a copy is COVERED when a spliced
transcript on its strand overlaps at least half of its exon bases (credits both genes of a read-through pair, row 1959).

## Claims (both anchored arms unless said)
- **C1, safe:** in the three half-silent sets the silent copies with a transcript are at most today's + 1.
- **C2, precise:** overall isoform precision is at least today's - 0.02 in all nine sets.
- **C3, gain on planted reads:** exact transcripts are at least 3% above today's in two of the three low-depth sets, for
  at least one score.
- **C4, gain on real reads:** covered annotated copies are at least today's + 2 on two of three samples for at least
  one score.
## Predictions (mine)
Step 0: most of the band lands on expressed copies in the E sets (all copies are expressed there) but few of those
alignments fit a true isoform of the landing copy; on ms the band holds more spliced alignments than on AS.
C1 holds. C2 fails in the all-expressed sets (every copy is anchored there, so the anchor filters nothing). C3 fails
(+0-2%). C4 fails. So I expect "safer than the open band, no real gain"; the user expects a gain can be found.

## Step 2, registered after step 0 and the planted part of step 1 were read, before any run of it
Step 0 showed that the anchor selects well (in the half-silent sets, 84% of the anchored alignments of the 0.9-0.98 band
fit an isoform truly expressed at the landing copy and 1% land on a silent copy, against 34% and 34% unanchored), and
step 1 that the anchored pool still adds wrong transcripts and no true one. The remaining wrong ones should come from the
anchored alignments that do NOT fit (an isoform of one copy carried to a copy that does not make it). So a stricter rule:
`janchor:S:R:K` = anchored, and every intron of the secondary is an intron of a primary alignment with mapping quality
>= 1 at that place. It cannot create a structure; it can only add reads to one the place's own reads show.
- **J1:** overall precision at least today's - 0.02 in all nine planted sets, both scores (C2, again).
- **J2:** silent copies with a transcript at most today's + 1 (C1, again).
- **J3:** exact transcripts at least 3% above today's in two of three low-depth sets for at least one score (C3, again).
Prediction (mine): J1 and J2 hold, J3 fails (0 to +2%): with the structure already shown by the place's reads the
assembler builds it without help.

## Result, step 0 (job 4507842; seed 1967, the other two seeds alike)
Secondary alignments away from the read's true copy, band 0.9-0.98 of the read's best:

| set, score | alignments | on an expressed copy | on a silent copy | no planted copy | fit a true isoform there | anchored | anchored and fit | anchored on a silent copy |
|---|---|---|---|---|---|---|---|---|
| all expressed, AS | 51,879 | 37,654 | 0 | 14,225 | 33,116 | 18,796 | 15,068 | 0 |
| all expressed, ms | 36,782 | 27,644 | 0 | 9,138 | 22,856 | 16,295 | 11,955 | 0 |
| half silent, AS | 28,159 | 10,691 | 9,681 | 7,787 | 9,555 | 5,827 | 4,897 | 61 |
| half silent, ms | 20,175 | 7,807 | 7,170 | 5,198 | 6,662 | 4,981 | 3,962 | 62 |
| 3 reads, AS | 2,496 | 1,839 | 0 | 657 | 1,562 | 731 | 611 | 0 |
| 3 reads, ms | 1,767 | 1,331 | 0 | 436 | 1,074 | 601 | 469 | 0 |
- **My prediction was wrong:** most alignments of the band that land on an expressed copy DO fit an isoform truly
  expressed there (88-89% in the half-silent set), because copies of a family share isoform structures.
- **The anchor selects well.** Half-silent set, AS: unanchored, 34% fit a true isoform and 34% land on a silent copy;
  anchored, 84% fit and 1% land on a silent copy. It keeps about half of the fitting alignments (4,897 of 9,555).
- **AS and ms differ in size, not in kind:** the band on ms holds about 30% fewer alignments (many of the AS band are
  at 0.98 or more on ms); the shares are the same within two points.
- Below 0.9 almost nothing fits (190 of 12,032; 2,560 of 27,295): 0.9 is a natural lower bound.

## Result, step 1, planted. Overall precision (exact transcripts); silent copies with a transcript in the S sets.
| set | today | open AS 0.90 | open ms 0.90 | anchored AS | anchored ms |
|---|---|---|---|---|---|
| 1967 E | 0.654 (197) | 0.547 (197) | 0.517 (196) | 0.561 (197) | 0.530 (196) |
| 1968 E | 0.654 (197) | 0.553 (197) | 0.528 (196) | 0.564 (197) | 0.540 (196) |
| 1969 E | 0.659 (197) | 0.540 (197) | 0.517 (196) | 0.553 (197) | 0.531 (196) |
| 1967 S | 0.452 (103), 45 | 0.377, 50 | 0.365, 51 | 0.426 (103), 45 | 0.424 (103), 45 |
| 1968 S | 0.417 (105), 48 | 0.342, 54 | 0.328, 54 | 0.382 (105), 49 | 0.367 (104), 49 |
| 1969 S | 0.488 (104), 40 | 0.419, 47 | 0.402, 51 | 0.460 (104), 40 | 0.454 (104), 40 |
| 1967 low E | 0.755 (166) | 0.697 (166) | 0.689 (162) | 0.719 (166) | 0.717 (162) |
| 1968 low E | 0.726 (162) | 0.750 (159) | 0.750 (156) | 0.711 (162) | 0.713 (159) |
| 1969 low E | 0.740 (168) | 0.671 (173) | 0.660 (171) | 0.695 (171) | 0.685 (170) |
- **C1 HOLDS** (silent copies 45, 49, 40 against 45, 48, 40; the open bands give 47-54). **C2 FAILS** (-0.09 to -0.13 in
  the all-expressed sets, -0.03 in the half-silent, -0.015 to -0.055 at low depth). **C3 FAILS** (exact 0, 0, +3 on AS;
  -4, -3, +2 on ms).
- The anchor does what it was built for (no extra silent copy), and the pool still adds wrong transcripts at real,
  expressed copies and no true one. They are in the `tied_only` layer (all-expressed sets: 105-108 transcripts on AS
  against 54-62 today, none exact): structures that no read placed at the copy shows, carried over from a relative.
- ms is slightly worse than AS in every set (more spliced alignments admitted, more carried-over structures).

## Result, step 2 (job 4507868): anchored and intron-supported. Overall precision (exact); silent copies in the S sets.
| set | today | anchored AS | anchored + introns, AS | anchored + introns, ms |
|---|---|---|---|---|
| 1967 E | 0.654 (197) | 0.561 (197) | 0.642 (197) | 0.643 (196) |
| 1968 E | 0.654 (197) | 0.564 (197) | 0.644 (197) | 0.645 (196) |
| 1969 E | 0.659 (197) | 0.553 (197) | 0.648 (197) | 0.649 (196) |
| 1967 S | 0.452 (103), 45 | 0.426, 45 | 0.452 (103), 45 | 0.452 (103), 45 |
| 1968 S | 0.417 (105), 48 | 0.382, 49 | 0.418 (105), 48 | 0.416 (104), 48 |
| 1969 S | 0.488 (104), 40 | 0.460, 40 | 0.475 (104), 40 | 0.475 (104), 40 |
| 1967 low E | 0.755 (166) | 0.719 (166) | 0.746 (167) | 0.747 (165) |
| 1968 low E | 0.726 (162) | 0.711 (162) | 0.726 (162) | 0.732 (161) |
| 1969 low E | 0.740 (168) | 0.695 (171) | 0.739 (170) | 0.741 (169) |
- **J1 HOLDS** (largest loss 0.013). **J2 HOLDS** (45, 48, 40 = today's). **J3 FAILS** (exact +1, 0, +2 on AS; -1, -1,
  +1 on ms).
- The rule keeps 54% of all secondary alignments against
  the anchored rule's 56%, and its `tied_only` layer is back at today's size (60, 61, 55 against 58, 62, 54).
- Reading: secondaries at 0.9 or more of the best score CAN be given to real loci without harm, when they are held to
  places with reads of their own and to introns those reads show. On planted reads they then change nothing that the
  scores see: a structure that the place's own reads show is assembled without them. Their possible use is as added
  support (read counts) for transcripts that already exist, not as seeds of new ones.

## Result, real samples (jobs 4507843 / 4507844 / 4507869 / 4507870; HG002 rerun, see below)
Covered = a spliced transcript on the copy's strand over at least half of its exon bases (`gained_copies.py`); in
brackets, the added copies by the reads placed there without a tie-break (none / 1-4 / 5 or more).

| | today | open AS 0.90 | open ms 0.90 | anchored AS | anchored ms | anchored + introns, AS | anchored + introns, ms |
|---|---|---|---|---|---|---|---|
| OR6737, covered of 49 | 26 | 40 (5 / 5 / 4) | 41 (6 / 5 / 4) | 31 (0 / 1 / 4) | 31 (0 / 1 / 4) | 28 (0 / 0 / 2) | 28 (0 / 0 / 2) |
| OR6737, junction precision / spliced transcripts | 0.67 / 176 | 0.66 / 435 | 0.66 / 449 | 0.61 / 363 | 0.61 / 383 | 0.65 / 269 | 0.65 / 269 |
| KB3781, covered of 49 | 23 | 29 (6 / 0 / 0) | 27 (4 / 0 / 0) | 24 (1 / 0 / 0) | 24 (1 / 0 / 0) | 23 | 23 |
| KB3781, junction precision / spliced transcripts | 0.63 / 108 | 0.67 / 276 | 0.64 / 275 | 0.60 / 235 | 0.59 / 245 | 0.63 / 183 | 0.63 / 186 |
C4 (covered + 2 on two of three samples) is met on OR6737 by both rules and not on KB3781; HG002 decides.

The five copies the anchored pool newly covers on OR6737:

| copy | own reads | today's best transcript covers | new transcript covers | its introns shown by own reads | strict rule too |
|---|---|---|---|---|---|
| LOC129523574 (LRPAP1-like, 7 exons) | 45 | 30% | 100% | 6 of 6 | yes |
| PLA2G4B (21 exons) | 23 | 0% | 100% | 20 of 25 | yes |
| LOC129527105 (6 exons) | 21 | 0% | 50% | 20 of 25 | no |
| SPTBN5 (67 exons) | 7 | 29% | 58% | 18 of 41 | no |
| LOC129534707 (57 exons) | 4 | 33% | 72% | 11 of 41 | no |
- **LOC129523574 is a rescue of the ANNOTATED isoform.** The copy's reads show two 3' ends: 40 splice out at
  11,294,250 (a transcript today's assembly has), 5 run on unspliced through the annotated 2.3 kb last exon. Today's
  assembly builds only the first. With the strict pool, 22 secondary alignments of relatives' reads that run through the
  same exon join the 5, and the annotated 7-exon transcript is assembled. Every intron of it is in the copy's own reads.
  This is the use the goal asks about: reads of relatives lifting a true, thinly supported isoform of a locus that has
  its own evidence.
- SPTBN5 and LOC129534707 (57-67 exons, 4-7 own reads that each cover under 30% of the gene) are assembled from introns
  that no read of the copy shows: possibly right, not shown by the sample; the strict rule does not build them.
- HG002: the first submission ran three arms at once with one pooling worker per contig, each holding the 9-10 GB score
  table, and stalled for 2.5 hours; cancelled by id and rerun with 10 workers, one arm at a time (jobs 4509313 / 4509314).

### HG002, the strict rule (jobs 4509313 / 4509314, 91 minutes an arm) and the open pools; anchored arms follow
| | today | open AS 0.90 | open ms 0.90 | anchored + introns, AS | anchored + introns, ms |
|---|---|---|---|---|---|
| copies with a spliced transcript (of 63) / junction precision / spliced transcripts | 37 / 0.45 / 442 | 44 / 0.41 / 615 | 45 / 0.42 / 656 | 37 / 0.45 / 430 | 37 / 0.45 / 451 |
| covered (of 133 members) | 22 | 27 | 29 | 22 | 23 |
| added to the covered: no own read / 1-4 / 5 or more | | 5 / 0 / 0 | 6 / 0 / 1 | none | 0 / 0 / 1 |
| added to the covered: with a truth-origin read | | 1 of 5 | 2 of 7 | | 1 of 1 |
The strict rule changes nothing on AS and adds one covered copy on ms, which has five or more reads of its own and a
read whose true origin is that copy; junction precision is unchanged. The open pools add 5 and 7 covered copies, of which
1 and 2 have a read that originates there.

### C4 and the reading over three real samples (strict rule)
Covered copies: OR6737 +2 (both confirmed by own reads; one is the annotated isoform of an LRPAP1-like copy), KB3781 0,
HG002 0 (AS) / +1 (ms, confirmed by the origin truth). **C4 FAILS as registered** (+2 on one sample only). No copy is
lost and junction precision moves by 0.02 at most. So: a small, real, safe gain on two of three samples (3 copies in
all on ms, 2 on AS), none of it unsupported; ms is not worse than AS here and finds the one HG002 copy.

### HG002, the anchored arms (the same two jobs)
Anchored AS: 36 copies with a spliced transcript (today 37), junction precision 0.41 (today 0.45), 485 spliced
transcripts, covered 22 (no change). Anchored ms: 37, 0.41, 514, covered 23 (the same confirmed copy the strict ms rule
adds). So on HG002 the anchored rule costs four points of junction precision for nothing the strict rule does not give.

## Summary of row 1961
- In the band 0.9-0.98 of the read's best score a third of the secondaries fit a truly expressed isoform at the place
  they land on and a third land on a silent copy; with an anchor (reads of the place's own) 84% and 1%. Under 0.9
  almost nothing fits. AS and ms give the same shares; the ms band is 30% smaller.
- Open band: wrong transcripts and copies without reads. Anchored: no extra silent copy, still structures the place's
  reads do not show. Anchored and intron-supported: as precise as today on planted reads, no copy lost on real ones,
  and 0-2 confirmed covered copies a sample (OR6737 +2, KB3781 0, HG002 0 on AS and +1 on ms).
- Claims: C1 holds; C2 fails for the anchored rule and holds for the strict one (J1); C3 / J3 and C4 fail (the gain is
  too small). The use of these secondaries is support for an isoform the locus's own reads already show.
