# PREREG: seed with uniquely placed primaries first, then bring in the tied reads (row 1963, 2026-11-25)

The user's idea: today's seeding takes every primary alignment, also the ones minimap2 placed by a tie-break (mapping
quality 0). Seed first with the primaries that have a real placement, then bring in the near-tied multimappers together
with their arbitrary primaries.

Known before this row: assembling without the MAPQ-0 primaries removes the made-up loci and half of the isoforms
(workspace `coin_toss_2026-11-25`); added as plain seeds, step 2 gives back today's assembly. So the idea changes
something in two forms only, tested as two arms.

## Definitions (fixed here)
- **Tied read:** a rival alignment within 10 `ms` of its primary (the band of the labels; mapping quality 0 would miss
  the near ties that are not exact). **Uniquely placed:** every other read.
- Tool `bench/o2_families/two_step.py`; workspace `r1963_two_step_seeding`.

## Arm A, the layers (no new assembly)
Today's layer `primary` = at least 2 primary alignments carry the chain, tie-broken ones included. New layer `unique` =
at least 2 UNIQUELY PLACED primaries carry it (relation `exact`, `exact_after_replacement` or `unspliced` in the
provenance table); everything else of today's assembly is layer `tied`.

## Arm B, a conditioned step 2 (filtered BAM, then the usual stages with `--seed-pool all`)
Step 1 = an assembly of the uniquely placed primaries alone (`--seed-pool primary`). Step 2 starts from today's seeds
(primaries and secondaries at AS >= 0.98 of the best) and only removes: for a tied read, among its seeds within 10 `ms`
of its best, when at least one lies on a step-1 transcript (50 exon bases), the ones that do not are dropped, the
primary included (a kept one is then made the primary record). When none lies on a step-1 transcript the group stays
whole, as today. A tied primary and a tied secondary are treated alike.

## Data
Nine planted sets of row 1952 (seeds 1967-1969: all expressed E and half silent S at 60 reads per isoform, low E at 3;
seen by rows 1952-1961, said here) scored with `score.py` and `units_score.py`; real OR6737 for arm B with
`toolcompare.py --species gorilla` against today's `r1937_shared_floor1/or/g.multicopy.gtf`.

## Claims
- **A1:** layer `unique` has at most half of layer `primary`'s false units in all three S sets.
- **A2:** its isoform precision is at least `primary`'s + 0.05 in all three S sets and not more than 0.02 below it in
  the other six. (Its exact isoforms are reported as the cost, no bar.)
- **B1, nothing lost:** exact isoforms at least today's - 2 in all nine sets.
- **B2, fewer made-up copies:** false units at most half of today's in all three S sets.
- **B3, precision:** isoform precision at least today's + 0.03 in six of the nine sets.
- **B4, real:** OR6737 annotated family copies at least today's - 1 and junction precision not below today's.

## Predictions (mine)
A1 and A2 hold; `unique` keeps 50-60% of the `primary` layer's exact isoforms at 60 reads and under half at 3. B1 holds.
B2 holds. B3 fails: the gain is in the S sets only (3 of 9), since in the E sets every copy has reads of its own. B4
holds, with fewer rescued loci. Against today's LABELS (not today's unlabelled assembly) arm B gains little: the labels
already mark a silent twin.

## Reading
A holds: the `unique` layer is the precise core, an opt-in product at most (defaults are the user's call). B1-B3 hold:
a case for a two-pass seeding inside the assembler, after one unseen sample. B2 holds and B3 fails: the conditioning is
a cleaner unlabelled assembly, not a better one than the labelled one.

## Result, planted (job 4509335; workspace `r1963_two_step_seeding/s*/out.txt`)
In these sets 64-66% of the reads are tied (5,074 of 7,980 at 60 reads per isoform), far more than in a real sample.

Arm A. Exact isoforms (isoform precision); false units in the S sets.
| set | today's layer `primary` | layer `unique` |
|---|---|---|
| 1967 E | 197 (0.824) | 65 (0.867) |
| 1968 E | 197 (0.835) | 65 (0.878) |
| 1969 E | 196 (0.807) | 65 (0.833) |
| 1967 S | 102 (0.650), 8 false units | 32 (0.865), 0 |
| 1968 S | 105 (0.644), 11 | 34 (0.810), 0 |
| 1969 S | 103 (0.652), 6 | 33 (0.917), 0 |
| 1967 low E | 82 (0.953) | 31 (0.969) |
| 1968 low E | 69 (0.932) | 25 (0.862) |
| 1969 low E | 72 (0.935) | 29 (0.906) |
- **A1 HOLDS** (0, 0, 0 false units against 8, 11, 6). **A2 FAILS** on its second half: +0.17 to +0.27 in the S sets as
  claimed, but 0.07 and 0.03 below `primary` in two low-depth sets (bar 0.02).
- Cost: `unique` keeps 31-38% of the `primary` layer's exact isoforms (I predicted 50-60%). The layer equals the step-1
  assembly within one or two transcripts, so it is the assembly of the uniquely placed primaries, nothing more.

Arm B. Exact isoforms (isoform precision); false units in the S sets.
| set | today | conditioned step 2 |
|---|---|---|
| 1967 E | 197 (0.654) | 192 (0.681) |
| 1968 E | 197 (0.654) | 194 (0.681) |
| 1969 E | 197 (0.659) | 194 (0.676) |
| 1967 S | 103 (0.452), 29 | 103 (0.466), 28 |
| 1968 S | 105 (0.417), 34 | 104 (0.444), 30 |
| 1969 S | 104 (0.488), 23 | 101 (0.518), 19 |
| 1967 low E | 166 (0.755) | 164 (0.759) |
| 1968 low E | 162 (0.726) | 158 (0.752) |
| 1969 low E | 168 (0.740) | 167 (0.746) |
- **B1 FAILS** (more than 2 exact isoforms lost in five of nine sets). **B2 FAILS** (false units 28, 30, 19 against 29,
  34, 23; the bar was half). **B3 FAILS** (+0.03 reached in one of nine; the gain is +0.004 to +0.030).
- Why B2 fails: a made-up copy is the identical twin of an expressed one, and at identical twins NO place has a uniquely
  placed read, so the group stays whole exactly where the problem is. Where one place does have unique reads, today's
  assembly was already mostly right.
- **B4 NOT RUN** (a deviation, decided after reading the planted part): with B1-B3 failed the real arm could not change
  the reading, and it costs about two node-hours.
- My predictions: A1 right, A2 wrong on the low-depth half, the cost too optimistic; B1 and B2 wrong, B3 right.

## Reading
Seeding first with the uniquely placed primaries gives a clean core (no made-up copy in three of three half-silent
sets) that holds a third of the isoforms. Bringing the tied reads in afterwards, conditioned on that core, changes
little and loses a few true isoforms: the tied reads that cause wrong copies are the ones for which no place has unique
evidence. Not adopted. What separates those places today is the evidence across reads (the EM labels), not the order of
seeding. Do not re-offer a two-pass seeding conditioned on unique primaries; the `unique` layer could be an opt-in
product (the user's call).
