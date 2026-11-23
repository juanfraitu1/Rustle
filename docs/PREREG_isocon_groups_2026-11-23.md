# PREREG: IsoCon on the reads of Rustle's tied groups, a hybrid count of expressed copies (row 1925, 2026-11-23)

Row 1924: on a planted sample that differs from the reference IsoCon reproduces the distinct expressed sequences of
reference-identical copies (66 of 67 units), my own counter does not. User: go ahead with the hybrid. Rustle supplies
the groups (loci joined by reads exactly tied on `ms`) and their tied reads; IsoCon runs on each group's reads. No
annotation and no truth decide the pools (every IsoCon arm so far was pooled per family by the truth).

## Built before the test
- `bench/o2_families/isocon_groups.sh ASSEMBLY.gtf BAM AS_TABLE READS.fq GENOME.fa OUT`: pools of >= 5 tied reads per
  group, IsoCon per pool, candidates aligned to the genome (`cand.sam`).
- Scorer `r1925_isocon_groups/score_iso.py`, a COUNT, stricter than row 1924's recall: per (unit of >= 2
  reference-identical copies, isoform) the number of candidates aligned there with exactly that intron chain, against
  the number of distinct sample sequences of that isoform among the unit's expressed copies. A candidate contained in a
  longer one of its pool (25 bases of slack per end) is dropped first as a truncated copy.

## Seen before registering (development, seed 1930 of row 1924, truth-pooled IsoCon only; said so the bars are honest)
Without the containment step IsoCon over-counts (exact 6 of 23 cells; 433 candidates for about 60 sequences): it emits
truncated candidates for the 5'-cut reads. With it: E 11 of 23, S 12 of 19, control E 18 of 23. So even with truth
pools IsoCon is right on about half of the cells as a count; row 1924's "66 of 67" was recall of sequences, not a
count. Constants on the same cells: one per cell 14 of 42, every copy of the unit 20 of 42.

## Test
Fresh seeds 1933, 1934, 1935 (`--private 0.001`, E and S) and a control (seed 1933, no variants). Arms: hybrid
(Rustle's groups), IsoCon with truth pools, the two constants. Same scorer for all.

## Claims (pooled over the six variant sets)
- H1: the hybrid's exact rate is not more than 0.10 under truth-pooled IsoCon's (Rustle's grouping costs little).
- H2: the hybrid's exact rate is above both constants.
- H3: control: the hybrid's exact rate >= 0.70 and it over-counts in no more than 2 cells.
H1 + H2 + H3 = worth offering as an opt-in stage (needs IsoCon installed); otherwise recorded with the counts.

## Predictions (mine)
H1 holds. H2 is a coin toss: on the development seed truth-pooled IsoCon beat "every copy of the unit" by 3 cells of
42. H3 uncertain. The honest expectation is a usable but noisy count (about half of the cells exact), limited by
IsoCon's truncated and spurious candidates rather than by Rustle's grouping.

## Result (job 4496802, workspace `r1925_isocon_groups`; cells with the count exact)

| set | cells | hybrid (Rustle's groups) | IsoCon, truth pools | one per cell | every copy of the unit |
|---|---|---|---|---|---|
| 1933 E | 23 | 13 | 9 | 1 | 18 |
| 1933 S | 23 | 16 | 13 | 14 | 4 |
| 1934 E | 23 | 13 | 9 | 2 | 14 |
| 1934 S | 15 | 9 | 11 | 5 | 4 |
| 1935 E | 23 | 13 | 11 | 0 | 21 |
| 1935 S | 19 | 8 | 7 | 10 | 2 |
| **sum, with variants** | 126 | **72 (0.571)** | 60 (0.476) | 32 (0.254) | 63 (0.500) |
| control E + S (truth = 1 everywhere) | 46 | 37 (0.804), over 0 | 36 (0.783), over 2 | 46 | 0 |

- **H1 HOLDS** (0.571 against 0.476: the hybrid is AHEAD of truth-pooled IsoCon, not just within 0.10). **H2 HOLDS**
  (above both constants, by 9 cells over the better one). **H3 HOLDS** (0.804, no cell over-counted).
- Why the hybrid is ahead of IsoCon given the truth's pools: Rustle's pools hold only the tied reads of one group, so
  IsoCon produces 100-110 candidates per set against 250-260 and fewer spurious ones (mean absolute error 0.39-0.89
  against 0.40-1.00).
- **How good it is, plainly:** the count is exact in 57% of cells, wrong in both directions (over 28, under 26 of 126),
  mean error 0.4-0.9 copies. Against "every copy of the unit" it wins clearly when half the copies are silent
  (33 against 10 of 57 cells) and LOSES when all are expressed (39 against 53 of 69): it is a noisy estimate that
  carries real information about silent copies, not a reliable count.
- Predictions: H1 held, H2 (coin toss) held narrowly, H3 held; "about half of the cells exact" was right.
- By the registered rule this is worth offering as an opt-in stage (needs IsoCon installed). Not wired; the user's
  call. Known before any real use: simulated uniform errors; no heterozygous variants (two alleles of one copy would
  count as two); a third of the candidates align on no planted isoform (intron chain differs) and are not counted.
