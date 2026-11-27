# PREREG: another test of "borrowing across the family": a thin copy beside a deep sibling (row 1980, 2026-11-26)

The user's goal: find another way to test the idea of borrowing information across the family for the family-wise
assembly. What the earlier tests were, and why this one is different:
- Sibling models confirming one-read chains (`siblings`, rows 1627 / 1628: quality passes, yield fails), sibling models
  as a guide (row 1632, fails), borrowed structure on KB3781 (row 1822, fails), pooling reads (closed), relatives' reads
  for thin isoforms genome-wide (row 1965: 0.2%). All were judged over WHOLE samples or families, where most copies are
  either well covered (nothing to gain) or have no well-covered relative (nothing to borrow).
- The planted low-depth sets (rows 1917, 1952) give EVERY copy 3 reads: there is no deep sibling to borrow from.
Borrowing can only pay in one situation: a thinly expressed copy next to a deeply expressed sibling. That situation was
never isolated. This test builds it and measures the thin copies alone.

## Part A, planted (exact truth)
New seeds 1981-1983, CHM13 families, ragged ends. Two sets per seed:
- **beside a deep sibling:** `plant.py --scenarios D --reads 9`: one copy per family gets 90 reads, the others 3 each;
- **control, no deep sibling:** `--scenarios E --reads 3`: every copy 3 reads.
THIN copies = those with 3 planted reads (a copy's reads are split 2 : 1 over two isoforms, so one isoform has 2 reads
and one has 1: the second can only be assembled with support from elsewhere).
Arms: Rustle on primaries only (no borrowing); today's mode (tied secondaries seed, 0.98 of the best AS); the strict
relatives rule of row 1961 with one own read required (`janchor:ms:0.9:1`); the open pool at 0.90 (`as:0.90`).
`bench/o2_truth/thin_copy_score.py`: at the thin copies only, planted isoforms found exactly (recall, also split into
the 2-read and the 1-read isoform), transcripts built there, and the share of them that are planted isoforms.

## Claims (thin copies only)
- **B1, borrowing happens:** with today's mode, isoform recall at thin copies beside a deep sibling is at least 10 points
  above the same arm's recall at thin copies of the control, in all three seeds.
- **B2, it is the tied reads:** beside a deep sibling, today's mode is at least 10 points above primaries-only, in all three.
- **B3, the one-read isoform:** beside a deep sibling at least 20% of the 1-read isoforms are found by today's mode
  (primaries-only cannot find one).
- **B4, cost:** the share of transcripts at thin copies that are planted isoforms is not more than 10 points lower
  beside a deep sibling than in the control, for today's mode.
## Predictions (mine)
B1 and B2 hold (a deep sibling's tied reads land on the thin copy and carry its shared isoforms over the floor).
B3 holds. B4 fails: the deep copy's isoforms that the thin copy does not make get built there too.
If B1-B3 hold, borrowing is real in exactly this corner, and Part B (real reads: a deep copy thinned to 2 reads beside
its untouched sibling, truth = its own full-depth assembly) is worth running.

## Result, Part A (jobs 4512518-4512520). Thin copies only (3 planted reads; 118 beside a deep sibling, 133 in the control).
Isoform recall at thin copies / of the isoforms with ONE planted read / share of the transcripts built there that are planted isoforms:

| arm | beside a deep sibling: seeds 1981 / 1982 / 1983 | control, every copy thin |
|---|---|---|
| primaries only | 0.371 / 0.337 / 0.381; 7% / 10% / 9%; 0.81 / 0.83 / 0.83 | 0.364 / 0.307 / 0.329; 7% / 6% / 7%; 0.94 / 0.90 / 0.93 |
| StringTie | 0.322 / 0.322 / 0.337; 7% / 11% / 10%; 0.73 / 0.72 / 0.72 | 0.320 / 0.307 / 0.325; 8% / 9% / 13%; 0.81 / 0.75 / 0.78 |
| today's mode (tied reads seed) | 0.728 / 0.717 / 0.703; 63% / 59% / 55%; 0.68 / 0.72 / 0.70 | 0.711 / 0.728 / 0.724; 56% / 60% / 61%; 0.77 / 0.70 / 0.75 |
| strict relatives rule, one own read | 0.738 / 0.707 / 0.698; 64% / 56% / 53%; 0.72 / 0.74 / 0.71 | 0.706 / 0.706 / 0.732; 55% / 56% / 62%; 0.81 / 0.74 / 0.78 |
| open pool 0.90 | 0.738 / 0.698 / 0.713; 63% / 52% / 56%; 0.65 / 0.69 / 0.70 | 0.706 / 0.715 / 0.737; 55% / 58% / 61%; 0.76 / 0.72 / 0.71 |

### Claims
- **B2 HOLDS:** at thin copies today's mode is 33 to 38 points above primaries-only (0.70-0.73 against 0.34-0.38).
- **B3 HOLDS:** 55-63% of the isoforms that have ONE read of their own are found (primaries-only 7-10%).
- **B1 FAILS:** the recall is the same with and without a deep sibling (0.70-0.73 against 0.71-0.73). I predicted a gap.
- **B4 HOLDS** (I predicted it would fail): the share of planted isoforms among the transcripts built is 9, -2 and 5
  points lower beside a deep sibling.

### Reading
- **Borrowing across the family is real and large at thin copies:** it doubles the isoforms recovered there and takes
  the isoforms with a single read of their own from under a tenth to about six in ten. This is the family-wise gain
  isolated on the copies where it can exist; over whole samples it was diluted by the copies that need nothing.
- **It does not need a deeply expressed sibling.** The help comes from the family as a whole: every copy's few reads
  tie with its relatives, and the tied reads together carry a shared isoform over the two-read floor at each copy.
- **Its price is at the same copies:** the share of built transcripts that are planted isoforms falls from 0.81-0.94 to
  0.68-0.77. The strict rule of row 1961 gives the same recall for 2-4 points less of that loss.
- Limits before this is quoted: planted reads, sample = reference, families chosen for being large and near-identical
  (the most favourable case for borrowing), isoforms shared by construction where the annotation has them at each copy.
  Part B (real reads: a deep copy thinned beside its untouched relatives, truth = its own full-depth assembly) is run next.

## Part B, real reads (registered after Part A, before any run of it)
OR6737 (gorilla). TARGETS = a random half (seed 1980) of the family copies of today's run (`r1937_shared_floor1/or`)
that have at least 20 reads placed there without a tie-break (primary, mapping quality >= 1, >= 50 exon bases).
The reads whose primary is at a target are thinned to 3 per target (the others are removed with all their alignments);
everything else is untouched, so the targets' relatives keep their depth (`bench/o2_truth/thin_targets.py`).
TRUTH for a target = the spliced intron chains a primaries-only assembly of the FULL reads builds there: what the
copy's own reads show when there are enough of them. No annotation, no borrowing in the truth.
Arms on the thinned reads: primaries only; today's mode; the strict relatives rule with one own read.
`bench/o2_truth/thin_real_score.py`: per arm, truth chains found at the targets, chains built there, and the share of
the built chains that are in the truth.
- **RB1:** today's mode finds at least 10 points more of the truth chains than primaries-only.
- **RB2:** the share of built chains that are in the truth is at most 15 points lower for today's mode than for
  primaries-only.
Predictions (mine): RB1 holds by a smaller margin than on planted reads (10-20 points: real copies differ more and
share fewer isoforms); RB2 fails.

## Result, Part B (jobs 4512538 / 4512539)
| | OR6737 | KB3781 |
|---|---|---|
| candidate copies / targets | 544 / 272 | 414 / 207 |
| reads placed at the targets, before / kept | 49,483 / 873 | 110,977 / 651 |
| truth chains at the targets (full reads, primaries only) | 1,439 | 1,175 |
| primaries only on the thinned reads: truth chains found / built / share in the truth | 126 (8.8%) / 139 / 0.906 | 116 (9.9%) / 124 / 0.935 |
| today's mode | 286 (19.9%) / 575 / 0.497 | 242 (20.6%) / 425 / 0.569 |
| strict relatives rule, one own read | 289 (20.1%) / 634 / 0.456 | 257 (21.9%) / 455 / 0.565 |
- **RB1 HOLDS:** today's mode finds 11.1 and 10.7 points more of the truth chains than primaries-only (2.3 and 2.1
  times as many: 160 and 126 chains that the copy's own three reads do not show and its full reads do).
- **RB2 FAILS** (as predicted): the share of built chains that are in the truth falls by 41 and 37 points.

### Reading of both parts
- **Borrowing across the family does what the advisor expects, at the copies where it can:** a thinly covered copy gets
  back twice as many of its own isoforms from its relatives' tied reads, on planted reads (0.35 -> 0.72 of the isoforms)
  and on real ones (0.09 -> 0.20 of what the full reads show). Real reads confirm the direction with a smaller size, as
  predicted: with three reads against a truth built from 180-540, most of the truth is out of reach of any method.
- **It costs precision at the same copies, and more on real reads:** about half of what is built at a thinned real copy
  is not in that copy's own full-depth assembly (planted: 0.68-0.77 against 0.81-0.94). Those chains are either a
  relative's isoforms that this copy does not make, or isoforms it makes below the truth assembly's own thresholds; the
  test does not separate the two.
- **So the effect is a trade, concentrated where depth is low,** which is why whole-sample tests saw almost nothing: the
  gain is real only for the thin copies, and they are a minority of a sample's transcripts. The layer tag marks exactly
  these transcripts (`tied_support`, `tied_only`).
- The strict rule of row 1961 does not help here: same recall, no better share.
- Limits: the truth is one assembler's primaries-only output at full depth, with its own floor and filters; thinning to
  three reads is one depth; targets are copies with at least 20 own reads (well-placed copies, not the hardest arrays).
