# Row 1879: near-miss junctions judged by the edit distance of the two transcripts (step 0, pre-registered 2026-11-21)

Row 1842 found that near misses differ by distance, not by side, and named two small cells: a weaker chain 1-2 bp from a stronger
one, and an intron moved at both sides. Both have one reading in common that row 1842 did not test: **the two spliced sequences
differ by very few bases.** A 1 bp shift at one site makes transcripts one base apart; an intron moved by the same amount at both
sides makes transcripts of equal length that differ by a few substitutions or by nothing. One read error at that base lets the
aligner choose the other placement. The 3 bp acceptor shift (NAGNAG) is a 3-base indel, which a read error rarely makes. The
junction replacement of rows 1819 / 1845 is blind to this: it asks for 20 x the support at any distance up to 10 bp.

This row measures only. Workspace `r1879_junction_edit` (`edit.py`). Registered before the script ran.

## Data
Today's default assemblies, no new run: row 1849 arm `dc` for the twelve datasets (fly, mouse, zebrafish simulated; mouse ES on
mouse contigs, heart, adrenal, cortex; OR6737, A119b, chimpanzee, bonobo, KB3781) and `r1860_polya_ends/wtc` (WTC11 r2). All are
development data for this lever (rows 1819, 1842, 1845 used the twelve).

## Measure, fixed now
- **Pair:** two chains of `x.gtf` on one strand whose intron lists are equal except at one position, where the two introns lie
  within 10 bp at each side. WEAKER = the chain with fewer `reads` (ties: no pair).
- **Edit distance D** of the pair: Levenshtein distance between the two spliced sequences over a common window
  (20 bp of exon before the leftmost donor to 20 bp after the rightmost acceptor, the intron of each chain cut out).
- Cells: D = 0, 1, 2, 3, 4-5, 6+. Per cell: weaker chains, correct, how many are also in `final.gtf`, and the weaker / stronger
  read ratio in four bins (< 0.05, 0.05-0.2, 0.2-0.5, >= 0.5). D = 3 is printed apart for "one site moved by 3 bp".
- A weaker chain with several stronger partners is counted once, at its smallest D.

## Reading, fixed now
- A cell is an **ERROR CELL** on a dataset when it holds at least 20 weaker chains and their share correct is under half the
  precision of that dataset's `x.gtf`.
- **Hypothesis H1:** the cells D <= 2 taken together are an error cell on at least 10 of the 13 datasets, in both groups
  (non-primate at least 5 of 7, primate at least 4 of 6), and stay so at a read ratio >= 0.05 (the pairs ratio 20 does not reach).
- **Hypothesis H2:** D = 3 is NOT an error cell on at least 10 of 13 (it must be spared).
- H1 and H2 both hold -> a rule is specified in its own row: replace the weaker chain's intron by the stronger one's when D <= 2,
  whatever the ratio above 1, as part of the junction replacement. Its size is stated from this table before it is built.
- H1 fails -> no rule; the 20 x replacement stays the only one.

## Limits, stated in advance
- Pairs are found among emitted chains only: a wrong chain whose true partner is not emitted is not seen (row 1842, table E:
  the partner is emitted in 81-94% on real libraries).
- Real-library truths are incomplete; a weaker chain one base from an annotated one is rarely a real isoform, but "wrong" there
  means "not annotated".
- The window is 20 bp: D saturates for pairs that differ more than the cells need.

## Result (2026-11-21): H1 FAILS (9 of 13), H2 holds; no rule, and the class is nearly empty

Job 4494583, thirteen runs exit 0; `r1879_junction_edit/edit.NAME.txt`. Weaker chains of `x.gtf` pairs, correct / all.
| dataset | `x.gtf` precision | D <= 2 | D <= 2 at ratio >= 0.05 | error cell | D = 3 | D = 4-5 |
|---|---|---|---|---|---|---|
| fly simulated | 0.941 | 0 / 8 | 0 / 8 | few | 12 / 17 | 9 / 13 |
| mouse simulated | 0.927 | 16 / 22 | 15 / 21 | no | 289 / 314 | 96 / 112 |
| zebrafish simulated | 0.839 | 2 / 14 | 2 / 14 | few | 16 / 21 | 8 / 21 |
| mouse ES | 0.389 | 5 / 73 | 5 / 71 | yes | 434 / 1,957 | 107 / 591 |
| mouse heart | 0.579 | 3 / 33 | 3 / 33 | yes | 251 / 698 | 43 / 183 |
| mouse adrenal | 0.576 | 2 / 37 | 2 / 37 | yes | 273 / 762 | 68 / 231 |
| mouse cortex | 0.497 | 3 / 48 | 3 / 48 | yes | 320 / 1,194 | 74 / 343 |
| OR6737 | 0.356 | 3 / 56 | 3 / 55 | yes | 1,037 / 1,883 | 237 / 634 |
| KB3781 | 0.392 | 8 / 62 | 8 / 61 | yes | 1,141 / 1,942 | 241 / 642 |
| A119b | 0.174 | 20 / 224 | 19 / 223 | no (0.089 against 0.087) | 1,200 / 3,865 | 349 / 1,668 |
| chimpanzee testis | 0.368 | 3 / 39 | 3 / 39 | yes | 652 / 1,070 | 165 / 393 |
| bonobo testis | 0.328 | 3 / 46 | 3 / 46 | yes | 683 / 1,264 | 137 / 448 |
| WTC11 r2 | 0.595 | 5 / 41 | 5 / 41 | yes | 670 / 983 | 174 / 346 |

- **H1 fails as registered:** D <= 2 is an error cell on 9 of 13 (non-primate 4 of 7, primate 5 of 6). The three simulated sets
  hold 8, 22 and 14 such chains, and on mouse simulated 16 of the 22 are true.
- **H2 holds:** D = 3 is an error cell on 0 of 13 (it is the one-site 3 bp shift in 99-100% of its chains).
- **The class is nearly empty.** After the default snap and the 20 x replacement, 33-73 weaker chains per real library are left
  at D <= 2 (224 on A119b): 0.1-0.3% of the wrong chains of `x.gtf`. On every real library they are 5-14% correct, so the
  reading of row 1842 stands (a one- or two-base difference is rarely a real isoform), but there is nothing left to gain: a rule
  would remove at most 28-68 wrong chains per library for 2-8 annotated ones (A119b 204 for 20).
- D = 0 (two placements of one sequence) occurs twice in thirteen datasets.
- Not built. **Near-miss junctions are closed as a lever** by rows 1819, 1835, 1842, 1845 and this one: what is left near a
  stronger junction is the 3 bp cell (0.22-0.36 correct on mouse libraries, 0.31-0.68 on primates, 0.71-0.92 simulated) and
  wider shifts at 0.18-0.50, which no table here separates from the rest of the assembly by a factor of two.
