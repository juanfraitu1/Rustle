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
