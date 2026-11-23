# PREREG: heterozygous copies, do alleles get counted as copies (row 1927, 2026-11-23)

Row 1926: IsoCon's test ported to Rustle's tied groups counts the distinct expressed sequences among reference-identical
copies (exact in 0.80 of cells). Stated limit: two alleles of one copy are two sequences. User: carry out the experiment
for the heterozygous case.

## Said before any run: what can and cannot be expected
A sequence count cannot tell an allele from a copy by itself. One copy with two alleles and two homozygous copies both
give two sequences. The only signal in the reads is balance: the two alleles of a copy are expressed about equally,
two copies need not be. So the experiment measures three things: whether the test still counts sequences correctly
when every copy has two alleles (more sequences, each with half the reads); how far "sequences = copies" is off; and
whether pairing sequences of balanced support recovers the copies.

## Built before any run
- `plant.py --het RATE`: every copy gets a second allele (its private variants plus heterozygous ones at RATE); each
  read is drawn from one of the two alleles at random; `sample.h2.fa`, `het.tsv`. RATE = 0.001 with `--private 0.001`.
- `recurrence_test.py` unchanged in its logic (an `isoform` column added to the output).
- Pairing estimate (scorer `r1927_het_alleles/score_het.py`, uses the tool's table only): inside one locus read and one
  isoform, sequences are paired greedily, closest read counts first, when a two-sided binomial test of their two counts
  against 1 : 1 gives p > 0.05; each sequence in at most one pair; copies = pairs + unpaired sequences.
- Truth per (unit, isoform) cell: haplotypes = distinct sequences over both alleles of the expressed copies; copies =
  distinct allele pairs.

## Data
Human named families, fresh seeds 1939, 1940, 1941; scenarios E (all expressed), S (half silent), D (one copy 10 x,
the others a third: copies of unequal expression, where balance could help).

## Claims (nine sets pooled)
- B1: the tool's sequence count equals the true number of HAPLOTYPES in >= 0.70 of cells.
- B2: taken as a count of COPIES the sequence count is exact in under half of the cells (the limit, measured).
- B3: the pairing estimate equals the true number of copies in >= 0.70 of cells and is above the raw count's rate.

## Predictions (mine)
B1 uncertain, 0.6-0.75: twice the sequences at half the support each. B2 holds. B3 FAILS: with equal expression two
homozygous copies are as balanced as two alleles, so pairing halves real copies; I expect it right in about half of the
cells, better in D than in E. If B3 fails the record says plainly that the count is of sequences and that telling
alleles from copies needs evidence the RNA reads do not hold (the sample's genome, or known heterozygous sites).

## Result (job 4496857, workspace `r1927_het_alleles`; cells with the count exact, three seeds)

| | E, all expressed (69) | S, half silent (61) | D, one copy 10 x (69) | all (199) |
|---|---|---|---|---|
| sequences against true HAPLOTYPES | 39 (0.57) | 42 (0.69) | 25 (0.36) | 106 (0.53); over 36, under 57 |
| sequences against true COPIES | 10 (0.14) | 25 (0.41) | 22 (0.32) | 57 (0.29); over 140, under 2 |
| paired estimate against COPIES | 39 (0.57) | 43 (0.70) | 26 (0.38) | 108 (0.54); over 28, under 63 |
| every copy of the unit, against COPIES | 64 | 2 | 64 | 130 (0.65) |
| one per cell, against COPIES | 0 | 44 | 0 | 44 (0.22) |

- **B1 FAILS** (0.53 against 0.70). With two alleles per copy the count of sequences itself degrades: twice the
  sequences at half the reads each. Worst where copies are weakly expressed (D: the copies at a third have about 20
  reads per isoform, 10 per allele; 34 of 69 cells under-counted). Row 1926's 0.80 holds for homozygous copies only.
- **B2 HOLDS** (0.29): read as copies the sequence count is too high in 140 of 199 cells. The limit, measured.
- **B3 FAILS** (0.54 against 0.70). Pairing balanced sequences does better than the raw count as a copy count (108
  against 57 cells) but is below "every copy of the unit" overall (130) and is right only as often as the sequence
  count is right about haplotypes. It beats that constant only where copies are silent (S: 43 against 2).
- Predictions: B2 and B3 as expected; B1 came in under my range (I said 0.6-0.75); "pairing better in D than in E" was
  wrong (0.38 against 0.57), because D's weak copies are lost before any pairing.
- **What stands.** The opt-in output is a count of SEQUENCES, reliable when copies are homozygous and well covered
  (about 30 reads per sequence), and it must not be read as a number of copies on a heterozygous sample. Telling an
  allele from a copy needs evidence the RNA reads do not hold: the sample's genome or known heterozygous sites. The
  pairing estimate is not built into the tool.
- Not tried, the user's call: support thresholds scaled to coverage; phasing-free use of a DNA call set for the sample.
