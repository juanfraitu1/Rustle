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
