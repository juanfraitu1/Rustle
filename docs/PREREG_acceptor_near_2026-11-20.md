# Row 1842: near misses at splice sites, acceptor against donor (pre-registered 2026-11-20)

User (2026-11-20): look at the rules that judge near misses at the intron acceptor sites. Measurement only, GTF level; no
assembler code, no default. Workspace `r1842_acceptor_near` (`near.py`). Committed before the script was written or run.

## What exists (code reading, 2026-11-20)
No rule treats the acceptor apart from the donor. Every near-miss rule uses one window of 10 bp on both sides:
| rule | state | what it does |
|---|---|---|
| `--polish-junction-snap reads` | default (row 1800) | junction K -> neighbour J (donor and acceptor within 10 bp, same locus) when fewer than 2 of K's reads align better at K than at J |
| `--polish-junction-ratio R` | opt-in (rows 1819, 1835; 20 proposed) | K -> J within 10 bp when J has R times K's support |
| `--polish-fuzzy-junction`, `--read-junction-snap` | opt-in, failed (rows 866, 1527, 1618) | merges within a window, no evidence test |
| chain score features (wobble class, non-GT-AG count) | in the `final.gtf` filter | no side |
| support profile flags | failed as a filter (row 1830) | shifts of 1-3 bp, no side |
Known: 18-27% of the wrong chains of `final.gtf` have an intron within 10 bp of a true one (row 1818); Rustle keeps 4-7 times
StringTie's near-junction wrong chains on real libraries (row 1838); weaker shifted-site chains were measured on apes only,
without side or distance (rows 1725-1729: frame-keeping shifts are more often annotated); splice-site sequence strength
failed as a score feature (row 1789). Not measured anywhere: the split by side and by distance, on non-primates.

## Measured (`x.final.gtf` of the row 1839 `base` arms, twelve datasets, truths of row 1839; `x.gtf` beside it)
A PAIR = two emitted chains of one contig and strand with the same number of introns, equal except in ONE intron whose
donor and acceptor each differ by at most 10 bp. The WEAKER chain is the one with fewer `reads` (ties: both left out).
Cell of a pair = side (acceptor only, donor only, both) x distance class (1-2, 3, 4-5, 6 or 9, 7-8 or 10 bp; for `both` the
larger of the two). Per cell: weaker chains, the share of them that is correct, their share with reads ratio < 0.05 to the
stronger chain, and whether the stronger chain is correct. For acceptor-only pairs also: is the weaker chain's acceptor an AG.
Beside it, table E: the wrong chains that differ from a TRUTH chain in one intron by <= 10 bp, by the same cells (the errors,
whether or not the true partner was emitted).

## Reading, fixed now
- **The sides differ** when, on at least 9 of 12 datasets, the weaker chains of acceptor-only pairs are at least twice or at
  most half as precise as those of donor-only pairs (all distances pooled). Otherwise a rule has no reason to treat them apart.
- **A cell is a rule candidate** when it holds at least 50 weaker chains and they are at most a quarter as precise as the
  dataset's `final.gtf`, on at least 6 of 7 non-primate and 4 of 5 primate datasets. A candidate goes to the user as a
  proposal for a registered rule row (drop or re-place the weaker chain of that cell); nothing is built here.
- **A cell is protected** (a rule must spare it) when its weaker chains are at least as precise as the dataset's `final.gtf`
  on at least 6 of 7 and 4 of 5.

## Expectation, written before the script
Acceptor-only pairs outnumber donor-only ones and the 3 bp acceptor cell (NAGNAG) is the most precise one, likely protected on
the simulated sets. Distances 1-2 and 4-5 at the acceptor are the artifact cells; I am unsure they reach "a quarter" on real
libraries, since rows 1725-1727 found the weaker shifted chains only mildly worse than the read count predicts on human.

## Limits
- Pairs differ in exactly one intron; chains with two shifted sites are not seen. The default snap has already acted, so the
  pairs are what it left.
- Real-library truths are incomplete: an unannotated true alternative acceptor counts as wrong.
- Simulated reads: junction errors come from the simulator and the aligner. Twelve reused datasets, said.
