# Row 1839: sub-chain collapse ratio at 0.01 (pre-registered 2026-11-20)

User (2026-11-20, after row 1838): try the sub-chain ratio at 1%. Row 1838 found that true chains StringTie / IsoQuant have and
Rustle lacks are more often contiguous parts of a longer Rustle chain on the non-primate datasets (13-27% of the missed chains
against 2-7% on apes). The ratio R keeps such a part when it has at least R times the longer chain's reads; default 0.3 (row
1810), tested before: 0.7, 0.1 and 0 (row 1800, at earlier defaults), 0.15 (row 1837, with other switches). Workspace
`r1839_subchain_ratio_001`. Committed before any arm ran.

## Arms (`assemble`, `RUSTLE_PTC=0`, binaries of `~/rustle_target` at the commit of this document)
`base` (ratio 0.3), `r0.1`, `r0.01` (`RUSTLE_POLISH_ISM_RATIO`), `r0` (the rule off, as the ceiling). Nothing else changes.

## Datasets (the twelve of row 1837, all used before; reuse per the user's ruling)
Non-primate: fly, mouse, zebrafish simulated; mouse ES, heart, adrenal, cortex. Primate: OR6737, A119b, chimpanzee and bonobo
testis, KB3781. Truths and scorer as row 1837 (`score.py`, exact intron chains; `x.gtf` and `x.final.gtf`).

## Reading, fixed now (arm against `base`, on `x.final.gtf`; `x.gtf` beside it)
- Lowering the ratio can only add chains (and through the filter move a few). Per dataset: correct chains added, wrong chains
  added, precision before and after, and the **marginal precision** = correct added / chains added.
- An arm *pays* on a dataset when its marginal precision is at least the dataset's `base` precision (the added chains are no
  worse than the output they join), or when precision falls by at most 0.005 while correct chains rise by at least 1%.
- **Candidate** = pays on at least 6 of 7 non-primate and 4 of 5 primate datasets. Otherwise a recall trade, reported with
  what one correct chain costs in wrong chains per dataset.
- Share of the ceiling: correct chains gained at 0.01 as a share of those gained at ratio 0.
- A candidate is not a default; the user decides.

## Expectation, written before the run
Row 1800 (earlier defaults, `x.gtf`): 0.3 -> 0.1 added 0.7-3.4% correct chains and cost 0.02-0.04 precision on real libraries.
I expect 0.01 to add little beyond 0.1 on the simulated sets and to cost precision on every real library: a recall trade, not a
candidate.

## Limits
Simulated reads lack degraded molecules, the main source of sub-chains, so they understate the cost. On real libraries the
truths are incomplete: a kept part can be a real unannotated transcript. The row-1838 pointer came from six of these datasets.
