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

## Result (2026-11-20): not a candidate; 0.01 is the rule switched off, and a recall trade that pays on 3 of 12

Jobs 4489087-4489098, 48 arms exit 0; `r1839_subchain_ratio_001/score.NAME.txt`. Repo `509e0a5f`, binaries of `~/rustle_target`
(rebuilt 2026-11-20 from `3a1005d1`). `base` reproduces row 1837's `base` on every dataset.

**0.01 against 0:** the same output (byte-identical `x.gtf` on mouse simulated; on A119b a different file with the same
counts). A two-read chain passes a ratio of 0.01 unless the longer chain has more than 200 reads, so 0.01 is the rule off.

`x.final.gtf`, arm `r0.01` against `base` (ratio 0.3):
| dataset | base correct / precision | correct added | wrong added | precision | marginal precision | pays |
|---|---|---|---|---|---|---|
| fly sim | 7,041 / 0.964 | +14 (0.2%) | 0 | 0.964 | 1.000 | yes |
| mouse sim | 27,159 / 0.953 | +515 (1.9%) | +4 | 0.954 | 0.992 | yes |
| zebrafish sim | 12,381 / 0.891 | +26 (0.2%) | +26 | 0.890 | 0.500 | no |
| mouse ES | 17,730 / 0.347 | +276 (1.6%) | +920 | 0.345 | 0.231 | yes (second clause) |
| mouse heart | 13,091 / 0.577 | +91 (0.7%) | +162 | 0.574 | 0.360 | no |
| mouse adrenal | 14,394 / 0.578 | +170 (1.2%) | +434 | 0.571 | 0.281 | no |
| mouse cortex | 15,679 / 0.487 | +137 (0.9%) | +257 | 0.485 | 0.348 | no |
| OR6737 | 25,635 / 0.377 | +43 (0.2%) | +1,532 | 0.370 | 0.027 | no |
| A119b | 38,091 / 0.219 | +111 (0.3%) | +1,730 | 0.217 | 0.060 | no |
| chimpanzee testis | 20,847 / 0.416 | +29 (0.1%) | +1,318 | 0.405 | 0.022 | no |
| bonobo testis | 17,884 / 0.337 | +37 (0.2%) | +1,397 | 0.329 | 0.026 | no |
| KB3781 | 25,745 / 0.373 | +39 (0.2%) | +2,328 | 0.361 | 0.016 | no |
Pays on 3 of 7 non-primate and 0 of 5 primate datasets (bar 6 and 4): **a recall trade, not a candidate.** Arm `r0.1` pays on
fly and mouse simulated only; it adds about 60-70% of what 0.01 adds on the mouse libraries (ES +159 / +443 wrong, adrenal
+108 / +247) and on the apes +20 to +108 correct for +856 to +1,603 wrong.

- **The rule does act differently on the two groups.** Without it a mouse library gains 0.7-1.9% correct chains at 2-3 wrong
  chains per correct one; an ape library gains 0.1-0.3% at 16-60 wrong per correct one. Row 1838's pointer was right in
  direction: the rule costs the non-primates more. It is still worth its cost on real mouse libraries by the registered reading.
- **Before the `final.gtf` filter the cost is much larger** (`x.gtf`, 0.01): mouse ES +1,347 correct / +24,399 wrong (0.331 ->
  0.242), mouse simulated +2,287 / +2,017 (0.903 -> 0.857), A119b +569 / +94,352. The filter removes most of what the ratio
  lets in, wrong and correct alike (ES keeps 276 of the 1,347).
- Expectation: held for the real libraries; wrong for mouse simulated, where 0.01 adds two thirds more than 0.1 (515 against 309).

Reading (mine): the read ratio does not tell a true contained chain from a truncated copy of the longer one; lowering it buys
recall at the rate the data allows. The rule from IsoQuant noted in row 1838 (keep a contained chain when its own end is a
read-end peak) is the test that could separate them; these arms are its reference.
