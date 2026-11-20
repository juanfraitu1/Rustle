# Row 1849: every read counted, re-run on today's defaults (pre-registered 2026-11-21)

User (2026-11-21): re-run the read-counting switch on top of today's defaults. Rows 1836 / 1837 measured `RUSTLE_DUP_COUNT=1`
(chains admitted on distinct coordinate records, `reads` = every record) before the junction replacement at 20 x support (row
1845), the spliced share (row 1848) and the end layers going opt-in (row 1846). Then: at isoform fraction 0.02 precision up on
most libraries and A119b -8.1% correct chains; at 0.01 about neutral, missing on mouse ES and mouse simulated. Workspace
`r1849_dup_count_rerun`. Committed before any run.

## Arms and data
`base` (today's defaults), `dc` (`RUSTLE_DUP_COUNT=1`), `dc01` (the same with `--polish-isoform-fraction 0.01`, by `sed` on the
driver snapshot as in row 1837: the driver has no switch for it) on the twelve datasets of rows 1837-1848. Stage `assemble`,
`RUSTLE_PTC=0`, `RUSTLE_ONEREAD_TIER=0`, binaries `~/rustle_target`. Scorer of row 1831. SIRV contigs of the mouse ES library
shown apart (development data). Reused datasets, said; no arm is new, so nothing is fitted here.

## Reading, fixed now (`x.final.gtf`; `x.gtf` beside it)
- An arm **qualifies** on a dataset when chain precision is higher than `base` and it keeps at least 98% of `base`'s correct
  chains (row 1837's rule).
- **Candidate** = an arm qualifies on at least 6 of 7 non-primate and 4 of 5 primate datasets, A119b among them (the library
  the collapse protects, row 1833). A candidate goes to the user; no default changes here.
- Reported beside it: wrong chains removed per correct chain lost, and single-exon models (the floor on distinct records
  also applies to them).

## Expectation
`dc`: precision up everywhere, more than 2% of the correct chains lost on A119b and several others: not a candidate.
`dc01`: close to neutral; I expect it to miss on mouse simulated (3.4% lost in row 1837) and to be decided by whether the two
new defaults changed mouse ES. On SIRV both should remove most of the 13 remaining shifted-site chains (row 1837: 0.877
precision in `x.gtf`).

## Limits
- The chain score filter and the 5' truncation rule were fitted on collapsed counts; `final.gtf` under real counts is outside
  that fit (this is why `final.gtf` is the judged file).
- Real-library truths are incomplete.
