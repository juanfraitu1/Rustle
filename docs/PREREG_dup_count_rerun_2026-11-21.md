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

## Result (2026-11-21): no candidate by the registered rule; on the spike-ins every false transcript is gone

Jobs 4491135-4491146, 36 runs exit 0; `r1849_dup_count_rerun/score.NAME.txt`.

`x.final.gtf`: correct chains / precision; share of base's correct chains kept; wrong chains removed per correct chain lost.
| dataset | base | `dc` (fraction 0.02) | | `dc01` (fraction 0.01) | |
|---|---|---|---|---|---|
| fly simulated | 6,973 / 0.968 | 6,753 / 0.973, kept 0.968, 0.2 | no | 6,894 / 0.972, kept 0.989, 0.4 | QUALIFIES |
| mouse simulated | 26,894 / 0.959 | 25,077 / 0.970, kept 0.932, 0.2 | no | 26,090 / 0.967, kept 0.970, 0.3 | no |
| zebrafish simulated | 12,318 / 0.896 | 12,129 / 0.912, kept 0.985, 1.4 | QUALIFIES | 12,245 / 0.910, kept 0.994, 3.0 | QUALIFIES |
| mouse ES | 17,663 / 0.371 | 17,074 / 0.410, kept 0.967, 9.0 | no | 17,614 / 0.363, kept 0.997, +945 wrong | no |
| mouse heart | 13,077 / 0.591 | 12,923 / 0.610, kept 0.988, 5.1 | QUALIFIES | 13,034 / 0.595, kept 0.997, 3.7 | QUALIFIES |
| mouse adrenal | 14,363 / 0.593 | 14,163 / 0.615, kept 0.986, 4.9 | QUALIFIES | 14,340 / 0.592, kept 0.998, +13 wrong | no |
| mouse cortex | 15,647 / 0.501 | 15,491 / 0.520, kept 0.990, 8.1 | QUALIFIES | 15,634 / 0.501, kept 0.999, +28 wrong | no |
| OR6737 | 25,530 / 0.384 | 24,315 / 0.421, kept 0.952, 6.1 | no | 25,432 / 0.388, kept 0.996, 8.6 | QUALIFIES |
| A119b | 37,804 / 0.223 | 34,882 / 0.236, kept 0.923, 6.3 | no | 37,521 / 0.222, kept 0.993, +408 wrong | no |
| chimpanzee | 20,789 / 0.419 | 20,309 / 0.432, kept 0.977, 4.5 | no | 20,738 / 0.421, kept 0.998, 6.3 | QUALIFIES |
| bonobo | 17,832 / 0.341 | 17,383 / 0.358, kept 0.975, 7.4 | no | 17,781 / 0.342, kept 0.997, 5.8 | QUALIFIES |
| KB3781 | 25,577 / 0.384 | 23,407 / 0.440, kept 0.915, 5.2 | no | 25,048 / 0.400, kept 0.979, 6.6 | no |
- **`dc`: 4 of 7 non-primate, 0 of 5 primate. `dc01`: 3 of 7 and 3 of 5 (not A119b). Neither is a candidate.**
- `dc` raises precision on all twelve (by 0.005-0.056) and keeps 91.5-99.0% of the correct chains; under 98% on eight. On the
  real libraries it removes 4.5-9 wrong chains per correct one lost; on fly and mouse simulated it loses five correct per wrong.
- `dc01` keeps 97.0-99.9% but its precision is LOWER than base on mouse ES, adrenal, cortex and A119b in `final.gtf` (it is
  higher there in `x.gtf` on all but ES): at fraction 0.01 the final filter, fitted on collapsed counts, lets more through.
- Single-exon models (the floor on distinct records applies to them too): truth found falls on every dataset, steeply on the
  simulated sets (fly 344 -> 97, mouse 89 -> 30, zebrafish 148 -> 52), by 4-17% on real libraries; the same in both arms.
- Expectation: `dc` as written. `dc01`: mouse simulated misses as expected (3.0% lost); mouse ES was not repaired by the new
  defaults, it got worse (precision -0.008).

**SIRV contigs of the mouse ES library (development data):**
| arm | `x.gtf`: reported / correct / wrong | `final.gtf`: reported / correct / wrong | all transcripts incl. single-exon, `final.gtf` |
|---|---|---|---|
| base | 75 / 55 / 20 | 62 / 49 / 13 | 77 |
| `dc` | 50 / 50 / 0 | 48 / 48 / 0 | 63 |
| `dc01` | 50 / 50 / 0 | 48 / 48 / 0 | 63 |
With every read counted, on today's defaults, **no false spliced SIRV transcript is left**, before or after the filter; one
correct one is lost in `final.gtf` (49 -> 48) and five before the filter (55 -> 50). The 15 single-exon spike-ins stay.
Rows 1836 / 1837 had 0.877 precision here; the junction replacement and the spliced share removed what was left.

**Decision:** no default proposed by the rule. The switch stays opt-in. For the user: `RUSTLE_DUP_COUNT=1` (fraction 0.02) is
the setting that gives a clean spike-in set and higher precision on every library, at 1-8.5% of the correct chains (most on
KB3781, A119b, mouse simulated) and a loss of single-exon transcripts. `dc01` is not recommended: it gives the precision back.
