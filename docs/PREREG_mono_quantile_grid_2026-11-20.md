# Row 1832: the read floor of single-exon transcripts, a grid of its quantile (pre-registered 2026-11-20)

User (2026-11-20), after row 1831 (`docs/PREREG_mono_end_rules_2026-11-20.md`): run the grid. Row 1831: the floor
(`--polish-mono-quantile 0.82`: a single-exon transcript needs the reads of the 0.82 quantile of the run's multi-exon
transcripts) costs 41-90% of the true single-exon transcripts found without it, on mouse, fly, zebrafish and the apes alike.
Registered before any arm ran. Workspace `r1832_mono_quantile_grid`. The floor acts on `PREFIX.gtf` only (`final.gtf` rebuilds
its single-exon models), so arms run with `RUSTLE_ASSEMBLY_LAYERS=off`; no default changes in this row.

## Step 1, the grid (the seven datasets of row 1831: diagnostic for this lever since that row)
Arms: quantile 0.82 (default), 0.75, 0.6, 0.5, 0.4, 0.25, 0.1, 0; stage `assemble`, `RUSTLE_PTC=0`. Scorer and measures of row
1831 (`r1831_mono_end_rules/score.py`): single-exon models, model precision, truth single-exon transcripts found; multi-exon
chains per arm (expected unchanged). Reported with them: the floor in reads each quantile gives per dataset (the smallest `reads`
of a single-exon model), since many quantiles map to the same floor.

**Choice, fixed now.** Non-primate datasets N = fly, mouse, zebrafish simulated, mouse ES, mouse heart.
- A quantile *pays* on a dataset against quantile 0 when single-exon model precision is higher and at most 5% of quantile 0's
  found truth transcripts are lost (row 1793's reading).
- q* = the largest grid value that pays on at least 4 of 5 N.
- If no value pays: q_b = the largest grid value that keeps at least 80% of quantile 0's found truth transcripts on at least 4
  of 5 N and has higher precision than quantile 0 on both real N libraries (mouse ES, heart). If none either: no candidate,
  the table goes to the user.
- OR6737 and A119b are reported beside the choice, not part of it (the lever was set on A119b).

## Step 2, check on a dataset not used for this lever
Mouse adrenal, ENCODE ENCSR852HRM (`reads/MMU/MMU_adrenal_mm.bam`, GENCODE vM27, mouse contigs). Its arms run with the grid but
are not scored until the choice of step 1 is written into this document.
- **C1:** the chosen quantile finds at least 10% more truth single-exon transcripts than 0.82.
- **C2:** its single-exon model precision is higher than quantile 0's.
Both pass = candidate default; the user decides. Either fails = reported, nothing proposed.

## Limits, stated in advance
- Simulated reads carry no intronic or antisense unspliced signal: on them a lower floor looks nearly free. The real libraries
  decide whether precision collapses, and there annotated single-exon genes are few, so unmatched models may be real.
- `PREFIX.gtf` feeds the family and assignment stages; this row measures the assembly only, not what more single-exon
  transcripts do there (the single-exon seeding trade of rows 1518 / 1520 is a separate lever).
- The floor is a quantile of each run's own support distribution: the same quantile is a different read count per library.
- Adrenal is the same species as two N datasets (another tissue and library).

## Result (2026-11-20): no quantile qualifies; the floor is a precision / recall dial, not a primate fit

Jobs 4484165-4484171 (grid), 4484172 (adrenal arms, run, NOT scored); every arm exit 0; tables
`r1832_mono_quantile_grid/score.NAME.txt`. Repo `197cdfeb`, binaries of 2026-11-20. Multi-exon chains are identical in all arms.

`PREFIX.gtf`, truth single-exon transcripts found / single-exon model precision:
| quantile | fly sim | mouse sim | zebrafish sim | mouse ES | mouse heart | OR6737 | A119b |
|---|---|---|---|---|---|---|---|
| 0 | 1,800 / 0.877 | 850 / 0.958 | 927 / 0.614 | 532 / 0.032 | 221 / 0.144 | 211 / 0.057 | 396 / 0.011 |
| 0.1 | 1,792 / 0.878 | 840 / 0.957 | 923 / 0.614 | 532 / 0.032 | 221 / 0.144 | 211 / 0.057 | 396 / 0.011 |
| 0.25 | 1,414 / 0.888 | 514 / 0.971 | 744 / 0.604 | 532 / 0.032 | 218 / 0.147 | 211 / 0.057 | 396 / 0.011 |
| 0.4 | 1,112 / 0.880 | 372 / 0.976 | 598 / 0.615 | 468 / 0.037 | 171 / 0.217 | 171 / 0.084 | 375 / 0.013 |
| 0.5 | 911 / 0.886 | 269 / 0.975 | 480 / 0.675 | 431 / 0.042 | 152 / 0.278 | 158 / 0.089 | 373 / 0.013 |
| 0.6 | 728 / 0.928 | 189 / 0.979 | 383 / 0.737 | 399 / 0.047 | 127 / 0.349 | 135 / 0.120 | 353 / 0.015 |
| 0.75 | 472 / 0.923 | 113 / 0.966 | 206 / 0.891 | 329 / 0.075 | 85 / 0.497 | 93 / 0.163 | 288 / 0.019 |
| 0.82 (default) | 343 / 0.921 | 85 / 0.966 | 144 / 0.861 | 275 / 0.121 | 59 / 0.557 | 74 / 0.216 | 234 / 0.021 |
Single-exon models at 0 / 0.82: mouse ES 15,900 / 2,218; heart 1,470 / 106; OR6737 3,694 / 342; A119b 37,272 / 10,848.
Smallest `reads` of a single-exon model per quantile (0 ... 0.82): OR6737 2, 2, 2, 3, 3, 4, 6, 8; A119b 2, 2, 2, 2, 3, 4, 7, 10. On
the other datasets that attribute stays at 2 while the models change, so it does not show the floor there (not examined).

**Choice by the registered rule: none.**
- q* (pays against quantile 0 on >= 4 of 5 N): 0.1 pays on fly and zebrafish only (on mouse ES and heart it changes nothing, on
  mouse simulated precision is 0.0005 lower); 0.25 pays on heart only; every higher value loses more than 5% everywhere.
- q_b (keeps >= 80% on >= 4 of 5 N, more precise than quantile 0 on both real N): 0.25 keeps 80% on 3 of 5 and is not more
  precise on mouse ES; 0.1 is not more precise on mouse ES or heart; higher values keep too little.
- Step 2 therefore has no value to check. The adrenal arms exist (`adrenal/q*/x.gtf`) and were not scored: adrenal stays unused
  for this lever.

**What the grid shows (my interpretation).**
- There is no knee. On the simulated sets precision hardly depends on the floor (fly 0.88-0.93, mouse 0.96-0.98) while the found
  truth grows 5-10 fold as it falls: the cost row 1831 measured there is the cost of a noise filter on data without noise.
- On every real library, mouse and ape alike, lowering the floor to 0 multiplies the found truth by 1.7-3.7 and divides
  precision by 2-4, with 3-14 times as many single-exon models. The floor trades the same way on mouse ES, mouse heart, OR6737
  and A119b: it is not fitted to primates.
- `final.gtf` does not use this floor: its single-exon models come from the read-end clusters (row 1831: mouse heart 128 found
  against 59 in `PREFIX.gtf`, A119b 315 against 234).
Nothing proposed. A lower floor is a recall choice for the user (0.6: heart 127 found at 0.35, mouse ES 399 at 0.047, A119b 353
with 23,428 models against 10,848); adrenal is the dataset to check it on.
