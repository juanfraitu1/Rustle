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
