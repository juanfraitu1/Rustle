# PREREG row 1697 (2026-11-13): quantification parity, re-measured on today's default assembly

User go 2026-11-13 ("Quantification parity"). Row 1088 (2026-09-23) read CLOSE on the assembler of that date. The
default assembly has changed since (polish, strand vote, sub-chain drop). This row repeats row 1088's measure on today's
default output with the comparators that exist on Fritz. (This file was written while job 4470931 was queued or
running and before its output was read; a shell slip skipped the first write.)

## Objects, join, quantity (row 1088's, verbatim where the inputs exist)
- Ours: default `PREFIX.gtf`, `reads` per transcript. A119b = workspace `r1617/default/a.gtf`; gorilla OR6737 =
  `r1628/orrun/or.gtf`.
- Tool: StringTie 3.0.1 `cov`, the lab runs on the same BAMs (`benchmark_collapse/{A119b,GGO}.stringtie.gtf.gz`).
- Join: strand-aware exact multi-exon intron chain; values summed over transcripts sharing a chain. Spearman rho
  (primary), Pearson r of log1p, shared chains, and the same for our `TPM` vs StringTie `TPM`.
- Not repeated: isoseq `count_fl` and FLAIR counts (those tables are not on Fritz), `--keep-coordinate-duplicates`.

## Reading rule (written before the output was read)
Per species: rho >= 0.85 on >= 1,000 shared chains = meets row 1088's bar for the StringTie cell. Both species meet it =
"StringTie cell COMPARABLE"; otherwise "CLOSE, unchanged" when both are >= 0.75, else "WORSE". Row 1088's full verdict
needs the isoseq cell too and is not re-issued by this row. Reference values: 0.898 (gorilla), 0.832 (human).

## Outcome (2026-11-13, job 4470934): CLOSE, unchanged
A119b: 86,456 shared chains, Spearman 0.827 (row 1088: 0.832), Pearson log1p 0.887; OR6737: 43,150 shared chains, Spearman 0.895 (0.898), Pearson 0.934. TPM vs TPM: 0.827 / 0.895.
Gorilla meets the 0.85 bar, human does not; both are within 0.005 of row 1088, so today's default assembly quantifies as the September one did. The systematic row 1088 named (de-duplicated counts sit below the tools' at depth) was not re-measured: `--keep-coordinate-duplicates`, isoseq and FLAIR counts are not in this row.
Disclosed: the first run (job 4470931) read StringTie's `longcov` in place of `cov` (the attribute pattern matched inside the longer name) and gave 0.472 / 0.569; the pattern was fixed and the job rerun. The registered quantity is `cov`; the wrong figures were seen before the fix.
