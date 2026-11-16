# PREREG: three assembly opt-ins on four unseen cell-line samples, with a promotion rule (2026-11-16, register 1760)

User (2026-11-16): "now that we have a lot of data are there any other opt ins that can be adopted as default for improving the
assembler?" ... "lets proceed with preregistration and run". Written before any rule output exists on these samples.

## Samples
Cell-line Iso-Seq of the T2T individuals (`docs/APE_SAMPLES_2026-11-15.md`): PPY AG05252 and PAB AG06213 (fibroblast), PPA PR00251
(fibroblast, fragmentary RNA), PTR AG18354 (lymphoblastoid). Default assemblies in workspace `cellline_asm/STEM/`. No rule has been
fitted or scored on them. Already seen (descriptive): default chains / correct / precision, precision by read stratum, the share of
chains that are a truncated piece of a RefSeq chain, read length, unspliced share. The same species' testis samples are spent.

## Tests, each sample on its own, bars copied unchanged from the rows that set them
- **A. `shiftdrop` (rows 1739 / 1747):** dropped chains in RefSeq <= 0.05; their moved intron annotated in RefSeq or CAT <= 0.15;
  correct chains lost <= 0.2%; precision of the filtered GTF not below the equal-size read cut. All four = PASS.
- **B. frozen logistic score (rows 1677 / 1747), fitted on OR6737 + A119b, not refitted:** score curve minus read-count curve at 0.90
  of the default correct chains >= +0.02.
- **C. `RUSTLE_EXTEND_TRANSFRAGS=1` (rows 1616-1618):** precision of the chains it adds (in the extended assembly, not in the default)
  >= the precision of the default's 2-read chains on the same sample. Correct chains lost are reported beside it.
- **Reported, no promotion question:** `trunc5` against the read-count curve (bar +0.02; row 1743 showed it drops CAGE-supported
  starts, so it is a tag and not a default filter whatever it scores here); the three r1677 chain classes.

## Promotion rule, stated now
A rule is a **default candidate** if it passes on at least three of the four samples AND on both orangutans. The bonobo cannot decide:
its fragmentary RNA favours any rule that removes truncated chains, so a result that needs the bonobo to reach three is not a candidate.
A candidate is put to the user with its history on the earlier samples (shiftdrop 5 of 6; score 4 of 5, the failure being fibroblast
KB3781; extension 1 of 2); switching a default stays the user's decision and is not made by this row. A rule that is not a candidate
stays opt-in and is not retested on these samples.

## Limits stated now
- RefSeq label only for A-C (CAT enters bar 2 of A); absolute precision is not compared across species.
- B is a score, not a stage: a candidate verdict means "worth building as a filter", with its cut still to be fixed on spent samples.
- Kinnex libraries are 3 to 8 times deeper than the testis samples; high read strata will be larger than in any earlier test.

Workspace `r1760_cellline_rules/` (`job.sbatch`, one job per sample).

## RESULT (register 1760; jobs 4473784-4473787; workspace `r1760_cellline_rules/job.*.out`, `frozen.*.txt`, `rules.*.txt`)
| test | PPY fib | PAB fib | PPA fib | PTR lcl | verdict |
|---|---|---|---|---|---|
| A `shiftdrop`: dropped in RefSeq <= 0.05 | 0.036 (277) | **0.061 (296)** | 0.031 (229) | 0.044 (319) | PASS / FAIL / PASS / PASS |
| A moved intron annotated <= 0.15 | 0.094 | 0.118 | 0.100 | 0.103 | PASS all |
| A correct lost <= 0.2% | 0.04% | 0.07% | 0.03% | 0.05% | PASS all |
| A precision vs equal-size read cut | 0.3040 vs 0.3034 | 0.3490 vs 0.3485 | 0.2132 vs 0.2130 | 0.3420 vs 0.3416 | PASS all |
| B frozen score at 0.90 (bar +0.02) | +0.072 | +0.023 | +0.096 | +0.026 | PASS all |
| C new chains vs default 2-read chains | 4,951 at 0.182 vs 0.127 | 3,054 at 0.199 vs 0.152 | 8,888 at 0.157 vs 0.089 | 4,835 at 0.181 vs 0.163 | PASS all |
| C net correct (correct lost) / precision | +875 (28) / 0.303 -> 0.302 | +580 (28) / 0.348 -> 0.346 | +1,352 (41) / 0.213 -> 0.217 | +804 (73) / 0.341 -> 0.337 | reported |
| `trunc5` vs read curve (bar +0.02, reported) | +0.027 | +0.018 | +0.032 | +0.014 | 2 of 4 |
**Verdict by the registered rule.** B (frozen score) and C (`RUSTLE_EXTEND_TRANSFRAGS=1`) are default candidates: four of four, both
orangutans. A (`shiftdrop`) is not: PAB misses bar 1 (18 of 296 dropped chains are in RefSeq), so it stays opt-in and is not retested
here. No default was switched by this row. The four samples are spent for these rules. Read with the verdict: C raises recall and
lowers whole-assembly precision slightly on three of four samples (its new chains are better than our weakest accepted chains, not
better than the average chain); B is a ranking, and a filter built on it still needs a cut fixed on spent samples.
