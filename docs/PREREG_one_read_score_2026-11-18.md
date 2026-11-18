# Row 1806: a score for one-read chains from the reads themselves (pre-registered 2026-11-18)

User (2026-11-18): build the score for one-read chains from read-level features. Follows row 1805 (evidence borrowed from
better-supported chains does not work). Registered before any run below. Datasets are reused (user ruling 2026-11-18); a rule
that comes out of this is checked on one unseen library in a later row. Nothing here changes a default.

**Runs.** Ours at today's defaults (polish junction snap on) with `RUSTLE_MIN_CHAIN_READS=1`, `RUSTLE_PTC=0`, on seven datasets:
fly, mouse and zebrafish simulated; mouse ES, mouse heart, OR6737, A119b. `PREFIX.gtf` only.
**One-read chain** = multi-exon chain with `reads "1"`. Label = exact intron chain in the truth of rows 1793 / 1800.

**Features** (existing extractor `bench/famcmp/chain_rules_feat.py`, all computed from the BAM and the genome, no annotation):
the 23 chain score features plus non-GT-AG introns and shortest intron, design as r1769 (within-sample standardisation over the
one-read chains, decile indicators). For a one-read chain these are: edit rate of the read (`de`), mapping quality, wobble around
the splice sites, poly-A signal and A-richness at the 3' end, intron retention, reads carrying a contiguous piece of the chain
(`partial`, `cont5`, `cont3`), support of its weakest intron, share of the locus, chain length and intron number.
Secondary model, shown: the same without the locus features (`locus_frac`, `locus_chains`, `min_isup`, `min_isup_rel`).

**Evaluation: leave one dataset out.** Weights fitted on the one-read chains of the other six (pooled logistic, `chain_rules.fit`),
applied to the held-out dataset.
**Bar, per held-out dataset.** Take the one-read chains in score order down to the last point where the accepted set is at least
as precise as that dataset's chains with two or more reads (so adding them does not lower precision). PASS = that set holds at
least 30% of the dataset's correct one-read chains. The cut is read on the held-out labels: this measures whether the weights
separate, not yet a usable threshold.
**Verdict.** SEPARABLE = PASS on at least 5 of 7, at least 2 of the 4 real libraries among them. Otherwise named per dataset.
Reported besides: the share kept at that point for the secondary model and for the read-count-free baselines (random order =
the class precision); a threshold carried over from the training datasets (largest of their own precision-matching cuts) and
what it accepts on the held-out one; the five largest weights of the all-data fit.

**Limits.** RefSeq / GENCODE are incomplete: a real lone one-read chain cannot be told from a wrong one on real libraries, so
labels there are noisy against us. Simulated reads lack the artefacts that make real one-read chains wrong (row 1805). A PASS
needs: a threshold rule, an implementation, the check on an unseen library.

## Result (2026-11-18): NOT SEPARABLE as registered; PASS on 1 of 7 (fly simulated)

Jobs 4482635-4482641, all exit 0, snap on (`x.assemble.polish` checked). Tables: workspace `r1806_one_read_score/lodo.txt`,
`diag.txt`, features `{dme,mmu,dre,es,heart,or,a1}.tsv`.

| held out | target (>= 2-read chains) | one-read chains / correct | accepted at target / correct | share kept | |
|---|---|---|---|---|---|
| fly sim | 0.855 | 547 / 371 | 306 / 262 | 0.706 | PASS |
| mouse sim | 0.895 | 2,520 / 1,801 | 48 / 43 | 0.024 | fail |
| zebrafish sim | 0.795 | 1,698 / 869 | 1 / 1 | 0.001 | fail |
| mouse ES | 0.326 | 19,789 / 1,005 | 0 | 0 | fail |
| mouse heart | 0.557 | 14,922 / 1,707 | 0 | 0 | fail |
| OR6737 | 0.351 | 55,259 / 3,262 | 0 | 0 | fail |
| A119b | 0.172 | 128,091 / 3,327 | 0 | 0 | fail |

Without the locus features: fly PASS (0.976), every other dataset at most 0.002. Carried-over cut: nothing usable (fly 290
accepted at 0.862; A119b 10,531 at 0.060; the others under 60 chains or far below target).
Largest contributions in the all-data fit: `r_partial`, `l_partial` (other reads carrying a piece of the chain), `pas`,
`retention`, `r_cont5`.

**Diagnostic, not registered** (`diag.txt`): is it the transfer of the weights, or the features? Precision of the best-scored 5%
of one-read chains, weights from the other six / weights fitted inside the dataset (two-fold):

| dataset | class precision | from the other six | inside the dataset | target |
|---|---|---|---|---|
| fly sim | 0.678 | 0.741 | 0.963 | 0.855 |
| mouse sim | 0.715 | 0.873 | 0.849 | 0.895 |
| zebrafish sim | 0.512 | 0.583 | 0.845 | 0.795 |
| mouse ES | 0.051 | 0.110 | 0.139 | 0.326 |
| mouse heart | 0.114 | 0.105 | 0.312 | 0.557 |
| OR6737 | 0.059 | 0.156 | 0.204 | 0.351 |
| A119b | 0.026 | 0.068 | 0.105 | 0.172 |

AUC 0.59-0.81 carried over, 0.71-0.93 inside. **Reading.** The features rank one-read chains better than chance everywhere, and
on simulated reads a score fitted in the dataset reaches the target. On the four real libraries even the best 5% under weights
fitted on the same library stay at about half the precision of the chains with two reads: with these features and an
annotation as truth, a one-read chain cannot be brought to the level of the rest of the output. Part of that gap is unannotated
truth, which this design cannot size.

**Follow-up (user, 2026-11-18): a separate ranked tier.** The opt-in driver stage `oneread` (register 1766, `PREFIX.oneread.gtf`
= the strict coding-frame tier) now also writes `PREFIX.oneread.ranked.gtf`: every one-read chain the default lacks, best
`oneread_score` first, with `oneread_rank` (`bench/famcmp/oneread_tier.py apply`, model `tools/data/oneread_r1806.json` = the
pooled fit on the seven feature tables of this row, 222,850 one-read chains). A ranking, not a filter; never appended to
`final.gtf`; `PREFIX.gtf` untouched. Test run on mouse heart (job 4482658, exit 0; in the fit, so not a transfer figure): 14,867
ranked chains, none in the default; GENCODE precision by rank: best 5% 0.301, next 15% 0.202, next 30% 0.128, lower half 0.060
(default chains 0.557); the strict tier of the same run holds 686 chains, 181 correct (0.264).

**Check of the ranked tier on a library outside the fit (registered 2026-11-18, before the run; user: one more of the tissues
we already have).** Chimpanzee testis (`reads/PTR/PTR_mm.bam`, truth RefSeq `cellline_asm/PTR.refseq.gtf`); not among the seven
fitted datasets (it was the development sample of the strict tier's rules, row 1766, which the ranking does not use).
`assemble` at today's defaults, then stage `oneread`. Rank slices: best 5%, next 15%, next 30%, lower half.
**The ranking TRANSFERS** when precision falls from each slice to the next and the best 5% is at least twice the precision of
all ranked chains. Shown: the strict tier of the same run. (Correction: earlier text named an "ONT" check; that was a misread
typo, there are no ONT reads here.)

**Result of the check (2026-11-18): the ranking TRANSFERS.** Job 4482667, both stages exit 0. Chimpanzee testis: default 56,754
chains, 21,897 correct, 0.386. Ranked tier 52,880 chains (none in the default), 4,950 correct, 0.094. By rank: best 5% 2,644 /
638 / 0.241 (2.6 times the tier); next 15% 7,932 / 1,287 / 0.162; next 30% 15,864 / 1,709 / 0.108; lower half 26,440 / 1,316 /
0.050. Strict tier of the same run (not independent here): 1,414 chains, 522 correct, 0.369.
The best slice stays below the default's 0.386, as row 1806 found: the tier is for browsing, not for the main output.

**Second check, a non-primate library outside the fit (registered 2026-11-18, before mapping; user: check a non-primate).**
ENCODE ENCSR862DMK, adult (2 months) B6CASTF1/J mouse left cerebral cortex, PacBio Sequel II, replicates ENCFF285NUG +
ENCFF977AFJ pooled (1,900,226 reads), truth GENCODE vM27. Same species as two fitted datasets (ES cells, heart), another tissue
and library; no non-primate library outside the fit was on disk. Shipped mapping, `assemble` at today's defaults, stage
`oneread`; the same slices and the same bar as the chimpanzee check.

**Result of the second check (2026-11-18): the ranking TRANSFERS.** Mapping 4482675 (1,900,226 reads in and out), job 4482676,
both stages exit 0. Mouse cortex: default 33,627 chains, 15,696 correct, 0.467. Ranked tier 20,778 chains (none in the
default), 1,796 correct, 0.086. By rank: best 5% 1,038 / 246 / 0.237 (2.8 times the tier); next 15% 3,117 / 479 / 0.154; next 30%
6,234 / 595 / 0.095; lower half 10,389 / 476 / 0.046. Strict tier of the same run (independent here): 935 chains, 207 correct, 0.221.
