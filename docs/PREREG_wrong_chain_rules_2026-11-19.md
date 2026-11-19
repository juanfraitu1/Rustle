# Row 1814: single-evidence rules against wrong chains, found on two non-primate samples (pre-registered 2026-11-19)

User (2026-11-19): pick two samples for new levers that are not overfit; take mouse ES and zebrafish simulated and see what can
remove wrong calls, on `PREFIX.gtf` only. Registered before any feature table of these runs was looked at for this purpose.
Nothing here changes a default.

**Development data.** `PREFIX.gtf` at today's defaults (row 1812 = arm `ext0_seed` of row 1811) of mouse ES (ENCSR172GXL; mouse
contigs against GENCODE vM27, SIRV contigs against the SIRV annotation) and zebrafish simulated. No primate sample.
**Evidence per chain.** The read features the chain score already computes (`bench/famcmp/chain_rules_feat.py`, tables of row
1813: `feat.es.ext0_seed.tsv`, `feat.dre.ext0_seed.tsv`), on their raw scale, plus the two genome extras. No annotation.

**Search, fixed now.** One rule = one feature, one direction, one threshold ("drop a chain when feature > t" or "< t"). For
every feature and direction, t runs over the pooled percentiles 1..99 of the two samples. A rule is a **candidate** when, on both
development samples (mouse contigs for ES), it loses at most 1% of the correct chains and removes at least 5% of the wrong
chains; per feature and direction the candidate removing the most wrong chains (summed over the two) is kept. The read count
itself is left out (that is the read floor). SIRV shown per candidate, not used for the choice.

**Check, bars per candidate, threshold frozen.** Datasets: fly simulated, mouse simulated, mouse heart, OR6737, A119b (arm
`ext0_seed`, tables of row 1813; used before for other rules, reused by the user's ruling) and mouse adrenal (ENCSR852HRM,
assembled now at today's defaults; never used to choose a rule).
- **W1.** On mouse adrenal: at most 1% of the correct chains lost, and the share of correct chains among the removed is at most
  half the precision of `PREFIX.gtf`.
- **W2.** The same two conditions on at least 4 of the other 5 datasets.
A candidate passing W1 and W2 TRANSFERS and is put to the user as an opt-in post-filter; correct and wrong removed are reported
separately for every dataset, and the overlap between passing rules.

**Limits.** GENCODE is incomplete: "wrong" on real mouse includes real chains, so a rule can look better than it is there; the
zebrafish and SIRV sides have complete truth. A post-filter only removes; inside the assembler the same rule could change later
steps. Thresholds on raw features may depend on depth. One threshold per feature: no combinations are searched.

## Result (2026-11-19): NO CANDIDATE

Development tables: mouse ES 54,530 chains on the mouse contigs (18,059 correct), zebrafish simulated 15,470 (12,578 correct).
No feature in either direction removes at least 5% of the wrong chains on both samples within 1% of the correct chains
(`search.txt`: 0 candidates). The check step did not run: there is nothing to check; the mouse adrenal table was made
(27,114 chains) and not opened.

Shown without a bar (`diag.txt`, chosen after the result): per feature the threshold with the largest wrong share removed on the
weaker of the two samples, at most 1% of correct chains lost on both (correct lost / wrong removed, share of the wrong):

| rule | mouse ES | zebrafish simulated |
|---|---|---|
| reads running past the chain's 5' end > 4 x the reads matching it (`ratio_cont5`) | 146 / 1,659 (4.5%) | 26 / 205 (7.1%) |
| weakest intron's share of its locus < 0.01 (`min_isup_rel`) | 174 / 3,109 (8.5%) | 100 / 79 (2.7%) |
| reads running past the 5' end > 122 (`cont5`) | 169 / 647 (1.8%) | 85 / 490 (16.9%) |
| shortest intron < 62 bp | 3 / 449 (1.2%) | 119 / 124 (4.3%) |
| mismatch rate > 0.03 (`de`) | 120 / 1,929 (5.3%) | 18 / 30 (1.0%) |

**Reading.** At a cost of 1% of the correct chains no single read feature removes much: the wrong chains that remain at today's
defaults do not differ from correct ones in any one feature. What the two samples lack differs (zebrafish: pieces of longer
chains; mouse ES: weakly supported introns and mismatched reads), so a rule good on one is weak on the other.
