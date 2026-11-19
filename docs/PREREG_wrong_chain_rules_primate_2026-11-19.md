# Row 1815: the single-evidence rule search of row 1814, with primate development samples (pre-registered 2026-11-19)

User (2026-11-19), after row 1814 found no candidate on mouse ES and zebrafish simulated: try it again with primate samples.
Registered before the search ran on these tables. Nothing here changes a default.

**Same as row 1814:** the evidence (read features of `chain_rules_feat.py` on their raw scale, two genome extras, no annotation,
read count left out), the search (one feature, one direction, one threshold over the pooled percentiles 1..99; candidate = at
most 1% of the correct chains lost and at least 5% of the wrong chains removed on both development samples; per feature and
direction the candidate removing the most wrong chains), the script (`r1814_wrong_chain_rules/rules.py`), `PREFIX.gtf` at today's
defaults.
**Development samples:** OR6737 (gorilla) and A119b (human), tables `feat.or.ext0_seed.tsv`, `feat.a1.ext0_seed.tsv` of row 1813.

**Check, thresholds frozen.** A dataset passes a rule when at most 1% of its correct chains are lost and the share of correct
chains among the removed is at most half the precision of its `PREFIX.gtf`.
- **P1, another primate.** Chimpanzee testis (`PTR_mm`), assembled now at today's defaults; never used to choose one of these rules
  (it is a fit sample of the chain score, which these rules do not use).
- **P2, beyond primates.** At least 4 of the 6 non-primate datasets: fly, mouse, zebrafish simulated, mouse ES, mouse heart
  (tables of row 1813) and mouse adrenal (table of row 1814, unopened so far).
Verdict per candidate: P1 and P2 = TRANSFERS; P1 only = PRIMATE RULE (said so, like the levers of rows 1793 / 1811); neither = no.
Correct and wrong removed are reported separately per dataset, and the overlap between passing rules.

**Limits.** As row 1814. RefSeq on the apes is thinner than GENCODE on mouse, so more real chains count as wrong there and a rule
looks better than it is. A rule found on the two samples every earlier threshold was tuned on is the pattern rows 1793 and 1811
call overfit; P2 is the guard.

## Result (2026-11-19): one candidate, and it TRANSFERS

Tables: workspace `r1815_wrong_chain_rules_primate` (`search.txt`, `rules.json`, `check.txt`); chimpanzee testis assembled by job
4482907 (59,700 chains).
**The candidate (the only one of 54 feature-direction pairs):** drop a chain when `(cont5 + 1) / (full + 1) > 3.22`, where `full`
= reads whose intron chain is exactly the chain and `cont5` = reads that carry the whole chain and splice on past its 5' end.
In words: more than about three reads continue upstream for every read that starts where the chain does, i.e. the chain looks
like a 5'-truncated copy of a longer molecule. Read evidence only.

| dataset | role | correct lost (share) | wrong removed (share of the wrong) | correct among the removed | precision of `PREFIX.gtf` | passes |
|---|---|---|---|---|---|---|
| OR6737 | development | 49 (0.2%) | 8,464 (15.5%) | | | |
| A119b | development | 368 (0.9%) | 46,787 (22.5%) | | | |
| chimpanzee testis | P1 | 34 (0.2%) | 5,720 (14.8%) | 0.006 | 0.351 | yes |
| fly simulated | P2 | 6 (0.1%) | 97 (17.0%) | 0.058 | 0.926 | yes |
| mouse simulated | P2 | 232 (0.8%) | 494 (16.5%) | 0.320 | 0.903 | yes |
| zebrafish simulated | P2 | 45 (0.4%) | 371 (12.8%) | 0.108 | 0.813 | yes |
| mouse ES | P2 | 194 (1.1%) | 1,989 (5.5%) | 0.089 | 0.331 | no (1.1% lost) |
| mouse heart | P2 | 62 (0.5%) | 665 (6.0%) | 0.085 | 0.546 | yes |
| mouse adrenal | P2 | 92 (0.6%) | 806 (6.5%) | 0.102 | 0.540 | yes |

P1 PASS, P2 PASS (5 of 6). SIRV (mouse ES): 3 correct and 5 wrong chains removed of 53 / 37.
Row 1814 had the same feature as its closest miss (threshold 4: 4.5% / 7.1% of the wrong chains on mouse ES / zebrafish), so the
rule is weaker on real mouse than on the apes but the same rule.
**Not measured:** how much of what it removes the chain score filter removes anyway; the rule inside the assembler; its relation
to the sub-chain collapse ratio (row 1810), which drops the same kind of chain by transcript counts instead of read counts.

## Follow-up (2026-11-19, user): implemented opt-in; the two overlaps (shown, no bar)

`RUSTLE_TRUNC5_RATIO=3.22` (unset = off): the chain score pass (`bench/famcmp/chainscore.py apply`, which already counts these
reads) also drops the rule's chains from `PREFIX.final.gtf`. Fly check: 751 chains dropped, 3 of them by the rule only.
Table `r1815_wrong_chain_rules_primate/overlap.txt` (`r1814_wrong_chain_rules/overlap.py`).

**A. On top of the chain score** (today's defaults; correct chains / precision of what is kept):

| dataset | r1769 alone | r1769 + rule | r1813 alone | r1813 + rule |
|---|---|---|---|---|
| A119b | 38,858 / 0.174 | 38,601 / 0.196 | 38,309 / 0.207 | 38,226 / 0.212 |
| OR6737 | 25,676 / 0.349 | 25,651 / 0.363 | 25,724 / 0.348 | 25,695 / 0.364 |
| chimpanzee testis | 20,926 / 0.373 | 20,904 / 0.391 | 20,896 / 0.383 | 20,880 / 0.397 |
| mouse ES | 16,162 / 0.424 | 16,139 / 0.427 | 17,852 / 0.340 | 17,792 / 0.344 |
| mouse heart | 12,937 / 0.586 | 12,923 / 0.589 | 13,211 / 0.557 | 13,189 / 0.565 |
| mouse adrenal | 14,163 / 0.584 | 14,154 / 0.587 | 14,534 / 0.552 | 14,509 / 0.560 |
| mouse simulated | 26,334 / 0.937 | 26,269 / 0.943 | 27,457 / 0.915 | 27,355 / 0.924 |
| zebrafish simulated | 11,443 / 0.913 | 11,439 / 0.913 | 12,490 / 0.834 | 12,471 / 0.845 |
| fly simulated | 6,658 / 0.964 | 6,657 / 0.965 | 7,068 / 0.938 | 7,066 / 0.944 |

The filter already removes half or more of the rule's chains. What is left is worth 0.014 to 0.022 precision on the three apes
for 22 to 257 correct chains (r1769), and 0.000 to 0.011 on the non-primate datasets.

**B. Against the sub-chain collapse ratio** (chains of the ratio-0.3 run; correct / wrong chains dropped):

| dataset | ratio 0.7 | the rule | both |
|---|---|---|---|
| A119b | 728 / 41,823 | 372 / 44,419 | 206 / 28,925 |
| OR6737 | 184 / 9,109 | 52 / 8,326 | 40 / 6,180 |
| mouse ES | 428 / 2,566 | 187 / 2,013 | 154 / 1,517 |
| mouse simulated | 647 / 993 | 231 / 518 | 209 / 484 |
| fly simulated | 37 / 198 | 5 / 98 | 4 / 96 |

The two drop largely the same chains. The rule loses a third to a half of the correct chains ratio 0.7 lost, for 0.5 to 1.06 of
the wrong chains: on the apes it gives back the precision that lowering the ratio to 0.3 cost (row 1810) at a smaller price.
