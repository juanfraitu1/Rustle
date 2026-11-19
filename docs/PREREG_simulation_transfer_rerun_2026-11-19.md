# Row 1817: the transfer test of rows 1779 / 1780 again, at the defaults of 2026-11-19 (pre-registered 2026-11-19)

Rows 1779 / 1780 (2026-11-18): DOES NOT TRANSFER as registered: never dominated, but below StringTie's curve on all four datasets
and the trained filter failed B3 on all three large ones. Since then five defaults changed (rows 1800, 1810, 1812, 1816):
junction snap on, sub-chain collapse ratio 0.3, transfrag extension off, `assemble` seeds with primary alignments, filter =
chain score r1813 + 5' truncation rule. Each was judged on its own; the combination has not been put beside the other tools.
Registered before any of these runs started. Nothing here changes a default.

**Runs.** Ours only, `assemble` at today's defaults (`RUSTLE_PTC=0` as in every run of this series), on the same three BAMs:
fly simulated, mouse simulated, mouse ES (with its SIRV contigs). StringTie, IsoQuant and FLAIR: the outputs of rows 1779 / 1780,
unchanged (`apetools/{DME_sim_mm,MMU_sim_mm,MMU_ES_mm}`). Same script (`r1779_simulated_reads/score/bars.py`), same truth files.

**Bars, the same three, per dataset** (fly, mouse simulated, mouse ES, SIRV; B3 not on SIRV):
- **B1** no other de novo tool has both more correct chains and higher precision than ours `PREFIX.gtf`.
- **B2** for every point of StringTie's coverage curve ours has a point on its read-count curve with at least as many correct
  chains and at least the same precision.
- **B3** `PREFIX.final.gtf` against `PREFIX.gtf`: precision higher and at least 0.93 of the correct chains kept.
**Verdict.** TRANSFERS = B1 and B2 on all four and B3 on the three large ones; otherwise the failed bars are named, as before.
**Shown without a bar:** the same table for the earlier defaults beside it; B1 and B2 read with `PREFIX.final.gtf` as ours;
`PREFIX.final.gtf` end to end on zebrafish simulated, mouse heart, OR6737 and A119b (correct chains / precision), which so far
exist only as feature-table estimates.

**Limits, stated in advance.** These datasets were used to judge the five changes (rows 1793-1816), so this is no longer a test
on unseen data: a pass says the changed defaults repair what the first test found, not that the assembler transfers to a new
species. The limits of rows 1779 / 1780 hold (simulated reads flatter every tool; GENCODE incomplete; SIRV loci simple).
The other tools ran once at their defaults and were not tuned.

## Result (2026-11-19): still DOES NOT TRANSFER as registered; the filter bar now passes, the curve bar still fails

Jobs 4482940-4482946, all exit 0. Tables: workspace `r1817_transfer_rerun` (`dme.txt`, `mmu.txt`, `mmu_es.txt`, `sirv.txt`,
`final_curve.txt`, `extra.txt`).

| dataset | B1 | B2 (2026-11-18 -> now) | B3 (2026-11-18 -> now) |
|---|---|---|---|
| fly simulated | PASS | FAIL every cut -> FAIL every cut | FAIL (kept 0.842) -> **PASS** 0.926 -> 0.964, kept 0.994 |
| mouse simulated | PASS | FAIL every cut -> FAIL every cut | FAIL (kept 0.842) -> **PASS** 0.903 -> 0.953, kept 0.979 |
| mouse ES | PASS | FAIL at cov >= 3, 4, 6, 11 -> FAIL at 3, 6, 11 | FAIL (kept 0.916) -> **PASS** 0.331 -> 0.347, kept 0.982 |
| SIRV | PASS | FAIL every cut -> **PASS** | not judged (0.589 -> 0.632, kept 0.906) |

Correct chains / precision, de novo (the other tools' outputs are those of rows 1779 / 1780):

| arm | fly (truth 9,379) | mouse sim (43,423) | mouse ES (114,584) | SIRV (60) |
|---|---|---|---|---|
| ours `PREFIX.gtf`, 2026-11-18 | 7,058 / 0.833 | 27,175 / 0.821 | 17,714 / 0.320 | 47 / 0.398 |
| ours `PREFIX.gtf`, now | 7,082 / 0.926 | 27,731 / 0.903 | 18,058 / 0.331 | 53 / 0.589 |
| ours `PREFIX.final.gtf`, 2026-11-18 | 5,945 / 0.974 | 22,892 / 0.916 | 16,233 / 0.411 | 41 / 0.547 |
| ours `PREFIX.final.gtf`, now | 7,041 / 0.964 | 27,159 / 0.953 | 17,730 / 0.347 | 48 / 0.632 |
| StringTie | 6,549 / 0.963 | 23,705 / 0.959 | 16,401 / 0.368 | 35 / 0.814 |
| IsoQuant | 6,592 / 0.972 | 24,098 / 0.967 | 13,443 / 0.472 | 39 / 0.848 |
| FLAIR | 6,614 / 0.720 | 27,048 / 0.535 | 16,686 / 0.096 | 48 / 0.081 |

**The limit that remains.** `PREFIX.gtf` at equal chain counts is still less precise than StringTie on the simulated sets
(fly 0.926 against 0.963 at the default point) and on mouse ES at middle coverage. The filter no longer costs recall.

**Shown without a bar: B2 read with `PREFIX.final.gtf` as ours.** Fly: covers every point of StringTie's curve (7,041 / 0.964
against 6,549 / 0.963 at the default point). SIRV: covers every point. Mouse simulated: below at six of seven cuts, by 0.001 to
0.006 precision with 13 to 15% more correct chains at the low cuts (27,159 / 0.953 against 23,705 / 0.959). Mouse ES: below at
cov >= 3 and 6 (16,679 / 0.400 against 15,243 / 0.450).
**`PREFIX.gtf` -> `PREFIX.final.gtf` end to end on the other four:** zebrafish simulated 12,574 / 0.813 -> 12,381 / 0.891; mouse
heart 13,287 / 0.546 -> 13,142 / 0.577; OR6737 25,792 / 0.321 -> 25,635 / 0.377; A119b 39,049 / 0.158 -> 38,091 / 0.219.
