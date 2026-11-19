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
