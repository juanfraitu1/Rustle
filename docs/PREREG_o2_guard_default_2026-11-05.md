# PREREG O2 PRIMARY-ON-FAMILY GUARD DEFAULT UNDER --tie-ms (register 1542, 2026-11-05, before any run)

Why (user 2026-11-05 "lets do that", the O2 ledger's one remaining next-item): r1514 measured the guard on KB3781 held-out
(decidable 129 right / 22 wrong vs the primary's 125 / 26; the 3,465 outside-catalog assignments of r1510 gone; panel all 1,275
kept plus 159) and FAILED two side bars written for identity, not for the decidable class. This row does not re-measure the
truth-based claim (KB3781 is the only sample with haplotype truth and both of its splits are spent); it registers the product
change and bars that a second sample can check without truth.
Change: under `--tie-ms W` the primary-on-family guard is ON; `--no-primary-on-family` is the escape; `--primary-on-family` keeps
its meaning without `--tie-ms` (exact-AS gate plus guard, as before). The default path (no `--tie-ms`) does not change.
Bars (all must hold):
1. Default path byte-identical: the r1509 KB3781 panel base run (`o2mstie/base.assignments.tsv`, default gate, no flags) rerun
   with the new binary -> every product file identical.
2. Adoption byte-identical: the r1514 run (`--tie-ms 10 --primary-on-family`, 494 families) rerun as `--tie-ms 10` alone with
   the new binary -> `assignments.tsv` identical to `pof10.assignments.tsv`; `--tie-ms 10 --no-primary-on-family` reproduces
   r1510's `all10.assignments.tsv` (the escape restores the measured pre-guard behaviour).
3. Second sample, no truth (human A119b, CHM13, catalog `~/runs/human_A119b/a1.fam.copies.tsv`): the 40 families with the most
   contested rows among those with 50-2,000 (`~/r1540/human_panel.txt`, avoids the chr13 hotspot), `--tie-ms 10`: (a) every
   assigned molecule has a primary record within its family's windows (the guard's defining property; a violation is a bug);
   (b) the share of `--no-primary-on-family` assignments removed by the guard is reported (r1514 held-out: 77%; no bar, no
   truth); (c) wall <= 1 h each arm on 64 threads.
Pass all -> the guard is part of `--tie-ms` (docs: definition, ledger, THESIS_OBJECTIVES); `--tie-ms` itself stays opt-in (the
r1510 held-out verdict stands: O2 = primary on correctness, the gain is not assigning reads the catalog cannot explain).
Any bar fails -> the change is reverted, the row recorded.

## RESULT r1542 (2026-11-05): PASS, adopted
Job 4442213 (`/anvme/workspace/v138da13-rustle/r1542`, `~/r1540/r1542.sbatch`), binary ~/rustle_target_test/release (commit
"copy_assign: primary-on-family guard on under --tie-ms").
1. Default path: KB3781 panel base rerun identical on assignments, families, quant, famcn_readonly, family_join (the job's
   "copy_conversions DIFFER" line is the script comparing a file neither run writes). PASS.
2. `--tie-ms 10` alone = r1514 `pof10.assignments.tsv` byte-identical (2,543 assigned, 310 s); `--tie-ms 10
   --no-primary-on-family` = r1510 `all10.assignments.tsv` byte-identical (5,804, 447 s). PASS.
3. Human A119b, 40 families (38 regions): guarded 5,800 assigned rows (5,702 molecules), escape 6,355 (6,246); removed share
   8.7% (KB3781 held-out was 77%: the human catalog's parent genes are in the catalog). Wall 689 s / 983 s (bar 1 h). PASS.
   (a) 5,800 of 5,800 assigned rows have a primary record in the family's read windows as the binary builds them (every copy of
   every family bound to the same swept region, +-50 kb). Under a narrower reading, the molecule's OWN family's copies +-50 kb,
   it is 5,719 of 5,800: 81 rows (1.4%) have their primary on a sibling family's copy of the same region, the class r1514
   reported as "of another family (4)". The bar is stated on the binary's windows; the narrower count is reported, not judged.
Adopted: under `--tie-ms` the guard is on; `--no-primary-on-family` restores r1510's gate; `--tie-ms` stays opt-in.
