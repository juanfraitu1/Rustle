# PREREG O2 FAMILY CHECK = SOLE-CANDIDATE ASSIGNMENTS UNDER --tie-ms (register 1656, 2026-11-12; first written as 1654, renumbered: another session took 1654-1655 the same hour)

Why: r1651-r1653: a few families carry O2's wrong near-tie assignments. Step 0 on the DISCOVERY panels (blind
per-family features vs truth; `bench/o2_contested/o2_row_classes.py`): the families are those whose assigned rows have
no deciding site. Row class = `sole_candidate == 1` (contested molecule, one candidate copy in the family, assigned with
n_decisive 0, p 1): its rival sits outside the family's candidates and is never examined (register 787 said this of the
default gate). Discovery, decidable assigned rows, right / wrong:
| | sole candidate | other rows |
|---|---|---|
| HG002 on CHM13 panel (r1651) | 0 / 121 (all MCL189) | 170 / 54 |
| HG002 on its own haplotype panel (r1644) | 29 / 37 (all MCL338) | 451 / 8 |
| KB3781 (r1514 runs) | 33 / 0 | 489 / 48 |
Share of all assigned rows: 22%, 47%, 3%. So "the family does not fit" is one known row class, and it is not always
wrong (KB3781).

Test, no new run: the held-out sets of r1653 (`r1653_o2_guard/{hgpat,hg13}`, 60 families each, truth already made; their
sole-candidate split has not been looked at). Arm = treat `sole_candidate == 1` rows as abstentions, which is what the
existing `--no-sole-candidate` flag does.
Bars:
0. Power: >= 20 decidable sole-candidate rows pooled, else UNDERPOWERED.
1. On each set with >= 10 such rows, their wrong rate is >= 3 x the wrong rate of the other assigned rows.
2. On each set the wrong rate of assigned decidable rows falls when they are dropped.
Readout: correct rows lost; families; share of assigned rows that are sole candidates on the held-out sets and on the
A119b default-gate base run.
Verdict: 0-2 pass = PASS: recommend `--no-sole-candidate` whenever `--tie-ms` is used (existing flag, no code; making
it the default under `--tie-ms` is the user's call). Else FAIL: MCL189 / MCL338 were two families.
Expected: PASS on the CHM13 set; uncertain on the own-haplotype set (discovery 29 / 37).

Outputs: `/anvme/workspace/v138da13-rustle/r1656_o2_sole/`.

## RESULT r1656 (2026-11-12, no new run; `r1656_o2_sole/{hgpat,hg13,a1}.txt`)
Held-out sets, decidable assigned rows, right / wrong:
| | sole candidate | other rows | wrong rate, all -> without sole | correct lost |
|---|---|---|---|---|
| HG002 on its own haplotype | 19 / 15 (44% wrong) | 953 / 48 (4.8%) | 6.1% -> 4.8% | 19 of 972 (2.0%) |
| HG002 on CHM13 | 13 / 71 (85% wrong) | 997 / 70 (6.6%) | 12.3% -> 6.6% | 13 of 1,010 (1.3%) |
- Bar 0 PASS (118 rows). Bar 1 PASS: 9.2 x and 12.8 x the other rows' wrong rate. Bar 2 PASS on both.
**PASS as registered.** Sole-candidate rows are 4.3% (423 of 9,935) and 19.9% (1,226 of 6,159) of assigned rows on the
held-out sets, and 13.0% (1,595 of 12,294) of the A119b default-gate base run. Families: MCL1256 (0 / 46), MCL338
(13 / 24) on CHM13; spread over MCL253, MCL141, MCL70 on the own haplotype.
Against r1653's guard on the same sets: dropping sole candidates removes 86 wrong / 32 right; the never-off-the-primary
guard removed 103 / 39 and would delete 83% of the A119b default output; this class is 13% of it.
Recommendation (not applied): run `--tie-ms` with `--no-sole-candidate` (existing flag, rows become `tied`). Making
that the default under `--tie-ms`, or everywhere, is the user's decision: register 787's ruling was already "a tied
read is never assigned on one candidate", while the shipped default still has `sole_candidate` on
(`docs/copy_assignment_definition.md` line 26). That gap needs the user's ruling, and a default change needs the
A119b / KB3781 identity and count rows.
Limit: on KB3781 (discovery) the class was 33 right / 0 wrong; it is not wrong everywhere.

## Follow-up checks (2026-11-12, readout, no bar; `r1656_o2_sole/kb.txt`)
1. KB3781 held-out split (`o2mstie/pof_heldout`, reads from the reference individual): sole-candidate rows 8 right /
   0 wrong (15 of 1,109 assigned rows); panel 25 / 0 (55 of 1,434). The class is clean on the reference individual and
   rare there (1-4% of assigned rows); other rows are 121 / 22 and 368 / 26.
2. Register 787 vs the shipped default: not a contradiction, a partial implementation. 09-06 (L3, ledger 6fj) ruled a
   certified sole candidate `assigned`; 787 (6gz) reversed it after 4,705 EIF3C reads leaked into NPIP, with the rule
   "candidates = ALL tied placements' loci, whichever family". What shipped is the guard `tie_outside_catalog`: a tied
   placement on NO family's unit forces `tied` (it covers `--tie-ms` rivals too). The sole-candidate promotion itself
   (`family.rs` `sole_candidate_promotes`) stayed on for molecules whose tied placements are all on catalog units.
3. What today's wrong sole rows are: the truth copy is in ANOTHER family on 46 of 71 (CHM13) and 12 of 15 (own
   haplotype): the tie partner is a unit of a different family, which the guard counts as inside and this family never
   scores as a candidate. That is the "whichever family" half of 787's rule, not shipped. The rest (25 and 3) have the
   truth copy in the same family with no chain there. Right rows: truth always in the row's family.

## Row 1695 (2026-11-13, registered before the run): sole-candidate rows under the DEFAULT gate, measurement only
The user's ruling of 2026-11-13 made sole-candidate rows `tied` under `--tie-ms`. The shipped pipeline runs the default
gate (no `--tie-ms`), where the rule is unchanged and row 1656 has no truth-scored figure. Measurement: the r1657 HG002
own-haplotype panel (same families, regions, truth, strict scorer at slack 0) without `--tie-ms`, arms D (default) and
DN (`--no-sole-candidate`); sole rows = D minus DN.
Reading rule: recommend extending the ruling to the default gate only when there are at least 20 decidable sole rows
and their wrong rate is at least twice that of the other assigned rows. Otherwise no recommendation. Limits: one
sample on its own haplotype (the friendly case); the truth table was built from the `--tie-ms 10` run's reads, so a
default-gate row outside it is undecidable. No default is changed by this row.
Outcome: **meets the reading rule as registered, on thin counts.** HG002 own-haplotype panel, default gate (job 4470903): 419 assigned rows with sole candidates (38 sole) vs 381 without. Decidable at slack 0: all rows 22 right / 4 wrong; without sole rows 5 / 0; so sole rows 17 right / 4 wrong (19% wrong, 21 decidable) and all 4 wrong rows of the default gate are sole rows. Both conditions hold (21 >= 20; 19% vs 0%), but the comparison class is 5 rows and the sample is the friendly case, so this is a weak recommendation to extend the ruling to the default gate; the cost on A119b would be up to 2,229 of 5,784 assigned rows. The default was the user's call: ruled "Extend to the default gate" the same day (see `LEDGER_O2.md`); A119b assigned rows 5,784 -> 3,555.
