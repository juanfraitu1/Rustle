# Row 1847: the two parked SIRV rules, tested on the twelve datasets (pre-registered 2026-11-20)

User (2026-11-20): the rules parked for the SIRVs were never run; see whether they help. This row runs the draft of row 1801
(`docs/PREREG_sirv_rules_2026-11-18.md`, bars never accepted, nothing run). Of its four rules two are settled since: every
read counted (rows 1833 / 1836, opt-in `RUSTLE_DUP_COUNT`) and extension off (default since row 1811). Left, both existing
opt-in switches, untested off the spike-ins:
- **S** `RUSTLE_POLISH_SPLICED_SHARE=0.05`: drop a transcript with an intron spliced by under 5% of the primary reads across it.
- **M** `RUSTLE_POLISH_MONO_SHADOW_RATIO=1`: the shadow rule keeps a single-exon transcript that has at least as many reads as
  the spliced transcript shadowing it.
Both values are the draft's, chosen on the SIRV contigs (development data). Workspace `r1847_sirv_rules`. Committed before any run.

## Arms and data
`base` (today's defaults: junction ratio 20, no end layer), `S`, `M`, `SM` on the twelve datasets of rows 1837-1846. Stage
`assemble`, `RUSTLE_PTC=0`, `RUSTLE_ONEREAD_TIER=0`, binaries `~/rustle_target` (`3e26f1e7`). Scorer of row 1831 (exact chains;
single-exon models matched at 50% of both lengths). SIRV contigs of the mouse ES library shown apart (development data, no bar).
The draft's check datasets (zebrafish, heart) are no longer unused; all twelve are reused, said (user ruling on reuse).

## Bars (the draft's D1 and D2, per dataset, on `x.gtf`; `x.final.gtf` shown)
- **D1, S against base:** chain precision not lower, and at most 0.5% of base's correct chains lost.
- **D2, M against base:** the multi-exon chains are the same set, single-exon model precision is not lower, and no fewer truth
  single-exon transcripts are found.
- A rule is a **candidate** when its bar holds on at least 6 of 7 non-primate and 4 of 5 primate datasets. `SM` is shown.
  A candidate goes to the user; no default changes here.

## Expectation
S removes few chains (on mouse ES the draft counted 1,015 dropped, 45 of them annotated: 0.25% of the correct) and should hold
D1 on real libraries; on simulated reads it can only cost, since no simulated read carries an artefact intron. M adds
single-exon models where a weak spliced artefact shadows a strong unspliced gene; I am unsure of D2 on the ape testis
libraries, where single-exon precision is 0.02-0.25. On SIRV: S removes the introns in unspliced spike-ins (4 or 5 of the 19
wrong in `final.gtf`), M restores up to 4 true single-exon spike-ins.

## Limits
- S re-reads the BAM per contig (wall time reported). Real-library truths are incomplete; an unannotated real intron counts
  as wrong and an unannotated single-exon gene as unmatched.
- The chain score filter was fitted without either rule.

## Result (2026-11-21): neither rule is a candidate; the spliced share misses by the two simulated sets, the shadow support fails everywhere

Jobs 4490962-4490973, 48 runs exit 0; `r1847_sirv_rules/score.NAME.txt`. The S arm costs 2-9% more wall time (A119b 783 -> 854 s).

**D1, spliced share 0.05 (`x.gtf`): correct chains lost, wrong chains removed, precision.**
| | fly sim | mouse sim | zebrafish sim | ES | heart | adrenal | cortex | OR6737 | A119b | chimp | bonobo | KB3781 |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| correct lost | 47 (0.67%) | 185 (0.67%) | 46 (0.37%) | 41 (0.23%) | 10 (0.08%) | 17 (0.12%) | 29 (0.18%) | 57 (0.22%) | 190 (0.49%) | 29 (0.14%) | 16 (0.09%) | 111 (0.43%) |
| wrong removed | 8 | 34 | 11 | 755 | 94 | 145 | 166 | 505 | 3,767 | 179 | 263 | 1,032 |
| precision | 0.9276 -> 0.9281 | 0.9065 -> 0.9070 | 0.8159 -> 0.8159 | 0.3485 -> 0.3531 | 0.5571 -> 0.5591 | 0.5511 -> 0.5539 | 0.4764 -> 0.4783 | 0.3254 -> 0.3270 | 0.1591 -> 0.1609 | 0.3535 -> 0.3542 | 0.3107 -> 0.3119 | 0.3393 -> 0.3430 |
| D1 | FAIL | FAIL | pass | pass | pass | pass | pass | pass | pass | pass | pass | pass |
**5 of 7 non-primate, 5 of 5 primate; bar 6 and 4: not a candidate.** Precision is not lower anywhere. The two failures are
fly and mouse simulated, where it removes five true chains per wrong one (0.67% of the correct), as expected: a simulated read
carries no artefact intron, so the rule only meets true minor introns there. On the nine real libraries it removes 7-20 wrong
chains per annotated one (mouse heart 9, ES 18, A119b 20) and under 0.5% of the correct. `final.gtf` moves the same way.

**D2, shadow support ratio 1 (`x.gtf`):** multi-exon chains identical on all twelve. Single-exon truth found rises everywhere
(fly 343 -> 399, mouse ES 275 -> 335, cortex 91 -> 127, A119b 234 -> 325) but model precision falls on all twelve: slightly
on the simulated sets (0.921 -> 0.912, 0.966 -> 0.964, 0.861 -> 0.805) and by half or more on real libraries (heart 0.557 ->
0.278, cortex 0.596 -> 0.280, OR6737 0.216 -> 0.079, A119b 0.021 -> 0.009 with 10,848 -> 35,795 models; A119b base precision
0.245 -> 0.145). **0 of 12: FAIL.**

**SIRV contigs (`final.gtf`, development data, no bar):**
| arm | multi-exon: reported / correct / wrong | shifted site | other new intron | combination | single-exon: models / matching, of 23 |
|---|---|---|---|---|---|
| base | 67 / 48 / 19 | 13 | 5 | 1 | 11 / 11 |
| S | 62 / 49 / 13 | 12 | 0 | 1 | 15 / 15 |
| M | 67 / 48 / 19 | 13 | 5 | 1 | 19 / 17 |
| S + M | 62 / 49 / 13 | 12 | 0 | 1 | 19 / 17 |
S alone does both jobs on the spike-ins: it removes all five introns reported in unspliced spike-ins and one shifted site,
gains a correct chain, and by removing the spliced artefact it lifts the shadow from four true single-exon spike-ins
(precision of the multi-exon chains 0.716 -> 0.790). M adds two more true single-exon spike-ins and two unmatched models.
Expectation: S as written (the single-exon gain through the shadow was not foreseen); M's doubt about the apes was right and
understated.

**Decision:** no candidate by the registered bars. S stays the opt-in it is (`RUSTLE_POLISH_SPLICED_SHARE=0.05`): harmless or
better on every real library, a 0.67% loss of true chains on clean simulated reads. Whether that is acceptable for a default
is the user's call. M (`RUSTLE_POLISH_MONO_SHADOW_RATIO`) should not be used as it stands.

## Row 1848 (2026-11-21, user decision): the spliced share is the driver default
User, after the result above: "yes, lets make it default". `tools/rustle_pipeline.sh`: `RUSTLE_POLISH_SPLICED_SHARE` unset =
0.05; `RUSTLE_POLISH_SPLICED_SHARE=0` = the command before. Decided against the registered miss (D1 fails on fly and mouse
simulated, 0.67% of the true chains lost) because every real library gains. Check (workspace `r1848_spliced_share_default`,
fly): the new default's `x.gtf` and `x.final.gtf` are byte-identical to row 1847's `S` arm, and `=0` to its `base`. Test
suite 1,585 passed, 0 failed. Today's `final.gtf` per dataset = the `S` arm of row 1847 (`r1847_sirv_rules/NAME/S`).
