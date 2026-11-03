# PREREG O2 MS-TIE PROOF (register 1509, 2026-11-03)

Why: the user's scoping rule (2026-11-03): O2 is for the reads minimap2 could not place; everything else is the aligner's
decision. Register 1508 sized that set: a rival within 10 on `ms` (the score minimap2 ranks placements by), any AS, minus
the K = 0 wall (both margins 0, same splice shape, same chromosome). On KB3781 (Jim, the assembled individual) that is
52,061 molecules, 43,289 informative. Assemblers (StringTie2, FLAIR) take the primary record only, so on this set they
commit to an ms coin toss; O2's claim is that the PSV certificate assigns the reads PSVs decide and abstains on the rest.
Registers 1403/1404 settled the AS-gap class (primary right 99.6%) and 1409-1413 showed no read-level score decides a
near tie; this row tests the certificate on the ms-near-tie set itself.

Change: `copy_assign --tie-ms W` (commit below): the gate admits a molecule whose runner-up is within W of the best on ms
(`gate_tied` = AS-tied OR ms-tied); the assignments table gains `ms_best`, `ms_second`, `ms_margin` under the flag only.
Default runs are unchanged (checked below).

Sample and panel (fast-feedback rule): KB3781 fibroblast on mGorGor1, base catalog `~/runs/gorilla_KB3781/kb.fam.copies.tsv`
(494 families). Panel = the 40 families with most ms-near-tie informative reads (`bench/o2_contested/panel.py`, W 10,
copy window +- 50 kb): 70,559 of 73,467 read-family pairs (96%); `~/o2contested/panel/{families,regions,panel}.{txt,tsv}`.
Note the three top families (MCL263-265, 2 copies each, ~15k reads each) overlap one region and will share molecules.

Arms (one job, wall 3 h, 64 threads, `RUSTLE_WINDOW_SWEEP=1`, `--tied-prefilter`):
- IDENTITY: slim A119b assign with the new binary and no `--tie-ms`: all six `x.assign.*` tables byte-identical to r1503.
- BASE: panel assign, default gate (exact AS ties): the rows O2 judges today.
- MS10: panel assign, `--tie-ms 10`.

Truth: `bench/o2_truth/haplotype_truth.py` (register 1402-era arm 4, unchanged): reads and copies aligned to Jim's mat and
pat haplotypes (indexes job 4421505); a read's truth = the copies whose haplotype loci overlap one of its near hits
(NM <= best + 2) given a near-perfect best hit (NM <= max(2, 0.003 x aligned)); scored only when the truth is one locus.
This replaces the T = 3 fewest-edit rule discussed before it was found (the haplotypes are stronger evidence than edit
counts against the mosaic reference). Applied to MS10's contested rows (and BASE's for reference).

Readout per class (near tie on both scores same shape; spliced primary vs unspliced rival; other), four-way per molecule:
assigned-right, assigned-wrong, abstained, truth outside the catalog (O3's input). Beside it: minimap2's primary
(what StringTie2/FLAIR use) right/wrong on the same decidable molecules, and the 10-AS rule (`as_margin > 10`).

Decision rule (pre-declared):
1. IDENTITY holds (else fix before reading anything else).
2. MS10 judges >= 5x BASE's contested molecules on the panel (the set exists and the gate reaches it).
3. On decidable molecules, O2's wrong-call rate (assigned-wrong / assigned) <= the primary's wrong rate / 2, in the
   near-tie-both class; and O2 assigns >= 20% of the decidable molecules there (otherwise it abstains on everything and
   proves nothing).
4. Spliced-vs-unspliced class: reported; no bar (1403/1404 cover its AS-gap sibling).
If 2 fails: the census and the gate disagree, investigate before any claim. If 3 fails on the wrong rate: O2 does not
beat the aligner on its own uncertainty, record as FAIL. If 3 fails on the assigned fraction only: O2 is safe but silent,
record as INCONCLUSIVE and report what the abstentions lacked (decisive sites vs certificate p).
