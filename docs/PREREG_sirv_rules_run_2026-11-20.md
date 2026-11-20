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
