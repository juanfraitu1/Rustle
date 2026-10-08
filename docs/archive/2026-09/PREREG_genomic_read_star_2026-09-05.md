# PREREG — genomic read-star: the read against each candidate's LOCUS, splice-aware (O2-9d, §6fd), 2026-09-05
Defect (rows 706–708): the origin certificate on the spliced UNIT cannot separate isoform structure from origin —
every edit counted rejects minority isoforms of the right copy (LCR16u: 1,865 of 1,874 ambiguities), substitutions
only admits a wrong call in 20 anchors.
Design: for each molecule, minimap2 `-x splice` of its sequence against every candidate copy's GENOMIC span (the
unit's extent from copies.tsv, forward strand, fetched from --fasta), all hits kept (`-p 0.3`), `--eqx`. A retained
intron aligns through the intron, an alternative exon aligns to the locus, introns appear as `N` (not edits).
Columns = the read's positions where the candidates' aligned GENOMIC bases differ (pairwise certificates as in
§6fc, on shared columns). Origin certificate = every edit except `N` (X, I, D) against the best candidate's locus ~
Binomial(aligned, error_rate): sequence the locus does not have still rejects; structure does not. No junction term.
`copy_assign --read-star-genomic` (opt-in for this test; the unit form stays the default until the predictions hold).
Predictions (paired sets, same reads):
- P1 NPIP's 62 audited anchors: 0 wrong.
- P2 paired 35: MAPQ-60 placement agreement >= 99.8 % (unit form 99.90 %).
- P3 LCR16u: origin-rejected unit reads fall by >= 50 % (449 ambiguous of 654 under the unit form) and assigned
  reads exceed the unit form's 42.
- P4 NPIP: assigned >= 316 and MAPQ<60 assigned >= 115 (the unit form's numbers), with P1.
Fail rule: P1 fails => the genomic form is not adopted; P2 fails by more than 0.3 points => report and keep the
unit form; P3 fails => the LCR16u rejections are not isoform structure and the diagnosis of row 706 is wrong.
