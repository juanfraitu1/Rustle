# PREREG O2 GUARD: NEVER ASSIGN OFF THE PRIMARY'S LOCUS (register 1653, 2026-11-12)

Why: r1651 / r1652: when `copy_assign --tie-ms 10` assigns a read to a copy whose span does not overlap the read's
primary record, it is right 5 times and wrong 145 over three sets (KB3781, HG002 on its own haplotype, HG002 on CHM13).
Those sets found the rule, so they cannot test it.

Guard G (prototype as a post-filter, `bench/o2_contested/o2_moves.py`): an assigned row whose copy span does not
overlap the molecule's primary record becomes an abstention. The binary already computes this per row
(`primary_on_assigned`); a port would be one opt-in flag and is NOT part of this row.

Held-out sets (families not in any earlier panel of this line): from the existing rankings
`r1644_hg002_unspliced/panel80/panel.tsv` (HG002 on its paternal haplotype) and `r1651_kallisto_hg13/panel80/panel.tsv`
(HG002 on CHM13), the next 60 families after the 40 used, same 2 Mb cost guard, at least 1 target read. Per set:
`copy_assign --tie-ms 10` (binary `target_r17`), primary spans from the BAM, `haplotype_truth.py` + `strict_truth.py`
slack 0 with the r1644 haplotype indexes. Rows = decidable assigned rows; primary judged by record overlap.

Bars:
0. Power: >= 20 removed (moved) decidable rows over the two sets, else UNDERPOWERED.
1. Removed rows: wrong >= 4 x right, pooled.
2. On each set the wrong rate of assigned decidable rows is lower with G than without.
3. On each set G loses <= 2% of the correct assignments.
4. Wall <= 3 h.
Readout, no bar: A119b default-gate base run (`a1.assign.assignments.tsv`, 12,294 assigned rows): share of assigned
rows G would remove (no truth there); families holding the removed rows.

Verdicts: 0-3 pass = PASS: propose the opt-in flag port; the default is the user's call after it. 1, 2 or 3 fails =
FAIL: the 5 / 145 was those three panels. Expected: PASS on bars 1-2; bar 0 is the risk on same-individual data (moves
were 4% of assignments there).

Outputs: `/anvme/workspace/v138da13-rustle/r1653_o2_guard/`.

## RESULT r1653 (2026-11-12, job 4467118, 50 min; `r1653_o2_guard/{hgpat,hg13}/moves.txt`, `a1.moves.txt`)
Held-out: 60 families each. HG002 on its own haplotype: 72 regions, 19.7 Mb, 9,935 assigned rows; HG002 on CHM13: 43
regions, 13.0 Mb, 6,159 assigned rows.
| decidable assigned rows | HG002 on own haplotype | HG002 on CHM13 |
|---|---|---|
| kept on the primary's locus: right / wrong | 956 / 40 | 987 / 61 |
| moved off it (what G removes): right / wrong | 16 / 23 | 23 / 80 |
| wrong rate without G -> with G | 6.1% -> 4.0% | 12.3% -> 5.8% |
| correct assignments lost to G | 16 of 972 (1.6%) | 23 of 1,010 (2.3%) |
- Bar 0 PASS (142 removed rows). Bar 1 FAIL: removed rows 103 wrong / 39 right = 2.6 to 1 (bar 4 to 1; discovery sets
  were 29 to 1). Bar 2 PASS on both. Bar 3 FAIL on CHM13 (2.3% > 2%), pass on the own haplotype. Bar 4 PASS.
**FAIL as registered.** Held out, a move off the primary is still wrong more often than right (0.41 and 0.22 right),
and G halves the wrong rate across individuals, but it is no longer near-free: it discards 39 correct relocations.
The discovery figure (5 / 145) was inflated by single families (MCL189, MCL338).
A119b readout (default gate, exact AS ties, no truth): G would remove 10,216 of 12,293 assigned rows (83%). Under the
default gate the primary is one of several identical placements, so assigning off it is what O2 is for; the haplotype
truth cannot judge those rows (unique best hit required). **G is therefore not a candidate default and not a general
guard**: its evidence covers only ms near ties that the truth can decide.
Standing: not adopted, not ported. What survives for O2: moved near-tie rows deserve a stricter certificate than kept
ones (wrong 59-78% vs 4-6%), and three to five families carry most wrong moves on every set (MCL1256, MCL338, MCL586,
MCL189): the per-family fit check is the better-aimed next row.
